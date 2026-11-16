"""The BAM source: one exact pass over a Dorado BAM or an aligned BAM.

What a pass gives, per run key and sample label found in the records: reads and bases, the length histogram with
bases, the mean-quality histogram, the length by quality grid, and for aligned files mapped and unmapped reads,
identity, mapping quality, and the secondary and supplementary records. **Primary records count reads once**:
records with flag 0x100 or 0x800 are only counted.

From the same records, in the same pass (``alnacc``): strand, the CIGAR totals and the mismatches of the primary
mapped records, aligned length, soft-clipped share, three 2-D grids, reads and aligned bases per contig, aligned bases
per window along a genome, read counts per transcript and coverage along the transcript body for a transcriptome
(``reference`` tells the two apart from the header), and from the basecaller's tags ``st``, ``du`` and ``pt``: output
per ten minutes of the clock, length, quality, aligned length, identity and measured speed per clock hour, and the
poly(A) tail length.

How it is read:

* an indexed BAM is cut into tasks of about equal compressed size at record starts taken from the index
  (``bamindex``); worker processes seek to a task's first record and stop at the next task's first record. Every
  aggregate is an integer, so the result does not depend on the number of workers;
* a BAM without an index is read from start to end by one process with htslib threads;
* the parent process never opens the file: the header, the index and every record are read in workers, and a pass
  whose workers make no progress for ``stall_timeout_s`` is abandoned.

Routing (per record, strongest first): the ``fn`` tag (name of the POD5 file: flow cell and run id), the ``RG`` tag
(through the header's ``@RG PU`` and ``DS runid=``, or its run-id prefix when the header has no read groups), the
single run all read groups of the header name. A record with none of these is *unrouted*: the pass keeps it under
no run and the link of the file decides where it belongs.

Tag hazards handled here: a tag wrapped by a FASTQ round trip as ``CO:Z:<tag>:<type>:<value>``; ``qs`` as float or
integer; read groups whose value ends in an alias the header does not list; read names prefixed by a cDNA
pre-processor; secondary records without a sequence. ``ts`` exists with two types (aligner and basecaller) and is
never read.
"""
from __future__ import annotations

import hashlib
import math
import os
import re
import struct
import time
import warnings
from datetime import datetime

import numpy as np

from flowcellar.model import runkey
from flowcellar.readdata import bamindex
from flowcellar.readdata._pool import beat, run_tasks
from flowcellar.readdata.acc import ReadAcc
from flowcellar.readdata.alnacc import (
    C_D,
    C_I,
    C_LEN,
    C_M,
    C_MAPPED,
    C_MAPQ,
    C_NM,
    C_Q,
    C_SA,
    NCOL,
    AlignAcc,
    TagAcc,
)
from flowcellar.readdata.reference import Layout, gene_index
from flowcellar.util.errors import InputError, MissingExtraError

BGZF_EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")
TASK_BYTES = 64 << 20
MAX_GROUPS = 1024
OTHER_LABEL = "__other__"
MAX_SPLIT_PARENTS = 2_000_000
BATCH = 4096 * NCOL * 8                 # bytes of the rows a group holds before it adds them up
FLIPS = 64
_ROW = struct.Struct(f"{NCOL}d").pack    # one row of alnacc.COLS; packing is three times faster than extending an array
#: per-transcript read counts are kept while groups x reference sequences stays below this (4 bytes each in a worker)
REF_COUNT_CELLS = 64_000_000
# A probe only decides which runs a file may touch (the pass routes every read). 500 records from each of 24 places
# find a second flow cell that holds one read in a thousand with a chance above 99.99 %, and cost about 20 MB of
# reading per file instead of 150 MB; the specification's bound is 100,000 records and 10 s.
PROBE_RECORDS = 12_000
PROBE_REGIONS = 24
PROBE_SECONDS = 10.0

_RG_RE = re.compile(r"^([0-9a-f]{8})[0-9a-f-]{24,32}_(.+?@v\d+(?:\.\d+)*)(?:_(.+))?$")
_RG_HASH_RE = re.compile(r"-[0-9A-F]{6,8}$")         # the suffix samtools adds to a read group id when it merges files
_RG_SHORT_RE = re.compile(r"^([0-9a-f]{8})[0-9a-f-]{24,32}(?:_|$)")
_BARCODE_RE = re.compile(r"(?:^|_)(barcode\d+[A-Za-z]?)$")
_PYCHOPPER_NAME_RE = re.compile(r"^\d+:\d+\|")
_MODEL_TIER_RE = re.compile(r"_(fast|hac|sup)(?:@|_|$)")
_NAN = float("nan")
_ERR_LUT = np.power(10.0, -np.arange(256) / 10.0)


def pysam_module():
    """pysam, or one clear line that names the extra to install."""
    try:
        import pysam
    except ImportError:
        raise MissingExtraError("pysam", "bam", "Reading BAM files") from None
    pysam.set_verbosity(0)
    return pysam


def has_eof_marker(path: str) -> bool:
    """True when the file ends with the empty BGZF block every complete BAM ends with."""
    try:
        with open(path, "rb") as fh:
            fh.seek(0, os.SEEK_END)
            if fh.tell() < len(BGZF_EOF):
                return False
            fh.seek(-len(BGZF_EOF), os.SEEK_END)
            return fh.read() == BGZF_EOF
    except OSError:
        return False


# ------------------------------------------------------------------------------------------ header
def norm_label(value: str | None) -> str | None:
    """A sample label as the records give it: ``SQK-NBD114-24_barcode07`` -> ``barcode07``; an alias stays;
    ``unclassified`` is no label."""
    if not value:
        return None
    v = value.strip()
    if not v or v.lower() == "unclassified" or v.lower().endswith("_unclassified"):
        return None
    m = _BARCODE_RE.search(v)
    return m.group(1) if m else v


def model_tier(model: str | None) -> str | None:
    m = _MODEL_TIER_RE.search(model or "")
    return m.group(1) if m else None


def _tokens(cl: str) -> list[str]:
    """A command line as Dorado and samtools write it into @PG (arguments may be double-quoted)."""
    return [t.strip('"') for t in re.findall(r'"[^"]*"|\S+', cl or "")]


def header_facts(hd: dict) -> dict:
    """What the header says about runs, samples and how the reads were made. Plain data, safe to store."""
    rg: dict[str, list] = {}
    runs: dict[str, int] = {}
    models: set[str] = set()
    experiments: set[str] = set()
    for g in hd.get("RG", []) or []:
        rid = str(g.get("ID") or "")
        ds = str(g.get("DS") or "")
        m = re.search(r"runid=([0-9a-fA-F-]{8,})", ds)
        runid = m.group(1).lower() if m else None
        pu = str(g.get("PU") or "")
        key = None
        if runid and runkey.FLOW_CELL_RE.match(pu):
            key = f"{pu}_{runid.replace('-', '')[:8]}"
        elif runid:
            key = "?_" + runid.replace("-", "")[:8]
        mm = re.search(r"basecall_model=(\S+)", ds)
        if mm:
            models.add(mm.group(1))
        me = re.search(r"experiment_id=(\S+)", ds)
        if me:
            experiments.add(me.group(1))
        label = norm_label(g.get("al") or g.get("SM"))
        rg[rid] = [key, label, mm.group(1) if mm else None]
        if key:
            runs[key] = runs.get(key, 0) + 1
    programs = []
    dorado: dict = {}
    aligner: dict = {}
    pychopper = False
    for p in hd.get("PG", []) or []:
        pn = str(p.get("PN") or p.get("ID") or "")
        cl = str(p.get("CL") or "")
        toks = _tokens(cl)
        programs.append({"name": pn, "version": p.get("VN"), "options": sorted({t for t in toks if t.startswith("--")})[:40]})
        low = pn.lower()
        if low == "dorado":
            sub = toks[1] if len(toks) > 1 else ""
            dorado.setdefault("version", p.get("VN"))
            dorado.setdefault("commands", []).append(sub)
            if "--min-qscore" in toks:
                try:
                    dorado["min_qscore"] = float(toks[toks.index("--min-qscore") + 1])
                except (ValueError, IndexError):
                    pass
            if sub in ("basecaller", "duplex"):
                dorado["no_trim"] = "--no-trim" in toks
                dorado["estimate_poly_a"] = "--estimate-poly-a" in toks
                dorado["modified_bases"] = any(t.startswith("--modified-bases") for t in toks)
            if sub == "demux":
                dorado["demux"] = True
                dorado["demux_no_trim"] = "--no-trim" in toks
            if "--barcode-both-ends" in toks:
                dorado["barcode_both_ends"] = True
            if "--kit-name" in toks:
                try:
                    dorado["kit_name"] = toks[toks.index("--kit-name") + 1]
                except IndexError:
                    pass
        elif low in ("minimap2", "winnowmap", "dorado_aligner", "ngmlr", "bwa", "vg"):
            aligner = {"name": pn, "version": p.get("VN"), "preset": next((toks[i + 1] for i, t in enumerate(toks[:-1]) if t in ("-ax", "-x")), None),
                       "copies_comments": "-y" in toks}
            if any("full_length" in os.path.basename(t) for t in toks if not t.startswith("-")):
                aligner["input_full_length"] = True
        if "pychopper" in low or "pychopper" in cl.lower():
            pychopper = True
    hdl = hd.get("HD") or {}
    return {"sort_order": hdl.get("SO"), "n_sq": len(hd.get("SQ", []) or []), "n_rg": len(rg), "rg": rg, "runs": runs,
            "models": sorted(models), "experiments": sorted(experiments), "programs": programs[:12], "dorado": dorado,
            "aligner": aligner, "pychopper_program": pychopper}


# ----------------------------------------------------------------------------------------- routing
class _Group:
    """One run key and one sample label: the accumulators, and the rows of the records not yet added to them
    (``alnacc.COLS`` per primary record, packed as float64 one row after the other)."""
    __slots__ = ("acc", "aln", "tags", "R", "parents", "ref_counts")

    def __init__(self):
        self.acc = ReadAcc()
        self.aln = AlignAcc()
        self.tags = TagAcc()
        self.parents: set[int] | None = set()
        self.ref_counts = True
        self.R = bytearray()

    def flush(self, layout: Layout | None = None):
        if not len(self.R):
            return
        R = np.ascontiguousarray(np.frombuffer(self.R, dtype=np.float64).reshape(-1, NCOL).T)     # columns by records
        mapped = R[C_MAPPED]
        cols = np.where(mapped > 0, R[C_M] + R[C_I] + R[C_D], -1.0)
        self.acc.add(R[C_LEN], R[C_Q], mapped, cols, R[C_NM], R[C_MAPQ], R[C_SA])
        self.tags.add(R)
        self.aln.add(R, layout, ref_counts=self.ref_counts)
        self.R = bytearray()


class _State:
    """Everything one task accumulates: groups by (run key, sample label), caches, counters."""

    def __init__(self, facts: dict, layout: Layout | None = None):
        self.layout = layout
        self.ref_cells = 0                  # reference sequences x groups that count reads per transcript
        self.clock: dict[str, float] = {}   # start time text -> ten-minute bin of the clock
        self.rg = facts["rg"]
        self.rg_ids = sorted(self.rg, key=len, reverse=True)
        self.default = next(iter(facts["runs"])) if len(facts["runs"]) == 1 else None
        self.groups: dict[tuple, _Group] = {}
        self.cache: dict[tuple, _Group] = {}
        self.fn_keys: dict[str, str | None] = {}
        self.rg_info: dict[str, tuple] = {}
        self.models: set[str] = set()
        self.tags: set[str] = set()
        self.records = 0
        self.primary = 0
        self.names_checked = 0
        self.pychopper_names = 0
        self.truncated = False
        self.error: str | None = None
        self.finished = False               # the range was read to its end
        self.last = 0

    def _rg(self, rg: str) -> tuple:
        info = self.rg_info.get(rg)
        if info is not None:
            return info
        key = label = model = None
        full = rg
        if rg not in self.rg:
            rg = _RG_HASH_RE.sub("", rg)    # a read group value can end in a hash suffix the header does not have
        rid = rg if rg in self.rg else next((r for r in self.rg_ids if rg.startswith(r + "_")), None)
        if rid is not None:
            key, label, model = self.rg[rid]
            tail = rg[len(rid):].lstrip("_")
            if tail:                        # the record's value carries an alias the header does not list
                label = norm_label(tail)
        m = _RG_RE.match(rg)
        if m:
            model = model or m.group(2)
            if rid is None:
                label = norm_label(m.group(3))
        if key is None:
            ms = _RG_SHORT_RE.match(rg)
            if ms:
                key = "?_" + ms.group(1)
        if model:
            self.models.add(model)
        info = self.rg_info[full] = (key, label)
        return info

    def resolve(self, fn, rg, bc) -> _Group:
        key = None
        if fn:
            if fn in self.fn_keys:
                key = self.fn_keys[fn]
            else:
                p = runkey.parse_minknow_name(fn)
                key = self.fn_keys[fn] = p.run_key if p else None
                if len(self.fn_keys) > 500_000:
                    self.fn_keys.clear()
        label = None
        if rg:
            rkey, label = self._rg(rg)
            if key is None:
                key = rkey
        if key is None:
            key = self.default
        if bc:
            label = norm_label(bc)
        gid = (key, label)
        g = self.groups.get(gid)
        if g is None:
            if len(self.groups) >= MAX_GROUPS:
                gid = (key, OTHER_LABEL)
                g = self.groups.get(gid)
            if g is None:
                g = self.groups[gid] = _Group()
                lay = self.layout
                if lay is not None and lay.kind == "transcriptome":
                    self.ref_cells += lay.n
                    g.ref_counts = self.ref_cells <= REF_COUNT_CELLS
        if len(self.cache) > 200_000:
            self.cache.clear()
        self.cache[(fn, rg, bc)] = g
        return g

    def clock_bin(self, text: str, key: str) -> float:
        """The ten-minute bin of the clock a start time falls in (floor(epoch seconds / 600)), -1 when unreadable."""
        try:
            tb = float(datetime.fromisoformat(text).timestamp() // 600)
        except (ValueError, TypeError, OverflowError):
            tb = -1.0
        if len(self.clock) > 100_000:
            self.clock.clear()
        self.clock[key] = tb
        return tb

    def result(self, seconds: float) -> dict:
        groups, parents, aln, tagacc = {}, {}, {}, {}
        for gid, g in self.groups.items():
            g.flush(self.layout)
            groups[gid] = g.acc.to_state()
            aln[gid] = g.aln.to_state() if g.aln.mapped else None
            tagacc[gid] = g.tags.to_state()
            parents[gid] = None if g.parents is None else np.fromiter(g.parents, dtype=np.uint64, count=len(g.parents))
        return {"groups": groups, "aln": aln, "tagacc": tagacc, "parents": parents, "records": self.records, "primary": self.primary,
                "names_checked": self.names_checked, "pychopper_names": self.pychopper_names, "tags": sorted(self.tags),
                "models": sorted(self.models), "truncated": self.truncated, "error": self.error,
                "finished": self.finished, "last": self.last, "seconds": seconds}


def _quality_from_bases(r) -> float:
    """Mean quality of a read from its base qualities (mean error probability, as the basecaller defines it)."""
    qual = r.query_qualities
    if qual is None or not len(qual):
        return _NAN
    p = float(_ERR_LUT[np.frombuffer(qual, dtype=np.uint8)].mean())
    return -10.0 * math.log10(p) if p > 0 else 60.0


def _scan(b, end, st: _State, max_primary: int | None = None, stop_wall: float | None = None) -> None:
    """Read records from the current position until the virtual offset ``end`` (None: end of file)."""
    tell = b.tell
    nxt = b.__next__
    cache = st.cache
    clock = st.clock
    row = _ROW
    layout = st.layout
    genome = layout is not None and layout.kind == "genome"
    win_first = layout.win_first_list if genome else None
    W = layout.window if genome else 1
    # Looking a tag up walks the tags of the record, and asking "is it there" before fetching it walks them twice.
    # A tag is on every record of a file or on none, nearly always: so each tag is fetched straight away while it
    # keeps turning up (one walk), and asked for first once it was missing (one walk again). A file that changes its
    # mind more than FLIPS times about a tag is asked first from then on.
    t_fn = t_rg = t_bc = t_co = t_qs = t_st = t_du = t_pt = True
    f_fn = f_rg = f_bc = f_co = f_qs = f_st = f_du = f_pt = 0
    n = primary = 0
    try:
        while True:
            if end is not None and tell() >= end:
                st.finished = True
                break
            try:
                r = nxt()
            except StopIteration:
                st.finished = True
                break
            n += 1
            flag = r.flag
            has = r.has_tag
            get = r.get_tag
            if t_fn:
                try:
                    fn = get("fn")
                except KeyError:
                    fn, t_fn = None, False
            elif has("fn"):
                fn = get("fn")
                f_fn += 1
                t_fn = f_fn < FLIPS
            else:
                fn = None
            if t_rg:
                try:
                    rg = get("RG")
                except KeyError:
                    rg, t_rg = None, False
            elif has("RG"):
                rg = get("RG")
                f_rg += 1
                t_rg = f_rg < FLIPS
            else:
                rg = None
            if t_bc:
                try:
                    bc = get("BC")
                except KeyError:
                    bc, t_bc = None, False
            elif has("BC"):
                bc = get("BC")
                f_bc += 1
                t_bc = f_bc < FLIPS
            else:
                bc = None
            if t_co:
                try:
                    co = get("CO")
                except KeyError:
                    co, t_co = None, False
            elif has("CO"):
                co = get("CO")
                f_co += 1
                t_co = f_co < FLIPS
            else:
                co = None
            wrapped_q = when = None
            if co is not None:
                if isinstance(co, str) and len(co) > 5 and co[2] == ":" and co[4] == ":":
                    t = co[:2]                       # the first tag of a FASTQ comment, wrapped by the aligner
                    if t == "BC":
                        if bc is None:
                            bc = co[5:]
                    elif t == "qs":
                        wrapped_q = co[5:]
                    elif t == "st":
                        when = co[5:]
                    elif t == "fn":
                        if fn is None:
                            fn = co[5:]
                    elif t == "RG":
                        if rg is None:
                            rg = co[5:]
            g = cache.get((fn, rg, bc))
            if g is None:
                g = st.resolve(fn, rg, bc)
            if flag & 0x900:                         # secondary or supplementary: counted, never read as a read
                if flag & 0x100:
                    g.acc.secondary += 1
                else:
                    g.acc.supplementary += 1
                continue
            primary += 1
            length = r.query_length
            if not length:
                length = r.infer_read_length() or 0
            if t_qs:
                try:
                    q = get("qs")
                except KeyError:
                    q, t_qs = None, False
            elif has("qs"):
                q = get("qs")
                f_qs += 1
                t_qs = f_qs < FLIPS
            else:
                q = None
            if q is not None:
                try:
                    q = float(q)
                except (TypeError, ValueError):
                    q = _NAN
            elif wrapped_q is not None:
                try:
                    q = float(wrapped_q)
                except ValueError:
                    q = _NAN
            else:
                q = _quality_from_bases(r)
                if q == q:
                    g.acc.q_computed += 1
            # the basecaller's own tags: when the read started, how long it took, its poly(A) tail
            if t_st:
                try:
                    when = get("st")
                except KeyError:
                    t_st = False
            elif has("st"):
                when = get("st")
                f_st += 1
                t_st = f_st < FLIPS
            if when:
                key = when[:15] + when[-6:]          # the ten minutes and the offset from UTC
                tb = clock.get(key)
                if tb is None:
                    tb = st.clock_bin(when, key)
            else:
                tb = -1.0
            if t_du:
                try:
                    du = get("du")
                except KeyError:
                    du, t_du = -1.0, False
            elif has("du"):
                du = get("du")
                f_du += 1
                t_du = f_du < FLIPS
            else:
                du = -1.0
            if t_pt:
                try:
                    pt = get("pt")
                except KeyError:
                    pt, t_pt = -2.0, False
            elif has("pt"):
                pt = get("pt")
                f_pt += 1
                t_pt = f_pt < FLIPS
            else:
                pt = -2.0
            if type(du) is not float:                # a tag of another type than the basecaller writes
                du = float(du) if isinstance(du, int) else -1.0
            if type(pt) is not int:
                pt = float(pt) if isinstance(pt, float) and pt == pt else -2.0
            if flag & 4:
                g.R += row(length, q, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0, 0.0, tb, du, pt, 0.0)
            else:
                cs = r.get_cigar_stats()[0]
                nm = cs[10]
                if nm == 0 and not has("NM"):
                    nm = -1
                m = cs[0] + cs[7] + cs[8]
                skip = cs[3]
                rid = r.reference_id
                start = r.reference_start
                walked = 0.0
                if skip and genome:                  # a spliced record: its aligned columns are not one stretch
                    wf = win_first[rid]
                    if wf >= 0 and start // W != (start + m + cs[2] + skip - 1) // W:
                        g.aln.walk(layout, wf, start, r.cigartuples)
                        walked = 1.0
                g.R += row(length, q, 1.0, m, cs[1], cs[2], skip, cs[4], cs[5], nm, r.mapping_quality,
                           1.0 if has("SA") else 0.0, 1.0 if flag & 16 else 0.0, rid, start, tb, du, pt, walked)
            if has("pi"):                            # a piece of a read the basecaller split: count the parent once
                g.acc.split_children += 1
                if g.parents is not None:
                    g.parents.add(int.from_bytes(hashlib.blake2b(str(get("pi")).encode(), digest_size=8).digest(), "little"))
                    if len(g.parents) > MAX_SPLIT_PARENTS:
                        g.parents = None
            if len(g.R) >= BATCH:
                g.flush(layout)
            if st.names_checked < 200:
                st.names_checked += 1
                if _PYCHOPPER_NAME_RE.match(r.query_name or ""):
                    st.pychopper_names += 1
                if st.names_checked <= 50:
                    st.tags.update(t for t, _ in r.get_tags())
            if max_primary is not None and primary >= max_primary:
                break
            if not n & 4095:
                beat(4096)
                if stop_wall is not None and time.time() >= stop_wall:
                    break
    except OSError as e:                             # a cut or damaged block: keep what was read
        st.truncated = True
        st.error = str(e)[:300]
    st.records += n
    st.primary += primary
    try:
        st.last = int(tell())
    except Exception:       # noqa: BLE001  (a closed or broken handle has no position)
        pass


# ----------------------------------------------------------------------------------------- workers
_OPEN: dict[tuple, tuple] = {}
_LAYOUT: dict[str, Layout] = {}


def _layout(path: str, b) -> Layout:
    """The reference layout of a file, built once per process from its header."""
    lay = _LAYOUT.get(path)
    if lay is None:
        lay = _LAYOUT[path] = Layout(b.get_reference_name, b.lengths)
    return lay


def _open(path: str, threads: int = 1):
    """One handle per file and process, with the facts of its header."""
    key = (path, threads)
    hit = _OPEN.get(key)
    if hit is None:
        pysam = pysam_module()
        try:
            if has_eof_marker(path):
                b = pysam.AlignmentFile(path, "rb", check_sq=False, threads=max(1, threads))
            else:               # cut short or still being written: htslib reads it only without extra threads
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")         # pysam warns; the pass records the state instead
                    b = pysam.AlignmentFile(path, "rb", check_sq=False, ignore_truncation=True)
        except (OSError, ValueError) as e:
            raise InputError(f"{path}: cannot be read as a BAM file ({e})") from None
        hit = _OPEN[key] = (b, header_facts(b.header.to_dict()), int(b.tell()))
    return hit


def _close_all() -> None:
    for b, _, _ in _OPEN.values():
        try:
            b.close()
        except Exception:       # noqa: BLE001
            pass
    _OPEN.clear()
    _LAYOUT.clear()


def plan_task(args: tuple) -> dict:
    return plan(*args)


def plan(path: str, task_bytes: int = TASK_BYTES, use_index: bool = True) -> dict:
    """The header facts, the index totals and the tasks of one BAM. Runs in a worker process."""
    b, facts, first = _open(path)
    size = os.path.getsize(path)
    eof_ok = has_eof_marker(path)
    out = {"facts": facts, "first": first, "size": size, "eof_ok": eof_ok, "index": None, "tasks": [(first, None)]}
    lay = _layout(path, b)
    out["reference"] = lay.describe()
    out["genes"] = gene_index(b.references) if lay.kind == "transcriptome" else None
    ipath = bamindex.find_index(path) if use_index else None
    if ipath and b.is_bam:
        info = {"path": ipath, "stale": bamindex.index_is_stale(path, ipath)}
        try:
            idx = bamindex.read_index(ipath)
            info.update(kind=idx.kind, n_ref=idx.n_ref, records=idx.records, mapped=idx.mapped,
                        unmapped_placed=idx.unmapped_placed, no_coordinate=idx.no_coordinate, cut_points=int(len(idx.offsets)))
            if eof_ok and not info["stale"]:
                out["tasks"] = bamindex.cut_tasks(first, idx.offsets, size, task_bytes)
        except (ValueError, OSError) as e:
            info["error"] = str(e)[:200]
        out["index"] = info
    return out


def task(args: tuple) -> dict:
    """One range of one file. ``args`` = (path, start, end, threads, max_primary, stop_wall)."""
    path, start, end, threads, max_primary, stop_wall = args
    t0 = time.time()
    if stop_wall is not None and t0 >= stop_wall:
        return {"groups": {}, "aln": {}, "tagacc": {}, "parents": {}, "records": 0, "primary": 0, "names_checked": 0, "pychopper_names": 0, "tags": [],
                "models": [], "truncated": False, "error": None, "finished": False, "last": start or 0, "seconds": 0.0,
                "range": (start, end), "skipped": True}
    b, facts, first = _open(path, threads)
    b.seek(start if start is not None else first)
    st = _State(facts, _layout(path, b))
    _scan(b, end, st, max_primary, stop_wall)
    res = st.result(time.time() - t0)
    res["range"] = (start, end)
    return res


# ------------------------------------------------------------------------------------------ merging
class Merged:
    """Results of tasks, added up."""

    def __init__(self):
        self.groups: dict[tuple, ReadAcc] = {}
        self.aln: dict[tuple, AlignAcc] = {}            # groups with at least one mapped primary record
        self.tagacc: dict[tuple, TagAcc] = {}
        self.parents: dict[tuple, set | None] = {}
        self.records = self.primary = 0
        self.names_checked = self.pychopper_names = 0
        self.tags: set[str] = set()
        self.models: set[str] = set()
        self.truncated = False
        self.errors: list[str] = []
        self.unfinished: list[tuple] = []
        self.bytes_done = 0
        self.seconds = 0.0

    def add(self, res: dict, size: int) -> None:
        for gid, state in res["groups"].items():
            acc = ReadAcc.from_state(state)
            if gid in self.groups:
                self.groups[gid].merge(acc)
            else:
                self.groups[gid] = acc
            for store, cls, state in ((self.aln, AlignAcc, res["aln"].get(gid)), (self.tagacc, TagAcc, res["tagacc"].get(gid))):
                if state is None:
                    continue
                part = cls.from_state(state)
                if gid in store:
                    store[gid].merge(part)
                else:
                    store[gid] = part
            new = res["parents"].get(gid)
            if gid not in self.parents:
                self.parents[gid] = set()
            if new is None or self.parents[gid] is None:
                self.parents[gid] = None
            else:
                self.parents[gid].update(int(x) for x in new)
                if len(self.parents[gid]) > MAX_SPLIT_PARENTS:
                    self.parents[gid] = None
        self.records += res["records"]
        self.primary += res["primary"]
        self.names_checked += res["names_checked"]
        self.pychopper_names += res["pychopper_names"]
        self.tags.update(res["tags"])
        self.models.update(res["models"])
        self.seconds += res["seconds"]
        start, end = res["range"]
        c0 = (start or 0) >> 16
        c1 = size if end is None else end >> 16
        if res["truncated"]:
            self.truncated = True
            self.errors.append(res["error"] or "the file ends early")
        if res["finished"]:
            self.bytes_done += max(c1 - c0, 0)
        else:
            self.unfinished.append((start, end))
            self.bytes_done += max(min(res["last"] >> 16, c1) - c0, 0)


def _stride_order(n: int) -> list[int]:
    """0..n-1 so that any prefix is spread over the whole range (a time-limited pass covers the file evenly)."""
    bits = max(1, (n - 1).bit_length())
    return sorted(range(n), key=lambda i: int(format(i, f"0{bits}b")[::-1], 2))


def scan_bam(path: str, *, jobs: int = 0, threads: int = 4, deadline: float | None = None, stall_timeout_s: float = 120.0,
             task_bytes: int = TASK_BYTES, inprocess: bool = False, use_index: bool = True, progress=None) -> dict:
    """One exact pass over one BAM. Returns ``{"plan", "merged": Merged, "state", "coverage", "errors", "timings"}``.

    ``deadline`` is a ``time.monotonic()`` value; a pass that reaches it stops and says what it covered.
    """
    pysam_module()
    t0 = time.time()
    stop_wall = None if deadline is None else time.time() + max(deadline - time.monotonic(), 0.0)
    workers = jobs if jobs and jobs > 0 else min(8, os.cpu_count() or 1)
    pl = run_tasks(plan_task, [(path, task_bytes, use_index)], workers=1, stall_timeout_s=stall_timeout_s,
                   inprocess=inprocess)["results"]
    if not pl:
        return {"plan": None, "merged": Merged(), "state": "error", "coverage": {"kind": "truncated", "records": 0, "records_total": None, "fraction": 0.0},
                "errors": [{"code": "stalled", "message": "the file did not answer in time"}], "timings": {"seconds": time.time() - t0, "workers": 0}}
    pl = pl[0]
    errors: list[dict] = []
    idx = pl["index"]
    if idx and idx.get("stale"):
        errors.append({"code": "index_stale", "message": "the index is older than the file and was not used"})
    if idx and idx.get("error"):
        errors.append({"code": "index_unreadable", "message": idx["error"]})
    if idx and pl["facts"]["n_sq"] > 100_000:
        workers = min(workers, 4)               # every worker holds the header: a transcriptome has 500,000 contigs

    def run(tasks: list[tuple], n_workers: int, th: int) -> tuple[Merged, bool]:
        order = _stride_order(len(tasks)) if len(tasks) > 1 else [0]
        args = [(path, tasks[i][0], tasks[i][1], th, None, stop_wall) for i in order]
        out = run_tasks(task, args, workers=n_workers, stall_timeout_s=stall_timeout_s, inprocess=inprocess, progress=progress)
        if inprocess:
            _close_all()
        m = Merged()
        for res in out["results"]:
            m.add(res, pl["size"])
        return m, out["stalled"]

    tasks = pl["tasks"]
    parallel = len(tasks) > 1
    merged, stalled = run(tasks, min(workers, len(tasks)) if parallel else 1, 1 if parallel else threads)
    total = idx.get("records") if idx else None
    timed_out = stop_wall is not None and time.time() >= stop_wall and bool(merged.unfinished)
    if parallel and not stalled and not timed_out and not merged.truncated and total is not None and merged.records != total:
        # the index does not describe this file: read it again from start to end, which needs no index
        errors.append({"code": "index_mismatch", "message": f"the index counts {total:,} records, a pass cut from it found "
                       f"{merged.records:,}; the file was read again without the index"})
        merged, stalled = run([(pl["first"], None)], 1, threads)
        total = None
    if stalled:
        errors.append({"code": "stalled", "message": f"no progress for {stall_timeout_s:.0f} s; the pass was abandoned"})
    if (pl.get("reference") or {}).get("kind") == "transcriptome" and len(merged.groups) * pl["facts"]["n_sq"] > REF_COUNT_CELLS:
        # the same rule whatever the tasks were: too many groups to count reads per transcript for each of them
        for a in merged.aln.values():
            a.ref_counts_dropped, a.ref_n = True, None
    fraction = min(merged.bytes_done / max(pl["size"] - (pl["first"] >> 16), 1), 1.0) if merged.unfinished else 1.0
    if stalled:
        state, kind = "error", "byte_ranges"
    elif merged.truncated or not pl["eof_ok"]:
        state, kind = "truncated", "truncated"
        if not merged.errors and not pl["eof_ok"]:
            merged.errors.append("the file has no end-of-file block: it was cut short or is still being written")
    elif merged.unfinished:
        state, kind = "partially_read", "byte_ranges"
    else:
        state, kind, fraction = "ok", "full", 1.0
    for msg in merged.errors[:3]:
        errors.append({"code": "truncated", "message": msg})
    seconds = time.time() - t0
    return {"plan": pl, "merged": merged, "state": state,
            "coverage": {"kind": kind, "records": int(merged.records), "records_total": total, "fraction": round(float(fraction), 6)},
            "errors": errors,
            "timings": {"seconds": round(seconds, 3), "mb_per_s": round(pl["size"] / 1e6 / max(seconds, 1e-6), 1),
                        "workers": min(workers, len(tasks)) if parallel else 1, "tasks": len(tasks)}}


def probe(path: str, max_primary: int = PROBE_RECORDS, regions: int = PROBE_REGIONS, seconds: float = PROBE_SECONDS) -> dict:
    """Identity evidence of one BAM: header facts and the run keys and sample labels of at most ``max_primary``
    primary records taken from up to ``regions`` places of the file (through the index when there is one).
    Runs in the calling process: the link step calls it inside a worker."""
    t0 = time.time()
    pl = plan(path, task_bytes=max(os.path.getsize(path) // max(regions, 1), 1 << 16))
    tasks = pl["tasks"]
    if len(tasks) > regions:
        step = len(tasks) / regions
        tasks = [tasks[int(i * step)] for i in range(regions)]
    per = max(max_primary // max(len(tasks), 1), 1)
    stop = time.time() + seconds
    m = Merged()
    for start, end in tasks:
        m.add(task((path, start, end, 1, per, stop)), pl["size"])
        if time.time() >= stop:
            break
    _close_all()
    runs: dict[str, int] = {}
    labels: dict[str, int] = {}
    unrouted = mapped = 0
    for (key, label), acc in m.groups.items():
        if key is None:
            unrouted += acc.reads
        else:
            runs[key] = runs.get(key, 0) + acc.reads
        labels[label or ""] = labels.get(label or "", 0) + acc.reads
        mapped += acc.mapped
    f = pl["facts"]
    return {"kind": "bam", "header_runs": f["runs"], "runs": runs, "unrouted": unrouted, "labels": labels, "primary": m.primary,
            "records": m.records, "mapped": mapped, "pychopper_names": m.pychopper_names, "names_checked": m.names_checked,
            "tags": sorted(m.tags), "models": sorted(set(f["models"]) | m.models), "sort_order": f["sort_order"], "n_sq": f["n_sq"],
            "n_rg": f["n_rg"], "dorado": f["dorado"], "aligner": f["aligner"], "pychopper_program": f["pychopper_program"],
            "experiments": f["experiments"], "index": pl["index"], "eof_ok": pl["eof_ok"], "regions": len(tasks),
            "seconds": round(time.time() - t0, 3)}
