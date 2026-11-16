"""Attaching read-data files to runs.

Identifiers inside a file always beat its path. The ladder, strongest first, as far as this build goes:

====================  ==========  ==========================================================================
method                confidence  evidence
====================  ==========  ==========================================================================
``run_dir``           certain     the file lies inside the run folder
``bam_header_rg``     certain     every read group of the header names one run (``PU`` + ``DS runid=``)
``read_tag``          certain     ``fn`` and ``RG`` tags of the probed records; gives the share of each run
``summary_column``    certain     file names and run ids in a read summary
``sample_sheet``      certain     ``flow_cell_id`` (and ``experiment_id``) of a basecaller sample sheet
``sample_name_chain`` high        same sample name as a file linked by content, inside one experiment folder
``path_rule``         medium      a ``[[link.rules]]`` entry with ``match``; shown, not trusted until confirmed
====================  ==========  ==========================================================================

A file is opened only to take its fingerprint and to probe it (header, 12,000 records from 24 places, 10 s), both in a
worker process; both results are cached under ``probes_dir``, so linking the same files again opens nothing.
Files with equal fingerprint and size are one input with several paths.
"""
from __future__ import annotations

import fnmatch
import os
import re

from flowcellar.contracts import Candidate, LinkedInput, RunIdentity
from flowcellar.model import partial as partial_model
from flowcellar.model import runkey
from flowcellar.readdata import record
from flowcellar.readdata._pool import run_tasks
from flowcellar.util import atomic
from flowcellar.util.fingerprint import fingerprint as file_fingerprint
from flowcellar.util.fingerprint import hash_obj

PROBE_VERSION = 4
#: the rules by which files are attached to runs. A bump makes a library link its files again at the next scan (no file
#: is opened for that: fingerprints and probes are cached). 2: a file without a variant name is the primary variant in
#: the sample-name chain.
LINK_VERSION = 2
BAM_ROLES = ("dorado.bam", "align.genome_bam", "align.transcriptome_bam")
TABLE_PASS_ROLES = ("polya.read_table",)              # large tables with a streamed pass of their own (readtable)
PASS_ROLES = BAM_ROLES + ("dorado.summary",) + TABLE_PASS_ROLES
COLLECTED_ROLES = ("pychopper.stats", "quant.salmon", "quant.oarfish", "quant.kallisto", "basecaller.sample_sheet")
KNOWN_ROLES = PASS_ROLES + COLLECTED_ROLES
#: versions of the section builders that read a pass: a bump changes the partial id, so only that pass is redone
#: (the one bump: the BAM pass now also collects strand, CIGAR totals, grids, depth windows, per-transcript
#: counts, the clock tables of st / du and the poly(A) tail of pt)
PASS_MODULES = {"yield.basecalled": "0.1.0", "len": "0.1.0", "qual": "0.1.0", "aln": "0.2.0", "aln.windows": "0.1.0", "aln.time": "0.1.0",
                "dyn.basecalled": "0.1.0", "rna.polya": "0.1.0", "rna.tx": "0.1.0", "rna.coverage": "0.1.0"}
CONTENT_METHODS = ("run_dir", "bam_header_rg", "read_tag", "summary_column")
INPROCESS_ENV = "FLOWCELLAR_READDATA_INPROCESS"


def inprocess() -> bool:
    """Tests and debugging: do the file reading of READ DATA in the calling process."""
    return os.environ.get(INPROCESS_ENV, "") not in ("", "0")


# ------------------------------------------------------------------------------------------- names
def canonical(name: str | None, config: dict | None) -> str | None:
    """A sample name after the library's alias map (recorded alias = canonical name)."""
    if not name:
        return None
    alias = ((config or {}).get("samples") or {}).get("alias") or {}
    return alias.get(name, name)


def is_unclassified(name: str | None, config: dict | None) -> bool:
    if not name:
        return False
    pats = ((config or {}).get("samples") or {}).get("unclassified") or ["unclassified"]
    return any(fnmatch.fnmatchcase(name, p) for p in pats)


def experiment_dir_of(c: Candidate) -> str:
    """The experiment folder of a file: what discovery found, else the folder above ``analysis`` or above a
    basecaller output tree, else the folder of the file."""
    if c.experiment_dir:
        return c.experiment_dir
    parts = c.path.split(os.sep)
    for i in range(len(parts) - 2, 0, -1):
        if parts[i] == "analysis" or parts[i] in ("bam_pass", "bam_fail", "fastq_pass", "fastq_fail"):
            return os.sep.join(parts[:i]) or os.sep
    return os.path.dirname(c.path)


# ------------------------------------------------------------------------------------------ probing
def _probe_path(probes_dir: str, fp: dict) -> str:
    v = fp["value"][:32]
    return os.path.join(probes_dir, v[:2], v + ".json")


def _probe_family(role: str) -> str:
    if role in BAM_ROLES:
        return "bam"
    if role == "basecaller.sample_sheet":
        return "sample_sheet"
    return "summary" if role == "dorado.summary" else "table"


def _probe_summary(path: str) -> dict:
    """Run keys named by the first rows of a read summary."""
    from flowcellar.readdata.tables import summary_columns
    cols = summary_columns(path)
    runs: dict[str, int] = {}
    fcol = next((c for c in ("filename_pod5", "filename", "filename_fast5", "filename_bam") if c in cols), None)
    rows = 0
    with open(path, "rb") as fh:
        fh.readline()
        for raw in fh.read(1 << 18).splitlines()[:-1][:2000]:
            p = raw.decode("utf-8", errors="replace").split("\t")
            if len(p) != len(cols):
                continue
            rows += 1
            key = None
            if fcol:
                pn = runkey.parse_minknow_name(p[cols.index(fcol)])
                key = pn.run_key if pn else None
            if key is None and "run_id" in cols:
                rid = p[cols.index("run_id")].replace("-", "").lower()
                if re.fullmatch(r"[0-9a-f]{8,}", rid):
                    key = "?_" + rid[:8]
            if key:
                runs[key] = runs.get(key, 0) + 1
    return {"kind": "summary", "columns": cols, "runs": runs, "rows": rows,
            "basecalled": "sequence_length_template" in cols and "mean_qscore_template" in cols}


def identify_task(args: tuple) -> dict:
    """Fingerprint and probe of one file, with both caches. ``args`` = (path, role, probes_dir). Worker process."""
    path, role, probes_dir = args
    out = {"path": path, "fingerprint": None, "probe": None, "error": None}
    try:
        fp = file_fingerprint(path)
        out["fingerprint"] = fp
        family = _probe_family(role)
        cache = _probe_path(probes_dir, fp) if probes_dir else None
        if cache and os.path.isfile(cache):
            try:
                d = atomic.read_json(cache)
                if d.get("version") == PROBE_VERSION and d.get("family") == family:
                    out["probe"] = d["probe"]
                    return out
            except (OSError, ValueError):
                pass
        if family == "bam":
            from flowcellar.readdata import bam
            probe = bam.probe(path)
        elif family == "sample_sheet":
            from flowcellar.readdata.sample_sheet import read_sample_sheet
            s = read_sample_sheet(path)
            probe = {"kind": "sample_sheet", "flow_cells": s["flow_cells"], "experiments": s["experiments"], "run_ids": s["run_ids"],
                     "rows": [{k: r.get(k) for k in ("flow_cell_id", "experiment_id", "barcode", "alias", "sample_id", "protocol_run_id")}
                              for r in s["rows"]]}
        elif family == "summary":
            probe = _probe_summary(path)
        else:
            probe = {"kind": "table"}
        out["probe"] = probe
        if cache:
            atomic.atomic_write_json(cache, {"format": "flowcellar.probe", "version": PROBE_VERSION, "family": family,
                                             "fingerprint": fp, "probe": probe})
    except Exception as e:      # noqa: BLE001  one unreadable file must not stop the linking of the others
        out["error"] = f"{type(e).__name__}: {e}"[:300]
    return out


def _stat_cache_path(probes_dir: str) -> str:
    return os.path.join(probes_dir, "stat.json.gz")


def identify(candidates: list[Candidate], probes_dir: str | None) -> dict[str, dict]:
    """path -> {"fingerprint", "probe", "error"} for every candidate, reading only what no cache answers."""
    stat_cache: dict = {}
    if probes_dir and os.path.isfile(_stat_cache_path(probes_dir)):
        try:
            stat_cache = (atomic.read_json(_stat_cache_path(probes_dir)) or {}).get("files", {})
        except (OSError, ValueError):
            stat_cache = {}
    out: dict[str, dict] = {}
    todo = []
    for c in candidates:
        hit = stat_cache.get(c.path)
        if hit and hit[0] == c.size and hit[1] == c.mtime_ns and probes_dir:
            fp = {"method": hit[2], "value": hit[3]}
            cache = _probe_path(probes_dir, fp)
            try:
                d = atomic.read_json(cache) if os.path.isfile(cache) else None
            except (OSError, ValueError):
                d = None
            if d and d.get("version") == PROBE_VERSION and d.get("family") == _probe_family(c.role):
                out[c.path] = {"path": c.path, "fingerprint": fp, "probe": d["probe"], "error": None}
                continue
        todo.append((c.path, c.role, probes_dir))
    if todo:
        res = run_tasks(identify_task, todo, workers=min(8, os.cpu_count() or 1), stall_timeout_s=120.0, inprocess=inprocess())
        for r in res["results"]:
            out[r["path"]] = r
        for path, _, _ in todo:
            out.setdefault(path, {"path": path, "fingerprint": None, "probe": None,
                                  "error": "the file did not answer in time" if res["stalled"] else "not read"})
        if probes_dir:
            changed = False
            by_path = {c.path: c for c in candidates}
            for path, _, _ in todo:
                fp = out[path].get("fingerprint")
                if fp:
                    c = by_path[path]
                    stat_cache[path] = [c.size, c.mtime_ns, fp["method"], fp["value"]]
                    changed = True
            if changed:
                atomic.atomic_write_json(_stat_cache_path(probes_dir), {"format": "flowcellar.stat_cache", "files": stat_cache})
    return out


# ----------------------------------------------------------------------------------------- the ladder
class _Runs:
    def __init__(self, runs):
        self.by_key: dict[str, RunIdentity] = {r.run_key: r for r in runs}
        self.by_short: dict[str, list[str]] = {}
        self.by_fc: dict[str, list[str]] = {}
        self.dirs: list[tuple[str, str]] = []
        for r in runs:
            self.by_short.setdefault(r.run_short[:8], []).append(r.run_key)
            self.by_fc.setdefault(r.flow_cell_id, []).append(r.run_key)
            if r.run_dir:
                self.dirs.append((os.path.join(os.path.normpath(r.run_dir), ""), r.run_key))

    def resolve(self, key: str) -> str | None:
        """A key found in a file -> a run key. A key that gives the run id only is completed when exactly one known
        run has that id; a full key the library does not know is kept (reads of a run it has no folder for)."""
        if key.startswith("?_"):
            hits = self.by_short.get(key[2:], [])
            return hits[0] if len(hits) == 1 else None
        return key if runkey.is_run_key(key) else None

    def inside(self, path: str) -> str | None:
        for d, k in self.dirs:
            if path.startswith(d):
                return k
        return None

    def under(self, folder: str) -> list[str]:
        base = os.path.join(os.path.normpath(folder), "")
        return sorted(k for d, k in self.dirs if d.startswith(base))


def _fmt_tally(t: dict[str, int]) -> str:
    return ", ".join(f"{k} {v:,}" for k, v in sorted(t.items(), key=lambda kv: (-kv[1], kv[0]))[:6])


def _by_content(c: Candidate, probe: dict | None, R: _Runs) -> tuple[list[str], dict, dict] | None:
    """(run keys, link, shares) from the identifiers inside a file or from its place in a run folder; None when the
    file names no run."""
    k = R.inside(c.path)
    if k:
        return [k], {"method": "run_dir", "confidence": "certain", "evidence": "the file lies inside the run folder", "rule": None}, {}
    if not probe:
        return None
    kind = probe.get("kind")
    if kind == "bam":
        tally: dict[str, int] = {}
        unknown: dict[str, int] = {}
        for key, n in (probe.get("runs") or {}).items():
            rk = R.resolve(key)
            if rk:
                tally[rk] = tally.get(rk, 0) + n
            else:
                unknown[key] = n
        header = sorted({R.resolve(k2) for k2 in (probe.get("header_runs") or {})} - {None})
        if len(header) == 1 and set(tally) <= set(header) and c.role == "dorado.bam":
            n_rg = sum((probe.get("header_runs") or {}).values())
            return header, {"method": "bam_header_rg", "confidence": "certain", "rule": None,
                            "evidence": f"every read group of the header names {header[0]} ({n_rg} read group{'s' if n_rg != 1 else ''}: PU and DS runid)"}, {}
        if tally:
            total = sum(tally.values())
            ev = f"fn and RG tags of {probe.get('primary', 0):,} probed records: {_fmt_tally(tally)}"
            if probe.get("unrouted"):
                ev += f"; {probe['unrouted']:,} records name no run"
            if unknown:
                ev += f"; run ids not in the library: {_fmt_tally({k2[2:]: v for k2, v in unknown.items()})}"
            return sorted(tally), {"method": "read_tag", "confidence": "certain", "evidence": ev, "rule": None}, \
                {k2: round(v / total, 6) for k2, v in tally.items()}
        if header:                                    # an aligned file that kept the read groups but no read tags
            return header, {"method": "bam_header_rg", "confidence": "certain", "rule": None,
                            "evidence": f"the read groups of the header name {', '.join(header)}"}, {}
        return None
    if kind == "summary":
        tally = {}
        for key, n in (probe.get("runs") or {}).items():
            rk = R.resolve(key)
            if rk:
                tally[rk] = tally.get(rk, 0) + n
        if tally:
            total = sum(tally.values())
            return sorted(tally), {"method": "summary_column", "confidence": "certain", "rule": None,
                                   "evidence": f"file names and run ids of the first {probe.get('rows', 0):,} rows: {_fmt_tally(tally)}"}, \
                {k2: round(v / total, 6) for k2, v in tally.items()}
        return None
    if kind == "sample_sheet":
        keys: list[str] = []
        how: list[str] = []
        weakest = "certain"
        exp_dir = experiment_dir_of(c)
        for rid in probe.get("run_ids") or []:
            short = rid.replace("-", "")[:8].lower()
            for rk in R.by_short.get(short, []):
                if rk not in keys:
                    keys.append(rk)
                    how.append(f"protocol_run_id column names {rk}")
        for fc in probe.get("flow_cells") or []:
            cands = list(R.by_fc.get(fc, []))
            if not cands or any(k2 in keys for k2 in cands):
                continue
            exps = {r.get("experiment_id") for r in probe.get("rows") or [] if r.get("flow_cell_id") == fc and r.get("experiment_id")}
            if len(cands) > 1:
                same_exp = [k2 for k2 in cands if (R.by_key[k2].facts or {}).get("experiment_id") in exps]
                near = [k2 for k2 in cands if k2 in R.under(exp_dir)]
                narrowed = same_exp or near
                if narrowed:
                    cands = narrowed
            if len(cands) == 1:
                how.append(f"flow_cell_id column names {fc}")
            else:
                # a flow cell id alone does not say which loading: shown, and trusted only when the experiment agrees
                agree = all((R.by_key[k2].facts or {}).get("experiment_id") in exps for k2 in cands)
                weakest = "high" if agree and weakest == "certain" else "medium" if not agree else weakest
                how.append(f"flow_cell_id column names {fc}, which the library holds {len(cands)} runs of"
                           + (" in the same experiment" if agree else "; the sheet does not say which"))
            keys.extend(k2 for k2 in cands if k2 not in keys)
        if keys:
            return keys, {"method": "sample_sheet", "confidence": weakest, "evidence": "; ".join(how), "rule": None}, {}
    return None


def _file_name_key(c: Candidate, R: _Runs) -> str | None:
    p = runkey.parse_minknow_name(os.path.basename(c.path))
    return p.run_key if p and p.run_key in R.by_key else None


def _path_matches(c: Candidate, pattern: str) -> bool:
    rel = c.relpath or c.path
    return fnmatch.fnmatchcase(rel, pattern) or fnmatch.fnmatchcase(c.path, pattern) or fnmatch.fnmatchcase(os.path.basename(c.path), pattern)


def link(candidates, runs, *, config=None, probes_dir=None) -> list[LinkedInput]:
    """See ``flowcellar.readdata.api.link``."""
    config = config or {}
    cands = [c for c in candidates if c.role in KNOWN_ROLES]
    info = identify(cands, probes_dir)
    R = _Runs(list(runs))
    # ---- equal fingerprints (and sizes) are one input with several paths
    uniq: dict[tuple, LinkedInput] = {}
    order: list[LinkedInput] = []
    probes: dict[int, dict | None] = {}
    for c in cands:
        i = info.get(c.path) or {}
        fp = i.get("fingerprint")
        if not fp:
            li = LinkedInput(candidate=c, fingerprint={"method": "none", "value": hash_obj([c.path, c.size, c.mtime_ns])}, run_keys=[],
                             link={"method": None, "confidence": None, "evidence": f"not readable: {i.get('error') or 'no answer'}", "rule": None})
            order.append(li)
            probes[id(li)] = None
            continue
        key = (fp["method"], fp["value"], c.size)
        if key in uniq:
            first = uniq[key]
            if c.path != first.candidate.path and c.path not in first.alt_paths:
                first.alt_paths.append(c.path)
            continue
        li = LinkedInput(candidate=c, fingerprint=fp, run_keys=[], link={"method": None, "confidence": None, "evidence": "", "rule": None})
        li.variant = c.variant
        uniq[key] = li
        order.append(li)
        probes[id(li)] = i.get("probe")
        if i.get("error"):
            li.link["evidence"] = f"not readable: {i['error']}"
    # ---- names: the sample of a file is the capture of its rule, else the one label its records carry
    for li in order:
        c, probe = li.candidate, probes[id(li)] or {}
        name = c.sample
        if name is None and probe.get("kind") == "bam":
            labels = [k for k, v in (probe.get("labels") or {}).items() if k and v]
            if len(labels) == 1 and not (probe.get("labels") or {}).get("", 0):
                name = labels[0]
        li.sample = canonical(name, config)
        li.needs = list(record.needs_of(c.role))
    # ---- ranks 1 to 5: location and identifiers inside the file
    for li in order:
        got = _by_content(li.candidate, probes[id(li)], R)
        if got:
            li.run_keys, li.link, li.shares = got
            li.link.setdefault("confirmed", False)
    # ---- rank 8: the sample-name chain inside one experiment folder
    anchors: dict[tuple, list[LinkedInput]] = {}
    for li in order:
        if li.run_keys and li.link.get("method") in CONTENT_METHODS and li.sample and li.candidate.role in PASS_ROLES:
            anchors.setdefault((experiment_dir_of(li.candidate), li.sample), []).append(li)
    sheet_alias: dict[tuple, list[tuple[LinkedInput, list[str]]]] = {}
    for li in order:
        probe = probes[id(li)] or {}
        if li.candidate.role == "basecaller.sample_sheet" and li.run_keys and li.link.get("confidence") in ("certain", "high"):
            for row in probe.get("rows") or []:
                for nm in {canonical(row.get("alias"), config), canonical(row.get("sample_id"), config)} - {None}:
                    keys = [k for k in li.run_keys if not row.get("flow_cell_id") or k.split("_")[0] == row["flow_cell_id"]]
                    sheet_alias.setdefault((experiment_dir_of(li.candidate), nm), []).append((li, keys))
    chain_rules = [r for r in ((config.get("link") or {}).get("rules") or []) if r.get("via") == "sample"]
    # a file whose name carries no variant is the primary variant (``[variants] primary``): the main alignment of a
    # sample chains to its main basecalled file, not only a named variant to the file of the same name
    primary_variant = (config.get("variants") or {}).get("primary")

    def version(v):
        return primary_variant if v is None else v
    for li in order:
        if li.run_keys or not li.sample or li.link.get("evidence", "").startswith("not readable"):
            continue
        c = li.candidate
        exp = experiment_dir_of(c)
        rule = next((r.get("name") for r in chain_rules if not r.get("roles") or c.role in r["roles"]), None)
        same = [a for a in anchors.get((exp, li.sample), []) if a is not li and (c.role not in BAM_ROLES or version(a.variant) == version(li.variant))]
        if same:
            keys = sorted({k for a in same for k in a.run_keys})
            best = sorted(same, key=lambda a: (record.READS_PREFERENCE.index(a.candidate.role) if a.candidate.role in record.READS_PREFERENCE else 9, a.candidate.path))[0]
            li.run_keys = keys
            li.link = {"method": "sample_name_chain", "confidence": "high", "rule": rule, "confirmed": False,
                       "evidence": f"same sample name '{li.sample}' as {os.path.basename(best.candidate.path)} ({best.link['method']}) in the same experiment folder"}
            continue
        sheets = sheet_alias.get((exp, li.sample), [])
        if sheets:
            keys = sorted({k for _, ks in sheets for k in ks})
            sheet = sheets[0][0]
            if keys:
                li.run_keys = keys
                li.link = {"method": "sample_name_chain", "confidence": "high", "rule": rule, "confirmed": False,
                           "evidence": f"sample name '{li.sample}' is listed in the sample sheet {os.path.basename(sheet.candidate.path)} of the same experiment folder"}
    # ---- rank 9: configured path rules, and the run key MinKNOW and Dorado write into file names
    path_rules = [r for r in ((config.get("link") or {}).get("rules") or []) if r.get("match")]
    for li in order:
        if li.run_keys or li.link.get("evidence", "").startswith("not readable"):
            continue
        c = li.candidate
        for r in path_rules:
            if r.get("roles") and c.role not in r["roles"]:
                continue
            if not _path_matches(c, r["match"]):
                continue
            if r.get("run") and r["run"] in R.by_key:
                keys = [r["run"]]
            elif r.get("flow_cell"):
                keys = sorted(R.by_fc.get(r["flow_cell"], []))
            else:
                keys = R.under(experiment_dir_of(c))
            if keys:
                li.run_keys = keys
                li.link = {"method": "path_rule", "confidence": "medium", "rule": r.get("name") or r["match"], "confirmed": bool(r.get("confirmed")),
                           "evidence": f"path matches the rule '{r.get('name') or r['match']}' ({r['match']}); "
                                       + (f"{len(keys)} run folders lie in its experiment folder" if len(keys) > 1 else f"links to {keys[0]}")}
                break
        if not li.run_keys:
            k = _file_name_key(c, R)
            if k:
                li.run_keys = [k]
                li.link = {"method": "path_rule", "confidence": "medium", "rule": "file name", "confirmed": False,
                           "evidence": f"the file name carries the flow cell and run id of {k}; nothing inside the file names a run"}
        if not li.run_keys and (c.rule or "").startswith("named on the command line") and len(R.by_key) == 1:
            li.run_keys = list(R.by_key)
            li.link = {"method": "manual", "confidence": "certain", "rule": None, "confirmed": True,
                       "evidence": "named on the command line for this run; nothing inside the file names a run"}
        if not li.run_keys and not li.link.get("evidence"):
            li.link["evidence"] = "nothing inside the file names a run, no file of the same sample is linked, and no link rule matches its path"
    # ---- scope, source of record, the pass each input needs
    for li in order:
        probe = probes[id(li)] or {}
        splittable = li.link.get("method") in ("read_tag", "summary_column") and li.candidate.role in PASS_ROLES
        if len(li.run_keys) > 1:
            li.scope = "shared" if splittable or li.candidate.role == "basecaller.sample_sheet" else "sample_all_runs"
        if li.candidate.role in TABLE_PASS_ROLES and li.fingerprint.get("method") != "none":
            # the table pass has its own module versions: a change there re-reads the tables only
            from flowcellar.readdata import readtable
            li.partial_id = partial_model.partial_id(li.fingerprint, li.candidate.role, readtable.MODULES, readtable.config_hash())
        elif li.candidate.role in PASS_ROLES and li.fingerprint.get("method") != "none":
            li.partial_id = partial_model.partial_id(li.fingerprint, li.candidate.role, PASS_MODULES, "")
        li.link["labels"] = sorted(k for k in (probe.get("labels") or {}) if k)[:50] if probe.get("kind") == "bam" else []
    record.choose(order, config)
    return order
