"""From an analysis to the flow cell runs it came from: the evidence ladder.

An analysis pools the reads of one or several runs, and most of what identifies a run is lost on the way. The
rungs are tried in this order and every link says which one made it:

1. ``read_tag``          the Dorado ``fn`` tag on the records of the workflow's tagged BAM: the name of the POD5
                         file, ``<flow cell>_<protocol run, 8 hex>_<acquisition, 8 hex>_<n>.pod5``. Tallied over at
                         most 100,000 primary records fetched from regions spread over the reference: runs and the
                         share of reads from each. High confidence.
2. ``bam_header_rg``     a Dorado ``@RG`` line that survived (``PU`` = flow cell, ``DS ... runid=``). High.
3. ``sibling_analysis``  a Space Ranger output takes the links of the workflow analysis made from the same reads.
4. ``read_id_sketch``    the read ids that start with a run of zeros, compared with the same sketch of each run.
5. ``path_rule``         a run folder beside the analysis (the same sample folder holds both) whose MinKNOW sample
                         id equals the sample name of the analysis. Medium.
   ``sample_name_chain`` the sample name alone, and only when exactly one run of the library carries it. Low.
6. ``manual``            stated by a person (``flowcellar analysis DIR --run FLOWCELL[@RUNID8]``). High.

Rungs 1, 2 and 6 need nothing but the analysis (``intrinsic_links``). Rungs 3 to 5 need the library
(``resolve``). A link whose run the library does not hold is kept, marked ``resolved: false``, and shown as
"run not in library"; it resolves by itself once that run is scanned.
"""
from __future__ import annotations

import os
import re
from dataclasses import dataclass, field

from flowcellar.util.errors import InputError, UsageError

# the name MinKNOW gives a POD5 file: flow cell, protocol run (8 hex), acquisition (8 hex), a counter
FN_RE = re.compile(r"^(?P<fc>[A-Z]{2,4}\d{4,6})_(?P<prot>[0-9a-f]{8})_(?P<acq>[0-9a-f]{8})_\d+\.pod5$")
_RG_RUNID = re.compile(r"\brunid=([0-9a-f]{8,64})")
_RUN_ARG = re.compile(r"^(?P<fc>[A-Za-z]{2,4}\d{4,6})(?:@(?P<run>[0-9a-fA-F]{8,12}))?$")
MAX_RECORDS = 100_000
N_REGIONS = 24
MIN_SKETCH_IDS = 8
CONFIDENCE_RANK = {"low": 0, "medium": 1, "high": 2, "certain": 3}


@dataclass(frozen=True)
class RunFact:
    """What the ladder needs to know about one run of the library."""
    run_key: str
    flow_cell_id: str
    run_short: str
    sample_id: str | None = None
    paths: tuple[str, ...] = ()                 # absolute paths of the run folder (it may lie under two roots)
    sketch: tuple[str, tuple[str, ...]] | None = None           # (prefix, ids) when the run bundle carries one
    acquisitions: tuple[str, ...] = ()


@dataclass
class Evidence:
    """What an analysis itself says about where its reads came from."""
    fn_runs: list[dict] = field(default_factory=list)       # from the fn tally
    rg_runs: list[dict] = field(default_factory=list)       # from @RG lines
    records: int = 0                                         # primary records looked at
    records_without_tag: int = 0
    note: str | None = None                                  # why the tally could not be made


def normalise_name(text: str | None) -> str:
    """A sample name for comparison: lower case, letters and digits only (``_Sample-1_A`` equals ``sample_1_a``)."""
    return re.sub(r"[^0-9a-z]+", "", (text or "").casefold())


# ---------------------------------------------------------------------------------------------------------
# rungs 1 and 2: what the BAM says
# ---------------------------------------------------------------------------------------------------------

def _regions(references: list[tuple[str, int]], n: int) -> list[tuple[str, int, int]]:
    """Up to ``n`` windows spread over the long references (all references when none is long)."""
    long_refs = [(r, ln) for r, ln in references if ln > 20_000_000] or sorted(references, key=lambda x: -x[1])[:n]
    out = []
    for r, ln in long_refs:
        for frac in (0.15, 0.4, 0.65, 0.9):
            start = int(ln * frac)
            out.append((r, start, min(ln, start + min(5_000_000, max(1, ln // 10)))))      # windows never overlap
    if len(out) > n:
        step = len(out) / n
        out = [out[int(i * step)] for i in range(n)]
    return out


def bam_evidence(bam: str, *, index: str | None = None, max_records: int = MAX_RECORDS, n_regions: int = N_REGIONS) -> Evidence:
    """Tally the ``fn`` tag of primary records fetched through the index, and read the ``@RG`` lines. Needs
    ``pysam``; without it, or without an index, the evidence says why it is empty."""
    ev = Evidence()
    try:
        import pysam  # noqa: PLC0415  the bam extra
    except ImportError:
        ev.note = 'the read tags of the BAM were not looked at: this needs pysam (pip install "flowcellar[bam]")'
        return ev
    index = index or bam + ".bai"
    if not os.path.isfile(index):
        ev.note = "the BAM has no index: its read tags were not looked at"
        return ev
    save = pysam.set_verbosity(0)
    try:
        with pysam.AlignmentFile(bam, "rb", index_filename=index, check_sq=False) as fh:
            for rg in fh.header.to_dict().get("RG", []):
                m = _RG_RUNID.search(rg.get("DS", "") or "")
                if rg.get("PU") and m and str(rg.get("PL", "")).upper() != "ILLUMINA":
                    ev.rg_runs.append({"flow_cell_id": rg["PU"], "acquisition_run_id_prefix": m.group(1)[:8]})
            refs = list(zip(fh.references, fh.lengths))
            if fh.mapped <= max_records:
                regions = [(r, 0, ln) for r, ln in refs]            # a small file: every record, so the shares are exact
                per = max_records
            else:
                regions = _regions(refs, n_regions)
                per = max(1, max_records // max(1, len(regions)))
            tally: dict[tuple, int] = {}
            st: dict[tuple, list] = {}
            for ref, a, b in regions:
                k = 0
                for rec in fh.fetch(ref, a, b):
                    if rec.is_secondary or rec.is_supplementary or rec.is_unmapped:
                        continue
                    k += 1
                    ev.records += 1
                    name = rec.get_tag("fn") if rec.has_tag("fn") else None
                    m = FN_RE.match(name) if isinstance(name, str) else None
                    if not m:
                        ev.records_without_tag += 1
                    else:
                        key = (m["fc"], m["prot"], m["acq"])
                        tally[key] = tally.get(key, 0) + 1
                        if rec.has_tag("st"):
                            s = str(rec.get_tag("st"))
                            lo_hi = st.setdefault(key, [s, s])
                            lo_hi[0], lo_hi[1] = min(lo_hi[0], s), max(lo_hi[1], s)
                    if k >= per:
                        break
    except (OSError, ValueError) as e:
        ev.note = f"the read tags of the BAM could not be read ({e})"
        return ev
    finally:
        pysam.set_verbosity(save)
    tagged = sum(tally.values())
    for key, n in sorted(tally.items(), key=lambda kv: (-kv[1], kv[0])):
        ev.fn_runs.append({"flow_cell_id": key[0], "protocol_run_id_prefix": key[1], "acquisition_run_id_prefix": key[2],
                           "reads": n, "share": n / tagged, "st_min": st.get(key, [None])[0], "st_max": st.get(key, [None, None])[1]})
    return ev


def _merge_acquisitions(fn_runs: list[dict]) -> list[dict]:
    """One entry per protocol run: a run with two acquisitions is still one run."""
    out: dict[tuple, dict] = {}
    for r in fn_runs:
        key = (r["flow_cell_id"], r["protocol_run_id_prefix"])
        e = out.get(key)
        if e is None:
            out[key] = dict(r)
        else:
            if r["reads"] > e["reads"]:
                e["acquisition_run_id_prefix"] = r["acquisition_run_id_prefix"]
            e["reads"] += r["reads"]
            e["share"] += r["share"]
    return sorted(out.values(), key=lambda r: -r["reads"])


def intrinsic_links(ev: Evidence) -> list[dict]:
    """The links an analysis carries by itself (rungs 1 and 2), in the layout of ``identity.analysis.links``."""
    links = []
    if ev.fn_runs:
        runs = _merge_acquisitions(ev.fn_runs)
        n = sum(r["reads"] for r in runs)
        for r in runs:
            links.append({"run_key": None, "flow_cell_id": r["flow_cell_id"], "protocol_run_id_prefix": r["protocol_run_id_prefix"],
                          "acquisition_run_id_prefix": r["acquisition_run_id_prefix"], "share": round(r["share"], 4),
                          "method": "read_tag", "confidence": "high", "resolved": False,
                          "evidence": f"file names stored on the reads: {r['reads']:,} of {n:,} records looked at name this run"})
    elif ev.rg_runs:
        for r in ev.rg_runs:
            links.append({"run_key": None, "flow_cell_id": r["flow_cell_id"], "protocol_run_id_prefix": None,
                          "acquisition_run_id_prefix": r["acquisition_run_id_prefix"], "share": None, "method": "bam_header_rg",
                          "confidence": "high", "resolved": False, "evidence": "read group of the basecaller in the BAM header"})
    return links


# ---------------------------------------------------------------------------------------------------------
# rung 6: stated by a person
# ---------------------------------------------------------------------------------------------------------

def parse_run_args(values) -> list[tuple[str, str | None]]:
    """``FLOWCELL`` or ``FLOWCELL@RUNID8`` -> (flow cell id, first 8 hex of the protocol run id or None)."""
    out = []
    for v in values or ():
        m = _RUN_ARG.match(v.strip())
        if not m:
            raise UsageError(f"--run takes FLOWCELL or FLOWCELL@RUNID8, got {v!r}", hint="for example --run TST00001@0a1b2c3d")
        out.append((m["fc"].upper(), m["run"].lower()[:8] if m["run"] else None))
    return out


def manual_links(stated: list[tuple[str, str | None]], runs: list[RunFact]) -> list[dict]:
    """Links as a person stated them. A flow cell with several runs in the library needs the run id as well."""
    links = []
    for fc, short in stated:
        if short is None:
            mine = [r for r in runs if r.flow_cell_id == fc]
            if len(mine) > 1:
                listed = "; ".join(f"{fc}@{r.run_short}" + (f" ({r.sample_id})" if r.sample_id else "") for r in mine)
                raise InputError(f"flow cell {fc} has {len(mine)} runs in this library: say which", hint="one of: " + listed)
            short = mine[0].run_short if mine else None
        links.append({"run_key": None, "flow_cell_id": fc, "protocol_run_id_prefix": short, "acquisition_run_id_prefix": None,
                      "share": None, "method": "manual", "confidence": "high", "resolved": False, "evidence": "stated by a person"})
    n = len(links)
    for ln in links:
        ln["share"] = round(1.0 / n, 4) if n == 1 else None
    return links


# ---------------------------------------------------------------------------------------------------------
# rungs 3 to 5 and the join with the library
# ---------------------------------------------------------------------------------------------------------

def _match(link: dict, runs: list[RunFact]) -> RunFact | None:
    mine = [r for r in runs if r.flow_cell_id == link.get("flow_cell_id")]
    short = link.get("protocol_run_id_prefix")
    if short:
        mine = [r for r in mine if r.run_short.startswith(short[:8]) or short.startswith(r.run_short)]
        return mine[0] if mine else None
    acq = link.get("acquisition_run_id_prefix")
    if acq:
        hit = [r for r in mine if any(a.startswith(acq) for a in r.acquisitions)]
        if len(hit) == 1:
            return hit[0]
    return mine[0] if len(mine) == 1 else None


def join(links: list[dict], runs: list[RunFact]) -> list[dict]:
    """Mark each link as resolved or not against the runs of the library, and give it its run key."""
    out = []
    for ln in links:
        ln = dict(ln)
        r = _match(ln, runs)
        ln["resolved"] = r is not None
        ln["run_key"] = r.run_key if r is not None else None
        if r is not None and not ln.get("protocol_run_id_prefix"):
            ln["protocol_run_id_prefix"] = r.run_short
        out.append(ln)
    return out


def sketch_links(sketch: dict | None, runs: list[RunFact]) -> list[dict]:
    """Rung 4. The share of the analysis' sketch found in each run's sketch gives the share of reads per run.
    Two sketches compare on the longer of their two prefixes."""
    ids = set((sketch or {}).get("ids") or ())
    if len(ids) < MIN_SKETCH_IDS:
        return []
    mine_prefix = sketch.get("prefix") or "0000"
    hits = []
    for r in runs:
        if not r.sketch:
            continue
        prefix, run_ids = r.sketch[0], set(r.sketch[1])
        common = mine_prefix if len(mine_prefix) >= len(prefix) else prefix
        a = {i for i in ids if i.startswith(common)}
        b = {i for i in run_ids if i.startswith(common)}
        found = len(a & b)
        if len(a) >= MIN_SKETCH_IDS and found >= max(3, 0.05 * len(a)):
            hits.append((r, found, len(a), len(b)))
    total = sum(h[1] for h in hits)
    return [{"run_key": r.run_key, "flow_cell_id": r.flow_cell_id, "protocol_run_id_prefix": r.run_short,
             "acquisition_run_id_prefix": None, "share": round(found / total, 4), "method": "read_id_sketch",
             "confidence": "high", "resolved": True,
             "evidence": f"read names: {found} of the {n_a} read ids of the analysis that start with {len(max(mine_prefix, r.sketch[0], key=len))} zeros are reads of this run"}
            for r, found, n_a, _ in sorted(hits, key=lambda h: -h[1])]


def name_links(names: list[str], analysis_path: str, runs: list[RunFact], *, max_levels: int = 6) -> list[dict]:
    """Rung 5. ``names`` are the sample names of the analysis (the workflow's alias, Space Ranger's sample id)."""
    wanted = {normalise_name(n) for n in names if n and normalise_name(n)}
    if not wanted:
        return []
    same = [r for r in runs if normalise_name(r.sample_id) in wanted]
    if not same:
        return []
    ancestors = []
    d = os.path.abspath(analysis_path)
    for _ in range(max_levels):
        d = os.path.dirname(d)
        if not d or d == os.path.dirname(d):
            break
        ancestors.append(d)
    beside = [r for r in same if any(os.path.dirname(p.rstrip("/")) in ancestors for p in r.paths)]
    if beside:
        chosen, method, conf = beside, "path_rule", "medium"
        words = "the run folder lies beside the analysis and MinKNOW's sample id equals the sample name of the analysis"
    elif len(same) == 1:
        chosen, method, conf = same, "sample_name_chain", "low"
        words = "by name only: one run of this library has a MinKNOW sample id equal to the sample name of the analysis"
    else:
        return []
    return [{"run_key": r.run_key, "flow_cell_id": r.flow_cell_id, "protocol_run_id_prefix": r.run_short,
             "acquisition_run_id_prefix": None, "share": round(1.0 / len(chosen), 4) if len(chosen) == 1 else None,
             "method": method, "confidence": conf, "resolved": True, "evidence": words} for r in chosen]


def resolve(*, intrinsic: list[dict], runs: list[RunFact], sibling_links: list[dict] | None = None,
            sketch: dict | None = None, names: list[str] | None = None, analysis_path: str = "",
            stated: list[tuple[str, str | None]] | None = None) -> list[dict]:
    """The links of one analysis inside a library: the first rung that yields any, joined with the runs."""
    if stated:
        return join(manual_links(stated, runs), runs)
    if intrinsic:
        return join(intrinsic, runs)
    if sibling_links:
        inherited = []
        for ln in sibling_links:
            ln = dict(ln)
            ln["evidence"] = "taken from the workflow analysis of the same reads (" + (ln.get("evidence") or ln.get("method") or "") + ")"
            ln["method"] = "sibling_analysis"
            inherited.append(ln)
        return join(inherited, runs)
    by_sketch = sketch_links(sketch, runs)
    if by_sketch:
        return by_sketch
    return name_links(names or [], analysis_path, runs)


def worst_confidence(links: list[dict]) -> str | None:
    if not links:
        return None
    return min((ln.get("confidence") or "low" for ln in links), key=lambda c: CONFIDENCE_RANK.get(c, 0))
