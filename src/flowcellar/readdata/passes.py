"""One pass over one input file -> one partial.

A partial of READ DATA holds, per run key and sample label found in the file, the fragments ``len``, ``qual`` and
(for a file that can hold alignments) ``aln``, each with one cell in the bundle's aggregate JSON. It is a pure
function of the file's content: nothing of the library's configuration or of the link is baked into its numbers.

``fragments[].run_key`` is a run key, ``"?_<8 hex>"`` when the records give a run id but no flow cell, or null for
records that name no run. ``fragments[].sample`` is the label the records carry (an alias or ``barcodeNN``), or null.
A read summary with a pass/fail column writes ``dims: ["filter"]`` and one cell per filter value; a BAM writes
``dims: []`` and its filter follows from where the basecaller put the file.
"""
from __future__ import annotations

import os
import time

from flowcellar.contracts import LinkedInput
from flowcellar.model import partial as partial_model
from flowcellar.readdata import link as L
from flowcellar.readdata._pool import run_tasks
from flowcellar.readdata.acc import ReadAcc
from flowcellar.util.timeutil import from_epoch


def _new(inp: LinkedInput) -> dict:
    c = inp.candidate
    return partial_model.new_partial(role=c.role, input_path=c.path, fp=inp.fingerprint, size=c.size,
                                     mtime=from_epoch(c.mtime_ns / 1e9) if c.mtime_ns else None, modules=dict(L.PASS_MODULES),
                                     config_hash="", sample=inp.sample, variant=inp.variant)


def _fragments(p: dict, key, label, cells: dict[str, dict], flt: str | None = None) -> None:
    for mod, cell in cells.items():
        existing = next((f for f in p["fragments"] if f["run_key"] == key and f["sample"] == label and f["module"] == mod), None)
        if flt is None:
            p["fragments"].append(partial_model.fragment(key, mod, sample=label, dims=[], cells=[{"at": [], **cell}]))
        elif existing is None:
            p["fragments"].append(partial_model.fragment(key, mod, sample=label, dims=["filter"], cells=[{"at": [flt], **cell}]))
        else:
            existing["cells"].append({"at": [flt], **cell})


def run_pass(inp: LinkedInput, *, jobs: int = 0, threads: int = 4, deadline: float | None = None, stall_timeout_s: float = 120.0,
             progress=None, inprocess: bool = False, task_bytes: int | None = None) -> dict:
    """The partial of one input (not written). ``state`` is ok, truncated, partially_read or error."""
    if inp.candidate.role == "dorado.summary":
        return _summary_pass(inp, deadline=deadline, stall_timeout_s=stall_timeout_s, inprocess=inprocess)
    if inp.candidate.role in L.TABLE_PASS_ROLES:        # an annotated read table: streamed, by readtable
        from flowcellar.readdata import readtable
        return readtable.run_pass(inp, deadline=deadline, stall_timeout_s=stall_timeout_s, inprocess=inprocess)
    return _bam_pass(inp, jobs=jobs, threads=threads, deadline=deadline, stall_timeout_s=stall_timeout_s, progress=progress,
                     inprocess=inprocess, task_bytes=task_bytes)


def _bam_pass(inp: LinkedInput, *, jobs, threads, deadline, stall_timeout_s, progress, inprocess, task_bytes) -> dict:
    from flowcellar.readdata import alnpass, bam
    kw = {"task_bytes": task_bytes} if task_bytes else {}
    res = bam.scan_bam(inp.candidate.path, jobs=jobs, threads=threads, deadline=deadline, stall_timeout_s=stall_timeout_s,
                       inprocess=inprocess, progress=progress, **kw)
    p = _new(inp)
    p["state"] = res["state"]
    p["coverage"] = res["coverage"]
    p["errors"] = res["errors"]
    p["timings"] = res["timings"]
    pl = res["plan"]
    if pl is None:
        return p
    m = res["merged"]
    facts = pl["facts"]
    aligned = facts["n_sq"] > 0 or any(a.mapped for a in m.groups.values())
    runs: dict[str, int] = {}
    labels: dict[str, int] = {}
    unrouted = secondary = supplementary = 0
    for (key, label), acc in sorted(m.groups.items(), key=lambda kv: (str(kv[0][0]), str(kv[0][1]))):
        cells = acc.cells(aligned=aligned)
        parents = m.parents.get((key, label))
        cells["len"]["scalars"]["split_parents"] = None if parents is None else len(parents)
        alnpass.extend(p, cells, m, (key, label), pl)         # strand, CIGAR, grids, depth, over time, poly(A)
        _fragments(p, key, label, cells)
        if key is None:
            unrouted += acc.reads
        else:
            runs[key] = runs.get(key, 0) + acc.reads
        labels[label or ""] = labels.get(label or "", 0) + acc.reads
        secondary += acc.secondary
        supplementary += acc.supplementary
    alnpass.finish(p, m, pl)
    p["routing"] = {"runs": runs, "unrouted": unrouted, "default_run": next(iter(facts["runs"])) if len(facts["runs"]) == 1 else None}
    idx = pl["index"] or {}
    p["facts"] = {"kind": "bam", "records": int(m.records), "primary": int(m.primary), "secondary": int(secondary),
                  "supplementary": int(supplementary), "labels": labels, "sort_order": facts["sort_order"], "n_sq": facts["n_sq"],
                  "n_rg": facts["n_rg"], "header_runs": facts["runs"], "experiments": facts["experiments"],
                  "models": sorted(set(facts["models"]) | m.models), "dorado": facts["dorado"], "aligner": facts["aligner"],
                  "programs": facts["programs"], "pychopper_program": facts["pychopper_program"],
                  "pychopper_names": int(m.pychopper_names), "names_checked": int(m.names_checked), "tags": sorted(m.tags),
                  "eof_ok": pl["eof_ok"], "reference": alnpass.reference_summary(pl.get("reference")),
                  "index": {k: idx.get(k) for k in ("kind", "records", "mapped", "unmapped_placed", "no_coordinate", "stale")} if idx else None}
    return p


def _summary_pass(inp: LinkedInput, *, deadline, stall_timeout_s, inprocess) -> dict:
    from flowcellar.readdata import tables
    t0 = time.time()
    stop_wall = None if deadline is None else time.time() + max(deadline - time.monotonic(), 0.0)
    p = _new(inp)
    try:
        out = run_tasks(tables.summary_task, [(inp.candidate.path, stop_wall)], workers=1, stall_timeout_s=stall_timeout_s, inprocess=inprocess)
    except Exception as e:      # noqa: BLE001  a table this reader cannot use is an input error, not a crash of the scan
        p["state"] = "error"
        p["errors"] = [{"code": "unreadable", "message": str(e)[:300]}]
        return p
    if out["stalled"] or not out["results"]:
        p["state"] = "error"
        p["errors"] = [{"code": "stalled", "message": f"no progress for {stall_timeout_s:.0f} s; the pass was abandoned"}]
        return p
    res = out["results"][0]
    runs: dict[str, int] = {}
    labels: dict[str, int] = {}
    unrouted = 0
    for (key, label, flt), state in sorted(res["groups"].items(), key=lambda kv: tuple(str(x) for x in kv[0])):
        acc = ReadAcc.from_state(state)
        cells = acc.cells(aligned=False)
        cells["len"]["scalars"]["split_parents"] = 0
        _fragments(p, key, label, cells, flt)
        if key is None:
            unrouted += acc.reads
        else:
            runs[key] = runs.get(key, 0) + acc.reads
        labels[label or ""] = labels.get(label or "", 0) + acc.reads
    p["routing"] = {"runs": runs, "unrouted": unrouted, "default_run": None}
    total = res["rows"] + res.get("rows_skipped", 0)
    if res["error"]:
        p["state"] = "truncated"
        p["errors"] = [{"code": "truncated", "message": res["error"]}]
        p["coverage"] = {"kind": "truncated", "records": res["rows"], "records_total": None, "fraction": None}
    elif not res["finished"]:
        p["state"] = "partially_read"
        p["coverage"] = {"kind": "truncated", "records": res["rows"], "records_total": None, "fraction": None}
    else:
        p["coverage"] = {"kind": "full", "records": res["rows"], "records_total": total, "fraction": 1.0}
        if total and (res["rows_bad"] + res.get("rows_skipped", 0)) / total > 0.01:
            p["state"] = "error"
            p["errors"] = [{"code": "bad_rows", "message": f"{res['rows_bad'] + res.get('rows_skipped', 0):,} of {total:,} rows could not be read"}]
    seconds = time.time() - t0
    p["timings"] = {"seconds": round(seconds, 3), "mb_per_s": round(os.path.getsize(inp.candidate.path) / 1e6 / max(seconds, 1e-6), 1), "workers": 1}
    p["facts"] = {"kind": "summary", "records": res["rows"], "primary": res["rows"] - res["rows_bad"], "rows_bad": res["rows_bad"],
                  "rows_skipped": res.get("rows_skipped", 0), "labels": labels, "columns": res["columns"], "dorado": {}, "aligner": {},
                  "pychopper_program": False, "pychopper_names": 0, "names_checked": 0, "models": []}
    return p
