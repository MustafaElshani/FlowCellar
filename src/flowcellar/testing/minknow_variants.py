"""More fake MinKNOW run folders: the older file layouts, and the folders that are not an ordinary finished run.
Everything is written from scratch with invented names (``flowcellar.testing.fixtures``).

    from flowcellar.testing import minknow_variants as V
    V.make_old_layout(tmp, generation="5.3")      # 15-column summary, run_id = acquisition id, instrument = device
    V.make_without_scan_table(tmp)                # a report without bream_info: pores come from the scan messages
    V.make_running(tmp); V.make_interrupted(tmp); V.make_raw_only(tmp); V.make_lookalike(tmp)
    V.make_reused_flow_cell(tmp)                  # a first loading, two restarts, a second loading
    V.truncate_summary(run)                       # cut the sequencing summary in the middle of a line

What differs between MinKNOW generations, as far as a reader can tell (learned from real run folders):

===========  =====================================================================================================
generation   differences from the newest layout (6.10, written by ``fixtures.make_run``)
===========  =====================================================================================================
6.2, 5.9     sequencing summary has 15 columns (no ``filename_bam``) and ``run_id`` holds the ACQUISITION id;
5.8, 5.7     final summary: ``instrument`` = device position, ``position`` empty; the report has no
             ``user_info.kit_info`` and no ``sample rate`` tag
5.5          as above, and ``software_versions`` is a top-level member of the report
5.3          as 5.5, and: no ``stop_origin``; no ``target_run_until_criteria`` (the planned time is the argument
             ``--experiment_time``); temperature samples carry no ``target_temperature``; no ``speed`` tag;
             the final summary counts fast5 and fastq files only
===========  =====================================================================================================
"""
from __future__ import annotations

import json
import os
import shutil
import time
from dataclasses import asdict, replace

from flowcellar.testing.fixtures import FixtureRun, RunSpec, make_run

GENERATIONS = {"5.3": "5.3.7", "5.5": "5.5.5", "5.7": "5.7.5", "5.8": "5.8.6", "5.9": "5.9.12", "6.2": "6.2.8"}


def _spec(spec: RunSpec | None, **overrides) -> RunSpec:
    spec = spec or RunSpec()
    return replace(spec, **overrides) if overrides else spec


def _read_report(run: FixtureRun) -> dict:
    with open(run.files["minknow.report_json"]) as fh:
        return json.load(fh)


def _write_report(run: FixtureRun, rep: dict) -> None:
    with open(run.files["minknow.report_json"], "w", newline="") as fh:
        json.dump(rep, fh, indent=1)


def _sequencing(rep: dict) -> dict:
    return next(a for a in rep["acquisitions"] if a["acquisition_run_info"]["config_summary"].get("purpose") == "SEQUENCING")


def settle_stop_sample(run: FixtureRun) -> FixtureRun:
    """Make the last temperature sample read like the others unless it falls on the end of the run.

    MinKNOW takes a sample as it stops a run that reaches its planned end, and that sample reads low; a run that is
    stopped between two samples has no such sample. The base generator writes a low last sample for every run.
    """
    rep = _read_report(run)
    seq = _sequencing(rep)
    bucket_s = float(seq.get("bucket_size") or 0) * 60
    run_s = run.expected["minutes"] * 60
    for blk in seq.get("temperature", []):
        samples = blk["temperatures"]
        if len(samples) >= 2 and bucket_s and abs((len(samples) - 1) * bucket_s - run_s) < 1:
            continue                                  # taken as the run stops: stays low
        if not samples:
            continue
        ref = samples[0] if len(samples) > 1 else None
        for k in ("promethion", "minion"):
            if k in samples[-1]:
                for name in ("flowcell_temperature", "heatsink_temperature"):
                    if name in samples[-1][k]:
                        tgt = samples[-1].get("target_temperature", {}).get("minimum")
                        samples[-1][k][name] = ref[k][name] if ref else (tgt if tgt is not None else samples[-1][k][name] + 2.1)
    _write_report(run, rep)
    return run


def make(out_dir, spec: RunSpec | None = None, **overrides) -> FixtureRun:
    """``fixtures.make_run`` with a temperature track as a real run of that length has it."""
    return settle_stop_sample(make_run(out_dir, _spec(spec, **overrides)))


# ------------------------------------------------------------------------------------------- older layouts
def make_old_layout(out_dir, spec: RunSpec | None = None, *, generation: str = "5.9", **overrides) -> FixtureRun:
    """A finished run as an older MinKNOW wrote it. ``generation``: 5.3, 5.5, 5.7, 5.8, 5.9 or 6.2."""
    if generation not in GENERATIONS:
        raise ValueError(f"generation must be one of {sorted(GENERATIONS)}")
    spec = _spec(spec, minknow=GENERATIONS[generation], **overrides)
    run = make(out_dir, spec)
    gen = tuple(int(x) for x in generation.split("."))
    # ---- sequencing summary: 15 columns, run_id = acquisition id
    path = run.files["minknow.sequencing_summary"]
    with open(path) as fh:
        lines = fh.read().split("\n")
    header = lines[0].split("\t")
    drop, rid = header.index("filename_bam"), header.index("run_id")
    out = []
    for i, line in enumerate(lines):
        if not line:
            out.append(line)
            continue
        f = line.split("\t")
        if i:
            f[rid] = run.acquisition_id
        del f[drop]
        out.append("\t".join(f))
    with open(path, "w", newline="") as fh:
        fh.write("\n".join(out))
    # ---- final summary: instrument = device, position empty
    fs_path = run.files["minknow.final_summary"]
    kv = dict(line.rstrip("\n").split("=", 1) for line in open(fs_path) if "=" in line)
    kv["instrument"], kv["position"] = kv["position"], ""
    if gen < (5, 9):
        for k in [k for k in kv if k.startswith(("bam_", "bai_"))]:
            del kv[k]
    if gen < (5, 7):
        for k in [k for k in kv if k.startswith("pod5_")]:
            del kv[k]
    if gen >= (5, 9) and gen < (6, 2):
        for k in [k for k in kv if k.startswith("bai_")]:
            kv.pop(k, None)
    with open(fs_path, "w", newline="") as fh:
        fh.write("".join(f"{k}={v}\n" for k, v in kv.items()))
    # ---- report
    rep = _read_report(run)
    pri = rep["protocol_run_info"]
    pri["user_info"].pop("kit_info", None)
    pri["meta_info"]["tags"].pop("sample rate", None)
    pri.pop("phase_history", None)
    seq = _sequencing(rep)
    if gen < (5, 7):
        sw = pri.pop("software_versions")
        sw["guppy_connected_version"] = sw.pop("basecaller_connected_version", None)
        rep = {"software_versions": sw, **rep}
        rep["host"] = {"serial": pri["device"]["device_id"], "product_name": "fake sequencer"}
    else:
        rep["host"] = {}
    if gen < (5, 5):
        pri.pop("stop_origin", None)
        pri["meta_info"]["tags"].pop("speed", None)
        planned = float(seq["acquisition_run_info"].pop("target_run_until_criteria")["stop_criteria"]["criteria"]["runtime"]["value"])
        pri["args"].append(f"--experiment_time={planned / 3600:g}")
        for blk in seq.get("temperature", []):
            for s in blk["temperatures"]:
                s.pop("target_temperature", None)
        seq["acquisition_run_info"]["writer_summary"]["bytes_to_write_failed"] = "0"
    elif pri.get("state") != "PROTOCOL_COMPLETED":
        pri["stop_origin"] = {"user": {"identity": "['']"}}
    _write_report(run, rep)
    run.expected.update(summary_columns=15, minknow=GENERATIONS[generation], device_position=kv["instrument"])
    return run


def make_without_scan_table(out_dir, spec: RunSpec | None = None, **overrides) -> FixtureRun:
    """A report without ``bream_info``: pores per scan are only in the ``mux_scan_result`` messages."""
    run = make(out_dir, _spec(spec, **overrides))
    rep = _read_report(run)
    _sequencing(rep)["acquisition_run_info"].pop("bream_info", None)
    _write_report(run, rep)
    run.expected["scan_table_in_report"] = False
    return run


def make_with_message(out_dir, spec: RunSpec | None = None, *, severity: str = "ERROR", identifier: str = "device_error",
                      text: str = "The flow cell lost its connection", **overrides) -> FixtureRun:
    """A run in which MinKNOW showed one message of the given severity (ERROR or WARNING)."""
    run = make(out_dir, _spec(spec, **overrides))
    rep = _read_report(run)
    first = rep["user_messages"][0]
    rep["user_messages"].append({"time": first["time"], "severity": f"MESSAGE_SEVERITY_{severity}", "identifier": identifier,
                                 "user_message": text})
    _write_report(run, rep)
    run.expected["messages_" + severity.lower()] = 1
    return run


def set_stop(run: FixtureRun, *, state: str = "PROTOCOL_FINISHED_WITH_ERROR", stop_reason: str = "STOPPED_DEVICE_ERROR") -> FixtureRun:
    """Make the run end the way a failing instrument ends it."""
    rep = _read_report(run)
    rep["protocol_run_info"]["state"] = state
    _sequencing(rep)["acquisition_run_info"]["stop_reason"] = stop_reason
    _write_report(run, rep)
    return run


def set_files_in_fallback(run: FixtureRun, n: int = 3) -> FixtureRun:
    path = run.files["minknow.final_summary"]
    text = open(path).read().replace("pod5_files_in_fallback=0", f"pod5_files_in_fallback={n}")
    with open(path, "w", newline="") as fh:
        fh.write(text)
    return run


def set_temperature(run: FixtureRun, offset_c: float, *, sample: int = 1) -> FixtureRun:
    """Move one temperature sample (not the one taken as the run stops) away from the target."""
    rep = _read_report(run)
    s = _sequencing(rep)["temperature"][0]["temperatures"][sample]
    for k in ("promethion", "minion"):
        if k in s:
            for name in ("flowcell_temperature", "heatsink_temperature"):
                if name in s[k]:
                    s[k][name] += offset_c
    _write_report(run, rep)
    return run


# --------------------------------------------------------------------------- folders that are not finished
def _touch_all(run_dir: str, when: float) -> None:
    for base, dirs, files in os.walk(run_dir):
        for n in files + dirs:
            os.utime(os.path.join(base, n), (when, when))
    os.utime(run_dir, (when, when))


def make_running(out_dir, spec: RunSpec | None = None, **overrides) -> FixtureRun:
    """A run in progress: sample sheet, pore scans so far and raw data; the report, throughput, pore activity,
    sequencing summary and final summary do not exist yet. Every file was changed just now."""
    run = make(out_dir, _spec(spec, **overrides))
    for role in ("minknow.final_summary", "minknow.report_json", "minknow.throughput", "minknow.pore_activity",
                 "minknow.sequencing_summary", "minknow.barcode_alignment"):
        os.unlink(run.files.pop(role))
    _touch_all(run.run_dir, time.time())
    run.expected["run_state"] = "running"
    return run


def make_interrupted(out_dir, spec: RunSpec | None = None, *, with_reports: bool = True, age_hours: float = 48.0, **overrides) -> FixtureRun:
    """A run without a final summary whose files have not changed for two days. ``with_reports=False``: MinKNOW
    never wrote its end-of-run files either."""
    run = make(out_dir, _spec(spec, **overrides))
    gone = ["minknow.final_summary"]
    if not with_reports:
        gone += ["minknow.report_json", "minknow.throughput", "minknow.pore_activity", "minknow.sequencing_summary", "minknow.barcode_alignment"]
    for role in gone:
        os.unlink(run.files.pop(role))
    _touch_all(run.run_dir, time.time() - age_hours * 3600)
    run.expected["run_state"] = "interrupted"
    return run


def make_raw_only(out_dir, spec: RunSpec | None = None, **overrides) -> FixtureRun:
    """A folder that holds raw data and ``other_reports`` only: every identity file is gone."""
    run = make(out_dir, _spec(spec, **overrides))
    for role in list(run.files):
        if role != "minknow.pore_scan":
            os.unlink(run.files.pop(role))
    run.expected["kind"] = "run_raw_only"
    return run


def make_lookalike(out_dir, spec: RunSpec | None = None, **overrides) -> str:
    """A folder named like a run folder that holds only basecaller output (``bam_pass``). Returns its path."""
    run = make(out_dir, _spec(spec, **overrides))
    look = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(run.run_dir))), "basecalled",
                        os.path.basename(run.run_dir))
    os.makedirs(os.path.join(look, "bam_pass"), exist_ok=True)
    with open(os.path.join(look, "bam_pass", f"{run.run_key}_0.bam"), "wb") as fh:
        fh.write(b"")
    shutil.rmtree(os.path.dirname(os.path.dirname(run.run_dir)))
    return look


def truncate_summary(run: FixtureRun, fraction: float = 0.6, *, role: str = "minknow.sequencing_summary") -> int:
    """Cut a table at ``fraction`` of its size, in the middle of a line. Returns the number of whole data lines left."""
    path = run.files[role]
    data = open(path, "rb").read()
    cut = int(len(data) * fraction)
    while cut < len(data) and data[cut - 1:cut] == b"\n":
        cut += 3
    with open(path, "wb") as fh:
        fh.write(data[:cut])
    return data[:cut].count(b"\n") - 1


# ------------------------------------------------------------------------------------------ sets of runs
def make_short(out_dir, *, minutes: int = 5, spec: RunSpec | None = None, **overrides) -> FixtureRun:
    """A run stopped after a few minutes (under 1 h: too short to judge) or after a few hours."""
    base = dict(hours=minutes / 60.0, planned_hours=72.0, state="PROTOCOL_STOPPED_BY_USER", reads=max(50, minutes * 8))
    base.update(overrides)
    return make(out_dir, _spec(spec, **base))


def make_promethion_72h(out_dir, spec: RunSpec | None = None, **overrides) -> FixtureRun:
    """A PromethION-like run of a little over 72 h (3,000 channels) that is small enough for tests: few reads, a
    pore scan every 8 h."""
    base = dict(family="promethion", hours=72.2, planned_hours=72.0, reads=6000, scan_period_h=8.0, flow_cell_id="TST00072")
    base.update(overrides)
    return make(out_dir, _spec(spec, **base))


def make_reused_flow_cell(out_dir, *, flow_cell_id: str = "TST00009", seed: int = 40) -> list[FixtureRun]:
    """One flow cell, four starts: a first loading (26 h), then on a later day a start of 5 minutes and one of
    1.5 h (restarts), and the second loading (26 h) with fewer pores."""
    common = dict(flow_cell_id=flow_cell_id, experiment="demo_reuse")
    a = make(out_dir, RunSpec(sample="load_1", start="2026-04-06T09:00:00Z", hours=26.1, reads=5000, seed=seed, **common))
    b = make(out_dir, RunSpec(sample="load_2_try_1", start="2026-04-08T09:00:00Z", hours=5 / 60, planned_hours=72, reads=40, seed=seed + 1,
                                  first_scan_pore_frac=0.40, state="PROTOCOL_STOPPED_BY_USER", **common))
    c = make(out_dir, RunSpec(sample="load_2_try_2", start="2026-04-08T09:20:00Z", hours=1.5, planned_hours=72, reads=300, seed=seed + 2,
                                  first_scan_pore_frac=0.36, state="PROTOCOL_STOPPED_BY_USER", **common))
    d = make(out_dir, RunSpec(sample="load_2", start="2026-04-08T11:10:00Z", hours=26.1, reads=3000, seed=seed + 3,
                                  first_scan_pore_frac=0.33, **common))
    return [a, b, c, d]


def make_series(out_dir, n: int = 9, *, hours: float = 24.2, reads: int = 2400, first_day: str = "2026-01-05", seed: int = 100,
                **overrides) -> list[FixtureRun]:
    """``n`` runs of one comparison group on ``n`` different days, each on its own flow cell and long enough to
    teach the lab-typical values: a history for the checks that are placed against earlier runs."""
    import datetime as dt
    d0 = dt.date.fromisoformat(first_day)
    out = []
    for i in range(n):
        day = d0 + dt.timedelta(days=3 * i)
        spec = RunSpec(flow_cell_id=f"TST{10100 + i:05d}", experiment="demo_series", sample=f"pool_{i + 1}",
                       start=f"{day.isoformat()}T10:00:00Z", hours=hours, reads=reads + 40 * (i % 5), seed=seed + i, scan_period_h=6.0)
        out.append(make(out_dir, replace(spec, **overrides) if overrides else spec))
    return out


def spec_dict(run: FixtureRun) -> dict:
    return asdict(run.spec)


KINDS = ("old-5.3", "old-5.5", "old-5.9", "old-6.2", "no-scan-table", "running", "interrupted", "interrupted-no-reports", "raw-only",
         "lookalike", "short-5min", "short-90min", "truncated-summary", "error-message", "promethion-72h", "reused-flow-cell", "series")


def make_kind(out_dir, kind: str) -> list[str]:
    """Write the fake folder(s) of one kind (see KINDS) and return their paths."""
    if kind.startswith("old-"):
        return [make_old_layout(out_dir, generation=kind[4:]).run_dir]
    if kind == "no-scan-table":
        return [make_without_scan_table(out_dir).run_dir]
    if kind == "running":
        return [make_running(out_dir).run_dir]
    if kind == "interrupted":
        return [make_interrupted(out_dir).run_dir]
    if kind == "interrupted-no-reports":
        return [make_interrupted(out_dir, with_reports=False).run_dir]
    if kind == "raw-only":
        return [make_raw_only(out_dir).run_dir]
    if kind == "lookalike":
        return [make_lookalike(out_dir)]
    if kind == "short-5min":
        return [make_short(out_dir, minutes=5).run_dir]
    if kind == "short-90min":
        return [make_short(out_dir, minutes=90).run_dir]
    if kind == "truncated-summary":
        run = make(out_dir)
        truncate_summary(run)
        return [run.run_dir]
    if kind == "error-message":
        return [make_with_message(out_dir).run_dir]
    if kind == "promethion-72h":
        return [make_promethion_72h(out_dir).run_dir]
    if kind == "reused-flow-cell":
        return [r.run_dir for r in make_reused_flow_cell(out_dir)]
    if kind == "series":
        return [r.run_dir for r in make_series(out_dir)]
    raise ValueError(f"kind must be one of {', '.join(KINDS)}")


def main(argv: list[str] | None = None) -> int:
    import argparse
    p = argparse.ArgumentParser(prog="python -m flowcellar.testing.minknow_variants",
                                description="Write fake MinKNOW run folders of the less ordinary kinds (invented names only).")
    p.add_argument("out_dir")
    p.add_argument("--kind", choices=KINDS, action="append", help="which folder to write (repeatable); default: every kind")
    a = p.parse_args(argv)
    for kind in a.kind or KINDS:
        for path in make_kind(os.path.join(a.out_dir, kind), kind):
            print(f"{kind}\t{path}")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
