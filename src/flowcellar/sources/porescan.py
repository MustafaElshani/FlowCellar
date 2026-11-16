"""Streaming pass over ``other_reports/pore_scan_data_*.csv`` (one row per channel, well and pore scan).

Reads eight of its 50 to 64 columns by name, so every MinKNOW generation is read by the same code. Gives the well
map, open-pore current and noise per scan, the running voltage per scan and the scan results (a cross-check of
the report, and the only source of pores per scan when the report has none).
"""
from __future__ import annotations

import os
import time

import numpy as np

from flowcellar.aggregates import build
from flowcellar.model import partial as P
from flowcellar.model.axes import bin_index, n_cells
from flowcellar.modules.common import GOOD_SCAN, SCAN_CATS, VERSIONS
from flowcellar.sources import workers
from flowcellar.util import fingerprint as _fp
from flowcellar.util.timeutil import from_epoch

ROLE = "minknow.pore_scan"
MODULES = {"fc": VERSIONS["fc"]}
REQUIRED = ("channel", "scan_number", "mux_scan_assessment")
WANT = ("channel", "well", "scan_number", "mux_scan_assessment", "pore_median", "pore_sd", "scan_running_voltage_mv",
        "seconds_since_start_of_run")
MAX_SCANS = 1024
PIECE = 16 << 20


def probe_task(args: dict, progress) -> dict:
    path = args["path"]
    st = os.stat(path)
    out = {"size": st.st_size, "mtime": st.st_mtime, "mtime_ns": st.st_mtime_ns, "fingerprint": _fp.fingerprint(path)}
    with open(path, "rb") as fh:
        line = fh.readline(1 << 20)
        out["header"] = [h.strip() for h in line.decode("utf-8", errors="replace").rstrip("\r\n").split(",")]
        out["body_start"] = fh.tell()
    return out


def _parse(data: bytes, header: list[str]):
    import pyarrow as pa
    import pyarrow.csv as pc
    types = {"channel": pa.int32(), "well": pa.int32(), "scan_number": pa.int32(),
             "mux_scan_assessment": pa.dictionary(pa.int32(), pa.string()), "pore_median": pa.float64(), "pore_sd": pa.float64(),
             "scan_running_voltage_mv": pa.float64(), "seconds_since_start_of_run": pa.float64()}
    cols = [c for c in WANT if c in header]
    skipped = [0]

    def on_invalid(row):
        skipped[0] += 1
        return "skip"
    tbl = pc.read_csv(pa.BufferReader(data),
                      read_options=pc.ReadOptions(column_names=header, use_threads=False, block_size=len(data) + 16),
                      parse_options=pc.ParseOptions(invalid_row_handler=on_invalid),
                      convert_options=pc.ConvertOptions(include_columns=cols, column_types={c: types[c] for c in cols},
                                                        null_values=["", "nan", "NaN", "None", "inf", "-inf"], strings_can_be_null=True))
    return tbl, cols, skipped[0]


def scan_task(args: dict, progress) -> dict:
    """Runs in a child: the whole table -> per-scan arrays."""
    path, channels, header = args["path"], args["channels"], args["header"]
    deadline = time.monotonic() + args["remaining_s"] if args.get("remaining_s") is not None else None
    C = channels + 1
    nc, nn = n_cells("current_2pa"), n_cells("noise_0p1pa")
    good = np.zeros((MAX_SCANS, C), np.int16)
    h_cur = np.zeros((MAX_SCANS, nc), np.int64)
    h_noise = np.zeros((MAX_SCANS, nn), np.int64)
    cats: dict[str, np.ndarray] = {}
    volt_sum, volt_n = np.zeros(MAX_SCANS), np.zeros(MAX_SCANS, np.int64)
    t_min = np.full(MAX_SCANS, np.inf)
    rows = bad = max_scan = 0
    pieces = workers.Pieces(path, args["body_start"], args["size"], args["body_start"], piece=PIECE, deadline=deadline)
    for data in pieces:
        try:
            batch, cols, skipped = _parse(data, header)
        except Exception:  # noqa: BLE001  a block that cannot be parsed is counted and left out
            bad += data.count(b"\n")
            progress(pieces.done)
            continue
        bad += skipped
        n = batch.num_rows
        if not n:
            continue
        ch_a, sc_a, ass = batch["channel"].combine_chunks(), batch["scan_number"].combine_chunks(), batch["mux_scan_assessment"].combine_chunks()
        valid = np.ones(n, bool)
        if ch_a.null_count or sc_a.null_count:
            valid &= ~(ch_a.is_null().to_numpy(zero_copy_only=False) | sc_a.is_null().to_numpy(zero_copy_only=False))
            bad += int((~valid).sum())
        ch = np.clip(ch_a.fill_null(0).to_numpy(zero_copy_only=False).astype(np.int64), 0, channels)
        sc = np.clip(sc_a.fill_null(0).to_numpy(zero_copy_only=False).astype(np.int64), 0, MAX_SCANS - 1)
        rows += int(valid.sum())
        if valid.any():
            max_scan = max(max_scan, int(sc[valid].max()))
        names = ass.dictionary.to_pylist()
        code = ass.indices.fill_null(-1).to_numpy(zero_copy_only=False).astype(np.int64)
        is_good = np.zeros(n, bool)
        for k, name in enumerate(names):
            m = (code == k) & valid
            if name is None or not m.any():
                continue
            cats.setdefault(name, np.zeros(MAX_SCANS, np.int64))
            cats[name] += np.bincount(sc[m], minlength=MAX_SCANS)
            if name in GOOD_SCAN:
                is_good |= m
        np.add.at(good, (sc[is_good], ch[is_good]), 1)
        if "pore_median" in cols:
            pm = batch["pore_median"].combine_chunks().fill_null(np.nan).to_numpy(zero_copy_only=False).astype(np.float64)
            m = is_good & np.isfinite(pm)
            h_cur += np.bincount(sc[m] * nc + bin_index(pm[m], "current_2pa"), minlength=MAX_SCANS * nc).reshape(MAX_SCANS, nc)
        if "pore_sd" in cols:
            sd = batch["pore_sd"].combine_chunks().fill_null(np.nan).to_numpy(zero_copy_only=False).astype(np.float64)
            m = is_good & np.isfinite(sd)
            h_noise += np.bincount(sc[m] * nn + bin_index(sd[m], "noise_0p1pa"), minlength=MAX_SCANS * nn).reshape(MAX_SCANS, nn)
        if "scan_running_voltage_mv" in cols:
            v = batch["scan_running_voltage_mv"].combine_chunks().fill_null(np.nan).to_numpy(zero_copy_only=False).astype(np.float64)
            m = valid & np.isfinite(v)
            volt_sum += np.bincount(sc[m], weights=v[m], minlength=MAX_SCANS)
            volt_n += np.bincount(sc[m], minlength=MAX_SCANS)
        if "seconds_since_start_of_run" in cols:
            s = batch["seconds_since_start_of_run"].combine_chunks().fill_null(np.nan).to_numpy(zero_copy_only=False).astype(np.float64)
            m = valid & np.isfinite(s)
            np.minimum.at(t_min, sc[m], s[m])
        progress(pieces.done)
    # scan numbers are 1-based in the file; index 0 of every output is scan 1
    idx = np.arange(1, max_scan + 1, dtype=np.int64)
    g_ = good[idx][:, 1:] if max_scan else np.zeros((0, channels), np.int16)
    last_good = np.zeros(channels, np.int64)
    for k in range(len(idx)):
        last_good[g_[k] > 0] = k + 1
    other = sum((v for name, v in cats.items() if name not in SCAN_CATS), np.zeros(MAX_SCANS, np.int64))
    by_cat = {c: [int(x) for x in (cats.get(c, np.zeros(MAX_SCANS, np.int64)) + (other if c == "other" else 0))[idx]] for c in SCAN_CATS}
    return {"rows": rows, "bad_rows": bad, "n_scans": int(max_scan), "counts": by_cat, "raw_categories": sorted(cats),
            "voltage_mv": [round(float(volt_sum[i] / volt_n[i]), 1) if volt_n[i] else None for i in idx],
            "t_s": [int(t_min[i]) if np.isfinite(t_min[i]) else None for i in idx],
            "h_current": h_cur[idx], "h_noise": h_noise[idx],
            "good_first": [int(x) for x in g_[0]] if max_scan else [], "good_last": [int(x) for x in g_[-1]] if max_scan else [],
            "last_good_scan": [int(x) for x in last_good], "done": pieces.done, "stopped": pieces.stopped,
            "cut_tail": pieces.cut_tail, "cut_bytes": pieces.cut_bytes}


def run_pass(path: str, probe: dict, *, channels: int, deadline: float | None = None, stall_timeout_s: float = 120.0,
             on_progress=None, gentle: bool = False) -> dict:
    t0 = time.monotonic()
    meta: dict = {"state": "error", "errors": [], "seconds": 0.0, "coverage_frac": 0.0, "rows": 0, "bytes_read": 0}
    missing = [c for c in REQUIRED if c not in probe["header"]]
    if missing:
        meta["errors"].append({"code": "columns_missing", "message": "The pore-scan table lacks the columns " + ", ".join(missing) + "."})
        return {"scan": None, "meta": meta}
    remaining = (deadline - time.monotonic()) if deadline is not None else None
    r = workers.run_tasks(f"{__name__}:scan_task", [{"path": path, "channels": channels, "header": probe["header"],
                                                     "body_start": probe["body_start"], "size": probe["size"], "remaining_s": remaining,
                                                     "_nice": 10 if gentle else 0}],
                          stall_timeout_s=stall_timeout_s, on_progress=on_progress)[0]
    meta["seconds"] = time.monotonic() - t0
    if r.state != "ok":
        code = {"stalled": "stalled", "died": "worker_died"}.get(r.state, "read_error")
        meta["errors"].append({"code": code, "message": (r.error or r.state).splitlines()[0][:300]})
        return {"scan": None, "meta": meta}
    v = r.value
    body = max(probe["size"] - probe["body_start"], 0)
    meta.update(rows=v["rows"], bytes_read=v["done"], coverage_frac=v["done"] / body if body else 1.0, bad_rows=v["bad_rows"])
    if v["bad_rows"] > 0.01 * max(v["rows"] + v["bad_rows"], 1):
        meta["errors"].append({"code": "too_many_bad_rows", "message": f"{v['bad_rows']:,} of {v['rows'] + v['bad_rows']:,} rows could not be read."})
    elif v["stopped"]:
        meta["state"] = "partially_read"
    elif v["cut_tail"]:
        meta["state"] = "truncated"
    else:
        meta["state"] = "ok"
    return {"scan": v, "meta": meta}


def fragment(scan: dict, run_key: str) -> dict:
    """Per-scan arrays -> the ``fc`` fragment (already in the units of the bundle: wells, pA, mV)."""
    n = scan["n_scans"]
    cell: dict = {"at": [], "scalars": {"pore_scan_rows": int(scan["rows"])}, "series": {}, "hists2d": {}, "channels": {}}
    if n:
        t = [0 if x is None else x for x in scan["t_s"]]
        cell["series"]["fc_pore_scan_table_series"] = build.series("mux_scan_index", {**scan["counts"], "t_s": t}, unit="wells",
                                                                  levels={"t_s": "same"}, n_bins=n)
        if any(v is not None for v in scan["voltage_mv"]):
            cell["series"]["fc_bias_voltage_series"] = build.series("mux_scan_index", {"voltage_mv": scan["voltage_mv"], "t_s": scan["t_s"]},
                                                                   unit="mV", levels={"voltage_mv": "same", "t_s": "same"})
        hc, hn = np.asarray(scan["h_current"]), np.asarray(scan["h_noise"])
        cell["hists2d"]["fc_open_pore_current"] = build.hist2d("mux_scan_index", "current_2pa", hc[:, 1:-1],
                                                               outside=int(hc[:, 0].sum() + hc[:, -1].sum()))
        cell["hists2d"]["fc_open_pore_noise"] = build.hist2d("mux_scan_index", "noise_0p1pa", hn[:, 1:-1],
                                                             outside=int(hn[:, 0].sum() + hn[:, -1].sum()))
        cell["channels"]["fc_well_map"] = {"layout": "channel",
                                           "values": {"good_wells_first_scan": scan["good_first"], "good_wells_last_scan": scan["good_last"],
                                                      "last_good_scan": scan["last_good_scan"]},
                                           "merge": {"good_wells_first_scan": "first", "good_wells_last_scan": "last", "last_good_scan": "max"}}
    return P.fragment(run_key, "fc", dims=[], cells=[cell], facts={"raw_categories": scan["raw_categories"]})


def make_partial(path: str, probe: dict, result: dict, run_key: str, *, channels: int) -> dict:
    scan, meta = result["scan"], result["meta"]
    p = P.new_partial(role=ROLE, input_path=path, fp=probe["fingerprint"], size=probe["size"], mtime=from_epoch(probe["mtime"]),
                      modules=MODULES)
    p["state"] = meta["state"]
    rows = meta.get("rows", 0)
    if meta["state"] == "ok":
        p["coverage"] = {"kind": "full", "records": rows, "records_total": rows, "fraction": 1.0}
    else:
        p["coverage"] = {"kind": "truncated" if meta["state"] == "truncated" else "byte_ranges", "records": rows, "records_total": None,
                         "fraction": round(meta.get("coverage_frac", 0.0), 6)}
    p["routing"] = {"runs": {run_key: rows}, "unrouted": 0, "default_run": run_key}
    # a table that was not read to its end lacks its last scans: nothing of it is used
    p["fragments"] = [fragment(scan, run_key)] if scan and meta["state"] in ("ok", "truncated") else []
    p["facts"] = {"rows": rows, "channel_count": channels, "columns": len(probe["header"]), "bad_rows": meta.get("bad_rows", 0),
                  "mtime_ns": probe.get("mtime_ns")}
    p["errors"] = list(meta["errors"])
    secs = meta["seconds"]
    p["timings"] = {"seconds": round(secs, 3), "mb_per_s": round(meta.get("bytes_read", 0) / 1e6 / secs, 1) if secs > 0 else None, "workers": 1}
    return p
