"""One streaming pass over a MinKNOW sequencing summary (tier 1, full detail).

The file is cut into one byte range per worker process; a worker owns every line that starts in its range, reads
it in pieces, parses the columns it needs with pyarrow and adds to integer numpy accumulators, which the parent
merges by element-wise addition (maximum for the time of the last read). Results are therefore identical for any
number of workers. The reader is header driven: every MinKNOW generation is read by the same code (15 or 16
columns, with or without ``filename_bam``; ``run_id`` is never used). The file is never held in memory: one piece
per worker at a time.

With a deadline every worker stops inside its own range, so what was read is spread evenly over the file; the
result then says which fraction it covers.
"""
from __future__ import annotations

import os
import time

import numpy as np

from flowcellar.aggregates import build
from flowcellar.model import partial as P
from flowcellar.model.axes import bin_index, n_cells
from flowcellar.modules.common import END_REASONS, TIME_AXIS_DAYS, VERSIONS, grid_bin_s
from flowcellar.sources import workers
from flowcellar.util import fingerprint as _fp
from flowcellar.util.timeutil import from_epoch

ROLE = "minknow.sequencing_summary"
MODULES = {k: VERSIONS[k] for k in ("dyn", "len.est", "fc")}
REQUIRED = ("channel", "minknow_events", "start_time", "duration")
OPTIONAL = ("mux", "end_reason", "pore_type")
HOURS = TIME_AXIS_DAYS * 24
N_ER = len(END_REASONS)
PIECE = 32 << 20
MIN_BYTES_PER_WORKER = 8 << 20
BAD_ROW_LIMIT = 0.01               # above this share of unreadable rows the input is an error
LEN_COUNT_CAP = 1 << 22            # the counting array of events per read ends here; longer reads share its last cell


# ------------------------------------------------------------------------------------------ accumulators
def new_acc(channels: int) -> dict:
    C = channels + 1
    G = (TIME_AXIS_DAYS * 86400) // grid_bin_s(channels)
    nd, nl, ns = n_cells("duration_log40"), n_cells("read_length_log40"), n_cells("speed_5bps")
    z = np.zeros
    return {
        "n": z(1, np.int64), "ev_total": z(1, np.int64), "bad_rows": z(1, np.int64), "late_reads": z(1, np.int64),
        "out_of_range_channel": z(1, np.int64),
        "h_dur": z((N_ER, nd), np.int64),
        "h_len_n": z((N_ER, nl), np.int64), "h_len_ev": z((N_ER, nl), np.int64),
        "h_spd": z((N_ER, ns), np.int64),                # reads of 1 s or longer
        "n_short_dur": z(N_ER, np.int64),                # reads under 1 s: no speed
        "h_spd_t": z((HOURS, ns), np.int64),             # reads that ended normally, 1 s or longer
        "ts_n": z((N_ER, HOURS), np.int64), "ts_ev": z((N_ER, HOURS), np.int64), "ts_dur_ms": z((N_ER, HOURS), np.int64),
        "ch_n": z(C, np.int64), "ch_ev": z(C, np.int64), "ch_dur_ms": z(C, np.int64), "max_ch_last_ms": z(C, np.int64),
        "chmux_n": z((C, 5), np.int64), "mux_n": z(5, np.int64), "mux_ev": z(5, np.int64),
        "grid_ms": z((C, G), np.int64),
        "max_events": z(1, np.int64), "max_dur_ms": z(1, np.int64),
        "ev_count": z(1, np.int64),                      # reads per number of events: grows to the longest read seen
        "pore_type": {},
    }


def _add_counts(acc: dict, c: np.ndarray) -> None:
    """Add a counting array of another length to ``acc["ev_count"]``, growing it when needed."""
    cur = acc["ev_count"]
    if len(c) > len(cur):
        grown = np.zeros(len(c), np.int64)
        grown[:len(cur)] = cur
        acc["ev_count"] = cur = grown
    cur[:len(c)] += c


def exact_lengths(counts) -> dict | None:
    """Exact N50 and median of the events per read, from the counting array (catalogue: ``len_est_n50``,
    ``len_est_median``). N50: the longest reads that together hold half of all events end at this length. Median:
    the middle read, or the mean of the two middle reads."""
    c = np.asarray(counts, np.int64)
    n = int(c.sum())
    if n <= 0:
        return None
    w = c * np.arange(len(c), dtype=np.int64)
    total = int(w.sum())
    n50 = None
    if total > 0:
        from_top = np.cumsum(w[::-1])
        n50 = len(c) - 1 - int(np.argmax(2 * from_top >= total))
    cum = np.cumsum(c)
    lo = int(np.searchsorted(cum, (n - 1) // 2 + 1, side="left"))
    hi = int(np.searchsorted(cum, n // 2 + 1, side="left"))
    return {"reads": n, "events_n50": n50, "events_median": (lo + hi) / 2.0,
            "over_cap": int(c[LEN_COUNT_CAP - 1]) if len(c) >= LEN_COUNT_CAP else 0}


def merge_acc(a: dict, b: dict) -> None:
    for k, v in b.items():
        if k == "pore_type":
            for name, n in v.items():
                a["pore_type"][name] = a["pore_type"].get(name, 0) + int(n)
        elif k == "ev_count":
            _add_counts(a, v)
        elif k.startswith("max_"):
            np.maximum(a[k], v, out=a[k])
        else:
            a[k] += v


def _h2(i, j, ni, nj, w=None):
    flat = i.astype(np.int64) * nj + j
    if w is None:
        return np.bincount(flat, minlength=ni * nj).reshape(ni, nj)
    return np.rint(np.bincount(flat, weights=w, minlength=ni * nj)).astype(np.int64).reshape(ni, nj)


def update(acc: dict, cols: dict, channels: int, ratio: float) -> None:
    """Add one batch. ``cols`` holds numpy arrays: channel, mux, minknow_events, start_time, duration (all of one
    length), ``end_reason`` as (codes, names) and optionally ``pore_type`` as (codes, names)."""
    ch = cols["channel"].astype(np.int64)
    n = len(ch)
    if n == 0:
        return
    mux = np.clip(cols["mux"].astype(np.int64), 0, 4)
    ev_i = np.maximum(cols["minknow_events"].astype(np.int64), 0)
    ev = ev_i.astype(np.float64)
    start = np.maximum(cols["start_time"].astype(np.float64), 0.0)
    dur = np.maximum(cols["duration"].astype(np.float64), 0.0)
    codes, names = cols["end_reason"]
    lut = np.array([END_REASONS.index(x) if x in END_REASONS else N_ER - 1 for x in names] or [N_ER - 1], np.int64)
    er = lut[codes]
    bad = (ch < 1) | (ch > channels)
    if bad.any():
        acc["out_of_range_channel"] += int(bad.sum())
        ch = np.where(bad, 0, ch)
    C = channels + 1
    acc["n"] += n
    acc["ev_total"] += int(ev_i.sum())
    acc["max_events"][0] = max(int(acc["max_events"][0]), int(ev_i.max()))
    _add_counts(acc, np.bincount(np.minimum(ev_i, LEN_COUNT_CAP - 1)))
    dms = np.rint(dur * 1000.0)
    acc["max_dur_ms"][0] = max(int(acc["max_dur_ms"][0]), int(dms.max()))
    acc["h_dur"] += _h2(er, bin_index(dur, "duration_log40"), N_ER, acc["h_dur"].shape[1])
    est = ev * ratio
    li = bin_index(est, "read_length_log40")
    nl = acc["h_len_n"].shape[1]
    acc["h_len_n"] += _h2(er, li, N_ER, nl)
    acc["h_len_ev"] += _h2(er, li, N_ER, nl, w=ev)
    ok = dur >= 1.0
    ns = acc["h_spd"].shape[1]
    acc["h_spd"] += _h2(er[ok], bin_index(est[ok] / dur[ok], "speed_5bps"), N_ER, ns)
    acc["n_short_dur"] += np.bincount(er[~ok], minlength=N_ER)
    hour_raw = (start / 3600.0).astype(np.int64)
    acc["late_reads"] += int((hour_raw >= HOURS).sum())
    hour = np.clip(hour_raw, 0, HOURS - 1)
    sp = ok & (er == 0)
    acc["h_spd_t"] += _h2(hour[sp], bin_index(est[sp] / dur[sp], "speed_5bps"), HOURS, ns)
    acc["ts_n"] += _h2(er, hour, N_ER, HOURS)
    acc["ts_ev"] += _h2(er, hour, N_ER, HOURS, w=ev)
    acc["ts_dur_ms"] += _h2(er, hour, N_ER, HOURS, w=dms)
    acc["ch_n"] += np.bincount(ch, minlength=C)
    acc["ch_ev"] += np.rint(np.bincount(ch, weights=ev, minlength=C)).astype(np.int64)
    acc["ch_dur_ms"] += np.rint(np.bincount(ch, weights=dms, minlength=C)).astype(np.int64)
    end_ms = (np.rint(start * 1000.0) + dms).astype(np.int64)
    np.maximum.at(acc["max_ch_last_ms"], ch, end_ms)
    acc["chmux_n"] += _h2(ch, mux, C, 5)
    acc["mux_n"] += np.bincount(mux, minlength=5)
    acc["mux_ev"] += np.rint(np.bincount(mux, weights=ev, minlength=5)).astype(np.int64)
    G = acc["grid_ms"].shape[1]
    gb = np.clip((start / grid_bin_s(channels)).astype(np.int64), 0, G - 1)
    acc["grid_ms"] += _h2(ch, gb, C, G, w=dms)
    if "pore_type" in cols:
        pcodes, pnames = cols["pore_type"]
        cnt = np.bincount(pcodes, minlength=len(pnames))
        for name, c in zip(pnames, cnt):
            if c:
                acc["pore_type"][str(name)] = acc["pore_type"].get(str(name), 0) + int(c)


# ------------------------------------------------------------------------------------------------ parsing
def _dict_col(col):
    import pyarrow as pa
    c = col.combine_chunks() if isinstance(col, pa.ChunkedArray) else col
    names = [("unknown" if x is None else x) for x in c.dictionary.to_pylist()]
    idx = c.indices
    codes = idx.fill_null(len(names)).to_numpy(zero_copy_only=False).astype(np.int64) if idx.null_count else \
        idx.to_numpy(zero_copy_only=False).astype(np.int64)
    if idx.null_count:
        names = names + ["unknown"]
    return codes, names


def parse_block(data: bytes, header: list[str]) -> tuple[dict, int]:
    """(columns of one block of whole lines, rows that could not be read)."""
    import pyarrow as pa
    import pyarrow.csv as pc
    want = [c for c in REQUIRED + OPTIONAL if c in header]
    types = {"channel": pa.int32(), "mux": pa.int32(), "minknow_events": pa.int64(), "start_time": pa.float64(),
             "duration": pa.float64(), "end_reason": pa.dictionary(pa.int32(), pa.string()),
             "pore_type": pa.dictionary(pa.int32(), pa.string())}
    skipped = [0]

    def on_invalid(row):
        skipped[0] += 1
        return "skip"
    if not data.isascii():
        # bytes that are not text (a damaged file): the fast reader would hand such a row to the handler above as
        # text and fail while doing so, printing a traceback; the line-by-line reader keeps what can be read
        return _parse_block_slow(data, header)
    try:
        tbl = pc.read_csv(pa.BufferReader(data),
                          read_options=pc.ReadOptions(column_names=header, use_threads=False, block_size=len(data) + 16),
                          parse_options=pc.ParseOptions(delimiter="\t", quote_char=False, invalid_row_handler=on_invalid),
                          convert_options=pc.ConvertOptions(include_columns=want, column_types={k: v for k, v in types.items() if k in want}))
        num = {k: tbl[k] for k in ("channel", "mux", "minknow_events", "start_time", "duration") if k in want}
        if any(c.null_count for c in num.values()):
            raise ValueError("empty numeric fields")
    except Exception:  # noqa: BLE001  a damaged block is read line by line instead
        return _parse_block_slow(data, header)
    n = tbl.num_rows
    cols: dict = {k: v.to_numpy() for k, v in num.items()}
    if "mux" not in cols:
        cols["mux"] = np.zeros(n, np.int64)
    cols["end_reason"] = _dict_col(tbl["end_reason"]) if "end_reason" in want else (np.zeros(n, np.int64), ["unknown"])
    if "pore_type" in want:
        cols["pore_type"] = _dict_col(tbl["pore_type"])
    return cols, skipped[0]


def _parse_block_slow(data: bytes, header: list[str]) -> tuple[dict, int]:
    """Line by line: keeps every row whose needed fields can be read, counts the others."""
    ix = {c: header.index(c) for c in REQUIRED + OPTIONAL if c in header}
    ncol = len(header)
    out = {k: [] for k in ("channel", "mux", "minknow_events", "start_time", "duration")}
    er_names: dict[str, int] = {}
    pt_names: dict[str, int] = {}
    er_codes: list[int] = []
    pt_codes: list[int] = []
    bad = 0
    for line in data.decode("utf-8", errors="replace").split("\n"):
        if not line:
            continue
        f = line.rstrip("\r").split("\t")
        if len(f) != ncol:
            bad += 1
            continue
        try:
            ch, ev = int(f[ix["channel"]]), int(f[ix["minknow_events"]])
            st, du = float(f[ix["start_time"]]), float(f[ix["duration"]])
            mx = int(f[ix["mux"]]) if "mux" in ix else 0
        except ValueError:
            bad += 1
            continue
        if not (np.isfinite(st) and np.isfinite(du)):
            bad += 1
            continue
        out["channel"].append(ch)
        out["mux"].append(mx)
        out["minknow_events"].append(ev)
        out["start_time"].append(st)
        out["duration"].append(du)
        er = f[ix["end_reason"]] if "end_reason" in ix else "unknown"
        er_codes.append(er_names.setdefault(er, len(er_names)))
        if "pore_type" in ix:
            pt = f[ix["pore_type"]]
            pt_codes.append(pt_names.setdefault(pt, len(pt_names)))
    cols: dict = {"channel": np.asarray(out["channel"], np.int64), "mux": np.asarray(out["mux"], np.int64),
                  "minknow_events": np.asarray(out["minknow_events"], np.int64),
                  "start_time": np.asarray(out["start_time"], np.float64), "duration": np.asarray(out["duration"], np.float64),
                  "end_reason": (np.asarray(er_codes, np.int64), list(er_names) or ["unknown"])}
    if "pore_type" in ix:
        cols["pore_type"] = (np.asarray(pt_codes, np.int64), list(pt_names))
    return cols, bad


# ------------------------------------------------------------------------------------------------- workers
def probe_task(args: dict, progress) -> dict:
    """Runs in a child: size, time, fingerprint, header and last byte of the summary. Reads about 2 MiB."""
    path = args["path"]
    st = os.stat(path)
    out = {"size": st.st_size, "mtime": st.st_mtime, "mtime_ns": st.st_mtime_ns, "fingerprint": _fp.fingerprint(path)}
    with open(path, "rb") as fh:
        line = fh.readline(1 << 20)
        out["header"] = line.decode("utf-8", errors="replace").rstrip("\r\n").split("\t")
        out["body_start"] = fh.tell()
        out["header_complete"] = line.endswith(b"\n")
    return out


def range_task(args: dict, progress) -> dict:
    """Runs in a child: one byte range of the summary -> accumulators."""
    channels, ratio, header = args["channels"], args["ratio"], args["header"]
    deadline = time.monotonic() + args["remaining_s"] if args.get("remaining_s") is not None else None
    acc = new_acc(channels)
    pieces = workers.Pieces(args["path"], args["lo"], args["hi"], args["body_start"], piece=args.get("piece", PIECE), deadline=deadline)
    for data in pieces:
        cols, bad = parse_block(data, header)
        acc["bad_rows"] += bad
        update(acc, cols, channels, ratio)
        progress(pieces.done)
    return {"acc": acc, "done": pieces.done, "span": max(args["hi"] - pieces.first, 0), "stopped": pieces.stopped,
            "cut_tail": pieces.cut_tail, "cut_bytes": pieces.cut_bytes}


def n_workers(jobs: int, body_bytes: int, gentle: bool) -> int:
    """How many worker processes: 2 for a gently read root, one below 1 GB, else ``min(8, cpus)``; never more
    than one per 8 MB of file."""
    if jobs and jobs > 0:
        n = jobs
    elif gentle:
        n = 2
    elif body_bytes < (1 << 30):
        n = 1
    else:
        n = min(8, os.cpu_count() or 1)
    return max(1, min(n, body_bytes // MIN_BYTES_PER_WORKER + 1))


def run_pass(path: str, probe: dict, *, channels: int, ratio: float, jobs: int = 0, gentle: bool = False,
             deadline: float | None = None, stall_timeout_s: float = 120.0, on_progress=None, piece: int = PIECE) -> dict:
    """The pass itself. Returns ``{"acc", "meta"}``; ``meta["state"]`` is ok, partially_read, truncated or error."""
    t0 = time.monotonic()
    header, body_start, size = probe["header"], probe["body_start"], probe["size"]
    meta: dict = {"size": size, "header": header, "bytes_read": 0, "coverage_frac": 0.0, "complete": False, "state": "error",
                  "errors": [], "workers": 0, "seconds": 0.0, "reads": 0, "cut_tail": False, "bad_rows": 0}
    missing = [c for c in REQUIRED if c not in header]
    if missing:
        meta["errors"].append({"code": "columns_missing", "message": "The sequencing summary lacks the columns " + ", ".join(missing) + "."})
        return {"acc": new_acc(channels), "meta": meta}
    body = max(size - body_start, 0)
    n = n_workers(jobs, body, gentle)
    remaining = (deadline - time.monotonic()) if deadline is not None else None
    tasks = [{"path": path, "lo": lo, "hi": hi, "body_start": body_start, "header": header, "channels": channels, "ratio": ratio,
              "remaining_s": remaining, "piece": piece, "_nice": 10 if gentle else 0}
             for lo, hi in workers.split_ranges(body_start, size, n)]
    results = workers.run_tasks(f"{__name__}:range_task", tasks, stall_timeout_s=stall_timeout_s, on_progress=on_progress)
    acc = new_acc(channels)
    done = cut_bytes = 0
    stopped = cut = False
    for r in results:
        if r.state != "ok":
            code = {"stalled": "stalled", "died": "worker_died"}.get(r.state, "read_error")
            meta["errors"].append({"code": code, "message": (r.error or r.state).splitlines()[0][:300]})
            continue
        v = r.value
        merge_acc(acc, v["acc"])
        done += v["done"]
        cut_bytes += v["cut_bytes"]
        stopped |= v["stopped"]
        cut |= v["cut_tail"]
    reads, bad = int(acc["n"][0]), int(acc["bad_rows"][0])
    meta.update(bytes_read=done, coverage_frac=(done / body if body else 1.0), workers=n, reads=reads, cut_tail=cut, bad_rows=bad,
                seconds=time.monotonic() - t0)
    if meta["errors"]:
        meta["state"] = "error"
    elif bad and bad > BAD_ROW_LIMIT * max(reads + bad, 1):
        meta["state"] = "error"
        meta["errors"].append({"code": "too_many_bad_rows", "message": f"{bad:,} of {reads + bad:,} rows could not be read."})
    elif stopped:
        meta["state"] = "partially_read"
    elif cut or not probe.get("header_complete", True):
        meta["state"] = "truncated"
    else:
        meta["state"] = "ok"
        meta["complete"] = done + cut_bytes == body
        if not meta["complete"]:
            meta["state"] = "truncated"
    return {"acc": acc, "meta": meta}


# ---------------------------------------------------------------------------------------------- fragments
def fragments(acc: dict, run_key: str, *, channels: int, ratio: float) -> list[dict]:
    """Accumulators -> cube fragments in raw units (events, milliseconds) for ``dyn``, ``len.est`` and ``fc``."""
    er_reads = acc["h_dur"].sum(1)
    total_by_hour = acc["ts_n"].sum(0)
    n_hours = max(int(np.nonzero(total_by_hour)[0].max()) + 1 if total_by_hour.sum() else 0, 1)
    dyn_cells, len_cells = [], []
    for k, name in enumerate(END_REASONS):
        if er_reads[k] == 0:
            continue
        cell = {"at": [name],
                "scalars": {"reads": int(er_reads[k]), "minknow_events": int(acc["ts_ev"][k].sum()),
                            "duration_ms": int(acc["ts_dur_ms"][k].sum()), "reads_under_1s": int(acc["n_short_dur"][k])},
                "hists": {"dyn_duration_hist": build.hist("duration_log40", acc["h_dur"][k]),
                          "dyn_speed_est_hist": build.hist("speed_5bps", acc["h_spd"][k], missing=int(acc["n_short_dur"][k]))},
                "series": {"dyn_end_reason_series": build.series(
                    "run_time_3600s", {"reads": acc["ts_n"][k][:n_hours], "minknow_events": acc["ts_ev"][k][:n_hours],
                                       "duration_ms": acc["ts_dur_ms"][k][:n_hours]}, n_bins=n_hours)}}
        if name == "signal_positive":
            m = acc["h_spd_t"][:n_hours]
            cell["hists2d"] = {"dyn_speed_est_series": build.hist2d("run_time_3600s", "speed_5bps", m[:, 1:-1],
                                                                    outside=int(m[:, 0].sum() + m[:, -1].sum()))}
        dyn_cells.append(cell)
        len_cells.append({"at": [name],
                          "scalars": {"reads": int(er_reads[k]), "minknow_events": int(acc["h_len_ev"][k].sum())},
                          "hists": {"len_est_hist": build.hist("read_length_log40", acc["h_len_n"][k], acc["h_len_ev"][k],
                                                               w_name="minknow_events", missing=0)}})
    pt = acc.get("pore_type") or {}
    run_level = {"counts": {"dyn_pore_type": build.counts(sorted(pt), [pt[k] for k in sorted(pt)], role="categories")}} if pt else None
    gs = grid_bin_s(channels)
    used = np.nonzero(acc["grid_ms"].sum(0))[0]
    ncol = int(used.max()) + 1 if len(used) else 1
    fc_cell = {
        "at": [],
        "scalars": {"reads": int(acc["n"][0]), "reads_in_channel_map": int(acc["ch_n"][1:].sum()), "reads_mux_unknown": int(acc["mux_n"][0]),
                    "reads_channel_out_of_range": int(acc["out_of_range_channel"][0])},
        "channels": {
            "fc_channel_map": {"layout": "channel", "values": {
                "reads": [int(v) for v in acc["ch_n"][1:]], "minknow_events": [int(v) for v in acc["ch_ev"][1:]],
                "duration_ms": [int(v) for v in acc["ch_dur_ms"][1:]], "last_read_end_ms": [int(v) for v in acc["max_ch_last_ms"][1:]]},
                "merge": {"last_read_end_ms": "max"}},
            "fc_channel_well_reads": {"layout": "channel", "wells": 4, "values": {"reads": [int(v) for v in acc["chmux_n"][1:, 1:5].ravel()]}}},
        "counts": {"fc_mux_usage": {"cats": ["1", "2", "3", "4"], "n": [int(v) for v in acc["mux_n"][1:5]],
                                    "w": [int(v) for v in acc["mux_ev"][1:5]], "w_name": "minknow_events", "role": "categories",
                                    "other": int(acc["mux_n"][0])}}}
    params = {"events_to_base_ratio": ratio, "channel_count": channels}
    return [
        P.fragment(run_key, "dyn", dims=["end_reason"], cells=dyn_cells, run_level=run_level, params=params,
                   facts={"reads": int(acc["n"][0]), "minknow_events": int(acc["ev_total"][0]), "max_events": int(acc["max_events"][0]),
                          "max_duration_ms": int(acc["max_dur_ms"][0]), "reads_after_time_axis": int(acc["late_reads"][0])}),
        P.fragment(run_key, "len.est", dims=["end_reason"], cells=len_cells, params=params),
        P.fragment(run_key, "fc", dims=[], cells=[fc_cell], params=params,
                   facts={"channel_time": {"bin_s": gs, "rows": channels, "cols": ncol,
                                           "ms": [int(v) for v in acc["grid_ms"][1:, :ncol].ravel()]}}),
    ]


def make_partial(path: str, probe: dict, result: dict, run_key: str, *, channels: int, ratio: float) -> dict:
    """The partial of one pass (``flowcellar.model.partial``), cube fragments in raw units."""
    acc, meta = result["acc"], result["meta"]
    p = P.new_partial(role=ROLE, input_path=path, fp=probe["fingerprint"], size=probe["size"], mtime=from_epoch(probe["mtime"]),
                      modules=MODULES)
    reads = meta["reads"]
    p["state"] = meta["state"]
    if meta["state"] == "ok":
        p["coverage"] = {"kind": "full", "records": reads, "records_total": reads, "fraction": 1.0}
    elif meta["state"] == "truncated":
        p["coverage"] = {"kind": "truncated", "records": reads, "records_total": None, "fraction": round(meta["coverage_frac"], 6)}
    else:
        p["coverage"] = {"kind": "byte_ranges", "records": reads, "records_total": None, "fraction": round(meta["coverage_frac"], 6)}
    p["routing"] = {"runs": {run_key: reads}, "unrouted": 0, "default_run": run_key}
    p["fragments"] = fragments(acc, run_key, channels=channels, ratio=ratio) if meta["state"] != "error" else []
    p["facts"] = {"reads": reads, "minknow_events": int(acc["ev_total"][0]), "events_to_base_ratio": ratio, "channel_count": channels,
                  "columns": len(meta["header"]), "bad_rows": meta["bad_rows"], "cut_tail": bool(meta["cut_tail"]),
                  "bytes_read": meta["bytes_read"], "mtime_ns": probe.get("mtime_ns")}
    if meta["state"] != "error":
        p["facts"]["exact"] = exact_lengths(acc.get("ev_count", ()))      # N50 and median without a histogram between
    p["errors"] = list(meta["errors"])
    secs = meta["seconds"]
    p["timings"] = {"seconds": round(secs, 3), "mb_per_s": round(meta["bytes_read"] / 1e6 / secs, 1) if secs > 0 else None,
                    "workers": meta["workers"]}
    return p
