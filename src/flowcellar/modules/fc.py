"""Section ``fc``: flow cell health.

Quick detail comes from the pore-activity file and the report (channel states, pore scans, temperature); full
detail adds what the sequencing summary and the pore-scan table give (channel map, channel x time grid, wells,
open-pore current and noise, voltage per scan).
"""
from __future__ import annotations

import math
import statistics

import numpy as np

from flowcellar.aggregates import build, derive
from flowcellar.modules.common import CHANNEL_STATES, SCAN_CATS, VERSIONS, grid_bin_s, note, r4
from flowcellar.sources.minknow import PoreActivity, g, to_float, to_int
from flowcellar.util.timeutil import parse_ts

STOP_SAMPLE_WINDOW_S = 120.0       # a temperature sample this close to the end of data reading was taken as the run stopped
_TEMP_AXES = {120: "run_time_7200s", 60: "run_time_3600s", 30: "run_time_1800s", 10: "run_time_600s"}


def _temperature(seq: dict, bream: dict, duration_s: float | None) -> dict:
    """Temperature samples of the sequencing acquisition, without the sample taken as the run stops."""
    tt = bream.get("target_temperature") if isinstance(bream.get("target_temperature"), dict) else {}
    bucket_min = to_float(seq.get("bucket_size")) or 0
    out = {"fc": [], "ch": [], "lo": [], "hi": [], "minute": [], "dev_max": 0.0, "n": 0, "dropped": 0, "blocks": 0,
           "bucket_min": bucket_min}
    blocks = [b for b in (seq.get("temperature") or []) if isinstance(b, dict)]
    out["blocks"] = len(blocks)
    for blk in blocks:
        samples = [s for s in (blk.get("temperatures") or []) if isinstance(s, dict)]
        for i, s in enumerate(samples):
            last = i == len(samples) - 1
            if last:
                # the sample MinKNOW takes while it stops the run reads low and says nothing about the run
                t_nominal = i * bucket_min * 60.0
                at_stop = (t_nominal >= duration_s - STOP_SAMPLE_WINDOW_S) if (duration_s and bucket_min) else True
                if at_stop and len(samples) > 1:
                    out["dropped"] += 1
                    continue
            v = c = None
            for k in ("promethion", "minion", "flongle"):
                if isinstance(s.get(k), dict):
                    v = to_float(s[k].get("flowcell_temperature"))
                    if v is None:
                        v = to_float(s[k].get("heatsink_temperature"))
                    c = to_float(s[k].get("chamber_temperature"))
            tg = tt or (s.get("target_temperature") if isinstance(s.get("target_temperature"), dict) else {})
            lo, hi = to_float(tg.get("minimum")), to_float(tg.get("maximum"))
            out["fc"].append(None if v is None else round(v, 2))
            out["ch"].append(None if c is None else round(c, 2))
            out["lo"].append(lo)
            out["hi"].append(hi)
            out["minute"].append(i * bucket_min)
            if v is None or lo is None or hi is None:
                continue
            tol = 0.0 if tt else 0.1          # a single target value gets the width MinKNOW's own band has
            out["dev_max"] = max(out["dev_max"], (lo - tol) - v, v - (hi + tol), 0.0)
            out["n"] += 1
    return out


def _scans(bream: dict, msgs: list, t0s: float | None) -> dict:
    """Pores per scan from ``bream_info.mux_scan_results``; from the scan messages when the report has no table."""
    scans = [s for s in (bream.get("mux_scan_results") or []) if isinstance(s, dict)]
    vals: dict[str, list[int]] = {c: [] for c in SCAN_CATS}
    t: list[int] = []
    for sc in scans:
        c = sc.get("counts") if isinstance(sc.get("counts"), dict) else {}
        other = sum(to_int(v) for k, v in c.items() if k not in SCAN_CATS)
        for k in SCAN_CATS:
            vals[k].append(to_int(c.get(k, 0)) + (other if k == "other" else 0))
        t.append(to_int(sc.get("mux_scan_timestamp", 0)))
    scan_msgs = [m for m in msgs if isinstance(m, dict) and m.get("identifier") == "mux_scan_result"]
    msg0_total = msg0_num = None
    if scan_msgs:
        ex = scan_msgs[0].get("extra_data") if isinstance(scan_msgs[0].get("extra_data"), dict) else {}
        msg0_total, msg0_num = to_int(ex.get("total_pores"), None), to_int(ex.get("num_pores"), None)
    from_messages = False
    if not scans and scan_msgs and msg0_total is not None:
        from_messages = True
        for m in scan_msgs:
            ex = m.get("extra_data") if isinstance(m.get("extra_data"), dict) else {}
            vals["single_pore"].append(to_int(ex.get("total_pores")))
            for k in SCAN_CATS[1:]:
                vals[k].append(0)
            ts = parse_ts(m.get("time"))
            t.append(max(int(ts - t0s), 0) if ts is not None and t0s is not None else 0)
    pores = [a + b for a, b in zip(vals["single_pore"], vals["reserved_pore"])]
    return {"vals": vals, "t": t, "pores": pores, "msg0_num": msg0_num, "msg0_total": msg0_total, "from_messages": from_messages}


def build_quick(pa: PoreActivity, rep: dict, seq: dict | None, *, channels: int | None, run_hours: float | None,
                in_rj: str | None, in_pa: str | None, waiting: bool) -> tuple[dict | None, dict]:
    """(the section in quick detail or None, numbers other sections and the status engine need)."""
    ari = (seq or {}).get("acquisition_run_info") or {}
    cfg = ari.get("config_summary") if isinstance(ari.get("config_summary"), dict) else {}
    bream = ari.get("bream_info") if isinstance(ari.get("bream_info"), dict) else {}
    msgs = [m for m in (rep.get("user_messages") or []) if isinstance(m, dict)]
    t0s = parse_ts(ari.get("data_read_start_time") or ari.get("start_time"))
    t1s = parse_ts(ari.get("data_read_end_time") or ari.get("end_time"))
    duration_s = (t1s - t0s) if t0s is not None and t1s is not None and t1s > t0s else None
    derived: dict = {"p_first": None, "p_last": None, "act_h": {}, "act10": None, "rate": None, "pa_hours": None}
    if not pa.usable and not seq:
        return None, derived
    notes: list[dict] = []
    series: dict = {}
    tables: dict = {}
    scalars: dict = {}
    summary: dict = {}
    params: dict = {}

    # ---------------------------------------------------------------- channel states
    if pa.usable:
        st = pa.states
        last = pa.last_minute
        rate = to_int(cfg.get("sample_rate"), 0)
        if not rate:
            per_min = float(np.median(sum(st.values()))) if st else 0.0
            rate = int(round(per_min / (channels * 60))) if channels and per_min else 5000
            notes.append(note("info", "sample_rate_inferred", f"The report gives no sample rate; {rate:,} per second is used, "
                              "from the samples per minute in the pore-activity file." if channels and per_min else
                              "The report gives no sample rate; 5,000 per second is assumed.", input=in_pa))
        derived["rate"], derived["pa_hours"] = rate, last / 60.0
        if run_hours is None:
            run_hours = last / 60.0
        states = list(CHANNEL_STATES) + sorted(s for s in st if s not in CHANNEL_STATES)
        zero = np.zeros(last + 1, np.int64)

        def S(*names):
            out = zero.copy()
            for n_ in names:
                if n_ in st:
                    out = out + st[n_]
            return out
        seq_m, pore_m, rec_m = S("strand", "adapter"), S("pore"), S("unavailable", "unblocking")
        tot_m = S(*st.keys())
        work_m = seq_m + pore_m + rec_m
        nb = last // 10 + 1

        def by10(a):
            pad = np.zeros(nb * 10, np.int64)
            pad[: len(a)] = a
            return pad.reshape(nb, 10).sum(1)
        per10 = {s: by10(st[s]) if s in st else np.zeros(nb, np.int64) for s in states}
        cs = {s: [0] * nb for s in states}
        for b in range(nb):
            exact = [per10[s][b] / rate for s in states]
            total = round(sum(int(per10[s][b]) for s in states) / rate)
            for s, v in zip(states, build.largest_remainder(exact, total)):
                cs[s][b] = v
        seq10, pore10 = by10(seq_m), by10(pore_m)
        with np.errstate(divide="ignore", invalid="ignore"):
            occ10 = np.where(seq10 + pore10 > 0, seq10 / np.maximum(seq10 + pore10, 1), np.nan)
        occ_valid = occ10[~np.isnan(occ10)]
        occ24 = occ10[:144][~np.isnan(occ10[:144])]
        sq, pr = int(seq_m[30:90].sum()), int(pore_m[30:90].sum())
        occ_first = sq / (sq + pr) if sq + pr > 0 else None
        tot_all, seq_all, work_all = int(tot_m.sum()), int(seq_m.sum()), int(work_m.sum())
        late = np.arange(last + 1) >= last * 2 / 3
        rec_last, work_last = int(rec_m[late].sum()), int(work_m[late].sum())
        act10 = by10(seq_m + pore_m) / (rate * 600.0)
        roll = np.array([act10[i:i + 6].mean() for i in range(max(len(act10) - 5, 0))])
        peak = half_life = None
        if len(roll):
            early = roll[:36]
            pk = int(np.argmax(early))
            peak = float(early[pk])
            for i in range(pk, len(roll)):
                if (roll[i:] < 0.5 * peak).all():
                    half_life = round((i - pk) / 6, 2)
                    break
        for h in (24, 48):
            derived["act_h"][h] = float((seq_m + pore_m)[: h * 60].sum()) / rate / 3600 if run_hours >= h else None
        derived["act10"] = act10
        adapter_all = int(S("adapter").sum())
        series["fc_channel_state_series"] = build.series("run_time_600s", {s: cs[s] for s in states}, unit="channel_seconds", n_bins=nb)
        scalars.update(strand_channel_seconds=int(sum(cs["strand"])), adapter_channel_seconds=int(sum(cs["adapter"])),
                       pore_channel_seconds=int(sum(cs["pore"])),
                       active_channel_seconds=int(round((int(seq_m.sum()) + int(pore_m.sum())) / rate)),
                       recovering_channel_seconds=int(sum(cs["unavailable"]) + sum(cs["unblocking"])),
                       total_channel_seconds=int(sum(sum(v) for v in cs.values())))
        summary.update(
            fc_occupancy_first_hour=r4(occ_first),
            fc_occupancy_median=r4(statistics.median(occ_valid)) if len(occ_valid) else None,
            fc_occupancy_median_24h=r4(np.median(occ24)) if len(occ24) else None,
            fc_sequencing_duty=r4(seq_all / tot_all) if tot_all else None,
            fc_adapter_frac=r4(adapter_all / seq_all) if seq_all else None,
            fc_recovering_frac=r4(int(rec_m.sum()) / work_all) if work_all else None,
            fc_recovering_last_third_frac=r4(rec_last / work_last) if work_last else None,
            fc_saturated_time_frac=r4(int(S("saturated").sum()) / tot_all) if tot_all else None,
            fc_active_channels_peak=round(peak) if peak is not None else None,
            fc_duty_half_life_h=half_life)
        params.update(sample_rate_hz=rate,
                      channel_seconds_rule="samples / sample rate, largest-remainder rounding per time bin so that every bin keeps its exact total",
                      occupancy_first_hour_window_min=[30, 90])
        if pa.state != "ok":
            notes.append(note("warning", "pore_activity_incomplete", "The pore-activity file is incomplete: " + "; ".join(pa.problems) + ".",
                              input=in_pa))

    # ---------------------------------------------------------------- pore scans
    if seq:
        sc = _scans(bream, msgs, t0s)
        pores, t = sc["pores"], sc["t"]
        p_first = pores[0] if pores else None
        p_last = pores[-1] if pores else None
        derived.update(p_first=p_first, p_last=p_last, msg0_num=sc["msg0_num"])
        pore_half = None
        survival: dict[int, float | None] = {}
        if pores and p_first:
            for i in range(len(pores) - 1):
                a, b = pores[i], pores[i + 1]
                if a >= 0.5 * p_first > b:
                    pore_half = round((t[i] + (t[i + 1] - t[i]) * (a - 0.5 * p_first) / (a - b)) / 3600, 1)
                    break
            for h in (24, 48, 72):
                j = min(range(len(pores)), key=lambda k: abs(t[k] - h * 3600))
                survival[h] = pores[j] / p_first if abs(t[j] - h * 3600) < 2 * 3600 else None
        if sc["from_messages"]:
            notes.append(note("info", "pore_scans_from_messages",
                              "The report holds no pore-scan table; pores per scan come from the scan messages and the other scan results are not known.",
                              input=in_rj))
        if pores:
            series["fc_pore_scan_series"] = build.series("mux_scan_index", {**sc["vals"], "t_s": t}, unit="wells",
                                                         levels={"t_s": "same"}, n_bins=len(pores))
        num = sc["msg0_num"]
        summary.update(
            fc_pores_first_scan=p_first, fc_pores_last_scan=p_last, fc_pore_half_life_h=pore_half,
            fc_channels_with_pore_first_scan=num,
            fc_channels_with_pore_first_scan_frac=r4(num / channels) if num is not None and channels else None,
            fc_pore_survival_24h_frac=r4(survival.get(24)), fc_pore_survival_48h_frac=r4(survival.get(48)),
            fc_pore_survival_72h_frac=r4(survival.get(72)))
        scalars["pore_scans"] = len(pores)

        # ------------------------------------------------------------ temperature
        tmp = _temperature(seq, bream, duration_s)
        scalars["temperature_samples"] = tmp["n"]
        summary["fc_temperature_max_dev_c"] = round(tmp["dev_max"], 2) if tmp["n"] else None
        axis = _TEMP_AXES.get(int(tmp["bucket_min"])) if tmp["bucket_min"] else None
        if tmp["fc"] and axis and tmp["blocks"] == 1:
            series["fc_temperature_series"] = build.series(
                axis, {"flowcell_c": tmp["fc"], "chamber_c": tmp["ch"], "target_min_c": tmp["lo"], "target_max_c": tmp["hi"]},
                unit="degC", levels={"flowcell_c": "mean", "chamber_c": "mean", "target_min_c": "same", "target_max_c": "same"})
        elif tmp["fc"]:
            tables["fc_temperature_table"] = {
                "columns": [{"id": "minute", "type": "float", "unit": "min"}, {"id": "flowcell_c", "type": "float", "unit": "degC"},
                            {"id": "chamber_c", "type": "float", "unit": "degC"}, {"id": "target_min_c", "type": "float"},
                            {"id": "target_max_c", "type": "float"}],
                "rows": [list(r) for r in zip(tmp["minute"], tmp["fc"], tmp["ch"], tmp["lo"], tmp["hi"])], "merge": "replace"}
        params["temperature_sample_period_min"] = tmp["bucket_min"] or None
        params["temperature_stop_sample_dropped"] = bool(tmp["dropped"])
        if tmp["dropped"]:
            notes.append(note("info", "temperature_stop_sample_dropped",
                              "The temperature sample MinKNOW takes while it stops the run is left out: it reads low and says nothing about the run.",
                              input=in_rj, value=tmp["dropped"]))
        ts_band = bream.get("target_translocation_speed") if isinstance(bream.get("target_translocation_speed"), dict) else {}
        params.update(mux_scan_period_h=to_float(g(bream, "mux_scan_metadata", "auto_mux_scan_period_hours")),
                      target_speed_min_bps=to_float(ts_band.get("minimum")), target_speed_max_bps=to_float(ts_band.get("maximum")))
    if not series and not scalars:
        return None, derived
    notes.insert(0, note("info", "quick_detail",
                         "Made without the sequencing summary and the pore-scan table: the channel map, the channel x time grid, "
                         "the well map and open-pore current are not in this section."))
    cell: dict = {"at": [], "scalars": scalars, "series": series}
    if tables:
        cell["tables"] = tables
    sec = {"version": VERSIONS["fc"], "tier": "run_health", "state": "partial",
           "partial_reason": "waiting" if waiting else "inputs_missing",
           "inputs": [i for i in (in_rj if seq else None, in_pa if pa.usable else None) if i], "dims": [], "cells": [cell],
           "params": params, "summary": summary, "notes": notes}
    return sec, derived


def _dense(h: dict, rows: int, cols: int) -> np.ndarray:
    m = np.zeros((rows, cols), np.int64)
    if h.get("i"):
        m[np.asarray(h["i"]), np.asarray(h["j"])] = np.asarray(h["n"])
    return m


def add_full(sec: dict, *, summary_frag: dict | None, summary_complete: bool, scan_frag: dict | None, ratio: float,
             channels: int | None, layout: str | None, hours: float, in_ss: str | None, in_ps: str | None,
             scan_missing: bool) -> None:
    """Add to a quick ``fc`` section what the two large passes found. ``summary_frag`` and ``scan_frag`` are the
    ``fc`` fragments of the partials, in raw units (events, milliseconds)."""
    cell = sec["cells"][0]
    summ = sec["summary"]
    notes = [n for n in sec.get("notes", []) if n["code"] != "quick_detail"]
    layout = layout or "unknown"
    if summary_frag and summary_frag.get("cells"):
        c = summary_frag["cells"][0]
        sc = c.get("scalars") or {}
        mux = (c.get("counts") or {}).get("fc_mux_usage")
        if mux:
            cell.setdefault("counts", {})["fc_mux_usage"] = build.counts(
                mux["cats"], mux["n"], w=[round(int(v) * ratio) for v in mux.get("w", [])], role="categories")
            cell["counts"]["fc_mux_usage"]["w_name"] = "estimated_bases"
            cell["counts"]["fc_mux_usage"]["other"] = int(mux.get("other", 0))
        cmap = (c.get("channels") or {}).get("fc_channel_map")
        if summary_complete and cmap and channels:
            v = cmap["values"]
            reads = np.asarray(v["reads"], np.int64)
            events = np.asarray(v["minknow_events"], np.int64)
            cell.setdefault("channels", {})["fc_channel_map"] = {
                "layout": layout,
                "values": {"reads": [int(x) for x in reads], "estimated_bases": [round(int(x) * ratio) for x in events],
                           "sequencing_seconds": [int(round(x / 1000)) for x in v["duration_ms"]],
                           "last_read_time": [int(round(x / 1000)) for x in v["last_read_end_ms"]]},
                "merge": {"last_read_time": "max"}}
            wells = (c.get("channels") or {}).get("fc_channel_well_reads")
            if wells:
                cell["channels"]["fc_channel_well_reads"] = {"layout": layout, "wells": wells.get("wells", 4),
                                                             "values": {"reads": [int(x) for x in wells["values"]["reads"]]}}
            ct = (summary_frag.get("facts") or {}).get("channel_time")
            if ct:
                gs = int(ct["bin_s"])
                ncol = int(max(1, min(math.ceil(hours * 3600 / gs), 4096)))
                ms = np.zeros((channels, ncol), np.int64)
                src = np.asarray(ct["ms"], np.int64).reshape(ct["rows"], ct["cols"])
                k = min(ncol, ct["cols"])
                ms[: min(channels, ct["rows"]), :k] = src[:channels, :k]
                secs = np.minimum(ms / 1000.0, gs)
                scale = math.sqrt(gs) / 63.0
                cell.setdefault("grids", {})["fc_channel_time_grid"] = {
                    "rows": channels, "cols": ncol, "rows_axis": "channel", "cols_axis": f"run_time_{gs}s", "layout": {"kind": "matrix"},
                    "channels": {"strand_seconds": {"enc": "u6", "data": build.u6_encode(np.rint(np.sqrt(secs) / scale)),
                                                    "scale": round(scale, 6), "transform": "sqrt", "levels": 64, "unit": "s", "nodata": False}}}
            active = int((reads > 0).sum())
            gini = derive.gini(events)
            cell["scalars"].update(active_channels=active, reads_in_channel_map=int(reads.sum()),
                                   reads_mux_unknown=int(sc.get("reads_mux_unknown", 0)))
            summ.update(fc_active_channels=active, fc_active_channels_frac=r4(active / channels),
                        fc_channel_yield_gini=round(gini, 3) if gini is not None else None)
        if sc.get("reads_channel_out_of_range"):
            notes.append(note("warning", "channel_out_of_range",
                              "Some reads carry a channel number outside the flow cell's channel count; they are counted in every total but not in the channel map.",
                              input=in_ss, value=int(sc["reads_channel_out_of_range"])))
        if in_ss and in_ss not in sec["inputs"]:
            sec["inputs"].append(in_ss)
    if scan_frag and scan_frag.get("cells"):
        c = scan_frag["cells"][0]
        for name in ("fc_open_pore_current", "fc_open_pore_noise"):
            h = (c.get("hists2d") or {}).get(name)
            if h:
                cell.setdefault("hists2d", {})[name] = h
        wm = (c.get("channels") or {}).get("fc_well_map")
        if wm:
            cell.setdefault("channels", {})["fc_well_map"] = {**wm, "layout": layout}
        bv = (c.get("series") or {}).get("fc_bias_voltage_series")
        if bv:
            cell["series"]["fc_bias_voltage_series"] = bv
            vv = [x for x in bv["values"]["voltage_mv"] if x is not None]
            if vv:
                summ.update(fc_bias_voltage_start_mv=vv[0], fc_bias_voltage_end_mv=vv[-1])
        cell["scalars"]["pore_scan_rows"] = int((c.get("scalars") or {}).get("pore_scan_rows", 0))
        tab = (c.get("series") or {}).get("fc_pore_scan_table_series")
        rep_scan = cell["series"].get("fc_pore_scan_series")
        if tab and rep_scan:
            a, b = rep_scan["values"], tab["values"]

            def pores(v):
                return [x + y for x, y in zip(v.get("single_pore", []), v.get("reserved_pore", []))]
            if pores(a) != pores(b):
                notes.append(note("warning", "pore_scan_sources_differ",
                                  "Pores per scan in the report and in the pore-scan table do not agree for every scan; the report is shown.",
                                  input=in_ps))
        elif tab and not rep_scan:
            # no scan results in the report at all: the table is the only source
            cell["series"]["fc_pore_scan_series"] = tab
            pores = [x + y for x, y in zip(tab["values"]["single_pore"], tab["values"]["reserved_pore"])]
            if pores:
                summ.update(fc_pores_first_scan=pores[0], fc_pores_last_scan=pores[-1])
                cell["scalars"]["pore_scans"] = len(pores)
        if in_ps and in_ps not in sec["inputs"]:
            sec["inputs"].append(in_ps)
    have_map = "fc_channel_map" in cell.get("channels", {})
    have_scan = bool(scan_frag and scan_frag.get("cells"))
    if have_map and (have_scan or scan_missing):
        sec["state"] = "ok"
        sec.pop("partial_reason", None)
    else:
        sec["state"] = "partial"
        sec["partial_reason"] = "partly_read" if (summary_frag and not summary_complete) else "waiting"
        if summary_frag and not summary_complete:
            notes.append(note("warning", "partly_read",
                              "The channel map, the channel x time grid and the count of channels with reads are left out: "
                              "they would be wrong from a part of the reads.", input=in_ss))
    sec["notes"] = notes
    sec.setdefault("params", {})["channel_time_grid_bin_s"] = grid_bin_s(channels) if have_map else None
