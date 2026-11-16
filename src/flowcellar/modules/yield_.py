"""Section ``yield``: output over time from MinKNOW's throughput file."""
from __future__ import annotations

import numpy as np

from flowcellar.aggregates import build
from flowcellar.modules.common import HORIZONS, VERSIONS, note
from flowcellar.sources.minknow import Throughput, to_float, to_int


def ratio_of(cfg: dict, tp: Throughput) -> tuple[float, str]:
    """(events_to_base_ratio, where it came from). It is not a constant (0.6 R10, 1.8 R9, 0.31 RNA004), so it is read
    from the run: the report's configuration, else the last row of the throughput file."""
    r = to_float((cfg or {}).get("events_to_base_ratio"))
    if r and r > 0:
        return r, "report"
    if tp.usable and len(tp.est_bases) and tp.events[-1] > 0 and tp.est_bases[-1] > 0:
        return round(float(tp.est_bases[-1]) / float(tp.events[-1]), 2), "throughput"
    return 0.6, "assumed"


def build_section(tp: Throughput, ari: dict, *, ratio: float, ratio_source: str, in_rj: str | None, in_tp: str | None,
                  p_first: int | None, act_h: dict) -> tuple[dict | None, dict]:
    """(section or None, derived numbers: run_hours, reads, estimated bases, per-horizon values, the output sparkline)."""
    ys = {k: to_int(v, None) for k, v in (ari.get("yield_summary") or {}).items()} if isinstance(ari.get("yield_summary"), dict) else {}
    derived: dict = {"run_hours": None, "reads": ys.get("read_count"), "selected_events": ys.get("selected_events"),
                     "est_bases": ys.get("estimated_selected_bases"), "bases_bin": None}
    if not tp.usable:
        if not ys:
            return None, derived
        # the report alone gives the totals, not their course over time
        sec = {"version": VERSIONS["yield"], "tier": "run_health", "state": "partial", "partial_reason": "inputs_missing",
               "inputs": [i for i in (in_rj,) if i], "dims": [],
               "cells": [{"at": [], "scalars": {k: v for k, v in (("reads", ys.get("read_count")), ("selected_events", ys.get("selected_events")),
                                                                  ("selected_raw_samples", ys.get("selected_raw_samples"))) if v is not None}}],
               "params": {"events_to_base_ratio": ratio},
               "summary": {"yield_reads": ys.get("read_count"), "yield_estimated_bases": ys.get("estimated_selected_bases")},
               "notes": [note("warning", "throughput_missing", "The throughput file is missing or unreadable"
                              + (": " + "; ".join(tp.problems) if tp.problems else "")
                              + ". Totals come from the report; output over time and at fixed hours is not known.", input=in_tp)]}
        return sec, derived
    notes: list[dict] = []
    last_min = int(tp.minute[-1])
    run_hours = last_min / 60.0
    n10 = last_min // 10 + 1
    # carry the last observation forward over any missing minute
    idx = np.searchsorted(tp.minute, np.arange(last_min + 1), side="right") - 1
    cum_r = tp.reads[np.clip(idx, 0, None)].copy()
    cum_r[idx < 0] = 0
    cum_e = tp.events[np.clip(idx, 0, None)].copy()
    cum_e[idx < 0] = 0
    ends = np.minimum(np.arange(n10) * 10 + 9, last_min)
    inc_r = np.diff(np.concatenate([[0], cum_r[ends]]))
    inc_e = np.diff(np.concatenate([[0], cum_e[ends]]))
    if (inc_r < 0).any() or (inc_e < 0).any():
        notes.append(note("warning", "throughput_not_monotonic",
                          "The throughput file is not cumulative everywhere; negative steps were set to zero.", input=in_tp))
        inc_r, inc_e = np.clip(inc_r, 0, None), np.clip(inc_e, 0, None)
    bases_bin = [round(int(e) * ratio) for e in inc_e]
    cum_b = np.cumsum(bases_bin)
    cum_rb = np.cumsum(inc_r)
    tot_reads, tot_bases_series = int(cum_rb[-1]), int(cum_b[-1])
    est_report = ys.get("estimated_selected_bases") or tot_bases_series

    def at_h(cum, h):
        b = h * 6
        return int(cum[b - 1]) if 0 < b <= len(cum) and run_hours >= h else None
    hourly = [int(cum_b[min(6 * (h + 1), len(cum_b)) - 1] - (cum_b[6 * h - 1] if h else 0)) for h in range((len(cum_b) + 5) // 6)]
    half_bin = int(np.searchsorted(cum_b, 0.5 * cum_b[-1])) if cum_b[-1] > 0 else 0
    summ: dict = {"yield_reads": tot_reads, "yield_estimated_bases": int(est_report)}
    for h in HORIZONS:
        summ[f"yield_reads_at_{h}h"] = at_h(cum_rb, h)
        summ[f"yield_estimated_bases_at_{h}h"] = at_h(cum_b, h)
    b24 = at_h(cum_b, 24)
    summ["yield_first_24h_frac"] = round(b24 / tot_bases_series, 4) if b24 and tot_bases_series else None
    summ["yield_half_time_h"] = round((half_bin + 1) / 6, 2) if cum_b[-1] > 0 else None
    summ["yield_rate_peak_gb_per_h"] = round(max(hourly) / 1e9, 3) if hourly else None
    summ["yield_per_pore_mb"] = round(est_report / p_first / 1e6, 2) if p_first else None
    for h in (24, 48, 72):
        v = summ.get(f"yield_estimated_bases_at_{h}h")
        summ[f"yield_per_pore_mb_{h}h"] = round(v / p_first / 1e6, 2) if v and p_first else None
    for h in (24, 48):
        v = summ.get(f"yield_estimated_bases_at_{h}h")
        summ[f"yield_kb_per_active_channel_hour_{h}h"] = round(v / act_h[h] / 1e3, 1) if v and act_h.get(h) else None
    if len(tp.est_bases):
        notes.append(note("info", "estimated_bases_sources_differ",
                          f"Three figures for estimated bases exist: {int(est_report):,} in the report (shown as the total), "
                          f"{tot_bases_series:,} from selected events x {ratio} (the stored series), and {int(tp.est_bases[-1]):,} "
                          "in the last row of the throughput file.", input=in_tp))
    if ratio_source != "report":
        notes.append(note("warning", "events_to_base_ratio_not_in_report",
                          f"The report does not give the number of bases per event; {ratio} is used, "
                          + ("from the last row of the throughput file." if ratio_source == "throughput" else "the value of the current chemistry."),
                          input=in_tp))
    if tp.state != "ok":
        notes.append(note("warning", "throughput_incomplete", "The throughput file is incomplete: " + "; ".join(tp.problems) + ".", input=in_tp))
    raw = ys.get("selected_raw_samples") or (int(tp.raw_samples[-1]) if len(tp.raw_samples) else None)
    scalars = {"reads": tot_reads, "selected_events": int(inc_e.sum())}
    if raw is not None:
        scalars["selected_raw_samples"] = int(raw)
    sec = {
        "version": VERSIONS["yield"], "tier": "run_health", "state": "ok" if tp.state == "ok" else "partial",
        "inputs": [i for i in (in_rj, in_tp) if i], "dims": [],
        "cells": [{"at": [], "scalars": scalars,
                   "series": {"yield_series": build.series("run_time_600s", {"reads": inc_r, "selected_events": inc_e,
                                                                             "estimated_bases": bases_bin}, n_bins=n10)}}],
        "params": {"events_to_base_ratio": ratio, "events_to_base_ratio_source": ratio_source,
                   "estimated_bases_rule": "round(selected_events x events_to_base_ratio) per stored cell",
                   "horizon_rule": "value at h hours = cumulative value in the throughput row of minute 60 h - 1"},
        "summary": summ,
        "value_sources": {"yield_estimated_bases": {"input": in_rj if ys.get("estimated_selected_bases") else in_tp,
                                                    "kind": "data" if ys.get("estimated_selected_bases") else "computed"}},
        "notes": notes,
    }
    if sec["state"] == "partial":
        sec["partial_reason"] = "partly_read"
    derived.update(run_hours=run_hours, reads=tot_reads, selected_events=int(inc_e.sum()), est_bases=int(est_report),
                   bases_bin=bases_bin, report_read_count=ys.get("read_count"), report_selected_events=ys.get("selected_events"),
                   throughput_reads=int(tp.reads[-1]), throughput_events=int(tp.events[-1]))
    return sec, derived
