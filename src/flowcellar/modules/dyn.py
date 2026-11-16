"""Section ``dyn``: read dynamics from the sequencing summary (speed through the pore, why reads ended, time in
the pore). Built from the ``dyn`` fragment of the summary partial, which is in raw units."""
from __future__ import annotations

import numpy as np

from flowcellar.aggregates import build, derive
from flowcellar.model.axes import EDGES, n_cells
from flowcellar.modules.common import VERSIONS, note, r4

VERSION = VERSIONS["dyn"]


def _sum_full(cells: list[dict], name: str, axis: str, key: str = "n") -> np.ndarray:
    out = np.zeros(n_cells(axis), np.int64)
    for c in cells:
        h = (c.get("hists") or {}).get(name)
        if h:
            out += build.hist_full(h, key)
    return out


def _quantile_full(axis: str, full: np.ndarray, q: float) -> float | None:
    """Quantile of an internal histogram array (under, bins, over): the one derivation of ``aggregates.derive``."""
    tot = int(full.sum())
    if tot <= 0:
        return None
    return derive.quantile(build.hist(axis, full), q)


def partly_read_note(state: str, coverage: dict, in_ss: str | None) -> dict | None:
    frac = coverage.get("fraction")
    pct = f"{frac * 100:.1f} % of the file" if isinstance(frac, (int, float)) else "a part of the file"
    if state == "partially_read":
        return note("warning", "partly_read",
                    f"The pass over the sequencing summary stopped at its time limit after {pct}, in ranges spread evenly over the "
                    "whole file. Counts in this section are for that part only; shares and medians are estimates.", input=in_ss, value=frac)
    if state == "truncated":
        return note("warning", "partly_read",
                    "The sequencing summary is cut off or incomplete. Counts in this section are for the reads it holds; "
                    "shares and medians are estimates.", input=in_ss, value=frac)
    return None


def build_section(frag: dict, *, state: str, coverage: dict, ratio: float, hours: float, speed_min: float | None, in_ss: str | None,
                  completeness: dict | None = None, extra_notes: list[dict] | None = None) -> dict:
    raw_cells = frag.get("cells") or []
    cells = []
    n = sum(int(c["scalars"]["reads"]) for c in raw_cells)
    by_reason = {c["at"][0]: c for c in raw_cells}
    for c in raw_cells:
        sc = c["scalars"]
        s = c["series"]["dyn_end_reason_series"]
        v = s["values"]
        cell = {"at": list(c["at"]),
                "scalars": {"reads": int(sc["reads"]), "minknow_events": int(sc["minknow_events"]),
                            "estimated_bases": round(int(sc["minknow_events"]) * ratio), "duration_ms": int(sc["duration_ms"]),
                            "reads_under_1s": int(sc["reads_under_1s"])},
                "hists": dict(c["hists"]),
                "series": {"dyn_end_reason_series": build.series(
                    s["axis"], {"reads": v["reads"], "estimated_bases": [round(int(x) * ratio) for x in v["minknow_events"]],
                                "duration_ms": v["duration_ms"]}, n_bins=len(v["reads"]))}}
        if c.get("hists2d"):
            cell["hists2d"] = dict(c["hists2d"])
        cells.append(cell)
    # ---- derived numbers, all from the stored aggregates
    sp_cell = by_reason.get("signal_positive")
    speed_med = speed_drop = None
    if sp_cell:
        speed_med = derive.median(sp_cell["hists"]["dyn_speed_est_hist"])
        h2 = (sp_cell.get("hists2d") or {}).get("dyn_speed_est_series")
        if h2 and hours >= 4:
            nb = len(EDGES["speed_5bps"]) - 1
            rows = max(int(max(h2["i"], default=0)) + 1, int(hours) + 1)
            m = np.zeros((rows, nb), np.int64)
            if h2["i"]:
                m[np.asarray(h2["i"]), np.asarray(h2["j"])] = np.asarray(h2["n"])
            q4 = max(1, int(hours // 4))
            first, last = m[:q4].sum(0), m[max(0, int(hours) - q4): int(hours)].sum(0)
            pad = np.zeros(1, np.int64)
            m_first = _quantile_full("speed_5bps", np.concatenate([pad, first, pad]), 0.5)
            m_last = _quantile_full("speed_5bps", np.concatenate([pad, last, pad]), 0.5)
            if m_first and m_last:
                speed_drop = 1 - m_last / m_first
    reads = {k: int(c["scalars"]["reads"]) for k, c in by_reason.items()}
    sp, ub, sn = reads.get("signal_positive", 0), reads.get("unblock_mux_change", 0), reads.get("signal_negative", 0)
    spd_all = _sum_full(raw_cells, "dyn_speed_est_hist", "speed_5bps")
    slow = None
    if speed_min and spd_all.sum():
        kb = int(np.searchsorted(EDGES["speed_5bps"], 0.6 * speed_min, side="right"))        # bins wholly below the edge
        slow = float(spd_all[:kb].sum() / spd_all.sum())
    dur_all = _sum_full(raw_cells, "dyn_duration_hist", "duration_log40")
    dur_med = _quantile_full("duration_log40", dur_all, 0.5)
    summary = {"dyn_speed_est_median": round(speed_med, 1) if speed_med is not None else None,
               "dyn_speed_drop_frac": r4(speed_drop),
               "dyn_signal_positive_frac": r4(sp / n) if n else None,
               "dyn_unblock_frac": r4(ub / n) if n else None,
               "dyn_signal_negative_frac": r4(sn / n) if n else None,
               "dyn_end_reason_other_frac": r4((n - sp - ub - sn) / n) if n else None,
               "dyn_duration_median_s": round(dur_med, 2) if dur_med is not None else None,
               "dyn_slow_read_frac": r4(slow),
               "dyn_reads": n}
    complete = state == "ok"
    params = {"events_to_base_ratio": ratio, "speed_rule": "minknow_events x events_to_base_ratio / duration, reads of 1 s or longer",
              "quantile_rule": "medians and N50 are read from the stored histograms (rank interpolation inside a bin)",
              "speed_drop_rule": "1 - median speed of the last quarter / median speed of the first quarter of the run, reads that ended normally",
              "complete": complete, "reads_read": n}
    if completeness is not None:
        params["completeness"] = completeness
    facts = frag.get("facts") or {}
    notes = list(extra_notes or [])
    pr = partly_read_note(state, coverage, in_ss)
    if pr:
        notes.append(pr)
    if facts.get("reads_after_time_axis"):
        notes.append(note("info", "reads_after_time_axis", "Some reads started later than three weeks into the run; they are counted in the last hour.",
                          input=in_ss, value=int(facts["reads_after_time_axis"])))
    sec = {"version": VERSION, "tier": "run_health", "state": "ok" if complete else "partial", "inputs": [in_ss] if in_ss else [],
           "dims": ["end_reason"], "cells": cells, "summary": summary, "params": params, "notes": notes, "coverage": coverage}
    if not complete:
        sec["partial_reason"] = "partly_read"
    if frag.get("run_level"):
        sec["run_level"] = frag["run_level"]
    return sec
