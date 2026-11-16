"""Section ``len.est``: read length estimated from the signal, before any basecalling.

Built from the ``len.est`` fragment of the summary partial: histograms of reads (``n``) and of MinKNOW events
(``w``) on the estimated-length axis. Events become estimated bases here, once.
"""
from __future__ import annotations

import numpy as np

from flowcellar.aggregates import build, derive
from flowcellar.model.axes import EDGES, n_cells
from flowcellar.modules.common import VERSIONS, note, r4
from flowcellar.modules.dyn import partly_read_note

VERSION = VERSIONS["len.est"]
AXIS = "read_length_log40"


def build_section(frag: dict, *, state: str, coverage: dict, ratio: float, in_ss: str | None, extra_notes: list[dict] | None = None,
                  exact: dict | None = None) -> dict:
    """``exact``: the pass's own N50 and median of events per read (``sources.seqsummary.exact_lengths``). A result
    cached by a build that did not keep them gives ``None``: the two numbers are then read from the histogram,
    and the section says so."""
    raw_cells = frag.get("cells") or []
    cells = []
    n_all = np.zeros(n_cells(AXIS), np.int64)
    ev_all = np.zeros(n_cells(AXIS), np.int64)
    events_total = 0
    for c in raw_cells:
        h = c["hists"]["len_est_hist"]
        n_full, ev_full = build.hist_full(h, "n"), build.hist_full(h, "w")
        n_all += n_full
        ev_all += ev_full
        ev = int(c["scalars"]["minknow_events"])
        events_total += ev
        w_full = np.asarray([round(int(v) * ratio) for v in ev_full], np.int64)
        cells.append({"at": list(c["at"]),
                      "scalars": {"reads": int(c["scalars"]["reads"]), "estimated_bases": round(ev * ratio)},
                      "hists": {"len_est_hist": build.hist(AXIS, n_full, w_full, w_name="estimated_bases", missing=0)}})
    n = int(n_all.sum())
    # N50 does not depend on the constant ratio between events and bases, so it is read from the event sums
    n50 = derive.nx(build.hist(AXIS, n_all, ev_all)) if ev_all.sum() > 0 else None
    median = derive.median(build.hist(AXIS, n_all)) if n else None
    source = "histogram" if n else None
    if n and exact and exact.get("reads") == n and not exact.get("over_cap") and exact.get("events_median") is not None:
        median = float(exact["events_median"]) * ratio
        if exact.get("events_n50") is not None:
            n50 = float(exact["events_n50"]) * ratio
        source = "exact"
    k200 = int(np.searchsorted(EDGES[AXIS], 200.0, side="right"))
    summary = {"len_est_n50": round(n50) if n50 else None,
               "len_est_median": round(median) if median is not None else None,
               "len_est_mean": round(events_total * ratio / n, 1) if n else None,
               "dyn_short_read_est_frac": r4(float(n_all[:k200].sum() / n)) if n else None}
    complete = state == "ok"
    notes = list(extra_notes or [])
    pr = partly_read_note(state, coverage, in_ss)
    if pr:
        notes.append(pr)
    if source == "histogram":
        notes.append(note("info", "quantiles_from_histogram",
                          "Estimated read N50 and median are read from a histogram with 40 steps per tenfold length: good to "
                          "a few per cent, not exact. Counting the reads of this run again makes them exact.", input=in_ss))
    sec = {"version": VERSION, "tier": "run_health", "state": "ok" if complete else "partial", "inputs": [in_ss] if in_ss else [],
           "dims": ["end_reason"], "cells": cells, "summary": summary,
           "params": {"events_to_base_ratio": ratio, "complete": complete, "reads_read": n, "quantile_source": source,
                      "quantile_rule": ("N50 and median are exact: counted per number of events in the pass over the sequencing summary"
                                        if source == "exact" else
                                        "N50 and median are read from the stored histogram by rank inside a bin (within a few per cent)"),
                      "estimated_bases_rule": "round(minknow_events x events_to_base_ratio) per stored cell"},
           "notes": notes, "coverage": coverage}
    if not complete:
        sec["partial_reason"] = "partly_read"
    return sec
