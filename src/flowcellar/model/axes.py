"""The axis registry: an axis id never changes meaning.

``AXES`` holds the definitions as they are copied into a bundle; ``EDGES`` the literal float64 edges used for
binning. There is one binning implementation: ``numpy.searchsorted(edges, x, side="right")`` on these tables.
Run-time axes are linear and open-ended, so they have a step and no edge table.
"""
from __future__ import annotations

import numpy as np

ORIGIN_TEXT = "start of data reading in the sequencing acquisition (identity.protocol_run.time_origin)"
CLOCK_TEXT = ("the whole step of the UTC clock that holds the run's time origin; the section's params.clock gives its "
              "time and its distance from the time origin")


def _log_edges(lo_exp: int, hi_exp: int, per_decade: int) -> np.ndarray:
    i = np.arange(lo_exp * per_decade, hi_exp * per_decade + 1)
    e = np.power(10.0, i / per_decade)
    whole = (i % per_decade) == 0
    e[whole] = np.power(10.0, (i[whole] // per_decade).astype(np.float64))     # whole decades are exact
    return e


def _lin_edges(lo: float, hi: float, inv_step: float) -> np.ndarray:
    return np.arange(int(round(lo * inv_step)), int(round(hi * inv_step)) + 1) / float(inv_step)


def _rt(step: int, coarser: list[str]) -> dict:
    return {"type": "linear", "start": 0, "step": step, "unit": "s", "family": "run_time", "coarsens_to": coarser,
            "origin": ORIGIN_TEXT}


def _lin(start, step, n, unit, family, coarser=None, **extra) -> dict:
    d = {"type": "linear", "start": start, "step": step, "n": n, "unit": unit, "family": family}
    if coarser:
        d["coarsens_to"] = coarser
    d.update(extra)
    return d


def _log(min_exp, per_decade, n, unit, family, coarser=None) -> dict:
    d = {"type": "log10", "min": min_exp, "bins_per_decade": per_decade, "n": n, "unit": unit, "family": family}
    if coarser:
        d["coarsens_to"] = coarser
    return d


AXES: dict[str, dict] = {
    # run time, from the time origin, open-ended
    "run_time_600s": _rt(600, ["run_time_1800s", "run_time_3600s", "run_time_7200s", "run_time_86400s"]),
    "run_time_1800s": _rt(1800, ["run_time_3600s", "run_time_7200s", "run_time_86400s"]),
    "run_time_3600s": _rt(3600, ["run_time_7200s", "run_time_86400s"]),
    "run_time_7200s": _rt(7200, ["run_time_86400s"]),
    "run_time_86400s": _rt(86400, []),
    # clock time, for what is binned from the start time each read carries (the st tag): bins begin on whole steps of
    # the UTC clock; bin 0 is the step that holds the run's time origin, so on the run-time axis every bin lies
    # ``params.clock.offset_s`` seconds (between minus one step and zero) from a multiple of the step
    "clock_600s": {"type": "linear", "start": 0, "step": 600, "unit": "s", "family": "clock_time", "coarsens_to": ["clock_3600s"],
                   "origin": CLOCK_TEXT},
    "clock_3600s": {"type": "linear", "start": 0, "step": 3600, "unit": "s", "family": "clock_time", "coarsens_to": [],
                    "origin": CLOCK_TEXT},
    # index axes
    "mux_scan_index": {"type": "index", "unit": "pore scan number, 1-based", "family": "index"},
    "channel": {"type": "index", "unit": "channel number, 1-based", "family": "channel"},
    "genome_window": {"type": "index", "unit": "window along the reference, 0-based; the section's params.reference "
                      "lists the contigs, their first window and the window length", "family": "genome_window"},
    # read length: 10 b to 10 Mb
    "read_length_log40": _log(1, 40, 240, "bases", "read_length", ["read_length_log20", "read_length_log4"]),
    "read_length_log20": _log(1, 20, 120, "bases", "read_length", ["read_length_log4"]),
    "read_length_log4": _log(1, 4, 24, "bases", "read_length"),
    # quality
    "qscore_0p25": _lin(0, 0.25, 240, "phred", "qscore", ["qscore_1", "qscore_5"]),
    "qscore_1": _lin(0, 1, 60, "phred", "qscore", ["qscore_5"]),
    "qscore_5": _lin(0, 5, 12, "phred", "qscore"),
    # duration: 0.1 s to 100,000 s
    "duration_log40": _log(-1, 40, 240, "s", "duration"),
    "speed_5bps": _lin(0, 5, 200, "bases/s", "speed"),
    # alignment identity on a Phred scale; identity 1.0 is the overflow cell
    "identity_q_0p1": _lin(0, 0.1, 500, "phred", "identity_q", ["identity_q_0p5", "identity_q_1"]),
    "identity_q_0p5": _lin(0, 0.5, 100, "phred", "identity_q", ["identity_q_1"]),
    "identity_q_1": _lin(0, 1, 50, "phred", "identity_q"),
    "mapq_1": _lin(0, 1, 61, "mapq", "mapq", note="MAPQ 255 is the overflow cell"),
    "polya_1nt": _lin(0, 1, 1000, "nt", "polya", ["polya_5nt"]),
    "polya_5nt": _lin(0, 5, 200, "nt", "polya"),
    "fraction_0p01": _lin(0, 0.01, 100, "fraction", "fraction", closed_last=True),
    "rel_pos_0p01": _lin(0, 0.01, 100, "relative position", "rel_pos", closed_last=True),
    "gc_1pct": _lin(0, 1, 100, "percent", "gc", closed_last=True),
    "current_2pa": _lin(0, 2, 250, "pA", "current"),
    "noise_0p1pa": _lin(0, 0.1, 100, "pA", "noise"),
    "trimmed_bases_2": _lin(0, 2, 150, "bases", "trimmed_bases"),
    "count_1": _lin(0, 1, 21, "count", "count"),
    "per_kb_1": _lin(0, 1, 100, "events per kb", "per_kb"),
}

EDGES: dict[str, np.ndarray] = {
    "read_length_log40": _log_edges(1, 7, 40),
    "read_length_log20": _log_edges(1, 7, 40)[::2],
    "read_length_log4": _log_edges(1, 7, 40)[::10],
    "qscore_0p25": _lin_edges(0, 60, 4),
    "qscore_1": _lin_edges(0, 60, 4)[::4],
    "qscore_5": _lin_edges(0, 60, 4)[::20],
    "duration_log40": _log_edges(-1, 5, 40),
    "speed_5bps": _lin_edges(0, 1000, 1)[::5],
    "identity_q_0p1": _lin_edges(0, 50, 10),
    "identity_q_0p5": _lin_edges(0, 50, 10)[::5],
    "identity_q_1": _lin_edges(0, 50, 10)[::10],
    "mapq_1": _lin_edges(0, 61, 1),
    "polya_1nt": _lin_edges(0, 1000, 1),
    "polya_5nt": _lin_edges(0, 1000, 1)[::5],
    "fraction_0p01": _lin_edges(0, 1, 100),
    "rel_pos_0p01": _lin_edges(0, 1, 100),
    "gc_1pct": _lin_edges(0, 100, 1),
    "current_2pa": _lin_edges(0, 500, 1)[::2],
    "noise_0p1pa": _lin_edges(0, 10, 10),
    "trimmed_bases_2": _lin_edges(0, 300, 1)[::2],
    "count_1": _lin_edges(0, 21, 1),
    "per_kb_1": _lin_edges(0, 100, 1),
}
for _e in EDGES.values():
    _e.setflags(write=False)


def axis_def(axis: str) -> dict:
    """The definition copied into ``bundle.axes``. KeyError for an axis that is not registered."""
    return dict(AXES[axis])


def has_edges(axis: str) -> bool:
    return axis in EDGES


def n_bins(axis: str) -> int:
    return len(EDGES[axis]) - 1


def n_cells(axis: str) -> int:
    """Cells of an internal histogram array: underflow, the bins, overflow."""
    return len(EDGES[axis]) + 1


def bin_index(x, axis: str) -> np.ndarray:
    """0 = underflow, k = edges[k-1] <= x < edges[k], len(edges) = overflow.

    On an axis marked ``closed_last`` a value equal to the last edge falls in the last bin.
    """
    e = EDGES[axis]
    x = np.asarray(x)
    idx = np.searchsorted(e, x, side="right")
    if AXES[axis].get("closed_last"):
        idx = np.where(x == e[-1], len(e) - 1, idx)
    return idx


def time_bin(seconds, axis: str) -> np.ndarray:
    """Bin of a run-time axis: floor(t / step); negative times go to bin 0."""
    step = AXES[axis]["step"]
    return np.maximum(np.floor_divide(np.asarray(seconds, dtype=np.float64), step).astype(np.int64), 0)


def coarsen_factor(fine: str, coarse: str) -> int:
    """How many bins of ``fine`` make one bin of ``coarse``; ValueError when they are not nested."""
    if fine == coarse:
        return 1
    if coarse not in AXES[fine].get("coarsens_to", []):
        raise ValueError(f"{fine} does not coarsen to {coarse}")
    if AXES[fine]["family"] == "run_time":
        return AXES[coarse]["step"] // AXES[fine]["step"]
    return n_bins(fine) // n_bins(coarse)


def coarsen(full: np.ndarray, factor: int) -> np.ndarray:
    """Internal histogram array (under, bins..., over) re-binned exactly onto an axis ``factor`` times coarser."""
    full = np.asarray(full)
    if factor == 1:
        return full.copy()
    inner = full[1:-1]
    return np.concatenate([full[:1], inner.reshape(-1, factor).sum(1), full[-1:]])


def used_axes(modules: dict) -> dict:
    """Definitions of every axis the sections of a bundle use (``bundle.axes``)."""
    used: set[str] = set()

    def from_cell(c: dict) -> None:
        for h in (c.get("hists") or {}).values():
            used.add(h["axis"])
        for h in (c.get("hists2d") or {}).values():
            used.update((h["x"], h["y"]))
        for s in (c.get("series") or {}).values():
            used.add(s["axis"])
        for g in (c.get("grids") or {}).values():
            for k in ("x", "y", "rows_axis", "cols_axis"):
                if isinstance(g.get(k), str) and g[k] in AXES:
                    used.add(g[k])
    for sec in modules.values():
        for c in sec.get("cells", []):
            from_cell(c)
        if sec.get("run_level"):
            from_cell(sec["run_level"])
    return {a: axis_def(a) for a in sorted(used) if a in AXES}
