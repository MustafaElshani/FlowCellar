"""The three fixed axes of the add-on, and its aggregates in their stored form.

The edges are the same for every analysis, so histograms of two analyses can be overlaid and added. The
definitions are copied into ``bundle.axes`` of every analysis bundle (a bundle describes itself); they are kept
here and not in the core's axis registry so that the core does not grow.

    log_counts_8_per_decade   8 bins per decade from 1 to 1,000,000; a count of 0 is the underflow cell
    fraction_50               50 bins of 0.02 from 0 to 1 (1.0 falls in the last bin)
    reads_per_molecule        1 to 64 one by one, then one bin per power of two up to 2^20
"""
from __future__ import annotations

import base64

import numpy as np

LOG_COUNTS = "log_counts_8_per_decade"
FRACTION_50 = "fraction_50"
READS_PER_MOLECULE = "reads_per_molecule"


def _log_edges() -> np.ndarray:
    i = np.arange(0, 6 * 8 + 1)
    e = np.power(10.0, i / 8.0)
    whole = (i % 8) == 0
    e[whole] = np.power(10.0, (i[whole] // 8).astype(np.float64))        # whole decades are exact
    return e


EDGES: dict[str, np.ndarray] = {
    LOG_COUNTS: _log_edges(),
    FRACTION_50: np.arange(0, 51) / 50.0,
    READS_PER_MOLECULE: np.concatenate([np.arange(1, 66), 2 ** np.arange(7, 21)]).astype(np.float64),
}
for _e in EDGES.values():
    _e.setflags(write=False)

AXES: dict[str, dict] = {
    LOG_COUNTS: {"type": "log10", "min": 0, "bins_per_decade": 8, "n": 48, "unit": "count", "family": "log_counts",
                 "note": "a unit with a count of 0 is in the underflow cell"},
    FRACTION_50: {"type": "linear", "start": 0, "step": 0.02, "n": 50, "unit": "fraction", "family": "fraction",
                  "closed_last": True},
    READS_PER_MOLECULE: {"type": "edges", "edges": [float(x) for x in EDGES[READS_PER_MOLECULE]], "unit": "reads",
                         "family": "reads_per_molecule",
                         "note": "1 to 64 reads one by one, then one bin per power of two; w holds the reads of each bin"},
}


def axis_defs(names) -> dict:
    return {n: dict(AXES[n]) for n in sorted(set(names))}


def bin_counts(values, axis: str, weights=None) -> np.ndarray:
    """Internal array (underflow, the bins, overflow) of ``values`` on ``axis``; with ``weights`` their sums."""
    e = EDGES[axis]
    v = np.asarray(values)
    idx = np.searchsorted(e, v, side="right")
    if AXES[axis].get("closed_last"):
        idx = np.where(v == e[-1], len(e) - 1, idx)
    if weights is None:
        return np.bincount(idx, minlength=len(e) + 1).astype(np.int64)
    return np.rint(np.bincount(idx, weights=np.asarray(weights, dtype=np.float64), minlength=len(e) + 1)).astype(np.int64)


def hist_from_full(axis: str, full, w_full=None, *, w_name: str | None = None, missing: int | None = None,
                   vmin=None, vmax=None) -> dict:
    """Internal array -> stored histogram with empty ends trimmed (the layout of the bundle schema)."""
    full = np.asarray(full)
    inner = full[1:-1]
    nz = np.nonzero(inner)[0]
    a, b = (int(nz[0]), int(nz[-1]) + 1) if len(nz) else (0, 0)
    h: dict = {"axis": axis, "start": a, "n": [int(x) for x in inner[a:b]]}
    if w_full is not None:
        w_full = np.asarray(w_full)
        h["w"] = [int(x) for x in w_full[1:-1][a:b]]
        if w_name:
            h["w_name"] = w_name
        h["w_under"], h["w_over"] = int(w_full[0]), int(w_full[-1])
    h["under"], h["over"] = int(full[0]), int(full[-1])
    if missing is not None:
        h["missing"] = int(missing)
    if vmin is not None:
        h["min"] = vmin
    if vmax is not None:
        h["max"] = vmax
    return h


def hist(axis: str, values, *, missing: int | None = None) -> dict:
    """Stored histogram of ``values`` (one per unit)."""
    v = np.asarray(values)
    vmin = float(v.min()) if len(v) else None
    vmax = float(v.max()) if len(v) else None
    return hist_from_full(axis, bin_counts(v, axis), missing=missing, vmin=_num(vmin), vmax=_num(vmax))


def _num(x):
    if x is None:
        return None
    return int(x) if float(x).is_integer() else float(x)


def reads_per_molecule_hist(counts_by_c) -> dict:
    """``counts_by_c[c]`` = molecules seen ``c`` times -> stored histogram with ``w`` = the reads of each bin.

    With ``w`` the mean number of reads of the molecules in a wide bin is known, and the total number of reads is
    exact, so the saturation curve drawn from the histogram ends exactly at the saturation of the data.
    """
    h = np.asarray(counts_by_c, dtype=np.int64)
    c = np.arange(len(h), dtype=np.int64)
    present = np.nonzero(h)[0]
    if len(present) == 0:
        return hist_from_full(READS_PER_MOLECULE, np.zeros(len(EDGES[READS_PER_MOLECULE]) + 1, np.int64),
                              np.zeros(len(EDGES[READS_PER_MOLECULE]) + 1, np.int64), w_name="reads")
    idx = np.searchsorted(EDGES[READS_PER_MOLECULE], c[present], side="right")
    n_cells = len(EDGES[READS_PER_MOLECULE]) + 1
    full = np.bincount(idx, weights=h[present], minlength=n_cells).astype(np.int64)
    w = np.bincount(idx, weights=(h[present] * c[present]).astype(np.float64), minlength=n_cells).astype(np.int64)
    return hist_from_full(READS_PER_MOLECULE, full, w, w_name="reads", vmin=int(c[present][0]), vmax=int(c[present][-1]))


def saturation_from_hist(h: dict, fraction: float = 1.0) -> float | None:
    """Saturation at a share ``fraction`` of the reads: ``1 - sum h[c] (1 - (1 - p)^c) / (p sum c h[c])``, the
    exact expectation under random subsampling. A wide bin stands for its mean number of reads per molecule."""
    n = np.asarray(h.get("n") or [], dtype=np.float64)
    w = np.asarray(h.get("w") or [], dtype=np.float64)
    if len(n) == 0 or w.sum() <= 0:
        return None
    ok = n > 0
    c = w[ok] / n[ok]
    p = float(fraction)
    unique = float((n[ok] * (1.0 - np.power(1.0 - p, c))).sum()) if p < 1 else float(n[ok].sum())
    return 1.0 - unique / (p * float(w.sum()))


def rank_curve(counts, *, points: int = 200) -> tuple[list[int], list[float]]:
    """Counts sorted from the largest down, kept at up to ``points`` log-spaced ranks (rank 1 and the last rank
    always kept): the barcode rank plot."""
    v = np.sort(np.asarray(counts))[::-1]
    v = v[v > 0]
    n = len(v)
    if n == 0:
        return [], []
    ranks = np.unique(np.concatenate([[1], np.round(np.logspace(0, np.log10(n), points)).astype(np.int64), [n]]))
    ranks = ranks[(ranks >= 1) & (ranks <= n)]
    return [int(r) for r in ranks], [_num(float(x)) for x in v[ranks - 1]]


def summary_stats(values) -> dict:
    """Median, mean and the 10th and 90th percentile of one number per unit (None when there is no unit)."""
    v = np.asarray(values, dtype=np.float64)
    if len(v) == 0:
        return {"median": None, "mean": None, "p10": None, "p90": None, "n": 0}
    p10, med, p90 = np.percentile(v, [10, 50, 90])
    return {"median": float(med), "mean": float(v.mean()), "p10": float(p10), "p90": float(p90), "n": int(len(v))}


# ---------------------------------------------------------------------------------------------------------
# grid channels (the spatial map): one byte per cell, row-major, base64 of the raw bytes ("u8-b64")
# ---------------------------------------------------------------------------------------------------------
LEVELS = 255            # levels 0..254 carry values; 255 means "no data" on a channel that says nodata
NODATA = 255


def encode_u8(q) -> str:
    return base64.b64encode(np.ascontiguousarray(q, dtype=np.uint8).tobytes()).decode("ascii")


def decode_u8(text: str) -> np.ndarray:
    return np.frombuffer(base64.b64decode(text), dtype=np.uint8)


def channel_log1p(a, *, unit: str | None = None, valid=None) -> dict:
    """value = expm1(level * scale). ``valid`` marks the cells that have a value (the others say "no data")."""
    a = np.asarray(a, dtype=np.float64)
    vmax = float(a.max()) if a.size else 0.0
    scale = float(np.log1p(vmax) / (LEVELS - 1)) if vmax > 0 else 1.0
    q = np.rint(np.log1p(np.maximum(a, 0)) / scale).astype(np.int64) if vmax > 0 else np.zeros(a.shape, np.int64)
    q = np.clip(q, 0, LEVELS - 1).astype(np.uint8)
    ch = {"enc": "u8-b64", "transform": "log1p", "scale": scale, "offset": 0, "levels": LEVELS}
    if valid is not None:
        q = np.where(np.asarray(valid, dtype=bool), q, NODATA).astype(np.uint8)
        ch["nodata"] = True
    if unit:
        ch["unit"] = unit
    ch["data"] = encode_u8(q)
    return ch


def channel_linear(a, *, lo: float = 0.0, hi: float = 1.0, unit: str | None = None, valid=None) -> dict:
    """value = offset + level * scale, between ``lo`` and ``hi``."""
    a = np.asarray(a, dtype=np.float64)
    scale = (hi - lo) / (LEVELS - 1) if hi > lo else 1.0
    q = np.clip(np.rint((np.clip(a, lo, hi) - lo) / scale), 0, LEVELS - 1).astype(np.uint8)
    ch = {"enc": "u8-b64", "transform": "linear", "scale": float(scale), "offset": float(lo), "levels": LEVELS}
    if valid is not None:
        q = np.where(np.asarray(valid, dtype=bool), q, NODATA).astype(np.uint8)
        ch["nodata"] = True
    if unit:
        ch["unit"] = unit
    ch["data"] = encode_u8(q)
    return ch


def channel_values(ch: dict) -> np.ndarray:
    """The values a channel stands for (NaN where it says "no data"): what the viewer shows, for the tests."""
    q = decode_u8(ch["data"]).astype(np.float64)
    v = np.expm1(q * ch["scale"]) if ch.get("transform") == "log1p" else ch.get("offset", 0.0) + q * ch["scale"]
    if ch.get("nodata"):
        v = np.where(q == NODATA, np.nan, v)
    return v
