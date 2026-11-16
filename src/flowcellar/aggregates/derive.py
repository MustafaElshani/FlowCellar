"""The fixed set of derivations: numbers read from stored aggregates.

Quantiles and N50 are read from histograms by rank interpolation inside a bin (geometric on log axes); they are
exact to bin resolution. The viewer implements the same functions; ``tests/vectors/derive`` holds them equal.
"""
from __future__ import annotations

import numpy as np

from flowcellar.aggregates.build import hist_full
from flowcellar.model.axes import AXES, EDGES


def total(h: dict, key: str = "n") -> int:
    """All values of a histogram, under- and overflow included (not ``missing``)."""
    if key == "n":
        return int(sum(h.get("n", [])) + h.get("under", 0) + h.get("over", 0))
    return int(sum(h.get("w", [])) + h.get("w_under", 0) + h.get("w_over", 0))


def ratio(num, den):
    return None if num is None or not den else num / den


def mean(h: dict) -> float | None:
    """Exact mean of the binned quantity: sum of weights over count (needs ``w``)."""
    n = total(h, "n")
    return None if "w" not in h or not n else total(h, "w") / n


def _interp(axis: str, b: int, frac: float) -> float:
    e = EDGES[axis]
    if b <= 0:
        return float(e[0])
    if b >= len(e):
        return float(e[-1])
    lo, hi = e[b - 1], e[b]
    if AXES[axis]["type"] == "log10":
        return float(lo * (hi / lo) ** frac)
    return float(lo + (hi - lo) * frac)


def quantile(h: dict, q: float) -> float | None:
    """The q-quantile of the binned values (0 < q < 1); None for an empty histogram."""
    c = hist_full(h, "n").astype(np.float64)
    tot = c.sum()
    if tot <= 0:
        return None
    cum = np.cumsum(c)
    b = int(np.searchsorted(cum, q * tot, side="left"))
    frac = (q * tot - (cum[b - 1] if b else 0.0)) / c[b] if c[b] else 0.0
    return _interp(h["axis"], b, frac)


def median(h: dict) -> float | None:
    return quantile(h, 0.5)


def nx(h: dict, x: float = 0.5) -> float | None:
    """N50-type statistic from the per-bin sums ``w``: the value at which the bins above it hold the share x of the
    weighted total (N50: x = 0.5)."""
    if "w" not in h:
        return None
    w = hist_full(h, "w").astype(np.float64)
    tot = w.sum()
    if tot <= 0:
        return None
    cum = np.cumsum(w)
    e = EDGES[h["axis"]]
    b = int(np.searchsorted(cum, (1.0 - x) * tot, side="left"))
    b = min(max(b, 1), len(e) - 1)
    frac = ((1.0 - x) * tot - cum[b - 1]) / w[b] if w[b] > 0 else 0.0
    return _interp(h["axis"], b, frac)


def frac_above(h: dict, x: float, key: str = "n") -> float | None:
    """Share of the values (or of the weight) at or above ``x``; ``x`` is snapped down to the edge of its bin."""
    c = hist_full(h, key).astype(np.float64)
    tot = c.sum()
    if tot <= 0:
        return None
    b = int(np.searchsorted(EDGES[h["axis"]], x, side="right"))
    return float(c[b:].sum() / tot)


def mode(h: dict) -> float | None:
    """Centre of the fullest bin (geometric on log axes)."""
    n = h.get("n") or []
    if not n or max(n) <= 0:
        return None
    return _interp(h["axis"], 1 + h["start"] + int(np.argmax(n)), 0.5)


def cumulative(values) -> list:
    return [int(v) if float(v).is_integer() else float(v) for v in np.cumsum(np.asarray(values, dtype=np.float64))]


def rate(values, step_s: float) -> list[float]:
    """Per-second rate of per-bin increments."""
    return [float(v) / step_s for v in values]


def gini(x) -> float | None:
    """Gini coefficient of per-channel output (0 = even)."""
    x = np.sort(np.asarray(x, dtype=np.float64))
    n, s = len(x), x.sum()
    if n == 0 or s <= 0:
        return None
    return float((2 * np.sum(np.arange(1, n + 1) * x) / (n * s)) - (n + 1) / n)
