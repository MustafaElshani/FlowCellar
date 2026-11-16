"""Reads per transcript, and how many transcripts fewer reads would have found: both from the read count of every
transcript, without a read identifier and without drawing a single random number.

**What the saturation curve measures.** Take the N reads that mapped to a transcript and the count c of each
transcript. If only n of those N reads had been sequenced (any n of them, each choice equally likely), a transcript
with c reads would be seen at least once with probability ``1 - C(N - c, n) / C(N, n)`` (the hypergeometric
distribution: drawing without replacement). The curve is the sum of that probability over all transcripts: the
*expected* number of transcripts detected at depth n. At n = N it is exactly the number of transcripts with a read.
It is not a replay of the run in time, and it says nothing about transcripts no read of this sample mapped to.

The same for "at least k reads" with the lower tail of the same distribution.
"""
from __future__ import annotations

import numpy as np

#: shares of the reads at which the curve is worked out
FRACTIONS = (0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
_HALF_LOG_2PI = 0.9189385332046727


def lgamma(x):
    """log Gamma(x) for an array of x >= 1, to better than 1e-9: the argument is moved up by ten and Stirling's
    series is used there."""
    x = np.asarray(x, dtype=np.float64)
    shift = np.zeros_like(x)
    for i in range(10):
        shift += np.log(x + i)
    y = x + 10.0
    inv = 1.0 / y
    inv2 = inv * inv
    series = inv * (1.0 / 12.0 - inv2 * (1.0 / 360.0 - inv2 * (1.0 / 1260.0 - inv2 / 1680.0)))
    return (y - 0.5) * np.log(y) - y + _HALF_LOG_2PI + series - shift


def _log_comb(a, b):
    """log C(a, b) for arrays with 0 <= b <= a."""
    return lgamma(a + 1.0) - lgamma(b + 1.0) - lgamma(a - b + 1.0)


def expected_detected(counts, depths, at_least: int = 1) -> list[float]:
    """For every depth n in ``depths``: the expected number of items with at least ``at_least`` reads when n of the
    N = sum(counts) reads are drawn without replacement."""
    c = np.asarray(counts, dtype=np.int64)
    c = c[c > 0]
    total = int(c.sum())
    if not len(c):
        return [0.0 for _ in depths]
    vals, mult = np.unique(c, return_counts=True)
    v = vals.astype(np.float64)
    out = []
    for n in depths:
        n = int(min(max(n, 0), total))
        if n == 0:
            out.append(0.0)
            continue
        if n == total:
            out.append(float(mult[vals >= at_least].sum()))
            continue
        base = _log_comb(np.float64(total), np.float64(n))
        below = np.zeros(len(vals))                         # P(fewer than at_least reads of the item are drawn)
        for j in range(at_least):
            ok = (v >= j) & (total - v >= n - j) & (n - j >= 0)
            if not ok.any():
                continue
            vv = v[ok]
            lp = _log_comb(vv, np.full(len(vv), float(j))) + _log_comb(total - vv, np.full(len(vv), float(n - j))) - base
            below[ok] += np.exp(lp)
        out.append(float((mult * (1.0 - np.minimum(below, 1.0))).sum()))
    return out


def depths_for(total: int) -> list[int]:
    """The depths of the curve for N reads: the fixed shares of N, each at least one read, without repeats."""
    out: list[int] = []
    for f in FRACTIONS:
        n = max(1, int(round(f * total)))
        if not out or n > out[-1]:
            out.append(n)
    return out if total > 0 else []


def curve(transcripts, genes=None) -> dict | None:
    """The saturation curve in the bundle's ``curve`` form: x = reads drawn, y = expected transcripts (and genes, and
    genes with ten reads or more) detected. None when there is no read."""
    t = np.asarray(transcripts, dtype=np.int64)
    total = int(t.sum())
    if total <= 0:
        return None
    x = depths_for(total)
    y = {"transcripts": [round(v, 1) for v in expected_detected(t, x)]}
    if genes is not None:
        g = np.asarray(genes, dtype=np.int64)
        y["genes"] = [round(v, 1) for v in expected_detected(g, x)]
        y["genes_10_reads"] = [round(v, 1) for v in expected_detected(g, x, at_least=10)]
    return {"x": x, "y": y, "x_unit": "reads", "y_unit": "detected", "x_scale": "linear"}


def reads_per_item(counts) -> dict | None:
    """How the reads are spread: items (transcripts) by their read count, in classes that double (1, 2 to 3, 4 to 7,
    ...). x = the smallest count of each class; y = how many items are in it and how many reads they hold."""
    c = np.asarray(counts, dtype=np.int64)
    c = c[c > 0]
    if not len(c):
        return None
    k = np.frexp(c.astype(np.float64))[1] - 1               # floor(log2 c), exact for whole numbers
    n = np.bincount(k)
    w = np.bincount(k, weights=c).astype(np.int64)
    return {"x": [int(2 ** i) for i in range(len(n))], "y": {"transcripts": [int(v) for v in n], "reads": [int(v) for v in w]},
            "x_unit": "reads per transcript, from", "y_unit": "count", "x_scale": "log10"}
