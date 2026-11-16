"""From numpy accumulators to the stored (JSON) form of each aggregate kind.

Internal histogram arrays are ``(under, bin 0 ... bin n-1, over)``: ``numpy.bincount`` of ``axes.bin_index``.
"""
from __future__ import annotations

import math

import numpy as np

from flowcellar.model.axes import EDGES

B64 = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/"
B64URL = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-_"


def hist(axis: str, full, w_full=None, *, w_name: str | None = None, missing: int | None = None,
         vmin=None, vmax=None) -> dict:
    """Internal array -> stored histogram with empty ends trimmed. ``w_full`` is the per-cell sum of the binned
    quantity (bases for length); its under- and overflow cells are kept as ``w_under`` and ``w_over``."""
    full = np.asarray(full)
    if len(full) != len(EDGES[axis]) + 1:
        raise ValueError(f"{axis}: expected {len(EDGES[axis]) + 1} cells, got {len(full)}")
    inner = full[1:-1]
    nz = np.nonzero(inner)[0]
    a, b = (int(nz[0]), int(nz[-1]) + 1) if len(nz) else (0, 0)
    h: dict = {"axis": axis, "start": a, "n": [int(v) for v in inner[a:b]]}
    if w_full is not None:
        w_full = np.asarray(w_full)
        h["w"] = [int(v) for v in w_full[1:-1][a:b]]
        if w_name:
            h["w_name"] = w_name
        h["w_under"] = int(w_full[0])
        h["w_over"] = int(w_full[-1])
    h["under"] = int(full[0])
    h["over"] = int(full[-1])
    if missing is not None:
        h["missing"] = int(missing)
    if vmin is not None:
        h["min"] = vmin
    if vmax is not None:
        h["max"] = vmax
    return h


def hist_full(h: dict, key: str = "n") -> np.ndarray:
    """Stored histogram -> internal array (under, bins..., over) of ``n`` or ``w``."""
    out = np.zeros(len(EDGES[h["axis"]]) + 1, np.int64)
    vals = h.get(key) or []
    out[1 + h["start"]: 1 + h["start"] + len(vals)] = vals
    if key == "n":
        out[0], out[-1] = h.get("under", 0), h.get("over", 0)
    else:
        out[0], out[-1] = h.get("w_under", 0), h.get("w_over", 0)
    return out


def hist2d(x: str, y: str, matrix, w_matrix=None, *, outside: int = 0, missing: int | None = None,
           w_name: str | None = None) -> dict:
    """Matrix of the INNER bins (rows = x bins, columns = y bins) -> stored sparse 2-D histogram."""
    m = np.asarray(matrix)
    i, j = np.nonzero(m)
    h: dict = {"x": x, "y": y, "enc": "coo", "i": [int(v) for v in i], "j": [int(v) for v in j],
               "n": [int(v) for v in m[i, j]]}
    if w_matrix is not None:
        h["w"] = [int(v) for v in np.asarray(w_matrix)[i, j]]
        if w_name:
            h["w_name"] = w_name
    h["outside"] = int(outside)
    if missing is not None:
        h["missing"] = int(missing)
    return h


def series(axis: str, values: dict, *, unit: str | None = None, levels: dict[str, str] | None = None,
           n_bins: int | None = None, start: int = 0) -> dict:
    """Named arrays on one axis. Increments are cut after the last bin in which any of them is non-zero (never
    before ``n_bins``); ``levels`` maps the names that are levels to their merge rule and may hold None values."""
    levels = levels or {}
    inc = {k: np.asarray(v) for k, v in values.items() if k not in levels}
    last = 0
    for v in inc.values():
        nz = np.nonzero(v)[0]
        if len(nz):
            last = max(last, int(nz[-1]) + 1)
    for k in levels:
        last = max(last, len(values[k]))
    if n_bins is not None:
        last = max(last, n_bins)
    out_vals: dict = {}
    for k, v in values.items():
        if k in levels:
            vv = [None if x is None else (x if isinstance(x, (int, float)) else float(x)) for x in v]
            out_vals[k] = vv + [None] * (last - len(vv))
        else:
            vv = [int(x) for x in np.asarray(v)[:last]]
            out_vals[k] = vv + [0] * (last - len(vv))
    s: dict = {"axis": axis, "start": start}
    if unit:
        s["unit"] = unit
    s["values"] = out_vals
    if levels:
        s["levels"] = sorted(levels)
        s["merge"] = dict(levels)
    return s


def counts(cats: list[str], n: list[int], *, w: list[int] | None = None, role: str | None = None,
           missing: list[str] | None = None, sources: dict | None = None) -> dict:
    c: dict = {"cats": list(cats), "n": [int(v) for v in n]}
    if w is not None:
        c["w"] = [int(v) for v in w]
    if role:
        c["role"] = role
    if missing:
        c["missing"] = list(missing)
    if sources:
        c["sources"] = sources
    return c


def largest_remainder(vals, total: int) -> list[int]:
    """Integers that sum to ``total``, each within 1 of its exact value."""
    fl = [int(math.floor(v)) for v in vals]
    rem = int(total - sum(fl))
    order = sorted(range(len(vals)), key=lambda i: vals[i] - fl[i], reverse=True)
    for i in order[: max(rem, 0)]:
        fl[i] += 1
    return fl


def u6_encode(levels, alphabet: str = B64) -> str:
    """One character per cell (levels 0..63), row-major: grids ``enc: u6`` and the index sparklines (B64URL)."""
    lut = np.frombuffer(alphabet.encode(), dtype=np.uint8)
    return lut[np.clip(np.asarray(levels), 0, 63).astype(np.int64).ravel()].tobytes().decode()


def u6_decode(text: str, alphabet: str = B64) -> np.ndarray:
    lut = np.full(256, 0, np.uint8)
    lut[np.frombuffer(alphabet.encode(), dtype=np.uint8)] = np.arange(64, dtype=np.uint8)
    return lut[np.frombuffer(text.encode(), dtype=np.uint8)]


def sparkline(values, *, points: int = 32) -> tuple[str, float]:
    """(32 characters of the base64url alphabet, the row's maximum): 6 bits per point, scaled to the row's maximum."""
    v = np.zeros(points)
    vals = np.asarray(values, dtype=np.float64)[:points]
    v[: len(vals)] = vals
    mx = float(v.max()) if len(v) else 0.0
    lv = np.rint(v / mx * 63).astype(np.int64) if mx > 0 else np.zeros(points, np.int64)
    return u6_encode(lv, B64URL), mx
