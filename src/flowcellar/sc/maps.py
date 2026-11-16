"""The spatial map of an analysis: a picture of the capture area at a stored resolution.

A map is a grid of channels, one byte per cell and channel (``flowcellar.sc.axes``). It is a picture: exact
totals, medians and under-tissue shares are computed from the full-resolution data and stored as numbers of
their own, never read off the map.

Visium HD      the 2 um grid of ``feature_slice.h5`` summed into square bins (64 um by default: 105 x 105 cells)
Visium         every spot of the slide, no binning; cells are the rows and columns of the staggered spot grid
ONT workflow   Visium HD without a Space Ranger folder: molecules, genes and reads per bin from the 8 um matrix
"""
from __future__ import annotations

import numpy as np

from flowcellar.sc import axes
from flowcellar.sc.spaceranger import rebin_sum

MAP_UM_CHOICES = (64, 32)
DEFAULT_MAP_UM = 64


def _cells_per_bin(n_fine: int, factor: int) -> np.ndarray:
    """How many fine rows (or columns) fall into each coarse one: the last one may hold fewer."""
    n = -(-n_fine // factor)
    out = np.full(n, factor, dtype=np.float64)
    if n_fine % factor:
        out[-1] = n_fine % factor
    return out


def hd_map(fs, *, map_um: int = DEFAULT_MAP_UM, genes8: tuple | None = None, tissue_image: bool = True,
           mito_prefix: str = "MT-") -> tuple[dict, dict]:
    """(grid, exact numbers) for a Visium HD output. ``fs`` is an open ``spaceranger.FeatureSlice``; ``genes8`` is
    (row, col, genes) of the 8 um bins under tissue, from the pointers of the 8 um matrix.

    The exact numbers are computed here because the layers are in memory anyway: molecules, reads and squares in
    all and under tissue, mitochondrial molecules, and per 8 um bin the molecules, reads and mitochondrial share.
    """
    pitch = fs.pitch_um or 2.0
    f = max(1, int(round(map_um / pitch)))
    nr, nc = fs.nrows, fs.ncols
    ur, uc, uv = fs.layer("umis/total")
    rr, rc, rv = fs.layer("reads/sequenced") if fs.has("reads/sequenced") else (ur[:0], uc[:0], uv[:0])
    mr, mc, mv = fs.layer("masks/filtered") if fs.has("masks/filtered") else (ur[:0], uc[:0], uv[:0])
    tissue2 = np.zeros((nr, nc), dtype=bool)
    tissue2[mr, mc] = mv > 0
    U = rebin_sum(ur, uc, uv, nr, nc, f)
    R = rebin_sum(rr, rc, rv, nr, nc, f)
    T = rebin_sum(mr, mc, mv > 0, nr, nc, f) / np.outer(_cells_per_bin(nr, f), _cells_per_bin(nc, f))
    tr_, tc_, tv_, n_mt = fs.mito_layers(mito_prefix)
    M = rebin_sum(tr_, tc_, tv_, nr, nc, f)
    channels: dict = {
        "umis": axes.channel_log1p(U, unit="umis"),
        "saturation": axes.channel_linear(1.0 - np.divide(U, R, out=np.zeros_like(U), where=R > 0), unit="fraction", valid=R > 0),
        "tissue": axes.channel_linear(T, unit="fraction"),
        "mito": _sqrt_channel(np.divide(M, U, out=np.zeros_like(U), where=U > 0), valid=U > 0),
    }
    f8 = max(1, int(round(8 / pitch)))
    exact: dict = {
        "umis_total": int(uv.sum(dtype=np.int64)), "reads_total": int(rv.sum(dtype=np.int64)),
        "squares_under_tissue": int((mv > 0).sum()),
        "umis_under_tissue": int(uv[tissue2[ur, uc]].sum(dtype=np.int64)),
        "reads_under_tissue": int(rv[tissue2[rr, rc]].sum(dtype=np.int64)) if len(rv) else None,
        "mito_umis": int(tv_.sum(dtype=np.int64)), "mito_features": int(n_mt),
        "grid_2um": [nr, nc], "pitch_um": pitch,
    }
    # per 8 um bin: molecules, reads, mitochondrial molecules, and which bins lie under tissue
    U8 = rebin_sum(ur, uc, uv, nr, nc, f8)
    R8 = rebin_sum(rr, rc, rv, nr, nc, f8)
    M8 = rebin_sum(tr_, tc_, tv_, nr, nc, f8)
    if fs.has("masks/square_008um"):
        r8, c8, v8 = fs.layer("masks/square_008um")
        T8 = np.zeros(U8.shape, dtype=bool)
        T8[r8, c8] = v8 > 0
    else:
        T8 = rebin_sum(mr, mc, mv > 0, nr, nc, f8) > 0
    exact.update(umis8=U8, reads8=R8, mito8=M8, tissue8=T8)
    if genes8 is not None:
        g_r, g_c, g_n = genes8
        k = max(1, f // f8)
        shape = U.shape
        flat = (g_r // k) * shape[1] + (g_c // k)
        ok = flat < shape[0] * shape[1]
        g_sum = np.bincount(flat[ok], weights=g_n[ok], minlength=shape[0] * shape[1]).reshape(shape)
        g_cnt = np.bincount(flat[ok], minlength=shape[0] * shape[1]).reshape(shape)
        channels["genes"] = axes.channel_log1p(np.divide(g_sum, g_cnt, out=np.zeros_like(g_sum), where=g_cnt > 0),
                                               unit="genes", valid=g_cnt > 0)
    if tissue_image:
        name = "images/microscope" if fs.has("images/microscope") else "images/cytassist" if fs.has("images/cytassist") else None
        if name:
            ir, ic, iv = fs.layer(name)
            s = rebin_sum(ir, ic, iv, nr, nc, f)
            n = rebin_sum(ir, ic, np.ones(len(iv)), nr, nc, f)
            img = np.divide(s, n, out=np.zeros_like(s), where=n > 0)
            channels["image"] = {"enc": "u8-b64", "transform": "linear", "scale": 1.0, "offset": 0.0, "levels": 256,
                                 "unit": "grey", "data": axes.encode_u8(np.clip(np.rint(img), 0, 255))}
    grid = {"rows": int(U.shape[0]), "cols": int(U.shape[1]), "layout": {"kind": "square", "bin_um": float(f * pitch)},
            "channels": channels}
    return grid, exact


def _sqrt_channel(a, *, valid=None) -> dict:
    """A share between 0 and 1 with fine steps near zero: value = (level * scale) squared."""
    a = np.clip(np.asarray(a, dtype=np.float64), 0.0, 1.0)
    scale = 1.0 / (axes.LEVELS - 1)
    q = np.clip(np.rint(np.sqrt(a) / scale), 0, axes.LEVELS - 1).astype(np.uint8)
    ch = {"enc": "u8-b64", "transform": "sqrt", "scale": scale, "offset": 0.0, "levels": axes.LEVELS, "unit": "fraction"}
    if valid is not None:
        q = np.where(np.asarray(valid, dtype=bool), q, axes.NODATA).astype(np.uint8)
        ch["nodata"] = True
    ch["data"] = axes.encode_u8(q)
    return ch


def sqrt_values(ch: dict) -> np.ndarray:
    q = axes.decode_u8(ch["data"]).astype(np.float64)
    v = np.square(q * ch["scale"]) + ch.get("offset", 0.0)
    return np.where(q == axes.NODATA, np.nan, v) if ch.get("nodata") else v


def visium_grid_shape(rows, cols) -> tuple[int, int]:
    """(rows, columns) of the stored grid for spots at (array_row, array_col): spots sit on every second column,
    shifted by one on odd rows, so the stored column is ``array_col // 2``."""
    return int(np.max(rows)) + 1, (int(np.max(cols)) + 2) // 2


def visium_map(rows, cols, in_tissue, *, umis=None, genes=None, mito=None, saturation=None, sat_valid=None,
               reads=None, pitch_um: float = 100.0) -> dict:
    """One value per spot of the slide, in (array_row, array_col // 2) cells. Every array has one entry per spot
    of the position list; a spot that has no data holds zero."""
    rows, cols = np.asarray(rows, dtype=np.int64), np.asarray(cols, dtype=np.int64)
    nr, nc = visium_grid_shape(rows, cols)
    flat = rows * nc + cols // 2

    def cells(values, fill=0.0):
        a = np.full(nr * nc, fill, dtype=np.float64)
        a[flat] = values
        return a.reshape(nr, nc)

    present = cells(np.ones(len(rows)), 0.0) > 0
    channels: dict = {"tissue": axes.channel_linear(cells(np.asarray(in_tissue, dtype=np.float64)), unit="fraction", valid=present)}
    if umis is not None:
        channels["umis"] = axes.channel_log1p(cells(umis), unit="umis", valid=present)
    if genes is not None:
        channels["genes"] = axes.channel_log1p(cells(genes), unit="genes", valid=present)
    if reads is not None:
        channels["reads"] = axes.channel_log1p(cells(reads), unit="reads", valid=present)
    if mito is not None and umis is not None:
        channels["mito"] = _sqrt_channel(cells(mito), valid=present & (cells(umis) > 0))
    if saturation is not None:
        ok = present & (cells(np.asarray(sat_valid, dtype=np.float64)) > 0) if sat_valid is not None else present
        channels["saturation"] = axes.channel_linear(cells(saturation), unit="fraction", valid=ok)
    return {"rows": nr, "cols": nc, "layout": {"kind": "hex", "bin_um": float(pitch_um)}, "channels": channels}


def bins_map(rows8, cols8, *, umis, genes, reads=None, bin_um: int = 8, map_um: int = DEFAULT_MAP_UM,
             grid_bins: int | None = None) -> dict:
    """Visium HD from the workflow's 8 um matrix alone: molecules (sum), genes (mean per 8 um bin with data) and
    reads (sum) on square bins of ``map_um``. No tissue mask and no image: the workflow does not know them."""
    k = max(1, map_um // bin_um)
    rows8, cols8 = np.asarray(rows8, dtype=np.int64), np.asarray(cols8, dtype=np.int64)
    n8 = max(int(rows8.max()) + 1 if len(rows8) else 1, int(cols8.max()) + 1 if len(cols8) else 1, grid_bins or 0)
    n = -(-n8 // k)
    flat = (rows8 // k) * n + (cols8 // k)
    U = np.bincount(flat, weights=np.asarray(umis, dtype=np.float64), minlength=n * n).reshape(n, n)
    G = np.bincount(flat, weights=np.asarray(genes, dtype=np.float64), minlength=n * n).reshape(n, n)
    C = np.bincount(flat, minlength=n * n).reshape(n, n)
    channels = {"umis": axes.channel_log1p(U, unit="umis"),
                "genes": axes.channel_log1p(np.divide(G, C, out=np.zeros_like(G), where=C > 0), unit="genes", valid=C > 0)}
    if reads is not None:
        r_r, r_c, r_n = reads
        fl = (np.asarray(r_r) // k) * n + (np.asarray(r_c) // k)
        ok = fl < n * n
        channels["reads"] = axes.channel_log1p(np.bincount(fl[ok], weights=np.asarray(r_n, dtype=np.float64)[ok], minlength=n * n).reshape(n, n),
                                               unit="reads")
    return {"rows": n, "cols": n, "layout": {"kind": "square", "bin_um": float(k * bin_um)}, "channels": channels}


def grid_bytes(grid: dict) -> int:
    """The JSON size of a grid: what it costs inside the bundle."""
    return sum(len(ch.get("data", "")) for ch in grid.get("channels", {}).values()) + 200 * len(grid.get("channels", {}))
