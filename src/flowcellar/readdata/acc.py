"""Integer accumulators of the read pass, and their stored (JSON) form.

One ``ReadAcc`` holds everything the read-level sections need about the primary records of one group (one run key
and one sample label inside one file). Every member is an integer or an integer array, so adding accumulators in any
order gives the same result: that is what makes a pass independent of the number of workers.
"""
from __future__ import annotations

import numpy as np

from flowcellar.aggregates import build
from flowcellar.model import axes

LEN_AXIS = "read_length_log40"
QUAL_AXIS = "qscore_0p25"
GRID_X = "read_length_log20"
GRID_Y = "qscore_1"
IDENT_AXIS = "identity_q_0p1"
MAPQ_AXIS = "mapq_1"

_SCALARS = ("reads", "bases", "q_missing", "q_computed", "mapped", "unmapped", "mapped_bases", "perfect", "aln_cols",
            "nm_sum", "ident_missing", "with_sa", "secondary", "supplementary", "split_children")
_ARRAYS = {"h_len": LEN_AXIS, "h_len_w": LEN_AXIS, "h_q": QUAL_AXIS, "h_id": IDENT_AXIS, "h_mapq": MAPQ_AXIS}
_GRID_SHAPE = (axes.n_cells(GRID_X), axes.n_cells(GRID_Y))


class ReadAcc:
    __slots__ = _SCALARS + tuple(_ARRAYS) + ("g_n", "g_w", "min_len", "max_len")

    def __init__(self):
        for k in _SCALARS:
            setattr(self, k, 0)
        for k, ax in _ARRAYS.items():
            setattr(self, k, np.zeros(axes.n_cells(ax), np.int64))
        self.g_n = None                     # length x quality grid, made on first use
        self.g_w = None
        self.min_len = None
        self.max_len = None

    # ------------------------------------------------------------------------------------ filling
    def add(self, length, q, mapped, cols, nm, mapq, sa) -> None:
        """A batch of primary records as parallel float64 arrays. ``q`` is NaN where a read has no quality;
        ``cols`` (aligned columns M + I + D) and ``nm`` are -1 where unknown."""
        n = len(length)
        if not n:
            return
        self.reads += n
        self.bases += int(length.sum())
        lo, hi = int(length.min()), int(length.max())
        self.min_len = lo if self.min_len is None else min(self.min_len, lo)
        self.max_len = hi if self.max_len is None else max(self.max_len, hi)
        li = axes.bin_index(length, LEN_AXIS)
        self.h_len += np.bincount(li, minlength=len(self.h_len))
        self.h_len_w += np.rint(np.bincount(li, weights=length, minlength=len(self.h_len))).astype(np.int64)
        ok = np.isfinite(q)
        n_ok = int(ok.sum())
        self.q_missing += n - n_ok
        if n_ok:
            qq, ll = q[ok], length[ok]
            self.h_q += np.bincount(axes.bin_index(qq, QUAL_AXIS), minlength=len(self.h_q))
            if self.g_n is None:
                self.g_n = np.zeros(_GRID_SHAPE, np.int64)
                self.g_w = np.zeros(_GRID_SHAPE, np.int64)
            flat = axes.bin_index(ll, GRID_X) * _GRID_SHAPE[1] + axes.bin_index(qq, GRID_Y)
            size = _GRID_SHAPE[0] * _GRID_SHAPE[1]
            self.g_n += np.bincount(flat, minlength=size).reshape(_GRID_SHAPE)
            self.g_w += np.rint(np.bincount(flat, weights=ll, minlength=size)).astype(np.int64).reshape(_GRID_SHAPE)
        m = mapped > 0
        n_m = int(m.sum())
        self.mapped += n_m
        self.unmapped += n - n_m
        if n_m:
            self.mapped_bases += int(length[m].sum())
            self.with_sa += int(sa[m].sum())
            self.h_mapq += np.bincount(axes.bin_index(mapq[m], MAPQ_AXIS), minlength=len(self.h_mapq))
            good = m & (cols > 0) & (nm >= 0)
            n_g = int(good.sum())
            self.ident_missing += n_m - n_g
            if n_g:
                c, e = cols[good], nm[good]
                self.aln_cols += int(c.sum())
                self.nm_sum += int(e.sum())
                self.perfect += int((e == 0).sum())
                with np.errstate(divide="ignore"):
                    idq = np.where(e > 0, -10.0 * np.log10(np.maximum(e, 1e-300) / c), np.inf)
                idq = np.maximum(idq, 0.0)          # more differences than columns cannot happen; guard all the same
                self.h_id += np.bincount(axes.bin_index(idq, IDENT_AXIS), minlength=len(self.h_id))

    # ------------------------------------------------------------------------------------ merging
    def merge(self, other: ReadAcc) -> None:
        for k in _SCALARS:
            setattr(self, k, getattr(self, k) + getattr(other, k))
        for k in _ARRAYS:
            getattr(self, k).__iadd__(getattr(other, k))
        if other.g_n is not None:
            if self.g_n is None:
                self.g_n, self.g_w = other.g_n.copy(), other.g_w.copy()
            else:
                self.g_n += other.g_n
                self.g_w += other.g_w
        for k, fn in (("min_len", min), ("max_len", max)):
            a, b = getattr(self, k), getattr(other, k)
            setattr(self, k, b if a is None else a if b is None else fn(a, b))

    def to_state(self) -> dict:
        """A picklable form that crosses the process boundary cheaply (sparse where it pays)."""
        d = {k: getattr(self, k) for k in _SCALARS}
        d.update({k: getattr(self, k) for k in _ARRAYS})
        d["min_len"], d["max_len"] = self.min_len, self.max_len
        if self.g_n is not None:
            idx = np.flatnonzero(self.g_n)
            d["g"] = (idx, self.g_n.ravel()[idx], self.g_w.ravel()[idx])
        return d

    @classmethod
    def from_state(cls, d: dict) -> ReadAcc:
        a = cls()
        for k in _SCALARS:
            setattr(a, k, d[k])
        for k in _ARRAYS:
            setattr(a, k, np.asarray(d[k], np.int64).copy())
        a.min_len, a.max_len = d["min_len"], d["max_len"]
        if "g" in d:
            idx, n, w = d["g"]
            a.g_n = np.zeros(_GRID_SHAPE, np.int64)
            a.g_w = np.zeros(_GRID_SHAPE, np.int64)
            a.g_n.ravel()[idx] = n
            a.g_w.ravel()[idx] = w
        return a

    # ------------------------------------------------------------------------------- stored cells
    def cells(self, *, aligned: bool) -> dict[str, dict]:
        """The cell members per module, in the bundle's aggregate JSON: ``{"len": {...}, "qual": {...}, "aln": {...}}``.
        ``aln`` is written only for a file that can hold alignments."""
        out: dict[str, dict] = {}
        out["len"] = {
            "scalars": {"reads": int(self.reads), "bases": int(self.bases), "split_children": int(self.split_children)},
            "hists": {"len_hist": build.hist(LEN_AXIS, self.h_len, self.h_len_w, w_name="bases", missing=0,
                                             vmin=self.min_len, vmax=self.max_len)},
        }
        qual: dict = {"scalars": {"reads": int(self.reads), "quality_from_base_qualities": int(self.q_computed)},
                      "hists": {"qual_hist": build.hist(QUAL_AXIS, self.h_q, missing=int(self.q_missing))}}
        if self.g_n is not None:
            inner_n, inner_w = self.g_n[1:-1, 1:-1], self.g_w[1:-1, 1:-1]
            h2 = build.hist2d(GRID_X, GRID_Y, inner_n, inner_w, outside=int(self.g_n.sum() - inner_n.sum()),
                              missing=int(self.q_missing), w_name="bases")
            h2["w_outside"] = int(self.g_w.sum() - inner_w.sum())
            qual["hists2d"] = {"qual_vs_length": h2}
        out["qual"] = qual
        if aligned:
            out["aln"] = {
                "scalars": {"reads": int(self.reads), "mapped_reads": int(self.mapped), "unmapped_reads": int(self.unmapped),
                            "mapped_bases": int(self.mapped_bases), "perfect_reads": int(self.perfect),
                            "alignment_columns": int(self.aln_cols), "edit_distance": int(self.nm_sum),
                            "reads_with_supplementary": int(self.with_sa),
                            "secondary_records": int(self.secondary), "supplementary_records": int(self.supplementary)},
                "hists": {"aln_identity_hist": build.hist(IDENT_AXIS, self.h_id,
                                                          missing=int(self.unmapped + self.ident_missing)),
                          "aln_mapq_hist": build.hist(MAPQ_AXIS, self.h_mapq, missing=int(self.unmapped))},
            }
        return out
