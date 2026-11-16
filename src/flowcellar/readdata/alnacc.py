"""Integer accumulators of the read pass, second part: what the alignment, over-time and poly(A) charts need.

``ReadAcc`` (acc.py) keeps what the first build counted. Here are, per group (one run key and one sample label in
one file):

* ``AlignAcc``: from the primary mapped records: strand, the CIGAR totals, mismatches from NM, aligned length,
  soft-clipped share per read, three 2-D grids, reads and bases per contig, aligned bases per window along a genome,
  per-transcript read counts and the coverage along the transcript body for a transcriptome, and aligned length and
  identity by clock hour;
* ``TagAcc``: from the tags the basecaller wrote on every read: reads and bases per ten minutes of the clock (``st``),
  read length, quality and measured speed by clock hour (``st``, ``qs``, ``du``), the speed histogram, and the
  poly(A) tail length (``pt``).

Every member is an integer or an integer array and every operation is an addition, so accumulators can be added in
any order: the result does not depend on how a file was cut into tasks or on the number of workers.

**Time is binned on the clock**, not on run time: a pass is a pure function of its file and does not know when the
run started. A ten-minute bin is ``floor(epoch seconds / 600)``, an hour bin ``floor(epoch seconds / 3600)``; the
sections turn them into time since the start of the run when the bundle is put together.

A group hands its records over as a batch of float64, one row per column of ``COLS`` and one column per primary
record (so that every quantity is one contiguous array).
"""
from __future__ import annotations

import numpy as np

from flowcellar.aggregates import build
from flowcellar.model import axes
from flowcellar.readdata.reference import BODY_BINS, OTHER_CONTIGS, TX_CLASSES, Layout

COLS = ("length", "q", "mapped", "m", "ins", "del", "skip", "soft", "hard", "nm", "mapq", "sa", "rev", "rid", "start",
        "t10", "du", "pt", "walked")
NCOL = len(COLS)
(C_LEN, C_Q, C_MAPPED, C_M, C_I, C_D, C_N, C_S, C_H, C_NM, C_MAPQ, C_SA, C_REV, C_RID, C_START, C_T10, C_DU, C_PT,
 C_WALKED) = range(NCOL)

LEN_AXIS = "read_length_log40"
GRID_LEN = "read_length_log20"
GRID_Q = "qscore_1"
GRID_ID = "identity_q_0p5"
CLIP_AXIS = "fraction_0p01"
SPEED_AXIS = "speed_5bps"
POLYA_AXIS = "polya_1nt"
MAX_TIME_ROWS = 24 * 366            # a group whose reads span more than a year of clock hours keeps no time tables
MIN_DURATION_S = 1.0                # speed is measured on reads that spent at least this long in the pore
PT_NO_TAG = -2.0

_NL = axes.n_cells(GRID_LEN)
_NQ = axes.n_cells(GRID_Q)
_NI = axes.n_cells(GRID_ID)
_NS = axes.n_cells(SPEED_AXIS)


# ------------------------------------------------------------------------------------------ time tables
class TimeTable:
    """Counts in ``width`` columns per time bin, over a range of bins that grows as reads arrive."""
    __slots__ = ("width", "lo", "a", "broken")

    def __init__(self, width: int):
        self.width = width
        self.lo: int | None = None
        self.a = np.zeros((0, width), np.int64)
        self.broken = False

    def _cover(self, lo: int, hi: int) -> bool:
        """Make the rows cover bins lo..hi. False (and broken) when the range is longer than any run."""
        if self.lo is None:
            new_lo, new_hi = lo, hi
        else:
            new_lo, new_hi = min(lo, self.lo), max(hi, self.lo + len(self.a) - 1)
        if new_hi - new_lo + 1 > MAX_TIME_ROWS * (6 if self.width <= 2 else 1):
            self.broken = True
            self.a = np.zeros((0, self.width), np.int64)
            self.lo = None
            return False
        if self.lo is None or new_lo != self.lo or new_hi - new_lo + 1 != len(self.a):
            a = np.zeros((new_hi - new_lo + 1, self.width), np.int64)
            if self.lo is not None and len(self.a):
                a[self.lo - new_lo: self.lo - new_lo + len(self.a)] = self.a
            self.a, self.lo = a, new_lo
        return True

    def add(self, t, v, w=None) -> None:
        """One count (or the weight ``w``) for every pair of time bin ``t`` and column ``v`` (integer arrays)."""
        if self.broken or not len(t):
            return
        t = t.astype(np.int64)
        if not self._cover(int(t.min()), int(t.max())):
            return
        flat = (t - self.lo) * self.width + v
        size = self.a.size
        if w is None:
            self.a += np.bincount(flat, minlength=size).reshape(self.a.shape)
        else:
            self.a += np.rint(np.bincount(flat, weights=w, minlength=size)).astype(np.int64).reshape(self.a.shape)

    def merge(self, other: TimeTable) -> None:
        if other.broken:
            self.broken = True
        if self.broken:
            self.a, self.lo = np.zeros((0, self.width), np.int64), None
            return
        if other.lo is None or not len(other.a):
            return
        if not self._cover(other.lo, other.lo + len(other.a) - 1):
            return
        self.a[other.lo - self.lo: other.lo - self.lo + len(other.a)] += other.a

    def trimmed(self) -> tuple[int | None, np.ndarray]:
        """(first bin, rows) without the empty rows at both ends."""
        if self.lo is None or not self.a.any():
            return None, np.zeros((0, self.width), np.int64)
        rows = np.flatnonzero(self.a.any(axis=1))
        return self.lo + int(rows[0]), self.a[rows[0]: rows[-1] + 1]

    def to_state(self):
        if self.broken:
            return "broken"
        lo, a = self.trimmed()
        if lo is None:
            return None
        idx = np.flatnonzero(a)
        return lo, len(a), idx, a.ravel()[idx]

    @classmethod
    def from_state(cls, width: int, state) -> TimeTable:
        t = cls(width)
        if state == "broken":
            t.broken = True
        elif state is not None:
            lo, rows, idx, val = state
            t.lo = int(lo)
            t.a = np.zeros((int(rows), width), np.int64)
            t.a.ravel()[idx] = val
        return t

    def to_json(self) -> dict | None:
        """The stored form (in the ``facts`` of a fragment): sparse, sorted by (row, column)."""
        if self.broken:
            return {"broken": True}
        lo, a = self.trimmed()
        if lo is None:
            return None
        i, j = np.nonzero(a)
        return {"lo": int(lo), "rows": int(len(a)), "i": [int(x) for x in i], "j": [int(x) for x in j], "n": [int(x) for x in a[i, j]]}

    @classmethod
    def from_json(cls, width: int, d: dict | None) -> TimeTable:
        t = cls(width)
        if not d:
            return t
        if d.get("broken"):
            t.broken = True
            return t
        t.lo = int(d["lo"])
        t.a = np.zeros((int(d["rows"]), width), np.int64)
        t.a[np.asarray(d["i"], np.int64), np.asarray(d["j"], np.int64)] = np.asarray(d["n"], np.int64)
        return t


def _grid_state(g):
    if g is None:
        return None
    idx = np.flatnonzero(g)
    return idx, g.ravel()[idx]


def _grid_from(shape, state):
    if state is None:
        return None
    g = np.zeros(shape, np.int64)
    g.ravel()[state[0]] = state[1]
    return g


def _add_opt(a, b):
    if b is None:
        return a
    return b.copy() if a is None else a + b


def _hist2d(x: str, y: str, full, *, missing: int | None = None) -> dict:
    """Internal grid with under- and overflow cells on both axes -> stored 2-D histogram."""
    inner = full[1:-1, 1:-1]
    return build.hist2d(x, y, inner, outside=int(full.sum() - inner.sum()), missing=missing)


def identity_q(nm, cols):
    """Identity on the Phred scale, -10 log10(NM / columns); no difference at all is +inf (the overflow cell)."""
    with np.errstate(divide="ignore", invalid="ignore"):
        idq = np.where(nm > 0, -10.0 * np.log10(np.maximum(nm, 1e-300) / np.maximum(cols, 1.0)), np.inf)
    return np.maximum(idq, 0.0)


# --------------------------------------------------------------------------------------------- AlignAcc
_A_SCALARS = ("fwd", "rev", "m_bases", "ins_bases", "del_bases", "skip_bases", "soft_bases", "hard_bases", "mism_bases",
              "m_bases_nm", "nm_negative", "grid_q_missing", "time_missing", "tx_reads")


class AlignAcc:
    """Primary mapped records of one group."""
    __slots__ = _A_SCALARS + ("h_alen", "h_alen_w", "h_clip", "g_len_aln", "g_aln_id", "g_q_id", "slot_n", "slot_w", "cov",
                              "ref_n", "body", "tx_class_n", "tt_alen", "tt_id", "ref_counts_dropped")

    def __init__(self):
        for k in _A_SCALARS:
            setattr(self, k, 0)
        self.h_alen = np.zeros(axes.n_cells(LEN_AXIS), np.int64)
        self.h_alen_w = np.zeros(axes.n_cells(LEN_AXIS), np.int64)
        self.h_clip = np.zeros(axes.n_cells(CLIP_AXIS), np.int64)
        self.g_len_aln = self.g_aln_id = self.g_q_id = None
        self.slot_n = self.slot_w = self.cov = None             # genome
        self.ref_n = self.body = self.tx_class_n = None         # transcriptome
        self.tt_alen = TimeTable(_NL)
        self.tt_id = TimeTable(_NI)
        self.ref_counts_dropped = False

    # ------------------------------------------------------------------------------------ filling
    def ensure_cov(self, layout: Layout) -> np.ndarray:
        if self.cov is None:
            self.cov = np.zeros(layout.n_windows + 1, np.int64)             # the last cell: contigs without windows
            self.slot_n = np.zeros(len(layout.kept) + 1, np.int64)
            self.slot_w = np.zeros(len(layout.kept) + 1, np.int64)
        return self.cov

    def add(self, R: np.ndarray, layout: Layout | None, *, ref_counts: bool = True) -> None:
        """A batch of primary records (``COLS`` by records). Unmapped records are ignored."""
        R = R[:, R[C_MAPPED] > 0]
        n = R.shape[1]
        if not n:
            return
        length, q = R[C_LEN], R[C_Q]
        M, INS, D, N, S, H, nm = R[C_M], R[C_I], R[C_D], R[C_N], R[C_S], R[C_H], R[C_NM]
        n_rev = int((R[C_REV] > 0).sum())
        self.rev += n_rev
        self.fwd += n - n_rev
        self.m_bases += int(M.sum())
        self.ins_bases += int(INS.sum())
        self.del_bases += int(D.sum())
        self.skip_bases += int(N.sum())
        self.soft_bases += int(S.sum())
        self.hard_bases += int(H.sum())
        known = nm >= 0
        mism = nm - INS - D
        self.nm_negative += int((known & (mism < 0)).sum())
        self.mism_bases += int(np.where(known, np.maximum(mism, 0.0), 0.0).sum())
        self.m_bases_nm += int(M[known].sum())
        alen = M + INS                                             # read bases inside the alignment
        ai = axes.bin_index(alen, LEN_AXIS)
        self.h_alen += np.bincount(ai, minlength=len(self.h_alen))
        self.h_alen_w += np.rint(np.bincount(ai, weights=alen, minlength=len(self.h_alen))).astype(np.int64)
        with np.errstate(divide="ignore", invalid="ignore"):
            clip = np.where(length > 0, S / np.maximum(length, 1.0), 0.0)
        self.h_clip += np.bincount(axes.bin_index(clip, CLIP_AXIS), minlength=len(self.h_clip))
        # ---- the 2-D grids
        if self.g_len_aln is None:
            self.g_len_aln = np.zeros((_NL, _NL), np.int64)
            self.g_aln_id = np.zeros((_NL, _NI), np.int64)
            self.g_q_id = np.zeros((_NQ, _NI), np.int64)
        bl = (ai + 1) // 2                                       # the 20-per-decade bin: its edges are every other edge
        self.g_len_aln += np.bincount(axes.bin_index(length, GRID_LEN) * _NL + bl, minlength=_NL * _NL).reshape(_NL, _NL)
        cols = M + INS + D
        good = known & (cols > 0)
        idq = identity_q(nm[good], cols[good])
        bi = axes.bin_index(idq, GRID_ID)
        self.g_aln_id += np.bincount(bl[good] * _NI + bi, minlength=_NL * _NI).reshape(_NL, _NI)
        qg = q[good]
        qok = np.isfinite(qg)
        self.grid_q_missing += int(len(qg) - qok.sum())
        self.g_q_id += np.bincount(axes.bin_index(qg[qok], GRID_Q) * _NI + bi[qok], minlength=_NQ * _NI).reshape(_NQ, _NI)
        # ---- by clock hour
        t10 = R[C_T10]
        timed = t10 >= 0
        self.time_missing += int(n - timed.sum())
        if timed.any():
            hour = np.floor_divide(t10, 6.0).astype(np.int64)
            self.tt_alen.add(hour[timed], bl[timed])
            tg = timed[good]
            self.tt_id.add(hour[good][tg], bi[tg])
        # ---- along the reference
        if layout is None or layout.kind is None:
            return
        rid = R[C_RID].astype(np.int64)
        start = R[C_START].astype(np.int64)
        if layout.kind == "genome":
            self._genome(layout, rid, start, M.astype(np.int64), D.astype(np.int64), N.astype(np.int64), R[C_WALKED] > 0)
        else:
            self._transcriptome(layout, rid, start, (M + D + N).astype(np.int64), ref_counts)

    def _genome(self, lay: Layout, rid, start, M, D, N, walked) -> None:
        cov = self.ensure_cov(lay)
        md = M + D                                               # reference positions under an aligned column
        slot = lay.slot[rid]
        self.slot_n += np.bincount(slot, minlength=len(self.slot_n))
        self.slot_w += np.bincount(slot, weights=md, minlength=len(self.slot_w)).astype(np.int64)
        wf = lay.win_first[rid]
        other = wf < 0
        if other.any():
            cov[-1] += int(md[other].sum())
        todo = ~other & ~walked                                  # a spliced record that crosses a window was walked in the loop
        if not todo.any():
            return
        W = lay.window
        wf, s, m = wf[todo], start[todo], md[todo]
        e = s + m + N[todo]
        w0 = s // W
        w1 = (e - 1) // W
        one = w0 == w1
        size = len(cov)
        if one.any():
            cov += np.bincount(wf[one] + w0[one], weights=m[one], minlength=size).astype(np.int64)
        if not one.all():                                        # no skipped region here: every position from start to end is aligned
            many = ~one
            wf, s, e, w0, w1 = wf[many], s[many], e[many], w0[many], w1[many]
            cov += np.bincount(wf + w0, weights=(w0 + 1) * W - s, minlength=size).astype(np.int64)
            cov += np.bincount(wf + w1, weights=e - w1 * W, minlength=size).astype(np.int64)
            diff = np.bincount(wf + w0 + 1, minlength=size + 1) - np.bincount(wf + w1, minlength=size + 1)
            cov += np.cumsum(diff)[:size] * W

    def walk(self, lay: Layout, first_window: int, start: int, cigar) -> None:
        """A record with a skipped region (N) that crosses a window boundary: add its aligned columns window by
        window. Called from the read loop for these records only."""
        cov = self.ensure_cov(lay)
        W = lay.window
        pos = start
        for op, ln in cigar:
            if op in (0, 2, 7, 8):
                end = pos + ln
                w = pos // W
                while pos < end:
                    stop = min(end, (w + 1) * W)
                    cov[first_window + w] += stop - pos
                    pos = stop
                    w += 1
            elif op == 3:
                pos += ln

    def _transcriptome(self, lay: Layout, rid, start, span, ref_counts: bool) -> None:
        self.tx_reads += len(rid)
        if self.body is None:
            self.body = np.zeros((len(TX_CLASSES), BODY_BINS), np.int64)
            self.tx_class_n = np.zeros(len(TX_CLASSES), np.int64)
        if ref_counts and not self.ref_counts_dropped:
            if self.ref_n is None:
                self.ref_n = np.zeros(lay.n, np.int32)
            u, c = np.unique(rid, return_counts=True)
            self.ref_n[u] += c.astype(np.int32)
        elif not ref_counts:
            self.ref_counts_dropped = True
            self.ref_n = None
        L = lay.lengths[rid]
        cls = lay.tx_class[rid]
        self.tx_class_n += np.bincount(cls, minlength=len(TX_CLASSES))
        # the read covers bin b of 100 when it covers the middle of the bin: 200 s <= (2 b + 1) L < 200 e
        b0 = np.clip(-((L - 200 * start) // (2 * L)), 0, BODY_BINS)
        b1 = np.clip(-((L - 200 * (start + span)) // (2 * L)), 0, BODY_BINS)
        width = BODY_BINS + 1
        diff = np.bincount(cls * width + b0, minlength=len(TX_CLASSES) * width) - np.bincount(cls * width + b1, minlength=len(TX_CLASSES) * width)
        self.body += np.cumsum(diff.reshape(len(TX_CLASSES), width), axis=1)[:, :BODY_BINS]

    # ------------------------------------------------------------------------------------ merging
    def merge(self, o: AlignAcc) -> None:
        for k in _A_SCALARS:
            setattr(self, k, getattr(self, k) + getattr(o, k))
        self.h_alen += o.h_alen
        self.h_alen_w += o.h_alen_w
        self.h_clip += o.h_clip
        for k in ("g_len_aln", "g_aln_id", "g_q_id", "slot_n", "slot_w", "cov", "body", "tx_class_n"):
            setattr(self, k, _add_opt(getattr(self, k), getattr(o, k)))
        self.ref_counts_dropped = self.ref_counts_dropped or o.ref_counts_dropped
        self.ref_n = None if self.ref_counts_dropped else _add_opt(self.ref_n, o.ref_n)
        self.tt_alen.merge(o.tt_alen)
        self.tt_id.merge(o.tt_id)

    def to_state(self) -> dict:
        d = {k: getattr(self, k) for k in _A_SCALARS}
        d.update(h_alen=self.h_alen, h_alen_w=self.h_alen_w, h_clip=self.h_clip, g_len_aln=_grid_state(self.g_len_aln),
                 g_aln_id=_grid_state(self.g_aln_id), g_q_id=_grid_state(self.g_q_id), slot_n=self.slot_n, slot_w=self.slot_w,
                 cov=self.cov, body=self.body, tx_class_n=self.tx_class_n, ref_n=_grid_state(self.ref_n),
                 ref_len=None if self.ref_n is None else len(self.ref_n), ref_counts_dropped=self.ref_counts_dropped,
                 tt_alen=self.tt_alen.to_state(), tt_id=self.tt_id.to_state())
        return d

    @classmethod
    def from_state(cls, d: dict) -> AlignAcc:
        a = cls()
        for k in _A_SCALARS:
            setattr(a, k, d[k])
        a.h_alen, a.h_alen_w, a.h_clip = (np.asarray(d[k], np.int64).copy() for k in ("h_alen", "h_alen_w", "h_clip"))
        a.g_len_aln = _grid_from((_NL, _NL), d["g_len_aln"])
        a.g_aln_id = _grid_from((_NL, _NI), d["g_aln_id"])
        a.g_q_id = _grid_from((_NQ, _NI), d["g_q_id"])
        for k in ("slot_n", "slot_w", "cov", "body", "tx_class_n"):
            setattr(a, k, None if d[k] is None else np.asarray(d[k], np.int64).copy())
        a.ref_counts_dropped = d["ref_counts_dropped"]
        if d["ref_n"] is not None and not a.ref_counts_dropped:
            a.ref_n = np.zeros(d["ref_len"], np.int32)
            a.ref_n[d["ref_n"][0]] = d["ref_n"][1]
        a.tt_alen = TimeTable.from_state(_NL, d["tt_alen"])
        a.tt_id = TimeTable.from_state(_NI, d["tt_id"])
        return a

    @property
    def mapped(self) -> int:
        return self.fwd + self.rev

    # ------------------------------------------------------------------------------- stored cells
    def cell(self, layout: Layout | None, *, unmapped: int) -> dict:
        """The members this accumulator adds to the ``aln`` cell of its group (the bundle's aggregate JSON)."""
        sc = {"forward_reads": int(self.fwd), "reverse_reads": int(self.rev),
              "match_or_mismatch_bases": int(self.m_bases), "inserted_bases": int(self.ins_bases),
              "deleted_bases": int(self.del_bases), "skipped_bases": int(self.skip_bases),
              "soft_clipped_bases": int(self.soft_bases), "hard_clipped_bases": int(self.hard_bases),
              "mismatched_bases": int(self.mism_bases), "match_or_mismatch_bases_with_nm": int(self.m_bases_nm),
              "nm_below_indels_reads": int(self.nm_negative)}
        out: dict = {"scalars": sc, "hists": {
            "aln_aligned_length_hist": build.hist(LEN_AXIS, self.h_alen, self.h_alen_w, w_name="aligned_bases", missing=int(unmapped)),
            "aln_soft_clip_hist": build.hist(CLIP_AXIS, self.h_clip, missing=int(unmapped))}}
        if self.g_len_aln is not None:
            out["hists2d"] = {
                "aln_length_vs_aligned": _hist2d(GRID_LEN, GRID_LEN, self.g_len_aln, missing=int(unmapped)),
                "aln_identity_vs_aligned": _hist2d(GRID_LEN, GRID_ID, self.g_aln_id, missing=int(unmapped + self.mapped - self.g_aln_id.sum())),
                "aln_q_calibration": _hist2d(GRID_Q, GRID_ID, self.g_q_id, missing=int(unmapped + self.mapped - self.g_q_id.sum()))}
        if layout is not None and layout.kind == "genome" and self.cov is not None:
            names = list(layout.names) + [OTHER_CONTIGS]
            out["counts"] = {"aln_contig_table": build.counts(names, [int(v) for v in self.slot_n], w=[int(v) for v in self.slot_w],
                                                              role="categories")}
            out["counts"]["aln_contig_table"]["w_name"] = "aligned_bases"
            out["series"] = {"aln_depth_windows": build.series("genome_window", {"aligned_bases": self.cov[:-1]}, unit="bases",
                                                               n_bins=layout.n_windows)}
        return out

    def time_json(self) -> dict | None:
        d = {"aligned_length": self.tt_alen.to_json(), "identity": self.tt_id.to_json()}
        return d if any(v for v in d.values()) else None

    def tx_cell(self) -> dict | None:
        """The cell of the ``tx`` fragment: coverage along the transcript body per transcript length class."""
        if self.body is None:
            return None
        hists = {}
        scalars = {"reads": int(self.tx_reads)}
        for k, (cid, _) in enumerate(TX_CLASSES):
            full = np.concatenate([[0], self.body[k], [0]])
            hists["rna_body_coverage_" + cid] = build.hist("rel_pos_0p01", full)
            scalars["reads_" + cid] = int(self.tx_class_n[k])
        return {"scalars": scalars, "hists": hists}


def sparse_counts(a: np.ndarray) -> dict:
    """Per-reference counts in the stored form: the indexes that hold a count, as differences, and the counts."""
    idx = np.flatnonzero(a)
    return {"size": int(len(a)), "d": [int(x) for x in np.diff(idx, prepend=0)], "n": [int(x) for x in a[idx]]}


def dense_counts(d: dict) -> np.ndarray:
    a = np.zeros(int(d["size"]), np.int64)
    if d["d"]:
        a[np.cumsum(np.asarray(d["d"], np.int64))] = np.asarray(d["n"], np.int64)
    return a


# ----------------------------------------------------------------------------------------------- TagAcc
_T_SCALARS = ("reads", "timed", "dur_reads", "dur_bases", "dur_ms", "pt_tag", "pt_found", "pt_zero", "pt_none")


class TagAcc:
    """Every primary record of one group: what the basecaller's own tags say (st, du, qs, pt)."""
    __slots__ = _T_SCALARS + ("t10", "tt_len", "tt_q", "tt_speed", "h_speed", "h_pt")

    def __init__(self):
        for k in _T_SCALARS:
            setattr(self, k, 0)
        self.t10 = TimeTable(2)                                  # column 0 reads, column 1 bases
        self.tt_len = TimeTable(_NL)
        self.tt_q = TimeTable(_NQ)
        self.tt_speed = TimeTable(_NS)
        self.h_speed = np.zeros(_NS, np.int64)
        self.h_pt = None

    def add(self, R: np.ndarray) -> None:
        """A batch of primary records (``COLS`` by records)."""
        n = R.shape[1]
        if not n:
            return
        self.reads += n
        length, q, t10, du, pt = R[C_LEN], R[C_Q], R[C_T10], R[C_DU], R[C_PT]
        timed = t10 >= 0
        n_t = int(timed.sum())
        self.timed += n_t
        fast = du >= MIN_DURATION_S
        hour = None
        if n_t:
            tb = t10[timed].astype(np.int64)
            zero = np.zeros(n_t, np.int64)
            self.t10.add(tb, zero)
            self.t10.add(tb, zero + 1, w=length[timed])
            hour = np.floor_divide(t10, 6.0).astype(np.int64)
            self.tt_len.add(hour[timed], axes.bin_index(length[timed], GRID_LEN))
            qok = timed & np.isfinite(q)
            self.tt_q.add(hour[qok], axes.bin_index(q[qok], GRID_Q))
        if fast.any():
            ln, d = length[fast], du[fast]
            self.dur_reads += int(fast.sum())
            self.dur_bases += int(ln.sum())
            self.dur_ms += int(np.rint(d * 1000.0).sum())
            sb = axes.bin_index(ln / d, SPEED_AXIS)
            self.h_speed += np.bincount(sb, minlength=_NS)
            if hour is not None:
                ft = timed[fast]
                self.tt_speed.add(hour[fast][ft], sb[ft])
        tagged = pt > PT_NO_TAG
        n_tag = int(tagged.sum())
        if n_tag:
            self.pt_tag += n_tag
            found = pt > 0
            self.pt_found += int(found.sum())
            self.pt_zero += int((pt == 0).sum())
            self.pt_none += int((tagged & (pt < 0)).sum())
            if self.h_pt is None:
                self.h_pt = np.zeros(axes.n_cells(POLYA_AXIS), np.int64)
            self.h_pt += np.bincount(axes.bin_index(pt[found], POLYA_AXIS), minlength=len(self.h_pt))

    def merge(self, o: TagAcc) -> None:
        for k in _T_SCALARS:
            setattr(self, k, getattr(self, k) + getattr(o, k))
        for k in ("t10", "tt_len", "tt_q", "tt_speed"):
            getattr(self, k).merge(getattr(o, k))
        self.h_speed += o.h_speed
        self.h_pt = _add_opt(self.h_pt, o.h_pt)

    def to_state(self) -> dict:
        d = {k: getattr(self, k) for k in _T_SCALARS}
        d.update(t10=self.t10.to_state(), tt_len=self.tt_len.to_state(), tt_q=self.tt_q.to_state(), tt_speed=self.tt_speed.to_state(),
                 h_speed=self.h_speed, h_pt=self.h_pt)
        return d

    @classmethod
    def from_state(cls, d: dict) -> TagAcc:
        a = cls()
        for k in _T_SCALARS:
            setattr(a, k, d[k])
        a.t10 = TimeTable.from_state(2, d["t10"])
        a.tt_len = TimeTable.from_state(_NL, d["tt_len"])
        a.tt_q = TimeTable.from_state(_NQ, d["tt_q"])
        a.tt_speed = TimeTable.from_state(_NS, d["tt_speed"])
        a.h_speed = np.asarray(d["h_speed"], np.int64).copy()
        a.h_pt = None if d["h_pt"] is None else np.asarray(d["h_pt"], np.int64).copy()
        return a

    # ------------------------------------------------------------------------------- stored cells
    def time_cell(self) -> dict | None:
        """The cell of the ``time`` fragment, or None when no read carries a start time or a duration."""
        if not self.timed and not self.dur_reads:
            return None
        return {"scalars": {"reads": int(self.reads), "reads_with_start_time": int(self.timed),
                            "reads_with_duration": int(self.dur_reads), "bases_with_duration": int(self.dur_bases),
                            "duration_ms": int(self.dur_ms)},
                "hists": {"dyn_speed_basecalled_hist": build.hist(SPEED_AXIS, self.h_speed, missing=int(self.reads - self.dur_reads))}}

    def time_json(self) -> dict | None:
        d = {"output": self.t10.to_json(), "length": self.tt_len.to_json(), "quality": self.tt_q.to_json(), "speed": self.tt_speed.to_json()}
        return d if any(v for v in d.values()) else None

    def polya_cell(self) -> dict | None:
        """The cell of the ``polya`` fragment, or None when no read carries the pt tag."""
        if not self.pt_tag or self.h_pt is None:
            return None
        return {"scalars": {"reads": int(self.reads), "reads_with_tag": int(self.pt_tag), "tail_found_reads": int(self.pt_found),
                            "tail_zero_reads": int(self.pt_zero), "tail_not_found_reads": int(self.pt_none)},
                "hists": {"rna_polya_hist": build.hist(POLYA_AXIS, self.h_pt, missing=int(self.reads - self.pt_found))}}
