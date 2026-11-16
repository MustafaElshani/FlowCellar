"""Merge laws and derivations of the aggregate kinds (section 16 "Merge laws")."""
import itertools

import numpy as np
import pytest

from flowcellar.aggregates import NotMergeable, build, derive, merge
from flowcellar.model import axes

RNG = np.random.default_rng(7)


def _len_hist(lengths, axis="read_length_log40"):
    idx = axes.bin_index(lengths, axis)
    n = np.bincount(idx, minlength=axes.n_cells(axis))
    w = np.bincount(idx, weights=lengths, minlength=axes.n_cells(axis)).astype(np.int64)
    return build.hist(axis, n, w, w_name="bases", vmin=int(lengths.min()), vmax=int(lengths.max()))


def _parts():
    return [np.rint(RNG.lognormal(7, 1, size=s)).astype(np.int64).clip(1) for s in (500, 800, 3)]


def test_hist_merge_equals_aggregate_of_the_union():
    a, b, c = _parts()
    whole = _len_hist(np.concatenate([a, b, c]))
    merged = merge.merge_hist(merge.merge_hist(_len_hist(a), _len_hist(b)), _len_hist(c))
    assert merged == whole
    assert derive.total(whole) == len(a) + len(b) + len(c)
    assert derive.total(whole, "w") == int(a.sum() + b.sum() + c.sum())      # weights make totals exact


def test_hist_merge_is_associative_and_commutative():
    hs = [_len_hist(p) for p in _parts()]
    results = []
    for order in itertools.permutations(range(3)):
        x = merge.merge_hist(merge.merge_hist(hs[order[0]], hs[order[1]]), hs[order[2]])
        y = merge.merge_hist(hs[order[0]], merge.merge_hist(hs[order[1]], hs[order[2]]))
        results += [x, y]
    assert all(r == results[0] for r in results)


def test_under_and_overflow_keep_their_weight():
    lengths = np.array([3, 5, 12, 400, 20_000_000])            # two below 10 b, one above 10 Mb
    h = _len_hist(lengths)
    assert (h["under"], h["over"]) == (2, 1)
    assert (h["w_under"], h["w_over"]) == (8, 20_000_000)
    assert derive.total(h, "w") == int(lengths.sum())
    assert derive.mean(h) == pytest.approx(lengths.mean())


def test_rebinning_then_merging_equals_merging_then_rebinning():
    a, b, _ = _parts()
    fine = merge.merge_hist(_len_hist(a), _len_hist(b))
    mixed = merge.merge_hist(_len_hist(a), _len_hist(b, "read_length_log20"))      # two levels of one family
    assert mixed["axis"] == "read_length_log20"
    assert mixed == merge.merge_hist(_len_hist(a, "read_length_log20"), _len_hist(b, "read_length_log20"))
    assert derive.total(mixed) == derive.total(fine) and derive.total(mixed, "w") == derive.total(fine, "w")
    with pytest.raises(NotMergeable):
        merge.merge_hist(_len_hist(a), build.hist("qscore_1", np.zeros(axes.n_cells("qscore_1"), int)))


def test_quantile_and_n50_are_close_to_the_exact_values():
    lengths = np.rint(RNG.lognormal(7.5, 0.8, size=50_000)).astype(np.int64).clip(10)
    h = _len_hist(lengths)
    assert derive.median(h) == pytest.approx(np.median(lengths), rel=0.02)
    s = np.sort(lengths)[::-1]
    n50 = s[np.searchsorted(np.cumsum(s), s.sum() / 2)]
    assert derive.nx(h, 0.5) == pytest.approx(n50, rel=0.03)
    assert derive.frac_above(h, 1000) == pytest.approx((lengths >= 1000).mean(), abs=1e-9)   # 1,000 b is an exact edge


def test_hist2d_merge():
    m1 = RNG.integers(0, 3, size=(6, 5)); m2 = RNG.integers(0, 3, size=(6, 5))
    a = build.hist2d("read_length_log20", "qscore_1", m1, outside=2)
    b = build.hist2d("read_length_log20", "qscore_1", m2, outside=1, missing=4)
    ab, ba = merge.merge_hist2d(a, b), merge.merge_hist2d(b, a)
    assert ab == ba == {**build.hist2d("read_length_log20", "qscore_1", m1 + m2, outside=3), "missing": 4}


def test_series_increments_add_and_levels_follow_their_rule():
    a = build.series("run_time_600s", {"reads": [1, 2, 3, 0, 0]})
    b = build.series("run_time_600s", {"reads": [10, 0, 0, 0, 7]})
    assert merge.merge_series(a, b)["values"]["reads"] == [11, 2, 3, 0, 7]
    assert merge.merge_series(a, b) == merge.merge_series(b, a)
    t1 = build.series("mux_scan_index", {"single_pore": [5, 4], "t_s": [0, 5400]}, levels={"t_s": "same"})
    t2 = build.series("mux_scan_index", {"single_pore": [1, 1], "t_s": [0, 5400]}, levels={"t_s": "same"})
    m = merge.merge_series(t1, t2)
    assert m["values"] == {"single_pore": [6, 5], "t_s": [0, 5400]} and m["levels"] == ["t_s"]
    lv = build.series("run_time_7200s", {"flowcell_c": [34.0, 34.1]}, levels={"flowcell_c": "mean"})
    with pytest.raises(NotMergeable):                 # a mean without a weight cannot be merged
        merge.merge_series(lv, lv)


def test_counts_merge_keeps_unknown_stages_unknown():
    a = build.counts(["sequenced", "basecalled", "given_to_a_sample"], [100, 0, 60], role="funnel", missing=["basecalled"])
    b = build.counts(["sequenced", "basecalled", "given_to_a_sample"], [50, 45, 30], role="funnel")
    m = merge.merge_counts(a, b)
    assert m["n"] == [150, 45, 90] and m["missing"] == ["basecalled"] and m["role"] == "funnel"
    assert merge.merge_counts(b, a)["missing"] == ["basecalled"]


def test_channels_and_scalars_rules():
    a = {"layout": "minion_512", "values": {"reads": [1, 2], "last_read_time": [10, 50]}, "merge": {"last_read_time": "max"}}
    b = {"layout": "minion_512", "values": {"reads": [3, 4], "last_read_time": [70, 20]}, "merge": {"last_read_time": "max"}}
    assert merge.merge_channels(a, b)["values"] == {"reads": [4, 6], "last_read_time": [70, 50]}
    s = merge.merge_scalars({"reads": 5, "max_events": 9, "x": None}, {"reads": 7, "max_events": 4, "x": 3}, {"max_events": "max"})
    assert s == {"reads": 12, "max_events": 9, "x": None}          # a sum with an unknown part is unknown


def test_rollup_over_samples_and_filters():
    dims = ["sample", "filter"]
    cells = []
    for sample, filt, n in (("s1", "pass", 10), ("s1", "fail", 2), ("s2", "pass", 20), ("s2", "fail", 5)):
        lengths = np.full(n, 500 if sample == "s1" else 2000)
        cells.append({"at": [sample, filt], "scalars": {"reads": n, "bases": int(lengths.sum())}, "hists": {"len_hist": _len_hist(lengths)}})
    total = merge.rollup(cells, dims, [])
    assert len(total) == 1 and total[0]["at"] == [] and total[0]["scalars"] == {"reads": 37, "bases": 12 * 500 + 25 * 2000}
    assert derive.total(total[0]["hists"]["len_hist"]) == 37
    per_sample = {c["at"][0]: c["scalars"]["reads"] for c in merge.rollup(cells, dims, ["sample"])}
    assert per_sample == {"s1": 12, "s2": 25}
    passed = merge.rollup(cells, dims, [], where={"filter": ["pass"]})
    assert passed[0]["scalars"]["reads"] == 30
    assert cells[0]["scalars"]["reads"] == 10                       # inputs are left untouched


def test_rollup_drops_curves_however_many_cells_hold_one():
    dims = ["sample", "filter"]
    curve = {"x": [0, 1], "y": {"reads": [1, 2]}}
    cells = [{"at": [s, "pass"], "scalars": {"reads": n}, "curves": {"demo_curve": dict(curve)}} for s, n in (("s1", 1), ("s2", 2), ("s3", 4), ("s4", 8))]
    total = merge.rollup(cells, dims, [])
    assert total[0]["scalars"] == {"reads": 15} and "curves" not in total[0], "a curve cannot be added, so the total has none"
    kept = merge.rollup(cells, dims, dims)
    assert all("curves" in c for c in kept) and len(kept) == 4, "nothing is rolled away: every cell keeps its curve"
    assert all("curves" in c for c in cells)                        # inputs are left untouched


def test_counts_keep_the_name_of_their_weight_and_add_what_was_folded_away():
    a = {"cats": ["x", "y"], "n": [1, 2], "w": [10.0, 20.0], "w_name": "aligned_bases", "other": 3}
    b = {"cats": ["y", "z"], "n": [5, 7], "w": [50.0, 70.0], "w_name": "aligned_bases"}
    m = merge.merge_counts(a, b)
    assert m == {"cats": ["x", "y", "z"], "n": [1, 7, 7], "w": [10.0, 70.0, 70.0], "w_name": "aligned_bases", "other": 3}
    assert merge.merge_counts(merge.merge_counts(a, b), a)["other"] == 6


def test_largest_remainder_and_sparkline():
    assert build.largest_remainder([0.5, 0.5, 2.0], 3) in ([1, 0, 2], [0, 1, 2])
    assert sum(build.largest_remainder([1.2, 3.3, 5.5], 10)) == 10
    text, mx = build.sparkline([0, 1, 2, 4])
    assert len(text) == 32 and mx == 4 and text[3] == "_" and text[0] == "A"
    assert list(build.u6_decode(build.u6_encode([0, 5, 63]))) == [0, 5, 63]
