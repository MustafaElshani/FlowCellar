"""History: judged as of the run's start and then frozen; comparison groups, the robust score, the 5 and 8 run
steps, the spread guard; flow cell use by loadings and restarts; re-use."""
import copy

import pytest
from engine_support import clone, history_of, judged

from flowcellar.contracts import Declared
from flowcellar.engine import api
from flowcellar.engine import baselines as BL
from flowcellar.engine import flowcell_use as FU
from flowcellar.testing import minknow_variants as V

READS = [2_600_000, 2_800_000, 3_000_000, 3_100_000, 3_200_000, 3_300_000, 3_500_000, 3_700_000, 3_900_000]


def _run(day_run, reads=3_000_000, **kw):
    """The run that is judged: 24.2 h, started on 1 June 2026."""
    return clone(day_run.bundle, yield__yield_reads_at_24h=reads, yield__yield_estimated_bases_at_24h=3_000_000_000, **kw)


def _hist(day_run, n, reads=READS, **kw):
    return history_of(day_run.bundle, n, yield__yield_reads_at_24h=reads[:n], yield__yield_estimated_bases_at_24h=[3_000_000_000] * n, **kw)


def _state(st, rid):
    r = st["_rules"][rid]
    return "pass" if r["outcome"] == "ok" else r["severity"]


def test_robust_score():
    b = BL.robust(50.0, [90, 100, 100, 110, 100, 95, 105, 100], 0.15)
    assert b["median"] == 100 and b["n"] == 8 and b["scale"] == 15.0            # the floor: 15 % of the median
    assert b["z"] == pytest.approx(-50 / 15) and b["rank_below"] == 0
    wide = BL.robust(100.0, [20, 60, 100, 140, 180, 30, 170, 100], 0.15)
    assert wide["scale"] == pytest.approx(1.4826 * 55) and wide["rel"] > 0.45           # 1.4826 x the median absolute deviation
    assert BL.robust(1.0, [], 0.1) is None and BL.robust(None, [1, 2, 3], 0.1) is None


def test_fewer_than_five_earlier_runs_is_no_history(day_run):
    st = judged(_run(day_run), _hist(day_run, 4))
    r = st["_rules"]["yield_reads_at_24h"]
    assert "baseline" not in r["threshold"] and r["threshold"]["applied"] == "fixed" and r["threshold"]["peers"] == 4
    assert "no_history" in [m["code"] for m in st["modifiers"]] and st["explanations"] == []
    assert "No lab history for this combination yet" in r["message"]


def test_five_to_seven_runs_is_an_early_comparison_that_can_warn_but_not_fail(day_run):
    hist = _hist(day_run, 6)
    st = judged(_run(day_run, reads=3_050_000), hist)
    r = st["_rules"]["yield_reads_at_24h"]
    bl = r["threshold"]["baseline"]
    assert (bl["trust"], bl["n"], bl["spread_guard"], bl["horizon_h"]) == ("provisional", 6, False, 24)
    assert bl["as_of"] == "2026-06-01T08:00:00Z" and bl["group"] == "MINION_MK1D | FLO-MIN114 | SQK-PCB114-24"
    assert r["outcome"] == "ok" and [m["code"] for m in st["modifiers"]] == ["baseline_provisional"]
    assert r["message"] == ("Reads at 24 h: 3.0 M. Typical here is 2.4 to 3.7 M (6 comparable runs). "
                            "Lab typical rests on 6 comparable runs: it can warn but not fail a run until there are 8.")
    assert r["threshold"]["warn"] == 2250000 and r["threshold"]["fail"] == 1125000 and r["threshold"]["applied"] == "fixed"
    low = judged(_run(day_run, reads=1_500_000), hist)                   # z far below -2.5, above the fixed fail level
    r = low["_rules"]["yield_reads_at_24h"]
    assert r["threshold"]["baseline"]["z"] <= -2.5 and (r["outcome"], r["severity"]) == ("triggered", "warn") and low["verdict"] == "warn"
    assert r["threshold"]["applied"] == "baseline" and "This lab's warning level at 24 h is" in r["message"]
    assert judged(_run(day_run, reads=1_000_000), hist)["_rules"]["yield_reads_at_24h"]["severity"] == "fail"       # the fixed level still fails


def test_eight_runs_is_a_trusted_history_that_can_fail(day_run):
    hist = _hist(day_run, 9)
    ok = judged(_run(day_run, reads=3_100_000), hist)
    r = ok["_rules"]["yield_reads_at_24h"]
    bl = r["threshold"]["baseline"]
    assert (bl["trust"], bl["n"]) == ("trusted", 9) and bl["median"] == 3_200_000 and r["threshold"]["applied"] == "baseline"
    assert r["threshold"]["warn"] is None and bl["warn_z"] == -1.5 and bl["fail_z"] == -2.5
    assert bl["warn_level"] == pytest.approx(bl["median"] - 1.5 * bl["scale"]) and bl["band"][0] == bl["warn_level"]
    assert ok["verdict"] == "pass" and not [m for m in ok["modifiers"] if m["code"].startswith(("baseline", "no_history"))]
    assert r["message"] == "Reads at 24 h: 3.1 M. Typical here is 2.5 to 3.9 M (9 comparable runs)."
    warn = judged(_run(day_run, reads=2_400_000), hist)
    assert _state(warn, "yield_reads_at_24h") == "warn" and -2.5 < warn["_rules"]["yield_reads_at_24h"]["threshold"]["baseline"]["z"] <= -1.5
    assert warn["_rules"]["yield_reads_at_24h"]["message"].startswith("Reads at 24 h: 2.4 M. This lab's warning level at 24 h is 2.5 M (9 comparable runs).")
    fail = judged(_run(day_run, reads=1_900_000), hist)
    assert _state(fail, "yield_reads_at_24h") == "fail" and fail["verdict"] == "fail"
    from flowcellar.model import schema
    assert schema.schema_errors(api_compose(_run(day_run, reads=1_900_000), fail)) == []
    # estimated bases beside the read count: placed against history, a warning at most (z <= -2.5)
    b = clone(_run(day_run, reads=3_100_000), yield__yield_estimated_bases_at_24h=1_000_000_000)
    assert _state(judged(b, hist), "yield_estimated_bases_at_24h") == "warn"


def api_compose(bundle, st):
    from flowcellar.model import bundle as B
    clean = {k: v for k, v in st.items() if not k.startswith("_")}
    return B.compose(bundle, status=clean, relations=st["_file"]["relations"])


def test_spread_guard_makes_a_wide_history_provisional(day_run):
    wide = [500_000, 900_000, 1_500_000, 2_500_000, 3_000_000, 4_500_000, 6_000_000, 7_000_000, 8_000_000]
    st = judged(_run(day_run, reads=600_000), _hist(day_run, 9, reads=wide))
    r = st["_rules"]["yield_reads_at_24h"]
    bl = r["threshold"]["baseline"]
    assert bl["spread"] > 0.45 and bl["spread_guard"] is True and bl["trust"] == "provisional" and bl["n"] == 9
    assert r["severity"] != "fail" or r["threshold"]["applied"] == "fixed"
    assert "baseline_spread_guard" in [m["code"] for m in st["modifiers"]]
    assert "Comparable runs differ widely from each other (9 runs)" in r["message"]


def test_peers_are_first_use_same_group_earlier_long_enough_and_not_the_same_day(day_run):
    run = _run(day_run)
    hist = _hist(day_run, 9)
    n = lambda cards: judged(run, cards)["_rules"]["yield_reads_at_24h"]["threshold"].get("baseline", {}).get("n", 0)  # noqa: E731
    assert n(hist) == 9
    later = copy.deepcopy(hist)
    later[0]["facts"]["started"] = "2026-06-02T09:00:00Z"                # started after the run
    later[0]["facts"]["day"] = "2026-06-02"
    assert n(later) == 8
    same_day = copy.deepcopy(hist)
    same_day[0]["facts"].update(started="2026-06-01T01:00:00Z", day="2026-06-01")
    assert n(same_day) == 8
    other_kit = copy.deepcopy(hist)
    other_kit[0]["facts"]["kit"] = "SQK-PCB111-24"
    assert n(other_kit) == 8
    short = copy.deepcopy(hist)
    short[0]["facts"]["run_hours"] = 20.0
    assert n(short) == 8
    quick_only = copy.deepcopy(hist)                                     # tier 1 not complete does not teach
    quick_only[0]["facts"]["t1_complete"] = False
    assert n(quick_only) == 8
    reused = copy.deepcopy(hist)                                         # a second loading of a flow cell never teaches
    reused[1]["facts"]["flow_cell_id"] = reused[0]["facts"]["flow_cell_id"]
    assert n(reused) == 8
    assert n(hist + [api.card(run, None)]) == 9                          # a run is never compared with itself


def test_window_is_the_last_forty_runs(day_run):
    cards = history_of(day_run.bundle, 28, month="2026-01", yield__yield_reads_at_24h=[1_000_000] * 28,
                       yield__yield_estimated_bases_at_24h=[3_000_000_000] * 28)
    cards += [c for c in history_of(day_run.bundle, 28, month="2026-02", yield__yield_reads_at_24h=[3_000_000] * 28,
                                    yield__yield_estimated_bases_at_24h=[3_000_000_000] * 28)]
    for i, c in enumerate(cards[28:]):
        c["facts"]["run_key"] = f"TST{30000 + i:05d}_{i + 1:08x}"
        c["facts"]["flow_cell_id"] = f"TST{30000 + i:05d}"
    bl = judged(_run(day_run), cards)["_rules"]["yield_reads_at_24h"]["threshold"]["baseline"]
    assert bl["n"] == 40 and bl["median"] == 3_000_000                   # 28 newer runs and the 12 newest of the older ones
    bl = judged(_run(day_run), cards, profile={"_baselines": {"window": 10}})["_rules"]["yield_reads_at_24h"]["threshold"]["baseline"]
    assert bl["n"] == 10


def test_group_label_splits_a_comparison_group(day_run):
    hist = _hist(day_run, 9)
    groups = [{"match": {"kit": "SQK-PCB114-24", "started_after": "2026-05-01"}, "label": "new protocol"}]
    st = judged(_run(day_run), hist, baseline_groups=groups)
    assert "baseline" not in st["_rules"]["yield_reads_at_24h"]["threshold"]             # the January runs are another group now
    assert st["_file"]["engine"]["group"] == "MINION_MK1D | FLO-MIN114 | SQK-PCB114-24 | new protocol"
    st = judged(_run(day_run), hist, declared=Declared(baseline_group="pilot"))
    assert "baseline" not in st["_rules"]["yield_reads_at_24h"]["threshold"]
    card = api.card(_run(day_run), st["_file"])
    assert card["columns"]["group"].endswith("| pilot") and card["facts"]["group"].endswith("| pilot")


def test_explanations_place_numbers_against_history_and_never_judge(day_run):
    n = 9
    hist = _hist(day_run, n, fc__fc_pores_first_scan=[1400, 1450, 1480, 1500, 1500, 1520, 1550, 1580, 1600],
                 fc__fc_occupancy_median_24h=[0.90] * n)
    st = judged(_run(day_run, fc__fc_pores_first_scan=900, fc__fc_occupancy_median_24h=0.70), hist)
    ex = {e["id"]: e for e in st["explanations"]}
    assert ex["fc_pores_first_scan"]["direction"] == "low" and ex["fc_pores_first_scan"]["aspect"] == "flow_cell"
    assert ex["fc_occupancy_median_24h"]["direction"] == "low" and ex["fc_occupancy_median_24h"]["aspect"] == "loading"
    assert ex["fc_pores_first_scan"]["baseline"]["trust"] == "trusted" and ex["fc_pores_first_scan"]["message"].startswith("Pores at the first scan: 900.")
    assert st["verdict"] == "pass"                                       # explanations never change a first-use verdict
    assert "Possibly re-used" in ex["fc_pores_first_scan"]["message"] and "flow_cell_use_unknown" in [m["code"] for m in st["modifiers"]]
    card = api.card(_run(day_run), st["_file"])["columns"]
    assert (card["mark_flow_cell"], card["mark_loading"]) == (-1, -1)
    assert card["outcome_metric"] == "yield_reads_at_24h" and card["outcome_n"] == 9 and card["outcome_trust"] == "trusted"
    assert card["outcome_level"] == st["_rules"]["yield_reads_at_24h"]["threshold"]["baseline"]["warn_level"]


def test_adapter_share_is_compared_within_one_minknow_version(day_run):
    hist = _hist(day_run, 9, fc__fc_adapter_frac=[0.07] * 9)
    st = judged(_run(day_run, fc__fc_adapter_frac=0.30), hist)
    assert {e["id"]: e["direction"] for e in st["explanations"]}["fc_adapter_frac"] == "high"
    for c in hist[:6]:
        c["facts"]["mk"] = "5.9"
    st = judged(_run(day_run, fc__fc_adapter_frac=0.30), hist)
    assert "fc_adapter_frac" not in {e["id"] for e in st["explanations"]}          # only three runs of this MinKNOW version are left


# --------------------------------------------------------------------------------------------- frozen
def test_a_later_run_never_changes_an_earlier_verdict(day_run):
    hist = _hist(day_run, 9)
    run = _run(day_run, reads=2_400_000)
    before = api.judge(run, hist)
    later = history_of(day_run.bundle, 20, month="2026-07", yield__yield_reads_at_24h=[9_000_000] * 20,
                       yield__yield_estimated_bases_at_24h=[9_000_000_000] * 20)
    for i, c in enumerate(later):
        c["facts"]["run_key"] = f"TST{40000 + i:05d}_{i + 1:08x}"
        c["facts"]["flow_cell_id"] = f"TST{40000 + i:05d}"
    after = api.judge(run, hist + later)

    def content(sf):
        st = copy.deepcopy(sf["status"])
        st.pop("evaluated")
        st["judged"].pop("at")
        return st
    assert content(before) == content(after)
    assert before["status"]["judged"]["history_hash"] == after["status"]["judged"]["history_hash"] is not None
    assert before["status"]["judged"]["as_of"] == "2026-06-01T08:00:00Z" and before["status"]["judged"]["frozen"] is True
    card = api.card(run, before)
    assert api.history_hash(card, hist + later) == before["status"]["judged"]["history_hash"]
    # an OLDER run added afterwards does change the hash: the library can flag the run without judging it again
    older = history_of(day_run.bundle, 1, month="2025-12", yield__yield_reads_at_24h=[3_000_000], yield__yield_estimated_bases_at_24h=[3_000_000_000])
    older[0]["facts"].update(run_key="TST50000_0000abcd", flow_cell_id="TST50000")
    assert api.history_hash(card, hist + older) != before["status"]["judged"]["history_hash"]
    # an earlier run that gets its read count later does not change it
    counted = copy.deepcopy(hist)
    counted[0]["facts"]["speed_med"] = 371.0
    assert api.history_hash(card, counted) == before["status"]["judged"]["history_hash"]


def test_judging_in_order_of_start_shows_what_a_library_kept_from_day_one_would_show(day_run):
    bundles = []
    for i in range(10):
        bundles.append(clone(day_run.bundle, key=f"TST{60000 + i:05d}_{i + 1:08x}", flow_cell=f"TST{60000 + i:05d}",
                             started=f"2026-03-{i + 1:02d}T09:00:00Z", yield__yield_reads_at_24h=READS[i % 9],
                             yield__yield_estimated_bases_at_24h=3_000_000_000))
    out = api.judge_all(list(reversed(bundles)))                         # given newest first: judged oldest first
    ns = [(out[b["bundle_id"]][0]["status"]["rules"], b["bundle_id"]) for b in bundles]
    peers = []
    for rules, _key in ns:
        r = next(x for x in rules if x["id"] == "yield_reads_at_24h")
        peers.append((r["threshold"].get("baseline") or {}).get("n", r["threshold"].get("peers", 0)))
    assert peers == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    trusts = [next(x for x in rules if x["id"] == "yield_reads_at_24h")["threshold"].get("baseline", {}).get("trust", "none") for rules, _ in ns]
    assert trusts == ["none"] * 5 + ["provisional"] * 3 + ["trusted"] * 2
    # one more run arrives: nothing that was judged changes when it is judged against the cards that exist
    known = [c for _sf, c in out.values()]
    extra = clone(day_run.bundle, key="TST60099_000000ff", flow_cell="TST60099", started="2026-03-20T09:00:00Z",
                  yield__yield_reads_at_24h=3_000_000, yield__yield_estimated_bases_at_24h=3_000_000_000)
    more = api.judge_all([extra], known=known)
    assert list(more) == ["TST60099_000000ff"]
    assert next(x for x in more["TST60099_000000ff"][0]["status"]["rules"] if x["id"] == "yield_reads_at_24h")["threshold"]["baseline"]["n"] == 10
    today = api.baselines_today(known + [more["TST60099_000000ff"][1]])
    g = today["groups"][0]
    assert g["n_runs"] == 11 and g["metrics"]["yield_reads_at_24h"]["trust"] == "trusted" and g["metrics"]["yield_reads_at_24h"]["n"] == 11
    assert len(g["curve_output_gb_per_3h"]["median"]) == 32


# ---------------------------------------------------------------------------------------- flow cell use
def _f(key, fc, started, hours, bases, **kw):
    return dict(run_key=key, flow_cell_id=fc, started=started, run_hours=hours, est_bases=bases, family="promethion", agroup="cdna",
                has_tier1=True, mux0_pores=5000, **kw)


def test_loadings_and_restarts():
    runs = [_f("A_1", "TST00001", "2026-01-01T10:00:00Z", 44.0, 38e9),
            _f("A_2", "TST00001", "2026-01-03T10:00:00Z", 0.1, 0),            # 5 minutes: a restart
            _f("A_3", "TST00001", "2026-01-03T10:20:00Z", 1.5, 0.5e9),        # 1.5 h, 0.5 Gb: a restart (under 0.54 Gb)
            _f("A_4", "TST00001", "2026-01-03T12:00:00Z", 67.0, 7e9),
            _f("B_1", "TST00002", "2026-01-02T10:00:00Z", 1.0, 0.6e9),        # short, but 0.6 Gb: a loading
            _f("B_2", "TST00002", "2026-01-04T10:00:00Z", 30.0, 9e9)]
    u = FU.annotate(runs)
    assert [u[k]["index"] for k in ("A_1", "A_2", "A_3", "A_4")] == [1, 2, 2, 2] and all(u[k]["of"] == 2 for k in ("A_1", "A_4"))
    assert [u[k]["restart"] for k in ("A_1", "A_2", "A_3", "A_4")] == [False, True, True, False]
    assert (u["A_2"]["restart_of"], u["A_3"]["restart_of"], u["A_4"]["restart_of"]) == ("A_4", "A_4", None)
    assert u["A_4"]["prior_runs"] == 3 and u["A_4"]["prior_hours"] == 45.6 and u["A_4"]["prior_estimated_bases"] == int(38.5e9)
    assert [u[k]["index"] for k in ("B_1", "B_2")] == [1, 2] and not u["B_1"]["restart"]
    assert [r["loading"] for r in u["A_1"]["runs"]] == [1, 2, 2, 2]
    use = FU.use_of(runs[3], runs[:3] + runs[4:])
    assert (use["index"], use["source"]) == (2, "library")
    assert FU.use_of(runs[0], [])["source"] == "unknown" and FU.use_of(runs[0], [], declared_use=3)["index"] == 3
    assert FU.use_of(runs[2], runs[:2], {"restart_max_hours": 1.0})["restart"] is False        # a configured restart length
    stub = dict(run_key="A_0", flow_cell_id="TST00001", started="2025-12-30T10:00:00Z", has_tier1=False, family="promethion", agroup="unknown")
    assert FU.annotate([stub] + runs)["A_1"]["index"] == 1               # a folder with raw data only is no loading


def test_reused_flow_cell_with_restarts_on_fake_runs(tmp_path, build):
    fixtures = V.make_reused_flow_cell(tmp_path)
    runs = [api.identify(f.run_dir) for f in fixtures]
    bundles = [api.quick_bundle(r) for r in runs]
    out = api.judge_all(bundles)
    sts = [out[f.run_key][0] for f in fixtures]
    use = [s["relations"]["flow_cell_use"] for s in sts]
    assert [u["index"] for u in use] == [1, 2, 2, 2] and [u["source"] for u in use] == ["unknown", "library", "library", "library"]
    assert use[1]["restart_of"] == use[2]["restart_of"] == fixtures[3].run_key and use[3]["restart_of"] is None
    assert len(use[3]["runs"]) == 4 and [r["restart"] for r in use[3]["runs"]] == [False, True, True, False]
    first, five_min, ninety_min, second = (s["status"] for s in sts)
    assert (five_min["verdict"], five_min["cause"], five_min["modifiers"]) == ("incomplete", "stopped_early", [])
    assert ninety_min["verdict"] == "incomplete" and ninety_min["cause"] == "output_not_judged"
    assert [m["code"] for m in ninety_min["modifiers"]][0] == "reused_flow_cell"
    assert second["modifiers"][0] == {"code": "reused_flow_cell", "message": "Re-used flow cell (use 2): judged on output per starting pore, "
                                                                           "occupancy, speed and read quality only."}
    ids = {r["id"] for r in second["rules"]}
    assert not (ids & {"fc_pores_first_scan", "fc_channels_with_pore_first_scan", "yield_reads_at_24h", "yield_estimated_bases_at_24h"})
    assert "fc_pores_first_scan" in {r["id"] for r in first["rules"]}
    cards = {f.run_key: out[f.run_key][1]["columns"] for f in fixtures}
    assert [cards[f.run_key]["use"] for f in fixtures] == [1, 2, 2, 2] and [cards[f.run_key]["restart"] for f in fixtures] == [False, True, True, False]
    cols = api.flow_cell_columns([out[f.run_key][1] for f in fixtures])
    assert cols[fixtures[0].run_key] == {"use": 1, "use_of": 2, "restart": False, "restart_of": None}        # the first loading learns of the second
    assert cols[fixtures[1].run_key]["restart_of"] == fixtures[3].run_key
    # adding the later runs did not change the status of the first one
    alone = api.judge(bundles[0], [])
    for sf in (alone, sts[0]):
        sf["status"].pop("evaluated")
        sf["status"]["judged"].pop("at")
    assert alone["status"] == sts[0]["status"]


def test_a_reuse_is_judged_on_output_per_starting_pore_against_first_use_history(day_run):
    hist = _hist(day_run, 9, fc__fc_pores_first_scan=[1500] * 9)
    earlier = clone(day_run.bundle, key="TST00030_0000aaaa", started="2026-05-20T08:00:00Z")           # the same flow cell, 12 days before
    cards = hist + [api.card(earlier, None)]
    run = _run(day_run, fc__fc_pores_first_scan=700)
    st = judged(run, cards)
    assert st["_file"]["relations"]["flow_cell_use"]["index"] == 2
    r = st["_rules"]["yield_per_pore_mb_24h"]
    assert r["role"] == "key" and r["threshold"]["applied"] == "baseline"
    assert r["threshold"]["baseline"]["n"] == 10                         # the nine runs and the first loading of this flow cell itself
    assert r["value"] == pytest.approx(3_000_000_000 / 700 / 1e6, abs=0.01) and r["outcome"] == "ok"
    assert "yield_reads_at_24h" not in st["_rules"] and "fc_pores_first_scan" not in st["_rules"]
    assert not any(e["id"] in ("fc_pores_first_scan", "fc_pore_survival_24h_frac", "fc_sequencing_duty") for e in st["explanations"])
    assert st["verdict"] == "pass" and st["modifiers"][0]["code"] == "reused_flow_cell"
    card = api.card(run, st["_file"])["columns"]
    assert (card["use"], card["use_source"], card["outcome_metric"]) == (2, "library", "yield_per_pore_mb_24h")
    low = judged(clone(run, yield__yield_estimated_bases_at_24h=300_000_000), cards)
    assert low["_rules"]["yield_per_pore_mb_24h"]["severity"] == "fail" and low["verdict"] == "fail"
    # declared with no library at all
    st = judged(run, (), declared=Declared(flow_cell_use=2), in_library=False)
    assert st["_file"]["relations"]["flow_cell_use"] == {**st["_file"]["relations"]["flow_cell_use"], "index": 2, "source": "declared"}
    assert st["modifiers"][0]["code"] == "reused_flow_cell" and "fc_pores_first_scan" not in st["_rules"] and st["judged"]["library"] is False
