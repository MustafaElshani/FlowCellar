"""Every tier 1 check on a run built to trigger it and on one built not to; the verdict table; quick
detail with its pending checks; the stored sentences."""
import copy
import re

import pytest
from engine_support import banned, clone, judged, problems, status_sentences

from flowcellar import catalogue
from flowcellar.contracts import Declared
from flowcellar.engine import api
from flowcellar.engine import copy_pending as CP
from flowcellar.engine.limits import PENDING, Limits
from flowcellar.model import status as S
from flowcellar.testing import minknow_variants as V
from flowcellar.testing.fixtures import RunSpec, make_run


def _state(st, rid):
    r = st["_rules"][rid]
    return "pass" if r["outcome"] == "ok" else r["severity"] if r["outcome"] == "triggered" else r["outcome"]


# ------------------------------------------------------------------- checks made on fake run folders
def test_a_clean_run_triggers_nothing_but_its_tiny_output(full8):
    st = full8.status
    trig = {r["id"] for r in st["rules"] if r["outcome"] == "triggered"}
    assert trig == {"yield_reads_at_6h", "yield_estimated_bases_at_6h"}          # 4,000 reads against a MinION's fixed level
    assert {r["id"] for r in st["rules"]} == {
        "run_stop_cause", "run_messages", "run_output_files", "fc_pores_first_scan", "fc_channels_with_pore_first_scan",
        "fc_occupancy_first_hour", "fc_recovering_frac", "fc_temperature_max_dev_c", "fc_active_channels", "dyn_speed_est_median",
        "dyn_signal_positive_frac", "dyn_unblock_frac", "len_est_n50", "yield_reads_at_6h", "yield_estimated_bases_at_6h"}
    assert st["by_tier"] == {"run_health": "fail", "read_level": "not_available"} and st["verdict"] == "fail"
    assert all(r["tier"] == "run_health" and r["panel"] for r in st["rules"])


def test_stop_cause(tmp_path, build):
    ok = build(make_run(tmp_path / "a", RunSpec(flow_cell_id="TST00070", seed=70, state="PROTOCOL_STOPPED_BY_USER")), full=False)
    assert ok.rules["run_stop_cause"]["outcome"] == "ok" and ok.rules["run_stop_cause"]["message"] == "The run was stopped by a person after 8.0 hours."
    bad = build(V.set_stop(make_run(tmp_path / "b", RunSpec(flow_cell_id="TST00071", seed=71))), full=False)
    r = bad.rules["run_stop_cause"]
    assert (r["outcome"], r["severity"], r["role"]) == ("triggered", "fail", "key") and r["message"] == "How the run ended: device error."
    assert bad.status["by_tier"]["run_health"] == "fail"


def test_messages(tmp_path, build):
    err = build(V.make_with_message(tmp_path / "a", flow_cell_id="TST00072", seed=72), full=False)
    assert err.rules["run_messages"]["severity"] == "fail" and err.rules["run_messages"]["message"] == "MinKNOW messages: 1 marked as an error."
    warn = build(V.make_with_message(tmp_path / "b", severity="WARNING", identifier="flow_cell_temperature", flow_cell_id="TST00073", seed=73), full=False)
    assert warn.rules["run_messages"]["severity"] == "warn" and warn.rules["run_messages"]["message"] == "MinKNOW messages: 1 warning."
    disk = build(make_run(tmp_path / "c", RunSpec(flow_cell_id="TST00074", seed=74, disk_warning=True)), full=False)
    r = disk.rules["run_messages"]                                  # a notice about the computer changes nothing
    assert r["outcome"] == "ok" and "1 notices about the computer's disk are in the run log" in r["message"]
    assert disk.bundle["modules"]["run"]["cells"][0]["scalars"]["messages_warning_computer"] == 1
    # the stored value is the count the sentence speaks of, and the metric names it
    assert (err.rules["run_messages"]["value"], err.rules["run_messages"]["metric"]) == (1, "run_messages_error")
    assert (warn.rules["run_messages"]["value"], warn.rules["run_messages"]["metric"]) == (1, "run_messages_warning")
    assert (r["value"], r["metric"]) == (0, "run_messages_error"), "nothing needs a look: the notice about the disk is not counted"


def test_output_files(tmp_path, build, quick8):
    assert quick8.rules["run_output_files"]["message"] == "Raw data written: every file reached the output folder."
    bad = build(V.set_files_in_fallback(make_run(tmp_path, RunSpec(flow_cell_id="TST00075", seed=75))), full=False)
    r = bad.rules["run_output_files"]
    assert (r["severity"], r["role"], r["value"]) == ("warn", "supporting", 3)
    assert r["message"] == "Raw data written: 3 files could not be written to the output folder."


def test_pores_and_channels_at_the_first_scan(tmp_path, build):
    low = build(make_run(tmp_path, RunSpec(flow_cell_id="TST00076", seed=76, first_scan_pore_frac=0.20)), full=False)
    r = low.rules["fc_pores_first_scan"]
    assert r["severity"] == "warn" and r["threshold"]["warn"] == 800 and r["threshold"]["fail"] is None
    assert r["message"] == catalogue.display("fc_pores_first_scan")["message_use_unknown"].replace("{value}", f"{r['value']:,}").replace("{limit}", "800")
    c = low.rules["fc_channels_with_pore_first_scan"]
    assert c["severity"] == "warn" and c["metric"] == "fc_channels_with_pore_first_scan_frac" and c["value"] < 0.65
    assert re.fullmatch(r"Channels with a working pore at the first scan: \d+ of 512 \(\d+(\.\d+)? %\)\. The limit is 65 % or more\.", c["message"])
    # a declared first use words the same count without "a used flow cell"
    st = judged(low.bundle, declared=Declared(flow_cell_use=1))
    assert "a used flow cell" not in st["_rules"]["fc_pores_first_scan"]["message"] and st["_file"]["relations"]["flow_cell_use"]["source"] == "declared"
    # the pore count of the flow cell check before loading replaces the first scan for the warranty figure
    st = judged(low.bundle, declared=Declared(flow_cell_check_pores=1200))
    assert st["_rules"]["fc_pores_first_scan"]["outcome"] == "ok" and "flow cell check before loading: 1,200" in st["_rules"]["fc_pores_first_scan"]["message"]


def test_occupancy(tmp_path, build):
    warn = build(make_run(tmp_path / "a", RunSpec(flow_cell_id="TST00077", seed=77, occupancy=0.55)), full=False)
    r = warn.rules["fc_occupancy_first_hour"]
    assert (r["severity"], r["threshold"]["warn"], r["threshold"]["fail"]) == ("warn", 0.7, 0.4)
    assert r["message"] == "Occupancy in the first hour: 55 %. ONT's guide is 70 % or more."
    fail = build(make_run(tmp_path / "b", RunSpec(flow_cell_id="TST00078", seed=78, occupancy=0.30)), full=False)
    assert fail.rules["fc_occupancy_first_hour"]["severity"] == "fail" and fail.status["verdict"] == "fail"


def test_temperature(tmp_path, build):
    f = V.set_temperature(make_run(tmp_path, RunSpec(flow_cell_id="TST00079", seed=79)), -1.4, sample=1)
    r = build(f, full=False).rules["fc_temperature_max_dev_c"]
    assert (r["severity"], r["role"]) == ("warn", "supporting") and r["message"] == "Temperature: up to 1.3 °C outside the target band. The limit is 0.5 °C or less."


# ------------------------------------------------------------- checks made on variations of one bundle
@pytest.mark.parametrize("check,summary,ok_value,bad_value,bad_level", [
    ("fc_recovering_frac", "fc__fc_recovering_frac", 0.07, 0.15, "warn"),
    ("fc_recovering_frac", "fc__fc_recovering_frac", 0.12, 0.21, "fail"),
    ("fc_occupancy_median", "fc__fc_occupancy_median", 0.70, 0.69, "warn"),
    ("fc_occupancy_median", "fc__fc_occupancy_median", 0.90, 0.30, "warn"),                 # supporting: a warning at most
    ("fc_recovering_last_third_frac", "fc__fc_recovering_last_third_frac", 0.25, 0.26, "warn"),
    ("fc_active_channels", "fc__fc_active_channels_frac", 0.65, 0.64, "warn"),
    ("dyn_signal_positive_frac", "dyn__dyn_signal_positive_frac", 0.90, 0.89, "warn"),
    ("dyn_signal_positive_frac", "dyn__dyn_signal_positive_frac", 0.95, 0.79, "fail"),
    ("dyn_unblock_frac", "dyn__dyn_unblock_frac", 0.08, 0.09, "warn"),
    ("dyn_speed_drop_frac", "dyn__dyn_speed_drop_frac", 0.15, 0.16, "warn"),
    ("len_est_n50", "len_est__len_est_n50", 500, 499, "warn"),
    ("dyn_speed_est_median", "dyn__dyn_speed_est_median", 332.5, 332.0, "warn"),             # the band 350 to 450, widened by 5 %
    ("dyn_speed_est_median", "dyn__dyn_speed_est_median", 472.5, 480.0, "warn"),
    ("dyn_speed_est_median", "dyn__dyn_speed_est_median", 400.0, 297.0, "fail"),             # more than 15 % below the band
    ("fc_temperature_max_dev_c", "fc__fc_temperature_max_dev_c", 0.5, 0.6, "warn"),
    ("fc_occupancy_first_hour", "fc__fc_occupancy_first_hour", 0.70, 0.6999, "warn"),
    ("fc_occupancy_first_hour", "fc__fc_occupancy_first_hour", 0.95, 0.39, "fail"),
    ("fc_channels_with_pore_first_scan", "fc__fc_channels_with_pore_first_scan_frac", 0.65, 0.6, "warn"),
])
def test_each_fixed_limit_on_both_sides(day_run, check, summary, ok_value, bad_value, bad_level):
    assert _state(judged(clone(day_run.bundle, **{summary: ok_value})), check) == "pass"
    st = judged(clone(day_run.bundle, **{summary: bad_value}))
    assert _state(st, check) == bad_level
    r = st["_rules"][check]
    assert r["severity"] == bad_level and (bad_level != "fail" or r["role"] == "key")
    assert not any(banned(t) for t in status_sentences(st))


def test_n50_limit_follows_the_assay(day_run):
    b = clone(day_run.bundle, len_est__len_est_n50=900)
    assert _state(judged(b), "len_est_n50") == "pass"                 # cDNA: 500
    g = copy.deepcopy(b)
    g["identity"]["assay"].update(type="gdna", group="gdna")
    r = judged(g)["_rules"]["len_est_n50"]                            # genomic DNA: 10,000
    assert r["severity"] == "warn" and r["message"] == "Estimated read N50: 900 bases. The limit is 10,000 bases or more."


def test_limits_come_from_the_catalogue():
    lim = Limits("promethion", "cdna")
    assert lim.warranty() == 5000 and Limits("minion", "cdna").warranty() == 800 and Limits("flongle", "gdna").warranty() == 50
    assert lim.fixed("fc_occupancy_first_hour") == {"kind": "min", "warn": 0.7, "fail": 0.4, "configured": False}
    assert lim.output_units() == [("reads", "key"), ("bases", "supporting")] and Limits("promethion", "gdna").output_units() == [("bases", "key")]
    assert lim.output_levels("reads", 72) == (50e6, 25e6) and lim.output_levels("reads", 24) == pytest.approx((22.5e6, 11.25e6))
    assert Limits("promethion", "drna").output_levels("reads", 12) == pytest.approx((8e6 * 0.35, 4e6 * 0.35))
    assert Limits("promethion", "gdna").output_levels("bases", 48) == pytest.approx((70e9 * 0.85, 40e9 * 0.85))
    assert Limits("promethion", "unknown").output_levels("bases", 72) is None
    for mid in ("fc_recovering_frac", "dyn_signal_positive_frac", "fc_occupancy_first_hour", "run_stop_cause", "run_messages"):
        assert lim.role(mid) == "key"
    for mid in ("fc_occupancy_median", "fc_active_channels", "dyn_unblock_frac", "len_est_n50", "run_output_files"):
        assert lim.role(mid) == "supporting"


def test_a_library_profile_changes_limits_roles_and_the_judged_unit(day_run):
    b = clone(day_run.bundle, fc__fc_recovering_frac=0.10, fc__fc_occupancy_median=0.3)
    base = judged(b)
    assert _state(base, "fc_recovering_frac") == "pass" and base["profile"]["source"] == "builtin"
    profile = {"cdna": {"rules": {"fc_recovering_frac": {"warn": 0.05, "fail": 0.08}, "fc_occupancy_median": {"enabled": False},
                                  "dyn_unblock_frac": {"role": "key", "warn": 0.001, "fail": 0.002}}}}
    st = judged(b, profile=profile)
    assert _state(st, "fc_recovering_frac") == "fail" and st["_rules"]["fc_recovering_frac"]["threshold"]["warn"] == 0.05
    assert "fc_occupancy_median" not in st["_rules"]
    assert _state(st, "dyn_unblock_frac") == "fail"                   # made a key check: it may fail now
    assert st["profile"]["source"] == "library_config" and st["profile"]["sha256"] != base["profile"]["sha256"]
    assert "custom_profile" in [m["code"] for m in st["modifiers"]]
    st = judged(b, profile={"cdna": {"outcome_unit": "bases"}})
    assert st["_rules"]["yield_estimated_bases_at_24h"]["role"] == "key" and "yield_reads_at_24h" not in st["_rules"]
    assert judged(b, profile={"gdna": {"rules": {"fc_recovering_frac": {"warn": 0.01}}}})["profile"]["source"] == "builtin"


# ------------------------------------------------------------------------------------ the verdict table
def test_under_one_hour_is_incomplete_stopped_early(tmp_path, build):
    built = build(V.make_short(tmp_path, minutes=6, flow_cell_id="TST00080", seed=80))
    st = built.status
    assert (st["verdict"], st["cause"], st["rules"]) == ("incomplete", "stopped_early", [])
    assert st["headline"] == "Stopped after 5 minutes. Too short to judge." and st["provisional"] is False
    assert problems(built.composed) == []


def test_one_to_six_hours_nothing_triggered_is_incomplete_output_not_judged(tmp_path, build):
    built = build(V.make_short(tmp_path, minutes=91, flow_cell_id="TST00081", seed=81, reads=600))
    st = built.status
    assert (st["verdict"], st["cause"]) == ("incomplete", "output_not_judged")
    assert st["headline"] == "Stopped after 1.5 hours. No check failed, but the run is too short to judge its output."
    assert not any(r["id"].startswith("yield_") for r in st["rules"]) and len(st["rules"]) >= 10
    assert problems(built.composed) == []


def test_one_to_six_hours_with_a_triggered_check_is_that_check(tmp_path, build):
    built = build(V.make_short(tmp_path, minutes=180, flow_cell_id="TST00082", seed=82, occupancy=0.2), full=False)
    st = built.status
    assert st["verdict"] == "fail" and st["cause"] is None and built.rules["fc_occupancy_first_hour"]["severity"] == "fail"
    assert st["headline"] == "Failed on occupancy in the first hour."


def test_six_hours_and_longer_is_judged_on_output_at_the_longest_horizon(day_run, full8):
    assert "yield_reads_at_6h" in full8.rules                           # an 8 h run: the 6 h horizon
    assert {"yield_reads_at_24h", "yield_estimated_bases_at_24h"} <= set(day_run.rules) and "yield_reads_at_6h" not in day_run.rules
    good = judged(clone(day_run.bundle, yield__yield_reads_at_24h=3_000_000, yield__yield_estimated_bases_at_24h=3_000_000_000))
    assert good["verdict"] == "pass" and good["by_tier"]["run_health"] == "pass" and good["cause"] is None
    assert good["headline"] == "The flow cell ran well (3.0 M reads by 24 h, 90 % of pores busy in the first hour)."
    r = good["_rules"]["yield_reads_at_24h"]
    assert r["threshold"] == {"kind": "min", "warn": 2250000, "fail": 1125000, "peers": 0, "applied": "fixed",
                              "basis": "fixed level for cdna on minion at 24 h: fewer than 5 comparable first-use runs started before this one"}
    assert r["message"] == "Reads at 24 h: 3.0 M. No lab history for this combination yet; the fixed level at 24 h is 2.2 M."
    assert [m["code"] for m in good["modifiers"]] == ["no_history"]
    warn = judged(clone(day_run.bundle, yield__yield_reads_at_24h=2_000_000, yield__yield_estimated_bases_at_24h=3_000_000_000))
    assert (warn["verdict"], _state(warn, "yield_reads_at_24h")) == ("warn", "warn")
    fail = judged(clone(day_run.bundle, yield__yield_reads_at_24h=1_000_000, yield__yield_estimated_bases_at_24h=100))
    assert _state(fail, "yield_reads_at_24h") == "fail" and _state(fail, "yield_estimated_bases_at_24h") == "warn"   # bases beside reads: a warning at most
    assert fail["verdict"] == "fail" and fail["headline"].startswith("Failed on reads at 24 h")


def test_genomic_dna_is_judged_on_estimated_bases(tmp_path, build):
    built = build(make_run(tmp_path, RunSpec(flow_cell_id="TST00083", seed=83, kit="SQK-NBD114-24")), full=False)
    assert built.bundle["identity"]["assay"]["group"] == "gdna"
    assert built.rules["yield_estimated_bases_at_6h"]["role"] == "key" and "yield_reads_at_6h" not in built.rules


def test_unknown_use_and_few_pores_cannot_fail_on_output(day_run):
    b = clone(day_run.bundle, yield__yield_reads_at_24h=1000, yield__yield_estimated_bases_at_24h=3_000_000_000, fc__fc_pores_first_scan=300)
    st = judged(b)
    assert _state(st, "yield_reads_at_24h") == "warn" and st["verdict"] == "warn"
    assert "flow_cell_use_unknown" in [m["code"] for m in st["modifiers"]]
    st = judged(b, declared=Declared(flow_cell_use=1))                   # first use is stated: the output check may fail
    assert _state(st, "yield_reads_at_24h") == "fail" and "flow_cell_use_unknown" not in [m["code"] for m in st["modifiers"]]


# ------------------------------------------------------------------------------------- quick detail
def test_quick_detail_is_provisional_with_pending_checks_and_says_so_far(day_run, quick8, tmp_path, build):
    st = quick8.status
    assert st["provisional"] is True and st["judged"]["frozen"] is False
    assert st["pending"] == [p for p in PENDING if p != "dyn_speed_drop_frac"]          # the fall in speed needs 24 h
    assert not (set(st["pending"]) & set(quick8.rules))
    assert st["detail"] == "Five checks wait for the read count. Read data: not available yet."
    assert "quick_detail" in [m["code"] for m in st["modifiers"]]
    quick_day = api.quick_bundle(day_run.run)
    good = judged(clone(quick_day, yield__yield_reads_at_24h=3_000_000, yield__yield_estimated_bases_at_24h=3_000_000_000))
    assert good["pending"] == list(PENDING) and good["detail"].startswith("Six checks wait for the read count.")
    assert good["verdict"] == "pass" and good["headline"] == "The flow cell ran well so far (3.0 M reads by 24 h, 90 % of pores busy in the first hour)."
    full = judged(clone(day_run.bundle, yield__yield_reads_at_24h=3_000_000, yield__yield_estimated_bases_at_24h=3_000_000_000))
    assert full["pending"] == [] and full["provisional"] is False and full["judged"]["frozen"] is True and "so far" not in full["headline"]
    card = api.card(quick_day, good["_file"])["columns"]
    assert card["provisional"] is True and card["pending"] == 6 and card["detail_t1"] == "quick" and card["speed"] is None


# ---------------------------------------------------------------------------------------- sentences
def test_sentences_of_checks_are_the_catalogues_copy(full8):
    cat = catalogue.metrics()
    for r in full8.status["rules"]:
        tpl = cat[r["id"]]["display"]["message"]
        head = tpl.split("{value}")[0]
        if r["outcome"] == "ok" and r["id"] in ("run_stop_cause", "run_output_files", "fc_temperature_max_dev_c"):
            continue                                                 # their all-clear sentences are in copy_pending
        assert r["message"].startswith(head), r["id"]
    lt = catalogue.limit_text()
    assert full8.rules["fc_recovering_frac"]["message"].endswith(lt["max"].replace("{limit}", "12 %"))
    assert full8.rules["dyn_speed_est_median"]["message"].endswith(lt["run_config"].replace("{low}", "350").replace("{high}", "450"))
    assert full8.rules["yield_reads_at_6h"]["message"].endswith(lt["history_none"].replace("{horizon}", "6 h").replace("{level}", "600,000"))


def test_no_stored_sentence_holds_a_word_the_interface_never_shows(full8, quick8, day_run):
    for built in (full8, quick8, day_run):
        for t in status_sentences(built.status):
            assert banned(t) == [], t
            assert "{" not in t and not re.search(r"\d%", t)
    for table in (CP.HEADLINE, CP.HEADLINE_FACT, CP.MODIFIER, CP.RULE, CP.REASON, CP.DETAIL):
        for text in table.values():
            assert banned(text) == [], text
    for text in catalogue.limit_text().values():
        assert banned(text.replace("{reason}", "x")) == []


def test_the_verdict_is_the_worst_triggered_check_and_only_key_checks_fail(day_run):
    b = clone(day_run.bundle, yield__yield_reads_at_24h=3_000_000, yield__yield_estimated_bases_at_24h=3_000_000_000,
              fc__fc_occupancy_median=0.2, dyn__dyn_unblock_frac=0.5)
    st = judged(b)
    assert st["verdict"] == "warn" == S.worst_of_rules(st["rules"], "run_health")
    assert all(r["severity"] != "fail" for r in st["rules"] if r["role"] != "key")
    assert st["headline"] == "Needs a look: occupancy over the run, reads ejected from blocked pores."
    assert S.reasons(st)[0]["check"] == "fc_occupancy_median"
    st = judged(clone(b, fc__fc_recovering_frac=0.3))
    # "Failed on" names only what failed; checks that only warned follow as "also needs a look"
    assert st["verdict"] == "fail"
    assert st["headline"] == "Failed on blocked pores. Also needs a look: occupancy over the run, reads ejected from blocked pores."
    failed = [r["id"] for r in st["rules"] if r["severity"] == "fail"]
    assert failed == ["fc_recovering_frac"]


def test_read_data_part_never_lowers_run_health(day_run, monkeypatch):
    b = clone(day_run.bundle, yield__yield_reads_at_24h=3_000_000, yield__yield_estimated_bases_at_24h=3_000_000_000)
    from flowcellar.readdata import api as readdata
    rule = {"id": "aln_mapped_frac", "module": "aln", "role": "key", "tier": "read_level", "severity": "fail", "outcome": "triggered",
            "metric": "aln_mapped_frac", "value": 0.2, "threshold": {"kind": "min", "warn": 0.8, "fail": 0.5, "applied": "fixed"},
            "message": "Reads that map: 20 %. The limit is 80 % or more."}
    monkeypatch.setattr(readdata, "read_checks", lambda bundle: [dict(rule)])
    st = judged(b)
    assert st["by_tier"] == {"run_health": "pass", "read_level": "fail"} and st["verdict"] == "fail"
    monkeypatch.setattr(readdata, "read_checks", lambda bundle: [dict(rule, outcome="not_evaluated", severity="none")])
    st = judged(b)
    assert st["by_tier"] == {"run_health": "pass", "read_level": "unknown"} and st["verdict"] == "pass"
    assert "None of the read-data checks could be assessed." in st["detail"]

    def broken(bundle):
        raise RuntimeError("read data is broken")
    monkeypatch.setattr(readdata, "read_checks", broken)
    assert judged(b)["by_tier"] == {"run_health": "pass", "read_level": "not_available"}


def test_a_failed_headline_never_lists_a_check_that_only_warned(day_run):
    """A supporting check can never fail: it must not stand under "Failed on" (Home, the Runs table and the scan
    output all show this sentence)."""
    b = clone(day_run.bundle, fc__fc_recovering_frac=0.3, len_est__len_est_n50=50, dyn__dyn_unblock_frac=0.5)
    st = judged(b)
    assert st["verdict"] == "fail"
    head = st["headline"]
    failed_part, _, also = head.partition(" Also needs a look: ")
    levels = {r["id"]: r["severity"] for r in st["rules"] if r["outcome"] == "triggered"}
    assert levels["fc_recovering_frac"] == "fail" and levels["len_est_n50"] == "warn" and levels["dyn_unblock_frac"] == "warn"
    assert failed_part.startswith("Failed on blocked pores")
    assert "estimated read N50" not in failed_part and "reads ejected" not in failed_part, head
    assert len([i for i, v in levels.items() if v == "fail"]) == failed_part.count(",") + 1, "every failed check, and only those"
    assert "estimated read N50" in also and "reads ejected from blocked pores" in also and len(head) <= 140
    # only warnings: the old sentence, unchanged
    st = judged(clone(day_run.bundle, yield__yield_reads_at_24h=3_000_000, yield__yield_estimated_bases_at_24h=3_000_000_000,
                      len_est__len_est_n50=50))
    assert st["verdict"] == "warn" and st["headline"].startswith("Needs a look: ") and "Failed" not in st["headline"]


def test_read_data_is_a_verdict_so_far_while_files_of_the_run_wait(day_run):
    """One check assessed from a small table while the BAMs wait must not read as a plain Pass."""
    b = clone(day_run.bundle)
    b["tiers"]["read_level"] = {"state": "partial", "samples_expected": 2, "samples_with_data": 2}
    b["provenance"]["inputs"] = list(b["provenance"]["inputs"]) + [
        {"id": "rd1", "role": "align.genome_bam", "path": "/x/a.bam", "state": "skipped", "notes": [{"level": "info", "code": "waiting", "message": "Not read yet."}]},
        {"id": "rd2", "role": "align.genome_bam", "path": "/x/b.bam", "state": "skipped", "notes": [{"level": "info", "code": "waiting", "message": "Not read yet."}]}]
    import flowcellar.readdata.api as rd_api
    rules = [{"id": "funnel_usable_frac", "module": "funnel", "role": "key", "tier": "read_level", "severity": "none", "outcome": "ok",
              "metric": "funnel_usable_frac", "value": 0.7, "threshold": {"kind": "min", "warn": 0.5, "fail": None, "basis": "built-in limit", "applied": "fixed"},
              "message": "Usable share: 70 %.", "panel": "funnel"},
             {"id": "len_n50", "module": "len", "role": "supporting", "tier": "read_level", "severity": "none", "outcome": "not_evaluated",
              "metric": "len_n50", "value": None, "threshold": {"kind": "rule", "basis": "waiting", "applied": "rule"},
              "message": "Read N50: not assessed, not every read has been read.", "not_evaluated_reason": "waiting", "panel": "len"}]
    import pytest as _pytest
    mp = _pytest.MonkeyPatch()
    try:
        mp.setattr(rd_api, "read_checks", lambda bundle: [dict(x) for x in rules])
        st = judged(b)
        mods = {m["code"]: m["message"] for m in st["modifiers"]}
        assert st["by_tier"]["read_level"] == "pass"
        assert mods["read_data_so_far"] == ("Read data so far: 2 files of this run have not been read yet; 1 of 2 read-data checks "
                                            "could be assessed.")
        b["provenance"]["inputs"] = b["provenance"]["inputs"][:-2]
        assert "read_data_so_far" not in [m["code"] for m in judged(b)["modifiers"]], "nothing waits: a plain verdict"
    finally:
        mp.undo()


def test_a_range_keeps_the_unit_of_its_lower_end_when_the_two_ends_differ_in_unit():
    # "Typical here is 727 Mb to 1.1 Gb", never "727b to 1.1 Gb"; one unit is written once, after the higher end
    from flowcellar.engine.status import _low_of_range
    from flowcellar.util.text import fmt_bases, fmt_count
    assert _low_of_range(fmt_bases(727e6), fmt_bases(1.1e9)) == "727 Mb"
    assert _low_of_range(fmt_bases(52.4e9), fmt_bases(142e9)) == "52.4"
    assert _low_of_range(fmt_bases(53e6), fmt_bases(84e6)) == "53"
    assert _low_of_range(fmt_count(63.2e6), fmt_count(140e6)) == "63.2"
    assert _low_of_range(fmt_count(525_747), fmt_count(831_020)) == "525,747"
    assert _low_of_range(fmt_count(950_000), fmt_count(1.2e6)) == "950,000"
