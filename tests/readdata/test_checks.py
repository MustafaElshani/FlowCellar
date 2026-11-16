"""The read-data checks: limits from the catalogue, evaluated only on a population they fit, per sample, in plain words."""
from __future__ import annotations

import copy

import pytest

from flowcellar.aggregates import build
from flowcellar.model import axes
from flowcellar.readdata import api, checks
from flowcellar.readdata import copy_pending as words

POP = {"stage": "basecalled_pass", "min_qscore": 7.0, "trimmed": False, "both_end_barcodes": None, "source_role": "dorado.bam", "evidence": ["test"]}


def len_cell(key, flt, lengths):
    import numpy as np
    arr = np.asarray(lengths, dtype=np.float64)
    idx = axes.bin_index(arr, "read_length_log40")
    n = np.bincount(idx, minlength=axes.n_cells("read_length_log40"))
    w = np.rint(np.bincount(idx, weights=arr, minlength=axes.n_cells("read_length_log40"))).astype(int)
    return {"at": [key, flt], "scalars": {"reads": len(lengths), "bases": int(arr.sum())}, "hists": {"len_hist": build.hist("read_length_log40", n, w, w_name="bases")}}


def qual_cell(key, flt, quals):
    import numpy as np
    n = np.bincount(axes.bin_index(np.asarray(quals, float), "qscore_0p25"), minlength=axes.n_cells("qscore_0p25"))
    return {"at": [key, flt], "scalars": {"reads": len(quals)}, "hists": {"qual_hist": build.hist("qscore_0p25", n, missing=0)}}


def section(cells, **kw):
    sec = {"version": "0.1.0", "tier": "read_level", "state": "ok", "dims": ["sample", "filter"], "cells": cells, "notes": [], "population": dict(POP)}
    sec.update(kw)
    return sec


def bundle(modules, *, group="cdna", samples=("one", "two", "three"), basecalling=None, family="promethion"):
    return {"identity": {"run_key": "TST00001_0a1b2c3d", "assay": {"group": group, "chemistry": "r10.4.1"}, "flow_cell": {"family": family},
                         "basecalling": basecalling if basecalling is not None else [{"model": "dna_r10.4.1_e8.2_400bps_sup@v5.2.0", "model_tier": "sup", "primary": True}]},
            "samples": [{"key": f"s{i + 1}", "name": n} for i, n in enumerate(samples)], "modules": modules}


def by_id(rules):
    return {r["id"]: r for r in rules}


def test_per_sample_check_is_one_rule_and_the_worst_sample_decides():
    cells = [len_cell("s1", "pass", [1500] * 60 + [900] * 40), len_cell("s2", "pass", [600] * 100), len_cell("s3", "pass", [1400] * 100)]
    r = by_id(api.read_checks(bundle({"len": section(cells)})))["len_n50"]
    # cDNA: warn below 700 bases, fail below 450. One sample is at about 600.
    assert r["outcome"] == "triggered" and r["severity"] == "warn" and r["role"] == "key" and r["tier"] == "read_level" and r["module"] == "len"
    assert r["scope"] == {"sample": "s2"} and 560 < r["value"] < 640 and r["threshold"]["warn"] == 700 and r["threshold"]["fail"] == 450
    assert r["threshold"]["kind"] == "min" and r["threshold"]["applied"] == "fixed" and "cDNA" in r["threshold"]["basis"]
    assert r["message"].startswith("Read N50: ") and "in sample two" in r["message"] and r["message"].endswith("The limit is 700 bases or more.")
    assert r["population"] == POP
    # two samples off: both are named, the rule names no single sample
    cells[2] = len_cell("s3", "pass", [300] * 100)
    r = by_id(api.read_checks(bundle({"len": section(cells)})))["len_n50"]
    assert r["severity"] == "fail" and r["scope"] == {"sample": None}
    assert "below the limit in 2 of 3 samples (three and two)" in r["message"] and "the lowest is" in r["message"]
    # all fine: the value over all samples and the lowest one
    ok = [len_cell("s1", "pass", [1500] * 100), len_cell("s2", "pass", [1200] * 100)]
    r = by_id(api.read_checks(bundle({"len": section(ok)})))["len_n50"]
    assert r["outcome"] == "ok" and r["severity"] == "none" and "over all samples; lowest sample" in r["message"]


def test_reads_without_a_barcode_and_failed_reads_do_not_decide_a_sample_check():
    cells = [len_cell("s1", "pass", [1500] * 100), len_cell("s1", "fail", [100] * 500), len_cell("unclassified", "pass", [90] * 900)]
    r = by_id(api.read_checks(bundle({"len": section(cells)}, samples=("one",))))
    assert r["len_n50"]["outcome"] == "ok" and r["len_short_read_frac"]["outcome"] == "ok" and r["len_short_read_frac"]["value"] == 0.0


def test_supporting_check_can_warn_but_not_fail():
    cells = [len_cell("s1", "pass", [100] * 90 + [1500] * 10)]              # 90 % of reads under 200 bases: beyond the fail level
    r = by_id(api.read_checks(bundle({"len": section(cells)}, samples=("one",))))["len_short_read_frac"]
    assert r["role"] == "supporting" and r["outcome"] == "triggered" and r["severity"] == "warn" and r["value"] == 0.9
    assert r["message"] == "Reads under 200 bases: 90.0 %. The limit is 5.0 % or less."


def test_limits_follow_chemistry_and_basecalling_model():
    cells = [qual_cell("s1", "pass", [14.5] * 100)]
    sup = by_id(api.read_checks(bundle({"qual": section(cells)}, samples=("one",))))["qual_median"]
    assert sup["outcome"] == "triggered" and sup["threshold"]["warn"] == 16 and "sup basecalling" in sup["threshold"]["basis"]
    hac = bundle({"qual": section(cells)}, samples=("one",), basecalling=[{"model": "dna_r10.4.1_e8.2_400bps_hac@v5.0.0", "model_tier": "hac", "primary": True}])
    r = by_id(api.read_checks(hac))["qual_median"]
    assert r["outcome"] == "ok" and r["threshold"]["warn"] == 13 and r["message"] == "Median read quality: Q 14.6. The limit is Q 13.0 or more."
    # no basecalling known and no chemistry: no clause matches, so the check is not assessed, never guessed
    none = bundle({"qual": section(cells)}, samples=("one",), basecalling=[])
    none["identity"]["assay"]["chemistry"] = None
    r = by_id(api.read_checks(none))["qual_median"]
    assert r["outcome"] == "not_evaluated" and r["not_evaluated_reason"] == words.REASONS["no_limit"] and r["severity"] == "none"
    assert r["message"] == "Not assessed: no limit is defined for this kind of run."


def test_check_is_not_assessed_on_a_population_it_does_not_fit():
    yb = section([{"at": ["s1", "pass"], "scalars": {"reads": 900, "bases": 900_000}}], summary={"yield_pass_frac": None, "yield_basecalled_reads": 900})
    b = bundle({"yield.basecalled": yb, "run": {"summary": {"run_completion_frac": 1.0}}}, samples=("one",))
    r = by_id(api.read_checks(b))
    assert r["yield_pass_frac"]["outcome"] == "not_evaluated" and r["yield_pass_frac"]["population"]["stage"] == "basecalled_pass"
    assert r["yield_pass_frac"]["message"] == ("Not assessed: it is defined for every read the basecaller wrote, and the reads on this machine are "
                                               "reads that passed the quality filter.")
    assert r["yield_basecalled_bases"]["outcome"] == "triggered"            # this one fits reads that passed
    full = copy.deepcopy(b)
    full["modules"]["yield.basecalled"]["population"]["stage"] = "full_length_trimmed"
    r2 = by_id(api.read_checks(full))
    assert r2["yield_basecalled_bases"]["outcome"] == "not_evaluated" and "full-length reads, trimmed" in r2["yield_basecalled_bases"]["message"]
    # every read the basecaller wrote: the share that passed is judged
    allr = copy.deepcopy(b)
    allr["modules"]["yield.basecalled"]["population"]["stage"] = "basecalled_all"
    allr["modules"]["yield.basecalled"]["summary"]["yield_pass_frac"] = 0.65
    r3 = by_id(api.read_checks(allr))["yield_pass_frac"]
    assert r3["outcome"] == "triggered" and r3["severity"] == "fail" and r3["value"] == 0.65
    assert r3["message"] == "Reads that passed the basecaller's quality filter: 65.0 %. The limit is 85.0 % or more."


def test_output_of_a_run_that_stopped_early_is_left_to_run_health():
    yb = section([{"at": ["s1", "pass"], "scalars": {"reads": 900, "bases": 900_000}}], summary={})
    for run_summary, expect in (({"run_completion_frac": 0.3, "run_duration_h": 5.0}, "not_evaluated"), ({"run_completion_frac": 0.3, "run_duration_h": 50.0}, "triggered"),
                                ({}, "not_evaluated")):
        r = by_id(api.read_checks(bundle({"yield.basecalled": yb, "run": {"summary": run_summary}}, samples=("one",))))["yield_basecalled_bases"]
        assert r["outcome"] == expect
        if expect == "not_evaluated":
            assert "planned length" in r["message"]


def test_cdna_full_length_share_is_not_judged_on_reads_with_a_barcode_at_both_ends():
    cd = {"version": "0.1.0", "tier": "assay", "state": "ok", "dims": ["sample"], "notes": [], "summary": {"cdna_full_length_frac": 0.978},
          "cells": [{"at": ["s1"], "scalars": {"pass_reads": 1000, "primers_found": 970, "rescued_reads": 8, "unusable": 22}},
                    {"at": ["s2"], "scalars": {"pass_reads": 1000, "primers_found": 500, "rescued_reads": 60, "unusable": 440}}]}
    r = by_id(api.read_checks(bundle({"cdna": cd}, samples=("one", "two"))))
    assert r["cdna_full_length_frac"]["outcome"] == "triggered" and r["cdna_full_length_frac"]["severity"] == "fail" and r["cdna_full_length_frac"]["scope"] == {"sample": "s2"}
    assert r["cdna_full_length_frac"]["value"] == 0.56 and "56.0 % in sample two" in r["cdna_full_length_frac"]["message"]
    assert r["cdna_rescued_frac"]["outcome"] == "triggered" and r["cdna_rescued_frac"]["severity"] == "warn" and r["cdna_rescued_frac"]["value"] == 0.06
    both = bundle({"cdna": cd}, samples=("one", "two"), basecalling=[{"model": "x_sup@v5", "model_tier": "sup", "primary": True, "barcode_both_ends": True}])
    r = by_id(api.read_checks(both))["cdna_full_length_frac"]
    assert r["outcome"] == "not_evaluated" and "barcode at both ends" in r["message"] and "97.8 %" in r["message"]
    assert by_id(api.read_checks(bundle({"cdna": cd}, group="gdna")))["cdna_full_length_frac"]["outcome"] == "not_evaluated"      # no limit outside cDNA


def test_read_fate_is_judged_only_when_every_stage_it_needs_is_known():
    fu = {"version": "0.1.0", "tier": "read_level", "state": "ok", "dims": [], "cells": [], "notes": [],
          "summary": {"funnel_usable_frac": 0.31, "funnel_usable_stage": "quantified", "funnel_basecall_loss_frac": 0.12}, "params": {"usable_not_assessed_reason": None}}
    r = by_id(api.read_checks(bundle({"funnel": fu})))
    assert r["funnel_usable_frac"]["outcome"] == "triggered" and r["funnel_usable_frac"]["severity"] == "fail" and r["funnel_usable_frac"]["role"] == "key"
    assert r["funnel_usable_frac"]["message"] == "Reads that ended up counted towards a transcript: 31.0 % of the reads sequenced. The limit is 55.0 % or more."
    assert r["funnel_basecall_loss_frac"]["outcome"] == "ok" and r["funnel_basecall_loss_frac"]["role"] == "supporting"
    fu["summary"].update(funnel_usable_frac=None, funnel_basecall_loss_frac=None)
    fu["params"]["usable_not_assessed_reason"] = "the reads counted towards a transcript are not known for every sample of this flow cell"
    r = by_id(api.read_checks(bundle({"funnel": fu})))
    assert r["funnel_usable_frac"]["outcome"] == "not_evaluated" and "not known for every sample" in r["funnel_usable_frac"]["message"]
    assert "funnel_basecall_loss_frac" not in r


@pytest.mark.parametrize("change, reason", [({"state": "partial", "partial_reason": "partly_read"}, "not_counted_yet"),
                                            ({"state": "partial", "partial_reason": "waiting"}, "not_counted_yet"),
                                            ({"notes": [{"level": "warning", "code": "weak_link", "message": "x"}]}, "weak_link"),
                                            ({"notes": [{"level": "warning", "code": "mixed_sources", "message": "x"}]}, "mixed_sources"),
                                            ({"cells": [], "notes": [{"level": "warning", "code": "two_candidate_sources", "message": "x"}]}, "two_sources")])
def test_nothing_is_judged_on_numbers_that_are_not_totals_or_not_trusted(change, reason):
    sec = section([len_cell("s1", "pass", [300] * 100)])
    sec.update(change)
    r = by_id(api.read_checks(bundle({"len": sec}, samples=("one",))))["len_n50"]
    assert r["outcome"] == "not_evaluated" and r["not_evaluated_reason"] == words.REASONS[reason] and r["value"] is None
    # samples missing is different: the samples that are there are judged
    sec = section([len_cell("s1", "pass", [300] * 100)], state="partial", partial_reason="samples_missing")
    assert by_id(api.read_checks(bundle({"len": sec}, samples=("one",))))["len_n50"]["outcome"] == "triggered"


def test_no_stored_sentence_uses_an_internal_word(project_data):
    banned = ("tier", "bundle", "partial", "population", "module", "fragment", "tier 2", "source of record", "aggregate", "sampled")
    sentences = [r["message"] for r in api.read_checks(project_data["bundle"])]
    sentences += [t for t in words.REASONS.values()] + [t for t in words.SENTENCES.values()]
    for s in sentences:
        low = s.lower()
        assert not any(b in low for b in banned), s
    assert all(m in checks.copy.CHECK_COPY for m in ("len_n50", "qual_median", "aln_identity_median", "funnel_usable_frac"))


def test_checks_of_a_bundle_without_read_data_are_none():
    assert api.read_checks({"identity": {}, "modules": {}}) == []
    assert api.read_checks({}) == []


def test_a_library_can_change_or_switch_off_the_limit_of_a_read_data_check():
    cells = [len_cell("s1", "pass", [1500] * 60 + [900] * 40), len_cell("s2", "pass", [600] * 100), len_cell("s3", "pass", [1400] * 100)]
    b = bundle({"len": section(cells)})
    # the built-in limit warns below 700 bases: one sample is at about 600
    assert by_id(api.read_checks(b))["len_n50"]["severity"] == "warn"
    # [profiles.cdna.rules.len_n50] warn = 500: the same reads pass, and the rule says where its limit comes from
    r = by_id(api.read_checks(b, {"cdna": {"rules": {"len_n50": {"warn": 500, "fail": 300}}}}))["len_n50"]
    assert r["outcome"] == "ok" and r["severity"] == "none" and r["threshold"]["warn"] == 500 and r["threshold"]["fail"] == 300
    assert "flowcellar.toml" in r["threshold"]["basis"] and r["message"].endswith("The limit is 500 bases or more.")
    # one level given: the other stays as built in
    r = by_id(api.read_checks(b, {"cdna": {"rules": {"len_n50": {"warn": 650}}}}))["len_n50"]
    assert r["threshold"]["warn"] == 650 and r["threshold"]["fail"] == 450 and r["severity"] == "warn"
    # a table for another assay group changes nothing; neither does a value that is not a number
    for other in ({"gdna": {"rules": {"len_n50": {"warn": 100}}}}, {"cdna": {"rules": {"len_n50": {"warn": "low"}}}}, {"cdna": {}}, {}):
        r = by_id(api.read_checks(b, other))["len_n50"]
        assert r["threshold"]["warn"] == 700 and "built-in" in r["threshold"]["basis"], other
    # enabled = false: stored as not assessed, with the reason
    r = by_id(api.read_checks(b, {"cdna": {"rules": {"len_n50": {"enabled": False}}}}))["len_n50"]
    assert r["outcome"] == "not_evaluated" and r["severity"] == "none" and "switched off" in r["message"]


def test_the_judge_hands_the_limits_of_the_library_to_the_read_data_checks():
    from flowcellar.engine import api as engine_api
    cells = [len_cell("s1", "pass", [600] * 100)]
    b = bundle({"len": section(cells)}, samples=("one",))
    b.update(bundle_id="TST00001_0a1b2c3d", tiers={"read_level": {"state": "complete"}}, provenance={"inputs": []})
    plain = {r["id"]: r for r in engine_api.judge(b, [], in_library=False)["status"]["rules"]}
    changed = {r["id"]: r for r in engine_api.judge(b, [], profile={"cdna": {"rules": {"len_n50": {"warn": 500}}}}, in_library=False)["status"]["rules"]}
    assert plain["len_n50"]["severity"] == "warn" and changed["len_n50"]["severity"] == "none" and changed["len_n50"]["threshold"]["warn"] == 500
