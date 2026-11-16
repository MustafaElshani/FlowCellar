"""The verdict of an analysis: each check on a number that triggers it and on one that does not."""
import copy

import pytest

from flowcellar.model import schema
from flowcellar.sc import catalogue, judge, store

from .conftest import composed, summary


def _with(bundle, **values):
    b = copy.deepcopy(bundle)
    summary(b).update(values)
    return b


def _rule(status, check_id):
    return next(r for r in status["rules"] if r["id"] == check_id)


def test_the_pipelines_alarm_is_the_first_reason_word_for_word(sr_visium):
    _, b = sr_visium
    st = judge.judge(b)
    assert st["verdict"] == "warn"
    r = _rule(st, "sp_pipeline_alarms")
    assert r["outcome"] == "triggered" and r["severity"] == "warn"
    assert r["message"] == ('Space Ranger alarm, in its own words: "High Fraction of Reads Mapped Antisense to Genes" (20.0%). '
                            "Ideal < 10%. This can indicate use of an unsupported chemistry type.")
    assert judge.reasons(st)[0]["check"] == "sp_pipeline_alarms"
    assert st["headline"].startswith("Space Ranger raised one alarm.")
    assert st["detail"] == "This verdict is about the analysis. It does not change the verdict of any flow cell."
    assert schema.schema_errors(composed(b)) == []


def test_an_error_alarm_fails(sr_visium):
    _, b = sr_visium
    b = copy.deepcopy(b)
    b["modules"]["sc"]["notes"].append({"level": "error", "code": "pipeline_alarm", "title": "No Tissue Detected", "message": "Check the image.",
                                        "source": "spaceranger"})
    st = judge.judge(b)
    assert st["verdict"] == "fail" and any(r["severity"] == "fail" and "Space Ranger error" in r["message"] for r in st["rules"])


@pytest.mark.parametrize("check, key, ok, warn, fail", [
    ("sr_antisense_fraction", "sr_antisense", 0.05, 0.2, 0.35),
    ("sr_usable_reads_fraction", "sr_usable_reads_fraction", 0.5, 0.2, 0.05),
    ("sr_valid_barcodes", "sr_valid_barcodes", 0.9, 0.7, 0.4),
    ("sp_fraction_reads_under_tissue", "spatial_fraction_reads_under_tissue", 0.9, 0.45, 0.2),
])
def test_key_checks_of_space_ranger(sr_hd, check, key, ok, warn, fail):
    _, b = sr_hd
    for value, outcome, sev in ((ok, "ok", "none"), (warn, "triggered", "warn"), (fail, "triggered", "fail")):
        st = judge.judge(_with(b, **{key: value}))
        r = _rule(st, check)
        assert (r["outcome"], r["severity"], r["role"]) == (outcome, sev, "key"), (check, value)
        assert r["value"] == value and r["threshold"]["warn"] == catalogue.limits(check)["warn"]
        if sev == "fail":
            assert st["verdict"] == "fail"


@pytest.mark.parametrize("check, key, ok, bad", [
    ("sr_valid_umis", "sr_valid_umis", 0.99, 0.6), ("sr_mapping", "sr_mapped_transcriptome", 0.8, 0.2),
    ("sp_segmentation_summary", "spatial_umis_in_cells", 0.64, 0.3), ("sc_mito_fraction_per_barcode", "mito_fraction_median", 0.02, 0.4),
])
def test_supporting_checks_warn_and_never_fail(sr_hd, check, key, ok, bad):
    _, b = sr_hd
    assert _rule(judge.judge(_with(b, **{key: ok})), check)["outcome"] == "ok"
    r = _rule(judge.judge(_with(b, **{key: bad})), check)
    assert (r["outcome"], r["severity"], r["role"]) == ("triggered", "warn", "supporting")


@pytest.mark.parametrize("check, key, ok, warn, fail", [
    ("sc_full_length_fraction", "reads_full_length_fraction", 0.95, 0.7, 0.5),
    ("sc_valid_barcode_fraction", "reads_valid_barcode_fraction", 0.85, 0.6, 0.4),
    ("sc_gene_assigned_fraction", "reads_gene_assigned_fraction", 0.8, 0.5, 0.3),
])
def test_key_checks_of_the_workflow(wf_v2, check, key, ok, warn, fail):
    _, b = wf_v2
    for value, sev in ((ok, "none"), (warn, "warn"), (fail, "fail")):
        r = _rule(judge.judge(_with(b, **{key: value})), check)
        assert r["severity"] == sev, (check, value)
    assert "matched to a spot" in _rule(judge.judge(b), "sc_valid_barcode_fraction")["message"]


def test_other_workflow_checks(wf_v2, wf_v3):
    _, b = wf_v2
    assert _rule(judge.judge(_with(b, reads_plus_fraction=0.9)), "sc_strand_plus_fraction")["severity"] == "warn"
    assert _rule(judge.judge(_with(b, reads_plus_fraction=0.5)), "sc_strand_plus_fraction")["threshold"]["pass_band"] == [0.35, 0.65]
    assert _rule(judge.judge(_with(b, reads_valid_umi_fraction=0.5)), "sc_valid_umi_fraction")["severity"] == "warn"
    assert _rule(judge.judge(_with(b, reads_transcript_assigned_fraction=0.1)), "sc_transcript_assigned_fraction")["severity"] == "warn"
    assert _rule(judge.judge(_with(b, reads_chrm_fraction=0.31)), "sc_reads_chrm_fraction")["message"] == \
        "31 % of the aligned reads map to the mitochondrial genome. Warn above 25 %."
    assert _rule(judge.judge(_with(b, workflow_failed_tasks=2)), "sc_wf_execution")["message"] == "2 of the 5 tasks of the workflow failed."
    # the 8 um bin takes "an"
    _, hd = wf_v3
    assert "matched to an 8 µm bin" in _rule(judge.judge(hd), "sc_valid_barcode_fraction")["message"]


def test_a_number_without_a_source_is_not_assessed_never_a_failure(wf_v3):
    _, b = wf_v3
    b = copy.deepcopy(b)
    for k in ("reads_full_length_fraction", "reads_valid_barcode_fraction", "reads_gene_assigned_fraction"):
        summary(b).pop(k)
    st = judge.judge(b)
    assert _rule(st, "sc_full_length_fraction")["outcome"] == "not_evaluated" and _rule(st, "sc_full_length_fraction")["severity"] == "none"
    assert st["verdict"] == "pass"


def test_spots_under_tissue(sr_visium, sr_hd):
    _, b = sr_visium
    assert _rule(judge.judge(_with(b, spatial_units_under_tissue=60)), "sp_units_under_tissue")["severity"] == "warn"
    assert _rule(judge.judge(_with(b, spatial_units_under_tissue=0)), "sp_units_under_tissue")["severity"] == "fail"
    assert _rule(judge.judge(_with(b, spatial_units_under_tissue=300)), "sp_units_under_tissue")["message"] == "300 spots lie under tissue."
    _, hd = sr_hd                        # the limit of 100 is about spots; a Visium HD slide with few bins only fails at none
    assert _rule(judge.judge(_with(hd, spatial_units_under_tissue=60)), "sp_units_under_tissue")["outcome"] == "ok"


def test_the_registration_finding_is_one_reason_not_two(sr_hd):
    _, b = sr_hd
    st = judge.judge(b)
    reg = _rule(st, "sp_registration")
    assert reg["outcome"] == "triggered" and reg["severity"] == "info" and reg["role"] == "info"     # the pipeline's alarm already says it
    assert [r["check"] for r in judge.reasons(st)] == ["sp_pipeline_alarms"] and st["verdict"] == "warn"
    quiet = copy.deepcopy(b)
    quiet["modules"]["sc"]["notes"] = []
    assert _rule(judge.judge(quiet), "sp_registration")["severity"] == "warn"


def test_cells_against_the_expected_number(wf_droplet):
    _, b = wf_droplet
    assert _rule(judge.judge(_with(b, barcodes_cells=3, expected_cells=12)), "sc_cells")["severity"] == "warn"
    assert _rule(judge.judge(_with(b, barcodes_cells=11, expected_cells=12)), "sc_cells")["outcome"] == "ok"


def test_links_never_lower_the_verdict(wf_v3):
    _, b = wf_v3
    st = judge.judge(b, in_library=True)                 # both runs are named and neither is in the library
    r = _rule(st, "sc_link_runs")
    assert (r["outcome"], r["severity"], r["role"]) == ("triggered", "info", "info") and "TST00001, TST00002" in r["message"]
    assert st["verdict"] == "pass" and judge.reasons(st) == []
    assert _rule(judge.judge(b, in_library=False), "sc_link_runs")["outcome"] == "not_evaluated"
    none = copy.deepcopy(b)
    none["identity"]["analysis"]["links"] = []
    assert "No flow cell run could be named" in _rule(judge.judge(none, in_library=True), "sc_link_runs")["message"]


def test_two_pipelines_counting_the_same_reads(wf_v3):
    _, b = wf_v3
    mol = summary(b)["barcodes_unique_molecules"]
    ok = _rule(judge.judge(b, sibling={"barcodes_unique_molecules": int(mol / 1.3)}), "sp_wf_sr_umi_ratio")
    assert ok["outcome"] == "ok" and ok["value"] == pytest.approx(1.3, abs=0.01)
    assert _rule(judge.judge(b, sibling={"barcodes_unique_molecules": int(mol / 2.5)}), "sp_wf_sr_umi_ratio")["severity"] == "warn"
    assert not any(r["id"] == "sp_wf_sr_umi_ratio" for r in judge.judge(b)["rules"])        # no sibling: nothing to compare


def test_a_library_can_move_or_switch_off_a_limit(sr_hd):
    _, b = sr_hd
    b = _with(b, sr_antisense=0.2)
    assert _rule(judge.judge(b), "sr_antisense_fraction")["severity"] == "warn"
    moved = judge.judge(b, profiles={"visium_hd": {"rules": {"sr_antisense_fraction": {"warn": 0.25}}}})
    assert _rule(moved, "sr_antisense_fraction")["outcome"] == "ok" and moved["profile"]["source"] == "library_config"
    off = judge.judge(b, profiles={"visium_hd": {"rules": {"sr_antisense_fraction": {"enabled": False}}}})
    assert not any(r["id"] == "sr_antisense_fraction" for r in off["rules"])
    other_assay = judge.judge(b, profiles={"visium": {"rules": {"sr_antisense_fraction": {"warn": 0.25}}}})
    assert _rule(other_assay, "sr_antisense_fraction")["severity"] == "warn"


def test_headlines(sr_hd, wf_v3):
    _, b = wf_v3
    st = judge.judge(b)
    assert st["verdict"] == "pass" and st["headline"].startswith("All ") and st["headline"].endswith("checks passed.")
    _, hd = sr_hd
    assert judge.judge(hd)["headline"] == "Space Ranger raised one alarm. Every other check is within its limit."
    two = judge.judge(_with(hd, sr_antisense=0.2, sr_valid_umis=0.5))
    assert two["headline"] == "Space Ranger raised one alarm. 2 more checks need a look."


def test_every_check_of_the_catalogue_has_its_sentences_and_a_metric():
    ids = set(catalogue.metrics())
    for c in catalogue.checks():
        assert c["id"] in ids and c["ok"] and c["role"] in ("key", "supporting", "info") and c["kind"] in ("min", "max", "range", "rule")
        if c["kind"] in ("min", "max") and c["role"] == "key":
            assert c.get("fail"), c["id"]
    # the one check that is shown only: a run that is not in the library never lowers the verdict of an analysis
    assert [c["id"] for c in catalogue.checks() if c["role"] == "info"] == ["sc_link_runs"]
    assert len(catalogue.metrics()) == 57


def test_card_row_is_small_and_holds_no_map(sr_hd):
    import json
    _, b = sr_hd
    st = judge.judge(b, in_library=True)
    card = store.card_for(b, {"status": st})
    row = card["row"]
    assert len(json.dumps(row)) < 2500 and "spatial_map" not in json.dumps(row)
    assert row["verdict"] == "warn" and row["alarms"] == 1 and row["assay"] == "visium_hd" and row["unit"] == "bin"
    assert row["spatial_units_under_tissue"] == summary(b)["spatial_units_under_tissue"]
