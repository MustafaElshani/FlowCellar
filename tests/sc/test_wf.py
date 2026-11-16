"""wf-single-cell output, both layouts: exact numbers, and the few that exist only in the workflow's report."""
import os

import numpy as np
import pytest

from flowcellar.library.maintain import check_bundle
from flowcellar.sc import axes, wfsc
from flowcellar.testing import sc as fx

from .conftest import cell, collected, composed, summary


def _check_summary(fake, b):
    s = summary(b)
    for key, want in fake.expected["summary"].items():
        got = s.get(key)
        assert got == (pytest.approx(want, rel=1e-12, abs=1e-12) if isinstance(want, float) else want), key


def test_v2_numbers_are_exact(wf_v2):
    fake, b = wf_v2
    _check_summary(fake, b)
    ana = b["identity"]["analysis"]
    assert (ana["label"], ana["assay"], ana["pipeline_version"], ana["reference"]) == ("demo_spatial_a", "visium", "2.2.0", "demo-2024-A")
    c = cell(b)
    adapter = dict(zip(c["counts"]["reads_adapter_config"]["cats"], c["counts"]["reads_adapter_config"]["n"]))
    assert adapter == fake.expected["adapter"]
    s = summary(b)
    assert c["counts"]["reads_funnel"]["n"] == [s["reads_input"], s["reads_full_length"], s["reads_valid_barcode"], s["reads_gene_assigned"],
                                                s["reads_transcript_assigned"], s["barcodes_unique_molecules"]]
    assert "missing" not in c["counts"]["reads_funnel"]
    src = b["modules"]["sc"].get("value_sources", {})
    assert src["pipeline_version"]["kind"] == "html_table" and src["workflow_parameters"]["kind"] == "html_table"
    assert "reads_input" not in src and "reads_valid_barcode" not in src           # these came from data files


def test_v2_molecules_counted_two_ways_agree(wf_v2):
    """The sum of the raw gene matrix is the number of distinct (barcode, gene, UMI) among gene-assigned reads."""
    fake, b = wf_v2
    assert b["modules"]["sc"]["params"]["cross_checks"]["molecules_from_read_table"] == summary(b)["barcodes_unique_molecules"]
    h = cell(b)["hists"]["saturation_reads_per_molecule"]
    want = np.array(fake.expected["molecule_hist"])
    assert sum(h["n"]) == want.sum() and sum(h["w"]) == int((want * np.arange(len(want))).sum())
    assert b["modules"]["sc"]["params"]["read_id_sketch"]["ids"] == fake.expected["sketch_ids"]


def test_v2_reads_per_unit(wf_v2):
    fake, b = wf_v2
    per = sorted(fake.expected["per_unit"].values(), reverse=True)
    knee = cell(b)["curves"]["barcodes_knee"]
    assert knee["y"]["all"][0] == per[0] and knee["x"][-1] == len(per) and knee["y_unit"] == "reads"
    h = cell(b)["hists"]["barcodes_reads_hist"]
    assert sum(h["n"]) + h["under"] + h["over"] == len(per)
    assert summary(b)["reads_chrm_fraction"] == pytest.approx(0.2)


def test_v3_numbers_are_exact_and_report_only_numbers_are_marked(wf_v3):
    fake, b = wf_v3
    _check_summary(fake, b)
    ana = b["identity"]["analysis"]
    assert (ana["assay"], ana["pipeline_version"], ana["unit"]) == ("visium_hd", "3.3.4", {"name": "bin", "size_um": 8.0})
    src = b["modules"]["sc"]["value_sources"]
    for key in ("reads_input", "reads_full_length", "reads_full_length_fraction", "reads_adapter_config", "saturation_reported"):
        assert src[key]["kind"] in ("html_table", "embedded_json"), key
    for key in ("reads_valid_barcode", "reads_gene_assigned", "reads_transcript_assigned"):
        assert key not in src, key                                # counted from the per-read table
    # a share of a report number and a counted number is marked too
    assert src["reads_gene_assigned_fraction"]["kind"] == "html_table"
    adapter = cell(b)["counts"]["reads_adapter_config"]
    assert adapter["n"][0] == round(0.987 * summary(b)["reads_input"]) and "rebuilt from the percentage" in adapter["sources"]["full_length"]
    assert b["modules"]["sc"]["params"]["sibling_hint"]["tail"] == "demo_analysis/spaceranger/outs"


def test_v3_units_are_8um_bins(wf_v3):
    fake, b = wf_v3
    e = fake.expected
    assert summary(b)["barcodes_detected"] == len(e["per_unit"])
    g = cell(b)["grids"]["spatial_map"]
    assert g["layout"] == {"kind": "square", "bin_um": 64.0} and set(g["channels"]) == {"umis", "genes", "reads"}
    umis = axes.channel_values(g["channels"]["umis"]).reshape(g["rows"], g["cols"])
    reads = axes.channel_values(g["channels"]["reads"]).reshape(g["rows"], g["cols"])
    assert umis.sum() == pytest.approx(sum(e["umis_by_unit"].values()), rel=0.03)
    assert reads.sum() == pytest.approx(sum(e["per_unit"].values()), rel=0.03)


def test_v3_links_from_the_read_tags(wf_v3):
    fake, b = wf_v3
    links = b["identity"]["analysis"]["links"]
    assert [(ln["flow_cell_id"], ln["protocol_run_id_prefix"], ln["acquisition_run_id_prefix"]) for ln in links] == \
        [("TST00001", "0a1b2c3d", "11112222"), ("TST00002", "4e5f6a7b", "33334444")]
    assert [ln["share"] for ln in links] == pytest.approx([0.6, 0.4], abs=0.02) and sum(ln["share"] for ln in links) == pytest.approx(1.0, abs=1e-3)
    assert all(ln["method"] == "read_tag" and ln["confidence"] == "high" and ln["resolved"] is False for ln in links)


def test_droplet_layout_with_its_own_adapter_statistics(wf_droplet):
    fake, b = wf_droplet
    _check_summary(fake, b)
    assert b["identity"]["analysis"]["assay"] == "singlecell" and b["identity"]["analysis"]["unit"] == {"name": "cell", "size_um": None}
    src = b["modules"]["sc"].get("value_sources", {})
    assert "reads_input" not in src and "reads_full_length" not in src      # config_stats.json was there: a data file
    assert summary(b)["reads_plus_fraction"] == pytest.approx(0.5, abs=0.01)
    per = np.sort(np.array(list(fake.expected["per_unit"].values())))[::-1]
    assert summary(b)["barcodes_cell_threshold"] == pytest.approx(np.percentile(per[:12], 95) / 20)
    assert "grids" not in cell(b)


def test_without_the_read_pass_counts_come_from_the_report_and_say_so(tmp_path):
    fake = fx.make_wf_v3(str(tmp_path), bam=False)
    b = collected(fake, read_pass=False)
    s, sec = summary(b), b["modules"]["sc"]
    assert s["reads_valid_barcode"] == fake.expected["summary"]["reads_valid_barcode"]
    assert sec["value_sources"]["reads_valid_barcode"]["kind"] == "html_table"
    assert sec["state"] == "partial" and any(n["code"] == "read_pass_skipped" for n in sec["notes"])
    assert "barcodes_reads_hist" not in cell(b)["hists"] and "saturation_reads_per_molecule" not in cell(b)["hists"]
    skipped = next(i for i in b["provenance"]["inputs"] if i["role"] == "wfsc.read_summary")
    assert skipped["state"] == "skipped"
    assert check_bundle(composed(b), where="fixture") == []


def test_without_the_report_the_numbers_it_alone_holds_are_not_available(tmp_path):
    fake = fx.make_wf_v3(str(tmp_path), bam=False)
    os.unlink(os.path.join(fake.path, fx.REPORT_NAME))
    b = collected(fx.FakeAnalysis(fake.files["sample"], "wf-single-cell", "wf_v3"))
    s = summary(b)
    assert "reads_input" not in s and "reads_full_length_fraction" not in s and b["identity"]["analysis"]["pipeline_version"] is None
    assert s["reads_valid_barcode"] == fake.expected["summary"]["reads_valid_barcode"]         # everything else is still reported
    assert cell(b)["counts"]["reads_funnel"]["missing"] == ["input", "full_length"]


def test_a_failed_workflow_task(tmp_path):
    fake = fx.make_wf_v2(str(tmp_path), failed_tasks=2)
    b = collected(fake)
    assert summary(b)["workflow_failed_tasks"] == 2 and summary(b)["workflow_tasks"] == 5
    rows = cell(b)["tables"]["workflow_execution"]["rows"]
    assert rows[0][:2] == ["assign_features", 2] and {r[0] for r in rows} == {"assign_features", "tag_bam", "call_paftools", "makeReport"}


def test_bundles_are_valid(wf_v2, wf_v3, wf_droplet):
    for _, b in (wf_v2, wf_v3, wf_droplet):
        assert check_bundle(composed(b), where="fixture") == []


def test_report_reader_ignores_scripts_and_reads_the_footer(tmp_path):
    fake = fx.make_wf_v3(str(tmp_path), bam=False)
    rep = wfsc.report(os.path.join(fake.path, fx.REPORT_NAME), "SAMPLE_S1_L001")
    assert rep["version"] == "3.3.4" and rep["adapter_pct"]["full_length"] == 98.7
    assert rep["saturation_reported"]["reads"] == fake.expected["summary"]["reads_valid_barcode"]
    assert wfsc.assay_of_kit("visium_hd:v1") == "visium_hd" and wfsc.assay_of_kit("visium:v1") == "visium" and wfsc.assay_of_kit("5prime:v2") == "singlecell"
