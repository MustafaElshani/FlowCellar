"""Space Ranger output: every number equals what the fixture wrote, and what the pipeline said is kept word for word."""
import os

import numpy as np
import pytest

from flowcellar.library.maintain import check_bundle
from flowcellar.sc import SECTION_BUDGET_BYTES, axes, maps
from flowcellar.sc import spaceranger as sr
from flowcellar.testing import sc as fx

from .conftest import cell, collected, composed, summary


def _same(got, want):
    if isinstance(want, float):
        return got == pytest.approx(want, rel=1e-12, abs=1e-12)
    return got == want


def test_visium_numbers_are_exact(sr_visium):
    fake, b = sr_visium
    s = summary(b)
    for key, want in fake.expected["summary"].items():
        if want is None or key.startswith("reads_bam"):
            continue
        assert _same(s.get(key), want), key
    ana = b["identity"]["analysis"]
    assert (ana["label"], ana["assay"], ana["pipeline_version"], ana["reference"]) == ("demo_spatial_a", "visium", "3.0.1", "demo-2024-A")
    assert ana["slide_serial_capture_area"] == ana["group_key"] == "V00T00-001-A1"
    assert ana["unit"] == {"name": "spot", "size_um": 55.0}
    assert b["bundle_id"].startswith("an_spaceranger_") and b["kind"] == "analysis"
    sc = cell(b)["scalars"]
    assert sc["reads_bam_mapped"] == 400 and sc["reads_bam_unmapped"] == 50 and sc["spatial_spots"] == 80


def test_visium_what_flowcellar_computes_agrees_with_the_pipeline(sr_visium):
    """The same numbers two ways: read from metrics_summary.csv, and recomputed from the matrix, the molecule table
    and the spot positions."""
    fake, b = sr_visium
    s, x = summary(b), b["modules"]["sc"]["params"]["cross_checks"]
    assert x["units_under_tissue"] == s["spatial_units_under_tissue"]
    assert x["umis_per_unit_median"] == s["spatial_umis_per_unit"] and x["genes_per_unit_median"] == s["spatial_genes_per_unit"]
    assert x["genes_detected_under_tissue"] == s["spatial_genes_detected"]
    assert x["fraction_reads_under_tissue"] == pytest.approx(s["spatial_fraction_reads_under_tissue"], abs=1e-15)
    assert x["saturation_from_molecules"] == pytest.approx(s["sr_saturation"], abs=1e-15)


def test_visium_histograms_and_curves(sr_visium):
    fake, b = sr_visium
    c, e = cell(b), fake.expected
    h = c["hists"]["saturation_reads_per_molecule"]
    want = np.array(e["count_hist"])
    assert sum(h["n"]) == want.sum() and sum(h["w"]) == int((want * np.arange(len(want))).sum())
    under = [u for u, t in zip(e["umis_by_spot"], e["under"]) if t]
    hu = c["hists"]["barcodes_umis_hist"]
    assert sum(hu["n"]) + hu["under"] + hu["over"] == len(under) and hu["max"] == max(under)
    knee = c["curves"]["barcodes_knee"]
    assert knee["y"]["under_tissue"][0] == max(under) and knee["x"][0] == 1 and knee["x"][-1] == len([u for u in under if u > 0])
    off = [u for u, t in zip(e["umis_by_spot"], e["under"]) if not t and u > 0]
    assert c["curves"]["barcodes_knee_off_tissue"]["x"][-1] == len(off)
    assert c["counts"]["sr_read_funnel"]["n"] == [e["summary"]["sr_reads"], e["summary"]["sr_reads_feature"], e["summary"]["sr_reads_usable"]]


def test_visium_map_holds_every_spot_at_its_place(sr_visium):
    fake, b = sr_visium
    g = cell(b)["grids"]["spatial_map"]
    e = fake.expected
    assert (g["rows"], g["cols"], g["layout"]) == (8, 10, {"kind": "hex", "bin_um": 100.0})
    tissue = axes.channel_values(g["channels"]["tissue"]).reshape(8, 10)
    umis = axes.channel_values(g["channels"]["umis"]).reshape(8, 10)
    genes = axes.channel_values(g["channels"]["genes"]).reshape(8, 10)
    for i, (r, ac) in enumerate(e["spots_rc"]):
        assert tissue[r, ac // 2] == e["under"][i]
        assert umis[r, ac // 2] == pytest.approx(e["umis_by_spot"][i], rel=0.03, abs=0.51)      # a picture: one byte per spot
        assert genes[r, ac // 2] == pytest.approx(e["genes_by_spot"][i], rel=0.03, abs=0.51)
    assert set(g["channels"]) == {"tissue", "umis", "genes", "mito", "saturation"}
    assert np.isnan(maps.sqrt_values(g["channels"]["mito"])).sum() == 2                         # the two spots without a molecule


def test_alarms_are_kept_word_for_word_and_marked_as_from_the_report(sr_visium):
    fake, b = sr_visium
    sec = b["modules"]["sc"]
    alarms = [n for n in sec["notes"] if n["code"] == "pipeline_alarm"]
    assert [a["title"] for a in alarms] == ["High Fraction of Reads Mapped Antisense to Genes"]
    assert alarms[0]["message"] == "Ideal < 10%. This can indicate use of an unsupported chemistry type."
    assert alarms[0]["value"] == "20.0%" and alarms[0]["source"] == "spaceranger" and alarms[0]["level"] == "warning"
    src = sec["value_sources"]
    assert src["spatial_alarms_warn"]["kind"] == "embedded_json" and src["spatial_sample"]["kind"] == "embedded_json"
    assert "sr_corrected_bc_frac" not in src                    # from the pipestance's own JSON: a data file
    web = next(i for i in b["provenance"]["inputs"] if i["role"] == "spaceranger.web_summary")
    assert src["spatial_alarms_warn"]["input"] == web["id"]


def test_a_stale_pipestance_is_not_used(tmp_path):
    fake = fx.make_spaceranger_visium(str(tmp_path), stale_pipestance=True)
    b = collected(fake)
    assert "sr_corrected_bc_frac" not in summary(b)
    stale = next(i for i in b["provenance"]["inputs"] if i["role"] == "spaceranger.pipestance_metrics")
    assert stale["state"] == "skipped" and "differs" in stale["notes"][0]["message"]
    assert summary(b)["sr_tso_frac"] == 0.0003 and b["modules"]["sc"]["value_sources"]["sr_tso_frac"]["kind"] == "embedded_json"


def test_hd_numbers_are_exact(sr_hd):
    fake, b = sr_hd
    s = summary(b)
    for key, want in fake.expected["summary"].items():
        assert _same(s.get(key), want), key
    for key, want in fake.expected["scalars"].items():
        assert cell(b)["scalars"][key] == want, key
    ana = b["identity"]["analysis"]
    assert (ana["assay"], ana["pipeline_version"], ana["slide_serial_capture_area"]) == ("visium_hd", "4.1.0", "H1-TESTSLD-A1")
    assert ana["unit"] == {"name": "bin", "size_um": 8.0}
    x = b["modules"]["sc"]["params"]["cross_checks"]
    assert x["squares_under_tissue"] == s["spatial_squares_under_tissue_2um"] and x["units_under_tissue"] == s["spatial_units_under_tissue"]
    assert x["umis_per_unit_mean"] == pytest.approx(s["spatial_umis_per_unit"], rel=1e-12)
    assert x["genes_per_unit_mean"] == pytest.approx(s["spatial_genes_per_unit"], rel=1e-12)
    assert x["saturation_from_molecules"] == pytest.approx(s["sr_saturation"], abs=1e-15)


def test_hd_tables(sr_hd):
    fake, b = sr_hd
    t = cell(b)["tables"]
    levels = t["spatial_bin_levels"]
    assert [c["id"] for c in levels["columns"]][:2] == ["bin_um", "units_under_tissue"] and [r[0] for r in levels["rows"]] == [2, 8, 16]
    assert levels["rows"][1][1] == fake.expected["summary"]["spatial_units_under_tissue"]
    seg = dict(t["spatial_segmentation"]["rows"])
    assert seg["cells"] == 412 and seg["umis_in_cells"] == 0.6429 and seg["median_cell_area_um2"] == 164.0      # two spellings of the micro sign
    reg = dict(t["spatial_registration"]["rows"])
    assert reg["alignment"] == "manual alignment file" and reg["is_gaussian_fit_outside_range"] == "true" and reg["gaussian_fit_x_offset"] == "14.96"
    assert dict(t["workflow_parameters"]["rows"])["transcriptome"] == "refs/refdata-demo-2024-A"       # a path is cut to two components


def test_hd_map_is_the_2um_grid_summed_into_bins(sr_hd):
    fake, b = sr_hd
    e = fake.expected
    g = cell(b)["grids"]["spatial_map"]
    assert (g["rows"], g["cols"], g["layout"]) == (4, 4, {"kind": "square", "bin_um": 64.0})        # 100 squares of 2 um, bins of 32 squares
    assert set(g["channels"]) == {"umis", "saturation", "tissue", "mito", "genes", "image"}
    f = 32
    want_u = np.add.reduceat(np.add.reduceat(e["umis"], np.arange(0, 100, f), axis=0), np.arange(0, 100, f), axis=1)
    got_u = axes.channel_values(g["channels"]["umis"]).reshape(4, 4)
    assert got_u == pytest.approx(want_u, rel=0.03, abs=0.51)
    cells = np.outer([32, 32, 32, 4], [32, 32, 32, 4])
    want_t = np.add.reduceat(np.add.reduceat(e["tissue"].astype(float), np.arange(0, 100, f), axis=0), np.arange(0, 100, f), axis=1) / cells
    assert axes.channel_values(g["channels"]["tissue"]).reshape(4, 4) == pytest.approx(want_t, abs=0.003)
    want_r = np.add.reduceat(np.add.reduceat(e["reads"], np.arange(0, 100, f), axis=0), np.arange(0, 100, f), axis=1)
    assert axes.channel_values(g["channels"]["saturation"]).reshape(4, 4) == pytest.approx(1 - want_u / want_r, abs=0.003)
    want_g = np.add.reduceat(np.add.reduceat(e["grey"].astype(float), np.arange(0, 100, f), axis=0), np.arange(0, 100, f), axis=1) / cells
    assert axes.channel_values(g["channels"]["image"]).reshape(4, 4) == pytest.approx(want_g, abs=0.51)
    # the exact totals are numbers of their own, never read off the picture
    assert cell(b)["scalars"]["total_umis"] == int(e["umis"].sum())


def test_hd_map_options(tmp_path):
    fake = fx.make_spaceranger_hd(str(tmp_path))
    b32 = collected(fake, map_um=32, tissue_image=False)
    g = cell(b32)["grids"]["spatial_map"]
    assert (g["rows"], g["layout"]["bin_um"]) == (7, 32.0) and "image" not in g["channels"]
    assert b32["modules"]["sc"]["params"]["map"]["tissue_image"] is False


def test_bundles_are_valid_and_inside_the_size_budget(sr_visium, sr_hd):
    for _, b in (sr_visium, sr_hd):
        assert check_bundle(composed(b), where="fixture") == []
        assert b["modules"]["sc"]["params"]["section_bytes"] < SECTION_BUDGET_BYTES
        assert set(b["axes"]) == {"log_counts_8_per_decade", "fraction_50", "reads_per_molecule"}


def test_without_h5py_the_summary_files_are_still_read(tmp_path, monkeypatch):
    fake = fx.make_spaceranger_visium(str(tmp_path))
    monkeypatch.setattr(sr, "h5py_or_none", lambda: None)
    b = collected(fake)
    s = summary(b)
    assert s["spatial_units_under_tissue"] == fake.expected["summary"]["spatial_units_under_tissue"] and s["sr_antisense"] == 0.04
    sec = b["modules"]["sc"]
    assert sec["state"] == "partial" and sec["partial_reason"] == "inputs_missing"
    note = next(n for n in sec["notes"] if n["code"] == "h5py_missing")
    assert 'pip install "flowcellar[sc]"' in note["message"]
    assert "grids" not in cell(b) and "barcodes_umis_median" not in s
    assert check_bundle(composed(b), where="fixture") == []


def test_an_unknown_report_layout_costs_the_alarms_only(tmp_path):
    fake = fx.make_spaceranger_visium(str(tmp_path))
    with open(os.path.join(fake.path, "web_summary.html"), "w") as fh:
        fh.write("<html><body>a page of another kind</body></html>")
    b = collected(fake)
    assert summary(b)["spatial_alarms_known"] is False and summary(b)["sr_reads"] == fake.expected["summary"]["sr_reads"]
    assert any(n["code"] == "report_unreadable" for n in b["modules"]["sc"]["notes"])


def test_copies_of_one_output_are_one_analysis_and_a_rerun_is_another(tmp_path):
    import shutil
    a = fx.make_spaceranger_visium(str(tmp_path / "a"))
    shutil.copytree(os.path.dirname(a.path), tmp_path / "copy" / "renamed_out")
    again = fx.make_spaceranger_visium(str(tmp_path / "b"), seed=9)
    ids = [collected(x)["bundle_id"] for x in (a, fx.FakeAnalysis(str(tmp_path / "copy" / "renamed_out" / "outs"), "spaceranger", "spaceranger"), again)]
    assert ids[0] == ids[1] != ids[2]


def test_bam_header_and_index_readers(sr_visium):
    fake, _ = sr_visium
    m = sr.bam_mapping(os.path.join(fake.path, "possorted_genome_bam.bam"))
    assert (m["mapped"], m["unmapped"], m["chrm_mapped"]) == (400, 50, 100)
    assert sr.slide_from_header(m["header_text"]) == "V00T00-001-A1"
