"""An analysis folder is recognised by the files in it, never by its name."""
import os
import shutil

from flowcellar.sc import ROLE_SPACERANGER, ROLE_WF, detect
from flowcellar.testing import sc as fx


def test_a_spaceranger_outs_folder_is_known_by_its_files_whatever_it_is_called(tmp_path):
    fake = fx.make_spaceranger_visium(str(tmp_path))
    c = detect.claims_at(fake.path)[0]
    assert (c.pipeline, c.layout, c.hd) == ("spaceranger", "spaceranger", False)
    assert c.variant == "spaceranger_out"                       # the folder that tells a re-run apart
    assert {"metrics_summary", "web_summary", "molecule_info", "raw_matrix", "tissue_positions", "pipestance_metrics", "cmdline"} <= set(c.files)
    renamed = tmp_path / "anything"
    shutil.copytree(fake.path, renamed)
    assert detect.claims_at(str(renamed))[0].layout == "spaceranger"
    assert detect.claims_at(os.path.dirname(fake.path))[0].path == fake.path      # the folder that holds outs/
    os.unlink(renamed / "web_summary.html")
    assert detect.claims_at(str(renamed)) == []


def test_visium_hd_is_told_by_its_files(tmp_path):
    c = detect.claims_at(fx.make_spaceranger_hd(str(tmp_path)).path)[0]
    assert c.hd and "feature_slice" in c.files and "filtered_matrix_008um" in c.files


def test_both_workflow_layouts(tmp_path):
    v2 = fx.make_wf_v2(str(tmp_path / "a"))
    v3 = fx.make_wf_v3(str(tmp_path / "b"))
    c2, c3 = detect.claims_at(v2.path), detect.claims_at(v3.path)
    assert [c.layout for c in c2] == ["wf_v2"] and [c.layout for c in c3] == ["wf_v3"]
    assert c2[0].alias == "demo_spatial_a" and not c2[0].hd and {"read_summary", "config_stats", "gene_matrix", "whitelist", "report", "trace"} <= set(c2[0].files)
    assert c3[0].alias == "SAMPLE_S1_L001" and c3[0].hd and {"read_summary", "params", "gene_matrix_8um", "gene_matrix_2um", "bam", "bam_index"} <= set(c3[0].files)
    assert "gene_matrix" not in c3[0].files                     # the processed (normalised) matrices are never taken
    assert detect.claims_at(v3.files["sample"])[0].path == v3.files["sample"]    # one sample folder is enough


def test_copies_and_scratch_are_skipped(tmp_path):
    v2 = fx.make_wf_v2(str(tmp_path))
    checkpoint = os.path.join(v2.path, ".ipynb_checkpoints", detect.REPORT_NAME)
    assert os.path.isfile(checkpoint)
    assert detect.claims_from_anchor(ROLE_WF, checkpoint) == []
    assert detect.skipped("/x/SPATIAL_RNA_COUNTER_CS/y/outs/metrics_summary.csv") and detect.skipped("/x/work/ab/outs")
    assert not detect.skipped("/x/analysis/spaceranger/outs/metrics_summary.csv")
    assert len(detect.claims_from_anchor(ROLE_WF, os.path.join(v2.path, detect.REPORT_NAME))) == 1


def test_a_scratch_name_above_the_scanned_folder_hides_nothing(tmp_path):
    # the data of a laboratory may lie under /work/... or in a folder somebody named "work": only the names BELOW the
    # scanned folder say that a path is a copy or scratch
    root = tmp_path / "work" / "store"
    v2 = fx.make_wf_v2(str(root / "proj" / "sample"))
    anchor = os.path.join(v2.path, detect.REPORT_NAME)
    assert detect.skipped(anchor)                                                   # the whole path holds "work" ...
    assert detect.claims_from_anchor(ROLE_WF, anchor) == []                          # ... which hides it when nothing says where the root is
    rel = os.path.relpath(anchor, root)
    assert len(detect.claims_from_anchor(ROLE_WF, anchor, rel)) == 1                 # below the scanned folder nothing is scratch
    inside = fx.make_wf_v2(str(root / "proj" / "work" / "ab"))
    anchor2 = os.path.join(inside.path, detect.REPORT_NAME)
    assert detect.claims_from_anchor(ROLE_WF, anchor2, os.path.relpath(anchor2, root)) == []     # a work folder below the root still is


def test_the_anchor_rules_name_files_not_folders_of_a_laboratory():
    rules = detect.input_rules()
    assert {r.role for r in rules} == {ROLE_SPACERANGER, ROLE_WF}
    assert sorted(r.glob for r in rules) == ["**/metrics_summary.csv", "**/wf-single-cell-report.html"], \
        "file names looked for in every folder: a rule that starts at a folder name would stop the walk at every 'outs'"


def test_a_folder_that_is_no_analysis(tmp_path):
    (tmp_path / "outs").mkdir()
    (tmp_path / "outs" / "metrics_summary.csv").write_text("a,b\n1,2\n")
    assert detect.claims_at(str(tmp_path)) == [] and detect.claims_at(str(tmp_path / "nothing")) == []
