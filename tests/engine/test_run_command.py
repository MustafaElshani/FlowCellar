"""flowcellar run: what it writes, what it prints, its exit codes."""
import json
import os
import re

import pytest

from flowcellar.cli import main
from flowcellar.model import bundle as B
from flowcellar.model import schema
from flowcellar.testing import minknow_variants as V
from flowcellar.testing.fixtures import RunSpec, make_run


def _viewer_built() -> bool:
    try:
        from flowcellar.report import api as report
        report.viewer_manifest()
        return True
    except Exception:  # noqa: BLE001
        return False


def run_cli(args, capsys):
    """(exit code, standard error text) of ``flowcellar <args>``."""
    with pytest.raises(SystemExit) as e:
        main([str(a) for a in args])
    cap = capsys.readouterr()
    return int(e.value.code or 0), cap.err


@pytest.fixture(scope="module")
def run_folder(tmp_path_factory):
    return make_run(tmp_path_factory.mktemp("cli_run"), RunSpec(flow_cell_id="TST00090", seed=90, reads=2500))


def test_run_writes_a_bundle_and_prints_its_lines(run_folder, tmp_path, capsys):
    out = tmp_path / "qc"
    code, err = run_cli(["run", run_folder.run_dir, "-o", out, "--no-report"], capsys)
    assert code == 0, err
    lines = err.splitlines()
    assert lines[0] == f"{run_folder.run_key}  MinION Mk1D · FLO-MIN114 · SQK-PCB114-24 · cDNA · 8.0 h"
    health = [ln for ln in lines if ln.startswith("  Run health")]
    assert re.search(r"so far · 10 checks done, 5 wait for the read count \(\d+ (kB|MB) to read\)$", health[0])
    assert re.search(r"Fail · 15 checks · reads at 6 h: [\d,]+", health[1])
    assert any(ln.startswith("  reading        sequencing_summary_") for ln in lines)
    assert lines[-1] == "  Read data      Not available yet"
    path = out / f"{run_folder.run_key}.fcb.json.gz"
    b = B.read_bundle(path)
    assert schema.schema_errors(b) == [] and b["tiers"]["run_health"]["detail"] == "full"
    assert b["status"]["by_tier"]["run_health"] == "fail" and b["status"]["judged"]["library"] is False
    assert b["relations"]["flow_cell_use"]["source"] == "unknown" and "user" in b
    assert (out / ".flowcellar-partials").is_dir() and not list(out.glob("*.html"))
    assert B.content_sha256(b) == b["content_sha256"]


def test_quick_stops_after_the_small_files(run_folder, tmp_path, capsys):
    code, err = run_cli(["run", run_folder.run_dir, "-o", tmp_path, "--quick", "--no-report", "--prefix", "mine"], capsys)
    assert code == 0
    assert "(--quick: the reads are not counted)" in err and "reading" not in err
    b = B.read_bundle(tmp_path / "mine.fcb.json.gz")
    assert b["tiers"]["run_health"]["detail"] == "quick" and b["status"]["provisional"] is True
    assert not (tmp_path / ".flowcellar-partials").exists()


@pytest.mark.skipif(not _viewer_built(), reason="the viewer is not built")
def test_run_writes_the_one_file_report(run_folder, tmp_path, capsys):
    code, err = run_cli(["run", run_folder.run_dir, "-o", tmp_path], capsys)
    assert code == 0, err
    html = (tmp_path / f"{run_folder.run_key}.flowcellar.html").read_text("utf-8")
    assert html.lstrip().lower().startswith("<!doctype html") and run_folder.run_key in html
    assert "  quick report   " in err and "  full report    " in err


def test_exit_codes(run_folder, tmp_path, capsys):
    code, err = run_cli(["run", run_folder.run_dir], capsys)
    assert code == 2 and "nowhere to write" in err
    code, err = run_cli(["run", tmp_path, "-o", tmp_path / "x"], capsys)
    assert code == 4 and "not a MinKNOW run folder" in err
    look = V.make_lookalike(tmp_path / "look")
    code, err = run_cli(["run", look, "-o", tmp_path / "x"], capsys)
    assert code == 4 and "basecaller output" in err
    code, err = run_cli(["run", run_folder.run_dir, "-o", tmp_path / "g", "--quick", "--no-report", "--fail-on", "warn"], capsys)
    assert code == 3 and "--fail-on warn" in err
    assert (tmp_path / "g" / f"{run_folder.run_key}.fcb.json.gz").exists()          # the files are written before the gate trips
    code, err = run_cli(["run", run_folder.run_dir, "-o", tmp_path / "h", "--reused", "--flow-cell-use", "3"], capsys)
    assert code == 2
    code, err = run_cli(["run", run_folder.run_dir, "-o", tmp_path / "h", "--input", "nonsense"], capsys)
    assert code == 2 and "ROLE=PATH" in err
    code, err = run_cli(["run", run_folder.run_dir, "-o", tmp_path / "h", "--prefix", "../up"], capsys)
    assert code == 2


def test_a_gate_does_not_trip_below_its_level(tmp_path, capsys):
    f = V.make_short(tmp_path / "run", minutes=5, flow_cell_id="TST00091", seed=91)
    code, _ = run_cli(["run", f.run_dir, "-o", tmp_path / "a", "--no-report", "--fail-on", "warn"], capsys)
    assert code == 0                                                 # incomplete is below a warning
    code, _ = run_cli(["run", f.run_dir, "-o", tmp_path / "b", "--no-report", "--fail-on", "incomplete"], capsys)
    assert code == 3


def test_deterministic_gives_identical_bytes(run_folder, tmp_path, capsys):
    for name in ("a", "b"):
        code, _ = run_cli(["run", run_folder.run_dir, "-o", tmp_path / name, "--no-report", "--deterministic"], capsys)
        assert code == 0
    name = f"{run_folder.run_key}.fcb.json.gz"
    assert (tmp_path / "a" / name).read_bytes() == (tmp_path / "b" / name).read_bytes()
    b = B.read_bundle(tmp_path / "a" / name)
    assert "timings" not in b["provenance"] and b["status"]["evaluated"] == b["provenance"]["created"]


def test_what_is_stated_on_the_command_line_reaches_the_verdict(run_folder, tmp_path, capsys):
    code, err = run_cli(["run", run_folder.run_dir, "-o", tmp_path, "--quick", "--no-report", "--reused", "--assay", "gdna",
                         "--baseline-group", "pilot"], capsys)
    assert code == 0, err
    b = B.read_bundle(tmp_path / f"{run_folder.run_key}.fcb.json.gz")
    assert b["relations"]["flow_cell_use"]["index"] == 2 and b["relations"]["flow_cell_use"]["source"] == "declared"
    assert b["status"]["modifiers"][0]["code"] == "reused_flow_cell"
    assert b["identity"]["assay"] == {**b["identity"]["assay"], "type": "gdna", "group": "gdna", "source": "user", "confidence": "certain"}
    assert "genomic DNA" in err.splitlines()[0]


def test_raw_only_and_running_folders(tmp_path, capsys):
    raw = V.make_raw_only(tmp_path / "raw", flow_cell_id="TST00092", seed=92)
    code, err = run_cli(["run", raw.run_dir, "-o", tmp_path / "o1", "--no-report"], capsys)
    assert code == 0 and err.splitlines()[0] == f"{raw.run_key}  raw data only"
    assert "Incomplete" in err and "Only raw data was found" in err
    b = B.read_bundle(tmp_path / "o1" / f"{raw.run_key}.fcb.json.gz")
    assert b["status"]["cause"] == "raw_only" and schema.schema_errors(b) == []
    live = V.make_running(tmp_path / "live", flow_cell_id="TST00093", seed=93)
    code, err = run_cli(["run", live.run_dir, "-o", tmp_path / "o2", "--no-report"], capsys)
    assert code == 0 and "Running" in err
    assert B.read_bundle(tmp_path / "o2" / f"{live.run_key}.fcb.json.gz")["status"]["verdict"] == "running"


def test_json_log_lines(run_folder, tmp_path, capsys):
    code, err = run_cli(["--log-format", "json", "run", run_folder.run_dir, "-o", tmp_path, "--quick", "--no-report"], capsys)
    assert code == 0
    events = [json.loads(ln) for ln in err.splitlines()]
    assert [e["event"] for e in events] == ["run.start", "run.health", "run.bundle", "run.data"]
    assert events[0]["run_key"] == run_folder.run_key and os.path.isabs(events[2]["path"])


def test_a_damaged_summary_ends_with_exit_code_5_and_a_bundle(tmp_path, capsys):
    f = make_run(tmp_path / "run", RunSpec(flow_cell_id="TST00094", seed=94))
    with open(f.files["minknow.sequencing_summary"], "w") as fh:
        fh.write("not\ta\tsummary\n1\t2\t3\n")
    code, err = run_cli(["run", f.run_dir, "-o", tmp_path / "o", "--no-report"], capsys)
    assert code == 5 and "not read" in err
    b = B.read_bundle(tmp_path / "o" / f"{f.run_key}.fcb.json.gz")
    assert b["modules"]["dyn"]["state"] == "error" and schema.schema_errors(b) == []


def _pysam_installed() -> bool:
    try:
        import pysam  # noqa: F401
        return True
    except ImportError:
        return False


@pytest.mark.bam
@pytest.mark.skipif(not _pysam_installed(), reason="needs pysam (the 'bam' extra)")
def test_run_with_read_data_goes_through_the_read_data_contract(tmp_path, capsys):
    from flowcellar.testing import readdata as R
    f = make_run(tmp_path / "runs", RunSpec(flow_cell_id="TST00095", seed=95, reads=2600))
    proj = R.make_project(str(tmp_path / "project"), R.run_ref(f), reads_per_sample=200, unclassified=40)
    code, err = run_cli(["run", f.run_dir, proj.project_dir, "--sample-sheet", proj.sample_sheet, "-o", tmp_path / "out", "--no-report"], capsys)
    assert code == 0, err
    b = B.read_bundle(tmp_path / "out" / f"{f.run_key}.fcb.json.gz")
    assert schema.schema_errors(b) == []
    assert list(b["modules"])[:5] == ["run", "fc", "yield", "dyn", "len.est"] and {"len", "qual", "funnel"} <= set(b["modules"])
    assert b["tiers"]["read_level"]["state"] in ("partial", "complete") and len(b["samples"]) >= 3
    assert b["status"]["by_tier"]["read_level"] != "not_available"
    assert any(i["id"].startswith("rd") for i in b["provenance"]["inputs"]) and len(b["provenance"]["partials"]) > 2
    assert err.splitlines()[-1].startswith("  Read data      ") and "Not available yet" not in err.splitlines()[-1]
    assert b["tiers"]["read_level"]["reads_expected"] == f.expected["reads"]
