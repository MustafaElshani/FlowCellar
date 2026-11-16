"""``flowcellar analysis``, and the promise that the core never loads the add-on's packages."""
from __future__ import annotations

import gzip
import json
import os
import subprocess
import sys

import pytest
from click.testing import CliRunner

from flowcellar.cli import cli, main
from flowcellar.library.maintain import validate_paths
from flowcellar.testing import sc as fx
from flowcellar.util.errors import FlowCellarError


@pytest.fixture(scope="module")
def viewer():
    from flowcellar.report import api as report_api
    try:
        return report_api.viewer_manifest()
    except FlowCellarError as e:
        pytest.skip(f"the viewer is not built: {e.message}")


def _bundle(path: str) -> dict:
    with gzip.open(path) as fh:
        return json.load(fh)


def test_one_folder_gives_a_bundle_and_a_one_file_report(tmp_path, viewer):
    fake = fx.make_spaceranger_hd(str(tmp_path / "data"))
    out = tmp_path / "qc"
    r = CliRunner().invoke(cli, ["analysis", fake.path, "-o", str(out)])
    assert r.exit_code == 0, r.output
    names = sorted(os.listdir(out))
    assert len(names) == 2 and names[0].endswith(".fcb.json.gz") and names[1].endswith(".flowcellar.html")
    aid = names[0][: -len(".fcb.json.gz")]
    b = _bundle(str(out / names[0]))
    assert b["bundle_id"] == aid and b["kind"] == "analysis" and b["status"]["verdict"] == "warn"
    assert [x for x in validate_paths([str(out / names[0])]) if not x["ok"]] == []
    html = (out / names[1]).read_text("utf-8")
    assert f'id="fcb:{aid}"' in html and 'id="fci"' in html and '"analysis":"' + aid + '"' in html
    assert "UMI Registration Warning" in html                     # the text summary, readable without scripts
    import html as html_lib
    head = html_lib.unescape(html.split("<title>")[0])
    assert "connect-src 'none'" in head and "default-src 'none'" in head       # the page asks the network for nothing
    assert "demo_hd_day2" in html and fake.path not in html.split('id="fc-boot"')[0]


def test_both_pipelines_of_one_sample_go_into_one_report(tmp_path, viewer):
    wf, sr = fx.make_sibling_pair(str(tmp_path / "data"))
    out = tmp_path / "qc"
    r = CliRunner().invoke(cli, ["analysis", wf.path, "-o", str(out)])
    assert r.exit_code == 0, r.output
    bundles = sorted(n for n in os.listdir(out) if n.endswith(".fcb.json.gz"))
    reports = [n for n in os.listdir(out) if n.endswith(".flowcellar.html")]
    assert len(bundles) == 2 and len(reports) == 1 and reports[0].startswith("an_wfsinglecell_")
    html = (out / reports[0]).read_text("utf-8")
    for n in bundles:
        assert f'id="fcb:{n[: -len(".fcb.json.gz")]}"' in html
    w = _bundle(str(out / [n for n in bundles if "wfsinglecell" in n][0]))
    s = _bundle(str(out / [n for n in bundles if "spaceranger" in n][0]))
    assert w["identity"]["analysis"]["sibling_analyses"] == [s["bundle_id"]] and s["identity"]["analysis"]["sibling_analyses"] == [w["bundle_id"]]
    assert w["identity"]["analysis"]["group_key"] == s["identity"]["analysis"]["group_key"] == "H1-TESTSLD-A1"
    assert [ln["method"] for ln in s["identity"]["analysis"]["links"]] == ["sibling_analysis", "sibling_analysis"]
    assert any(r["id"] == "sp_wf_sr_umi_ratio" for r in w["status"]["rules"])
    # --no-sibling: the one folder that was named, and nothing else
    out2 = tmp_path / "qc2"
    assert CliRunner().invoke(cli, ["analysis", wf.path, "-o", str(out2), "--no-sibling", "--no-report"]).exit_code == 0
    assert len(os.listdir(out2)) == 1


def test_stated_runs_labels_and_gates(tmp_path, viewer):
    fake = fx.make_wf_v2(str(tmp_path / "data"))
    out = tmp_path / "qc"
    r = CliRunner().invoke(cli, ["analysis", fake.path, "-o", str(out), "--run", "TST00005@aabbccdd", "--label", "my sample", "--no-read-pass",
                                 "--prefix", "one"])
    assert r.exit_code == 0, r.output
    b = _bundle(str(out / "one.fcb.json.gz"))
    assert b["identity"]["analysis"]["label"] == "my sample"
    assert [(ln["flow_cell_id"], ln["protocol_run_id_prefix"], ln["method"]) for ln in b["identity"]["analysis"]["links"]] == [("TST00005", "aabbccdd", "manual")]
    assert os.path.isfile(out / "one.flowcellar.html")
    sr = fx.make_spaceranger_visium(str(tmp_path / "sr"), antisense=0.2)
    with pytest.raises(SystemExit) as e:
        main(["analysis", sr.path, "-o", str(tmp_path / "g"), "--fail-on", "warn"])
    assert e.value.code == 3


@pytest.mark.parametrize("argv, code", [
    (["analysis", "{tmp}/nothing", "-o", "{tmp}/o"], 4),            # no such folder
    (["analysis", "{tmp}", "-o", "{tmp}/o"], 4),                    # not an analysis folder
    (["analysis", "{tmp}"], 2),                                     # nowhere to write
    (["analysis", "{tmp}", "-o", "{tmp}/o", "--run", "what"], 2),
])
def test_exit_codes(tmp_path, capsys, argv, code):
    with pytest.raises(SystemExit) as e:
        main([a.format(tmp=tmp_path) for a in argv])
    assert e.value.code == code
    assert "flowcellar" in capsys.readouterr().err or code == 2


def test_adding_one_analysis_to_a_library(tmp_path, viewer, fake_run):
    lib_dir = str(tmp_path / "library")
    assert CliRunner().invoke(cli, ["run", fake_run.run_dir, "-L", lib_dir, "--quick"]).exit_code == 0      # a library with one run
    wf, sr = fx.make_sibling_pair(str(tmp_path / "data"))
    r = CliRunner().invoke(cli, ["analysis", sr.path, "-L", lib_dir, "--run", fake_run.flow_cell_id])
    assert r.exit_code == 0, r.output
    with gzip.open(os.path.join(lib_dir, "index", "index.json.gz")) as fh:
        idx = json.load(fh)
    assert idx["runs"]["n"] == 1 and len(idx["analyses"]) == 2     # the workflow analysis of the same reads came with it
    row = next(a for a in idx["analyses"] if a["pipeline"] == "spaceranger")
    assert [(ln["run"], ln["method"], ln["resolved"]) for ln in row["links"]] == [(fake_run.run_key, "manual", True)]
    other = next(a for a in idx["analyses"] if a["pipeline"] == "wf-single-cell")
    assert [ln["resolved"] for ln in other["links"]] == [False, False] and other["siblings"] == [row["id"]]     # its reads name other runs
    assert all(x["ok"] for x in validate_paths([lib_dir]))
    assert os.path.isfile(os.path.join(lib_dir, "hub", "data", row["id"] + ".js"))


CORE_IMPORTS = ("import flowcellar, flowcellar.cli, flowcellar.engine.api, flowcellar.library.api, flowcellar.readdata.api, "
                "flowcellar.report.api, flowcellar.library.scan, flowcellar.library.index, flowcellar.report.hub, flowcellar.serve.server\n")
CHECK = ("import sys\n"
         "heavy = [m for m in sys.modules if m == 'h5py' or m.startswith('h5py.')]\n"
         "assert not heavy, heavy\n")


def _python(code: str) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)


def test_importing_the_core_loads_nothing_of_the_add_on():
    done = _python(CORE_IMPORTS + CHECK + "assert not [m for m in sys.modules if m.startswith('flowcellar.sc')], 'the add-on was imported'\n")
    assert done.returncode == 0, done.stderr[-2000:]


def test_a_run_report_loads_nothing_of_the_add_on(tmp_path, fake_run, viewer):
    code = ("from flowcellar.cli import main\n"
            "try:\n"
            f"    main(['-q', 'run', {fake_run.run_dir!r}, '-o', {str(tmp_path / 'qc')!r}])\n"
            "except SystemExit as e:\n"
            "    assert e.code == 0, e.code\n"
            + CHECK + "assert not [m for m in sys.modules if m.startswith('flowcellar.sc')], 'the add-on was imported'\n")
    done = _python(code)
    assert done.returncode == 0, done.stderr[-2000:]
    assert any(n.endswith(".flowcellar.html") for n in os.listdir(tmp_path / "qc"))


def test_the_add_on_itself_loads_h5py_only_to_open_a_file():
    done = _python("import flowcellar.sc, flowcellar.sc.detect, flowcellar.sc.collect, flowcellar.sc.collect_sr, flowcellar.sc.collect_wf, "
                   "flowcellar.sc.store, flowcellar.sc.judge, flowcellar.sc.link, flowcellar.sc.report, flowcellar.sc.scan, flowcellar.cli.analysis\n"
                   + CHECK + "assert 'pysam' not in sys.modules\n")
    assert done.returncode == 0, done.stderr[-2000:]


def test_the_extra_is_declared():
    import tomllib
    from pathlib import Path
    meta = tomllib.loads((Path(__file__).resolve().parents[2] / "pyproject.toml").read_text("utf-8"))
    extras = meta["project"]["optional-dependencies"]
    assert any(d.startswith("h5py") for d in extras["sc"])
    assert not any("h5py" in d for d in meta["project"]["dependencies"])
