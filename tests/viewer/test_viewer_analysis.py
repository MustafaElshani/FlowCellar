"""The screens of the single-cell and spatial add-on, in a real browser (viewer/verify/analysis.mjs).

A library is made by the engine itself from the synthetic spatial tree of ``flowcellar.testing.sc`` (one run, a
Visium sample and a Visium HD sample, each analysed by both pipelines; invented names only) and its hub is opened in
headless Firefox at 1,366 and 2,560 px in both themes: the Analyses tab, the page of an analysis with its spatial
map (pointer, keyboard, the switch, the table), the panel on the run page, the block on Home and a one-file report::

    FLOWCELLAR_BROWSER=1 .venv/bin/python -m pytest tests/viewer/test_viewer_analysis.py -q

The first tests need no browser: they check what the hub and the report hand to the viewer.
"""
from __future__ import annotations

import base64
import gzip
import json
import os
import re
import shutil
import subprocess
from pathlib import Path

import pytest

pytest.importorskip("h5py", reason="the single-cell and spatial add-on reads HDF5 with h5py (the 'sc' extra)")
pytest.importorskip("pysam", reason="the fixtures write BAM files with pysam (the 'bam' extra)")

from flowcellar.library import scan as scan_mod  # noqa: E402
from flowcellar.sc import report as sc_report  # noqa: E402
from flowcellar.sc.store import AnalysisStore  # noqa: E402
from flowcellar.testing import sc as fx  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
pytestmark = [pytest.mark.viewer]


@pytest.fixture(scope="module")
def spatial_library(tmp_path_factory, viewer_built):
    base = tmp_path_factory.mktemp("viewer_spatial")
    tree = fx.make_spatial_tree(str(base / "data"))
    lib_dir = base / "library"
    res = scan_mod.scan(str(lib_dir), scan_mod.ScanOptions(roots=(tree.roots["seq"], tree.roots["store"]), wait=False))
    assert res.exit_code == 0
    return {"dir": lib_dir, "tree": tree}


def hub_index(lib_dir: Path) -> dict:
    text = (lib_dir / "hub" / "index.js").read_text("utf-8")
    m = re.fullmatch(r'__FC_INDEX__\("([A-Za-z0-9+/=]*)"\);\n', text)
    assert m, text[:80]
    return json.loads(gzip.decompress(base64.b64decode(m.group(1))))


def test_the_hub_of_a_library_with_analyses_carries_them(spatial_library):
    lib_dir = spatial_library["dir"]
    ix = hub_index(lib_dir)
    rows = ix["analyses"]
    assert len(rows) == 4 and ix["runs"]["n"] == 1
    for r in rows:
        # the table of the Analyses tab is drawn from these alone
        assert {"id", "label", "pipeline", "assay", "verdict", "links", "siblings", "folder", "reasons", "reasons_n"} <= set(r)
        assert len(r["reasons"]) <= 2 and all(len(x) <= 200 for x in r["reasons"])
        side = lib_dir / "hub" / "data" / (r["id"] + ".js")
        assert side.is_file()                                  # the page of the analysis reads this, from file:// too
        head = side.read_text("utf-8")[:120]
        assert head.startswith("__FCB_SIDECAR__(") and r["id"] in head
    page = (lib_dir / "hub" / "index.html").read_text("utf-8")
    assert not any(r["id"] in page or (r["label"] or "\0") in page for r in rows)       # the page holds no data


def test_a_one_file_report_of_a_pair_opens_on_its_analysis(spatial_library, tmp_path):
    lib_dir = spatial_library["dir"]
    store = AnalysisStore(str(lib_dir))
    rows = {(r["pipeline"], r["assay"]): r for r in hub_index(lib_dir)["analyses"]}
    wf, sr = rows[("wf-single-cell", "visium_hd")], rows[("spaceranger", "visium_hd")]
    out = sc_report.write_report([store.composed(wf["id"]), store.composed(sr["id"])], str(tmp_path / "pair.flowcellar.html"))
    page = Path(out).read_text("utf-8")
    assert "fcb:" + wf["id"] in page and "fcb:" + sr["id"] in page
    assert '"analysis":"' + wf["id"] + '"' in page.replace("&quot;", '"').replace(" ", "")
    assert len(page) < 3_000_000


def test_hub_reports_says_that_analyses_get_no_report_file(spatial_library):
    """``hub --reports`` writes one report per run. The analyses of the library get none from it: the command says
    so, and names the command that writes the report of one analysis, instead of leaving them out silently."""
    from click.testing import CliRunner

    from flowcellar.cli import cli
    lib_dir = spatial_library["dir"]
    res = CliRunner().invoke(cli, ["hub", str(lib_dir), "--reports"])
    assert res.exit_code == 0, res.output
    assert "1 report written" in res.output
    assert "no report file was written for the 4 single-cell and spatial analyses" in res.output and "flowcellar analysis ANALYSIS_DIR -o OUT" in res.output
    assert len([n for n in os.listdir(lib_dir / "reports") if n.endswith(".flowcellar.html")]) == 1
    res = CliRunner().invoke(cli, ["hub", str(lib_dir)])
    assert res.exit_code == 0 and "no report file" not in res.output


@pytest.mark.slow
@pytest.mark.skipif(os.environ.get("FLOWCELLAR_BROWSER") != "1", reason="set FLOWCELLAR_BROWSER=1 to run the browser checks")
def test_the_analysis_screens_in_a_real_browser(spatial_library, tmp_path):
    node = shutil.which("node")
    if not node or not (REPO / "viewer" / "node_modules" / "playwright").is_dir():
        pytest.skip("Node or viewer/node_modules/playwright is missing")
    lib_dir = spatial_library["dir"]
    store = AnalysisStore(str(lib_dir))
    rows = {(r["pipeline"], r["assay"]): r for r in hub_index(lib_dir)["analyses"]}
    report = sc_report.write_report([store.composed(rows[("wf-single-cell", "visium_hd")]["id"]), store.composed(rows[("spaceranger", "visium_hd")]["id"])],
                                    str(tmp_path / "pair.flowcellar.html"))
    r = subprocess.run([node, str(REPO / "viewer" / "verify" / "analysis.mjs"), "--hub", str(lib_dir / "hub"), "--out", str(tmp_path / "shots"),
                        "--widths", "1366,2560", "--report", report], cwd=REPO / "viewer", capture_output=True, text=True, timeout=1800)
    if r.returncode == 2 and "No Playwright Firefox found" in r.stderr:
        pytest.skip("no Playwright Firefox on this machine")
    tail = "\n".join((r.stdout + r.stderr).splitlines()[-40:])
    assert r.returncode == 0, tail
