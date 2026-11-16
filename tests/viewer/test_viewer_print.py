"""Save as PDF, in real browsers (viewer/verify/print.mjs).

A run page with every section (3 and 13 samples), an analysis page with its spatial map, a comparison of three runs
and the history of a flow cell are made by the engine from invented data (tests/viewer/viewer_print_fixture.py), laid
out for paper in Chromium and in Firefox, and saved as PDFs on A4 and on Letter: Chromium through Playwright, Firefox
through its own print pipeline. The check holds what a PDF must be (the head and the selection line, no control, no
blank sheet, no heading alone at the foot of a sheet, charts painted again at the width of the sheet, under 5 MB) and
that the screen is afterwards exactly as it was::

    FLOWCELLAR_BROWSER=1 .venv/bin/python -m pytest tests/viewer/test_viewer_print.py -q

It needs a Playwright Chromium and a Playwright Firefox (PLAYWRIGHT_BROWSERS_PATH or FLOWCELLAR_PW_BROWSERS names
their folder; ``cd viewer && PLAYWRIGHT_BROWSERS_PATH=pw-browsers npx playwright install chromium firefox``) and
poppler's ``pdfinfo``, ``pdftotext``, ``pdfimages`` and ``pdftoppm``. The first tests need no browser.
"""
from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
pytestmark = [pytest.mark.viewer]
POPPLER = ("pdfinfo", "pdftotext", "pdfimages", "pdftoppm")


def test_the_built_viewer_carries_the_control_and_the_layout_for_paper(viewer_built):
    dist = REPO / "src" / "flowcellar" / "viewer_dist"
    script, styles = (dist / "viewer.core.js").read_text("utf-8"), (dist / "viewer.css").read_text("utf-8")
    assert "Save as PDF" in script and "PDF with every chart" in script
    assert "beforeprint" in script and "afterprint" in script
    assert "html.fc-print" in styles and "@page" in styles and "--fc-paper-w" in styles
    # the layout for paper hangs on a class that the script sets, never on the width of the window: no rule of it
    # stands inside a media query (the one `@media print` block after it is the fall-back for a page printed without the script)
    paper = re.sub(r"/\*.*?\*/", "", styles[styles.index("/* ==== Paper"):], flags=re.S)
    assert "@media" not in paper.split("@media print")[0] and paper.count("@media") == 1


def test_the_fixture_run_has_every_section_a_run_page_can_show(tmp_path, viewer_built):
    pytest.importorskip("pysam", reason="the fixture writes BAM files (extra 'bam')")
    sys.path.insert(0, str(Path(__file__).parent))
    import viewer_print_fixture

    from flowcellar.model import bundle as bundle_model

    report = Path(viewer_print_fixture.make_run_report(tmp_path / "run3", 3, 300))
    bundle = bundle_model.read_bundle(str(next(report.parent.glob("*.fcb.json.gz"))))
    assert {"funnel", "len", "qual", "aln", "aln.time", "dyn.basecalled", "rna.tx", "rna.polya", "rna.coverage", "cdna.pychopper", "rna.biotype", "rna.genes"} <= set(bundle["modules"])
    assert len([s for s in bundle["samples"] if s["key"] != "unclassified"]) == 3
    page = report.read_text("utf-8")
    assert "fcb:" + bundle["bundle_id"] in page and len(page) < 3_000_000


@pytest.mark.slow
@pytest.mark.skipif(os.environ.get("FLOWCELLAR_BROWSER") != "1", reason="set FLOWCELLAR_BROWSER=1 to run the browser checks")
def test_pages_go_to_paper_and_come_back_in_chromium_and_firefox(tmp_path, viewer_built):
    pytest.importorskip("pysam", reason="the fixture writes BAM files (extra 'bam')")
    node = shutil.which("node")
    if not node or not (REPO / "viewer" / "node_modules" / "playwright").is_dir():
        pytest.skip("Node or viewer/node_modules/playwright is missing")
    missing = [tool for tool in POPPLER if not shutil.which(tool)]
    if missing:
        pytest.skip("poppler's " + ", ".join(missing) + " not found: the PDFs cannot be read")
    sys.path.insert(0, str(Path(__file__).parent))
    import viewer_print_fixture

    made = viewer_print_fixture.make(str(tmp_path / "fx"))
    fixture = tmp_path / "fixture.json"
    fixture.write_text(json.dumps(made), "utf-8")
    env = {k: v for k, v in os.environ.items() if k != "FLOWCELLAR_BROWSER_ENGINE"}
    r = subprocess.run([node, str(REPO / "viewer" / "verify" / "print.mjs"), "--fixture", str(fixture), "--out", str(tmp_path / "out")], cwd=REPO / "viewer",
                       capture_output=True, text=True, timeout=3000, env=env)
    if r.returncode == 2 and "No Playwright Chromium found" in r.stderr:
        pytest.skip("no Playwright Chromium on this machine (cd viewer && PLAYWRIGHT_BROWSERS_PATH=pw-browsers npx playwright install chromium)")
    if r.returncode == 2 and "No Playwright Firefox found" in r.stderr:
        pytest.skip("no Playwright Firefox on this machine")
    tail = "\n".join((r.stdout + r.stderr).splitlines()[-60:])
    assert r.returncode == 0, tail
