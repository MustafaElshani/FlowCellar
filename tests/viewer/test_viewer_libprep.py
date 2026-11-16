"""The Pychopper and RNA composition sections of the run page, in a real browser (viewer/verify/libprep.mjs).

Two fixture runs are made by the engine itself (a fake run folder, a sample sheet, complete pychopper statistics
tables and annotated read tables for 4 and for 24 samples), written as one-run reports, and opened in headless
Firefox at 1,366 and 2,560 px in both themes, with all samples and with three samples chosen::

    FLOWCELLAR_BROWSER=1 .venv/bin/python -m pytest tests/viewer/test_viewer_libprep.py -q

The first test needs no browser: it checks that the report of such a run carries the three sections.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from flowcellar.model import bundle as bundle_model
from flowcellar.testing import libprep as F
from flowcellar.testing import readdata as T
from flowcellar.testing.fixtures import RunSpec, make_run

REPO = Path(__file__).resolve().parents[2]
pytestmark = [pytest.mark.viewer]
NAMES = ("ctrl", "treated", "knockdown", "rescue")


def make_report(base: Path, n: int) -> tuple[Path, dict]:
    """A fixture run with ``n`` samples -> (the one-run report, its bundle)."""
    fr = make_run(base / "project", RunSpec(experiment="demo_experiment", sample=f"pool_{n}", seed=100 + n))
    ref = T.run_ref(fr)
    aliases = [f"{NAMES[i % 4]}_rep{i // 4 + 1}" for i in range(n)]
    rows = [{"experiment_id": ref.experiment, "kit": T.KIT, "flow_cell_id": ref.flow_cell_id, "barcode": f"barcode{i + 1:02d}", "sample_id": a, "alias": a}
            for i, a in enumerate(aliases)]
    sheet = T.write_sample_sheet(str(base / "project" / "samples.csv"), rows)
    ana = base / "project" / "analysis"
    for i, a in enumerate(aliases):
        spec = F.PychopperSpec(primers_found=40000 + 3100 * (i % 7), plus=19000 + 1700 * (i % 5), rescued={2: 60 + 9 * i, 3: 5 + i % 4, 6: i % 2},
                               unusable_hits={0: 300 + 40 * (i % 6), 2: 60, 3: 500 + 90 * (i % 4), 4: 20, 7: 2, 30: 1}, qc_fail=400 + 35 * i, len_fail=90 + 7 * i,
                               umi=(38000 + 2900 * (i % 7), 37500 + 2900 * (i % 7)), cutoff_index=1 if i != 2 else 0, seed=200 + i)
        F.write_pychopper_stats_full(str(ana / "pychopper_reports" / f"{a}_pychopper_stats.tsv"), spec)
        F.write_annotated_table(str(ana / "dorado_polyA" / "annotated_readid" / f"{a}_readid_annotated.csv"),
                                F.make_table_rows(3000 + 400 * (i % 5), seed=300 + i, not_primary=3))
    out = base / "out"
    env = {k: v for k, v in os.environ.items() if not k.startswith("FLOWCELLAR_")}
    r = subprocess.run([sys.executable, "-m", "flowcellar", "run", fr.run_dir, str(ana / "pychopper_reports"), str(ana / "dorado_polyA"),
                        "--sample-sheet", sheet, "-o", str(out)], capture_output=True, text=True, timeout=600, env=env)
    reports = sorted(out.glob("*.flowcellar.html")) if out.is_dir() else []
    assert reports, r.stdout + r.stderr
    return reports[0], bundle_model.read_bundle(str(next(out.glob("*.fcb.json.gz"))))


def test_report_of_a_run_with_library_preparation_tables(tmp_path, viewer_built):
    report, bundle = make_report(tmp_path / "four", 4)
    mods = bundle["modules"]
    assert {"cdna.pychopper", "rna.biotype", "rna.genes"} <= set(mods)
    assert len(mods["cdna.pychopper"]["cells"]) == 4 and mods["rna.biotype"]["summary"]["rna_biotype_samples"] == 4
    page = report.read_text("utf-8")
    assert "fcb:" + bundle["bundle_id"] in page and len(page) < 3_000_000


@pytest.mark.slow
@pytest.mark.skipif(os.environ.get("FLOWCELLAR_BROWSER") != "1", reason="set FLOWCELLAR_BROWSER=1 to run the browser checks")
def test_sections_in_a_real_browser_for_4_and_24_samples(tmp_path, viewer_built):
    node = shutil.which("node")
    if not node or not (REPO / "viewer" / "node_modules" / "playwright").is_dir():
        pytest.skip("Node or viewer/node_modules/playwright is missing")
    reports = [str(make_report(tmp_path / f"n{n}", n)[0]) for n in (4, 24)]
    r = subprocess.run([node, str(REPO / "viewer" / "verify" / "libprep.mjs"), *reports], cwd=REPO / "viewer", capture_output=True, text=True, timeout=1800)
    if r.returncode == 2 and "No Playwright Firefox found" in r.stderr:
        pytest.skip("no Playwright Firefox on this machine")
    tail = "\n".join((r.stdout + r.stderr).splitlines()[-40:])
    assert r.returncode == 0, tail
