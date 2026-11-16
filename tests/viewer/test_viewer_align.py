"""The alignment, over-time and transcript charts of the run page, in a real browser (viewer/verify/align.mjs).

Two fixture runs are made by the engine itself (a fake run folder, a sample sheet, Dorado BAMs, genome BAMs and
transcriptome BAMs with known content for 3 and for 13 samples), written as one-run reports, and opened in headless
Firefox at 1,366 and 2,560 px in both themes, with all samples, with one sample and with two samples chosen::

    FLOWCELLAR_BROWSER=1 .venv/bin/python -m pytest tests/viewer/test_viewer_align.py -q

The first test needs no browser: it checks that the report of such a run carries the sections and their charts' data.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from flowcellar.model import bundle as bundle_model
from flowcellar.model import schema
from flowcellar.testing import readdata as T
from flowcellar.testing.fixtures import RunSpec, make_run

REPO = Path(__file__).resolve().parents[2]
pytestmark = [pytest.mark.viewer]
pytest.importorskip("pysam", reason="the fixture writes BAM files (extra 'bam')")


def make_report(base: Path, n: int, reads: int = 500) -> tuple[Path, dict]:
    """A fixture run with ``n`` samples, each with a Dorado BAM, a genome BAM and a transcriptome BAM -> (report, bundle)."""
    fr = make_run(base / "project", RunSpec(experiment="demo_experiment", sample=f"pool_{n}", seed=300 + n))
    ref = T.run_ref(fr)
    aliases = tuple(f"sample_{chr(97 + i)}" for i in range(n))
    proj = T.make_project(str(base / "project"), ref, aliases=aliases, reads_per_sample=reads, seed=61, quant=False)
    args = []
    for i, alias in enumerate(aliases):
        rs = T.make_reads(reads, ref, seed=80 + i, alias=alias, barcode=proj.barcodes[i])
        p = base / "project" / "analysis" / "alignments" / f"{alias}_transcriptome_sorted.bam"
        T.write_transcriptome_bam(str(p), rs, seed=90 + i)
        args += ["--input", f"align.transcriptome_bam={p}"]
    out = base / "out"
    env = {k: v for k, v in os.environ.items() if not k.startswith("FLOWCELLAR_")}
    r = subprocess.run([sys.executable, "-m", "flowcellar", "run", fr.run_dir, str(base / "project"), *args, "--sample-sheet", proj.sample_sheet, "-o", str(out)],
                       capture_output=True, text=True, timeout=900, env=env)
    reports = sorted(out.glob("*.flowcellar.html")) if out.is_dir() else []
    assert reports, r.stdout + r.stderr
    return reports[0], bundle_model.read_bundle(str(next(out.glob("*.fcb.json.gz"))))


def test_report_of_a_run_with_alignments_carries_the_new_sections(tmp_path, viewer_built):
    report, bundle = make_report(tmp_path / "three", 3)
    mods = bundle["modules"]
    assert {"aln", "aln.time", "aln.windows", "dyn.basecalled", "rna.polya", "rna.tx", "rna.coverage"} <= set(mods)
    assert schema.schema_errors(bundle) == []
    assert len(mods["aln.windows"]["cells"]) == 3 and len(mods["rna.coverage"]["cells"]) == 12
    clock = mods["dyn.basecalled"]["params"]["clock"]
    assert clock["origin_known"] and -600 < clock["offset_600s_s"] <= 0 and -3600 < clock["offset_3600s_s"] <= 0
    assert {"clock_600s", "clock_3600s", "genome_window"} <= set(bundle["axes"])
    page = report.read_text("utf-8")
    assert "fcb:" + bundle["bundle_id"] in page and len(page) < 3_000_000


@pytest.mark.slow
@pytest.mark.skipif(os.environ.get("FLOWCELLAR_BROWSER") != "1", reason="set FLOWCELLAR_BROWSER=1 to run the browser checks")
def test_charts_in_a_real_browser_for_3_and_13_samples(tmp_path, viewer_built):
    node = shutil.which("node")
    if not node or not (REPO / "viewer" / "node_modules" / "playwright").is_dir():
        pytest.skip("Node or viewer/node_modules/playwright is missing")
    reports = [str(make_report(tmp_path / f"n{n}", n, reads)[0]) for n, reads in ((3, 1500), (13, 200))]
    r = subprocess.run([node, str(REPO / "viewer" / "verify" / "align.mjs"), *reports], cwd=REPO / "viewer", capture_output=True, text=True, timeout=2400)
    if r.returncode == 2 and "No Playwright Firefox found" in r.stderr:
        pytest.skip("no Playwright Firefox on this machine")
    tail = "\n".join((r.stdout + r.stderr).splitlines()[-60:])
    assert r.returncode == 0, tail
