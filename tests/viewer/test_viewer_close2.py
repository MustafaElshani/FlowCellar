"""Regressions held as checks in a real browser (viewer/verify/close2.mjs).

One library is made by the engine itself (invented names only) from two trees:

* a project whose three samples were sequenced on **two flow cells**: one alignment file, one pychopper statistics
  table and one annotated read table per sample, each covering both flow cells; the tables are read, the alignment
  files are left waiting (``scan --read tables``), so the hub holds a table shared by two runs, sections whose files
  wait, and samples that only the tables have numbers for;
* the spatial tree of ``flowcellar.testing.sc`` (a Visium and a Visium HD sample, each analysed by both pipelines).

Its hub is opened in headless Firefox::

    FLOWCELLAR_BROWSER=1 .venv/bin/python -m pytest tests/viewer/test_viewer_close2.py -q

The first test needs no browser: it checks that the bundles carry what the page must say.
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

pytest.importorskip("pysam", reason="the fixtures write BAM files with pysam (the 'bam' extra)")

from flowcellar.library import scan as scan_mod  # noqa: E402
from flowcellar.model import bundle as bundle_model  # noqa: E402
from flowcellar.readdata import link as link_mod  # noqa: E402
from flowcellar.testing import libprep as F  # noqa: E402
from flowcellar.testing import readdata as T  # noqa: E402
from flowcellar.testing.fixtures import RunSpec, make_run  # noqa: E402
from flowcellar.testing.trees import age_tree  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
pytestmark = [pytest.mark.viewer]
# names that differ in the middle only, as the names of one experiment do: a cut in the middle would make them equal
ALIASES = ("demo_batch_101_plate_ctrl_day1", "demo_batch_101_plate_trt_day1", "demo_batch_102_plate_ctrl_day1")
RUN_READS = 2500

CONFIG = """\
[library]
name = "closing demo"

[discovery]
min_age_seconds = 0

[[inputs]]
role = "align.genome_bam"
glob = "**/analysis/alignments/{sample}_genome_sorted.bam"
[[inputs]]
role = "pychopper.stats"
glob = "**/analysis/pychopper_reports/{sample}_pychopper_stats.tsv"
[[inputs]]
role = "polya.read_table"
glob = "**/analysis/dorado_polyA/annotated_readid/{sample}_readid_annotated.csv"
[[inputs]]
role = "basecaller.sample_sheet"
sniff = "csv_columns:barcode,alias,flow_cell_id|position_id"

[[link.rules]]
name = "experiment-folder"
roles = ["pychopper.stats", "polya.read_table"]
via = "sample"
within = "experiment"
"""


def make_shared_project(root: Path) -> dict:
    """``root/project_shared``: two run folders, a sample sheet for each flow cell, and under ``analysis/`` one
    alignment file, one pychopper table and one annotated read table per sample, each for both flow cells."""
    proj = root / "project_shared"
    fa = make_run(proj, RunSpec(flow_cell_id="TST00021", experiment="demo_shared", sample="pool_day1", reads=RUN_READS, seed=71))
    fb = make_run(proj, RunSpec(flow_cell_id="TST00022", experiment="demo_shared", sample="pool_day2", reads=RUN_READS, seed=72, start="2026-03-04T09:00:00Z"))
    A, B = T.run_ref(fa), T.run_ref(fb)
    for ref, name in ((A, "sheet_day1.csv"), (B, "sheet_day2.csv")):
        rows = [{"experiment_id": ref.experiment, "kit": T.KIT, "flow_cell_id": ref.flow_cell_id, "barcode": f"barcode{i + 1:02d}", "sample_id": a, "alias": a}
                for i, a in enumerate(ALIASES)]
        T.write_sample_sheet(str(proj / name), rows)
    ana = proj / "analysis"
    (ana / "alignments").mkdir(parents=True)
    for i, a in enumerate(ALIASES):
        reads = T.make_reads(260 + 20 * i, A, seed=10 + i, alias=a, barcode=f"barcode{i + 1:02d}") + T.make_reads(200 + 20 * i, B, seed=20 + i, alias=a, barcode=f"barcode{i + 1:02d}")
        T.write_aligned_bam(str(ana / "alignments" / f"{a}_genome_sorted.bam"), reads, seed=30 + i)
        spec = F.PychopperSpec(primers_found=3000 + 400 * i, plus=1400 + 200 * i, rescued={2: 30 + 5 * i, 3: 4}, unusable_hits={0: 80 + 10 * i, 2: 20, 3: 160, 4: 6},
                               qc_fail=120 + 15 * i, len_fail=40 + 5 * i, umi=(2800 + 380 * i, 2750 + 380 * i), cutoff_index=1, seed=40 + i)
        F.write_pychopper_stats_full(str(ana / "pychopper_reports" / f"{a}_pychopper_stats.tsv"), spec)
        # more rows than one flow cell sequenced: the total of the section can only be said as "both flow cells"
        F.write_annotated_table(str(ana / "dorado_polyA" / "annotated_readid" / f"{a}_readid_annotated.csv"), F.make_table_rows(1100 + 150 * i, seed=50 + i, not_primary=3))
    return {"root": str(root), "runs": [fa.run_key, fb.run_key]}


@pytest.fixture(scope="module")
def closing_library(tmp_path_factory, viewer_built):
    before = os.environ.get(link_mod.INPROCESS_ENV)
    os.environ[link_mod.INPROCESS_ENV] = "1"
    try:
        base = tmp_path_factory.mktemp("viewer_close2")
        shared = make_shared_project(base / "store_lr")
        roots = [shared["root"]]
        try:
            import h5py  # noqa: F401

            from flowcellar.testing import sc as fx
            tree = fx.make_spatial_tree(str(base / "spatial"))
            roots += [tree.roots["seq"], tree.roots["store"]]
        except ImportError:                       # without the 'sc' extra the library holds the runs only
            pass
        for r in roots:
            age_tree(r)
        lib = base / "library"
        lib.mkdir()
        (lib / "flowcellar.toml").write_text(CONFIG)
        res = scan_mod.scan(str(lib), scan_mod.ScanOptions(roots=tuple(roots), wait=False, read=("tables", "analyses")))
        assert res.exit_code == 0
        return {"dir": lib, "runs": shared["runs"]}
    finally:
        if before is None:
            os.environ.pop(link_mod.INPROCESS_ENV, None)
        else:
            os.environ[link_mod.INPROCESS_ENV] = before


def hub_index(lib_dir: Path) -> dict:
    m = re.fullmatch(r'__FC_INDEX__\("([A-Za-z0-9+/=]*)"\);\n', (lib_dir / "hub" / "index.js").read_text("utf-8"))
    return json.loads(gzip.decompress(base64.b64decode(m.group(1))))


def test_the_bundles_carry_what_the_page_must_say(closing_library):
    lib = closing_library["dir"]
    totals = []
    for key in closing_library["runs"]:
        b = bundle_model.read_bundle(str(lib / "bundles" / bundle_model.bundle_filename(key)))
        mods = b["modules"]
        assert [s["name"] for s in b["samples"]] == list(ALIASES), "the samples of the sheet, and no second sample named after a file"
        for sec in ("rna.biotype", "cdna.pychopper"):
            notes = [n for n in mods[sec]["notes"] if n["code"] == "covers_several_runs"]
            assert len(notes) == len(ALIASES) and all(n["value"] == 2 and "merged before" in n["message"] for n in notes), sec
        # the alignment files hold the reads of both flow cells and wait: no cell, and the reason is "waiting"
        for sec in ("len", "aln"):
            assert mods[sec]["state"] == "partial" and mods[sec]["partial_reason"] == "waiting" and mods[sec]["cells"] == [], sec
        table_reads = sum(c["scalars"]["reads"] for c in mods["rna.biotype"]["cells"])
        totals.append(table_reads)
        assert table_reads > mods["yield"]["summary"]["yield_reads"], "the fixture must hold more annotated reads than one flow cell sequenced"
    assert totals[0] == totals[1], "the same table stands on both runs: its total is the reads of both flow cells"


@pytest.mark.slow
@pytest.mark.skipif(os.environ.get("FLOWCELLAR_BROWSER") != "1", reason="set FLOWCELLAR_BROWSER=1 to run the browser checks")
def test_the_findings_of_the_second_test_round_in_a_real_browser(closing_library, tmp_path):
    node = shutil.which("node")
    if not node or not (REPO / "viewer" / "node_modules" / "playwright").is_dir():
        pytest.skip("Node or viewer/node_modules/playwright is missing")
    lib = closing_library["dir"]
    r = subprocess.run([node, str(REPO / "viewer" / "verify" / "close2.mjs"), "--hub", str(lib / "hub"), "--out", str(tmp_path / "shots"),
                        "--need", "shared,waiting,control,names"], cwd=REPO / "viewer", capture_output=True, text=True, timeout=1800)
    if r.returncode == 2 and "No Playwright Firefox found" in r.stderr:
        pytest.skip("no Playwright Firefox on this machine")
    tail = "\n".join((r.stdout + r.stderr).splitlines()[-40:])
    assert r.returncode == 0, tail
    rep = json.loads((tmp_path / "shots" / "close2_report.json").read_text())
    cases = rep["summary"]["cases"]
    assert rep["summary"]["failed"] == 0 and rep["summary"]["errors"] == 0 and rep["summary"]["requests_to_another_origin"] == 0
    assert cases.get("shared", 0) >= 4 and cases.get("waiting", 0) >= 4 and cases.get("control", 0) >= 2 and cases.get("names", 0) >= 2, cases
