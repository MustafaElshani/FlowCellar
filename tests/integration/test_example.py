"""The example of the README (``python -m flowcellar.testing.demo``), written small and scanned once.

What the README and ``docs/quick-start.md`` say about the example is held here: how many runs, flow cells and
projects, which verdicts, which runs have read data, the sample that went on two flow cells, the two basecalling
passes, the analyses and the flow cells they are linked to. A change of the generator or of the engine that makes
one of those sentences untrue fails a test.
"""
from __future__ import annotations

import gzip
import json
import os
import re
import shutil
import subprocess
import sys
import time
import tomllib
from pathlib import Path

import pytest

from flowcellar.config import loader
from flowcellar.model import bundle as bundle_model
from flowcellar.testing import demo

REPO = Path(__file__).resolve().parents[2]
SCALE = 0.002
pytestmark = pytest.mark.slow


def flowcellar(*args: str, cwd: str | None = None) -> subprocess.CompletedProcess:
    env = {k: v for k, v in os.environ.items() if not k.startswith("FLOWCELLAR_")}
    return subprocess.run([sys.executable, "-m", "flowcellar", *args], capture_output=True, text=True, timeout=900, env=env, cwd=cwd)


@pytest.fixture(scope="module")
def example(tmp_path_factory):
    pytest.importorskip("pysam", reason="the example writes BAM files (extra 'bam')")
    pytest.importorskip("h5py", reason="the example writes Space Ranger folders (extra 'sc')")
    base = tmp_path_factory.mktemp("example")
    tree = demo.make_demo(base / "demo", scale=SCALE, jobs=4)
    before = {p: os.stat(p).st_mtime_ns for d, _s, files in os.walk(tree.roots["minknow"]) for p in (os.path.join(d, f) for f in files)}
    scan = flowcellar("scan", "demo/minknow", "demo/projects", "-L", "demo/library", cwd=str(base))
    lib = Path(tree.library)
    with gzip.open(lib / "index" / "index.json.gz") as fh:
        index = json.load(fh)
    cols = index["runs"]["columns"]

    def col(name: str) -> list:
        c = cols[name]
        return [c["dict"][i] for i in c["codes"]] if isinstance(c, dict) else list(c)
    yield {"base": base, "tree": tree, "lib": lib, "scan": scan, "index": index, "col": col, "before": before}
    shutil.rmtree(base, ignore_errors=True)


def _bundle(example, run_key: str) -> dict:
    return bundle_model.read_bundle(str(example["lib"] / "bundles" / bundle_model.bundle_filename(run_key)))


def _status(example, run_key: str) -> dict:
    return json.loads((example["lib"] / "status" / f"{run_key}.status.json").read_text("utf-8"))["status"]


def test_the_example_is_what_its_description_says(example):
    tree, runs = example["tree"], demo.demo_runs()
    assert len(runs) == 15 and len({r.flow_cell_id for r in runs}) == 13 and len({r.experiment for r in runs}) == 4
    assert all(re.fullmatch(r"TST\d{5}", r.flow_cell_id) for r in runs), "invented flow cell ids only"
    assert set(tree.runs) == {r.label for r in runs} and not tree.skipped
    assert sorted(os.listdir(tree.library)) == ["flowcellar.toml"] or (Path(tree.library) / "LIBRARY").exists()   # before the scan: the configuration only
    # the configuration it writes is accepted as it is, and says that its output limits are those of small runs
    text = demo.config_text(SCALE)
    loaded = loader.merge(tomllib.loads(text), source="the example")
    assert loaded.warnings == [] and loaded.notices == []
    assert "0.2 %" in text and "[profiles.cdna.rules.yield_basecalled_bases]" in text and "[profiles." not in demo.config_text(1.0)
    # equal arguments give the same runs: the names the README types are those of every copy of the example
    assert demo.run_key({r.label: r for r in runs}["g3"]) == "TST00202_5cf5f173"
    readme = (REPO / "README.md").read_text("utf-8")
    typed = set(re.findall(r"\bdemo/(?:minknow|projects)/[\w./-]+", readme))
    assert len(typed) >= 3
    for path in typed:
        assert (example["base"] / path).exists(), f"README names {path}, which the example does not hold"
    # nothing in it is newer than two minutes: the first scan reads all of it
    assert max(example["before"].values()) / 1e9 < time.time() - 300


def test_the_first_scan_gives_the_library_the_readme_shows(example):
    scan, col, index = example["scan"], example["col"], example["index"]
    assert scan.returncode == 0, scan.stdout[-2000:] + scan.stderr[-2000:]
    said = scan.stdout + scan.stderr
    assert "found 15 runs (15 with MinKNOW reports, 0 with raw data only)" in said
    assert "done: 15 runs · 12 Pass, 2 Warning, 0 Fail, 1 Incomplete" in said and "83 of 83 files read" in said
    assert said.count("analysis  visium_") >= 4 and "hub ready: " in said
    assert index["runs"]["n"] == 15 and sorted(p["id"] for p in index["projects"] if p["id"] != "_none") == [
        "demo_cdna_timecourse", "demo_direct_rna", "demo_genome", "demo_spatial"]
    verdict = dict(zip(col("id"), col("verdict")))
    by_fc = {}
    for key, fc in zip(col("id"), col("fc")):
        by_fc.setdefault(fc, []).append(key)
    assert sorted(v for v in verdict.values()) == ["incomplete"] + ["pass"] * 12 + ["warn"] * 2
    assert {fc for fc, keys in by_fc.items() if any(verdict[k] == "warn" for k in keys)} == {"TST00104", "TST00403"}
    # nothing was written where the data is
    after = {p: os.stat(p).st_mtime_ns for p in example["before"]}
    assert after == example["before"] and not [d for d in os.listdir(example["tree"].roots["minknow"]) if d.startswith(".")]
    # the run that was stopped and started again; the flow cell that was used twice
    use, restart = dict(zip(col("id"), col("use"))), dict(zip(col("id"), col("restart")))
    assert len(by_fc["TST00103"]) == 2 and sorted(verdict[k] for k in by_fc["TST00103"]) == ["incomplete", "pass"]
    assert sum(1 for k in by_fc["TST00103"] if restart[k]) == 1
    assert sorted(use[k] for k in by_fc["TST00201"]) == [1, 2]
    # read data: three cDNA runs, one genomic, one direct RNA; the others are judged on run health alone
    data = dict(zip(col("fc"), col("read_data")))
    assert {fc for fc, state in data.items() if state == "complete"} == {"TST00105", "TST00106", "TST00107", "TST00202", "TST00301"}
    # two folders were scanned, and one project lies under both
    assert len(index["roots"]) == 2
    cdna = next(p for p in index["projects"] if p["id"] == "demo_cdna_timecourse")
    assert len(cdna.get("roots") or cdna.get("root_ids") or []) == 2, cdna


def test_a_sample_on_two_flow_cells_has_a_page_per_flow_cell_and_one_merged_analysis(example):
    first, second = _bundle(example, "TST00105_8f6031c7"), _bundle(example, "TST00106_1ae7d84b")
    for b in (first, second):
        named = [(s["name"], s.get("barcode")) for s in b["samples"] if s.get("barcode")]
        assert named == [(n, f"barcode{i + 1:02d}") for i, n in enumerate(demo.SAMPLES)], "names and barcodes from the sample sheet of the basecaller"
        assert any(n.get("code") == "covers_several_runs" for n in b["modules"]["cdna.pychopper"]["notes"])
    # each flow cell has its own reads; the pychopper analysis is one, the same on both pages
    reads = [b["modules"]["yield.basecalled"]["summary"] for b in (first, second)]
    assert reads[0] != reads[1]
    assert first["modules"]["cdna.pychopper"]["summary"] == second["modules"]["cdna.pychopper"]["summary"]
    # the pool that went on one flow cell has every stage, so its usable share is known
    usable = next(r for r in _status(example, "TST00107_0dfa2105")["rules"] if r["id"] == "funnel_usable_frac")
    assert usable["outcome"] == "ok" and 0.55 < usable["value"] < 0.85


def test_of_two_basecalling_passes_one_is_counted_and_the_other_is_listed(example):
    b = _bundle(example, "TST00202_5cf5f173")
    passes = {i["variant"]: i for i in b["provenance"]["inputs"] if i["role"] == "dorado.bam"}
    assert set(passes) == {"sup", "sup_5mC"}
    assert passes["sup"]["counted_for"] == ["reads"] and passes["sup"]["state"] == "ok"
    assert passes["sup_5mC"]["counted_for"] == [] and passes["sup_5mC"]["state"] == "skipped"
    assert [s["name"] for s in b["samples"]] == ["isolate_3"]


def test_analyses_are_linked_to_every_flow_cell_their_reads_came_from(example):
    rows = {(a["label"], a["pipeline"]): a for a in example["index"]["analyses"]}
    assert set(rows) == {("visium_a", "spaceranger"), ("visium_a", "wf-single-cell"), ("visium_hd_b", "spaceranger"), ("visium_hd_b", "wf-single-cell")}
    assert all(a["verdict"] == "pass" for a in rows.values()), {k: (a["verdict"], a["headline"]) for k, a in rows.items()}
    # the Visium HD slide: two flow cells, told by the reads of the workflow's BAM, with unequal shares
    for pipeline, method in (("wf-single-cell", "read_tag"), ("spaceranger", "sibling_analysis")):
        links = rows[("visium_hd_b", pipeline)]["links"]
        assert [(ln["fc"], ln["method"], ln["resolved"]) for ln in links] == [("TST00402", method, True), ("TST00403", method, True)]
        assert 0.55 < links[0]["share"] < 0.70 and abs(links[0]["share"] + links[1]["share"] - 1.0) < 1e-6
    # the Visium slide: the files name no run; the sample name does, because one run of the library carries it
    links = rows[("visium_a", "spaceranger")]["links"]
    assert [(ln["fc"], ln["method"]) for ln in links] == [("TST00401", "sample_name_chain")]
    # an analysis never changes the verdict of a flow cell: the flow cell with few pores stays a warning beside a passing analysis
    assert dict(zip(example["col"]("fc"), example["col"]("verdict")))["TST00403"] == "warn"


def test_one_run_without_a_library_as_the_readme_types_it(example):
    out = example["base"] / "qc"
    r = flowcellar("--config", "demo/library/flowcellar.toml", "run",
                   "demo/minknow/demo_genome/isolate_3/20261013_1000_MD-000001_TST00202_5cf5f173", "demo/projects/demo_genome", "-o", "qc",
                   cwd=str(example["base"]))
    assert r.returncode == 0, r.stdout + r.stderr
    assert (out / "TST00202_5cf5f173.fcb.json.gz").is_file() and (out / "TST00202_5cf5f173.flowcellar.html").is_file()
    said = r.stdout + r.stderr
    assert re.search(r"Run health\s+Pass", said) and re.search(r"Read data\s+Pass", said), said[-800:]
