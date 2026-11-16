"""Analyses in a library: found by a scan, stored as bundles of their own, linked to the runs, listed in the index."""
from __future__ import annotations

import gzip
import json
import os

import pytest
from click.testing import CliRunner

from flowcellar.cli import cli
from flowcellar.library import scan as scan_mod
from flowcellar.library.maintain import validate_paths
from flowcellar.library.store import LibraryStore
from flowcellar.sc import collect as collect_mod
from flowcellar.sc import store as sc_store
from flowcellar.testing import sc as fx
from flowcellar.util.errors import FlowCellarError


def run_scan(lib_dir: str, *roots: str, **options):
    options.setdefault("wait", False)
    return scan_mod.scan(lib_dir, scan_mod.ScanOptions(roots=tuple(roots), **options))


def index_of(lib_dir: str) -> dict:
    with gzip.open(os.path.join(lib_dir, "index", "index.json.gz")) as fh:
        return json.load(fh)


def by_pipeline(index: dict) -> dict:
    out = {}
    for a in index["analyses"]:
        out[(a["pipeline"], a["assay"])] = a
    return out


@pytest.fixture(scope="module")
def lib(tmp_path_factory):
    """A library scanned once from the two roots of the spatial tree (the run on the second flow cell is not in it)."""
    from flowcellar.report import api as report_api
    try:
        report_api.viewer_manifest()
    except FlowCellarError as e:
        pytest.skip(f"the viewer is not built: {e.message}")
    base = tmp_path_factory.mktemp("spatial")
    tree = fx.make_spatial_tree(str(base / "data"))
    lib_dir = str(base / "library")
    res = run_scan(lib_dir, tree.roots["seq"], tree.roots["store"])
    return {"tree": tree, "dir": lib_dir, "first": res, "index": index_of(lib_dir)}


def test_a_scan_finds_the_analyses_and_stores_a_bundle_for_each(lib):
    idx, tree = lib["index"], lib["tree"]
    assert lib["first"].exit_code == 0 and idx["runs"]["n"] == 1
    rows = idx["analyses"]
    assert len(rows) == 4                                   # the copy of the Space Ranger folder is the same analysis
    assert {(r["pipeline"], r["assay"]) for r in rows} == {("spaceranger", "visium"), ("wf-single-cell", "visium"),
                                                           ("spaceranger", "visium_hd"), ("wf-single-cell", "visium_hd")}
    store = sc_store.AnalysisStore(lib["dir"])
    assert sorted(store.ids()) == sorted(r["id"] for r in rows)
    for r in rows:
        b = store.read_measured(r["id"])
        assert b["kind"] == "analysis" and b["bundle_id"] == r["id"] and "status" not in b      # measured content only
        assert store.read_status_file(r["id"])["status"]["verdict"] == r["verdict"]
        assert len(json.dumps(r)) < 3000                     # scalars only: no map, no curve
    hd = by_pipeline(idx)[("spaceranger", "visium_hd")]
    assert hd["copies"] == 2 and hd["spatial_units_under_tissue"] == tree.analyses["sr_hd"].expected["summary"]["spatial_units_under_tissue"]
    assert LibraryStore.open(lib["dir"]).run_keys() == [tree.runs["a"].run_key]       # the run code does not see the analyses


def test_links_resolve_against_the_library_and_a_missing_run_stays_named(lib):
    idx, tree = lib["index"], lib["tree"]
    rows = by_pipeline(idx)
    key_a, key_b = tree.runs["a"].run_key, tree.runs["b"].run_key
    wf = rows[("wf-single-cell", "visium_hd")]
    assert [(ln["run"], ln["resolved"], ln["method"]) for ln in wf["links"]] == [(key_a, True, "read_tag"), (key_b, False, "read_tag")]
    assert [ln["share"] for ln in wf["links"]] == pytest.approx([0.6, 0.4], abs=0.02)
    # the Space Ranger output of the same reads has no identifier of its own: it takes the links of its sibling
    sr = rows[("spaceranger", "visium_hd")]
    assert [(ln["run"], ln["resolved"], ln["method"]) for ln in sr["links"]] == [(key_a, True, "sibling_analysis"), (key_b, False, "sibling_analysis")]
    assert sr["siblings"] == [wf["id"]] and wf["siblings"] == [sr["id"]] and wf["group"] == sr["group"] == "H1-TESTSLD-A1"
    # the older workflow keeps no identifier: the run folder beside the analysis with the same sample name
    wf2 = rows[("wf-single-cell", "visium")]
    assert [(ln["run"], ln["method"], ln["confidence"], ln["resolved"]) for ln in wf2["links"]] == [(key_a, "path_rule", "medium", True)]
    sr2 = rows[("spaceranger", "visium")]
    assert [(ln["run"], ln["method"]) for ln in sr2["links"]] == [(key_a, "sibling_analysis")] and sr2["group"] == wf2["group"]
    b = sc_store.AnalysisStore(lib["dir"]).composed(wf["id"])
    links = b["identity"]["analysis"]["links"]
    assert links[1]["run_key"] is None and links[1]["flow_cell_id"] == "TST00077" and links[1]["resolved"] is False
    rule = next(r for r in b["status"]["rules"] if r["id"] == "sc_link_runs")
    assert rule["severity"] == "info" and "TST00077" in rule["message"] and b["status"]["verdict"] == "pass"


def test_bundles_of_a_library_are_valid_with_their_verdicts(lib):
    results = validate_paths([lib["dir"]])
    assert [r for r in results if not r["ok"]] == [] and len(results) == 5       # one run and four analyses
    r = CliRunner().invoke(cli, ["validate", "-L", lib["dir"]])
    assert r.exit_code == 0 and r.output.count("ok  ") == 5


def test_a_second_scan_opens_no_file_of_an_analysis(lib, monkeypatch):
    calls = []
    monkeypatch.setattr(collect_mod, "collect", lambda *a, **k: calls.append("collect") or (_ for _ in ()).throw(AssertionError("read again")))
    monkeypatch.setattr(collect_mod, "quick_id", lambda *a, **k: calls.append("id") or (_ for _ in ()).throw(AssertionError("opened again")))
    store = sc_store.AnalysisStore(lib["dir"])
    before = {aid: os.stat(store.bundle_path(aid)).st_mtime_ns for aid in store.ids()}
    res = run_scan(lib["dir"])
    assert res.exit_code == 0 and calls == []
    assert {aid: os.stat(store.bundle_path(aid)).st_mtime_ns for aid in store.ids()} == before
    assert [a["id"] for a in index_of(lib["dir"])["analyses"]] == [a["id"] for a in lib["index"]["analyses"]]


def test_ls_and_show_work_for_analyses(lib):
    runner = CliRunner()
    r = runner.invoke(cli, ["ls", "-L", lib["dir"], "--analyses", "--format", "json", "--columns", "id,label,assay,verdict"])
    rows = json.loads(r.output)
    assert r.exit_code == 0 and len(rows) == 4 and {x["assay"] for x in rows} == {"visium", "visium_hd"}
    r = runner.invoke(cli, ["ls", "-L", lib["dir"], "--analyses"])
    assert r.exit_code == 0 and "Visium HD" in r.output and "demo_hd_day2" in r.output
    aid = by_pipeline(lib["index"])[("spaceranger", "visium_hd")]["id"]
    r = runner.invoke(cli, ["show", aid, "-L", lib["dir"]])
    assert r.exit_code == 0 and "UMI Registration Warning" in r.output and "run not in library" in r.output and "same reads" in r.output
    r = runner.invoke(cli, ["show", aid, "-L", lib["dir"], "--status"])
    assert r.exit_code == 0 and "sp_pipeline_alarms" in r.output and "sr_antisense_fraction" in r.output
    r = runner.invoke(cli, ["show", aid, "-L", lib["dir"], "--section", "sc"])
    assert r.exit_code == 0 and "spatial_units_under_tissue" in r.output and "(from report)" in r.output
    r = runner.invoke(cli, ["show", aid, "-L", lib["dir"], "--inputs"])
    assert r.exit_code == 0 and "spaceranger.metrics_summary" in r.output and "spaceranger.feature_slice" in r.output
    bundle = json.loads(runner.invoke(cli, ["show", aid, "-L", lib["dir"], "--json"]).output)
    assert bundle["bundle_id"] == aid and bundle["status"]["verdict"] == "warn"


def test_the_hub_and_the_server_carry_the_analyses(lib):
    from flowcellar.report import api as report_api
    store = LibraryStore.open(lib["dir"])
    for a in lib["index"]["analyses"]:
        path = os.path.join(lib["dir"], "hub", "data", a["id"] + ".js")
        text = open(path).read()
        assert text.startswith(f'__FCB_SIDECAR__("{a["id"]}","') and text.endswith('");\n')
        served = report_api.serve_file(store, f"/data/{a['id']}.js")
        assert served is not None and served.body.decode() == text and "Content-Security-Policy" in served.headers
    assert report_api.serve_file(store, "/data/an_spaceranger_000000000000.js") is None
    assert report_api.serve_file(store, "/data/an_../../LIBRARY.js") is None
    hub_index = open(os.path.join(lib["dir"], "hub", "index.js")).read()
    assert hub_index.startswith('__FC_INDEX__("')


def test_a_run_that_arrives_later_resolves_the_links_that_named_it(tmp_path):
    tree = fx.make_spatial_tree(str(tmp_path / "data"))
    lib_dir = str(tmp_path / "library")
    run_scan(lib_dir, tree.roots["store"], no_hub=True)          # the analyses first, and none of their runs
    rows = by_pipeline(index_of(lib_dir))
    wf = rows[("wf-single-cell", "visium_hd")]
    assert [ln["resolved"] for ln in wf["links"]] == [False, False] and [ln["fc"] for ln in wf["links"]] == ["TST00001", "TST00077"]
    assert wf["verdict"] == "pass"                               # a run that is not in the library is not a finding about the analysis
    assert validate_paths([lib_dir]) and all(r["ok"] for r in validate_paths([lib_dir]))
    res = run_scan(lib_dir, tree.late_root, no_hub=True)
    assert res.exit_code == 0
    rows = by_pipeline(index_of(lib_dir))
    key_b = tree.runs["b"].run_key
    for pipe in ("wf-single-cell", "spaceranger"):
        links = rows[(pipe, "visium_hd")]["links"]
        assert [(ln["run"], ln["resolved"]) for ln in links] == [(tree.runs["a"].run_key, False), (key_b, True)], pipe
    stored = sc_store.AnalysisStore(lib_dir).read_measured(wf["id"])["identity"]["analysis"]["links"]
    assert [ln["run_key"] for ln in stored] == [None, key_b]


def test_adding_one_run_by_hand_keeps_the_analyses_in_the_index(tmp_path):
    tree = fx.make_spatial_tree(str(tmp_path / "data"))
    lib_dir = str(tmp_path / "library")
    run_scan(lib_dir, tree.roots["store"], no_hub=True)
    r = CliRunner().invoke(cli, ["run", tree.runs["b"].run_dir, "-L", lib_dir])
    assert r.exit_code == 0, r.output
    idx = index_of(lib_dir)
    assert idx["runs"]["n"] == 1 and len(idx["analyses"]) == 2
    wf = by_pipeline(idx)[("wf-single-cell", "visium_hd")]
    assert [ln["resolved"] for ln in wf["links"]] == [False, True]     # the index joins on the run key


def test_a_library_without_analyses_has_an_empty_table_and_does_not_load_the_reader(tmp_path):
    from flowcellar.testing.fixtures import make_library_tree
    from flowcellar.testing.trees import age_tree
    root = tmp_path / "data"
    make_library_tree(root)
    age_tree(root)
    lib_dir = str(tmp_path / "library")
    import subprocess
    import sys
    code = ("import sys, gzip, json\n"
            "from flowcellar.library import scan as s\n"
            f"s.scan({lib_dir!r}, s.ScanOptions(roots=({str(root)!r},), wait=False, no_hub=True, quick=True))\n"
            f"idx = json.load(gzip.open({os.path.join(lib_dir, 'index', 'index.json.gz')!r}))\n"
            "assert idx['analyses'] == [] and idx['runs']['n'] == 3, idx['runs']['n']\n"
            "bad = [m for m in ('h5py', 'flowcellar.sc.collect', 'flowcellar.sc.store', 'flowcellar.sc.spaceranger') if m in sys.modules]\n"
            "assert not bad, bad\n")
    done = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert done.returncode == 0, done.stderr[-2000:]


def test_an_analysis_folder_that_cannot_be_read_costs_a_line_not_the_scan(tmp_path):
    tree = fx.make_spatial_tree(str(tmp_path / "data"))
    outs = tree.analyses["sr_visium"].path
    with open(os.path.join(outs, "metrics_summary.csv"), "w") as fh:
        fh.write("only one line\n")
    os.utime(os.path.join(outs, "metrics_summary.csv"), (1e9, 1e9))
    lib_dir = str(tmp_path / "library")
    res = run_scan(lib_dir, tree.roots["seq"], tree.roots["store"], no_hub=True)
    idx = index_of(lib_dir)
    assert idx["runs"]["n"] == 1 and len(idx["analyses"]) == 3 and ("spaceranger", "visium") not in by_pipeline(idx)
    assert res.exit_code in (0, 5)
