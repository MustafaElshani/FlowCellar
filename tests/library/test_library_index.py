"""The index: the one payload Home, Runs and Projects are drawn from."""
from __future__ import annotations

import gzip
import json
import os

from library_helpers import run_scan, set_root

from flowcellar.library import api
from flowcellar.library.index import COLUMN_ORDER, build_index, column, rows
from flowcellar.library.store import LibraryStore


def index_of(lib_dir: str) -> dict:
    with gzip.open(os.path.join(lib_dir, "index", "index.json.gz")) as fh:
        return json.load(fh)


def test_the_tables_home_needs(scanned, tree, lib_dir):
    idx = index_of(lib_dir)
    assert list(idx) == ["format", "schema", "written", "meta", "library", "runs", "samples", "projects", "roots", "activity",
                         "baselines", "facets", "notices", "progress", "analyses"]
    assert idx["analyses"] == []                          # the table of the single-cell and spatial analyses: none here
    assert idx["format"] == "flowcellar.index" and idx["schema"] == "0.1.1"
    runs = tree.runs
    total_reads = sum(r.expected["reads"] for r in runs.values())

    # ---- library: name, counts, totals, when it was written, and its own path for the commands shown on Home
    lib = idx["library"]
    assert lib["name"] == "library" and lib["path"] == lib_dir and lib["id"].startswith("lib_")
    assert lib["created"] and lib["scanned"] and lib["scanned"] >= lib["created"]
    assert lib["counts"] == {"runs": 6, "flow_cells": 5, "projects": 2, "samples": 0,
                             "by_verdict": {"pass": 5, "warn": 0, "fail": 0, "incomplete": 1, "running": 0, "unknown": 0},
                             "by_detail": {"full": 5, "quick": 0, "none": 1}, "with_read_data": 0}
    assert lib["reads"] == total_reads and lib["estimated_bases"] == total_reads * 900
    assert lib["first_run"] == "2026-02-01T07:00:00Z" and lib["last_run"] == "2026-04-10T11:00:00Z"
    assert idx["meta"]["library"] == "library" and idx["meta"]["runs"] == 6 and idx["meta"]["generator"]["name"] == "flowcellar"

    # ---- projects: newest first; alpha_project lies under two roots and is ONE project
    assert [p["id"] for p in idx["projects"]] == ["alpha_project", "beta_project", "_none"]
    alpha, beta, none = idx["projects"]
    assert alpha == {"id": "alpha_project", "name": "alpha_project", "runs": 3, "flow_cells": 2, "samples": 0,
                     "first": "2026-03-02", "last": "2026-04-10",
                     "verdicts": {"pass": 3, "warn": 0, "fail": 0, "incomplete": 0, "running": 0, "unknown": 0},
                     "estimated_bases": sum(runs[k].expected["reads"] for k in ("a1", "a2", "a3")) * 900,
                     "reads": sum(runs[k].expected["reads"] for k in ("a1", "a2", "a3")), "roots": ["root_seq", "root_store"]}
    assert beta["runs"] == 2 and beta["roots"] == ["root_seq"] and beta["first"] == "2026-02-01" and beta["last"] == "2026-03-05"
    assert none["name"] is None and none["runs"] == 1 and none["roots"] == ["root_arch"] and none["verdicts"]["incomplete"] == 1

    # ---- roots: one row per [[roots]] entry, in the order of the configuration
    assert [(r["id"], r["path"], r["runs"], r["state"]) for r in idx["roots"]] == [
        ("root_seq", tree.roots["seq"], 4, "reachable"), ("root_store", tree.roots["store"], 1, "reachable"),
        ("root_arch", tree.roots["arch"], 1, "reachable")]
    assert [r["detail"] for r in idx["roots"]] == [{"full": 4, "quick": 0, "raw_only": 0}, {"full": 1, "quick": 0, "raw_only": 0},
                                                   {"full": 0, "quick": 0, "raw_only": 1}]
    assert [r["project_depth"] for r in idx["roots"]] == [1, 2, 1] and all(r["last_scanned"] for r in idx["roots"])
    assert all(r["offline_since"] is None and r["live"] is False for r in idx["roots"])

    # ---- activity: every month from the first run to the last
    assert idx["activity"]["months"] == ["2026-02", "2026-03", "2026-04"] and idx["activity"]["runs"] == [1, 3, 1]
    assert sum(idx["activity"]["reads"]) == total_reads and sum(idx["activity"]["estimated_bases"]) == lib["estimated_bases"]

    # ---- nothing to say, nothing in progress
    assert idx["notices"] == []
    assert idx["progress"]["state"] == "idle" and idx["progress"]["tier1"]["runs_waiting"] == 0 and idx["progress"]["message"] is None
    assert idx["facets"]["root"] == ["root_arch", "root_seq", "root_store"]
    assert idx["facets"]["project"] == ["alpha_project", "beta_project"] and idx["facets"]["verdict"] == ["incomplete", "pass"]
    assert idx["baselines"] == {"groups": []}


def test_runs_are_columns_newest_first_with_stored_sentences(scanned, tree, lib_dir):
    idx = index_of(lib_dir)
    cols = idx["runs"]["columns"]
    assert idx["runs"]["n"] == 6 and all(len(column(idx, c)) == 6 for c in cols)
    started = column(idx, "started")
    assert started[:5] == sorted(started[:5], reverse=True) and started[5] is None, "newest first; no start time last"
    assert column(idx, "id")[0] == tree.runs["a3"].run_key
    # low-cardinality text is dictionary-encoded, identifiers and sentences are not; null is an entry of the dictionary
    assert isinstance(cols["kit"], dict) and set(cols["kit"]) == {"dict", "codes"} and None in cols["kit"]["dict"]
    assert all(isinstance(cols[c], list) for c in ("id", "fc", "started", "headline", "roots", "project_ids", "reasons"))
    assert -1 not in cols["kit"]["codes"]
    order = [c for c in COLUMN_ORDER if c in cols]
    assert list(cols)[:len(order)] == order, "columns come in the order of the plan"
    by_id = {r["id"]: r for r in rows(idx)}
    a3 = by_id[tree.runs["a3"].run_key]
    assert a3["roots"] == ["root_store"] and a3["project"] == "alpha_project" and a3["project_ids"] == ["alpha_project"]
    assert a3["folder"] == os.path.basename(tree.runs["a3"].run_dir) and a3["source_state"] == "ok"
    assert a3["bundle"] == a3["id"] + ".fcb.json.gz" and a3["bundle_gzip_bytes"] == os.path.getsize(
        os.path.join(lib_dir, "bundles", a3["bundle"])) and a3["bundle_bytes"] > a3["bundle_gzip_bytes"]
    status = json.load(open(os.path.join(lib_dir, "status", a3["id"] + ".status.json")))
    assert a3["headline"] == status["status"]["headline"], "word for word as in the status file"
    raw = by_id["TST00009_0f0e0d0c"]
    assert raw["project"] is None and raw["project_ids"] == [] and raw["verdict"] == "incomplete" and raw["cause"] == "raw_only"


def test_a_project_under_two_roots_follows_the_configuration(parts, tree, lib_dir):
    run_scan(lib_dir, *tree.root_list())                               # every root: project = first folder
    idx = index_of(lib_dir)
    assert [p["id"] for p in idx["projects"]] == ["group_one", "beta_project", "alpha_project", "_none"]
    before = os.stat(os.path.join(lib_dir, "bundles", tree.runs["a3"].run_key + ".fcb.json.gz")).st_mtime_ns
    judged = parts.engine.calls["judge"]
    set_root(lib_dir, "root_store", project_depth=2)                   # the person corrects one line
    run_scan(lib_dir)
    idx = index_of(lib_dir)
    assert [(p["id"], p["roots"]) for p in idx["projects"]] == [("alpha_project", ["root_seq", "root_store"]),
                                                                ("beta_project", ["root_seq"]), ("_none", ["root_arch"])]
    assert os.stat(os.path.join(lib_dir, "bundles", tree.runs["a3"].run_key + ".fcb.json.gz")).st_mtime_ns == before
    assert parts.engine.calls["judge"] == judged, "a project name is library knowledge: nothing is judged again"
    set_root(lib_dir, "root_store", project_depth=0)                   # 0: no project from this root
    run_scan(lib_dir)
    assert {r["id"]: r for r in rows(index_of(lib_dir))}[tree.runs["a3"].run_key]["project"] is None


def test_project_pattern(parts, tree, lib_dir):
    run_scan(lib_dir, tree.roots["store"])
    path = os.path.join(lib_dir, "flowcellar.toml")
    text = open(path).read().replace("project_depth = 1", 'project_pattern = "group_*/{project}_project"')
    open(path, "w").write(text)
    run_scan(lib_dir)
    assert [(p["id"], p["name"]) for p in index_of(lib_dir)["projects"]] == [("alpha", "alpha")]


def test_the_index_of_an_empty_library_and_of_one_run_without_a_library(parts, tree, tmp_path):
    empty = build_index([])
    assert empty["runs"]["n"] == 0 and empty["projects"] == [] and empty["roots"] == []
    assert empty["activity"] == {"months": [], "runs": [], "estimated_bases": [], "reads": []}
    assert empty["library"]["counts"]["runs"] == 0 and empty["library"]["path"] is None
    assert empty["progress"]["state"] == "idle" and all(v == [] for v in empty["runs"]["columns"].values())
    # a one-run report builds its index the same way: build_index([card_for(bundle, status_file)])
    from flowcellar.engine import api as engine_api
    ident = engine_api.identify(tree.runs["a1"].run_dir)
    bundle = engine_api.quick_bundle(ident)
    sf = engine_api.judge(bundle, [], in_library=False)
    card = api.card_for(bundle, sf)
    assert card["run_key"] == ident.run_key and card["columns"]["roots"] == [] and card["columns"]["project"] is None
    assert card["columns"]["folder"] == os.path.basename(tree.runs["a1"].run_dir)
    assert card["columns"]["bundle_gzip_bytes"] > 0 and "facts" in card and card["samples"] == []
    one = api.build_index([card])
    assert one["runs"]["n"] == 1 and column(one, "id") == [ident.run_key] and one["library"]["counts"]["flow_cells"] == 1
    assert one["projects"] == [{"id": "_none", "name": None, "runs": 1, "flow_cells": 1, "samples": 0,
                                "first": "2026-03-02", "last": "2026-03-02",
                                "verdicts": {"pass": 1, "warn": 0, "fail": 0, "incomplete": 0, "running": 0, "unknown": 0},
                                "estimated_bases": one["library"]["estimated_bases"], "reads": 0, "roots": []}]


def test_the_library_protocol_other_parts_read_through(scanned, tree, lib_dir):
    lib = api.open_library(lib_dir)
    key = tree.runs["a1"].run_key
    assert lib.path == lib_dir and lib.id.startswith("lib_") and lib.name == "library" and lib.config["scan"]["threads"] == 4
    assert len(lib.run_keys()) == 6 and key in lib.run_keys()
    measured = lib.read_measured(key)
    assert "status" not in measured and "relations" not in measured and measured["content_sha256"]
    sf = lib.read_status_file(key)
    assert sf["format"] == "flowcellar.status" and sf["bundle_id"] == key and sf["content_sha256"] == measured["content_sha256"]
    assert [p["root"] for p in sf["relations"]["paths"]] == ["root_seq"]
    composed = lib.composed(key)
    assert composed["status"] == sf["status"] and composed["relations"]["projects"][0]["id"] == "alpha_project"
    assert composed["user"]["tags"] == [] and list(composed)[:2] == ["format", "schema_version"]
    assert json.loads(gzip.decompress(lib.composed_gzip(key))) == json.loads(json.dumps(composed))
    assert lib.card(key)["columns"]["id"] == key and lib.card("TST00099_00000000") is None
    assert lib.index()["runs"]["n"] == 6 and gzip.decompress(lib.index_gzip())[:1] == b"{"
    assert lib.hub_dir() == os.path.join(lib_dir, "hub")
    assert isinstance(lib, LibraryStore)


def test_open_library_creates_only_when_asked(tmp_path):
    import pytest

    from flowcellar.util.errors import InputError, LibraryError
    with pytest.raises(LibraryError) as e:
        api.open_library(str(tmp_path / "nope"))
    assert e.value.exit_code == 6
    lib = api.open_library(str(tmp_path / "new"), create=True, name="fresh", roots=[str(tmp_path / "data")])
    assert lib.name == "fresh" and [r["path"] for r in lib.config["roots"]] == [str(tmp_path / "data")]
    assert api.open_library(str(tmp_path / "new")).id == lib.id
    marker = tmp_path / "new" / "LIBRARY"
    marker.write_text(json.dumps({**json.loads(marker.read_text()), "layout_version": 99}))
    with pytest.raises(LibraryError) as e:
        api.open_library(str(tmp_path / "new"))
    assert "newer FlowCellar" in e.value.message
    with pytest.raises(InputError):
        lib.read_measured("../../etc/passwd")
