"""The library together with the real engine, the real read-data code and the real hub writer, on fixtures.

The other library tests use stand-ins so that they test the library alone. These run the whole chain the way
``flowcellar scan`` does and check what the LIBRARY promises about it. They are skipped while a part is not built.
"""
from __future__ import annotations

import gzip
import json
import os
import time

import pytest
from library_helpers import log_events, run_scan, set_root

from flowcellar.library.index import column, rows
from flowcellar.library.maintain import validate_paths
from flowcellar.library.store import LibraryStore
from flowcellar.testing import trees
from flowcellar.util.errors import FlowCellarError, NotBuiltError


@pytest.fixture(scope="module")
def real(tmp_path_factory):
    """One library scanned once from the three-root tree, with nothing replaced."""
    from flowcellar.report import api as report_api
    try:
        report_api.viewer_manifest()
    except FlowCellarError as e:
        pytest.skip(f"the viewer is not built: {e.message}")
    base = tmp_path_factory.mktemp("real_parts")
    tree = trees.make_three_roots(base / "data", reads=1500, hours=8.0)
    lib_dir = str(base / "library")
    try:
        t0 = time.monotonic()
        res = run_scan(lib_dir, *tree.root_list())
        seconds = time.monotonic() - t0
    except NotBuiltError as e:
        pytest.skip(e.message)
    return {"tree": tree, "lib": lib_dir, "first": res, "seconds": seconds}


def index_of(lib_dir: str) -> dict:
    with gzip.open(os.path.join(lib_dir, "index", "index.json.gz")) as fh:
        return json.load(fh)


def test_the_first_scan_builds_a_library_whose_bundles_are_valid(real):
    tree, lib, res = real["tree"], real["lib"], real["first"]
    assert res.exit_code == 0 and len(res.new) == 6
    idx = index_of(lib)
    assert idx["runs"]["n"] == 6 and sorted(column(idx, "detail_t1"), key=str) == ["full"] * 5 + ["none"]
    by_id = {r["id"]: r for r in rows(idx)}
    for run in tree.runs.values():
        assert by_id[run.run_key]["reads"] == run.expected["reads"], "the read count of the full pass is in the index"
        assert by_id[run.run_key]["fc"] == run.flow_cell_id and by_id[run.run_key]["headline"]
    assert by_id["TST00009_0f0e0d0c"]["verdict"] == "incomplete"
    assert by_id[tree.runs["a2"].run_key]["use"] == 2, "the second loading of a flow cell is known from the library"
    results = validate_paths([lib])
    assert [r["problems"] for r in results if not r["ok"]] == [] and len(results) == 6
    for name in ("index.html", "index.js"):
        assert os.path.getsize(os.path.join(lib, "hub", name)) > 0
    assert sorted(os.listdir(os.path.join(lib, "hub", "data"))) == sorted(k + ".js" for k in by_id)
    assert idx["library"]["counts"]["runs"] == 6 and idx["library"]["path"] == lib
    assert len(idx["projects"]) >= 3 and [r["id"] for r in idx["roots"]] == ["root_seq", "root_store", "root_arch"]


def test_the_second_scan_reads_no_input(real):
    lib = real["lib"]
    before = {sub: {n: os.stat(os.path.join(lib, sub, n)).st_mtime_ns for n in os.listdir(os.path.join(lib, sub))}
              for sub in ("bundles", "status", "cards")}
    t0 = time.monotonic()
    res = run_scan(lib)
    seconds = time.monotonic() - t0
    assert res.exit_code == 0 and res.new == [] and res.changed == []
    events = {e["event"] for e in log_events(res)}
    assert not {"identify", "quick", "pass.start", "read_pass.start", "bundle.written", "judged"} & events
    counts = log_events(res, "scan.done")[0]["counts"]
    assert counts["identify"] == counts["quick"] == counts["full_pass"] == counts["read_pass"] == counts["bundles_written"] == 0
    after = {sub: {n: os.stat(os.path.join(lib, sub, n)).st_mtime_ns for n in os.listdir(os.path.join(lib, sub))}
             for sub in ("bundles", "status", "cards")}
    assert after == before
    assert seconds < max(5.0, real["seconds"])


def test_projects_moves_and_a_lost_root(real):
    tree, lib = real["tree"], real["lib"]
    set_root(lib, "root_store", project_depth=2)
    run_scan(lib)
    idx = index_of(lib)
    alpha = next(p for p in idx["projects"] if p["id"] == "alpha_project")
    assert alpha["roots"] == ["root_seq", "root_store"] and alpha["runs"] == 3
    b2 = tree.runs["b2"]
    bundle = os.path.join(lib, "bundles", b2.run_key + ".fcb.json.gz")
    stamp = os.stat(bundle).st_mtime_ns
    new_dir = trees.move_run(b2, os.path.join(tree.roots["arch"], "beta_project", "archived"))
    res = run_scan(lib)
    assert res.moved == [b2.run_key] and res.changed == [] and os.stat(bundle).st_mtime_ns == stamp
    row = {r["id"]: r for r in rows(index_of(lib))}[b2.run_key]
    assert row["roots"] == ["root_arch"] and row["folder"] == os.path.basename(new_dir)
    beta = next(p for p in index_of(lib)["projects"] if p["id"] == "beta_project")
    assert beta["roots"] == ["root_seq", "root_arch"], "one project under two roots is one project"
    trees.take_offline(tree.roots["seq"], "empty")
    try:
        res = run_scan(lib)
        assert res.exit_code == 5 and res.offline_roots == ["root_seq"] and res.missing == []
        idx = index_of(lib)
        assert [n["code"] for n in idx["notices"]] == ["root_offline"] and idx["runs"]["n"] == 6
    finally:
        trees.bring_back(tree.roots["seq"])
    assert run_scan(lib).exit_code == 0 and index_of(lib)["notices"] == []


def test_serve_hands_out_the_real_hub(real):
    import http.client

    from flowcellar.serve import HubServer
    srv = HubServer(LibraryStore.open(real["lib"]), "127.0.0.1", 0, quiet=True)
    thread = srv.start_in_thread()
    try:
        port = srv.server_address[1]
        key = real["tree"].runs["a1"].run_key
        for path, start in (("/", b"<!doctype html"), ("/index.js", b"__FC_INDEX__("), (f"/data/{key}.js", b"__FCB_SIDECAR__(")):
            conn = http.client.HTTPConnection("127.0.0.1", port, timeout=10)
            conn.request("GET", path)
            resp = conn.getresponse()
            body = resp.read()
            conn.close()
            assert resp.status == 200 and body[:40].lower().startswith(start.lower()), path
            assert "script-src" in resp.getheader("Content-Security-Policy") and resp.getheader("Cache-Control") == "no-store"
        conn = http.client.HTTPConnection("127.0.0.1", port, timeout=10)
        conn.request("GET", "/flowcellar.toml")
        assert conn.getresponse().status == 404
        conn.close()
    finally:
        srv.shutdown()
        srv.server_close()
        thread.join(5)
