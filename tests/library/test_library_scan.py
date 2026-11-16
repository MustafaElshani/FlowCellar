"""The scan: what the lab asked about. "Do people have a library built, and when they import a flow cell or
project the library is built? How are we going to handle multiple projects located in multiple areas?"

Everything here runs on the fixture tree of ``flowcellar.testing.trees`` (three roots, invented names) with the
stand-ins of ``fakes.py`` for the engine, read data and the hub writer.
"""
from __future__ import annotations

import gzip
import json
import os
import shutil
import time

import pytest
from library_helpers import log_events, run_scan, set_root

from flowcellar.library import scan as scan_mod
from flowcellar.library.index import column, rows
from flowcellar.library.lease import LeaseBusy
from flowcellar.library.store import LibraryStore
from flowcellar.testing import trees
from flowcellar.util.errors import LibraryError, UsageError


def index_of(lib_dir: str) -> dict:
    with gzip.open(os.path.join(lib_dir, "index", "index.json.gz")) as fh:
        return json.load(fh)


def files_under(folder: str) -> dict[str, int]:
    out = {}
    for d, _dirs, names in os.walk(folder):
        for n in names:
            p = os.path.join(d, n)
            out[os.path.relpath(p, folder)] = os.stat(p).st_mtime_ns
    return out


def snapshot(lib_dir: str) -> dict[str, int]:
    """Modification times of everything a scan must not touch when nothing changed."""
    out = {}
    for sub in ("bundles", "status", "cards", "state", "partials"):
        for name, m in files_under(os.path.join(lib_dir, sub)).items():
            if name not in ("roots.json", "notices.json"):
                out[sub + "/" + name] = m
    return out


# ------------------------------------------------------------------------------------------------ creation and roots
def test_scan_creates_the_library_and_remembers_the_roots(parts, tree, lib_dir):
    assert not os.path.exists(lib_dir)
    res = run_scan(lib_dir, *tree.root_list(), name="demo library")
    assert res.exit_code == 0 and res.created
    assert sorted(res.new) == sorted([r.run_key for r in tree.runs.values()] + list(tree.raw_only))
    store = LibraryStore.open(lib_dir)
    assert store.name == "demo library"
    assert [r["id"] for r in store.config["roots"]] == ["root_seq", "root_store", "root_arch"]
    assert [r["path"] for r in store.config["roots"]] == tree.root_list()
    marker = json.load(open(os.path.join(lib_dir, "LIBRARY")))
    assert marker["format"] == "flowcellar.library" and marker["layout_version"] == 1 and marker["library_id"].startswith("lib_")
    for key in res.new:
        for sub, suffix in (("bundles", ".fcb.json.gz"), ("status", ".status.json"), ("cards", ".card.json"), ("state", ".state.json")):
            assert os.path.isfile(os.path.join(lib_dir, sub, key + suffix)), (sub, key)
    assert os.path.isfile(os.path.join(lib_dir, "hub", "index.html"))
    assert not os.path.exists(os.path.join(lib_dir, "claims", "scan.lease")), "the lease is released"
    # the run that lies under a folder pruned by name was not found, the look-alike is not a run
    assert all("TST00008" not in k for k in res.new)
    assert len(store.run_keys()) == 6


def test_nothing_is_written_outside_the_library(parts, tree, lib_dir):
    before = files_under(tree.base)
    run_scan(lib_dir, *tree.root_list())
    run_scan(lib_dir)
    assert files_under(tree.base) == before


def test_a_scan_without_roots_uses_the_remembered_ones_and_a_new_root_is_added(parts, tree, lib_dir, tmp_path):
    run_scan(lib_dir, tree.roots["seq"])
    assert len(LibraryStore.open(lib_dir).run_keys()) == 4
    text_before = open(os.path.join(lib_dir, "flowcellar.toml")).read()
    res = run_scan(lib_dir)                                           # no roots: the remembered one
    assert res.new == [] and res.exit_code == 0
    assert open(os.path.join(lib_dir, "flowcellar.toml")).read() == text_before
    res = run_scan(lib_dir, tree.roots["store"])                      # a new root later: added, all are scanned
    text = open(os.path.join(lib_dir, "flowcellar.toml")).read()
    assert text.startswith(text_before) and text.count("[[roots]]") == text_before.count("[[roots]]") + 1
    assert "added by flowcellar scan" in text
    assert res.new == [tree.runs["a3"].run_key]
    assert [r["id"] for r in LibraryStore.open(lib_dir).config["roots"]] == ["root_seq", "root_store"]
    res = run_scan(lib_dir, tree.roots["store"], tree.roots["seq"])   # known roots given again change nothing
    assert open(os.path.join(lib_dir, "flowcellar.toml")).read() == text and res.new == []
    res = run_scan(lib_dir, os.path.join(tree.roots["seq"], "alpha_project"))     # inside a known root: not added
    assert open(os.path.join(lib_dir, "flowcellar.toml")).read() == text


def test_two_roots_with_the_same_last_name_get_distinct_ids(parts, tmp_path, lib_dir):
    a, b = tmp_path / "site_a" / "projects", tmp_path / "site_b" / "projects"
    a.mkdir(parents=True)
    b.mkdir(parents=True)
    run_scan(lib_dir, str(a), str(b))
    assert [r["id"] for r in LibraryStore.open(lib_dir).config["roots"]] == ["site_a_projects", "site_b_projects"]


def test_scan_refuses_a_folder_that_is_not_a_library_and_a_library_without_roots(parts, tree, tmp_path):
    other = tmp_path / "documents"
    other.mkdir()
    (other / "thesis.txt").write_text("not a library")
    with pytest.raises(LibraryError) as e:
        run_scan(str(other), tree.roots["seq"])
    assert e.value.exit_code == 6 and (other / "thesis.txt").read_text() == "not a library"
    assert sorted(os.listdir(other)) == ["thesis.txt"]
    with pytest.raises(UsageError):
        run_scan(str(tmp_path / "new_library"))                       # no library there and no root given


# ------------------------------------------------------------------------------------------------ order of work
def test_quick_verdicts_for_every_run_come_before_any_pass_and_the_hub_is_written_in_between(parts, tree, lib_dir):
    order = []
    eng = parts.engine
    real_quick, real_pass = eng.quick_bundle, eng.full_pass
    hub_real = parts.report.write_hub

    def quick(run, ctx=None):
        order.append("quick")
        return real_quick(run, ctx)

    def full(run, partials_dir, opts=None):
        order.append("pass:" + run.run_key)
        return real_pass(run, partials_dir, opts)

    def hub(library, **kw):
        order.append("hub")
        return hub_real(library, **kw)
    from flowcellar.engine import api as engine_api
    from flowcellar.report import api as report_api
    eng.quick_bundle, eng.full_pass = quick, full
    engine_api.quick_bundle, engine_api.full_pass, report_api.write_hub = quick, full, hub
    try:
        res = run_scan(lib_dir, *tree.root_list())
    finally:
        engine_api.quick_bundle, engine_api.full_pass, report_api.write_hub = real_quick, real_pass, hub_real
    first_pass = next(i for i, x in enumerate(order) if x.startswith("pass:"))
    assert order[:first_pass].count("quick") == 6, "every new run has its quick bundle before the first pass"
    assert "hub" in order[:first_pass], "the hub is written before the counting starts"
    # passes run newest run first
    started = {r.run_key: r.spec.start for r in tree.runs.values()}
    passed = [x.split(":")[1] for x in order if x.startswith("pass:")]
    assert passed == sorted(passed, key=lambda k: started[k], reverse=True)
    # after the scan every run with reports has full detail
    idx = index_of(lib_dir)
    assert sorted(column(idx, "detail_t1")) == ["full"] * 5 + ["none"]
    assert res.waiting_tier1 == 0


def test_runs_are_judged_in_order_of_run_start(parts, tree, lib_dir):
    run_scan(lib_dir, *tree.root_list())
    started = {r.run_key: r.spec.start for r in tree.runs.values()}
    first_round = parts.engine.judged[:6]                             # the quick round: oldest first
    with_start = [k for k, _ in first_round if k in started]
    assert with_start == sorted(with_start, key=lambda k: started[k])
    # each run saw the earlier ones as history when it was judged for the last time
    last = {}
    for key, history in parts.engine.judged:
        last[key] = history
    newest = max(started, key=started.get)
    assert set(last[newest]) >= set(started) - {newest}
    # the flow cell that was loaded twice: the later run knows it is the second use
    idx = index_of(lib_dir)
    by_id = {r["id"]: r for r in rows(idx)}
    assert by_id[tree.runs["a1"].run_key]["use"] == 1 and by_id[tree.runs["a2"].run_key]["use"] == 2


def test_the_terminal_lines_say_what_was_found_and_what_is_still_counted(parts, tree, lib_dir, capsys):
    from flowcellar.util import log
    log.setup(0)
    run_scan(lib_dir, *tree.root_list())
    err = capsys.readouterr().err
    lines = err.strip().split("\n")
    assert lines[0].startswith("created library ") and "(3 roots; edit flowcellar.toml to add file patterns)" in lines[0]
    assert lines[1].startswith("found 6 runs (5 with MinKNOW reports, 1 with raw data only), ") and "of read tables to count" in lines[1]
    run_lines = lines[2:8]
    assert run_lines[0].startswith(tree.runs["a3"].run_key) and "Pass so far" in run_lines[0], "newest run first"
    assert "wait for the read count" in run_lines[0]
    assert any(x.startswith("TST00009_0f0e0d0c") and "Incomplete" in x for x in run_lines)
    assert lines[8].startswith("hub ready: ") and lines[8].split()[2].endswith("hub/index.html")
    assert lines[-1].startswith("done: 6 runs · 5 Pass, 0 Warning, 0 Fail, 1 Incomplete · read data not linked yet (no file patterns)")
    run_scan(lib_dir)
    assert capsys.readouterr().err.strip() == \
        "done: nothing new · 6 runs · 5 Pass, 0 Warning, 0 Fail, 1 Incomplete · read data not linked yet (no file patterns)"


# ------------------------------------------------------------------------------------------------ incremental
def test_a_second_scan_reads_no_input_and_writes_no_run_file(parts, tree, lib_dir):
    run_scan(lib_dir, *tree.root_list())
    eng = parts.engine
    calls, reads = dict(eng.calls), list(eng.reads)
    assert len(reads) == 10, "five runs: one pass over the read table and one over the pore-scan table each"
    before = snapshot(lib_dir)
    opened = []
    import builtins
    real_open = builtins.open

    def counting_open(file, *a, **k):
        if isinstance(file, (str, bytes, os.PathLike)) and os.fspath(file).startswith(tree.base):
            opened.append(os.fspath(file))
        return real_open(file, *a, **k)
    builtins.open = counting_open
    try:
        t0 = time.monotonic()
        res = run_scan(lib_dir)
        seconds = time.monotonic() - t0
    finally:
        builtins.open = real_open
    assert res.exit_code == 0 and res.new == [] and res.changed == [] and res.moved == []
    assert opened == [], "no file under any root was opened"
    for name in ("identify", "quick_bundle", "full_pass", "assemble", "judge"):
        assert eng.calls[name] == calls[name], name
    assert eng.reads == reads
    assert snapshot(lib_dir) == before, "no bundle, status file, card, state file or partial was rewritten"
    assert seconds < 5
    events = [e["event"] for e in log_events(res)]
    assert not {"identify", "quick", "pass.start", "read_pass.start", "bundle.written", "judged", "link"} & set(events)
    done = log_events(res, "scan.done")[0]
    assert done["counts"]["identify"] == 0 and done["counts"]["bundles_written"] == 0 and done["counts"]["status_written"] == 0


def test_a_later_run_changes_no_stored_verdict(parts, tree, lib_dir):
    run_scan(lib_dir, tree.roots["seq"])
    before = snapshot(lib_dir)
    res = run_scan(lib_dir, tree.roots["store"])                      # brings a3, the newest run
    assert res.new == [tree.runs["a3"].run_key]
    after = snapshot(lib_dir)
    for name, m in before.items():
        if name.startswith(("bundles/", "status/", "cards/")):
            assert after[name] == m, f"{name} was rewritten because a later run arrived"
    assert not [n for n in index_of(lib_dir)["notices"] if n["code"] == "judged_before_older_runs"]


def test_an_older_run_added_later_is_reported_and_verdicts_are_kept(parts, tree, lib_dir, tmp_path):
    # hold the oldest run back, scan, then let it appear
    held = str(tmp_path / "held_back")
    b2 = tree.runs["b2"]
    shutil.move(b2.run_dir, held)
    run_scan(lib_dir, tree.roots["seq"])
    before = snapshot(lib_dir)
    shutil.move(held, b2.run_dir)
    res = run_scan(lib_dir)
    assert res.new == [b2.run_key]
    after = snapshot(lib_dir)
    assert all(after[n] == m for n, m in before.items() if n.startswith("status/")), "verdicts are frozen"
    notice = [n for n in index_of(lib_dir)["notices"] if n["code"] == "judged_before_older_runs"]
    assert len(notice) == 1 and "3 runs were judged before older runs were added" in notice[0]["message"]
    assert "recompute --status-only --rebaseline" in notice[0]["message"]


def test_a_changed_input_is_counted_again_and_only_that_run(parts, tree, lib_dir):
    run_scan(lib_dir, *tree.root_list())
    eng = parts.engine
    a1 = tree.runs["a1"]
    summary = a1.files["minknow.sequencing_summary"]
    with open(summary, "a") as fh:
        fh.write("\t".join(["-"] * 16) + "\n")                        # one more read
    trees.age_tree(a1.run_dir, 1800)
    reads_before = len(eng.reads)
    before = snapshot(lib_dir)
    res = run_scan(lib_dir)
    assert eng.reads[reads_before:] == [summary], "only the changed file was read again"
    assert res.changed == [a1.run_key]
    changed = {n for n, m in snapshot(lib_dir).items() if before.get(n) != m and n.startswith("bundles/")}
    assert changed == {f"bundles/{a1.run_key}.fcb.json.gz"}
    by_id = {r["id"]: r for r in rows(index_of(lib_dir))}
    assert by_id[a1.run_key]["reads"] == a1.expected["reads"] + 1


def test_fresh_files_are_left_for_the_next_scan(parts, tmp_path, lib_dir):
    t = trees.make_three_roots(tmp_path / "fresh", aged=False)        # written a moment ago: still being copied?
    res = run_scan(lib_dir, t.roots["seq"])
    assert res.new == [] and res.exit_code == 0
    assert len(log_events(res, "run.settling")) == 4
    trees.age_tree(t.base)
    assert len(run_scan(lib_dir).new) == 4


# ------------------------------------------------------------------------------------------------ moved, mirrored, gone
def test_a_moved_run_keeps_its_bundle_and_nothing_is_computed_again(parts, tree, lib_dir):
    run_scan(lib_dir, *tree.root_list())
    eng = parts.engine
    a1 = tree.runs["a1"]
    calls = dict(eng.calls)
    bundles = files_under(os.path.join(lib_dir, "bundles"))
    new_dir = trees.move_run(a1, os.path.join(tree.roots["store"], "group_one", "alpha_project", "moved_here"))
    res = run_scan(lib_dir)
    assert res.moved == [a1.run_key] and res.new == [] and res.changed == [] and res.exit_code == 0
    assert eng.calls["identify"] == calls["identify"], "the files have the same names, sizes and times: not even opened"
    assert eng.calls["full_pass"] == calls["full_pass"] and eng.calls["quick_bundle"] == calls["quick_bundle"]
    assert files_under(os.path.join(lib_dir, "bundles")) == bundles
    state = json.load(open(os.path.join(lib_dir, "state", a1.run_key + ".state.json")))
    assert [p["path"] for p in state["paths"]] == [new_dir] and state["paths"][0]["root"] == "root_store"
    row = {r["id"]: r for r in rows(index_of(lib_dir))}[a1.run_key]
    assert row["roots"] == ["root_store"] and row["source_state"] == "ok"
    status = json.load(open(os.path.join(lib_dir, "status", a1.run_key + ".status.json")))
    assert [p["root"] for p in status["relations"]["paths"]] == ["root_store"]
    assert len(LibraryStore.open(lib_dir).run_keys()) == 6, "no second run appeared"


def test_a_moved_run_whose_file_times_changed_is_recognised_by_its_content(parts, tree, lib_dir):
    run_scan(lib_dir, *tree.root_list())
    eng = parts.engine
    a1 = tree.runs["a1"]
    reads = len(eng.reads)
    new_dir = trees.move_run(a1, os.path.join(tree.roots["arch"], "copied"))
    trees.age_tree(new_dir, 7200)                                     # a copy that did not keep the times
    res = run_scan(lib_dir)
    assert res.moved == [a1.run_key] and res.changed == []
    assert len(eng.reads) == reads, "the passes were found in the cache by the fingerprint of the files"
    assert len(LibraryStore.open(lib_dir).run_keys()) == 6


def test_the_same_run_under_two_roots_is_one_run_with_two_paths(parts, tree, lib_dir):
    b1 = tree.runs["b1"]
    copy = trees.mirror_run(b1, os.path.join(tree.roots["arch"], "mirror_of_beta", "pool_b"))
    res = run_scan(lib_dir, *tree.root_list())
    assert len(res.new) == 6 and len(LibraryStore.open(lib_dir).run_keys()) == 6
    state = json.load(open(os.path.join(lib_dir, "state", b1.run_key + ".state.json")))
    assert [p["path"] for p in state["paths"]] == [b1.run_dir, copy]
    idx = index_of(lib_dir)
    row = {r["id"]: r for r in rows(idx)}[b1.run_key]
    assert row["roots"] == ["root_seq", "root_arch"]
    assert parts.engine.reads.count(b1.files["minknow.sequencing_summary"]) == 1, "read once, from the first root"
    assert idx["library"]["counts"]["runs"] == 6 and idx["library"]["reads"] == sum(r.expected["reads"] for r in tree.runs.values())
    # the mirror appearing later is "also found at", and nothing is computed
    lib2 = lib_dir + "_2"
    shutil.rmtree(copy)
    run_scan(lib2, *tree.root_list())
    calls = dict(parts.engine.calls)
    trees.mirror_run(b1, os.path.join(tree.roots["arch"], "mirror_of_beta", "pool_b"))
    res = run_scan(lib2)
    assert res.moved == [b1.run_key] and res.new == [] and res.changed == []
    assert parts.engine.calls["full_pass"] == calls["full_pass"]
    assert {r["id"]: r for r in rows(index_of(lib2))}[b1.run_key]["roots"] == ["root_seq", "root_arch"]


def test_a_run_whose_folder_disappeared_keeps_its_bundle(parts, tree, lib_dir):
    run_scan(lib_dir, *tree.root_list())
    b1 = tree.runs["b1"]
    bundle = os.path.join(lib_dir, "bundles", b1.run_key + ".fcb.json.gz")
    before = os.stat(bundle).st_mtime_ns
    shutil.rmtree(b1.run_dir)
    res = run_scan(lib_dir)
    assert res.missing == [b1.run_key] and res.exit_code == 0
    assert os.stat(bundle).st_mtime_ns == before
    idx = index_of(lib_dir)
    row = {r["id"]: r for r in rows(idx)}[b1.run_key]
    assert row["source_state"] == "missing" and row["verdict"] == "pass" and row["reads"] == b1.expected["reads"]
    assert idx["library"]["counts"]["runs"] == 6
    assert run_scan(lib_dir).missing == [], "said once, not at every scan"


# ------------------------------------------------------------------------------------------------ roots that do not answer
@pytest.mark.parametrize("how", ["empty", "missing"])
def test_an_unreachable_root_is_skipped_and_nothing_under_it_changes(parts, tree, lib_dir, how):
    run_scan(lib_dir, *tree.root_list())
    before = snapshot(lib_dir)
    calls = dict(parts.engine.calls)
    trees.take_offline(tree.roots["seq"], how)                        # four of the six runs live there
    res = run_scan(lib_dir)
    assert res.exit_code == 5 and res.offline_roots == ["root_seq"]
    assert res.missing == [] and res.new == [] and res.changed == []
    assert snapshot(lib_dir) == before, "no run under the lost root changed in any way"
    assert parts.engine.calls["identify"] == calls["identify"]
    idx = index_of(lib_dir)
    root = {r["id"]: r for r in idx["roots"]}["root_seq"]
    assert root["state"] == "offline" and root["offline_since"] and root["runs"] == 4
    notices = [n for n in idx["notices"] if n["code"] == "root_offline"]
    assert len(notices) == 1 and notices[0]["root"] == "root_seq" and "Nothing known under it was changed" in notices[0]["message"]
    by_id = {r["id"]: r for r in rows(idx)}
    assert by_id[tree.runs["a1"].run_key]["source_state"] == "offline"
    assert by_id[tree.runs["a3"].run_key]["source_state"] == "ok"
    assert idx["library"]["counts"]["runs"] == 6
    since = root["offline_since"]
    assert {r["id"]: r for r in index_of(lib_dir)["roots"]}["root_seq"]["offline_since"] == since
    # it comes back: everything is as before, the notice is gone
    trees.bring_back(tree.roots["seq"])
    res = run_scan(lib_dir)
    assert res.exit_code == 0 and res.missing == [] and res.new == []
    idx = index_of(lib_dir)
    assert {r["id"]: r for r in idx["roots"]}["root_seq"]["state"] == "reachable"
    assert not [n for n in idx["notices"] if n["code"] == "root_offline"]
    assert snapshot(lib_dir) == before


def test_a_root_that_stalls_is_skipped(parts, tree, lib_dir, monkeypatch):
    run_scan(lib_dir, *tree.root_list())
    before = snapshot(lib_dir)
    from flowcellar.discovery import guard
    real = guard.check_roots

    def slow(roots, **kw):
        kw.update(target=trees.stalling_child, stall_seconds=0.5)
        return real(roots, **kw)
    monkeypatch.setattr("flowcellar.discovery.check_roots", slow)
    t0 = time.monotonic()
    res = run_scan(lib_dir)
    assert time.monotonic() - t0 < 15
    assert res.exit_code == 5 and sorted(res.offline_roots) == ["root_arch", "root_seq", "root_store"]
    assert snapshot(lib_dir) == before
    states = {r["id"]: r["state"] for r in index_of(lib_dir)["roots"]}
    assert set(states.values()) == {"skipped"}


def test_a_subtree_over_the_listing_budget_is_named_and_runs_below_it_are_not_called_missing(parts, tree, lib_dir):
    run_scan(lib_dir, *tree.root_list())
    trees.make_deep_tree(os.path.join(tree.roots["seq"], "alpha_project"), folders=60)
    trees.age_tree(tree.roots["seq"])
    path = os.path.join(lib_dir, "flowcellar.toml")
    open(path, "a").write("\n[discovery]\nsubtree_listing_budget = 8\n")
    res = run_scan(lib_dir)
    notices = [n for n in index_of(lib_dir)["notices"] if n["code"] == "subtree_over_budget"]
    assert len(notices) == 1 and os.path.join(tree.roots["seq"], "alpha_project") in notices[0]["message"]
    assert "add it to prune" in notices[0]["message"]
    assert res.missing == [], "runs under a subtree that was given up are not reported as gone"


# ------------------------------------------------------------------------------------------------ lease
def test_one_scan_at_a_time_and_a_stale_lease_is_taken_over(parts, tree, lib_dir):
    run_scan(lib_dir, tree.roots["seq"])
    lease = os.path.join(lib_dir, "claims", "scan.lease")
    os.makedirs(os.path.dirname(lease), exist_ok=True)
    import socket
    json.dump({"host": socket.gethostname(), "pid": os.getpid(), "started": "2026-01-01T00:00:00Z"}, open(lease, "w"))
    with pytest.raises(LeaseBusy) as e:                               # held by a live process on this host
        run_scan(lib_dir, wait=False)
    assert e.value.exit_code == 5 and str(os.getpid()) in e.value.message
    json.dump({"host": socket.gethostname(), "pid": 2 ** 22 + 12345, "started": "2026-01-01T00:00:00Z"}, open(lease, "w"))
    res = run_scan(lib_dir, wait=False)                               # that process is gone: stale
    assert res.exit_code == 0
    assert [n["code"] for n in index_of(lib_dir)["notices"]] == ["stale_lease_taken"]
    assert not os.path.exists(lease)
    # held by another host: stale only when its time stamp stands still
    json.dump({"host": "another-host", "pid": 1, "started": "2026-01-01T00:00:00Z"}, open(lease, "w"))
    with pytest.raises(LeaseBusy):
        run_scan(lib_dir, wait=False)
    res = run_scan(lib_dir, wait=True, lease_poll_seconds=0.05, lease_stale_seconds=0.3)
    assert res.exit_code == 0 and [n["code"] for n in index_of(lib_dir)["notices"]] == ["stale_lease_taken"]
    assert run_scan(lib_dir).exit_code == 0 and index_of(lib_dir)["notices"] == []


# ------------------------------------------------------------------------------------------------ time box
def test_a_time_box_stops_passes_and_a_later_scan_finishes_them(parts, tree, lib_dir, monkeypatch):
    clock = [1000.0]
    monkeypatch.setattr(scan_mod, "_monotonic", lambda: clock[0])
    seen = set()

    def tick(run, role):                                              # every run's counting takes 25 minutes
        if run.run_key not in seen:
            seen.add(run.run_key)
            clock[0] += 25 * 60
    parts.engine.on_pass = tick
    res = run_scan(lib_dir, *tree.root_list(), time_box_minutes=60)
    started = {r.run_key: r.spec.start for r in tree.runs.values()}
    newest_first = sorted(started, key=started.get, reverse=True)
    assert sorted(seen) == sorted(newest_first[:3]), "three passes of 25 minutes start inside one hour, newest runs first"
    assert res.waiting_tier1 == 2 and res.exit_code == 0, "work left for later is not a failure"
    idx = index_of(lib_dir)
    by_id = {r["id"]: r for r in rows(idx)}
    assert [by_id[k]["detail_t1"] for k in newest_first] == ["full", "full", "full", "quick", "quick"]
    assert all(by_id[k]["verdict"] == "pass" for k in newest_first), "the runs not counted yet keep their quick verdict"
    assert by_id[newest_first[3]]["provisional"] is True
    assert idx["progress"]["state"] == "idle" and idx["progress"]["tier1"]["runs_waiting"] == 2
    assert "3 of 5 runs done" in idx["progress"]["message"]
    assert len(log_events(res, "scan.time_box")) == 1
    # the next scan goes on where this one stopped, and reads only what is left
    reads = len(parts.engine.reads)
    res = run_scan(lib_dir, time_box_minutes=60)
    assert sorted(seen) == sorted(newest_first) and res.waiting_tier1 == 0
    assert len(parts.engine.reads) == reads + 4
    idx = index_of(lib_dir)
    assert sorted(column(idx, "detail_t1")) == ["full"] * 5 + ["none"]
    assert idx["progress"]["tier1"]["runs_waiting"] == 0 and idx["progress"]["message"] is None
    assert sorted(res.changed) == sorted(newest_first[3:])


def test_time_box_zero_means_no_limit(parts, tree, lib_dir, monkeypatch):
    clock = [0.0]
    monkeypatch.setattr(scan_mod, "_monotonic", lambda: clock[0])
    parts.engine.on_pass = lambda run, role: clock.__setitem__(0, clock[0] + 100000)
    res = run_scan(lib_dir, *tree.root_list(), time_box_minutes=0)
    assert res.waiting_tier1 == 0


# ------------------------------------------------------------------------------------------------ options
def test_tier_only_since_force_and_dry_run(parts, tree, lib_dir):
    eng = parts.engine
    a1, a3, b2 = tree.runs["a1"], tree.runs["a3"], tree.runs["b2"]
    res = run_scan(lib_dir, *tree.root_list(), dry_run=True)
    assert not os.path.exists(lib_dir), "--dry-run writes nothing, not even the library"
    assert len(res.planned["new"]) == 6 and res.planned["bytes_to_count"] > 0
    res = run_scan(lib_dir, *tree.root_list(), tier="2")              # read data only: quick verdicts, no counting
    assert len(res.new) == 6 and eng.calls["full_pass"] == 0 and res.waiting_tier1 == 5
    res = run_scan(lib_dir, only=(a1.flow_cell_id,))                  # both runs of that flow cell
    assert eng.calls["full_pass"] == 2 and res.waiting_tier1 == 3
    res = run_scan(lib_dir, since="2026-04-01")                       # only a3 started after that day
    assert eng.calls["full_pass"] == 3 and res.waiting_tier1 == 2
    assert res.changed == [a3.run_key]
    res = run_scan(lib_dir, tier="1")
    assert eng.calls["full_pass"] == 5 and res.waiting_tier1 == 0
    reads = len(eng.reads)
    res = run_scan(lib_dir, only=(b2.run_key,), force=True)           # recompute although nothing changed
    assert len(eng.reads) == reads + 2 and eng.calls["full_pass"] == 6
    planned = run_scan(lib_dir, dry_run=True).planned
    assert planned["new"] == [] and len(planned["known"]) == 6
    before = files_under(lib_dir)
    run_scan(lib_dir, dry_run=True)
    assert files_under(lib_dir) == before


def test_no_hub_and_no_index(parts, tree, lib_dir):
    run_scan(lib_dir, tree.roots["seq"], no_hub=True)
    assert os.path.isfile(os.path.join(lib_dir, "index", "index.json.gz")) and not os.path.exists(os.path.join(lib_dir, "hub"))
    lib2 = lib_dir + "_2"
    run_scan(lib2, tree.roots["seq"], no_index=True)
    assert not os.path.exists(os.path.join(lib2, "index")) and not os.path.exists(os.path.join(lib2, "hub"))
    assert len(LibraryStore.open(lib2).run_keys()) == 4
    assert LibraryStore.open(lib2).index()["runs"]["n"] == 4, "the index can always be built from the cards"


def test_one_bad_run_folder_does_not_stop_the_scan(parts, tree, lib_dir):
    bad = tree.runs["b1"]
    parts.engine.fail_identify.add(os.path.basename(bad.run_dir))
    parts.engine.fail_pass.add(tree.runs["a1"].run_key)
    res = run_scan(lib_dir, *tree.root_list())
    assert res.exit_code == 5
    assert bad.run_key not in LibraryStore.open(lib_dir).run_keys() and len(LibraryStore.open(lib_dir).run_keys()) == 5
    by_id = {r["id"]: r for r in rows(index_of(lib_dir))}
    assert by_id[tree.runs["a1"].run_key]["detail_t1"] == "quick", "its quick verdict stays"
    assert by_id[tree.runs["a2"].run_key]["detail_t1"] == "full"
    assert [n["code"] for n in index_of(lib_dir)["notices"]] == ["inputs_failed"]
    # a failed pass is not tried at every scan; --retry-failed tries it
    parts.engine.fail_pass.clear()
    parts.engine.fail_identify.clear()
    calls = parts.engine.calls["full_pass"]
    res = run_scan(lib_dir)
    assert parts.engine.calls["full_pass"] == calls + 1 and res.new == [bad.run_key], "only the new run was counted"
    assert {r["id"]: r for r in rows(index_of(lib_dir))}[tree.runs["a1"].run_key]["detail_t1"] == "quick"
    run_scan(lib_dir, retry_failed=True)
    assert {r["id"]: r for r in rows(index_of(lib_dir))}[tree.runs["a1"].run_key]["detail_t1"] == "full"


# ------------------------------------------------------------------------------------------------ read data
def _add_read_files(tree, lib_dir) -> dict[str, str]:
    """Fake read-data files beside the runs, and the rules that find them."""
    a1, a3 = tree.runs["a1"], tree.runs["a3"]
    out = {}
    for length, (run, root, sub) in enumerate(((a1, "seq", "alpha_project"), (a3, "store", "group_one/alpha_project")), 900):
        d = os.path.join(tree.roots[root], sub, "analysis", "alignments")
        os.makedirs(d, exist_ok=True)
        for sample, n in (("s_one", 30), ("s_two", 50)):
            p = os.path.join(d, f"{run.run_key}__{sample}.fakebam")
            open(p, "w").write("\n".join([str(length)] * n) + "\n")
            out[f"{run.run_key}/{sample}"] = p
        trees.age_tree(os.path.join(tree.roots[root], sub, "analysis"))
    open(os.path.join(lib_dir, "flowcellar.toml"), "a").write(
        '\n[[inputs]]\nrole = "align.genome_bam"\nglob = "**/analysis/alignments/*__{sample}.fakebam"\n')
    return out


def test_read_data_comes_after_run_health_is_cached_and_survives_a_deleted_file(parts, tree, lib_dir):
    run_scan(lib_dir, *tree.root_list())
    set_root(lib_dir, "root_store", project_depth=2)
    files = _add_read_files(tree, lib_dir)
    a1, a3 = tree.runs["a1"], tree.runs["a3"]
    order = []
    parts.readdata.on_pass = lambda inp: order.append(os.path.basename(inp.candidate.path))
    res = run_scan(lib_dir)
    assert parts.readdata.calls["link"] == 1 and len(parts.readdata.reads) == 4
    assert order[0].startswith(a3.run_key) and order[-1].startswith(a1.run_key), "newest run first"
    assert sorted(res.changed) == sorted([a1.run_key, a3.run_key])
    idx = index_of(lib_dir)
    by_id = {r["id"]: r for r in rows(idx)}
    assert by_id[a1.run_key]["read_data"] == "complete" and by_id[a1.run_key]["t2_reads"] == 80
    assert by_id[a1.run_key]["data"] == "pass" and by_id[tree.runs["b1"].run_key]["data"] == "not_available"
    assert idx["library"]["counts"]["with_read_data"] == 2 and idx["library"]["counts"]["samples"] == 2
    assert {r["id"]: r["read_data_files"] for r in idx["roots"]} == {
        "root_seq": {"found": 2, "read": 2}, "root_store": {"found": 2, "read": 2}, "root_arch": {"found": 0, "read": 0}}
    assert idx["progress"]["read_data"] == {"files_total": 4, "files_read": 4, "files_waiting": 0, "bytes_waiting": 0, "eta_minutes": None}
    state = json.load(open(os.path.join(lib_dir, "state", a1.run_key + ".state.json")))
    assert sorted(state["partials_planned"]) == sorted(state["partials_merged"]) and len(state["partials_merged"]) == 4
    # nothing new: no link, no pass
    before = snapshot(lib_dir)
    res = run_scan(lib_dir)
    assert parts.readdata.calls["link"] == 1 and len(parts.readdata.reads) == 4 and res.changed == []
    assert snapshot(lib_dir) == before
    assert "read data: 4 of 4 files read" in log_events(res, "scan.done")[0]["message"]
    # a file that was read is deleted (pipelines do that): the bundle keeps what was read from it
    os.unlink(files[f"{a1.run_key}/s_one"])
    res = run_scan(lib_dir)
    by_id = {r["id"]: r for r in rows(index_of(lib_dir))}
    assert by_id[a1.run_key]["t2_reads"] == 80 and by_id[a1.run_key]["read_data"] == "complete"
    assert res.changed == [] and len(parts.readdata.reads) == 4
    assert log_events(res, "inputs.gone")[0]["n"] == 1


def test_read_data_obeys_tier_time_box_and_retry(parts, tree, lib_dir, monkeypatch):
    run_scan(lib_dir, *tree.root_list())
    files = _add_read_files(tree, lib_dir)
    a1, a3 = tree.runs["a1"], tree.runs["a3"]
    res = run_scan(lib_dir, tier="1")                                  # run health only: nothing is linked or read
    assert parts.readdata.calls["link"] == 0 and parts.readdata.reads == []
    clock = [0.0]
    monkeypatch.setattr(scan_mod, "_monotonic", lambda: clock[0])
    parts.readdata.on_pass = lambda inp: clock.__setitem__(0, clock[0] + 40 * 60)
    parts.readdata.fail.add(os.path.basename(files[f"{a3.run_key}/s_two"]))
    res = run_scan(lib_dir, time_box_minutes=60)
    # newest run first: its first file is read (40 min), its second fails; then one file of the older run fits
    # before the hour is over, and the last one is left for the next scan
    assert [os.path.basename(p) for p in parts.readdata.reads] == [f"{a3.run_key}__s_one.fakebam", f"{a1.run_key}__s_one.fakebam"]
    assert res.exit_code == 5 and res.waiting_read_data == 1
    idx = index_of(lib_dir)
    assert idx["progress"]["read_data"]["files_waiting"] == 1 and "Read data: 2 of 4 files read" in idx["progress"]["message"]
    by_id = {r["id"]: r for r in rows(idx)}
    assert by_id[a3.run_key]["read_data"] == "partial" and by_id[a3.run_key]["t2_files_waiting"] == 1
    assert by_id[a1.run_key]["read_data"] == "partial" and by_id[a1.run_key]["t2_reads"] == 30
    assert [n["code"] for n in idx["notices"]] == ["inputs_failed"]
    parts.readdata.on_pass = None
    res = run_scan(lib_dir)                                            # the rest; the failed file is left alone
    assert len(parts.readdata.reads) == 3 and res.waiting_read_data == 0
    assert {r["id"]: r for r in rows(index_of(lib_dir))}[a1.run_key]["read_data"] == "complete"
    parts.readdata.fail.clear()
    res = run_scan(lib_dir, retry_failed=True)
    assert len(parts.readdata.reads) == 4 and res.exit_code == 0
    assert {r["id"]: r for r in rows(index_of(lib_dir))}[a3.run_key]["read_data"] == "complete"


def test_the_run_folder_and_its_files_give_a_run_two_projects(parts, tree, lib_dir):
    run_scan(lib_dir, *tree.root_list())
    set_root(lib_dir, "root_store", project_depth=2)
    b1 = tree.runs["b1"]                                               # sequenced for beta, analysed under alpha
    d = os.path.join(tree.roots["store"], "group_one", "alpha_project", "analysis", "alignments")
    os.makedirs(d)
    open(os.path.join(d, f"{b1.run_key}__s_x.fakebam"), "w").write("700\n700\n")
    trees.age_tree(os.path.dirname(d))
    open(os.path.join(lib_dir, "flowcellar.toml"), "a").write(
        '\n[[inputs]]\nrole = "align.genome_bam"\nglob = "**/analysis/alignments/*__{sample}.fakebam"\n')
    run_scan(lib_dir)
    status = json.load(open(os.path.join(lib_dir, "status", b1.run_key + ".status.json")))
    assert [(p["id"], p["roots"]) for p in status["relations"]["projects"]] == [("beta_project", ["root_seq"]),
                                                                                ("alpha_project", ["root_store"])]
    row = {r["id"]: r for r in rows(index_of(lib_dir))}[b1.run_key]
    assert row["project"] == "beta_project" and row["project_ids"] == ["beta_project", "alpha_project"]


# ------------------------------------------------------------------------------------------------ runs in progress, one run by hand
def test_a_run_in_progress_is_looked_at_again_at_every_scan_until_it_has_finished(parts, tmp_path, lib_dir):
    t = trees.make_three_roots(tmp_path / "live", aged=False)
    a1 = t.runs["a1"]
    final = a1.files["minknow.final_summary"]
    held = str(tmp_path / "final_summary.held")
    os.rename(final, held)                                             # MinKNOW writes it at the very end
    res = run_scan(lib_dir, os.path.dirname(os.path.dirname(a1.run_dir)))
    assert res.new == [], "a root that is not marked live: files this fresh are left for the next scan"
    set_root(lib_dir, "alpha_project", live=True)
    res = run_scan(lib_dir)
    assert sorted(res.new) == sorted([a1.run_key, t.runs["a2"].run_key])
    state = json.load(open(os.path.join(lib_dir, "state", a1.run_key + ".state.json")))
    assert state["run_state"] == "running"
    quick = parts.engine.calls["quick_bundle"]
    run_scan(lib_dir)
    assert parts.engine.calls["quick_bundle"] == quick + 1, "only the running run is computed again"
    os.rename(held, final)
    res = run_scan(lib_dir)
    state = json.load(open(os.path.join(lib_dir, "state", a1.run_key + ".state.json")))
    assert state["run_state"] == "finished" and state["t1"] == "done"
    quick = parts.engine.calls["quick_bundle"]
    trees.age_tree(t.base)
    run_scan(lib_dir)
    run_scan(lib_dir)
    assert parts.engine.calls["quick_bundle"] == quick, "finished: left alone from now on"


def test_a_raw_only_folder_becomes_a_run_when_its_reports_arrive(parts, tmp_path, lib_dir):
    t = trees.make_three_roots(tmp_path / "data")
    b1 = t.runs["b1"]
    away = tmp_path / "reports_not_copied_yet"
    away.mkdir()
    names = [n for n in os.listdir(b1.run_dir) if os.path.isfile(os.path.join(b1.run_dir, n))]
    for n in names:
        shutil.move(os.path.join(b1.run_dir, n), away / n)
    shutil.move(os.path.join(b1.run_dir, "other_reports"), away / "other_reports")
    os.makedirs(os.path.join(b1.run_dir, "other_reports"))
    trees.age_tree(b1.run_dir)
    res = run_scan(lib_dir, t.roots["seq"])
    by_id = {r["id"]: r for r in rows(index_of(lib_dir))}
    assert b1.run_key in res.new and by_id[b1.run_key]["verdict"] == "incomplete" and by_id[b1.run_key]["cause"] == "raw_only"
    os.rmdir(os.path.join(b1.run_dir, "other_reports"))
    shutil.move(str(away / "other_reports"), os.path.join(b1.run_dir, "other_reports"))
    for n in names:
        shutil.move(str(away / n), os.path.join(b1.run_dir, n))
    trees.age_tree(b1.run_dir)
    res = run_scan(lib_dir)
    by_id = {r["id"]: r for r in rows(index_of(lib_dir))}
    assert res.new == [] and res.changed == [b1.run_key] and len(by_id) == 4
    assert by_id[b1.run_key]["verdict"] == "pass" and by_id[b1.run_key]["detail_t1"] == "full"
    assert by_id[b1.run_key]["reads"] == b1.expected["reads"]


def test_one_run_added_by_hand(parts, tree, lib_dir, tmp_path):
    from flowcellar.library import api
    from flowcellar.testing.fixtures import RunSpec, make_run
    run_scan(lib_dir, tree.roots["store"])
    lib = api.open_library(lib_dir)
    outside = make_run(tmp_path / "elsewhere", RunSpec(reads=120, hours=1.0, flow_cell_id="TST00042", seed=42,
                                                         start="2026-05-01T08:00:00Z"))
    before = snapshot(lib_dir)
    key = api.add_run(lib, outside.run_dir)                            # a folder under no root of the library
    assert key == outside.run_key
    idx = index_of(lib_dir)
    row = {r["id"]: r for r in rows(idx)}[key]
    assert row["detail_t1"] == "full" and row["roots"] == [] and row["project"] is None and row["reads"] == 120
    assert idx["runs"]["n"] == 2 and os.path.isfile(os.path.join(lib_dir, "hub", "data", key + ".js"))
    assert all(snapshot(lib_dir)[n] == m for n, m in before.items() if n.startswith(("bundles/", "status/"))), "nothing else changed"
    res = run_scan(lib_dir)                                            # the next scan leaves it alone: its folder is still there
    assert res.missing == [] and res.new == [] and index_of(lib_dir)["runs"]["n"] == 2
    a1 = tree.runs["a1"]                                               # a folder under a root the library does not scan yet
    quick_only = api.add_run(lib, a1.run_dir, quick=True)
    row = {r["id"]: r for r in rows(index_of(lib_dir))}[quick_only]
    assert row["detail_t1"] == "quick" and row["provisional"] is True
    assert api.add_run(lib, a1.run_dir) == a1.run_key                  # again, now with the read count
    assert {r["id"]: r for r in rows(index_of(lib_dir))}[a1.run_key]["detail_t1"] == "full"
    shutil.rmtree(outside.run_dir)
    assert run_scan(lib_dir).missing == [key]
    with pytest.raises(Exception) as e:
        api.add_run(lib, str(tmp_path))
    assert getattr(e.value, "exit_code", None) == 4


def test_a_run_that_could_not_be_judged_is_judged_by_the_next_scan(parts, tree, lib_dir, monkeypatch):
    from flowcellar.engine import api as engine_api
    a3 = tree.runs["a3"]
    real = engine_api.judge

    def broken(bundle, history, **kw):
        if bundle["bundle_id"] == a3.run_key:
            raise RuntimeError("limits table unreadable (test)")
        return real(bundle, history, **kw)
    monkeypatch.setattr(engine_api, "judge", broken)
    res = run_scan(lib_dir, *tree.root_list())
    assert res.exit_code == 5 and index_of(lib_dir)["runs"]["n"] == 5, "the run has a bundle but no verdict yet"
    assert os.path.isfile(os.path.join(lib_dir, "bundles", a3.run_key + ".fcb.json.gz"))
    monkeypatch.setattr(engine_api, "judge", real)
    res = run_scan(lib_dir)
    assert res.exit_code == 0 and index_of(lib_dir)["runs"]["n"] == 6
    assert {r["id"]: r for r in rows(index_of(lib_dir))}[a3.run_key]["verdict"] == "pass"


def test_only_with_an_id_nobody_has_says_so(parts, tree, lib_dir):
    res = run_scan(lib_dir, tree.roots["seq"], only=("TST00099",))
    assert res.new == [] and [e["message"] for e in log_events(res, "scan.only")] == [
        "--only TST00099: no run of this library or under its roots has this run key or flow cell id"]


def test_a_root_taken_out_of_the_configuration_does_not_make_its_runs_disappear(parts, tree, lib_dir):
    run_scan(lib_dir, tree.roots["seq"], tree.roots["store"])
    path = os.path.join(lib_dir, "flowcellar.toml")
    text = open(path).read()
    cut = text.index('[[roots]]\nid = "root_store"')
    end = text.index("project_depth = 1", cut) + len("project_depth = 1")
    open(path, "w").write(text[:cut] + text[end:])
    res = run_scan(lib_dir)
    assert res.missing == [] and index_of(lib_dir)["runs"]["n"] == 5
    assert [r["id"] for r in index_of(lib_dir)["roots"]] == ["root_seq"]


def test_of_two_copies_of_a_run_the_complete_one_is_read(parts, tree, lib_dir):
    b1 = tree.runs["b1"]
    partial_copy = trees.mirror_run(b1, os.path.join(tree.roots["seq"], "a_first_in_order", "pool_b"))
    os.unlink(os.path.join(partial_copy, os.path.basename(b1.files["minknow.sequencing_summary"])))   # a copy still in transit
    trees.age_tree(tree.roots["seq"])
    run_scan(lib_dir, tree.roots["seq"])
    state = json.load(open(os.path.join(lib_dir, "state", b1.run_key + ".state.json")))
    assert [p["path"] for p in state["paths"]] == [b1.run_dir, partial_copy]
    assert {r["id"]: r for r in rows(index_of(lib_dir))}[b1.run_key]["reads"] == b1.expected["reads"]


def test_the_hub_is_rewritten_as_results_land_but_not_more_than_once_a_minute(parts, tree, lib_dir, monkeypatch):
    clock = [0.0]
    monkeypatch.setattr(scan_mod, "_monotonic", lambda: clock[0])
    seen = set()

    def tick(run, role):                                              # every run's counting takes ten minutes
        if run.run_key not in seen:
            seen.add(run.run_key)
            clock[0] += 600
    parts.engine.on_pass = tick
    run_scan(lib_dir, *tree.root_list())
    started = {r.run_key: r.spec.start for r in tree.runs.values()}
    newest_first = sorted(started, key=started.get, reverse=True)
    writes = parts.report.hub_writes
    assert writes[0] is None, "first the whole hub, before any counting"
    assert writes[1:-1] == [[k] for k in newest_first], "then the data file of each run as its count lands"
    assert writes[-1] is None, "and the whole hub once more at the end"
    lib2 = lib_dir + "_fast"
    parts.report.hub_writes.clear()
    parts.engine.on_pass = lambda run, role: clock.__setitem__(0, clock[0] + 1)      # two passes a run, one second each
    run_scan(lib2, *tree.root_list())
    assert parts.report.hub_writes == [None, None], "fast passes: the index is not rewritten after every run"
    idx = index_of(lib2)
    assert idx["progress"]["state"] == "idle" and sorted(column(idx, "detail_t1")) == ["full"] * 5 + ["none"]
