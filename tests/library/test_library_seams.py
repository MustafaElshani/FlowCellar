"""What the integration added to the scan: a quick-only scan, linking after the first hub, loadings counted again
at every index build, and the ``[baselines]`` settings handed to the status engine. Stand-ins as in the other
library tests; the real chain is in tests/integration."""
from __future__ import annotations

import gzip
import json
import os

from library_helpers import log_events, run_scan

from flowcellar.library import cards as cards_mod
from flowcellar.library.index import build_index, column, rows


def index_of(lib_dir: str) -> dict:
    with gzip.open(os.path.join(lib_dir, "index", "index.json.gz")) as fh:
        return json.load(fh)


def test_a_quick_scan_judges_every_run_and_starts_no_pass(parts, tree, lib_dir, capsys):
    res = run_scan(lib_dir, *tree.root_list(), quick=True)
    assert res.exit_code == 0 and len(res.new) == 6 and res.waiting_tier1 == 5
    assert parts.engine.calls["quick_bundle"] == 6 and parts.engine.calls["full_pass"] == 0 and parts.engine.reads == []
    assert parts.readdata.calls["link"] == 0 and parts.readdata.calls["read_pass"] == 0
    idx = index_of(lib_dir)
    assert sorted(column(idx, "detail_t1")) == ["none"] + ["quick"] * 5 and all(column(idx, "headline"))
    assert os.path.isfile(os.path.join(lib_dir, "hub", "index.html"))
    done = log_events(res, "scan.done")[0]["message"]
    assert "the reads of 5 runs are still to be counted" in done and "read data not looked at (--quick)" in done
    # the next scan without --quick counts them
    res = run_scan(lib_dir)
    assert res.waiting_tier1 == 0 and parts.engine.calls["full_pass"] == 5
    assert sorted(column(index_of(lib_dir), "detail_t1")) == ["full"] * 5 + ["none"]


def test_files_are_linked_only_after_the_first_hub_is_written(parts, tree, lib_dir):
    """Linking opens every read-data file it has not seen; the verdicts from the small files must not wait for it."""
    order = []
    from flowcellar.readdata import api as readdata_api
    from flowcellar.report import api as report_api
    real_link, real_hub = readdata_api.link, report_api.write_hub

    def link(*a, **kw):
        order.append("link")
        return real_link(*a, **kw)

    def hub(*a, **kw):
        order.append("hub")
        return real_hub(*a, **kw)

    a3 = tree.runs["a3"]
    d = os.path.join(tree.roots["store"], "group_one", "alpha_project", "analysis", "alignments")
    os.makedirs(d)
    with open(os.path.join(d, f"{a3.run_key}__s_one.fakebam"), "w") as fh:
        fh.write("100\n200\n")
    from flowcellar.testing import trees
    trees.age_tree(tree.base)
    # a configuration a person wrote before the first scan: the roots given on the command line are added to it
    os.makedirs(lib_dir)
    with open(os.path.join(lib_dir, "flowcellar.toml"), "w") as fh:
        fh.write('[[inputs]]\nrole = "align.genome_bam"\nglob = "**/analysis/alignments/*__{sample}.fakebam"\n')
    readdata_api.link, report_api.write_hub = link, hub
    try:
        run_scan(lib_dir, *tree.root_list())
    finally:
        readdata_api.link, report_api.write_hub = real_link, real_hub
    assert "link" in order and order.index("hub") < order.index("link")
    idx = index_of(lib_dir)
    assert [r["id"] for r in idx["roots"]] == ["root_seq", "root_store", "root_arch"], "the roots were added to the file"
    by_id = {r["id"]: r for r in rows(idx)}
    assert by_id[a3.run_key]["read_data"] == "complete" and by_id[a3.run_key]["t2_reads"] == 2


def _card(key, fc, started, use, of):
    return {"run_key": key, "columns": {"id": key, "fc": fc, "started": started, "use": use, "use_of": of, "use_source": "library",
                                         "restart": False, "restart_of": None, "roots": [], "project_ids": []},
            "facts": {"run_key": key, "flow_cell_id": fc, "started": started, "run_hours": 30.0, "estimated_bases": 5e9}}


def test_the_index_counts_the_loadings_of_a_flow_cell_as_of_today():
    """A status file is frozen when its run is judged: the first loading of a flow cell says "1 of 1" for ever.
    The index is library knowledge and says how many loadings there are now."""
    first = _card("TST00001_0a0a0a0a", "TST00001", "2026-03-02T09:00:00Z", 1, 1)
    second = _card("TST00001_0b0b0b0b", "TST00001", "2026-03-04T09:00:00Z", 2, 2)
    other = _card("TST00002_0c0c0c0c", "TST00002", "2026-03-05T09:00:00Z", 1, 1)
    by_id = {r["id"]: r for r in rows(build_index([first, second, other]))}
    assert (by_id[first["run_key"]]["use"], by_id[first["run_key"]]["use_of"]) == (1, 2)
    assert (by_id[second["run_key"]]["use"], by_id[second["run_key"]]["use_of"]) == (2, 2)
    assert (by_id[other["run_key"]]["use"], by_id[other["run_key"]]["use_of"]) == (1, 1)
    # what a person stated about a run is not counted over
    stated = _card("TST00003_0d0d0d0d", "TST00003", "2026-03-06T09:00:00Z", 3, 3)
    stated["columns"]["use_source"] = "declared"
    row = rows(build_index([stated]))[0]
    assert (row["use"], row["use_of"]) == (3, 3)


def test_the_baseline_settings_of_the_configuration_reach_the_status_engine():
    cfg = {"profiles": {"cdna": {"rules": {"fc_occupancy_first_hour": {"warn": 0.6}}}},
           "baselines": {"window": 12, "min_runs": 3, "trusted_runs": 6, "spread_guard": 0.3, "groups": [{"label": "g", "match": {"kit": "K"}}]},
           "flow_cell": {"restart_max_hours": 1}}
    kw = cards_mod.judge_settings(cfg)
    assert kw["profile"]["_baselines"] == {"window": 12, "min_runs": 3, "trusted_runs": 6, "spread_guard": 0.3}
    assert kw["profile"]["cdna"] == cfg["profiles"]["cdna"] and kw["baseline_groups"] == ({"label": "g", "match": {"kit": "K"}},)
    assert kw["flow_cell"] == {"restart_max_hours": 1}
    assert cards_mod.judge_settings({})["profile"] == {"_baselines": {}}


def test_a_read_table_that_is_cut_off_is_counted_once_and_does_not_wait_for_ever(parts, tree, lib_dir):
    """Found on real data: one sequencing summary on disk ended in the middle of a line. The pass says so
    (state "truncated"); the run must not be announced as "still to be counted" at every later scan."""
    from dataclasses import replace

    from flowcellar.engine import api as engine_api
    cut = tree.runs["b1"].run_key
    real_pass = engine_api.full_pass

    def full_pass(run, partials_dir, opts=None):
        refs = real_pass(run, partials_dir, opts)
        if run.run_key != cut:
            return refs
        return [replace(r, state="truncated", complete=False) if r.role == "minknow.sequencing_summary" else r for r in refs]

    engine_api.full_pass = full_pass
    try:
        res = run_scan(lib_dir, *tree.root_list())
        assert res.exit_code == 0 and res.waiting_tier1 == 0, "nothing is left to count: the file is what it is"
        done = log_events(res, "scan.done")[0]["message"]
        assert "1 run has a read table that is cut off" in done and "still to be counted" not in done
        assert len(log_events(res, "run.truncated")) == 1
        with open(os.path.join(lib_dir, "state", cut + ".state.json")) as fh:
            assert json.load(fh)["t1"] == "truncated"
        assert index_of(lib_dir)["progress"]["tier1"]["runs_waiting"] == 0
        calls = parts.engine.calls["full_pass"]
        res = run_scan(lib_dir)
        assert parts.engine.calls["full_pass"] == calls and res.changed == [], "the next scan does not try again"
    finally:
        engine_api.full_pass = real_pass


def test_adding_one_run_by_hand_keeps_what_the_last_scan_said_about_the_roots(parts, tree, lib_dir):
    from flowcellar.library import scan as scan_mod
    from flowcellar.testing import trees
    run_scan(lib_dir, *tree.root_list())
    trees.take_offline(tree.roots["arch"], "empty")
    try:
        assert run_scan(lib_dir).exit_code == 5
        assert [n["code"] for n in index_of(lib_dir)["notices"]] == ["root_offline"]
        scan_mod.Scan(lib_dir, scan_mod.ScanOptions(wait=False)).add_single(tree.runs["a1"].run_dir)
        assert [n["code"] for n in index_of(lib_dir)["notices"]] == ["root_offline"], "no root was looked at: the notice stands"
    finally:
        trees.bring_back(tree.roots["arch"])


def test_a_recompute_takes_back_the_notice_that_asked_for_it(parts, tree, lib_dir):
    from flowcellar.library import maintain
    from flowcellar.library.store import LibraryStore
    run_scan(lib_dir, *tree.root_list())
    store = LibraryStore.open(lib_dir)
    notices = store.read_aux(store.layout.notices, {})
    notices["notices"] = [{"code": "judged_before_older_runs", "level": "info", "root": None, "since": None, "message": "2 runs ..."},
                          {"code": "root_offline", "level": "warn", "root": "root_arch", "since": None, "message": "kept"}]
    store.write_aux(store.layout.notices, notices)
    out = maintain.recompute_status(store, store.run_keys(), rebaseline=True, wait=False)
    assert out["failed"] == {} and len(out["judged"]) == 6
    assert [n["code"] for n in index_of(lib_dir)["notices"]] == ["root_offline"]


def test_a_run_gets_its_read_data_as_soon_as_its_own_files_are_read(parts, tree, lib_dir):
    """Read data of the newest run must be on its page while the files of older runs are still being read, not
    only when the whole scan has ended."""
    from library_helpers import set_root
    from test_library_scan import _add_read_files
    run_scan(lib_dir, *tree.root_list())
    set_root(lib_dir, "root_store", project_depth=2)
    _add_read_files(tree, lib_dir)
    res = run_scan(lib_dir)
    a1, a3 = tree.runs["a1"].run_key, tree.runs["a3"].run_key
    events = [(e["event"], e.get("run_key") or os.path.basename(e.get("path", ""))) for e in log_events(res)
              if e["event"] in ("read_pass.done", "bundle.written")]
    newest_written = events.index(("bundle.written", a3))
    first_read_of_older = next(i for i, (ev, what) in enumerate(events) if ev == "read_pass.done" and what.startswith(a1))
    assert newest_written < first_read_of_older, events
    assert log_events(res, "scan.linking"), "the person is told that files are being linked"
