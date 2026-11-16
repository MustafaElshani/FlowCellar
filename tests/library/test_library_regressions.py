"""Regressions, each held by a test on fixtures.

These run the real engine and the real hub writer (no stand-ins): the faults were in how the parts meet.
"""
from __future__ import annotations

import base64
import gzip
import json
import logging
import os
import signal
import socket
import threading
import time

import pytest
from library_helpers import log_events, run_scan

from flowcellar.cli import main
from flowcellar.library import lease as lease_mod
from flowcellar.library import scan as scan_mod
from flowcellar.library.index import build_index, rows
from flowcellar.library.lease import Lease, LeaseBusy, read_holder
from flowcellar.library.store import LibraryStore
from flowcellar.testing import minknow_variants as mv
from flowcellar.testing import trees
from flowcellar.testing.fixtures import make_run
from flowcellar.util.errors import FlowCellarError

needs_permissions = pytest.mark.skipif(hasattr(os, "geteuid") and os.geteuid() == 0, reason="root reads every file")


@pytest.fixture(autouse=True)
def _viewer_built():
    from flowcellar.report import api as report_api
    try:
        report_api.viewer_manifest()
    except FlowCellarError as e:
        pytest.skip(f"the viewer is not built: {e.message}")


def index_of(lib_dir: str) -> dict:
    with gzip.open(os.path.join(lib_dir, "index", "index.json.gz")) as fh:
        return json.load(fh)


def hub_index_of(lib_dir: str) -> dict:
    """What the page loads: ``hub/index.js`` is ``__FC_INDEX__("<base64 of gzip>");``."""
    text = open(os.path.join(lib_dir, "hub", "index.js")).read()
    b64 = text[text.index('("') + 2:text.rindex('")')]
    return json.loads(gzip.decompress(base64.b64decode(b64)))


def state_of(lib_dir: str, key: str) -> dict:
    return json.load(open(os.path.join(lib_dir, "state", f"{key}.state.json")))


def one_root(tmp_path, n: int = 1, **spec):
    """A root with ``n`` runs, each in a project folder of its own: ``project_1/pool/<run>``."""
    root = tmp_path / "root"
    runs = []
    for i in range(1, n + 1):
        runs.append(make_run(root / f"project_{i}" / "pool", flow_cell_id=f"TST0010{i}", seed=10 + i, reads=600, hours=3.0, **spec))
    trees.age_tree(root)
    return str(root), runs


def summary_of(run) -> str:
    return next(os.path.join(run.run_dir, n) for n in os.listdir(run.run_dir) if n.startswith("sequencing_summary"))


# ------------------------------------------------------------------------------------- a read table that cannot be opened
@needs_permissions
def test_a_run_whose_read_table_cannot_be_opened_is_not_counted_and_is_tried_again(tmp_path, caplog):
    root, (run,) = one_root(tmp_path)
    lib = str(tmp_path / "library")
    path = summary_of(run)
    os.chmod(path, 0)
    try:
        with caplog.at_level(logging.WARNING, logger="flowcellar"):
            res = run_scan(lib, root)
        assert res.exit_code == 5, "a file that was skipped is exit code 5, never a clean finish"
        assert any("could not be opened" in r.getMessage() and run.run_key in r.getMessage() for r in caplog.records)
        st = state_of(lib, run.run_key)
        assert st["t1"] == "waiting", "not counted is not 'done': the next scan must try again"
        by_role = {i["role"]: i for i in st["inputs"]}
        assert by_role["minknow.sequencing_summary"]["pass"] == "error"
        idx = index_of(lib)
        assert idx["progress"]["tier1"]["runs_waiting"] == 1 and "The next scan goes on from here." in idx["progress"]["message"]
        row = rows(idx)[0]
        assert row["detail_t1"] == "quick" and row["provisional"] is True and res.waiting_tier1 == 1
    finally:
        os.chmod(path, 0o644)
    again = run_scan(lib)                                  # readable again: counted without --force or --retry-failed
    assert again.exit_code == 0 and again.waiting_tier1 == 0
    assert state_of(lib, run.run_key)["t1"] == "done"
    row = rows(index_of(lib))[0]
    assert row["detail_t1"] == "full" and row["reads"] == run.expected["reads"] and row["provisional"] is False


def test_full_pass_answers_for_an_input_it_could_not_open(tmp_path):
    from flowcellar.engine import api as engine_api
    run = make_run(tmp_path / "runs", reads=300, hours=2.0)
    ident = engine_api.identify(run.run_dir)
    os.rename(summary_of(run), summary_of(run) + ".moved")          # the identity still names the file
    refs = engine_api.full_pass(ident, str(tmp_path / "partials"))
    by_role = {r.role: r for r in refs}
    bad = by_role["minknow.sequencing_summary"]
    assert bad.state == "error" and not bad.complete and not bad.partial_id, "an error, with no cached result behind it"
    assert by_role["minknow.pore_scan"].complete


# ------------------------------------------------------------------------------------- folders that cannot be read
@needs_permissions
def test_a_folder_that_cannot_be_read_is_reported_and_what_is_known_below_it_stays(tmp_path, caplog):
    root, (r1, r2) = one_root(tmp_path, n=2)
    lib = str(tmp_path / "library")
    locked = os.path.join(root, "project_2")
    os.chmod(locked, 0)
    try:
        with caplog.at_level(logging.WARNING, logger="flowcellar"):
            res = run_scan(lib, root, quick=True)
        assert res.exit_code == 5 and res.new == [r1.run_key]
        assert any("could not be read" in r.getMessage() and locked in r.getMessage() for r in caplog.records)
        notices = index_of(lib)["notices"]
        assert [n["code"] for n in notices] == ["folders_unreadable"] and locked in notices[0]["message"]
    finally:
        os.chmod(locked, 0o755)
    res = run_scan(lib, quick=True)                         # readable: the second run arrives, the notice goes
    assert res.exit_code == 0 and res.new == [r2.run_key] and index_of(lib)["notices"] == []
    os.chmod(locked, 0)
    try:
        res = run_scan(lib, quick=True)
        assert res.missing == [] and res.exit_code == 5, "a run under a folder that cannot be read is not 'gone'"
        by_id = {r["id"]: r for r in rows(index_of(lib))}
        assert by_id[r2.run_key]["source_state"] == "ok"
        assert not log_events(res, "run.missing")
    finally:
        os.chmod(locked, 0o755)


def test_a_run_added_by_hand_under_a_pruned_folder_is_not_declared_gone(tmp_path):
    root, (r1,) = one_root(tmp_path)
    hidden = make_run(os.path.join(root, "project_1", "work", "pool"), flow_cell_id="TST00199", seed=77, reads=300, hours=2.0)
    trees.age_tree(root)
    lib = str(tmp_path / "library")
    first = run_scan(lib, root, quick=True)
    assert first.new == [r1.run_key], "'work' is a name the walk does not enter"
    key = scan_mod.Scan(lib, scan_mod.ScanOptions(wait=False)).add_single(hidden.run_dir, quick=True)
    assert key == hidden.run_key
    res = run_scan(lib, quick=True)
    assert res.missing == [] and res.exit_code == 0
    idx = index_of(lib)
    by_id = {r["id"]: r for r in rows(idx)}
    assert by_id[key]["source_state"] == "ok"
    note = next(n for n in idx["notices"] if n["code"] == "runs_not_walked")
    assert key in note["message"] and "'work'" in note["message"] and note["level"] == "info"
    # and a folder that really went away is still noticed
    import shutil
    shutil.rmtree(hidden.run_dir)
    res = run_scan(lib, quick=True)
    assert res.missing == [key]


# ------------------------------------------------------------------------------------- interrupt
def _interrupt_during_first_pass(monkeypatch, how):
    from flowcellar.engine import api as engine_api
    real = engine_api.full_pass
    calls = []

    def full_pass(*a, **k):
        calls.append(1)
        if len(calls) == 1:
            return real(*a, **k)
        how()
        time.sleep(5)                                       # a SIGTERM arrives while the pass waits
        raise AssertionError("the signal did not interrupt the scan")
    monkeypatch.setattr(engine_api, "full_pass", full_pass)


def _raise_interrupt():
    raise KeyboardInterrupt


@pytest.mark.parametrize("how, error", [(_raise_interrupt, KeyboardInterrupt),
                                        (lambda: os.kill(os.getpid(), signal.SIGTERM), scan_mod.Terminated)])
def test_after_an_interrupt_the_hub_no_longer_says_a_scan_is_running(tmp_path, monkeypatch, how, error):
    root, runs = one_root(tmp_path, n=2)
    lib = str(tmp_path / "library")
    before = signal.getsignal(signal.SIGTERM)
    with monkeypatch.context() as patched:
        _interrupt_during_first_pass(patched, how)
        with pytest.raises(error):
            run_scan(lib, root)
    assert signal.getsignal(signal.SIGTERM) is before, "the handler of the caller is put back"
    assert not os.path.exists(os.path.join(lib, "claims", "scan.lease")), "the lease is released"
    assert index_of(lib)["progress"]["state"] == "idle"
    hub = hub_index_of(lib)
    assert hub["progress"]["state"] == "idle", "the page would reload a 'scanning' index every minute for ever"
    assert hub["written"] == index_of(lib)["written"]
    assert hub["progress"]["tier1"]["runs_waiting"] == 1 and sorted(r["detail_t1"] for r in rows(hub)) == ["full", "quick"]
    res = run_scan(lib)                                     # and the next scan simply goes on
    assert res.exit_code == 0 and res.waiting_tier1 == 0


# ------------------------------------------------------------------------------------- lease
def test_a_lease_of_this_host_whose_pid_now_is_another_process_is_taken_over(tmp_path):
    path = str(tmp_path / "claims" / "scan.lease")
    os.makedirs(os.path.dirname(path))
    # pid 1 is always alive and never a scan; the lease was last touched two days ago
    json.dump({"host": socket.gethostname(), "pid": 1, "started": "2026-01-01T00:00:00Z"}, open(path, "w"))
    old = time.time() - 2 * 86400
    os.utime(path, (old, old))
    with Lease(path, wait=False) as lease:
        assert lease.taken_over is not None and lease.taken_over.pid == 1
    # the same file, fresh: a live scan touches it every minute, so it is busy, and the message says how to recover
    json.dump({"host": socket.gethostname(), "pid": 1, "started": "2026-01-01T00:00:00Z"}, open(path, "w"))
    with pytest.raises(LeaseBusy) as e:
        Lease(path, wait=False).acquire()
    assert path in e.value.hint and "last touched" in e.value.hint and "delete that file" in e.value.hint
    os.unlink(path)


def test_a_lease_that_records_another_start_time_for_its_pid_is_stale_at_once(tmp_path):
    if lease_mod._pid_started(os.getpid()) is None:
        pytest.skip("this system does not tell when a process started")
    path = str(tmp_path / "scan.lease")
    born = lease_mod._pid_started(os.getpid())
    json.dump({"host": socket.gethostname(), "pid": os.getpid(), "started": "2026-01-01T00:00:00Z", "pid_started": born - 3600},
              open(path, "w"))
    with Lease(path, wait=False) as lease:                  # fresh mtime, live pid, but not the process that wrote it
        assert lease.taken_over.pid == os.getpid()
    json.dump({"host": socket.gethostname(), "pid": os.getpid(), "started": "2026-01-01T00:00:00Z", "pid_started": born},
              open(path, "w"))
    with pytest.raises(LeaseBusy):
        Lease(path, wait=False).acquire()


def test_an_unreadable_lease_is_named_as_such(tmp_path):
    path = str(tmp_path / "scan.lease")
    open(path, "w").close()
    old = time.time() - 3600
    os.utime(path, (old, old))
    with Lease(path, wait=False) as lease:
        assert lease.taken_over.words() == "a lease file that could not be read"
        assert "process 0" not in lease.taken_over.words()
    assert read_holder(path) is None


def test_run_into_a_library_says_that_it_waits_for_a_scan(tmp_path, caplog):
    root, (run,) = one_root(tmp_path)
    lib = str(tmp_path / "library")
    run_scan(lib, root, quick=True)
    other = make_run(tmp_path / "elsewhere", flow_cell_id="TST00177", seed=5, reads=300, hours=2.0)
    held = Lease(os.path.join(lib, "claims", "scan.lease")).acquire()
    threading.Timer(0.4, held.release).start()
    with caplog.at_level(logging.INFO, logger="flowcellar"):
        key = scan_mod.Scan(lib, scan_mod.ScanOptions(lease_poll_seconds=0.05)).add_single(other.run_dir, quick=True)
    assert key == other.run_key
    assert any("another scan of this library is running" in r.getMessage() and "waiting for it" in r.getMessage()
               for r in caplog.records)


def test_the_time_box_starts_when_the_lease_is_held(tmp_path, monkeypatch):
    root, (run,) = one_root(tmp_path)
    lib = str(tmp_path / "library")
    clock = {"t": 1000.0}
    monkeypatch.setattr(scan_mod, "_monotonic", lambda: clock["t"])
    real = Lease.acquire

    def slow_acquire(self, heartbeat=True):
        clock["t"] += 3600.0                                # an hour behind another scan
        return real(self, heartbeat)
    monkeypatch.setattr(Lease, "acquire", slow_acquire)
    res = run_scan(lib, root, time_box_minutes=1)
    assert res.waiting_tier1 == 0 and state_of(lib, run.run_key)["t1"] == "done", "the wait is not taken out of the box"
    assert not log_events(res, "scan.time_box")


# ------------------------------------------------------------------------------------- an empty library
def test_ls_on_a_library_without_runs(tmp_path, capsys):
    lib = str(tmp_path / "library")
    for argv in (["init", lib], ["ls", "-L", lib]):
        try:
            main(argv)
        except SystemExit as e:
            assert e.code == 0, capsys.readouterr().err
    out, err = capsys.readouterr()
    assert out.split() == ["id", "started", "verdict", "headline"] and "0 runs" in err
    with pytest.raises(SystemExit) as e:
        main(["ls", "-L", lib, "--format", "json"])
    assert e.value.code == 0 and json.loads(capsys.readouterr().out) == []


def test_an_index_without_runs_has_every_column():
    from flowcellar.library.index import COLUMN_ORDER
    idx = build_index([])
    assert idx["runs"]["n"] == 0 and list(idx["runs"]["columns"]) == list(COLUMN_ORDER)
    assert all(v == [] for v in idx["runs"]["columns"].values())


# ------------------------------------------------------------------------------------- a damaged read table
def test_a_read_table_of_bytes_that_are_not_text_gives_no_traceback_and_is_named_on_the_run(tmp_path, capfd):
    root, (run,) = one_root(tmp_path)
    path = summary_of(run)
    size = os.path.getsize(path)
    with open(path, "wb") as fh:
        fh.write(bytes(range(256)) * (size // 256 + 1))
    trees.age_tree(root)
    lib = str(tmp_path / "library")
    res = run_scan(lib, root)
    err = capfd.readouterr().err
    assert "Traceback" not in err and "Exception ignored" not in err and "UnicodeDecodeError" not in err
    assert res.exit_code == 5
    store = LibraryStore.open(lib)
    bundle = store.composed(run.run_key)
    assert bundle["tiers"]["run_health"]["state"] == "partial", "a read table that could not be read is not a complete picture"
    codes = [m["code"] for m in bundle["status"]["modifiers"]]
    assert "reads_not_counted" in codes
    assert "reads_not_counted" in rows(index_of(lib))[0]["modifiers"]


def test_the_fast_reader_leaves_bytes_that_are_not_text_to_the_slow_one():
    from flowcellar.sources import seqsummary
    header = ["read_id", "channel", "mux", "start_time", "duration", "minknow_events", "end_reason"]
    good = b"r1\t1\t1\t0.5\t1.0\t100\tsignal_positive\n"
    cols, bad = seqsummary.parse_block(good + b"\xff\xfe\x00junk\tno\n" + good, header)
    assert len(cols["channel"]) == 2 and bad == 1


# ------------------------------------------------------------------------------------- words
def test_plural_of_the_run_that_still_waits(tmp_path, caplog):
    root, (run,) = one_root(tmp_path)
    lib = str(tmp_path / "library")
    with caplog.at_level(logging.INFO, logger="flowcellar"):
        run_scan(lib, root, quick=True)
    text = " ".join(r.getMessage() for r in caplog.records)
    assert "the reads of 1 run are" in text and "1 run is" not in text


def test_validate_takes_the_library_option(tmp_path, capsys):
    root, _ = one_root(tmp_path)
    lib = str(tmp_path / "library")
    run_scan(lib, root, quick=True)
    with pytest.raises(SystemExit) as e:
        main(["validate", "-L", lib])
    out, err = capsys.readouterr()
    assert e.value.code == 0 and "1 of 1 valid" in err
    assert mv is not None


# ------------------------------------------------------------------------------------- projects
def test_a_pipeline_folder_at_project_depth_is_not_a_project():
    from flowcellar.config.loader import not_project_names
    from flowcellar.library.projects import project_name
    cfg = {"inputs": [{"role": "align.genome_bam", "glob": "**/analysis/alignments/{sample}_sorted.bam"},
                      {"role": "dorado.bam", "glob": "**/*_basecalled/{sample}.bam"}],
           "project": {"exclude": ["Scratch"]}}
    skip = not_project_names(cfg)
    assert skip == {"analysis", "alignments", "scratch"}
    root = {"id": "store", "path": "/data", "project_depth": 2}
    assert project_name(root, "group/alpha_project/analysis/alignments/s1_sorted.bam", skip) == "alpha_project"
    # the same pipeline one level higher: "analysis" stands at project depth and names no project
    assert project_name(root, "alpha_project/analysis/alignments/s1_sorted.bam", skip) is None
    assert project_name(root, "alpha_project/analysis/alignments/s1_sorted.bam") == "analysis", "without the rule it would"
    assert project_name(root, "group/scratch/x/run_folder", skip) is None
    assert project_name({"id": "r", "path": "/d", "project_pattern": "*/{project}/**"}, "g/Analysis/x/y.bam", skip) is None


def test_ls_names_the_clock_of_its_times(tmp_path, capsys):
    root, (run,) = one_root(tmp_path)
    lib = str(tmp_path / "library")
    run_scan(lib, root, quick=True)
    capsys.readouterr()
    with pytest.raises(SystemExit) as e:
        main(["ls", "-L", lib])
    out = capsys.readouterr().out
    assert e.value.code == 0 and run.run_key in out
    import re
    assert re.search(r"\d{4}-\d{2}-\d{2} \d{2}:\d{2} UTC", out), out
