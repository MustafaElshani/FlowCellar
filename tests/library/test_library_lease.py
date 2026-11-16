"""The scan lease: one scan at a time per library."""
from __future__ import annotations

import json
import os
import socket
import threading
import time

import pytest

from flowcellar.library.lease import Lease, LeaseBusy, read_holder


def test_take_hold_release(tmp_path):
    path = str(tmp_path / "claims" / "scan.lease")
    with Lease(path, heartbeat_seconds=0.05) as lease:
        held = json.load(open(path))
        assert held == {"host": socket.gethostname(), "pid": os.getpid(), "started": lease.me.started,
                        "pid_started": lease.me.pid_started}
        first = os.stat(path).st_mtime_ns
        time.sleep(0.3)
        assert os.stat(path).st_mtime_ns > first, "the time stamp is kept fresh while the lease is held"
        with pytest.raises(LeaseBusy) as e:
            Lease(path, wait=False).acquire()
        assert e.value.exit_code == 5 and socket.gethostname() in e.value.message
    assert not os.path.exists(path) and lease.taken_over is None


def test_a_second_scan_waits_for_the_first(tmp_path):
    path = str(tmp_path / "scan.lease")
    first = Lease(path).acquire()
    waited = []
    got = []

    def second():
        with Lease(path, wait=True, poll_seconds=0.05, on_wait=waited.append) as lease:
            got.append(lease.taken_over)
    t = threading.Thread(target=second)
    t.start()
    time.sleep(0.3)
    assert got == [] and len(waited) == 1 and waited[0].pid == os.getpid()
    first.release()
    t.join(5)
    assert got == [None] and not os.path.exists(path)


def test_stale_leases(tmp_path):
    path = str(tmp_path / "scan.lease")
    json.dump({"host": socket.gethostname(), "pid": 2 ** 22 + 4242, "started": "2026-01-01T00:00:00Z"}, open(path, "w"))
    with Lease(path, wait=False) as lease:                            # this host, process gone: stale at once
        assert lease.taken_over.pid == 2 ** 22 + 4242 and read_holder(path) == lease.me
    json.dump({"host": "elsewhere", "pid": 7, "started": "2026-01-01T00:00:00Z"}, open(path, "w"))
    t0 = time.monotonic()
    with Lease(path, wait=True, poll_seconds=0.05, stale_seconds=0.4) as lease:   # another host, time stamp stands still
        assert lease.taken_over.host == "elsewhere" and time.monotonic() - t0 >= 0.4
    # another host that is alive keeps touching its lease: never stale
    json.dump({"host": "elsewhere", "pid": 7, "started": "2026-01-01T00:00:00Z"}, open(path, "w"))
    stop = threading.Event()

    def alive():
        while not stop.wait(0.05):
            os.utime(path)
    t = threading.Thread(target=alive)
    t.start()
    result = []

    def waiter():
        try:
            lease = Lease(path, wait=True, poll_seconds=0.05, stale_seconds=0.4).acquire()
            result.append(lease.taken_over)
            lease.release()
        except Exception as e:  # noqa: BLE001
            result.append(e)
    w = threading.Thread(target=waiter, daemon=True)
    w.start()
    time.sleep(1.2)
    assert result == [], "a lease that is kept fresh is not taken over"
    stop.set()
    t.join(2)
    w.join(5)
    assert len(result) == 1 and result[0].host == "elsewhere"


def test_release_never_removes_somebody_elses_lease(tmp_path):
    path = str(tmp_path / "scan.lease")
    lease = Lease(path).acquire()
    json.dump({"host": "elsewhere", "pid": 7, "started": "2026-01-02T00:00:00Z"}, open(path, "w"))     # taken over from us
    lease.release()
    assert read_holder(path).host == "elsewhere"
