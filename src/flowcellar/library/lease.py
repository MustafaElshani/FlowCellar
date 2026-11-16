"""One scan at a time per library: ``claims/scan.lease``.

The lease is a small JSON file made with ``O_CREAT | O_EXCL`` that holds host, pid and start time; its mtime is
touched every minute while the scan runs. It is **stale** when

* it was taken on this host and that process is gone, or
* it was taken on this host and the process that now has its pid is another one (the pid was used again, as
  after a reboot): the lease records when its process started and that differs, or the lease has not been
  touched for ``stale_seconds`` (3 minutes) although a scan touches it every minute, or
* it was taken on another host and its mtime did not change across two looks ``stale_seconds`` apart.
  Clocks of two hosts are never compared.

A stale lease is taken over, and the caller is told (it writes the notice ``stale_lease_taken``).
"""
from __future__ import annotations

import json
import os
import socket
import threading
import time
from dataclasses import dataclass

from flowcellar.util.errors import ExitCode, FlowCellarError
from flowcellar.util.timeutil import now_z


class LeaseBusy(FlowCellarError):
    """Another scan holds the library and the caller did not want to wait (exit code 5)."""
    exit_code = ExitCode.PARTIAL


@dataclass(frozen=True)
class Holder:
    host: str
    pid: int
    started: str
    pid_started: float | None = None        # when that process itself started (this host's clock); None: not known

    def words(self) -> str:
        if self.pid <= 0 and not self.started:
            return "a lease file that could not be read"
        return f"process {self.pid} on {self.host} (since {self.started})"


def _age_words(seconds: float) -> str:
    if seconds < 90:
        return f"{int(max(seconds, 0))} s"
    if seconds < 5400:
        return f"{seconds / 60:.0f} min"
    if seconds < 36 * 3600:
        return f"{seconds / 3600:.0f} h"
    return f"{seconds / 86400:.0f} days"


def _pid_started(pid: int) -> float | None:
    """When the process with this pid started, as a time of this host's clock; None where that cannot be told
    (no /proc)."""
    try:
        with open(f"/proc/{pid}/stat", "rb") as fh:
            ticks = int(fh.read().rsplit(b")", 1)[1].split()[19])      # field 22: start time in clock ticks since boot
        with open("/proc/stat", "rb") as fh:
            boot = next(int(line.split()[1]) for line in fh if line.startswith(b"btime"))
        return boot + ticks / os.sysconf("SC_CLK_TCK")
    except (OSError, ValueError, IndexError, StopIteration, AttributeError):
        return None


def _pid_alive(pid: int) -> bool:
    if pid <= 0:
        return False
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    except OSError:
        return False
    return True


def read_holder(path: str) -> Holder | None:
    """Who holds the lease; None when there is no lease or it cannot be read (yet)."""
    try:
        with open(path, "rb") as fh:
            d = json.loads(fh.read() or b"{}")
        born = d.get("pid_started")
        return Holder(str(d["host"]), int(d["pid"]), str(d.get("started", "")),
                      float(born) if isinstance(born, (int, float)) and not isinstance(born, bool) else None)
    except (OSError, ValueError, KeyError, TypeError):
        return None


class Lease:
    """``with Lease(path, wait=True) as lease: ...``; ``lease.taken_over`` names a stale holder that was replaced."""

    def __init__(self, path: str, *, wait: bool = True, poll_seconds: float = 5.0, stale_seconds: float = 180.0,
                 heartbeat_seconds: float = 60.0, on_wait=None, host: str | None = None):
        self.path = path
        self.wait = wait
        self.poll_seconds, self.stale_seconds, self.heartbeat_seconds = poll_seconds, stale_seconds, heartbeat_seconds
        self.on_wait = on_wait
        born = _pid_started(os.getpid())
        self.me = Holder(host or socket.gethostname(), os.getpid(), now_z(), round(born, 2) if born is not None else None)
        self.taken_over: Holder | None = None
        self._stop = threading.Event()
        self._thread: threading.Thread | None = None
        self._touched = 0.0
        self.held = False

    # ------------------------------------------------------------------------------------------ taking
    def _create(self) -> bool:
        os.makedirs(os.path.dirname(self.path), exist_ok=True)
        try:
            fd = os.open(self.path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o664)
        except FileExistsError:
            return False
        with os.fdopen(fd, "w") as fh:
            json.dump({"host": self.me.host, "pid": self.me.pid, "started": self.me.started,
                       "pid_started": self.me.pid_started}, fh)
            fh.flush()
            os.fsync(fh.fileno())
        return True

    def _remove_stale(self, seen: Holder | None) -> None:
        """Move the stale file away, and put it back if it turned out to be somebody's fresh lease."""
        aside = f"{self.path}.stale.{self.me.pid}"
        try:
            os.rename(self.path, aside)
        except OSError:
            return
        now = read_holder(aside)
        if now != seen and now is not None:
            try:
                os.link(aside, self.path)                  # fails when a third process already made a new one
            except OSError:
                pass
        try:
            os.unlink(aside)
        except OSError:
            pass

    def acquire(self, heartbeat: bool = True) -> Lease:
        """Take the lease. With ``heartbeat=False`` the caller keeps it fresh with ``touch_if_due`` until it calls
        ``start_heartbeat`` (the scan does this while it still starts child processes)."""
        watched: tuple[int, float] | None = None           # (mtime of a foreign lease, when we first saw it)
        told = False
        while True:
            if self._create():
                self.held = True
                self._touched = time.monotonic()
                if heartbeat:
                    self.start_heartbeat()
                return self
            holder = read_holder(self.path)
            stale = False
            if holder is None:
                try:                                       # unreadable: stale once it is old enough
                    age = time.time() - os.stat(self.path).st_mtime
                except OSError:
                    continue                               # it went away: try again at once
                stale = age > self.stale_seconds
            elif holder.host == self.me.host:
                stale = not _pid_alive(holder.pid) or self._pid_is_another(holder)
            else:
                try:
                    mtime = os.stat(self.path).st_mtime_ns
                except OSError:
                    continue
                now = time.monotonic()
                if watched is None or watched[0] != mtime:
                    watched = (mtime, now)
                stale = now - watched[1] >= self.stale_seconds
            if stale:
                self.taken_over = holder or Holder("unknown", 0, "")
                self._remove_stale(holder)
                watched = None
                continue
            if not self.wait:
                who = holder.words() if holder else "another process"
                try:
                    age = f", last touched {_age_words(time.time() - os.stat(self.path).st_mtime)} ago"
                except OSError:
                    age = ""
                raise LeaseBusy(f"another scan of this library is running: {who}",
                                hint=f"wait for it, or run scan with --wait. Its lease is {self.path}{age} (a running scan "
                                     "touches it every minute); if no scan is running, delete that file")
            if not told and self.on_wait:
                self.on_wait(holder)
                told = True
            time.sleep(self.poll_seconds)

    def _pid_is_another(self, holder: Holder) -> bool:
        """A lease of this host whose pid is alive: True when that process cannot be the scan that wrote it."""
        try:
            mtime = os.stat(self.path).st_mtime
        except OSError:
            return False
        born = _pid_started(holder.pid)
        if born is not None and holder.pid_started is not None and abs(born - holder.pid_started) > 2.0:
            return True                                    # the pid now belongs to a process born at another moment
        # the heartbeat touches the file every minute; a file left alone for the stale time belongs to nobody
        return time.time() - mtime > max(self.stale_seconds, 3 * self.heartbeat_seconds)

    # ------------------------------------------------------------------------------------------ holding
    def start_heartbeat(self) -> None:
        if self._thread is not None or not self.held:
            return

        def beat() -> None:
            while not self._stop.wait(self.heartbeat_seconds):
                try:
                    os.utime(self.path)
                except OSError:
                    pass
        self._thread = threading.Thread(target=beat, name="flowcellar-lease", daemon=True)
        self._thread.start()

    def touch(self) -> None:
        try:
            os.utime(self.path)
        except OSError:
            pass
        self._touched = time.monotonic()

    def touch_if_due(self) -> None:
        if self.held and time.monotonic() - self._touched >= self.heartbeat_seconds / 2:
            self.touch()

    def release(self) -> None:
        if not self.held:
            return
        self._stop.set()
        if self._thread:
            self._thread.join(2)
        if read_holder(self.path) == self.me:              # never remove a lease that was taken over from us
            try:
                os.unlink(self.path)
            except OSError:
                pass
        self.held = False

    def __enter__(self) -> Lease:
        return self.acquire()

    def __exit__(self, *exc) -> None:
        self.release()
