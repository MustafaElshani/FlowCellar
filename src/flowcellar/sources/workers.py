"""Worker processes for reading large inputs.

All input I/O of a pass happens in child processes: a read that never returns (a hard network mount whose server
died) then costs one abandoned child and never the program. The parent only waits, with a clock: a child that
reports no progress for ``stall_timeout_s`` is killed where possible and abandoned otherwise.

A task function has the signature ``fn(args: dict, progress: Callable[[int], None]) -> picklable`` and lives at
module level; it is named by ``"package.module:function"`` so that every start method can find it.
"""
from __future__ import annotations

import importlib
import multiprocessing as mp
import os
import threading
import time
import traceback
from dataclasses import dataclass
from multiprocessing import spawn as _spawn
from multiprocessing.connection import wait as _wait
from typing import Callable

_PRELOAD = ["flowcellar.sources.seqsummary", "flowcellar.sources.porescan"]
_ctx = None


def context():
    """The start method: a fork server where the platform has one (children start in milliseconds from a clean,
    single-threaded process that has numpy and pyarrow loaded), else spawn. FLOWCELLAR_MP_START overrides it."""
    global _ctx
    if _ctx is None:
        want = os.environ.get("FLOWCELLAR_MP_START", "")
        methods = mp.get_all_start_methods()
        method = want if want in methods else "forkserver" if "forkserver" in methods else "spawn"
        _ctx = mp.get_context(method)
        if method == "forkserver":
            try:
                _ctx.set_forkserver_preload(_PRELOAD)
            except Exception:  # noqa: BLE001  preloading is an optimisation only
                pass
    return _ctx


_start_lock = threading.Lock()


def _start(process) -> None:
    """Start a child that does not import the caller's ``__main__``.

    The spawn and fork-server start methods re-import the main module in every child, which runs a script again
    when it has no ``if __name__ == "__main__"`` guard and fails outright for a program read from standard input.
    A worker here needs nothing from the main module (its task function is named by module path), so that step is
    left out of the data the child is prepared with.
    """
    with _start_lock:
        original = _spawn.get_preparation_data

        def without_main(name):
            d = original(name)
            d.pop("init_main_from_name", None)
            d.pop("init_main_from_path", None)
            return d
        _spawn.get_preparation_data = without_main
        try:
            process.start()
        finally:
            _spawn.get_preparation_data = original


@dataclass
class TaskResult:
    state: str                  # ok, error, stalled, died
    value: object = None
    error: str | None = None
    seconds: float = 0.0


def _child(fn_name: str, args: dict, conn) -> None:
    try:
        try:
            os.nice(int(args.get("_nice", 0)))
        except (OSError, AttributeError, ValueError):
            pass
        mod, _, name = fn_name.partition(":")
        fn = getattr(importlib.import_module(mod), name)
        last = [0.0]

        def progress(n: int) -> None:
            now = time.monotonic()
            if now - last[0] >= 0.25:
                last[0] = now
                conn.send(("p", int(n)))
        conn.send(("p", 0))
        value = fn(args, progress)
        conn.send(("ok", value))
    except BaseException as e:  # noqa: BLE001  everything is reported to the parent
        try:
            conn.send(("error", f"{type(e).__name__}: {e}\n{traceback.format_exc(limit=6)}"))
        except Exception:  # noqa: BLE001
            pass
    finally:
        try:
            conn.close()
        except Exception:  # noqa: BLE001
            pass


def run_tasks(fn_name: str, tasks: list[dict], *, stall_timeout_s: float = 120.0,
              on_progress: Callable[[int], None] | None = None) -> list[TaskResult]:
    """Run every task in its own child process, all at once, and wait for them. Never raises for a task: the
    outcome of each is in its TaskResult."""
    ctx = context()
    n = len(tasks)
    results: list[TaskResult | None] = [None] * n
    procs, conns = [], []
    t_start = time.monotonic()
    for t in tasks:
        parent, child = ctx.Pipe(duplex=False)
        p = ctx.Process(target=_child, args=(fn_name, t, child), daemon=True)
        _start(p)
        child.close()
        procs.append(p)
        conns.append(parent)
    last_seen = [time.monotonic()] * n
    done_bytes = [0] * n
    pending = set(range(n))
    while pending:
        ready = _wait([conns[i] for i in pending], timeout=1.0)
        now = time.monotonic()
        for conn in ready:
            i = conns.index(conn)
            try:
                kind, payload = conn.recv()
            except (EOFError, OSError):
                code = procs[i].exitcode
                results[i] = TaskResult("died", error=f"the worker process ended without a result (exit code {code})", seconds=now - t_start)
                pending.discard(i)
                continue
            last_seen[i] = now
            if kind == "p":
                done_bytes[i] = payload
            elif kind == "ok":
                results[i] = TaskResult("ok", value=payload, seconds=now - t_start)
                pending.discard(i)
            else:
                results[i] = TaskResult("error", error=str(payload), seconds=now - t_start)
                pending.discard(i)
        if on_progress and ready:
            try:
                on_progress(sum(done_bytes))
            except Exception:  # noqa: BLE001  a display problem never stops a pass
                pass
        for i in list(pending):
            if now - last_seen[i] > stall_timeout_s:
                results[i] = TaskResult("stalled", error=f"no progress for {stall_timeout_s:.0f} s", seconds=now - t_start)
                pending.discard(i)
                try:
                    procs[i].kill()
                except Exception:  # noqa: BLE001
                    pass
    for i, p in enumerate(procs):
        p.join(timeout=0.2 if results[i] and results[i].state == "stalled" else 10)
        if p.is_alive():
            try:
                p.kill()
            except Exception:  # noqa: BLE001
                pass
        try:
            conns[i].close()
        except Exception:  # noqa: BLE001
            pass
    return [r if r is not None else TaskResult("died", error="no result") for r in results]


def run_one(fn_name: str, args: dict, *, stall_timeout_s: float = 120.0) -> TaskResult:
    return run_tasks(fn_name, [args], stall_timeout_s=stall_timeout_s)[0]


# ------------------------------------------------------------------------------------- reading by pieces
class Pieces:
    """Blocks of whole lines from one byte range of a text table.

    A range owns every line that starts in it: the line that straddles ``lo`` belongs to the range before, the line
    that straddles ``hi`` is finished here. After the loop: ``done`` (bytes of whole lines handed out), ``stopped``
    (the deadline was reached), ``cut_tail`` (the file ends inside a line; that line was not handed out).
    """

    def __init__(self, path: str, lo: int, hi: int, body_start: int, *, piece: int = 32 << 20, deadline: float | None = None):
        self.path, self.lo, self.hi, self.body_start = path, lo, hi, body_start
        self.piece, self.deadline = piece, deadline
        self.done = 0
        self.first = lo
        self.stopped = False
        self.cut_tail = False
        self.cut_bytes = 0

    def __iter__(self):
        with open(self.path, "rb", buffering=0) as fh:
            if self.lo > self.body_start:
                fh.seek(self.lo - 1)
                if fh.read(1) != b"\n":
                    fh.readline()                 # the line that straddles lo belongs to the range before
            else:
                fh.seek(self.body_start)
            pos = self.first = fh.tell()
            carry = b""
            at_eof = False
            while pos < self.hi:
                if self.deadline is not None and time.monotonic() > self.deadline:
                    self.stopped = True
                    break
                buf = fh.read(min(self.piece, self.hi - pos))
                if not buf:
                    at_eof = True
                    break
                pos += len(buf)
                if pos >= self.hi and not buf.endswith(b"\n"):
                    tail = fh.readline()          # finish the line that straddles hi
                    buf += tail
                    if not tail.endswith(b"\n"):
                        at_eof = True
                data = carry + buf
                cut = data.rfind(b"\n") + 1
                carry = data[cut:]
                data = data[:cut]
                if data:
                    self.done += len(data)
                    yield data
            if carry and not self.stopped and (at_eof or pos >= self.hi):
                self.cut_tail = True
                self.cut_bytes = len(carry)


def split_ranges(body_start: int, size: int, n: int) -> list[tuple[int, int]]:
    body = max(size - body_start, 0)
    n = max(1, n)
    step = body / n
    bounds = [body_start + int(round(i * step)) for i in range(n)] + [size]
    return [(bounds[i], bounds[i + 1]) for i in range(n)]


def _sleep_task(args: dict, progress) -> None:
    """A task that makes no progress: used to test that a stalled worker is abandoned."""
    time.sleep(float(args.get("seconds", 1)))
