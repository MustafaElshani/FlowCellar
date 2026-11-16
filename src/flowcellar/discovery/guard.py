"""Checking and walking roots in child processes, so that a root that does not answer cannot hang a scan.

A hard network mount blocks for ever when its server is gone, and a process blocked there cannot even be killed.
So every root is listed and walked in its own child process. The parent only waits on pipes:

* the child's first message is the list of first-level names. No message within ``root_check_seconds`` -> the
  root is **offline**. Fewer than half of the names recorded at the last good scan -> offline as well (an
  automount that answers with an empty folder fails this);
* afterwards the child reports progress. Silence for ``stall_timeout_seconds`` -> the root is **skipped**;
* a child that cannot be stopped is abandoned, not waited for.

The caller treats "offline" and "skipped" alike: nothing it knows under that root changes.
"""
from __future__ import annotations

import contextlib
import multiprocessing
import multiprocessing.connection
import os
import sys
import threading
import time
from dataclasses import dataclass, field

from flowcellar.contracts import InputRule
from flowcellar.discovery.walk import RootWalk, WalkSettings, walk_root

REACHABLE, OFFLINE, SKIPPED = "reachable", "offline", "skipped"


@dataclass
class RootOutcome:
    root_id: str
    path: str
    state: str = REACHABLE                  # reachable, offline, skipped
    reason: str | None = None               # why it is not reachable, in words
    first_level: list[str] = field(default_factory=list)
    walk: RootWalk | None = None            # None unless reachable and a walk was asked for
    check_seconds: float = 0.0
    seconds: float = 0.0


def _child(conn, root_id: str, path: str, rules: tuple, settings: WalkSettings, do_walk: bool) -> None:
    """Runs in the child: list the root, say so, walk it, send the result."""
    try:
        try:
            with os.scandir(path) as it:
                names = sorted(e.name for e in it if not e.name.startswith("."))
        except FileNotFoundError:
            conn.send(("gone", "the folder does not exist"))
            return
        except NotADirectoryError:
            conn.send(("gone", "not a folder"))
            return
        except PermissionError:
            conn.send(("gone", "no permission to list the folder"))
            return
        conn.send(("root", names))
        if do_walk:
            result = walk_root(root_id, path, rules, settings, beat=lambda: conn.send(("beat", None)))
            conn.send(("done", result))
        else:
            conn.send(("done", None))
    except BaseException as e:  # noqa: BLE001  reported to the parent, which decides
        try:
            conn.send(("error", f"{type(e).__name__}: {e}"))
        except Exception:  # noqa: BLE001
            pass
    finally:
        conn.close()


def holds_enough(now: list[str], recorded: list[str] | None) -> bool:
    """At least half of the first-level names recorded for a root are still there (always true for a new root)."""
    if not recorded:
        return True
    have = set(now)
    return 2 * sum(1 for n in recorded if n in have) >= len(recorded)


def _abandon(proc) -> None:
    """Stop a child; if it cannot die (blocked on a dead mount) forget it so that exit does not wait for it."""
    try:
        proc.kill()
        proc.join(0.5)
        if proc.is_alive():
            multiprocessing.process._children.discard(proc)       # type: ignore[attr-defined]
    except Exception:  # noqa: BLE001
        pass


def _start_method() -> str:
    """``fork`` on Linux while this process has one thread (instant, and nothing is imported again); otherwise
    ``spawn``, which is safe with threads and on every platform."""
    if sys.platform == "linux" and _os_threads() == 1 and "fork" in multiprocessing.get_all_start_methods():
        return "fork"
    return "spawn"


def _os_threads() -> int:
    """Threads of this process as the system counts them (numeric libraries start their own, unseen by Python)."""
    try:
        with open("/proc/self/status") as fh:
            for line in fh:
                if line.startswith("Threads:"):
                    return int(line.split()[1])
    except (OSError, ValueError, IndexError):
        pass
    return max(threading.active_count(), 2)


@contextlib.contextmanager
def _importable_main():
    """``spawn`` would run the parent's main script again in every child. The child needs nothing from it (its
    function is found by module name), and running it again fails for code piped into the interpreter and for a
    script without the ``if __name__ == "__main__"`` guard. So the main module is hidden for the moment of starting."""
    main = sys.modules.get("__main__")
    saved = (getattr(main, "__file__", None), getattr(main, "__spec__", None)) if main is not None else None
    try:
        if main is not None:
            main.__file__, main.__spec__ = None, None
        yield
    finally:
        if main is not None:
            main.__file__, main.__spec__ = saved


def check_roots(roots: list[dict], *, rules: tuple[InputRule, ...] = (), settings: dict[str, WalkSettings] | None = None,
                recorded: dict[str, list[str]] | None = None, check_seconds: float = 20.0,
                stall_seconds: float = 120.0, walk: bool = True, start_method: str | None = None,
                target=None, tick=None) -> dict[str, RootOutcome]:
    """Check every root (``{"id", "path"}``) and, when ``walk`` is true, walk the ones that answer.

    All roots are worked on at the same time, one child each; the call returns when every root has finished,
    failed its check or stalled. ``recorded`` maps a root id to the first-level names of its last good scan.
    ``target`` replaces the child's function (tests use it to act out a mount that hangs). ``tick`` is called
    about five times a second while waiting (the scan keeps its lease fresh with it).
    """
    ctx = multiprocessing.get_context(start_method or _start_method())
    recorded = recorded or {}
    settings = settings or {}
    out: dict[str, RootOutcome] = {}
    live: dict = {}                                   # connection -> [root id, process, checked, last word, started]
    t0 = time.monotonic()
    for r in roots:
        out[r["id"]] = RootOutcome(r["id"], r["path"])
        parent, child = ctx.Pipe(duplex=False)
        proc = ctx.Process(target=target or _child, name=f"flowcellar-root-{r['id']}", daemon=True,
                           args=(child, r["id"], r["path"], tuple(rules), settings.get(r["id"], WalkSettings()), walk))
        with _importable_main():
            proc.start()
        child.close()
        live[parent] = [r["id"], proc, False, time.monotonic(), time.monotonic()]

    def finish(conn, state: str | None = None, reason: str | None = None) -> None:
        rid, proc, _, _, started = live.pop(conn)
        o = out[rid]
        if state:
            o.state, o.reason, o.walk = state, reason, None
        o.seconds = time.monotonic() - started
        try:
            conn.close()
        except OSError:
            pass
        if state:
            _abandon(proc)
        else:
            proc.join(5)
            if proc.is_alive():
                _abandon(proc)

    while live:
        ready = multiprocessing.connection.wait(list(live), timeout=0.2)
        now = time.monotonic()
        if tick:
            tick()
        for conn in ready:
            rid = live[conn][0]
            try:
                kind, payload = conn.recv()
            except (EOFError, OSError):
                finish(conn, OFFLINE if not live[conn][2] else SKIPPED, "the walk stopped unexpectedly")
                continue
            live[conn][3] = now
            if kind == "root":
                live[conn][2] = True
                o = out[rid]
                o.first_level, o.check_seconds = payload, now - live[conn][4]
                if not holds_enough(payload, recorded.get(rid)):
                    had = len(recorded[rid])
                    finish(conn, OFFLINE, f"it holds {sum(1 for n in recorded[rid] if n in set(payload))} of the {had} "
                                          "first-level folders it had at the last scan (an empty mount point?)")
            elif kind == "done":
                out[rid].walk = payload
                finish(conn)
            elif kind == "gone":
                finish(conn, OFFLINE, payload)
            elif kind == "error":
                finish(conn, OFFLINE if not live[conn][2] else SKIPPED, payload)
        for conn in list(live):
            rid, _, checked, last, started = live[conn]
            if not checked and now - started > check_seconds:
                finish(conn, OFFLINE, f"it did not answer within {check_seconds:g} s")
            elif checked and now - last > stall_seconds:
                finish(conn, SKIPPED, f"the walk made no progress for {stall_seconds:g} s")
    for o in out.values():
        o.seconds = o.seconds or (time.monotonic() - t0)
    return out
