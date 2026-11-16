"""Worker processes for the read pass.

All input I/O of READ DATA happens in worker processes: a file on a network mount that stops answering
blocks a worker, never the parent. Workers advance a shared counter as they read; a pool whose counter stands still
and that returns nothing for ``stall_timeout_s`` is terminated and the caller records the file as stalled.

Workers are started from a fork server (or spawned where the platform has none), never forked from the caller: the
caller may hold threads and open handles.
"""
from __future__ import annotations

import multiprocessing as mp
import time

_BEAT = None            # set in each worker: the shared progress counter


def _init(counter) -> None:
    global _BEAT
    _BEAT = counter


def beat(n: int = 1) -> None:
    """Called by a worker while it reads. Not exact (no lock): it only has to change while work goes on."""
    c = _BEAT
    if c is not None:
        c.value += n


def context():
    if "forkserver" in mp.get_all_start_methods():
        ctx = mp.get_context("forkserver")
        try:
            ctx.set_forkserver_preload(["flowcellar.readdata.bam", "flowcellar.readdata.tables"])
        except Exception:       # noqa: BLE001  (a server that is already running keeps its own list)
            pass
        return ctx
    return mp.get_context("spawn")


def run_tasks(fn, args: list, *, workers: int = 1, stall_timeout_s: float = 120.0, inprocess: bool = False,
              progress=None) -> dict:
    """``fn(a)`` for every ``a`` in ``args``, in worker processes. Returns ``{"results": [...], "stalled": bool}``;
    results come in completion order. An exception raised by ``fn`` is raised here.

    ``progress(done, total)`` is called after every finished task. ``inprocess=True`` runs everything in the calling
    process (tests, debugging).
    """
    if not args:
        return {"results": [], "stalled": False}
    if inprocess:
        out = []
        for i, a in enumerate(args):
            out.append(fn(a))
            if progress:
                progress(i + 1, len(args))
        return {"results": out, "stalled": False}
    ctx = context()
    counter = ctx.RawValue("q", 0)
    results: list = []
    stalled = False
    pool = ctx.Pool(processes=max(1, min(int(workers), len(args))), initializer=_init, initargs=(counter,))
    try:
        it = pool.imap_unordered(fn, args)
        seen, changed = 0, time.monotonic()
        poll = max(0.05, min(1.0, stall_timeout_s / 4))
        while True:
            try:
                res = it.next(timeout=poll)
            except StopIteration:
                break
            except mp.TimeoutError:
                v = counter.value
                if v != seen:
                    seen, changed = v, time.monotonic()
                elif time.monotonic() - changed > stall_timeout_s:
                    stalled = True
                    break
                continue
            results.append(res)
            changed = time.monotonic()
            if progress:
                progress(len(results), len(args))
    finally:
        pool.terminate()
        pool.join()
    return {"results": results, "stalled": stalled}
