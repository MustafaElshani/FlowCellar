"""The two large passes of tier 1 and their cache.

``full_pass`` reads the sequencing summary and the pore-scan table of one run, each once, and stores each result
as a partial keyed by the fingerprint of its input. A second call finds the partials and reads nothing. Every
byte of input is read in worker processes (``flowcellar.sources.workers``); this process only waits.
"""
from __future__ import annotations

import os
import time

from flowcellar.contracts import PartialRef, PassOptions, ProgressEvent, RunIdentity
from flowcellar.model import partial as P
from flowcellar.sources import minknow as mk
from flowcellar.sources import porescan, seqsummary, workers
from flowcellar.util import atomic
from flowcellar.util import log as _log
from flowcellar.util.fingerprint import hash_obj

#: errors of a pass that reading the same bytes again would repeat: such a partial is final
_FINAL_ERRORS = {"columns_missing", "too_many_bad_rows"}
_DEFAULT_CHANNELS = {"promethion": 3000, "minion": 512, "flongle": 126}


# ---------------------------------------------------------------------------- fingerprints without reading
def _stat_cache_path(partials_dir: str, path: str) -> str:
    return os.path.join(partials_dir, "stat", hash_obj({"path": os.path.abspath(path)}, 32) + ".json")


def cached_fingerprint(partials_dir: str, path: str) -> dict | None:
    """The fingerprint taken at an earlier pass, when size and modification time are still the same."""
    try:
        d = atomic.read_json(_stat_cache_path(partials_dir, path), retries=1)
        st = os.stat(path)
    except (OSError, ValueError):
        return None
    if isinstance(d, dict) and d.get("size") == st.st_size and d.get("mtime_ns") == st.st_mtime_ns and isinstance(d.get("fingerprint"), dict):
        return d["fingerprint"]
    return None


def _remember_fingerprint(partials_dir: str, path: str, probe: dict) -> None:
    try:
        atomic.atomic_write_json(_stat_cache_path(partials_dir, path),
                                 {"path": os.path.abspath(path), "size": probe["size"], "mtime_ns": probe["mtime_ns"],
                                  "fingerprint": probe["fingerprint"]})
    except OSError:
        pass


# --------------------------------------------------------------------------------------- the prepare task
def prep_task(args: dict, progress) -> dict:
    """Runs in a child: what the passes need to know before they start. Probes each large input (size, time,
    fingerprint, header) and reads the number of bases per event and the channel count from the report, or from
    the throughput file when the report cannot be read."""
    out: dict = {"probes": {}, "errors": {}}
    for role, path in args["inputs"].items():
        try:
            out["probes"][role] = (seqsummary if role == seqsummary.ROLE else porescan).probe_task({"path": path}, progress)
        except OSError as e:
            out["errors"][role] = f"{e.strerror or e}"
        progress(0)
    ratio = channels = None
    rj = args.get("report")
    if rj:
        rep, _problem = mk.read_report(rj)
        seq, _n = mk.sequencing_acquisition(rep)
        cfg = mk.g(seq or {}, "acquisition_run_info", "config_summary", default={}) or {}
        ratio = mk.to_float(cfg.get("events_to_base_ratio"))
        fc = mk.g(rep, "protocol_run_info", "flow_cell", default={}) or {}
        dev = mk.g(rep, "protocol_run_info", "device", default={}) or {}
        channels = mk.to_int(fc.get("channel_count"), 0) or mk.to_int(dev.get("max_channel_count"), 0) or mk.to_int(cfg.get("channel_count"), 0) or None
        del rep
    progress(0)
    source = "report"
    if not ratio or ratio <= 0:
        tp = mk.read_throughput(args.get("throughput"))
        from flowcellar.modules.yield_ import ratio_of
        ratio, source = ratio_of({}, tp)
    out.update(ratio=ratio, ratio_source=source, channels=channels)
    return out


def pass_parameters(run: RunIdentity, prep: dict) -> tuple[int, float]:
    channels = prep.get("channels") or run.facts.get("channel_count") or _DEFAULT_CHANNELS.get(run.facts.get("family") or "", 3000)
    return int(channels), float(prep.get("ratio") or 0.6)


# ------------------------------------------------------------------------------------------------ the pass
def _ref(p: dict, path: str, *, cached: bool, seconds: float = 0.0) -> PartialRef:
    return PartialRef(partial_id=p["partial_id"], path=path, role=p["role"], input_path=p["input"]["path"], state=p["state"],
                      complete=P.is_complete(p), cached=cached, seconds=seconds)


def _unopened(role: str, path: str) -> PartialRef:
    """The answer for an input that could not be opened: an error with no partial behind it (``partial_id`` is
    empty). The caller must not take the run for counted; nothing is cached, so the next call tries again."""
    return PartialRef(partial_id="", path="", role=role, input_path=path, state="error", complete=False, cached=False)


def _is_final(p: dict) -> bool:
    """A cached partial that a new pass over the same bytes could not improve."""
    if P.is_complete(p) or p.get("state") == "truncated":
        return True
    return p.get("state") == "error" and bool(p.get("errors")) and all(e.get("code") in _FINAL_ERRORS for e in p["errors"])


def _find(partials_dir: str, fp: dict, role: str, modules: dict) -> tuple[str, dict | None, str | None]:
    pid = P.partial_id(fp, role, modules)
    path = P.find_partial(partials_dir, pid)
    if not path:
        return pid, None, None
    try:
        return pid, P.read_partial(path), path
    except Exception:  # noqa: BLE001  an unreadable cache file is no cache
        return pid, None, None


def full_pass(run: RunIdentity, partials_dir: str, opts: PassOptions = PassOptions()) -> list[PartialRef]:
    if run.kind != "run" or run.run_state == "running":
        return []                       # nothing large to read, or files that still change: never cached
    plan = [(seqsummary.ROLE, seqsummary.MODULES), (porescan.ROLE, porescan.MODULES)]
    plan = [(role, mods) for role, mods in plan if run.files.get(role)]
    if not plan:
        return []
    refs: dict[str, PartialRef] = {}
    earlier: dict[str, PartialRef] = {}          # a time-limited result of an earlier call: better than nothing
    todo: list[tuple[str, dict]] = []
    for role, mods in plan:
        fp = cached_fingerprint(partials_dir, run.files[role])
        if fp:
            _pid, p, path = _find(partials_dir, fp, role, mods)
            if p is not None and _is_final(p):
                refs[role] = _ref(p, path, cached=True)
                continue
            if p is not None:
                earlier[role] = _ref(p, path, cached=True)
        todo.append((role, mods))

    def emit(what, path=None, fraction=None, mb_per_s=None, message=None):
        if opts.progress:
            try:
                opts.progress(ProgressEvent(what=what, path=path, fraction=fraction, mb_per_s=mb_per_s, message=message))
            except Exception:  # noqa: BLE001
                pass

    def out_of_time() -> bool:
        return opts.deadline is not None and time.monotonic() >= opts.deadline

    if todo and not out_of_time():
        prep_r = workers.run_one(f"{__name__}:prep_task",
                                 {"inputs": {role: run.files[role] for role, _ in todo}, "report": run.files.get("minknow.report_json"),
                                  "throughput": run.files.get("minknow.throughput")}, stall_timeout_s=opts.stall_timeout_s)
        if prep_r.state != "ok":
            _log.event("pass.error", f"{run.run_key}: the input files could not be opened ({(prep_r.error or prep_r.state).splitlines()[0]})",
                       level=30, run_key=run.run_key)
            emit("error", message=(prep_r.error or prep_r.state).splitlines()[0])
            for role, _mods in todo:                 # nothing was read: say so for every file, never "nothing to do"
                if role not in earlier:
                    refs[role] = _unopened(role, run.files[role])
            todo = []
            prep = {}
        else:
            prep = prep_r.value
        channels, ratio = pass_parameters(run, prep) if prep else (0, 0.0)
        for role, mods in todo:
            path = run.files[role]
            probe = prep["probes"].get(role)
            if probe is None:
                emit("error", path=path, message=prep["errors"].get(role, "could not be opened"))
                if role not in earlier:
                    refs[role] = _unopened(role, path)
                continue
            _remember_fingerprint(partials_dir, path, probe)
            _pid, old, old_path = _find(partials_dir, probe["fingerprint"], role, mods)
            if old is not None and _is_final(old) and (old.get("facts") or {}).get("channel_count") in (None, channels) \
                    and (role != seqsummary.ROLE or (old.get("facts") or {}).get("events_to_base_ratio") == ratio):
                refs[role] = _ref(old, old_path, cached=True)
                continue
            if out_of_time():
                if old is not None:
                    refs[role] = _ref(old, old_path, cached=True)
                emit("skipped", path=path, message="time limit reached before this file")
                continue
            t0 = time.monotonic()
            body = max(probe["size"] - probe["body_start"], 1)

            def on_progress(done, _t0=t0, _body=body, _path=path):
                dt = time.monotonic() - _t0
                emit("reading", path=_path, fraction=min(done / _body, 1.0), mb_per_s=(done / 1e6 / dt) if dt > 0 else None)
            if role == seqsummary.ROLE:
                res = seqsummary.run_pass(path, probe, channels=channels, ratio=ratio, jobs=opts.jobs, gentle=opts.gentle,
                                          deadline=opts.deadline, stall_timeout_s=opts.stall_timeout_s, on_progress=on_progress)
                part = seqsummary.make_partial(path, probe, res, run.run_key, channels=channels, ratio=ratio)
            else:
                res = porescan.run_pass(path, probe, channels=channels, deadline=opts.deadline, stall_timeout_s=opts.stall_timeout_s,
                                        on_progress=on_progress, gentle=opts.gentle)
                part = porescan.make_partial(path, probe, res, run.run_key, channels=channels)
            secs = time.monotonic() - t0
            if old is not None and not P.is_complete(part) and part["state"] != "truncated" and \
                    ((old.get("coverage") or {}).get("fraction") or 0) >= ((part.get("coverage") or {}).get("fraction") or 0) \
                    and old.get("state") != "error":
                refs[role] = _ref(old, old_path, cached=True, seconds=secs)       # the earlier, larger part stays
            else:
                written = P.write_partial(partials_dir, part)
                refs[role] = _ref(part, written, cached=False, seconds=secs)
            st = refs[role].state
            emit("done" if st in ("ok", "truncated", "partially_read") else "error", path=path, fraction=(part["coverage"] or {}).get("fraction"),
                 mb_per_s=(part.get("timings") or {}).get("mb_per_s"),
                 message=None if st == "ok" else "; ".join(e.get("message", "") for e in part.get("errors", [])) or st)
    for role, ref in earlier.items():
        refs.setdefault(role, ref)
    return [refs[role] for role, _ in plan if role in refs]
