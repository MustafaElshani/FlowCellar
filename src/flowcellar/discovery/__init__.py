"""Discovery: which run folders and read-data files lie under the roots of a library.

    found = discover(config, recorded=first_level_names_of_the_last_scan)
    found.roots["minknow"].state          "reachable", "offline" or "skipped"
    found.runs                            [(root id, FoundFolder), ...]   run folders and raw-only folders
    found.candidates                      [Candidate, ...]                read-data files with their role

Nothing is opened except small tables that are told by their content (through ``readdata.api.sniff``), and
nothing is written. A root that does not answer is reported and left alone.
"""
from __future__ import annotations

import dataclasses
import os
import time
from dataclasses import dataclass, field

from flowcellar.contracts import Candidate, InputRule
from flowcellar.discovery.guard import OFFLINE, REACHABLE, SKIPPED, RootOutcome, check_roots
from flowcellar.discovery.patterns import rules_from_config
from flowcellar.discovery.walk import FoundFolder, RootWalk, WalkSettings, experiment_dir, walk_root
from flowcellar.readdata import api as readdata_api

__all__ = ["Discovery", "discover", "input_rules", "walk_root", "check_roots", "RootOutcome", "RootWalk", "FoundFolder",
           "WalkSettings", "REACHABLE", "OFFLINE", "SKIPPED"]

SNIFF_MAX_BYTES = 1 << 20


@dataclass
class Discovery:
    roots: dict[str, RootOutcome]
    runs: list[tuple[str, FoundFolder]] = field(default_factory=list)
    lookalikes: list[tuple[str, FoundFolder]] = field(default_factory=list)
    candidates: list[Candidate] = field(default_factory=list)
    analyses: list[Candidate] = field(default_factory=list)                # anchor files of single-cell and spatial analyses
    problems: list[str] = field(default_factory=list)
    over_budget: list[tuple[str, str]] = field(default_factory=list)       # (root id, subtree path)
    unreadable: list[tuple[str, str]] = field(default_factory=list)        # (root id, folder that could not be listed)
    seconds: float = 0.0

    def counts(self) -> dict[str, int]:
        c = {"listings": 0, "entries": 0, "stats": 0, "pruned": 0, "anchors": 0}
        for o in self.roots.values():
            if o.walk:
                for k in ("listings", "entries", "stats", "pruned"):
                    c[k] += o.walk.counts.get(k, 0)
                c["anchors"] += len(o.walk.anchors)
        c["runs"] = sum(1 for _, f in self.runs if f.kind == "run")
        c["raw_only"] = sum(1 for _, f in self.runs if f.kind == "run_raw_only")
        c["lookalikes"] = len(self.lookalikes)
        c["files"] = len(self.candidates)
        c["analyses"] = len(self.analyses)
        return c


def input_rules(config: dict) -> list[InputRule]:
    """The configured ``[[inputs]]`` rules first, then the generic ones the package ships: those of read data,
    and the two anchor files by which analyses of the single-cell and spatial add-on are found (names only:
    nothing of the add-on but its list of patterns is loaded here)."""
    from flowcellar.sc.detect import input_rules as analysis_rules
    return rules_from_config(config, readdata_api.default_input_rules()) + analysis_rules()


def _nested(roots: list[dict]) -> dict[str, tuple[str, ...]]:
    """For each root, the other roots that lie inside it (each folder is walked under one root only)."""
    out = {}
    for r in roots:
        inner = [o["path"] for o in roots if o["id"] != r["id"] and o["path"].startswith(r["path"].rstrip("/") + "/")]
        out[r["id"]] = tuple(inner)
    return out


def discover(config: dict, *, roots: list[dict] | None = None, recorded: dict[str, list[str]] | None = None,
             sniff_cache: dict | None = None, walk: bool = True, start_method: str | None = None,
             tick=None) -> Discovery:
    """Check and walk the roots of ``config`` (or the given ones). ``sniff_cache`` (path -> [size, mtime_ns, role])
    is read and updated, so a table that did not change is not opened again."""
    t0 = time.monotonic()
    roots = list(config["roots"] if roots is None else roots)
    rules = tuple(input_rules(config))
    nested = _nested(roots)
    settings = {r["id"]: WalkSettings.from_config(config, nested[r["id"]]) for r in roots}
    d = config["discovery"]
    outcomes = check_roots(roots, rules=rules, settings=settings, recorded=recorded,
                           check_seconds=float(d["root_check_seconds"]), stall_seconds=float(d["stall_timeout_seconds"]),
                           walk=walk, start_method=start_method, tick=tick)
    found = Discovery(roots=outcomes)
    seen_dirs: set[str] = set()
    for r in roots:                                   # configuration order: the first root that holds a path wins
        o = outcomes[r["id"]]
        if o.state != REACHABLE or o.walk is None:
            continue
        w = o.walk
        found.problems.extend(p for p in w.problems if p not in found.problems)
        found.over_budget.extend((r["id"], p) for p in w.over_budget)
        found.unreadable.extend((r["id"], p) for p in getattr(w, "unreadable", ()))
        for f in w.runs:
            if f.path not in seen_dirs:
                seen_dirs.add(f.path)
                found.runs.append((r["id"], f))
        found.lookalikes.extend((r["id"], f) for f in w.lookalikes)
        holders = set(w.holders)
        sheets = _sniff(w, r, sniff_cache, float(d["min_age_seconds"]))
        holders.update(os.path.dirname(c.path) for c in sheets)
        for c in list(w.candidates) + sheets:
            c = dataclasses.replace(c, experiment_dir=experiment_dir(c.path, w.path, holders))
            # an analysis is not a read-data file of a run: it is kept apart and never linked as one
            (found.analyses if c.role.startswith("analysis.") else found.candidates).append(c)
    found.seconds = time.monotonic() - t0
    return found


def _sniff(w: RootWalk, root: dict, cache: dict | None, min_age: float) -> list[Candidate]:
    """Tables with free names in the folders that hold a run folder or an anchor, told by their content."""
    out = []
    now = time.time()
    known = {c.path for c in w.candidates}
    for path in w.sniffable:
        if path in known:
            continue
        try:
            st = os.lstat(path)
        except OSError:
            continue
        if st.st_size > SNIFF_MAX_BYTES or st.st_size == 0 or now - st.st_mtime < min_age:
            continue
        hit = cache.get(path) if cache is not None else None
        if hit and hit[0] == st.st_size and hit[1] == st.st_mtime_ns:
            role = hit[2]
        else:
            try:
                role = readdata_api.sniff(path)
            except OSError:
                role = None
            if cache is not None:
                cache[path] = [st.st_size, st.st_mtime_ns, role]
        if role:
            out.append(Candidate(role=role, path=path, size=st.st_size, mtime_ns=st.st_mtime_ns, root_id=root["id"],
                                 relpath=os.path.relpath(path, w.path), rule="told by its content"))
    return out
