"""The discovery walk over one root. Listings and ``stat`` calls only: no file is opened.

Part 1 is a generic pruned walk. It recognises run folders by the files in them (``model.runfolder``), never by
their name, does not enter them, and records the *anchors* of the input patterns without entering those either.
Part 2 resolves every input pattern from its anchors by directed descent (``discovery.patterns``).

Rules: symlinks are never followed; hidden names are ignored; a read-data file changed within
``min_age_seconds`` is left for the next scan; part 1 stops at 12 levels and gives up a first-level subtree that
needs more than ``budget`` listings, naming it.
"""
from __future__ import annotations

import fnmatch
import os
import stat as _stat
import time
from dataclasses import dataclass, field
from typing import Callable

from flowcellar.contracts import Candidate, InputRule
from flowcellar.discovery.patterns import Rule, compile_rules, descend
from flowcellar.model import runfolder

MAX_DEPTH = 12
_ROLE_PATTERNS = tuple(runfolder.TIER1_ROLES.values())
SNIFF_SUFFIXES = (".csv",)
SNIFF_MAX_PER_FOLDER = 50


@dataclass(frozen=True)
class FoundFolder:
    """A folder the walk recognised: ``kind`` is run, run_raw_only or lookalike."""
    kind: str
    path: str
    relpath: str                     # under its root; "" when the root itself is the folder
    evidence: tuple[str, ...] = ()   # the names in its listing that decided
    files: tuple[tuple[str, int, int], ...] = ()   # run folders: (name, size, mtime in ns) of the MinKNOW files


@dataclass
class RootWalk:
    """What one root holds. Plain data: it crosses a process boundary."""
    root_id: str
    path: str
    first_level: list[str] = field(default_factory=list)
    runs: list[FoundFolder] = field(default_factory=list)
    lookalikes: list[FoundFolder] = field(default_factory=list)
    anchors: list[str] = field(default_factory=list)
    candidates: list[Candidate] = field(default_factory=list)
    holders: list[str] = field(default_factory=list)        # folders that hold a run folder or an anchor
    sniffable: list[str] = field(default_factory=list)      # small tables in those folders, to be told by content
    over_budget: list[str] = field(default_factory=list)    # first-level subtrees given up
    unreadable: list[str] = field(default_factory=list)     # folders that could not be listed (no permission, I/O error)
    problems: list[str] = field(default_factory=list)
    counts: dict[str, int] = field(default_factory=dict)
    seconds_runs: float = 0.0
    seconds_files: float = 0.0


@dataclass(frozen=True)
class WalkSettings:
    prune: tuple[str, ...] = ()              # the configured list, on top of the built-in one
    budget: int = 20000
    min_age_seconds: float = 120.0
    skip_paths: tuple[str, ...] = ()         # other roots that lie inside this one: walked on their own
    now: float | None = None                 # for tests

    @classmethod
    def from_config(cls, config: dict, skip_paths: tuple[str, ...] = ()) -> WalkSettings:
        d = config["discovery"]
        return cls(prune=tuple(d["prune"]), budget=int(d["subtree_listing_budget"]),
                   min_age_seconds=float(d["min_age_seconds"]), skip_paths=tuple(skip_paths))


class _Fs:
    """Listings and stat calls, counted. Symlinks answer as neither folder nor file."""

    def __init__(self, counts: dict[str, int], beat: Callable[[], None] | None):
        self.counts = counts
        self.beat = beat
        self.cache: dict[str, list[tuple[str, bool, bool]]] = {}
        self.caching = False             # part 2 lists a folder once even when two patterns pass through it
        self.on_list: Callable[[str, list], None] | None = None
        self.unreadable: list[str] = []  # folders that could not be listed: what lies below them is not known

    def listdir(self, path: str) -> list[tuple[str, bool, bool]]:
        hit = self.cache.get(path)
        if hit is not None:
            return hit
        out: list[tuple[str, bool, bool]] = []
        try:
            with os.scandir(path) as it:
                for e in it:
                    try:
                        if e.is_symlink():
                            self.counts["symlinks"] += 1
                            continue
                        out.append((e.name, e.is_dir(follow_symlinks=False), e.is_file(follow_symlinks=False)))
                    except OSError:
                        self.counts["unreadable"] += 1
        except FileNotFoundError:        # it went away between two listings: nothing to report
            return []
        except OSError:
            self.counts["unreadable"] += 1
            if path not in self.unreadable:
                self.unreadable.append(path)
            return []
        out.sort()
        self.counts["listings"] += 1
        self.counts["entries"] += len(out)
        if self.beat and self.counts["listings"] % 100 == 0:
            self.beat()
        if self.caching:
            if len(self.cache) >= 512:
                self.cache.clear()
            self.cache[path] = out
        if self.on_list:
            self.on_list(path, out)
        return out

    def lstat(self, path: str):
        self.counts["stats"] += 1
        try:
            return os.lstat(path)
        except OSError:
            return None

    def is_dir(self, path: str) -> bool:
        st = self.lstat(path)
        return st is not None and _stat.S_ISDIR(st.st_mode)

    def is_file(self, path: str) -> bool:
        st = self.lstat(path)
        return st is not None and _stat.S_ISREG(st.st_mode)


def classify(dirname: str, entries: list[tuple[str, bool, bool]]) -> str | None:
    """``model.runfolder.classify_listing`` on a listing that knows files from folders."""
    return runfolder.classify_listing(dirname, [n for n, _, _ in entries], [n for n, is_dir, _ in entries if is_dir])


def _evidence(kind: str, names: list[str]) -> tuple[str, ...]:
    pats = runfolder.IDENTITY_PATTERNS if kind == "run" else runfolder.RAW_PATTERNS if kind == "run_raw_only" \
        else runfolder.BASECALLER_OUTPUT
    hit = [n for n in names if any(fnmatch.fnmatchcase(n, p) for p in pats)]
    return tuple(sorted(hit)[:8])


def walk_root(root_id: str, root_path: str, rules: list[InputRule] | tuple[InputRule, ...] = (),
              settings: WalkSettings = WalkSettings(), beat: Callable[[], None] | None = None) -> RootWalk:
    """Walk one root. Raises OSError when the root itself cannot be listed."""
    compiled, problems = compile_rules(rules)
    w = _Walker(root_id, os.path.abspath(root_path), compiled, settings, beat)
    w.out.problems.extend(problems)
    return w.run()


class _Walker:
    def __init__(self, root_id: str, root: str, rules: list[Rule], settings: WalkSettings, beat):
        self.root, self.rules, self.s = root, rules, settings
        self.now = settings.now if settings.now is not None else time.time()
        self.counts = {k: 0 for k in ("listings", "entries", "stats", "pruned", "depth_capped", "unreadable", "symlinks",
                                      "too_new", "excluded")}
        self.fs = _Fs(self.counts, beat)
        self.out = RootWalk(root_id=root_id, path=root, counts=self.counts)
        self.heads = [r for r in rules if r.floating and not r.free_file]
        self.free = [r for r in rules if r.free_file]
        self.fixed = [r for r in rules if not r.floating]
        self.anchor_hits: list[tuple[str, Rule, dict]] = []
        self.found: dict[str, Candidate] = {}
        self.holders: set[str] = set()
        self.tables: dict[str, list[str]] = {}
        self.skip = set(settings.skip_paths)
        self._sub = 0
        self._over = False

    # ------------------------------------------------------------------------------------------ part 1
    def run(self) -> RootWalk:
        t0 = time.monotonic()
        with os.scandir(self.root):                      # raises when the root cannot be listed
            pass
        top = self.fs.listdir(self.root)
        self.out.first_level = [n for n, _, _ in top if not n.startswith(".")]
        if not self._folder(self.root, top):
            for name, is_dir, _ in top:
                if not is_dir:
                    continue
                self._sub, self._over = 0, False
                self._child(self.root, name, depth=0)
                if self._over:
                    self.out.over_budget.append(os.path.join(self.root, name))
        self.out.seconds_runs = time.monotonic() - t0
        t1 = time.monotonic()
        self._part2()
        self.out.seconds_files = time.monotonic() - t1
        self.out.anchors = sorted({a for a, _, _ in self.anchor_hits})
        self.out.holders = sorted(self.holders)
        for d in self.out.holders:
            self.out.sniffable.extend(os.path.join(d, n) for n in self.tables.get(d, [])[:SNIFF_MAX_PER_FOLDER])
        self.out.candidates = sorted(self.found.values(), key=lambda c: c.path)
        self.out.unreadable = sorted(self.fs.unreadable)
        return self.out

    def _rel(self, path: str) -> str:
        return "" if path == self.root else os.path.relpath(path, self.root)

    def _folder(self, d: str, entries) -> bool:
        """Look at one listed folder. True when it is a run folder or a look-alike (the walk does not go below)."""
        names = [n for n, _, _ in entries if not n.startswith(".")]
        if self.free:
            self._free_files(d, entries)
        tables = [n for n, _, is_file in entries if is_file and n.lower().endswith(SNIFF_SUFFIXES) and not n.startswith(".")]
        if tables:
            self.tables[d] = tables
        kind = classify(os.path.basename(d), entries)
        if kind is None:
            return False
        found = FoundFolder(kind, d, self._rel(d), _evidence(kind, names),
                            self._run_files(d, entries) if kind == "run" else ())
        (self.out.lookalikes if kind == "lookalike" else self.out.runs).append(found)
        self.holders.add(os.path.dirname(d) if d != self.root else d)
        for name, is_dir, _ in entries:                 # basecaller output inside it may be the anchor of a rule
            if is_dir and not name.startswith("."):
                self._anchor(d, name, hold=False)
        return True

    def _run_files(self, d: str, entries) -> tuple[tuple[str, int, int], ...]:
        """One stat per MinKNOW file of a run folder: enough to know, without opening anything, that a run the
        library already holds did not change."""
        names = [n for n, _, is_file in entries if is_file and not n.startswith(".")]
        if any(n == "other_reports" and is_dir for n, is_dir, _ in entries):
            names += ["other_reports/" + n for n, _, is_file in self.fs.listdir(os.path.join(d, "other_reports")) if is_file]
        out = []
        for n in names:
            if any(fnmatch.fnmatchcase(n, p) for p in _ROLE_PATTERNS):
                st = self.fs.lstat(os.path.join(d, n))
                if st is not None and _stat.S_ISREG(st.st_mode):
                    out.append((n, st.st_size, st.st_mtime_ns))
        return tuple(sorted(out))

    def _anchor(self, parent: str, name: str, hold: bool = True) -> bool:
        hit = False
        for r in self.heads:
            caps = r.head.match(name, {})
            if caps is not None:
                self.anchor_hits.append((os.path.join(parent, name), r, caps))
                hit = True
        if hit and hold:
            self.holders.add(parent)
        return hit

    def _child(self, parent: str, name: str, depth: int) -> None:
        path = os.path.join(parent, name)
        if name.startswith(".") or path in self.skip:
            return
        if any(fnmatch.fnmatchcase(name, p) for p in self.s.prune):
            self.counts["pruned"] += 1
            return
        if self._anchor(parent, name):
            return
        if runfolder.is_pruned(name):
            self.counts["pruned"] += 1
            return
        self._visit(path, depth + 1)

    def _visit(self, d: str, depth: int) -> None:
        if self._over:
            return
        if self._sub >= self.s.budget:
            self._over = True
            return
        self._sub += 1
        entries = self.fs.listdir(d)
        if self._folder(d, entries):
            return
        if depth >= MAX_DEPTH:
            self.counts["depth_capped"] += 1
            return
        for name, is_dir, _ in entries:
            if is_dir:
                self._child(d, name, depth)
                if self._over:
                    return

    # ------------------------------------------------------------------------------------------ part 2
    def _emit(self, rule: Rule, path: str, caps: dict) -> None:
        if path in self.found:
            return
        rel = self._rel(path)
        name = os.path.basename(path)
        if name.startswith("."):
            return
        if rule.excluded(name, rel):
            self.counts["excluded"] += 1
            return
        st = self.fs.lstat(path)
        if st is None or not _stat.S_ISREG(st.st_mode):
            return
        if self.now - st.st_mtime < self.s.min_age_seconds:
            self.counts["too_new"] += 1
            return
        self.found[path] = Candidate(role=rule.rule.role, path=path, size=st.st_size, mtime_ns=st.st_mtime_ns,
                                     root_id=self.out.root_id, relpath=rel, sample=caps.get("sample"),
                                     variant=caps.get("variant") or rule.rule.variant, primary=rule.rule.primary,
                                     rule=rule.words)

    def _seen_below_anchor(self, d: str, entries) -> None:
        """Basecaller output named like a run, met while a pattern is resolved: reported, never taken for a run."""
        if classify(os.path.basename(d), entries) == "lookalike" and not any(f.path == d for f in self.out.lookalikes):
            names = [n for n, _, _ in entries if not n.startswith(".")]
            self.out.lookalikes.append(FoundFolder("lookalike", d, self._rel(d), _evidence("lookalike", names)))

    def _free_files(self, d: str, entries) -> None:
        for name, _, is_file in entries:
            if not is_file or name.startswith("."):
                continue
            for r in self.free:
                caps = r.head.match(name, {})
                if caps is not None:
                    self._emit(r, os.path.join(d, name), caps)
                    break

    def _part2(self) -> None:
        def pruned(name: str) -> bool:
            return runfolder.is_pruned(name, self.s.prune)
        self.fs.caching = True
        self.fs.on_list = self._seen_below_anchor
        jobs: list[tuple[int, str, Rule, dict, tuple]] = []
        for path, rule, caps in self.anchor_hits:
            jobs.append((rule.index, path, rule, caps, rule.parts[1:]))
        for rule in self.fixed:
            jobs.append((rule.index, self.root, rule, {}, rule.parts))
        for _, base, rule, caps, parts in sorted(jobs, key=lambda j: (j[0], j[1])):   # the first rule wins a file
            descend(self.fs, base, parts, caps, 0, lambda p, c, r=rule: self._emit(r, p, c), pruned)


def experiment_dir(path: str, root: str, holders: set[str]) -> str | None:
    """The nearest folder above ``path`` (inside ``root``) that holds a run folder, a sample sheet or an anchor."""
    d = os.path.dirname(path)
    while True:
        if d in holders:
            return d
        if d == root or len(d) <= len(root):
            return None
        parent = os.path.dirname(d)
        if parent == d:
            return None
        d = parent
