"""Finding read-data files by pattern, for ``flowcellar run`` without a library.

The same glob language as the ``[[inputs]]`` rules of a library: ``**`` crosses folders, ``*`` and ``?`` stay inside
one path component, ``{sample}`` captures a fragment of one component and ``{variant:a|b}`` one of the listed words.
A pattern that starts with ``**`` is resolved from its **anchors**: the folders whose name matches the component
after ``**``. Anchors are found by a pruned walk that does not enter them; below an anchor the descent is directed:
a folder is listed only if the pattern can still match below it, a literal tail costs one ``stat``, and ``**``
inside a pattern goes at most six levels below the anchor. Symbolic links are never followed.

(The scan of a library has its own walk over whole roots; this module serves the paths given on a command line.)
"""
from __future__ import annotations

import fnmatch
import os
import re
import stat as _stat
import time
from dataclasses import dataclass

from flowcellar.contracts import Candidate, InputRule
from flowcellar.model import runfolder

MAX_BELOW_ANCHOR = 6
MAX_WALK_DEPTH = 12
MAX_LISTINGS = 20_000
_CAPTURE = re.compile(r"\{([A-Za-z_][A-Za-z0-9_]*)(?::([^{}]*))?\}")


class PatternError(ValueError):
    pass


def _component(text: str) -> tuple[str, re.Pattern | None]:
    """("deep" | "literal" | "pattern", regex)."""
    if text == "**":
        return "deep", None
    if not re.search(r"[*?{]", text):
        return "literal", None
    out, i, seen = [], 0, set()
    while i < len(text):
        ch = text[i]
        if ch == "*":
            out.append("[^/]*")
        elif ch == "?":
            out.append("[^/]")
        elif ch == "{":
            m = _CAPTURE.match(text, i)
            if not m:
                raise PatternError(f"unreadable capture at {text[i:]!r}; write {{sample}} or {{variant:a|b}}")
            name, alts = m.group(1), m.group(2)
            if name in seen:
                out.append(f"(?P={name})")
            elif alts is not None:
                out.append(f"(?P<{name}>" + "|".join(re.escape(w) for w in alts.split("|") if w) + ")")
            else:
                out.append(f"(?P<{name}>[^/]+?)")
            seen.add(name)
            i = m.end() - 1
        else:
            out.append(re.escape(ch))
        i += 1
    return "pattern", re.compile("".join(out))


@dataclass(frozen=True)
class Compiled:
    rule: InputRule
    parts: tuple            # ((kind, text, regex), ...) without the leading **
    floating: bool

    @property
    def words(self) -> str:
        where = "configured" if self.rule.source == "config" else "built-in"
        return f"{where} rule {self.rule.role}: {self.rule.glob}" + (f", variant {self.rule.variant}" if self.rule.variant else "")


def compile_rule(rule: InputRule) -> Compiled:
    comps = [c for c in (rule.glob or "").strip().split("/") if c and c != "."]
    if not comps or ".." in comps or comps[-1] == "**":
        raise PatternError(f"{rule.glob!r}: not a usable pattern")
    parts = []
    for c in comps:
        kind, rx = _component(c)
        if kind == "deep" and parts and parts[-1][0] == "deep":
            continue
        parts.append((kind, c, rx))
    floating = parts[0][0] == "deep"
    return Compiled(rule, tuple(parts[1:] if floating else parts), floating)


def _match(part, name: str, caps: dict) -> dict | None:
    kind, text, rx = part
    if kind == "literal":
        return caps if name == text else None
    m = rx.fullmatch(name)
    if not m:
        return None
    new = m.groupdict()
    if any(k in caps and caps[k] != v for k, v in new.items()):
        return None
    return {**caps, **new} if new else caps


class _Finder:
    def __init__(self, min_age_seconds: float = 0.0, prune: tuple[str, ...] = ()):
        self.listings = 0
        self.min_age = min_age_seconds
        self.prune = prune
        self.now = time.time()
        self.found: dict[str, tuple] = {}          # path -> (compiled, caps, stat)

    def listdir(self, path: str) -> list[tuple[str, bool]]:
        """(name, is a folder) of the entries that are plain files or folders; links are left out."""
        self.listings += 1
        out = []
        try:
            with os.scandir(path) as it:
                for e in it:
                    try:
                        if e.is_symlink():
                            continue
                        out.append((e.name, e.is_dir(follow_symlinks=False)))
                    except OSError:
                        continue
        except OSError:
            return []
        out.sort()
        return out

    def emit(self, c: Compiled, path: str, caps: dict) -> None:
        if path in self.found:
            return                                   # the first rule that matches a file wins
        name = os.path.basename(path)
        if name.startswith(".") or any(fnmatch.fnmatchcase(name, x) or fnmatch.fnmatchcase(path, x) for x in c.rule.exclude):
            return
        try:
            st = os.lstat(path)
        except OSError:
            return
        if not _stat.S_ISREG(st.st_mode) or self.now - st.st_mtime < self.min_age:
            return
        self.found[path] = (c, caps, st)

    def descend(self, c: Compiled, base: str, parts: tuple, caps: dict, depth: int) -> None:
        if not parts or self.listings > MAX_LISTINGS:
            return
        part, last = parts[0], len(parts) == 1
        if part[0] == "deep":
            self.descend(c, base, parts[1:], caps, depth)
            if depth < MAX_BELOW_ANCHOR:
                for name, is_dir in self.listdir(base):
                    if is_dir and not name.startswith("."):
                        self.descend(c, os.path.join(base, name), parts, caps, depth + 1)
            return
        if part[0] == "literal":
            if all(p[0] == "literal" for p in parts):            # a literal tail: one stat, no listing
                self.emit(c, os.path.join(base, *[p[1] for p in parts]), caps)
                return
            nxt = os.path.join(base, part[1])
            if os.path.isdir(nxt) and not os.path.islink(nxt):
                self.descend(c, nxt, parts[1:], caps, depth + 1)
            return
        for name, is_dir in self.listdir(base):
            if name.startswith("."):
                continue
            got = _match(part, name, caps)
            if got is None:
                continue
            full = os.path.join(base, name)
            if last:
                if not is_dir:
                    self.emit(c, full, got)
            elif is_dir:
                self.descend(c, full, parts[1:], got, depth + 1)

    def walk(self, top: str, compiled: list[Compiled]) -> None:
        """Find the anchors of the floating patterns below ``top`` (pruned, anchors are not entered) and resolve
        every pattern. A pattern that does not float is resolved from ``top``."""
        floating = [c for c in compiled if c.floating]
        for c in compiled:
            if not c.floating:
                self.descend(c, top, c.parts, {}, 0)
        if not floating:
            return
        stack = [(top, 0)]
        while stack and self.listings <= MAX_LISTINGS:
            d, depth = stack.pop()
            for name, is_dir in self.listdir(d):
                if name.startswith("."):
                    continue
                full = os.path.join(d, name)
                hit = False
                for c in floating:
                    got = _match(c.parts[0], name, {})
                    if got is None:
                        continue
                    if len(c.parts) == 1:
                        if not is_dir:
                            self.emit(c, full, got)
                    elif is_dir:
                        hit = True
                        self.descend(c, full, c.parts[1:], got, 1)
                if is_dir and not hit and depth < MAX_WALK_DEPTH and not runfolder.is_pruned(name, self.prune):
                    stack.append((full, depth + 1))


def find(top: str, rules: list[InputRule], *, min_age_seconds: float = 0.0, prune: tuple[str, ...] = (),
         root_id: str | None = None, experiment_dir: str | None = None) -> tuple[list[Candidate], list[str]]:
    """The files below ``top`` that the rules with a glob name, as candidates, and the rules that could not be read."""
    compiled, problems = [], []
    for r in rules:
        if not r.glob:
            continue
        try:
            compiled.append(compile_rule(r))
        except PatternError as e:
            problems.append(str(e))
    f = _Finder(min_age_seconds, prune)
    f.walk(top, compiled)
    out = []
    for path in sorted(f.found):
        c, caps, st = f.found[path]
        out.append(Candidate(role=c.rule.role, path=path, size=st.st_size, mtime_ns=st.st_mtime_ns, root_id=root_id,
                             relpath=os.path.relpath(path, top), sample=caps.get("sample"),
                             variant=caps.get("variant") or c.rule.variant, primary=c.rule.primary, rule=c.words,
                             experiment_dir=experiment_dir))
    if f.listings > MAX_LISTINGS:
        problems.append(f"{top}: more than {MAX_LISTINGS:,} folders were listed; the search stopped there")
    return out, problems
