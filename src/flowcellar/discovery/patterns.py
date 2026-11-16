"""Input patterns: globs with named captures, and directed descent from an anchor.

A pattern is a path of components relative to a root::

    **/analysis/alignments/{sample}_sorted.bam
    **/*_basecalled/{variant:main|alt}/{sample}.bam
    reads/{sample}/*.fastq.gz

``**`` crosses folders, ``*`` and ``?`` stay inside one component, ``{name}`` captures a fragment of one component
and ``{name:a|b}`` one of the listed words. A pattern that starts with ``**`` *floats*: the component after it (its
*head*) is looked for as a folder name during the generic walk, and every folder that matches is an **anchor** the
pattern is resolved from. Resolution is a directed descent: a folder is listed only when the pattern can still
match below it, a literal tail is tested with one ``stat``, and ``**`` goes at most six levels below the anchor.
"""
from __future__ import annotations

import fnmatch
import os
import re
from dataclasses import dataclass
from typing import Callable, Iterable, Protocol

from flowcellar.contracts import InputRule

MAX_BELOW_ANCHOR = 6
_CAPTURE = re.compile(r"\{([A-Za-z_][A-Za-z0-9_]*)(?::([^{}]*))?\}")
_WILD = re.compile(r"[*?\[{]")


class PatternError(ValueError):
    """A glob that cannot be read; the message says which and why."""


@dataclass(frozen=True)
class Part:
    kind: str                       # "deep" (**), "literal" or "pattern"
    text: str
    regex: re.Pattern | None = None

    def match(self, name: str, caps: dict) -> dict | None:
        """The captures after matching ``name``, or None. A capture seen before must repeat its value."""
        if self.kind == "literal":
            return caps if name == self.text else None
        m = self.regex.fullmatch(name)
        if not m:
            return None
        new = m.groupdict()
        for k, v in new.items():
            if k in caps and caps[k] != v:
                return None
        return {**caps, **new} if new else caps


def _component(text: str, glob: str) -> Part:
    if text == "**":
        return Part("deep", text)
    if "**" in text:
        raise PatternError(f"{glob!r}: '**' must be a whole path component")
    if not _WILD.search(text):
        return Part("literal", text)
    out, i, seen = [], 0, set()
    while i < len(text):
        ch = text[i]
        if ch == "*":
            out.append("[^/]*")
        elif ch == "?":
            out.append("[^/]")
        elif ch == "[":
            j = text.find("]", i + 2)
            if j < 0:
                out.append(re.escape(ch))
            else:
                body = text[i + 1:j]
                body = "^" + body[1:] if body.startswith("!") else body
                out.append("[" + body.replace("\\", "\\\\") + "]")
                i = j
        elif ch == "{":
            m = _CAPTURE.match(text, i)
            if not m:
                raise PatternError(f"{glob!r}: unreadable capture at {text[i:]!r}; write {{sample}} or {{variant:a|b}}")
            name, alts = m.group(1), m.group(2)
            if name in seen:
                out.append(f"(?P={name})")
            elif alts is not None:
                words = [w for w in alts.split("|") if w]
                if not words:
                    raise PatternError(f"{glob!r}: the capture {{{name}:...}} lists no word")
                out.append(f"(?P<{name}>" + "|".join(re.escape(w) for w in words) + ")")
            else:
                out.append(f"(?P<{name}>[^/]+?)")
            seen.add(name)
            i = m.end() - 1
        else:
            out.append(re.escape(ch))
        i += 1
    return Part("pattern", text, re.compile("".join(out)))


def parse_glob(glob: str) -> tuple[Part, ...]:
    """The components of a glob. Empty components and a leading ``./`` are dropped; ``**/**`` is one ``**``."""
    comps = [c for c in glob.strip().split("/") if c and c != "."]
    if not comps or ".." in comps:
        raise PatternError(f"{glob!r}: not a usable pattern")
    parts: list[Part] = []
    for c in comps:
        p = _component(c, glob)
        if p.kind == "deep" and parts and parts[-1].kind == "deep":
            continue
        parts.append(p)
    if parts[-1].kind == "deep":
        raise PatternError(f"{glob!r}: a pattern must end in a file name, not in '**'")
    return tuple(parts)


@dataclass(frozen=True)
class Rule:
    """One compiled ``[[inputs]]`` rule (or built-in default) with a glob."""
    index: int                      # position in the rule list: the first rule that matches a file wins
    rule: InputRule
    parts: tuple[Part, ...]         # without the leading ** of a floating pattern
    floating: bool

    @property
    def head(self) -> Part:
        return self.parts[0]

    @property
    def free_file(self) -> bool:
        """``**/name``: a file looked for in every folder the generic walk lists."""
        return self.floating and len(self.parts) == 1

    @property
    def words(self) -> str:
        where = "flowcellar.toml" if self.rule.source == "config" else "built-in"
        extra = f", variant {self.rule.variant}" if self.rule.variant else ""
        return f"{where} rule {self.rule.role}: {self.rule.glob}{extra}"

    def excluded(self, name: str, relpath: str) -> bool:
        return any(fnmatch.fnmatchcase(name, x) or fnmatch.fnmatchcase(relpath, x) for x in self.rule.exclude)


def compile_rules(rules: Iterable[InputRule]) -> tuple[list[Rule], list[str]]:
    """(compiled rules that have a glob, problems). A rule that cannot be read is left out and named."""
    out, problems = [], []
    for i, r in enumerate(rules):
        if not r.glob:
            continue
        try:
            parts = parse_glob(r.glob)
        except PatternError as e:
            problems.append(str(e))
            continue
        floating = parts[0].kind == "deep"
        out.append(Rule(i, r, parts[1:] if floating else parts, floating))
    return out, problems


def rules_from_config(config: dict, default_rules: Iterable[InputRule] = ()) -> list[InputRule]:
    """The ``[[inputs]]`` entries of a configuration as InputRule, followed by the shipped defaults."""
    out = []
    for t in config.get("inputs", []):
        out.append(InputRule(role=t["role"], glob=t.get("glob"), sniff=t.get("sniff"), variant=t.get("variant"),
                             primary=t.get("primary"), exclude=tuple(t.get("exclude", ())), source="config"))
    out.extend(default_rules)
    return out


class Lister(Protocol):
    def listdir(self, path: str) -> list[tuple[str, bool, bool]]: ...        # (name, is a folder, is a file)
    def is_dir(self, path: str) -> bool: ...
    def is_file(self, path: str) -> bool: ...


def descend(fs: Lister, base: str, parts: tuple[Part, ...], caps: dict, depth: int,
            emit: Callable[[str, dict], None], pruned: Callable[[str], bool]) -> None:
    """Resolve ``parts`` below ``base``; ``depth`` counts levels below the anchor. Symlinks are never followed
    (the lister answers False for them)."""
    if not parts:
        return
    p, last = parts[0], len(parts) == 1
    if p.kind == "deep":
        descend(fs, base, parts[1:], caps, depth, emit, pruned)            # ** matches no level at all
        if depth >= MAX_BELOW_ANCHOR:
            return
        for name, is_dir, _ in fs.listdir(base):
            if is_dir and not name.startswith(".") and not pruned(name):
                descend(fs, os.path.join(base, name), parts, caps, depth + 1, emit, pruned)
        return
    if p.kind == "literal":
        if all(q.kind == "literal" for q in parts):                         # a literal tail: one stat
            full = os.path.join(base, *[q.text for q in parts])
            if fs.is_file(full):
                emit(full, caps)
            return
        nxt = os.path.join(base, p.text)
        if fs.is_dir(nxt):
            descend(fs, nxt, parts[1:], caps, depth + 1, emit, pruned)
        return
    for name, is_dir, is_file in fs.listdir(base):
        if name.startswith("."):
            continue
        c = p.match(name, caps)
        if c is None:
            continue
        if last:
            if is_file:
                emit(os.path.join(base, name), c)
        elif is_dir:
            descend(fs, os.path.join(base, name), parts[1:], c, depth + 1, emit, pruned)
