"""Projects: named from the folder structure of each root.

For a path under a root, the project is the path component at ``project_depth`` (1 = the first folder under
the root) or the capture ``{project}`` of ``project_pattern``. The component must lie ABOVE the run folder or
the file: a run directly under its root has no project. The project id is the name in lower case with every
run of characters outside a-z and 0-9 replaced by ``_``; equal ids under two roots are one project.
"""
from __future__ import annotations

from functools import lru_cache

from flowcellar.config.loader import slug
from flowcellar.discovery.patterns import Part, PatternError, parse_glob

NO_PROJECT = "_none"


def project_id(name: str) -> str:
    """Starts and ends with a letter or a digit (the schema's ``safe_id``); ``_none`` stays reserved for "no project"."""
    return slug(name)[:80].strip("_") or "project"


@lru_cache(maxsize=64)
def _pattern(glob: str) -> tuple[Part, ...] | None:
    """The components of a project pattern. It describes the START of a path: a trailing ``/**`` is implied."""
    g = glob.strip().rstrip("/")
    while g.endswith("/**") or g.endswith("/*"):
        g = g.rsplit("/", 1)[0]
    try:
        return parse_glob(g)
    except PatternError:
        return None


def _match(parts: tuple[Part, ...], comps: list[str], caps: dict) -> dict | None:
    """Match ``parts`` against the leading components; at least one component (the run folder or file) must be left."""
    if not parts:
        return caps if comps else None
    p = parts[0]
    if p.kind == "deep":
        for skip in range(len(comps) + 1):
            hit = _match(parts[1:], comps[skip:], caps)
            if hit is not None:
                return hit
        return None
    if not comps:
        return None
    c = p.match(comps[0], caps)
    return None if c is None else _match(parts[1:], comps[1:], c)


def project_name(root: dict | None, relpath: str | None, not_names: frozenset[str] = frozenset()) -> str | None:
    """The project a path under ``root`` belongs to, or None. ``not_names`` (lower case): folder names that are
    never a project (``config.loader.not_project_names``): a pipeline's ``analysis`` folder at project depth
    names no project."""
    if not root or not relpath:
        return None
    comps = [c for c in relpath.split("/") if c and c != "."]
    if root.get("project_pattern"):
        parts = _pattern(root["project_pattern"])
        hit = _match(parts, comps, {}) if parts else None
        name = (hit or {}).get("project") or None
    else:
        depth = int(root.get("project_depth", 1) or 0)
        if depth < 1 or len(comps) <= depth:
            return None
        name = comps[depth - 1]
    return None if name is not None and name.lower() in not_names else name


def merge_projects(found: list[tuple[str | None, str | None]], root_order: list[str]) -> list[dict]:
    """``[(project name, root id), ...]`` in order of importance -> the ``relations.projects`` list: one entry per
    project id, roots in configuration order."""
    out: dict[str, dict] = {}
    for name, root_id in found:
        if not name:
            continue
        pid = project_id(name)
        p = out.setdefault(pid, {"id": pid, "name": name[:256], "source": "path", "roots": []})
        if root_id and root_id not in p["roots"]:
            p["roots"].append(root_id)
    for p in out.values():
        p["roots"].sort(key=lambda r: root_order.index(r) if r in root_order else len(root_order))
    return list(out.values())
