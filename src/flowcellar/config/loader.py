"""``flowcellar.toml``: the defaults, a library's own file, and the few places FlowCellar writes it.

People own the file. FlowCellar writes it when ``scan`` or ``init`` creates a library, and afterwards only ever
ADDS a ``[[roots]]`` entry at its end (``append_roots``); nothing else in it is touched.

    loaded = load("/path/library/flowcellar.toml")
    loaded.config["roots"]      -> [{"id", "path", "live", "project_depth", "project_pattern"}, ...]
    loaded.config["discovery"]  -> {"min_age_seconds": 120, ...}        every default filled in
    loaded.warnings             -> ["unknown key [scan] job (did you mean jobs?)", ...]
    loaded.notices              -> [{"code": "config_not_used", ...}]   keys read and not used by this build

An unknown key is a warning, never an error. A value of the wrong type is an error with the key named.
"""
from __future__ import annotations

import copy
import difflib
import os
import re
import tomllib
from dataclasses import dataclass, field
from importlib import resources

from flowcellar.util import atomic
from flowcellar.util.errors import LibraryError
from flowcellar.util.fingerprint import hash_obj

CONFIG_NAME = "flowcellar.toml"
DEFAULTS_FILE = "defaults.toml"
ROOT_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,79}$")

# section -> key -> accepted types. Tables that hold free keys (aliases, profiles) are checked by hand below.
_SCHEMA: dict[str, dict[str, tuple]] = {
    "library": {"name": (str,)},
    "discovery": {"min_age_seconds": (int, float), "stale_after_hours": (int, float), "prune": (list,),
                  "subtree_listing_budget": (int,), "root_check_seconds": (int, float),
                  "stall_timeout_seconds": (int, float)},
    "variants": {"primary": (str,)},
    "samples": {"unclassified": (list,), "alias": (dict,)},
    "project": {"from": (str,), "exclude": (list,)},
    "baselines": {"window": (int,), "min_runs": (int,), "trusted_runs": (int,), "spread_guard": (int, float),
                  "groups": (list,)},
    "flow_cell": {"restart_max_hours": (int, float), "restart_max_output_frac": (int, float)},
    "scan": {"jobs": (int,), "threads": (int,), "quick_first": (bool,), "time_box_minutes": (int, float),
             "live_workers": (int,)},
    "hub": {"sidecars": (bool,), "reports": (bool,)},
    "serve": {"host": (str,), "port": (int,)},
}
_ROOT_KEYS = {"id": (str,), "path": (str,), "live": (bool,), "project_depth": (int,), "project_pattern": (str,)}
_INPUT_KEYS = {"role": (str,), "glob": (str,), "sniff": (str,), "variant": (str,), "primary": (bool,), "exclude": (list,)}
_LINK_RULE_KEYS = {"name": (str,), "roles": (list,), "via": (str,), "within": (str,), "match": (str,),
                   "confirmed": (bool,), "run": (str,), "flow_cell": (str,)}
_TOP_LEVEL = set(_SCHEMA) | {"roots", "inputs", "link", "profiles", "projects"}
# read, kept out of the way, and reported once each with the notice config_not_used
_NOT_USED = {("hub", "views"), ("hub", "default_group"), ("roots", "minknow_tmp"), ("roots", "max_read_mb_s"),
             ("projects", None)}


def not_project_names(config: dict | None) -> frozenset[str]:
    """Folder names that are never the name of a project, in lower case: the names listed under ``[project]
    exclude``, and every folder name an ``[[inputs]]`` pattern spells out (``**/analysis/alignments/{sample}.bam``
    names ``analysis`` and ``alignments``). By the configuration's own word those are the folders of a pipeline;
    when such a name stands at project depth, the path is simply one level shallower than the project rule
    expects, and the file names no project (the run folder's own project stands)."""
    cfg = config or {}
    names = {str(n).lower() for n in (cfg.get("project") or {}).get("exclude") or [] if isinstance(n, str)}
    for rule in cfg.get("inputs") or []:
        glob = rule.get("glob") if isinstance(rule, dict) else getattr(rule, "glob", None)
        for comp in str(glob or "").split("/")[:-1]:
            if comp and not any(ch in comp for ch in "*?[]{}"):
                names.add(comp.lower())
    return frozenset(names)


class ConfigError(LibraryError):
    """flowcellar.toml cannot be used as it is (exit code 6): the message names the file and the key."""


@dataclass
class Loaded:
    config: dict
    path: str | None = None                 # the file that was read; None = defaults only
    warnings: list[str] = field(default_factory=list)
    notices: list[dict] = field(default_factory=list)


def defaults() -> dict:
    """The shipped defaults as a dict, with the list-valued sections present and empty."""
    text = resources.files("flowcellar.config").joinpath(DEFAULTS_FILE).read_text("utf-8")
    d = tomllib.loads(text)
    d.setdefault("library", {})
    d.setdefault("variants", {})
    d["samples"].setdefault("alias", {})
    d["baselines"].setdefault("groups", [])
    d.update(roots=[], inputs=[], link={"rules": []}, profiles={})
    return d


def _type_name(types: tuple) -> str:
    names = {str: "text", int: "a whole number", float: "a number", bool: "true or false", list: "a list", dict: "a table"}
    return " or ".join(dict.fromkeys(names[t] for t in types))


def _check(value, types: tuple, where: str, source: str) -> None:
    ok = isinstance(value, types) and not (isinstance(value, bool) and bool not in types)
    if not ok:
        raise ConfigError(f"{source}: {where} must be {_type_name(types)}, not {value!r}")


def _closest(key: str, known) -> str:
    m = difflib.get_close_matches(key, list(known), n=1, cutoff=0.7)
    return f" (did you mean {m[0]}?)" if m else ""


def slug(text: str) -> str:
    """Lower case, every run of characters outside a-z and 0-9 replaced by one underscore."""
    return re.sub(r"[^a-z0-9]+", "_", text.lower()).strip("_")


def root_ids_for(paths: list[str], taken: set[str]) -> list[str]:
    """An id for each new root path: its last path component, made unique. Two new roots that end in the same name
    both get their parent folder in front (``/a/projects``, ``/b/projects`` -> ``a_projects``, ``b_projects``)."""
    base = [slug(os.path.basename(p.rstrip("/"))) or "root" for p in paths]
    out: list[str] = []
    used = set(taken)
    for p, b in zip(paths, base):
        cand = b
        if base.count(b) > 1 or b in used:
            parent = slug(os.path.basename(os.path.dirname(p.rstrip("/"))))
            cand = f"{parent}_{b}" if parent else b
        n, unique = 2, cand
        while unique in used:
            unique, n = f"{cand}_{n}", n + 1
        used.add(unique)
        out.append(unique[:80])
    return out


def normalise_root_path(path: str) -> str:
    return os.path.abspath(os.path.expanduser(path))


def _roots(raw, source: str, loaded: Loaded) -> list[dict]:
    if not isinstance(raw, list):
        raise ConfigError(f"{source}: roots must be written as [[roots]] tables")
    out: list[dict] = []
    for i, r in enumerate(raw, 1):
        where = f"[[roots]] entry {i}"
        if not isinstance(r, dict) or "path" not in r:
            raise ConfigError(f"{source}: {where} needs a path")
        for k, v in r.items():
            if ("roots", k) in _NOT_USED:
                _not_used(loaded, f"roots.{k}", f"{where}: '{k}' is read but not used by this version of FlowCellar")
            elif k not in _ROOT_KEYS:
                loaded.warnings.append(f"{source}: unknown key '{k}' in {where}{_closest(k, _ROOT_KEYS)}")
            else:
                _check(v, _ROOT_KEYS[k], f"{where}: {k}", source)
        depth = r.get("project_depth", 1)
        if depth < 0:
            raise ConfigError(f"{source}: {where}: project_depth must be 0 or more (0 = no project from this root)")
        pattern = r.get("project_pattern")
        if pattern is not None and "{project}" not in pattern:
            raise ConfigError(f"{source}: {where}: project_pattern must hold the capture {{project}}")
        out.append({"id": r.get("id"), "path": normalise_root_path(r["path"]), "live": bool(r.get("live", False)),
                    "project_depth": depth, "project_pattern": pattern})
    taken = {r["id"] for r in out if r["id"]}
    missing = [r for r in out if not r["id"]]
    for r, rid in zip(missing, root_ids_for([r["path"] for r in missing], taken)):
        r["id"] = rid
    seen: set[str] = set()
    for i, r in enumerate(out, 1):
        if not ROOT_ID_RE.match(r["id"]):
            raise ConfigError(f"{source}: [[roots]] entry {i}: id {r['id']!r} may hold letters, digits, '_', '-' and '.' only")
        if r["id"] in seen:
            raise ConfigError(f"{source}: two [[roots]] entries have the id {r['id']!r}")
        seen.add(r["id"])
    return out


def _not_used(loaded: Loaded, key: str, message: str) -> None:
    if not any(n.get("key") == key for n in loaded.notices):
        loaded.notices.append({"code": "config_not_used", "level": "info", "message": message, "root": None,
                               "since": None, "key": key})


def _list_of_tables(raw, name: str, keys: dict, source: str, loaded: Loaded) -> list[dict]:
    if not isinstance(raw, list):
        raise ConfigError(f"{source}: {name} must be written as [[{name}]] tables")
    out = []
    for i, t in enumerate(raw, 1):
        where = f"[[{name}]] entry {i}"
        if not isinstance(t, dict):
            raise ConfigError(f"{source}: {where} must be a table")
        for k, v in t.items():
            if k not in keys:
                loaded.warnings.append(f"{source}: unknown key '{k}' in {where}{_closest(k, keys)}")
            else:
                _check(v, keys[k], f"{where}: {k}", source)
        out.append(dict(t))
    return out


def merge(raw: dict, source: str = CONFIG_NAME) -> Loaded:
    """The defaults with ``raw`` (a parsed flowcellar.toml) laid over them, checked."""
    loaded = Loaded(config=defaults())
    cfg = loaded.config
    for top, value in raw.items():
        if top not in _TOP_LEVEL:
            loaded.warnings.append(f"{source}: unknown section [{top}]{_closest(top, _TOP_LEVEL)}")
            continue
        if top == "roots":
            cfg["roots"] = _roots(value, source, loaded)
        elif top == "inputs":
            cfg["inputs"] = _list_of_tables(value, "inputs", _INPUT_KEYS, source, loaded)
            for i, t in enumerate(cfg["inputs"], 1):
                if "role" not in t or ("glob" not in t and "sniff" not in t):
                    raise ConfigError(f"{source}: [[inputs]] entry {i} needs a role and a glob or a sniff rule")
        elif top == "link":
            if not isinstance(value, dict):
                raise ConfigError(f"{source}: [link] must be a table")
            for k in value:
                if k != "rules":
                    loaded.warnings.append(f"{source}: unknown key '{k}' in [link]")
            cfg["link"] = {"rules": _list_of_tables(value.get("rules", []), "link.rules", _LINK_RULE_KEYS, source, loaded)}
        elif top == "profiles":
            if not isinstance(value, dict):
                raise ConfigError(f"{source}: [profiles] must be a table of assay groups")
            cfg["profiles"] = copy.deepcopy(value)
        elif top == "projects":
            _not_used(loaded, "projects", "[[projects]] targets are read but not used by this version of FlowCellar")
        else:
            if not isinstance(value, dict):
                raise ConfigError(f"{source}: [{top}] must be a table")
            for k, v in value.items():
                if (top, k) in _NOT_USED:
                    _not_used(loaded, f"{top}.{k}", f"[{top}] {k} is read but not used by this version of FlowCellar")
                    continue
                if k not in _SCHEMA[top]:
                    loaded.warnings.append(f"{source}: unknown key '{k}' in [{top}]{_closest(k, _SCHEMA[top])}")
                    continue
                _check(v, _SCHEMA[top][k], f"[{top}] {k}", source)
                if isinstance(v, dict) and isinstance(cfg[top].get(k), dict):
                    cfg[top][k] = {**cfg[top][k], **v}
                else:
                    cfg[top][k] = copy.deepcopy(v)
    if cfg["project"].get("from", "path") != "path":
        _not_used(loaded, "project.from", f"[project] from = {cfg['project']['from']!r} is not used: projects are named from the "
                                          "folder structure of each root (project_depth or project_pattern)")
    if cfg["scan"]["quick_first"] is False:
        _not_used(loaded, "scan.quick_first", "[scan] quick_first = false is not used: verdicts from the small files always come first")
    for k in ("min_age_seconds", "stale_after_hours", "subtree_listing_budget", "root_check_seconds", "stall_timeout_seconds"):
        if cfg["discovery"][k] < 0:
            raise ConfigError(f"{source}: [discovery] {k} must not be negative")
    if cfg["scan"]["time_box_minutes"] < 0 or cfg["scan"]["jobs"] < 0 or cfg["scan"]["threads"] < 1:
        raise ConfigError(f"{source}: [scan] jobs and time_box_minutes must not be negative, threads must be 1 or more")
    if not all(isinstance(p, str) for p in cfg["discovery"]["prune"]):
        raise ConfigError(f"{source}: [discovery] prune must be a list of folder names")
    for g in cfg["baselines"]["groups"]:
        if not isinstance(g, dict):
            raise ConfigError(f"{source}: [[baselines.groups]] entries must be tables")
    return loaded


def load(path: str | os.PathLike | None) -> Loaded:
    """Read one flowcellar.toml over the defaults. ``None`` or a file that does not exist gives the defaults."""
    if path is None or not os.path.isfile(path):
        return merge({})
    path = os.fspath(path)
    try:
        raw = tomllib.loads(atomic.read_bytes(path).decode("utf-8"))
    except tomllib.TOMLDecodeError as e:
        raise ConfigError(f"{path}: not valid TOML: {e}", hint="fix the line named above; FlowCellar changed nothing") from None
    except UnicodeDecodeError:
        raise ConfigError(f"{path}: not a UTF-8 text file") from None
    loaded = merge(raw, source=path)
    loaded.path = path
    return loaded


def shaping_hash(config: dict) -> str:
    """Hash of the settings that decide which files belong to which run: a change means linking is done again."""
    return hash_obj({k: config.get(k) for k in ("inputs", "variants", "link", "samples", "roots")}, 16)


# ---------------------------------------------------------------------------------------------- writing
def toml_str(s: str) -> str:
    """A TOML basic string."""
    out = s.replace("\\", "\\\\").replace('"', '\\"')
    out = "".join(ch if ch >= " " and ch != "\x7f" else "\\u%04x" % ord(ch) for ch in out)
    return f'"{out}"'


def root_block(root: dict, comment: str | None = None) -> str:
    lines = [f"# {comment}"] if comment else []
    lines += ["[[roots]]", f"id = {toml_str(root['id'])}", f"path = {toml_str(root['path'])}"]
    lines.append(f"live = {'true' if root.get('live') else 'false'}")
    if root.get("project_pattern"):
        lines.append(f"project_pattern = {toml_str(root['project_pattern'])}")
    else:
        lines.append(f"project_depth = {int(root.get('project_depth', 1))}")
    return "\n".join(lines) + "\n"


_HEAD = """\
# FlowCellar library configuration.
#
# This file is yours: edit it freely. FlowCellar reads it at every scan and only ever ADDS a [[roots]]
# entry at the end, when "flowcellar scan NEW_FOLDER -L <this library>" is given a folder it does not know.
# A key that is left out keeps its default. Nothing is ever copied from, moved in or written to a root.

[library]
name = {name}

# ---- Where the library comes from ---------------------------------------------------------------------
# One [[roots]] entry per folder that is scanned.
#   id             short name shown on the Home page and stored with every run; do not change it later
#   live           true for the folder a sequencer is writing to: it is read gently
#   project_depth  which folder level under the root names the project (1 = the first folder under the
#                  root, 2 = the second, 0 = no project from this root); or instead
#   project_pattern = "*/{{project}}/**"   a pattern with the capture {{project}}
"""

_TAIL = """\

# ---- Read data (optional) -----------------------------------------------------------------------------
# FlowCellar finds by itself what MinKNOW and Dorado write inside a run folder. Files that your own
# pipeline writes elsewhere need one [[inputs]] rule each. "{sample}" captures the sample name,
# "{variant:a|b}" one of several versions of a file, "**" crosses folders. Patterns are matched under
# every root. See what a rule finds, before anything is computed, with:
#     flowcellar discover -L <this library> --explain
#
# [[inputs]]
# role = "align.genome_bam"
# glob = "**/analysis/alignments/{sample}_sorted.bam"
#
# [discovery]
# prune = ["notebooks", "figures"]        # folder names the walk never enters, added to the built-in list
#
# [scan]
# time_box_minutes = 240                  # stop starting new passes after this long; 0 = no limit
"""


def render_new(name: str, roots: list[dict]) -> str:
    """The text of the flowcellar.toml that ``scan`` writes when it creates a library."""
    parts = [_HEAD.format(name=toml_str(name))]
    for r in roots:
        parts.append("\n" + root_block(r))
    parts.append(_TAIL)
    return "".join(parts)


def write_new(path: str | os.PathLike, name: str, roots: list[dict]) -> None:
    atomic.atomic_write_text(path, render_new(name, roots))


def append_roots(path: str | os.PathLike, roots: list[dict], when: str) -> None:
    """Add ``[[roots]]`` entries at the end of an existing file (whole-file replace; every other byte is kept)."""
    text = atomic.read_bytes(path).decode("utf-8") if os.path.isfile(path) else ""
    if text and not text.endswith("\n"):
        text += "\n"
    for r in roots:
        text += "\n" + root_block(r, comment=f"added by flowcellar scan on {when[:10]}")
    try:
        tomllib.loads(text)                  # never leave a file behind that cannot be read
    except tomllib.TOMLDecodeError as e:
        raise ConfigError(f"{os.fspath(path)}: a [[roots]] entry cannot be added to this file: {e}",
                          hint="write the roots as [[roots]] tables, then scan again") from None
    atomic.atomic_write_text(path, text)
