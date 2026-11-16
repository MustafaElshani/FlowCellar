"""Helpers shared by the library tests (kept out of conftest.py so that they can be imported by name)."""
from __future__ import annotations

import os

from flowcellar.library import scan as scan_mod


def run_scan(lib_dir: str, *roots: str, **options):
    """One scan through the Python entry point (the command line is tested separately)."""
    options.setdefault("wait", False)
    return scan_mod.scan(lib_dir, scan_mod.ScanOptions(roots=tuple(roots), **options))


def log_events(result, name: str | None = None) -> list[dict]:
    """The lines of a scan's log file, optionally of one event name."""
    import json
    with open(result.log) as fh:
        lines = [json.loads(x) for x in fh if x.strip()]
    return [x for x in lines if name is None or x["event"] == name]


def set_root(lib_dir: str, root_id: str, **changes) -> None:
    """Edit one ``[[roots]]`` entry of a library's flowcellar.toml the way a person would."""
    path = os.path.join(lib_dir, "flowcellar.toml")
    lines = open(path).read().split("\n")
    out, inside = [], False
    for line in lines:
        if line.strip() == "[[roots]]":
            inside = False
        if line.strip() == f'id = "{root_id}"':
            inside = True
        key = line.split("=")[0].strip()
        if inside and key in changes:
            v = changes.pop(key)
            line = f"{key} = {str(v).lower() if isinstance(v, bool) else v if isinstance(v, int) else chr(34) + str(v) + chr(34)}"
        out.append(line)
        if inside and line.strip().startswith("path =") and changes:
            pass
    assert not changes, f"keys not found in the root entry: {changes}"
    open(path, "w").write("\n".join(out))
