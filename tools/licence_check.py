#!/usr/bin/env python3
"""Check the licences of everything FlowCellar ships or needs at run time.

    python tools/licence_check.py            exit 1 on a licence outside the allow-list or a missing notice

Two things are looked at:

* **What is bundled.** The packages of ``viewer/package-lock.json`` that are not development tools end up
  inside the viewer script, and therefore inside every wheel and every report. Each must carry a licence of
  the allow-list and must be named, with its exact version, in ``THIRD_PARTY_NOTICES.md``.
* **What is installed beside it.** The distributions that ``flowcellar[all]`` pulls in (numpy, pyarrow, click,
  pysam, jsonschema and what those need), read from the environment this runs in. Each must carry a licence
  of the allow-list. A distribution that is not installed here is reported and skipped.

FlowCellar is GPL-3.0-or-later: only permissive licences are allowed in, all of which may be combined with it,
and nothing that would put further conditions on the whole. Standard library only.
"""
from __future__ import annotations

import json
import re
import sys
from importlib import metadata
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ALLOWED = {
    "MIT", "MIT-0", "BSD-2-Clause", "BSD-3-Clause", "Apache-2.0", "ISC", "0BSD", "Zlib", "PSF-2.0", "Python-2.0",
    "CC0-1.0", "Unlicense", "BSL-1.0", "LLVM-exception", "GCC-exception-3.1", "CC-BY-4.0", "BSD-3-Clause-Clear",
}
# Older packages state their licence as free text or as a classifier: what those say, in the terms above.
ALIASES = {
    "mit license": "MIT", "mit": "MIT", "the mit license (mit)": "MIT", "bsd license": "BSD-3-Clause", "bsd": "BSD-3-Clause",
    "bsd-3-clause": "BSD-3-Clause", "bsd 3-clause": "BSD-3-Clause", "new bsd": "BSD-3-Clause",
    "apache software license": "Apache-2.0", "apache license, version 2.0": "Apache-2.0", "apache 2.0": "Apache-2.0",
    "apache-2.0": "Apache-2.0", "python software foundation license": "PSF-2.0", "psf-2.0": "PSF-2.0", "isc license (iscl)": "ISC",
}


def bundled_packages() -> list[dict]:
    lock = json.loads((REPO / "viewer" / "package-lock.json").read_text("utf-8"))
    out = []
    for path, entry in lock.get("packages", {}).items():
        if path and not entry.get("dev"):
            out.append({"name": path.split("node_modules/")[-1], "version": entry.get("version"), "licence": entry.get("license")})
    return sorted(out, key=lambda p: p["name"])


def _expression_ids(expression: str) -> list[str]:
    return [t for t in re.split(r"\s+(?:AND|OR|WITH)\s+|[()\s]+", expression) if t and t not in ("AND", "OR", "WITH")]


def licence_ids(dist: metadata.Distribution) -> list[str]:
    meta = dist.metadata
    if meta.get("License-Expression"):
        return _expression_ids(meta["License-Expression"])
    names = [c.rsplit(" :: ", 1)[-1] for c in meta.get_all("Classifier", []) if c.startswith("License ::")]
    if not names and meta.get("License"):
        first = meta["License"].strip().splitlines()[0].strip()
        names = [first] if len(first) < 60 else []
    out = []
    for n in names:
        key = n.strip().lower()
        out.append(ALIASES.get(key) or (n if n in ALLOWED else next((v for k, v in ALIASES.items() if key.startswith(k)), n)))
    return out


def closure(root: str = "flowcellar", extras: tuple[str, ...] = ("bam", "validate", "sc")) -> tuple[list[str], list[str]]:
    """(installed distributions that flowcellar with these extras needs, names that are not installed here)."""
    try:
        requires = metadata.requires(root) or []
    except metadata.PackageNotFoundError:
        raise SystemExit("flowcellar is not installed in this environment (pip install -e .)") from None
    seen: dict[str, bool] = {}
    todo = []
    for req in requires:
        m = re.search(r"extra\s*==\s*['\"]([^'\"]+)['\"]", req)
        if m is None or m.group(1) in extras:
            todo.append(req)
    while todo:
        req = todo.pop()
        name = re.split(r"[ ;<>=!~\[(]", req.strip(), maxsplit=1)[0]
        key = name.lower().replace("_", "-")
        if key in seen or key == root or not key:
            continue
        try:
            sub = metadata.requires(name) or []
            seen[key] = True
        except metadata.PackageNotFoundError:
            seen[key] = False
            continue
        todo += [r for r in sub if "extra ==" not in r]
    return sorted(k for k, v in seen.items() if v), sorted(k for k, v in seen.items() if not v)


def main() -> int:
    problems = []
    notices = (REPO / "THIRD_PARTY_NOTICES.md").read_text("utf-8")
    print("bundled in the viewer (viewer/package-lock.json, not development tools):")
    for p in bundled_packages():
        named = re.search(r"\[" + re.escape(p["name"]) + r"\]\([^)]*\)\s*\|\s*" + re.escape(p["version"]) + r"\s*\|", notices)
        ok = p["licence"] in ALLOWED and bool(named)
        print(f"  {'ok  ' if ok else 'FAIL'}  {p['name']} {p['version']}  {p['licence']}" + ("" if named else "  (not in THIRD_PARTY_NOTICES.md with this version)"))
        if not ok:
            problems.append(p["name"])
    installed, absent = closure()
    print("installed beside FlowCellar (flowcellar[all] and what it needs):")
    for name in installed:
        dist = metadata.distribution(name)
        ids = licence_ids(dist)
        ok = bool(ids) and all(i in ALLOWED for i in ids)
        print(f"  {'ok  ' if ok else 'FAIL'}  {dist.metadata['Name']} {dist.version}  {' + '.join(ids) or 'no licence stated'}")
        if not ok:
            problems.append(name)
    for name in absent:
        print(f"  skip  {name}  (not installed in this environment)")
    print("FAILED: " + ", ".join(problems) if problems else "all licences are on the allow-list")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
