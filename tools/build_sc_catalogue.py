#!/usr/bin/env python3
"""Build the JSON the single-cell and spatial add-on loads at run time from its catalogue source.

    python tools/build_sc_catalogue.py            # catalogue/metrics_singlecell_spatial.yaml -> src/flowcellar/sc/catalogue.json
    python tools/build_sc_catalogue.py --check    # exit 1 when the JSON is not up to date

PyYAML is a development dependency only; the package never imports it.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parents[1]
SRC = REPO / "catalogue" / "metrics_singlecell_spatial.yaml"
DST = REPO / "src" / "flowcellar" / "sc" / "catalogue.json"
# No id that looks like a real flow cell belongs in the catalogue.
FORBIDDEN = re.compile(r"(?<![A-Za-z0-9])[PF][A-Z]{2}[0-9]{5}(?![0-9])")
# Working remarks do not ship in the package.
WORKING_REMARKS = re.compile(r"(?i)\bobserved\b|\bverified\b|lab baseline|on this machine")
CHECK_KINDS = {"min", "max", "range", "rule"}


def build() -> str:
    cat = yaml.safe_load(SRC.read_text("utf-8"))
    ids = [m["id"] for m in cat["metrics"]]
    dup = sorted({i for i in ids if ids.count(i) > 1})
    if dup:
        raise SystemExit(f"duplicate metric ids: {', '.join(dup)}")
    for m in cat["metrics"]:
        d = m.get("display") or {}
        if not d.get("label") or not d.get("short_label"):
            raise SystemExit(f"{m['id']}: display.label and display.short_label are needed")
        if len(d["short_label"]) > 18:
            raise SystemExit(f"{m['id']}: short_label is longer than 18 characters")
    for c in cat.get("checks", []):
        if c["id"] not in ids:
            raise SystemExit(f"check {c['id']}: no such metric")
        if c.get("kind") not in CHECK_KINDS or c.get("role") not in ("key", "supporting", "info"):
            raise SystemExit(f"check {c['id']}: kind or role not known")
        if not c.get("ok"):
            raise SystemExit(f"check {c['id']}: the sentence for a check that passed is missing")
    text = json.dumps(cat, ensure_ascii=False, separators=(",", ":"), sort_keys=False) + "\n"
    bad = FORBIDDEN.search(text)
    if bad:
        raise SystemExit(f"forbidden text in the catalogue: {bad.group(0)!r}")
    lab = WORKING_REMARKS.search(text)
    if lab:
        raise SystemExit(f"a working remark is in the catalogue: {lab.group(0)!r}")
    return text


def main(argv: list[str]) -> int:
    text = build()
    if "--check" in argv:
        if not DST.exists() or DST.read_text("utf-8") != text:
            print(f"{DST.relative_to(REPO)} is not up to date: run tools/build_sc_catalogue.py", file=sys.stderr)
            return 1
        return 0
    DST.parent.mkdir(parents=True, exist_ok=True)
    DST.write_text(text, "utf-8")
    print(f"wrote {DST.relative_to(REPO)} ({len(text):,} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
