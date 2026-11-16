#!/usr/bin/env python3
"""Build the JSON the package loads at run time from the catalogue source.

    python tools/build_catalogue.py            # catalogue/metrics_nanopore.yaml -> src/flowcellar/catalogue/metrics_nanopore.json
    python tools/build_catalogue.py --check    # exit 1 when the JSON is not up to date

PyYAML is a development dependency only; the package never imports it.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parents[1]
SRC = REPO / "catalogue" / "metrics_nanopore.yaml"
DST = REPO / "src" / "flowcellar" / "catalogue" / "metrics_nanopore.json"


def build() -> str:
    cat = yaml.safe_load(SRC.read_text("utf-8"))
    ids = [m["id"] for m in cat["metrics"]] + [d["id"] for d in cat.get("derived", [])]
    dup = sorted({i for i in ids if ids.count(i) > 1})
    if dup:
        raise SystemExit(f"duplicate metric ids: {', '.join(dup)}")
    for m in cat["metrics"] + cat.get("derived", []):
        short = (m.get("display") or {}).get("short_label")
        if short and len(short) > 18:
            raise SystemExit(f"{m['id']}: short_label is longer than 18 characters")
    # A free-text note on a limit is a working note of whoever tunes it and is never packaged.
    for m in cat["metrics"] + cat.get("derived", []):
        th = m.get("thresholds")
        if isinstance(th, dict):
            th.pop("note", None)
    text = json.dumps(cat, ensure_ascii=False, separators=(",", ":"), sort_keys=False) + "\n"
    return text


def main(argv: list[str]) -> int:
    text = build()
    if "--check" in argv:
        if not DST.exists() or DST.read_text("utf-8") != text:
            print(f"{DST.relative_to(REPO)} is not up to date: run tools/build_catalogue.py", file=sys.stderr)
            return 1
        return 0
    DST.parent.mkdir(parents=True, exist_ok=True)
    DST.write_text(text, "utf-8")
    print(f"wrote {DST.relative_to(REPO)} ({len(text):,} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
