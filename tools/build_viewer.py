#!/usr/bin/env python3
"""Build the viewer: viewer/ -> src/flowcellar/viewer_dist/ (script, stylesheet, template, manifest with hashes).

    python tools/build_viewer.py            minified, as shipped
    python tools/build_viewer.py --dev      readable script
    python tools/build_viewer.py --check    build into a temporary folder and say whether viewer_dist/ is up to date

This is a thin wrapper around ``node viewer/build.mjs`` for people who work from Python. Node (22 or newer) and the
packages of viewer/package-lock.json (``npm ci`` in viewer/) are needed to build; nothing of them is needed to run
FlowCellar, which ships the built files.
"""
from __future__ import annotations

import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
VIEWER = REPO / "viewer"
DIST = REPO / "src" / "flowcellar" / "viewer_dist"


def main(argv: list[str]) -> int:
    node = shutil.which("node")
    if not node:
        print("build_viewer: Node is not installed (needed only to build the viewer)", file=sys.stderr)
        return 2
    if not (VIEWER / "node_modules" / "esbuild").is_dir():
        print("build_viewer: viewer/node_modules is missing; run 'npm ci' in viewer/ first", file=sys.stderr)
        return 2
    if "--check" in argv:
        if not (DIST / "manifest.json").exists():
            print("viewer_dist/ does not exist: run python tools/build_viewer.py")
            return 1
        have = json.loads((DIST / "manifest.json").read_text("utf-8"))
        with tempfile.TemporaryDirectory() as tmp:
            flags = [] if have.get("minified", True) else ["--dev"]
            r = subprocess.run([node, str(VIEWER / "build.mjs"), "--out", tmp, *flags], cwd=VIEWER, capture_output=True, text=True)
            if r.returncode != 0:
                print(r.stdout + r.stderr, file=sys.stderr)
                return r.returncode
            fresh = json.loads((Path(tmp) / "manifest.json").read_text("utf-8"))
        same = {k: v["sha256"] for k, v in have["files"].items()} == {k: v["sha256"] for k, v in fresh["files"].items()}
        print("viewer_dist/ is up to date" if same else "viewer_dist/ is older than viewer/: run python tools/build_viewer.py")
        return 0 if same else 1
    return subprocess.run([node, str(VIEWER / "build.mjs"), *[a for a in argv if a in ("--dev",)]], cwd=VIEWER).returncode


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
