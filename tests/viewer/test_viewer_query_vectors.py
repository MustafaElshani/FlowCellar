"""The Runs table's search grammar (viewer/src/hub/query.js) against the shared vectors of tests/vectors/query.

The same vectors are run by the Python side (``flowcellar ls --where``); together the two tests hold the two
implementations of the one grammar equal. Needs Node and viewer/node_modules; skipped where they are missing.
"""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
pytestmark = pytest.mark.viewer


def test_the_viewer_answers_every_shared_query_vector():
    node = shutil.which("node")
    if not node or not (REPO / "viewer" / "node_modules" / "esbuild").is_dir():
        pytest.skip("Node or viewer/node_modules is missing (run 'npm ci' in viewer/)")
    if not list((REPO / "tests" / "vectors" / "query").glob("*.json")):
        pytest.skip("no query vectors yet")
    r = subprocess.run([node, str(REPO / "viewer" / "verify" / "query_vectors.mjs")], cwd=REPO / "viewer", capture_output=True, text=True, timeout=120)
    assert r.returncode == 0, (r.stdout + r.stderr)[-3000:]
