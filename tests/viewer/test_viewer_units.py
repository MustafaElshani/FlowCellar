"""The viewer's plain functions, checked without a browser (viewer/verify/units.mjs): the unit of a number follows
the end of its id (every catalogue number in hours is written in hours), a time is shown on the clock of the
computer the page is read on with its zone said, and the columns of the Runs table grow with the text."""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
pytestmark = pytest.mark.viewer


def test_units_times_and_columns():
    node = shutil.which("node")
    if not node or not (REPO / "viewer" / "node_modules" / "esbuild").is_dir():
        pytest.skip("Node or viewer/node_modules is missing (run 'npm ci' in viewer/)")
    r = subprocess.run([node, str(REPO / "viewer" / "verify" / "units.mjs")], cwd=REPO / "viewer", capture_output=True, text=True, timeout=180)
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-2000:]
    assert "all checks passed in 5 time zones" in r.stdout
