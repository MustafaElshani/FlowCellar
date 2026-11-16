"""The browser checks (viewer/verify/fixture.mjs) as one pytest test.

They open the synthetic library in headless Firefox: every screen in both themes at three widths, from file://
and over HTTP, with no page error, no console error, no policy violation and no request to another origin.
They take a few minutes, so they run only when asked for::

    FLOWCELLAR_BROWSER=1 .venv/bin/python -m pytest tests/viewer/test_viewer_browser.py -q

A Playwright Firefox must be on this machine (PLAYWRIGHT_BROWSERS_PATH or FLOWCELLAR_PW_BROWSERS names its folder).
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
pytestmark = [pytest.mark.viewer, pytest.mark.slow]


@pytest.mark.skipif(os.environ.get("FLOWCELLAR_BROWSER") != "1", reason="set FLOWCELLAR_BROWSER=1 to run the browser checks")
def test_every_screen_in_a_real_browser(tmp_path):
    node = shutil.which("node")
    if not node or not (REPO / "viewer" / "node_modules" / "playwright").is_dir():
        pytest.skip("Node or viewer/node_modules/playwright is missing")
    env = dict(os.environ, FLOWCELLAR_PYTHON=sys.executable)
    r = subprocess.run([node, str(REPO / "viewer" / "verify" / "fixture.mjs"), "--out", str(tmp_path / "verify")], cwd=REPO / "viewer",
                       capture_output=True, text=True, timeout=1800, env=env)
    if r.returncode == 2 and "No Playwright Firefox found" in r.stderr:
        pytest.skip("no Playwright Firefox on this machine")
    tail = "\n".join((r.stdout + r.stderr).splitlines()[-40:])
    assert r.returncode == 0, tail
