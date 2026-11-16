"""The built wheel, installed away from the source tree, carries the viewer and runs ``scan`` and ``run``.

Needs a wheel, so it runs only on request::

    python -m build --wheel
    FLOWCELLAR_WHEEL=dist/flowcellar-<version>-py3-none-any.whl python -m pytest tests/integration/test_wheel.py

The wheel is installed alone into a new virtual environment; its dependencies (numpy, pyarrow, click, pysam) are
taken from the environment the tests run in, so nothing is downloaded.
"""
from __future__ import annotations

import os
import site
import subprocess
import sys
import venv

import pytest

WHEEL = os.environ.get("FLOWCELLAR_WHEEL")
pytestmark = [pytest.mark.slow, pytest.mark.skipif(not WHEEL, reason="set FLOWCELLAR_WHEEL to a built wheel to run this")]


def _run(cmd: list[str], cwd: str) -> subprocess.CompletedProcess:
    env = {k: v for k, v in os.environ.items() if not k.startswith(("FLOWCELLAR_", "PYTHON"))}
    return subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, timeout=600, env=env)


def test_the_wheel_runs_in_an_environment_of_its_own(tmp_path):
    wheel = os.path.abspath(WHEEL)
    assert os.path.isfile(wheel), wheel
    env_dir = tmp_path / "env"
    venv.EnvBuilder(with_pip=True).create(env_dir)
    py = str(env_dir / "bin" / "python")
    where = _run([py, "-c", "import site; print(site.getsitepackages()[0])"], str(tmp_path)).stdout.strip()
    # the dependencies, not the source tree: a path named in a .pth file is added as it is, so the editable
    # install of the development environment does not come along
    with open(os.path.join(where, "zz_dependencies.pth"), "w") as fh:
        fh.write("\n".join(site.getsitepackages()) + "\n")
    out = _run([py, "-m", "pip", "install", "--no-deps", "--no-index", "--quiet", wheel], str(tmp_path))
    assert out.returncode == 0, out.stderr
    work = tmp_path / "work"
    work.mkdir()
    flowcellar = str(env_dir / "bin" / "flowcellar")
    out = _run([py, "-c", "import flowcellar; from flowcellar.report import api; print(flowcellar.__file__); "
                          "print(sorted(api.viewer_manifest()['files']))"], str(work))
    assert out.returncode == 0, out.stderr
    assert str(env_dir) in out.stdout.splitlines()[0], "the package is imported from the wheel, not from the source tree"
    assert "viewer.core.js" in out.stdout and "viewer.css" in out.stdout and "template.html" in out.stdout
    assert _run([py, "-m", "flowcellar.testing.fixtures", "tree", "--tree"], str(work)).returncode == 0
    assert _run([py, "-c", "from flowcellar.testing.trees import age_tree; age_tree('tree')"], str(work)).returncode == 0
    out = _run([flowcellar, "scan", "tree", "-L", "library"], str(work))
    assert out.returncode == 0 and "done: 3 runs" in out.stderr, out.stderr
    for name in ("index.html", "index.js"):
        assert os.path.getsize(work / "library" / "hub" / name) > 0
    assert len(os.listdir(work / "library" / "hub" / "data")) == 3
    assert _run([flowcellar, "validate", "library"], str(work)).returncode == 0
    run_dir = next(os.path.dirname(os.path.join(d, n)) for d, _, names in os.walk(work / "tree") for n in names
                   if n.startswith("final_summary_"))
    out = _run([flowcellar, "run", run_dir, "-o", "qc"], str(work))
    assert out.returncode == 0, out.stderr
    made = sorted(os.listdir(work / "qc"))
    assert [n.rsplit(".", 2)[-2:] for n in made if not n.startswith(".")] == [["json", "gz"], ["flowcellar", "html"]], made
    assert sys.version_info >= (3, 11)
