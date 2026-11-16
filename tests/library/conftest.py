"""Fixtures of the library tests: the three stand-ins and the fixture trees. Synthetic data only."""
from __future__ import annotations

import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(__file__))

from library_fakes import FakeEngine, FakeReadData, FakeReport  # noqa: E402
from library_helpers import run_scan, set_root  # noqa: E402

from flowcellar.testing import trees  # noqa: E402


@pytest.fixture(autouse=True)
def _quiet_logger():
    """Each test starts with the package logger in its untouched state (a test may point it at a captured stream)."""
    import logging
    lg = logging.getLogger("flowcellar")
    lg.handlers[:] = []
    lg.propagate = True
    lg.setLevel(logging.INFO)
    yield
    lg.handlers[:] = []
    lg.propagate = True


@pytest.fixture
def engine(monkeypatch):
    return FakeEngine().install(monkeypatch)


@pytest.fixture
def readdata(monkeypatch):
    return FakeReadData().install(monkeypatch)


@pytest.fixture
def report(monkeypatch):
    return FakeReport().install(monkeypatch)


@pytest.fixture
def parts(engine, readdata, report):
    """All three stand-ins installed: ``parts.engine``, ``parts.readdata``, ``parts.report``."""
    class Parts:
        pass
    p = Parts()
    p.engine, p.readdata, p.report = engine, readdata, report
    return p


@pytest.fixture
def tree(tmp_path):
    """Three roots, five runs with reports, one raw-only folder; alpha_project lies under two roots."""
    return trees.make_three_roots(tmp_path / "data")


@pytest.fixture
def lib_dir(tmp_path):
    return str(tmp_path / "library")


@pytest.fixture
def scanned(parts, tree, lib_dir):
    """A library made by one scan of the three roots of ``tree`` (store: project = second folder)."""
    res = run_scan(lib_dir, *tree.root_list())
    set_root(lib_dir, "root_store", project_depth=2)
    res = run_scan(lib_dir)
    return res
