"""Fixtures of the viewer's tests: the synthetic library (invented names only) and its read-only stand-in."""
from __future__ import annotations

import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import viewer_synth  # noqa: E402
from viewer_fake_library import FakeLibrary  # noqa: E402

from flowcellar.util.errors import InputError  # noqa: E402


@pytest.fixture(scope="session")
def viewer_built():
    """The built viewer. Tests that need it are skipped, with the command to run, when it is not there."""
    from flowcellar.report import html as page

    try:
        return page.load_viewer()
    except InputError as e:
        pytest.skip(f"{e.message} ({e.hint})")


@pytest.fixture(scope="session")
def synth(tmp_path_factory):
    """The synthetic library on disk: ``{"path", "specs", "bundles", "index"}``. Read-only for the tests."""
    return viewer_synth.make_library(str(tmp_path_factory.mktemp("synth_library") / "library"))


@pytest.fixture()
def library(tmp_path):
    """A fresh copy of the synthetic library that a test may write a hub into and change."""
    made = viewer_synth.make_library(str(tmp_path / "library"))
    return FakeLibrary(made["path"]), made
