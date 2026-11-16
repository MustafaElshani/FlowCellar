"""Shared pytest fixtures. Everything here is synthetic: no real run, sample or project name may appear in tests."""
from __future__ import annotations

import pytest

from flowcellar.testing.fixtures import RunSpec, make_library_tree, make_run


@pytest.fixture(scope="session")
def fake_run(tmp_path_factory):
    """One fake MinION-like run folder: 8 h, 4,000 reads, basecalling off, newest file layout."""
    return make_run(tmp_path_factory.mktemp("fake_run"), RunSpec())


@pytest.fixture(scope="session")
def fake_tree(tmp_path_factory):
    """Two roots, three runs on two flow cells (see flowcellar.testing.fixtures.make_library_tree)."""
    root = tmp_path_factory.mktemp("fake_tree")
    return root, make_library_tree(root)
