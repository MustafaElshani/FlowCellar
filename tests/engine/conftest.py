"""Fixtures of the ENGINE tests. Everything is synthetic: fake run folders written by flowcellar.testing."""
from __future__ import annotations

import pytest

from flowcellar.engine import api
from flowcellar.model import bundle as B
from flowcellar.testing import minknow_variants as V
from flowcellar.testing.fixtures import RunSpec


class Built:
    """One fake run taken through the engine."""

    def __init__(self, fixture, partials_dir, *, full=True, history=(), **judge_kw):
        self.fixture = fixture
        self.run = api.identify(fixture.run_dir)
        self.refs = api.full_pass(self.run, str(partials_dir)) if full else []
        self.bundle = api.assemble(self.run, self.refs)
        self.status_file = api.judge(self.bundle, history, in_library=bool(history), **judge_kw)
        self.status = self.status_file["status"]
        self.rules = {r["id"]: r for r in self.status["rules"]}

    @property
    def composed(self) -> dict:
        return B.compose(self.bundle, status=self.status, relations=self.status_file["relations"])


@pytest.fixture(scope="session")
def partials_dir(tmp_path_factory):
    return tmp_path_factory.mktemp("partials")


@pytest.fixture(scope="session")
def build(partials_dir):
    def _build(fixture, **kw):
        return Built(fixture, partials_dir, **kw)
    return _build


@pytest.fixture(scope="session")
def full8(fake_run, build):
    """The standard 8 h fake run in full detail."""
    return build(fake_run)


@pytest.fixture(scope="session")
def quick8(fake_run, build):
    return build(fake_run, full=False)


@pytest.fixture(scope="session")
def day_run(tmp_path_factory, build):
    """A 24.2 h MinION-like cDNA run in full detail: long enough to be judged at 24 h and to teach."""
    f = V.make(tmp_path_factory.mktemp("day_run"), RunSpec(flow_cell_id="TST00030", hours=24.2, planned_hours=24.0, reads=3000,
                                                           scan_period_h=6.0, seed=30, start="2026-06-01T08:00:00Z"))
    return build(f)
