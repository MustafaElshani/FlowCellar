"""Fixtures of the read-data tests. Everything is synthetic (flowcellar.testing.readdata)."""
from __future__ import annotations

import os

import pytest

pytest.importorskip("pysam", reason="the read-data tests write and read BAM files (extra 'bam')")

from flowcellar.readdata import api, link  # noqa: E402
from flowcellar.testing import readdata as T  # noqa: E402
from flowcellar.testing.fixtures import RunSpec, make_run  # noqa: E402


@pytest.fixture(autouse=True)
def _inprocess(request, monkeypatch):
    """Read in the test process, unless the test has 'pool' in its name: those go through worker processes."""
    if "pool" in request.node.name:
        monkeypatch.delenv(link.INPROCESS_ENV, raising=False)
    else:
        monkeypatch.setenv(link.INPROCESS_ENV, "1")


from rd_helpers import bundle_of, do_passes, ident, project_candidates, tier1_of  # noqa: E402


@pytest.fixture(scope="session")
def project(tmp_path_factory):
    """One fake run folder with one fake project beside it: (FixtureRun, ProjectFixture)."""
    root = tmp_path_factory.mktemp("rd_project")
    fr = make_run(root / "project_one", RunSpec(experiment="demo_experiment", sample="pool_1", seed=21))
    proj = T.make_project(str(root / "project_one"), T.run_ref(fr), seed=31)
    return fr, proj


@pytest.fixture(scope="session")
def project_data(project, tmp_path_factory):
    """The fixture project linked, read and assembled once: a dict with everything the section tests look at."""
    before = os.environ.get(link.INPROCESS_ENV)
    os.environ[link.INPROCESS_ENV] = "1"                 # a session fixture is made outside the per-test setting
    try:
        fr, proj = project
        lib = tmp_path_factory.mktemp("rd_library")
        run = ident(fr)
        linked = api.link(project_candidates(proj), [run], config={}, probes_dir=str(lib / "probes"))
        refs = do_passes(linked, lib / "partials", [run.run_key])
        tier1 = tier1_of(fr)
        rd = api.read_data(run, linked, str(lib / "partials"), tier1=tier1, config={})
    finally:
        if before is None:
            os.environ.pop(link.INPROCESS_ENV, None)
        else:
            os.environ[link.INPROCESS_ENV] = before
    return {"fr": fr, "proj": proj, "run": run, "linked": linked, "refs": refs, "tier1": tier1, "rd": rd, "lib": lib,
            "bundle": bundle_of(tier1, rd)}
