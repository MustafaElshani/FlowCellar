"""Fixtures of the single-cell and spatial add-on. Everything is synthetic (``flowcellar.testing.sc``)."""
from __future__ import annotations

import pytest

from flowcellar.model import bundle as bundle_model

h5py = pytest.importorskip("h5py", reason="the single-cell and spatial add-on reads HDF5 with h5py (the 'sc' extra)")
pysam = pytest.importorskip("pysam", reason="the fixtures write BAM files with pysam (the 'bam' extra)")

from flowcellar.sc import collect, detect, judge, store  # noqa: E402
from flowcellar.testing import sc as fx  # noqa: E402


def collected(fake, **options) -> dict:
    """The measured bundle of a fake analysis folder, with its content hash."""
    claims = detect.claims_at(fake.path)
    assert len(claims) == 1, claims
    b = collect.collect(claims[0], collect.Options(workers=2, **options))
    b["content_sha256"] = bundle_model.content_sha256(b)
    return b


def composed(bundle: dict, **kw) -> dict:
    return store.compose(bundle, judge.judge(bundle, **kw))


def summary(bundle: dict) -> dict:
    return bundle["modules"]["sc"]["summary"]


def cell(bundle: dict) -> dict:
    return bundle["modules"]["sc"]["cells"][0]


@pytest.fixture(scope="session")
def sr_visium(tmp_path_factory):
    fake = fx.make_spaceranger_visium(str(tmp_path_factory.mktemp("sr_visium")), bam=True, antisense=0.2)
    return fake, collected(fake)


@pytest.fixture(scope="session")
def sr_hd(tmp_path_factory):
    fake = fx.make_spaceranger_hd(str(tmp_path_factory.mktemp("sr_hd")))
    return fake, collected(fake)


@pytest.fixture(scope="session")
def wf_v2(tmp_path_factory):
    fake = fx.make_wf_v2(str(tmp_path_factory.mktemp("wf_v2")), bam=True)
    return fake, collected(fake)


@pytest.fixture(scope="session")
def wf_v3(tmp_path_factory):
    fake = fx.make_wf_v3(str(tmp_path_factory.mktemp("wf_v3")))
    return fake, collected(fake)


@pytest.fixture(scope="session")
def wf_droplet(tmp_path_factory):
    fake = fx.make_wf_v3(str(tmp_path_factory.mktemp("wf_droplet")), kit="3prime:v3", config_stats=True, expected_cells=12,
                         spaceranger_bam=None, alias="demo_cells")
    return fake, collected(fake)
