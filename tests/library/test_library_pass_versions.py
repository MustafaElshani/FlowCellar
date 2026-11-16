"""A new version of a pass, and ``scan --read``, with the real engine and the real read-data passes.

What the library showed before an upgrade must come back: a file read by an earlier version of the pass is read
again (once), before anything that was never read, and only the kind of file whose pass changed. Synthetic data only.
"""
from __future__ import annotations

import os

import pytest
from library_helpers import log_events, run_scan

from flowcellar.library import scan as scan_mod
from flowcellar.model import bundle as bundle_model
from flowcellar.model import schema
from flowcellar.readdata import link as link_mod
from flowcellar.testing import libprep as F
from flowcellar.testing import readdata as T
from flowcellar.testing.fixtures import RunSpec, make_run

CONFIG = """\
[library]
name = "pass version demo"

[discovery]
min_age_seconds = 0

[[inputs]]
role = "align.genome_bam"
glob = "**/analysis/alignments/{sample}_genome_sorted.bam"
[[inputs]]
role = "polya.read_table"
glob = "**/analysis/dorado_polyA/annotated_readid/{sample}_readid_annotated.csv"
[[inputs]]
role = "basecaller.sample_sheet"
sniff = "csv_columns:barcode,alias,flow_cell_id|position_id"

[[link.rules]]
name = "experiment-folder"
roles = ["polya.read_table"]
via = "sample"
within = "experiment"
"""


@pytest.fixture
def project(tmp_path, monkeypatch):
    pytest.importorskip("pysam", reason="the fixture project writes BAM files (extra 'bam')")
    monkeypatch.setenv(link_mod.INPROCESS_ENV, "1")
    root = tmp_path / "store"
    fr = make_run(root / "project_one", RunSpec(experiment="demo_experiment", sample="pool_1", seed=31))
    proj = T.make_project(str(root / "project_one"), T.run_ref(fr), seed=33, quant=False, dorado=False, reads_per_sample=120)
    for i, alias in enumerate(proj.aliases):
        path = os.path.join(proj.project_dir, "analysis", "dorado_polyA", "annotated_readid", f"{alias}_readid_annotated.csv")
        F.write_annotated_table(path, F.make_table_rows(300 + 20 * i, seed=40 + i))
    lib = tmp_path / "library"
    lib.mkdir()
    (lib / "flowcellar.toml").write_text(CONFIG)
    return {"lib": str(lib), "root": str(root), "run_key": fr.run_key, "aliases": list(proj.aliases)}


def _reads(res) -> list[tuple[str, str]]:
    """(role, file name) of every pass a scan really made, in the order it made them."""
    done = {e["path"]: e for e in log_events(res, "read_pass.done")}
    return [(e["role"], os.path.basename(e["path"])) for e in log_events(res, "read_pass.start") if not done[e["path"]]["cached"]]


def _partials(lib: str) -> set[str]:
    return {n for _, _, names in os.walk(os.path.join(lib, "partials")) for n in names if n.endswith(".fcp.json.gz")}


def _sections(lib: str, run_key: str) -> dict:
    b = bundle_model.read_bundle(os.path.join(lib, "bundles", bundle_model.bundle_filename(run_key)))
    assert schema.schema_errors(b) == []
    return b["modules"]


def test_a_new_pass_version_reads_its_files_again_once_and_only_those(project, monkeypatch):
    lib, n = project["lib"], len(project["aliases"])
    first = run_scan(lib, project["root"], no_hub=True)
    assert first.exit_code == 0
    assert sorted(r for r, _ in _reads(first)) == ["align.genome_bam"] * n + ["polya.read_table"] * n
    assert _sections(lib, project["run_key"])["aln"]["state"] == "ok"
    before = _partials(lib)
    assert not _reads(run_scan(lib, no_hub=True)), "a second scan reads nothing"

    monkeypatch.setitem(link_mod.PASS_MODULES, "aln", "9.9.9")        # the BAM pass of a later version
    again = run_scan(lib, no_hub=True)
    assert again.exit_code == 0
    assert sorted(_reads(again)) == sorted(("align.genome_bam", f"{a}_genome_sorted.bam") for a in project["aliases"]), \
        "every BAM is read again, once; the tables are not, because their pass did not change"
    assert len(log_events(again, "scan.read_again")) == 1
    after = _partials(lib)
    assert len(after) == len(before) and len(after - before) == n, "the partials of the earlier pass are gone, not kept beside the new ones"
    assert _sections(lib, project["run_key"])["aln"]["state"] == "ok"
    assert again.waiting_read_data == 0
    assert not _reads(run_scan(lib, no_hub=True)), "and the scan after that reads nothing"


def test_files_read_by_an_earlier_pass_come_before_files_never_read(project, monkeypatch):
    lib, aliases = project["lib"], project["aliases"]
    clock = [0.0]
    monkeypatch.setattr(scan_mod, "_monotonic", lambda: clock[0])
    real = scan_mod.Scan._read_one

    def slow(self, li, known):                                        # every pass takes 25 minutes
        real(self, li, known)
        clock[0] += 25 * 60
    monkeypatch.setattr(scan_mod.Scan, "_read_one", slow)
    first = run_scan(lib, project["root"], no_hub=True, read=("bam",), time_box_minutes=40)
    read_first = [name for _, name in _reads(first)]
    assert len(read_first) == 2 and first.waiting_read_data == len(aliases) * 2 - 2, "two passes start inside 40 minutes"

    monkeypatch.setitem(link_mod.PASS_MODULES, "aln", "9.9.9")
    clock[0] = 0.0
    second = run_scan(lib, no_hub=True, read=("bam",), time_box_minutes=0)
    names = [name for _, name in _reads(second)]
    assert sorted(names[:2]) == sorted(read_first), "what was read before is read again first"
    assert len(names) == len(aliases) == len(set(names)), "then the rest, and no file twice"
    waiting = _sections(lib, project["run_key"])
    assert waiting["aln"]["state"] == "ok"


def test_read_chooses_the_kind_of_file_a_scan_reads(project):
    lib, n = project["lib"], len(project["aliases"])
    tables = run_scan(lib, project["root"], no_hub=True, read=("tables",))
    assert {r for r, _ in _reads(tables)} == {"polya.read_table"} and len(_reads(tables)) == n
    assert tables.waiting_read_data == n, "the BAMs wait for a scan that reads them"
    mods = _sections(lib, project["run_key"])
    assert mods["rna.biotype"]["state"] == "ok"
    bams = run_scan(lib, no_hub=True, read=("bam",))
    assert {r for r, _ in _reads(bams)} == {"align.genome_bam"} and bams.waiting_read_data == 0
    assert _sections(lib, project["run_key"])["aln"]["state"] == "ok"


def test_new_section_or_linking_rules_put_the_bundle_together_again_without_reading(project, monkeypatch):
    """A later version that attaches files differently, or builds the sections differently from the same partials:
    the next scan links again and assembles again, and reads no file for it."""
    from flowcellar.readdata import sections as sections_mod
    lib = project["lib"]
    assert run_scan(lib, project["root"], no_hub=True).exit_code == 0
    path = os.path.join(lib, "bundles", bundle_model.bundle_filename(project["run_key"]))
    before = os.stat(path).st_mtime_ns
    same = run_scan(lib, no_hub=True)
    assert not _reads(same) and os.stat(path).st_mtime_ns == before and not log_events(same, "link"), "nothing changed: nothing is linked or written again"

    monkeypatch.setattr(link_mod, "LINK_VERSION", link_mod.LINK_VERSION + 1)
    relinked = run_scan(lib, no_hub=True)
    assert relinked.exit_code == 0 and len(log_events(relinked, "link")) == 1 and not _reads(relinked), "linked again, from the cached probes"

    monkeypatch.setattr(sections_mod, "VERSION", "9.9.9")
    again = run_scan(lib, no_hub=True)
    assert again.exit_code == 0 and not _reads(again), "no pass is made for a new way of putting sections together"
    assert os.stat(path).st_mtime_ns != before, "the bundle was assembled again"
    assert _sections(lib, project["run_key"])["aln"]["version"] == "9.9.9"
    settled = os.stat(path).st_mtime_ns
    assert not _reads(run_scan(lib, no_hub=True)) and os.stat(path).st_mtime_ns == settled, "and the scan after that leaves it alone"
