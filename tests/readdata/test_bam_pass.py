"""The BAM source: exact counts against what the generator wrote, whatever the number of workers, and every hazard
of real files on a fake one."""
from __future__ import annotations

import os
import shutil

import numpy as np
import pytest

from flowcellar.aggregates import build
from flowcellar.model import axes
from flowcellar.readdata import bam, bamindex
from flowcellar.testing import readdata as T

A = T.fake_run_ref("TST00001", 1)
B = T.fake_run_ref("TST00002", 2)


def scan(path, **kw):
    kw.setdefault("inprocess", True)
    return bam.scan_bam(str(path), **kw)


def assert_groups(merged, expected):
    """Every aggregate of every (run key, sample label) equals what the generator counted."""
    assert set(merged.groups) == set(expected["groups"])
    for gid, acc in merged.groups.items():
        x = expected["groups"][gid]
        got = dict(reads=acc.reads, bases=acc.bases, mapped=acc.mapped, unmapped=acc.unmapped, secondary=acc.secondary,
                   supplementary=acc.supplementary, nm=acc.nm_sum, cols=acc.aln_cols, perfect=acc.perfect, with_sa=acc.with_sa,
                   q_missing=acc.q_missing, split_children=acc.split_children)
        assert got == {k: x[k] for k in got}, gid
        parents = merged.parents[gid]
        assert (len(parents) if parents is not None else None) == x["split_parents"]
        # the histograms, bin by bin, from the generator's own lists
        lengths = np.asarray(x["lengths"], dtype=np.float64)
        want = np.bincount(axes.bin_index(lengths, "read_length_log40"), minlength=axes.n_cells("read_length_log40"))
        assert (acc.h_len == want).all()
        assert int(acc.h_len_w.sum()) == x["bases"] == int(lengths.sum())
        quals = np.asarray(x["quals"], dtype=np.float64)
        assert (acc.h_q == np.bincount(axes.bin_index(quals, "qscore_0p25"), minlength=axes.n_cells("qscore_0p25"))).all()
        if x["mapq"]:
            assert (acc.h_mapq == np.bincount(axes.bin_index(np.asarray(x["mapq"], float), "mapq_1"), minlength=axes.n_cells("mapq_1"))).all()
            idq = np.asarray(x["identity_q"], dtype=np.float64)
            assert (acc.h_id == np.bincount(axes.bin_index(idq, "identity_q_0p1"), minlength=axes.n_cells("identity_q_0p1"))).all()
            assert int(acc.h_id[-1]) == x["perfect"]            # identity 1 is the overflow cell
        if acc.g_n is not None and len(lengths) == len(quals):
            assert int(acc.g_n.sum()) == len(quals) and int(acc.g_w.sum()) == x["bases"]


@pytest.fixture(scope="module")
def files(tmp_path_factory):
    d = tmp_path_factory.mktemp("bams")
    ra = T.make_reads(1500, A, seed=1, alias="alpha", barcode="barcode01", split_every=50)
    rb = T.make_reads(1300, B, seed=2, alias="alpha", barcode="barcode01")
    out = {"dir": d, "reads_a": ra, "reads_b": rb}
    out["dorado"] = (d / "dorado.bam", T.write_dorado_bam(str(d / "dorado.bam"), ra))
    out["aligned"] = (d / "two_cells_genome_sorted.bam", T.write_aligned_bam(str(d / "two_cells_genome_sorted.bam"), ra + rb, seed=5))
    return out


def test_dorado_bam_counts_are_exact(files):
    path, exp = files["dorado"]
    r = scan(path)
    m = r["merged"]
    assert r["state"] == "ok" and r["coverage"] == {"kind": "full", "records": 1500, "records_total": None, "fraction": 1.0}
    assert m.records == exp["records"] == m.primary == 1500
    assert_groups(m, exp)
    (gid, acc), = m.groups.items()
    assert gid == (A.run_key, "alpha") and acc.split_children == 60 and len(m.parents[gid]) == 30
    f = r["plan"]["facts"]
    assert f["dorado"]["min_qscore"] == 7.0 and f["dorado"]["no_trim"] and f["dorado"]["barcode_both_ends"] and f["dorado"]["demux"]
    assert f["runs"] == {A.run_key: 2} and f["models"] == [T.MODEL]


def test_aligned_bam_counts_are_exact_and_primary_records_count_reads_once(files):
    path, exp = files["aligned"]
    r = scan(path)
    m = r["merged"]
    assert r["state"] == "ok"
    assert m.records == exp["records"] and m.primary == exp["primary"] == 2800
    assert exp["secondary"] > 500 and exp["supplementary"] > 20            # the file really holds what the test is about
    assert sum(a.secondary for a in m.groups.values()) == exp["secondary"]
    assert sum(a.supplementary for a in m.groups.values()) == exp["supplementary"]
    assert m.primary + exp["secondary"] + exp["supplementary"] == m.records
    assert_groups(m, exp)
    # the index totals are reproduced: a complete pass finds every record the index counts
    idx = r["plan"]["index"]
    assert idx["records"] == m.records == r["coverage"]["records_total"]
    assert idx["no_coordinate"] == sum(a.unmapped for a in m.groups.values())


def test_a_file_with_reads_of_two_flow_cells_is_split_per_read(files):
    path, exp = files["aligned"]
    m = scan(path)["merged"]
    by_run = {key: acc.reads for (key, _), acc in m.groups.items()}
    assert by_run == {A.run_key: 1500, B.run_key: 1300}
    assert sum(by_run.values()) == m.primary                              # routed reads add up to the primary records
    assert {k: g["bases"] for (k, _), g in exp["groups"].items()} == {k: a.bases for (k, _), a in m.groups.items()}


@pytest.mark.parametrize("task_bytes", [1 << 30, 60_000, 9_000])
def test_result_does_not_depend_on_how_the_file_is_cut(files, task_bytes):
    path, exp = files["aligned"]
    r = scan(path, task_bytes=task_bytes)
    assert r["state"] == "ok" and r["errors"] == []
    if task_bytes < 1 << 20:
        assert r["timings"]["tasks"] > 5
    assert_groups(r["merged"], exp)
    assert r["merged"].records == exp["records"]


@pytest.mark.slow
@pytest.mark.parametrize("jobs", [1, 2, 5])
def test_worker_pool_gives_identical_partials(files, jobs, tmp_path):
    """The same file through 1, 2 and 5 worker processes: equal, cell for cell."""
    path, exp = files["aligned"]
    r = bam.scan_bam(str(path), jobs=jobs, task_bytes=20_000, inprocess=False, stall_timeout_s=60)
    assert r["state"] == "ok"
    assert r["timings"]["workers"] == jobs and r["timings"]["tasks"] > jobs
    assert_groups(r["merged"], exp)
    cells = {gid: acc.cells(aligned=True) for gid, acc in r["merged"].groups.items()}
    ref = {gid: acc.cells(aligned=True) for gid, acc in scan(path)["merged"].groups.items()}
    assert cells == ref


def test_unindexed_and_indexed_passes_agree(files, tmp_path):
    path, exp = files["aligned"]
    plain = tmp_path / "no_index.bam"
    shutil.copy(path, plain)
    a = scan(plain)
    b = scan(path, task_bytes=15_000)
    assert a["plan"]["index"] is None and a["timings"]["tasks"] == 1 and b["timings"]["tasks"] > 1
    assert a["coverage"]["records_total"] is None and b["coverage"]["records_total"] == exp["records"]
    assert {g: x.cells(aligned=True) for g, x in a["merged"].groups.items()} == {g: x.cells(aligned=True) for g, x in b["merged"].groups.items()}


def test_csi_index_is_read(tmp_path):
    reads = T.make_reads(900, A, seed=4, alias="alpha", barcode="barcode01")
    p = tmp_path / "csi.bam"
    exp = T.write_aligned_bam(str(p), reads, index="csi")
    idx = bamindex.read_index(str(p) + ".csi")
    assert idx.kind == "csi" and idx.records == exp["records"] and len(idx.offsets) > 3
    r = scan(p, task_bytes=10_000)
    assert r["timings"]["tasks"] > 1 and r["merged"].records == exp["records"]
    assert_groups(r["merged"], exp)


def test_a_stale_index_is_not_used(tmp_path):
    reads = T.make_reads(400, A, seed=4, alias="alpha", barcode="barcode01")
    p = tmp_path / "stale.bam"
    exp = T.write_aligned_bam(str(p), reads)
    old = os.stat(str(p) + ".bai").st_mtime - 3600
    os.utime(str(p) + ".bai", (old, old))
    r = scan(p, task_bytes=5_000)
    assert r["timings"]["tasks"] == 1 and r["state"] == "ok"
    assert [e["code"] for e in r["errors"]] == ["index_stale"]
    assert_groups(r["merged"], exp)


def test_an_index_of_another_file_is_noticed_and_the_file_read_again(tmp_path):
    a = tmp_path / "a.bam"
    b = tmp_path / "b.bam"
    T.write_aligned_bam(str(a), T.make_reads(700, A, seed=4, alias="alpha"), seed=1)
    exp = T.write_aligned_bam(str(b), T.make_reads(500, A, seed=5, alias="alpha"), seed=2, index=None)
    shutil.copy(str(a) + ".bai", str(b) + ".bai")
    r = scan(b, task_bytes=5_000)
    assert "index_mismatch" in [e["code"] for e in r["errors"]] or r["state"] == "truncated"
    if r["state"] == "ok":
        assert_groups(r["merged"], exp)


@pytest.mark.parametrize("wrap", ["BC", "qs", "mixed", None])
def test_wrapped_first_tag_is_unwrapped(tmp_path, wrap):
    """minimap2 -y after samtools fastq -T: the first tag arrives as CO:Z:<tag>:<type>:<value>."""
    reads = T.make_reads(300, A, seed=6, alias="alpha_one", barcode="barcode03")
    p = tmp_path / f"wrap_{wrap}.bam"
    exp = T.write_aligned_bam(str(p), reads, wrap=wrap)
    m = scan(p)["merged"]
    assert set(m.groups) == {(A.run_key, "alpha_one")}                     # the alias, wrapped or not, never the barcode of the read group
    assert_groups(m, exp)
    assert sum(a.q_missing for a in m.groups.values()) == 0 and sum(a.q_computed for a in m.groups.values()) == 0


def test_duplicate_ts_tag_and_prefixed_names_do_not_disturb(tmp_path):
    reads = T.make_reads(250, A, seed=7, alias="alpha", barcode="barcode01")
    p = tmp_path / "ts.bam"
    exp = T.write_aligned_bam(str(p), reads, duplicate_ts=True, pychopper_names=True)
    import pysam
    with pysam.AlignmentFile(str(p)) as fh:
        tags = [t for t, _ in next(r for r in fh if not r.is_unmapped).get_tags()]
    assert tags.count("ts") == 2                                           # the fixture really carries both
    r = scan(p)
    assert_groups(r["merged"], exp)
    assert r["merged"].pychopper_names == r["merged"].names_checked == 200


def test_read_without_qs_takes_quality_from_its_bases_and_without_bases_is_missing(tmp_path):
    reads = T.make_reads(200, A, seed=8, alias="alpha")
    p1, p2 = tmp_path / "noqs.bam", tmp_path / "noqual.bam"
    e1 = T.write_aligned_bam(str(p1), reads, dorado_tags=False, wrap=None)
    e2 = T.write_aligned_bam(str(p2), reads, dorado_tags=False, wrap=None, base_qualities=False)
    m1, m2 = scan(p1)["merged"], scan(p2)["merged"]
    (g1, a1), = m1.groups.items()
    (g2, a2), = m2.groups.items()
    assert g1 == g2 == (None, None)                                        # no identifier at all: unrouted, unlabelled
    assert a1.q_computed == 200 and a1.q_missing == 0 and int(a1.h_q.sum()) == 200
    assert a2.q_missing == 200 and int(a2.h_q.sum()) == 0 and a2.g_n is None
    assert_groups(m1, e1)
    assert e2["groups"][(None, None)]["q_missing"] == 200
    cell = a2.cells(aligned=True)["qual"]
    assert cell["hists"]["qual_hist"]["missing"] == 200 and "hists2d" not in cell     # counted as missing, never binned, never zero


def test_read_group_without_header_gives_the_run_id_only(tmp_path):
    reads = T.make_reads(120, A, seed=9, alias="alpha", barcode="barcode05")
    p = tmp_path / "rgonly.bam"
    exp = T.write_dorado_bam(str(p), reads, header_rg=False, fn_tag=False, bc_tag=False)
    m = scan(p)["merged"]
    assert set(m.groups) == {("?_" + A.run8, "alpha")} == set(exp["groups"])
    assert_groups(m, exp)


def test_header_alone_routes_a_file_whose_reads_carry_no_tag(tmp_path):
    reads = T.make_reads(90, A, seed=10)
    p = tmp_path / "header_only.bam"
    exp = T.write_dorado_bam(str(p), reads, fn_tag=False, rg_tag=False, bc_tag=False)
    m = scan(p)["merged"]
    assert set(m.groups) == {(A.run_key, None)}
    assert_groups(m, exp)


def test_truncated_file_keeps_what_was_read_and_says_so(tmp_path):
    reads = T.make_reads(3000, A, seed=11, alias="alpha")
    p = tmp_path / "cut.bam"
    exp = T.write_dorado_bam(str(p), reads, truncate_frac=0.55)
    assert not bam.has_eof_marker(str(p))
    r = scan(p)
    assert r["state"] == "truncated" and r["coverage"]["kind"] == "truncated"
    assert 0 < r["merged"].primary < 3000 and exp["truncated"]
    assert any(e["code"] == "truncated" for e in r["errors"])


def test_time_limit_stops_the_pass_and_reports_coverage(files):
    import time
    path, exp = files["aligned"]
    r = scan(path, task_bytes=5_000, deadline=time.monotonic() - 1.0)
    assert r["state"] == "partially_read" and r["coverage"]["kind"] == "byte_ranges"
    assert r["coverage"]["fraction"] < 1.0 and r["merged"].records < exp["records"]


def test_probe_reads_a_bounded_part_and_names_the_runs(files):
    path, exp = files["aligned"]
    pr = bam.probe(str(path), max_primary=400, regions=8)
    assert set(pr["runs"]) == {A.run_key, B.run_key} and pr["unrouted"] == 0
    assert 0 < pr["primary"] <= 400 + 8 and pr["regions"] <= 8
    assert pr["labels"] == {"alpha": pr["primary"]}
    assert pr["index"]["records"] == exp["records"] and pr["aligner"]["input_full_length"]
    pd = bam.probe(str(files["dorado"][0]))
    assert pd["header_runs"] == {A.run_key: 2} and pd["runs"] == {A.run_key: 1500} and pd["n_sq"] == 0


def test_cut_tasks_cover_the_file_without_overlap():
    offs = np.array([(c << 16) | u for c in range(1000, 500_000, 3_000) for u in (0, 777)], dtype=np.uint64)
    first = 1000 << 16
    tasks = bamindex.cut_tasks(first, offs, 500_000, 40_000)
    assert tasks[0][0] == first and tasks[-1][1] is None
    assert all(tasks[i][1] == tasks[i + 1][0] for i in range(len(tasks) - 1))
    assert all(a < b for a, b in tasks[:-1])
    sizes = [(b >> 16) - (a >> 16) for a, b in tasks[:-1]]
    assert min(sizes) >= 40_000 and max(sizes) < 40_000 + 3_000 + 1
    assert bamindex.cut_tasks(first, offs[:0], 500_000, 40_000) == [(first, None)]


def test_stored_cells_hold_exact_totals(files):
    path, exp = files["aligned"]
    acc = scan(path)["merged"].groups[(A.run_key, "alpha")]
    x = exp["groups"][(A.run_key, "alpha")]
    c = acc.cells(aligned=True)
    h = c["len"]["hists"]["len_hist"]
    assert sum(h["n"]) + h["under"] + h["over"] == x["reads"]
    assert sum(h["w"]) + h["w_under"] + h["w_over"] == x["bases"] and h["min"] == min(x["lengths"]) and h["max"] == max(x["lengths"])
    g = c["qual"]["hists2d"]["qual_vs_length"]
    assert sum(g["n"]) + g["outside"] == x["reads"] and sum(g["w"]) + g["w_outside"] == x["bases"]
    a = c["aln"]
    assert a["scalars"]["mapped_reads"] + a["scalars"]["unmapped_reads"] == a["scalars"]["reads"] == x["reads"]
    ih = a["hists"]["aln_identity_hist"]
    assert sum(ih["n"]) + ih["under"] + ih["over"] + ih["missing"] == x["reads"] and ih["over"] == x["perfect"]
    assert (build.hist_full(a["hists"]["aln_mapq_hist"]).sum()) == x["mapped"]


def test_read_group_value_with_a_hash_suffix_or_an_alias_the_header_does_not_list():
    hd = {"RG": [{"ID": f"{A.protocol_run_id}_{T.MODEL}", "PU": A.flow_cell_id, "DS": f"runid={A.protocol_run_id} basecall_model={T.MODEL}"},
                 {"ID": f"{A.protocol_run_id}_{T.MODEL}_{T.KIT}_barcode03", "PU": A.flow_cell_id, "DS": f"runid={A.protocol_run_id} basecall_model={T.MODEL}",
                  "SM": "barcode03", "al": "barcode03"}]}
    st = bam._State(bam.header_facts(hd))
    base = f"{A.protocol_run_id}_{T.MODEL}"
    assert st._rg(base) == (A.run_key, None)
    assert st._rg(base + "-1A2B3C4D") == (A.run_key, None)                       # merged by samtools: the suffix is not a sample
    assert st._rg(base + f"_{T.KIT}_barcode03") == (A.run_key, "barcode03")
    assert st._rg(base + f"_{T.KIT}_barcode03-0F0F0F0F") == (A.run_key, "barcode03")
    assert st._rg(base + "_my_sample_name") == (A.run_key, "my_sample_name")     # renamed by a sample sheet
    assert st._rg(base + "_unclassified") == (A.run_key, None)
    # no header at all: the run id of the value, the alias after the model name
    bare = bam._State(bam.header_facts({}))
    assert bare._rg(base + "_my_sample_name") == ("?_" + A.run8, "my_sample_name")
    assert bare._rg(base + f"_{T.KIT}_barcode11-AB12CD34") == ("?_" + A.run8, "barcode11")
    assert bare._rg("not a dorado read group") == (None, None)
    assert bam.norm_label("SQK-NBD114-24_barcode07") == "barcode07" and bam.norm_label("unclassified") is None and bam.norm_label("alias_x") == "alias_x"
    assert bam.model_tier(T.MODEL) == "sup" and bam.model_tier("dna_r9.4.1_e8_hac@v3.3") == "hac" and bam.model_tier(None) is None


@pytest.mark.slow
def test_pool_that_makes_no_progress_is_abandoned():
    """A worker blocked on a file that never answers must not hold the parent: the pool is terminated and says so."""
    import time

    from flowcellar.readdata import _pool
    t0 = time.monotonic()
    out = _pool.run_tasks(time.sleep, [30, 30], workers=2, stall_timeout_s=0.5)
    assert out == {"results": [], "stalled": True} and time.monotonic() - t0 < 15
    done = _pool.run_tasks(abs, [-1, -2, -3], workers=2, stall_timeout_s=30)
    assert sorted(done["results"]) == [1, 2, 3] and not done["stalled"]
    seen = []
    assert _pool.run_tasks(abs, [-4], inprocess=True, progress=lambda d, n: seen.append((d, n)))["results"] == [4] and seen == [(1, 1)]
