"""The second part of the BAM pass: strand, CIGAR totals, grids, depth along a genome, transcripts, start times,
measured speed and poly(A). Every number is compared with what the generator of the fake file wrote down while it
wrote the records, or with a count done here in the plainest possible way (one array cell per reference base)."""
from __future__ import annotations

import math
from array import array
from fractions import Fraction

import numpy as np
import pytest

from flowcellar.aggregates import build, derive
from flowcellar.model import axes
from flowcellar.readdata import alnacc, bam, reference, saturation
from flowcellar.readdata.alnacc import AlignAcc, TagAcc, TimeTable
from flowcellar.readdata.sections_align import table_quantiles
from flowcellar.testing import readdata as T

A = T.fake_run_ref("TST00001", 1)
B = T.fake_run_ref("TST00002", 2)
GID = (A.run_key, "alpha")


def scan(path, **kw):
    kw.setdefault("inprocess", True)
    return bam.scan_bam(str(path), **kw)


def bincount(values, axis):
    return np.bincount(axes.bin_index(np.asarray(values, dtype=np.float64), axis), minlength=axes.n_cells(axis))


def grid(xs, x_axis, ys, y_axis):
    nx, ny = axes.n_cells(x_axis), axes.n_cells(y_axis)
    flat = axes.bin_index(np.asarray(xs, float), x_axis) * ny + axes.bin_index(np.asarray(ys, float), y_axis)
    return np.bincount(flat, minlength=nx * ny).reshape(nx, ny)


def table(pairs, width):
    """{time bin: counts over ``width`` columns} from (time bin, column) pairs."""
    out: dict[int, np.ndarray] = {}
    for t, v in pairs:
        out.setdefault(int(t), np.zeros(width, np.int64))[int(v)] += 1
    return out


def as_dict(tt: TimeTable):
    lo, a = tt.trimmed()
    return {} if lo is None else {lo + i: row for i, row in enumerate(a) if row.any()}


def same_tables(got: dict, want: dict):
    assert set(got) == set(want)
    for k in want:
        assert (got[k] == want[k]).all(), k


def depth_by_base(aligned, contigs):
    """Aligned columns per reference base, counted one base at a time: the independent count."""
    per = [np.zeros(ln + 1, np.int64) for _, ln in contigs]
    for a in aligned:
        for s, e in a["blocks"]:
            per[a["contig"]][s] += 1
            per[a["contig"]][e] -= 1
    return [np.cumsum(p)[:-1] for p in per]


def window_sums(per_base, window):
    return np.concatenate([np.add.reduceat(p, np.arange(0, len(p), window)) for p in per_base])


@pytest.fixture(scope="module")
def genome(tmp_path_factory):
    d = tmp_path_factory.mktemp("align")
    ra = T.make_reads(1500, A, seed=1, alias="alpha", barcode="barcode01", split_every=50)
    rb = T.make_reads(1300, B, seed=2, alias="alpha", barcode="barcode01")
    path = d / "two_cells_genome_sorted.bam"
    exp = T.write_aligned_bam(str(path), ra + rb, seed=5)
    return {"dir": d, "path": path, "exp": exp, "reads": ra + rb, "scan": scan(path, task_bytes=30_000)}


@pytest.fixture(scope="module")
def txome(tmp_path_factory):
    d = tmp_path_factory.mktemp("tx")
    reads = T.make_reads(2500, A, seed=3, alias="alpha", barcode="barcode01")
    path = d / "alpha_txome_sorted.bam"
    exp = T.write_transcriptome_bam(str(path), reads, seed=9)
    return {"dir": d, "path": path, "exp": exp, "reads": reads, "scan": scan(path, task_bytes=30_000)}


# ------------------------------------------------------------------------------------------- the header
def test_genome_or_transcriptome_is_told_from_the_header_alone():
    names = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
    lengths = [100_000_000 - i * 1_000_000 for i in range(24)] + [16_569]
    lay = reference.Layout(lambda i: names[i], lengths)
    assert lay.kind == "genome" and lay.names == names and lay.mito == 24 and lay.describe()["mitochondrial"] == "chrM"
    # a reference of transcripts by its names, whatever its size
    tx = T.fake_transcripts(30)
    lt = reference.Layout(lambda i: tx[i][0], [ln for _, ln in tx])
    assert lt.kind == "transcriptome" and "transcript identifiers" in lt.evidence and lt.n_windows == 0
    # ... and by its shape when the names say nothing: thousands of short sequences
    shape = reference.Layout(lambda i: f"seq{i}", [1500] * 5000)
    assert shape.kind == "transcriptome" and "too many and too short" in shape.evidence
    # a draft assembly of thousands of scaffolds with megabase ones is a genome
    draft = reference.Layout(lambda i: f"scaffold{i}", [3_000_000] * 5 + [20_000] * 4000)
    assert draft.kind == "genome" and len(draft.kept) == reference.MAX_CONTIGS and draft.describe()["other"]["sequences"] == 4005 - 64
    assert reference.Layout(lambda i: "", []).kind is None
    # one small genome
    assert reference.Layout(lambda i: "NC_000000.0", [29_903]).kind == "genome"


def test_window_length_gives_about_three_thousand_windows():
    assert reference.window_length(3_099_750_718) == 1_100_000              # a human genome: 2,818 windows and the contig ends
    assert reference.window_length(3_016_500) == 1100 and reference.window_length(30_000) == 10 and reference.window_length(1) == 1
    for total in (5_000, 123_456, 4_600_000, 12_000_000, 2_700_000_000, 16_000_000_000):
        w = reference.window_length(total)
        assert 2400 <= total / w <= 3000, (total, w)
    lay = reference.Layout(lambda i: T.CONTIGS[i][0], [ln for _, ln in T.CONTIGS])
    d = lay.describe()
    assert d["window_bases"] == 1100 and [c["windows"] for c in d["contigs"]] == [1819, 910, 15] and d["windows"] == 2744
    assert [c["first_window"] for c in d["contigs"]] == [0, 1819, 2729] and d["mitochondrial"] == "chrM"
    # two files with the same header get the same layout id; another reference another id
    assert lay.ref_id == reference.Layout(lambda i: T.CONTIGS[i][0], [ln for _, ln in T.CONTIGS]).ref_id
    assert lay.ref_id != reference.Layout(lambda i: T.CONTIGS[i][0], [ln + 1 for _, ln in T.CONTIGS]).ref_id


# ---------------------------------------------------------------------------------------- genome: exact
def test_strand_cigar_totals_and_mismatches_are_exact(genome):
    m, exp = genome["scan"]["merged"], genome["exp"]
    assert genome["scan"]["timings"]["tasks"] > 5 and set(m.aln) == set(exp["groups"])
    for gid, x in exp["groups"].items():
        al, a = x["aligned"], m.aln[gid]
        assert len(al) == x["mapped"] == a.mapped
        assert (a.fwd, a.rev) == (sum(1 for r in al if not r["reverse"]), sum(1 for r in al if r["reverse"]))
        got = (a.m_bases, a.ins_bases, a.del_bases, a.skip_bases, a.soft_bases, a.hard_bases, a.mism_bases, a.m_bases_nm, a.nm_negative)
        want = tuple(sum(r[k] for r in al) for k in ("m", "ins", "dele", "skip", "soft", "hard", "mism", "m")) + (0,)
        assert got == want, gid
        assert a.skip_bases > 0 and a.hard_bases > 0 and a.del_bases > 0 and a.rev > 100       # the file holds what the test is about
        # match + mismatch + inserted + soft-clipped = the read bases of the mapped reads
        assert a.m_bases + a.ins_bases + a.soft_bases == sum(r["length"] for r in al)
        assert sum(r["nm"] for r in al) == x["nm"] and a.mism_bases == x["nm"] - a.ins_bases - a.del_bases


def test_aligned_length_soft_clip_share_and_the_grids_are_exact(genome):
    m, exp = genome["scan"]["merged"], genome["exp"]
    for gid, x in exp["groups"].items():
        al, a = x["aligned"], m.aln[gid]
        alen = [r["m"] + r["ins"] for r in al]
        assert (a.h_alen == bincount(alen, "read_length_log40")).all() and int(a.h_alen_w.sum()) == sum(alen)
        assert (a.h_clip == bincount([r["soft"] / r["length"] for r in al], "fraction_0p01")).all()
        assert (a.g_len_aln == grid([r["length"] for r in al], "read_length_log20", alen, "read_length_log20")).all()
        idq = [math.inf if r["nm"] == 0 else -10 * math.log10(r["nm"] / r["cols"]) for r in al]
        assert (a.g_aln_id == grid(alen, "read_length_log20", idq, "identity_q_0p5")).all()
        assert (a.g_q_id == grid([r["q"] for r in al], "qscore_1", idq, "identity_q_0p5")).all()
        assert int(a.g_aln_id[:, -1].sum()) == x["perfect"]                    # no difference at all: the overflow cell
        cell = a.cell(None, unmapped=x["unmapped"])
        h = cell["hists"]["aln_aligned_length_hist"]
        assert derive.total(h) == x["mapped"] and derive.total(h, "w") == sum(alen) and h["missing"] == x["unmapped"]
        for name in ("aln_length_vs_aligned", "aln_identity_vs_aligned", "aln_q_calibration"):
            g = cell["hists2d"][name]
            assert sum(g["n"]) + g["outside"] + g["missing"] == x["reads"], name
        assert cell["hists2d"]["aln_identity_vs_aligned"]["outside"] >= x["perfect"]


def test_reads_and_bases_per_contig_and_depth_per_window_are_exact(genome):
    m, exp, ref = genome["scan"]["merged"], genome["exp"], genome["scan"]["plan"]["reference"]
    assert ref["kind"] == "genome" and ref["window_bases"] == 1100 and ref["windows"] == 2744
    lay = alnpass_layout(ref)
    crossing = spliced = 0
    for gid, x in exp["groups"].items():
        al, a = x["aligned"], m.aln[gid]
        per_base = depth_by_base(al, T.CONTIGS)
        want = window_sums(per_base, 1100)
        assert len(want) == 2744 and (a.cov[:-1] == want).all() and a.cov[-1] == 0
        assert int(a.cov.sum()) == sum(r["m"] + r["dele"] for r in al) == sum(int(p.sum()) for p in per_base)
        assert [int(v) for v in a.slot_n] == [sum(1 for r in al if r["contig"] == k) for k in range(3)] + [0]
        assert [int(v) for v in a.slot_w] == [sum(r["m"] + r["dele"] for r in al if r["contig"] == k) for k in range(3)] + [0]
        cell = a.cell(lay, unmapped=x["unmapped"])
        t = cell["counts"]["aln_contig_table"]
        assert t["cats"] == ["chrT1", "chrT2", "chrM", "other contigs"] and sum(t["n"]) == x["mapped"]
        s = cell["series"]["aln_depth_windows"]
        assert s["axis"] == "genome_window" and len(s["values"]["aligned_bases"]) == 2744
        for r in al:                                                            # both ways through the code are taken
            end = r["blocks"][-1][1]
            crossing += r["start"] // 1100 != (end - 1) // 1100 and r["skip"] == 0
            spliced += r["start"] // 1100 != (end - 1) // 1100 and r["skip"] > 0
    assert crossing > 50 and spliced > 500


def alnpass_layout(ref):
    from flowcellar.readdata.alnpass import _Described
    return _Described(ref)


def test_spliced_and_long_records_on_a_made_up_reference(tmp_path):
    """Seventy contigs, 73-base windows: records that cross many windows with and without a skipped region, a record
    that ends on the last base of its contig, a contig that has no window of its own."""
    import pysam
    contigs = [("big1", 50_000), ("big2", 50_000), ("big3", 50_011)] + [(f"small{i}", 1000 + (i % 7)) for i in range(67)]
    hd = {"HD": {"VN": "1.6", "SO": "unsorted"}, "SQ": [{"SN": n, "LN": ln} for n, ln in contigs]}
    lay = reference.Layout(lambda i: contigs[i][0], [ln for _, ln in contigs])
    assert lay.kind == "genome" and lay.window == 73 and len(lay.kept) == 64 and lay.describe()["other"]["sequences"] == 6
    other = [i for i in range(70) if i not in lay.kept]
    recs = [  # (contig, start, cigar)
        (0, 10, [(4, 5), (0, 400), (2, 7), (0, 300), (1, 4), (0, 200)]),            # long, unspliced, many windows
        (0, 70, [(0, 20), (3, 5000), (0, 30), (2, 3), (0, 50), (3, 1234), (0, 400)]),   # spliced over many windows
        (1, 146, [(0, 30), (3, 10), (0, 20)]),                                      # spliced inside one window
        (1, 100, [(0, 46)]),                                                        # ends on a window boundary
        (2, 50_011 - 300, [(0, 300)]),                                              # ends on the last base of its contig
        (2, 0, [(5, 10), (0, 73), (4, 3)]),                                         # exactly the first window
        (other[0], 5, [(0, 500)]), (other[1], 0, [(0, 100), (3, 300), (0, 100)]),   # contigs without a window
        (5, 0, [(0, 1000)]),                                                        # a whole small contig
    ]
    path = tmp_path / "made_up.bam"
    aligned = []
    with pysam.AlignmentFile(str(path), "wb", header=hd) as out:
        for k, (rid, start, cigar) in enumerate(recs * 3):
            a = pysam.AlignedSegment(out.header)
            qlen = sum(n for op, n in cigar if op in (0, 1, 4))
            a.query_name, a.query_sequence, a.flag = f"r{k}", "A" * qlen, 16 if k % 2 else 0
            a.query_qualities = array("B", [30]) * qlen
            a.reference_id, a.reference_start, a.mapping_quality, a.cigartuples = rid, start, 60, cigar
            a.set_tag("NM", sum(n for op, n in cigar if op in (1, 2)), "i")
            out.write(a)
            pos, blocks = start, []
            for op, n in cigar:
                if op in (0, 2):
                    blocks.append((pos, pos + n))
                if op in (0, 2, 3):
                    pos += n
            aligned.append({"contig": rid, "blocks": blocks})
    r = scan(path)
    (acc,) = r["merged"].aln.values()
    per_base = depth_by_base(aligned, contigs)
    kept = np.concatenate([np.add.reduceat(per_base[i], np.arange(0, contigs[i][1], 73)) for i in lay.kept])
    assert len(kept) == lay.n_windows == len(acc.cov) - 1 and (acc.cov[:-1] == kept).all()
    assert acc.cov[-1] == sum(int(per_base[i].sum()) for i in other) == 3 * 700
    assert int(acc.slot_n[-1]) == 6 and int(acc.slot_n.sum()) == 27 and int(acc.slot_w.sum()) == int(acc.cov.sum())
    assert acc.skip_bases == 3 * (5000 + 1234 + 10 + 300) and acc.hard_bases == 30 and acc.soft_bases == 3 * 8
    assert acc.mism_bases == 0 and acc.ins_bases == 12 and acc.del_bases == 30


# ------------------------------------------------------------------------------- start time, speed, tail
def test_output_length_quality_and_speed_by_clock_time_are_exact(genome):
    m, exp = genome["scan"]["merged"], genome["exp"]
    nl, nq, ns = axes.n_cells("read_length_log20"), axes.n_cells("qscore_1"), axes.n_cells("speed_5bps")
    for gid, x in exp["groups"].items():
        tg, t = x["tagged"], m.tagacc[gid]
        assert t.reads == x["reads"] == t.timed == len(tg)
        out: dict[int, np.ndarray] = {}
        for t10, length, *_ in tg:
            out.setdefault(t10, np.zeros(2, np.int64))
            out[t10] += (1, length)
        same_tables(as_dict(t.t10), out)
        assert len(out) > 150                                                   # start times really spread over the run
        same_tables(as_dict(t.tt_len), table([(t10 // 6, axes.bin_index(ln, "read_length_log20")) for t10, ln, *_ in tg], nl))
        same_tables(as_dict(t.tt_q), table([(t10 // 6, axes.bin_index(q, "qscore_1")) for t10, _, q, *_ in tg], nq))
        fast = [(t10, ln, du) for t10, ln, _, du, _ in tg if du >= 1.0]
        assert 0 < len(fast) < len(tg) and t.dur_reads == len(fast) and t.dur_bases == sum(f[1] for f in fast)
        assert t.dur_ms == sum(int(np.rint(f[2] * 1000.0)) for f in fast)
        same_tables(as_dict(t.tt_speed), table([(t10 // 6, axes.bin_index(ln / du, "speed_5bps")) for t10, ln, du in fast], ns))
        assert (t.h_speed == bincount([ln / du for _, ln, du in fast], "speed_5bps")).all()
        # the aligned reads by clock hour: aligned length and identity
        al, a = x["aligned"], m.aln[gid]
        same_tables(as_dict(a.tt_alen), table([(r["t10"] // 6, axes.bin_index(r["m"] + r["ins"], "read_length_log20")) for r in al], nl))
        idq = [math.inf if r["nm"] == 0 else -10 * math.log10(r["nm"] / r["cols"]) for r in al]
        same_tables(as_dict(a.tt_id), table([(r["t10"] // 6, axes.bin_index(v, "identity_q_0p5")) for r, v in zip(al, idq)], axes.n_cells("identity_q_0p5")))


def test_poly_a_tail_from_the_pt_tag(genome):
    m, exp = genome["scan"]["merged"], genome["exp"]
    for gid, x in exp["groups"].items():
        t, pts = m.tagacc[gid], [v[4] for v in x["tagged"]]
        found = [p for p in pts if p is not None and p > 0]
        assert (t.pt_tag, t.pt_found, t.pt_zero, t.pt_none) == (sum(p is not None for p in pts), len(found), pts.count(0), pts.count(-1))
        assert t.pt_zero > 20 and t.pt_none > 50 and t.pt_tag < t.reads         # every class is in the file
        assert (t.h_pt == bincount(found, "polya_1nt")).all()
        cell = t.polya_cell()
        h = cell["hists"]["rna_polya_hist"]
        assert derive.total(h) == len(found) and h["missing"] == t.reads - len(found) and h["axis"] == "polya_1nt"
        assert cell["scalars"] == {"reads": t.reads, "reads_with_tag": t.pt_tag, "tail_found_reads": len(found),
                                   "tail_zero_reads": pts.count(0), "tail_not_found_reads": pts.count(-1)}
        assert derive.median(h) == pytest.approx(float(np.median(found)), abs=1.0)


def test_start_time_is_read_in_every_form_and_a_file_without_tags_has_no_time(tmp_path):
    st = bam._State(bam.header_facts({}))
    want = float(1_772_442_252 // 600)                                         # 2026-03-02T09:04:12Z
    for text in ("2026-03-02T09:04:12.345000+00:00", "2026-03-02T09:04:12Z", "2026-03-02T09:04:12.345+00:00",
                 "2026-03-02T10:04:12.345000+01:00", "2026-03-02T03:34:12-05:30"):
        assert st.clock_bin(text, text[:15] + text[-6:]) == want, text
    assert st.clock_bin("yesterday", "x") == -1.0 and st.clock_bin("", "") == -1.0
    # a file whose reads carry no tag of the basecaller: no time table, no speed, no tail, and nothing breaks
    reads = T.make_reads(200, A, seed=8, alias="alpha")
    p = tmp_path / "bare.bam"
    T.write_aligned_bam(str(p), reads, dorado_tags=False, wrap=None)
    m = scan(p)["merged"]
    (tags,), (aln,) = m.tagacc.values(), m.aln.values()
    assert tags.timed == 0 and tags.dur_reads == 0 and tags.pt_tag == 0 and tags.time_cell() is None and tags.polya_cell() is None
    assert aln.time_missing == aln.mapped and aln.time_json() is None and aln.mapped > 150


def test_time_table_grows_merges_and_refuses_a_span_no_run_has():
    a, b = TimeTable(3), TimeTable(3)
    a.add(np.array([10, 12, 10]), np.array([0, 2, 0]))
    b.add(np.array([7, 12]), np.array([1, 2]), w=np.array([5.0, 7.0]))
    a.merge(b)
    assert a.lo == 7 and a.a.tolist() == [[0, 5, 0], [0, 0, 0], [0, 0, 0], [2, 0, 0], [0, 0, 0], [0, 0, 8]]
    back = TimeTable.from_json(3, a.to_json())
    assert back.lo == 7 and (back.a == a.a).all() and (TimeTable.from_state(3, a.to_state()).a == a.a).all()
    assert TimeTable(3).to_json() is None and TimeTable.from_json(3, None).lo is None
    far = TimeTable(3)
    far.add(np.array([0, alnacc.MAX_TIME_ROWS + 5]), np.array([0, 0]))
    assert far.broken and far.to_json() == {"broken": True}
    a.merge(far)
    assert a.broken and a.trimmed()[0] is None                                # one broken part spoils the sum, whatever the order


# ------------------------------------------------------------------------------------------ transcripts
def test_transcriptome_counts_per_transcript_and_gene(txome):
    r, exp = txome["scan"], txome["exp"]
    ref = r["plan"]["reference"]
    assert ref["kind"] == "transcriptome" and "windows" not in ref and r["timings"]["tasks"] > 5
    assert [c["transcripts"] for c in ref["length_classes"]] == [sum(1 for _, ln in exp["transcripts"] if lo <= ln < hi)
                                                                 for lo, hi in ((0, 1000), (1000, 2000), (2000, 4000), (4000, 10**9))]
    gene_of, n_genes = r["plan"]["genes"]
    assert n_genes == 140 and gene_of.tolist() == exp["genes"]
    x, a = exp["groups"][GID], r["merged"].aln[GID]
    want = np.bincount([t["transcript"] for t in x["tx"]], minlength=len(exp["transcripts"]))
    assert a.tx_reads == x["mapped"] == len(x["tx"]) and (a.ref_n == want).all() and a.cov is None
    assert (a.fwd, a.rev) == (sum(not t["reverse"] for t in x["tx"]), sum(t["reverse"] for t in x["tx"])) and a.rev > 50
    assert exp["secondary"] > 1000 and sum(v.secondary for v in r["merged"].groups.values()) == exp["secondary"]
    assert 50 < int((want > 0).sum()) < len(want) - 50                         # detected, and many that are not


def test_coverage_along_the_transcript_body_per_length_class(txome):
    exp, a = txome["exp"]["groups"][GID], txome["scan"]["merged"].aln[GID]
    want = np.zeros((4, 100), np.int64)
    per_class = np.zeros(4, np.int64)
    for t in exp["tx"]:
        k = 0 if t["length"] < 1000 else 1 if t["length"] < 2000 else 2 if t["length"] < 4000 else 3
        per_class[k] += 1
        for b in range(100):                                                    # the read covers the middle of step b
            if 200 * t["start"] <= (2 * b + 1) * t["length"] < 200 * t["end"]:
                want[k, b] += 1
    assert (a.body == want).all() and (a.tx_class_n == per_class).all() and per_class.min() > 20
    assert (a.body[:, -1] > a.body[:, 0]).all()                                # reads end at the 3' end: coverage rises along the body
    cell = a.tx_cell()
    assert cell["scalars"] == {"reads": int(per_class.sum()), "reads_under_1kb": int(per_class[0]), "reads_1_to_2kb": int(per_class[1]),
                               "reads_2_to_4kb": int(per_class[2]), "reads_4kb_and_longer": int(per_class[3])}
    h = cell["hists"]["rna_body_coverage_2_to_4kb"]
    assert h["axis"] == "rel_pos_0p01" and build.hist_full(h)[1:-1].tolist() == want[2].tolist()


def test_reads_per_transcript_and_the_saturation_curve(txome):
    a = txome["scan"]["merged"].aln[GID]
    counts = a.ref_n.astype(np.int64)
    dist = saturation.reads_per_item(counts)
    pos = counts[counts > 0]
    assert sum(dist["y"]["transcripts"]) == len(pos) and sum(dist["y"]["reads"]) == int(pos.sum())
    for lo, n, w in zip(dist["x"], dist["y"]["transcripts"], dist["y"]["reads"]):
        inside = pos[(pos >= lo) & (pos < 2 * lo)]
        assert (n, w) == (len(inside), int(inside.sum()))
    gene_of, n_genes = txome["scan"]["plan"]["genes"]
    genes = np.bincount(gene_of, weights=counts, minlength=n_genes).astype(np.int64)
    cv = saturation.curve(counts, genes)
    total = int(counts.sum())
    assert cv["x"][-1] == total and cv["y"]["transcripts"][-1] == len(pos) and cv["y"]["genes"][-1] == int((genes > 0).sum())
    assert cv["y"]["genes_10_reads"][-1] == int((genes >= 10).sum())
    assert all(b >= a for a, b in zip(cv["y"]["transcripts"], cv["y"]["transcripts"][1:]))             # more reads never find fewer

    def exact(c, n, k=1):                                                      # the same expectation in whole numbers
        big = int(sum(c))
        tot = Fraction(0)
        for ci in (int(v) for v in c if v):
            tot += 1 - sum(Fraction(math.comb(ci, j) * math.comb(big - ci, n - j), math.comb(big, n)) for j in range(k) if n - j >= 0)
        return float(tot)
    for i in (0, 3, 8, 12):
        assert cv["y"]["transcripts"][i] == pytest.approx(exact(counts, cv["x"][i]), abs=0.06)
        assert cv["y"]["genes"][i] == pytest.approx(exact(genes, cv["x"][i]), abs=0.06)
        assert cv["y"]["genes_10_reads"][i] == pytest.approx(exact(genes, cv["x"][i], 10), abs=0.06)
    assert saturation.curve(np.zeros(5, np.int64)) is None and saturation.reads_per_item([0, 0]) is None
    x = np.array([1.0, 1.5, 2, 3, 7.3, 10, 55.5, 1e3, 1e7, 3e9])
    assert np.allclose(saturation.lgamma(x), [math.lgamma(v) for v in x], rtol=1e-12, atol=1e-10)


def test_too_many_groups_for_per_transcript_counts_says_so(txome, monkeypatch):
    monkeypatch.setattr(bam, "REF_COUNT_CELLS", 100)                          # 420 transcripts x 1 group is already too many
    r = scan(txome["path"], task_bytes=30_000)
    a = r["merged"].aln[GID]
    assert a.ref_n is None and a.ref_counts_dropped and a.tx_reads == txome["scan"]["merged"].aln[GID].tx_reads
    assert (a.body == txome["scan"]["merged"].aln[GID].body).all()            # everything else is still counted


# ----------------------------------------------------------------- workers, cuts, merging: one result
def _states(m):
    return ({g: a.cell(None, unmapped=0) for g, a in m.aln.items()}, {g: a.time_json() for g, a in m.aln.items()},
            {g: (t.time_cell(), t.time_json(), t.polya_cell()) for g, t in m.tagacc.items()},
            {g: (a.tx_cell(), None if a.ref_n is None else a.ref_n.tolist()) for g, a in m.aln.items()})


@pytest.mark.parametrize("task_bytes", [1 << 30, 60_000, 9_000])
def test_result_does_not_depend_on_how_the_file_is_cut(genome, txome, task_bytes):
    for fx in (genome, txome):
        r = scan(fx["path"], task_bytes=task_bytes)
        assert r["state"] == "ok" and (task_bytes > 1 << 20 or r["timings"]["tasks"] > 5)
        assert _states(r["merged"]) == _states(fx["scan"]["merged"])


@pytest.mark.slow
@pytest.mark.parametrize("jobs", [1, 3])
def test_worker_pool_gives_the_same_alignment_numbers(genome, txome, jobs):
    for fx in (genome, txome):
        r = bam.scan_bam(str(fx["path"]), jobs=jobs, task_bytes=20_000, inprocess=False, stall_timeout_s=60)
        assert r["state"] == "ok" and r["timings"]["workers"] == jobs
        assert _states(r["merged"]) == _states(fx["scan"]["merged"])


def test_two_files_add_up_to_the_file_that_holds_both(genome, tmp_path):
    """Accumulators are sums: the reads of one file split over two files give the same numbers once merged."""
    reads = [r for r in genome["reads"] if r.run is A]
    whole, first, second = tmp_path / "whole.bam", tmp_path / "first.bam", tmp_path / "second.bam"
    T.write_aligned_bam(str(whole), reads, seed=5, secondary_rate=0, supplementary_frac=0)
    # the same records, in two files: every record keeps its position in the random stream by writing the whole file
    # and filtering it
    import pysam
    with pysam.AlignmentFile(str(whole)) as src, pysam.AlignmentFile(str(first), "wb", header=src.header) as o1, \
            pysam.AlignmentFile(str(second), "wb", header=src.header) as o2:
        for k, rec in enumerate(src):
            (o1 if k % 3 else o2).write(rec)
    w, a, b = (scan(p)["merged"] for p in (whole, first, second))
    for store, cls in (("aln", AlignAcc), ("tagacc", TagAcc)):
        left, right = getattr(a, store)[GID], getattr(b, store)[GID]
        merged_ab = cls.from_state(left.to_state())
        merged_ab.merge(right)
        merged_ba = cls.from_state(right.to_state())
        merged_ba.merge(left)
        one = getattr(w, store)[GID]
        for got in (merged_ab, merged_ba):
            if cls is AlignAcc:
                assert got.cell(alnpass_layout(scan(whole)["plan"]["reference"]), unmapped=0) == one.cell(alnpass_layout(scan(whole)["plan"]["reference"]), unmapped=0)
                assert got.time_json() == one.time_json()
            else:
                assert (got.time_cell(), got.time_json(), got.polya_cell()) == (one.time_cell(), one.time_json(), one.polya_cell())


def test_quartiles_of_a_table_equal_the_quartiles_of_each_histogram():
    rng = np.random.default_rng(3)
    for axis in ("read_length_log20", "qscore_1", "identity_q_0p5", "speed_5bps"):
        rows = rng.integers(0, 40, size=(30, axes.n_cells(axis)))
        rows[3] = 0                                                             # an hour without a read
        rows[5, :] = 0
        rows[5, 17] = 9                                                         # every read of the hour in one bin
        rows[7, 1:-1] = 0                                                       # only values outside the axis
        got = table_quantiles(rows, axis)
        for i, row in enumerate(rows):
            h = build.hist(axis, row)
            for k, q in enumerate((0.25, 0.5, 0.75)):
                want = derive.quantile(h, q)
                assert (np.isnan(got[i, k]) and want is None) or got[i, k] == pytest.approx(want, rel=1e-12), (axis, i, q)
