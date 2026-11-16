"""The alignment, over-time, poly(A) and transcript sections of a run (sections_align): what the bundle holds equals
what the generator wrote, per sample and added up; clock bins become bins counted from the run's time origin; and
what cannot be added up (quartiles, reads per transcript, the saturation curve) is worked out from the full tables."""
from __future__ import annotations

import os

import jsonschema
import numpy as np
import pytest
from rd_helpers import cand, do_passes, ident, project_candidates, tier1_of

from flowcellar import catalogue
from flowcellar.aggregates import build, derive
from flowcellar.model import axes, schema
from flowcellar.readdata import api, link, record, saturation, sections_align
from flowcellar.testing import readdata as T
from flowcellar.testing.fixtures import RunSpec, make_run

NEW = ("aln.time", "aln.windows", "dyn.basecalled", "rna.polya", "rna.tx", "rna.coverage")
ORIGIN = "2026-03-02T09:04:12.000000000Z"            # 252 s after the protocol start of the fixture run
BIN0_600 = 1_772_442_000 // 600                      # 09:00:00, the ten minutes that hold the origin
BIN0_3600 = 1_772_442_000 // 3600


def section_validator():
    sch = schema.load_schema()
    return jsonschema.Draft202012Validator({"$ref": "#/$defs/section", "$defs": sch["$defs"]})


def build_project(root, *, samples=("alpha_ctrl", "alpha_trt", "beta_ctrl"), reads=600, origin=ORIGIN, tx=True, seed=31):
    fr = make_run(root / "project_one", RunSpec(experiment="demo_experiment", sample="pool_1", seed=21))
    proj = T.make_project(str(root / "project_one"), T.run_ref(fr), aliases=samples, reads_per_sample=reads, seed=seed, quant=False)
    cands = project_candidates(proj)
    tx_exp = {}
    if tx:
        for i, alias in enumerate(proj.aliases):
            rs = T.make_reads(400, T.run_ref(fr), seed=40 + i, alias=alias, barcode=proj.barcodes[i])
            p = os.path.join(proj.project_dir, "analysis", "alignments", f"{alias}_txome_sorted.bam")
            tx_exp[alias] = T.write_transcriptome_bam(p, rs, seed=50 + i)
            cands.append(cand("align.transcriptome_bam", p, sample=alias, rule="configured rule", exp=proj.project_dir))
    run = ident(fr)
    linked = api.link(cands, [run], config={}, probes_dir=str(root / "probes"))
    do_passes(linked, root / "partials", [run.run_key])
    tier1 = tier1_of(fr)
    if origin:
        tier1["identity"]["protocol_run"] = {"time_origin": origin}
    rd = api.read_data(run, linked, str(root / "partials"), tier1=tier1, config={})
    return {"fr": fr, "proj": proj, "run": run, "linked": linked, "rd": rd, "tx": tx_exp, "tier1": tier1, "root": root}


@pytest.fixture(scope="module")
def data(tmp_path_factory):
    before = os.environ.get(link.INPROCESS_ENV)
    os.environ[link.INPROCESS_ENV] = "1"
    try:
        return build_project(tmp_path_factory.mktemp("align_sections"))
    finally:
        if before is None:
            os.environ.pop(link.INPROCESS_ENV, None)
        else:
            os.environ[link.INPROCESS_ENV] = before


def keys_of(d):
    return {s["name"]: s["key"] for s in d["rd"].samples}


def cell(sec, *at):
    return next(c for c in sec["cells"] if c["at"] == list(at))


def genome_exp(d, alias):
    return d["proj"].expected["files"][d["proj"].files["align.genome_bam"][alias]]["groups"][(d["proj"].run.run_key, alias)]


def dorado_exp(d, alias):
    return d["proj"].expected["files"][d["proj"].files["dorado.bam"][alias]]["groups"][(d["proj"].run.run_key, alias)]


def test_the_alignment_to_the_transcripts_is_read_now():
    assert record.INACTIVE_NEEDS == () and record.needs_of("align.transcriptome_bam") == ("reads", "transcriptome")
    assert sections_align.need_of("aln", ["reads", "genome"]) == "genome" and sections_align.need_of("aln", ["transcriptome"]) == "transcriptome"
    assert sections_align.need_of("tx", ["genome"]) == "transcriptome" and sections_align.need_of("time", ["reads"]) == "reads"
    assert sections_align.need_of("polya", []) == "reads" and sections_align.need_of("len", ["genome"]) == "reads"


def test_every_new_section_is_valid_complete_and_names_its_reads(data):
    rd, v = data["rd"], section_validator()
    for name in NEW + ("aln",):
        sec = rd.modules[name]
        assert [e.message for e in v.iter_errors(sec)] == [], name
        assert sec["state"] == "ok" and sec["coverage"]["kind"] == "full" and sec["detail"] == "standard", name
        assert sec["population"] and sec["population"]["stage"] and sec["tier"] == "read_level", name
        for key in sec["summary"]:                                 # every summary number has a name people can read
            assert catalogue.label(key) != key, (name, key)
    assert rd.modules["rna.coverage"]["dims"] == ["sample", "read_class"]
    assert rd.modules["dyn.basecalled"]["population"]["source_role"] == "dorado.bam"
    assert rd.modules["aln.windows"]["population"]["source_role"] == "align.genome_bam"
    assert rd.modules["rna.tx"]["population"]["source_role"] == "align.transcriptome_bam"
    used = axes.used_axes(rd.modules)
    assert {"clock_600s", "clock_3600s", "genome_window", "identity_q_1", "rel_pos_0p01", "polya_1nt", "speed_5bps"} <= set(used)
    by_role = {}
    for i in rd.inputs:
        by_role.setdefault(i["role"], i)
    assert by_role["align.transcriptome_bam"]["counted_for"] == ["transcriptome"] and by_role["align.transcriptome_bam"]["state"] == "ok"
    assert by_role["align.transcriptome_bam"]["used_by"] == ["rna.coverage", "rna.tx"]


def test_strand_cigar_and_contigs_per_sample_and_the_run_summary(data):
    rd, keys = data["rd"], keys_of(data)
    tot = {k: 0 for k in ("m", "ins", "dele", "soft", "mism", "fwd", "mapped", "chrm", "cols")}
    for alias in data["proj"].aliases:
        x = genome_exp(data, alias)
        al = x["aligned"]
        c = cell(rd.modules["aln"], keys[alias], "pass")
        sc = c["scalars"]
        assert (sc["forward_reads"], sc["reverse_reads"]) == (sum(not r["reverse"] for r in al), sum(r["reverse"] for r in al))
        want = {"match_or_mismatch_bases": "m", "inserted_bases": "ins", "deleted_bases": "dele", "skipped_bases": "skip",
                "soft_clipped_bases": "soft", "hard_clipped_bases": "hard", "mismatched_bases": "mism", "match_or_mismatch_bases_with_nm": "m"}
        assert {k: sc[k] for k in want} == {k: sum(r[f] for r in al) for k, f in want.items()}
        t = c["counts"]["aln_contig_table"]
        assert t["cats"] == ["chrT1", "chrT2", "chrM", "other contigs"] and t["w_name"] == "aligned_bases"
        assert t["n"] == [sum(r["contig"] == k for r in al) for k in range(3)] + [0]
        assert t["w"] == [sum(r["m"] + r["dele"] for r in al if r["contig"] == k) for k in range(3)] + [0]
        h = c["hists"]["aln_aligned_length_hist"]
        assert derive.total(h) == len(al) and derive.total(h, "w") == sum(r["m"] + r["ins"] for r in al) and h["missing"] == x["unmapped"]
        assert derive.total(c["hists"]["aln_soft_clip_hist"]) == len(al)
        for name, (gx, gy) in {"aln_length_vs_aligned": ("read_length_log20", "read_length_log20"), "aln_identity_vs_aligned": ("read_length_log20", "identity_q_1"),
                               "aln_q_calibration": ("qscore_1", "identity_q_1")}.items():
            g = c["hists2d"][name]
            assert (g["x"], g["y"]) == (gx, gy) and sum(g["n"]) + g["outside"] + g["missing"] == x["reads"], name
        for k, f in (("m", "m"), ("ins", "ins"), ("dele", "dele"), ("soft", "soft"), ("mism", "mism")):
            tot[k] += sum(r[f] for r in al)
        tot["fwd"] += sum(not r["reverse"] for r in al)
        tot["mapped"] += len(al)
        tot["chrm"] += sum(r["contig"] == 2 for r in al)
        tot["cols"] += sum(r["cols"] for r in al)
    s = rd.modules["aln"]["summary"]
    columns = tot["m"] + tot["ins"] + tot["dele"]
    assert s["aln_strand_frac"] == round(tot["fwd"] / tot["mapped"], 6) and s["aln_chrm_frac"] == round(tot["chrm"] / tot["mapped"], 6)
    assert s["aln_mismatch_frac"] == round(tot["mism"] / tot["cols"], 6) and s["aln_insertion_frac"] == round(tot["ins"] / columns, 6)
    assert s["aln_deletion_frac"] == round(tot["dele"] / columns, 6)
    assert s["aln_soft_clipped_frac"] == round(tot["soft"] / (tot["m"] + tot["ins"] + tot["soft"]), 6)
    assert s["aln_mean_depth"] == round((tot["m"] + tot["dele"]) / 3_016_500, 4)
    ref = rd.modules["aln"]["params"]["reference"]
    assert ref["kind"] == "genome" and ref["mitochondrial"] == "chrM" and ref["sequences"] == 3 and ref["total_length"] == 3_016_500
    assert "curves" not in rd.modules["aln"]["cells"][0]                           # the checks add these cells up: nothing here that cannot be added


def test_depth_windows_per_sample(data):
    rd, keys = data["rd"], keys_of(data)
    sec = rd.modules["aln.windows"]
    ref = sec["params"]["reference"]
    assert ref["window_bases"] == 1100 and ref["windows"] == 2744 and [c["name"] for c in ref["contigs"]] == ["chrT1", "chrT2", "chrM"]
    assert [(c["first_window"], c["windows"]) for c in ref["contigs"]] == [(0, 1819), (1819, 910), (2729, 15)]
    for alias in data["proj"].aliases:
        al = genome_exp(data, alias)["aligned"]
        per = [np.zeros(ln + 1, np.int64) for _, ln in T.CONTIGS]
        for r in al:
            for s, e in r["blocks"]:
                per[r["contig"]][s] += 1
                per[r["contig"]][e] -= 1
        want = np.concatenate([np.add.reduceat(np.cumsum(p)[:-1], np.arange(0, len(p) - 1, 1100)) for p in per])
        c = cell(sec, keys[alias], "pass")
        got = c["series"]["aln_depth_windows"]
        assert got["axis"] == "genome_window" and got["start"] == 0 and got["values"]["aligned_bases"] == want.tolist()
        assert c["scalars"]["aligned_bases"] == int(want.sum())


def test_clock_bins_become_bins_from_the_time_origin(data):
    rd, keys = data["rd"], keys_of(data)
    sec = rd.modules["dyn.basecalled"]
    clock = sec["params"]["clock"]
    assert clock["origin_known"] and clock["time_origin"] == "2026-03-02T09:04:12Z"
    assert (clock["bin0_600s"], clock["offset_600s_s"]) == ("2026-03-02T09:00:00Z", -252.0)
    assert (clock["bin0_3600s"], clock["offset_3600s_s"]) == ("2026-03-02T09:00:00Z", -252.0)
    # all passed reads of all samples, ten minutes by ten minutes
    names = list(data["proj"].aliases) + ["unclassified"]
    tagged = {a: (dorado_exp(data, a) if a != "unclassified" else
                  data["proj"].expected["files"][data["proj"].files["dorado.bam"]["unclassified"]]["groups"][(data["proj"].run.run_key, None)])["tagged"] for a in names}
    n10 = max(t[0] for ts in tagged.values() for t in ts) - BIN0_600 + 1
    reads, bases = np.zeros(n10, np.int64), np.zeros(n10, np.int64)
    for ts in tagged.values():
        for t10, length, *_ in ts:
            reads[t10 - BIN0_600] += 1
            bases[t10 - BIN0_600] += length
    s = sec["run_level"]["series"]["dyn_output_series_600s"]
    assert s["axis"] == "clock_600s" and s["start"] == 0 and set(s["values"]) == {"reads_pass", "bases_pass"}
    assert s["values"]["reads_pass"] == reads.tolist() and s["values"]["bases_pass"] == bases.tolist()
    assert sum(s["values"]["reads_pass"]) == rd.modules["yield.basecalled"]["summary"]["yield_basecalled_reads"]
    # one sample, hour by hour: output, and the quartiles of length, quality and speed
    for alias in data["proj"].aliases:
        ts = tagged[alias]
        c = cell(sec, keys[alias], "pass")
        hours = max(t[0] // 6 for t in ts) - BIN0_3600 + 1
        hr = np.zeros(hours, np.int64)
        hb = np.zeros(hours, np.int64)
        for t10, length, *_ in ts:
            hr[t10 // 6 - BIN0_3600] += 1
            hb[t10 // 6 - BIN0_3600] += length
        got = c["series"]["dyn_output_series"]
        assert got["axis"] == "clock_3600s" and got["values"]["reads"] == hr.tolist() and got["values"]["bases"] == hb.tolist()
        assert c["scalars"]["reads"] == c["scalars"]["reads_with_start_time"] == len(ts)
        fast = [(t10, ln, du) for t10, ln, _, du, _ in ts if du >= 1.0]
        assert c["scalars"]["reads_with_duration"] == len(fast) and c["scalars"]["bases_with_duration"] == sum(f[1] for f in fast)
        for cid, axis, values, digits in (("len_over_time", "read_length_log20", [(t[0] // 6, t[1]) for t in ts], 0),
                                          ("qual_over_time", "qscore_1", [(t[0] // 6, t[2]) for t in ts], 2),
                                          ("dyn_speed_over_time", "speed_5bps", [(t10 // 6, ln / du) for t10, ln, du in fast], 1)):
            cv = c["curves"][cid]
            by_hour = {}
            for h, val in values:
                by_hour.setdefault(h - BIN0_3600, []).append(val)
            assert cv["x"] == sorted(by_hour) and cv["y"]["reads"] == [len(by_hour[h]) for h in sorted(by_hour)], cid
            for k, h in enumerate(sorted(by_hour)):
                hist = build.hist(axis, np.bincount(axes.bin_index(np.asarray(by_hour[h], float), axis), minlength=axes.n_cells(axis)))
                for name, q in (("q25", 0.25), ("median", 0.5), ("q75", 0.75)):
                    want = derive.quantile(hist, q)
                    assert cv["y"][name][k] == (int(round(want)) if digits == 0 else round(want, digits)), (cid, h, name)
            # a quartile read from a 20-per-decade histogram lies within one bin of the true one
            if cid == "len_over_time":
                big = max(by_hour, key=lambda h: len(by_hour[h]))
                true = float(np.median(by_hour[big]))
                assert abs(np.log10(cv["y"]["median"][cv["x"].index(big)] / true)) < 0.05
    # the grids under run_level hold every passed read once
    g = sec["run_level"]["hists2d"]["len_over_time_pass"]
    assert (g["x"], g["y"]) == ("clock_3600s", "read_length_log20") and sum(g["n"]) + g["outside"] == int(reads.sum())
    assert sec["summary"]["dyn_reads_with_start_time_frac"] == 1.0 and 250 < sec["summary"]["dyn_speed_basecalled"] < 550


def test_aligned_length_and_identity_over_time(data):
    rd, keys = data["rd"], keys_of(data)
    sec = rd.modules["aln.time"]
    assert sec["params"]["clock"]["offset_3600s_s"] == -252.0
    total = 0
    for alias in data["proj"].aliases:
        al = genome_exp(data, alias)["aligned"]
        c = cell(sec, keys[alias], "pass")
        assert c["scalars"]["reads_with_start_time"] == len(al)
        total += len(al)
        by_hour = {}
        for r in al:
            by_hour.setdefault(r["t10"] // 6 - BIN0_3600, []).append(r)
        for cid, axis, f in (("aln_aligned_length_over_time", "read_length_log20", lambda r: r["m"] + r["ins"]),
                             ("aln_identity_over_time", "identity_q_0p5", lambda r: np.inf if r["nm"] == 0 else -10 * np.log10(r["nm"] / r["cols"]))):
            cv = c["curves"][cid]
            assert cv["x"] == sorted(by_hour) and cv["y"]["reads"] == [len(by_hour[h]) for h in sorted(by_hour)]
            h0 = sorted(by_hour)[len(by_hour) // 2]
            vals = np.asarray([f(r) for r in by_hour[h0]], float)
            hist = build.hist(axis, np.bincount(axes.bin_index(vals, axis), minlength=axes.n_cells(axis)))
            k = cv["x"].index(h0)
            assert cv["y"]["median"][k] == (int(round(derive.median(hist))) if cid.startswith("aln_aligned") else round(derive.median(hist), 2))
    g = sec["run_level"]["hists2d"]["aln_identity_over_time_pass"]
    assert (g["x"], g["y"]) == ("clock_3600s", "identity_q_0p5") and sum(g["n"]) + g["outside"] == total


def test_poly_a_tail_section(data):
    rd, keys = data["rd"], keys_of(data)
    sec = rd.modules["rna.polya"]
    found_all, reads_all = [], 0
    for alias in data["proj"].aliases:
        pts = [t[4] for t in dorado_exp(data, alias)["tagged"]]
        found = [p for p in pts if p is not None and p > 0]
        c = cell(sec, keys[alias], "pass")
        assert c["scalars"] == {"reads": len(pts), "reads_with_tag": sum(p is not None for p in pts), "tail_found_reads": len(found),
                                "tail_zero_reads": pts.count(0), "tail_not_found_reads": pts.count(-1)}
        h = c["hists"]["rna_polya_hist"]
        assert h["axis"] == "polya_1nt" and derive.total(h) == len(found) and h["missing"] == len(pts) - len(found)
        full = build.hist_full(h)
        assert full[1:-1].tolist() == np.bincount(found, minlength=1000).tolist()
    for c in sec["cells"]:
        if c["at"][1] == "pass":
            found_all.append(c["scalars"]["tail_found_reads"])
            reads_all += c["scalars"]["reads"]
    assert sec["summary"]["rna_polya_detect_frac"] == round(sum(found_all) / reads_all, 6) and 60 < sec["summary"]["rna_polya_median"] < 130


def test_transcripts_per_sample_and_added_up(data):
    rd, keys = data["rd"], keys_of(data)
    sec = rd.modules["rna.tx"]
    pooled = np.zeros(420, np.int64)
    for alias in data["proj"].aliases:
        e = data["tx"][alias]
        x = e["groups"][(data["proj"].run.run_key, alias)]
        counts = np.bincount([t["transcript"] for t in x["tx"]], minlength=420)
        genes = np.bincount(e["genes"], weights=counts, minlength=140).astype(np.int64)
        pooled += counts
        c = cell(sec, keys[alias], "pass")
        sc = c["scalars"]
        assert (sc["reads"], sc["mapped_reads"], sc["unmapped_reads"], sc["secondary_records"]) == (x["reads"], x["mapped"], x["unmapped"], x["secondary"])
        assert (sc["forward_reads"], sc["reverse_reads"]) == (sum(not t["reverse"] for t in x["tx"]), sum(t["reverse"] for t in x["tx"]))
        assert c["curves"]["rna_reads_per_transcript"] == saturation.reads_per_item(counts)
        cv = c["curves"]["rna_saturation_curve"]
        assert cv == saturation.curve(counts, genes) and cv["x"][-1] == x["mapped"] and cv["y"]["transcripts"][-1] == int((counts > 0).sum())
        assert "hists2d" not in c and "series" not in c and "counts" not in c
    genes = np.bincount(data["tx"][data["proj"].aliases[0]]["genes"], weights=pooled, minlength=140).astype(np.int64)
    s = sec["summary"]
    assert s["rna_tx_transcripts_detected"] == int((pooled > 0).sum()) and s["rna_tx_genes_detected"] == int((genes > 0).sum())
    assert sec["run_level"]["curves"]["rna_saturation_curve"] == saturation.curve(pooled, genes)
    assert sec["run_level"]["curves"]["rna_reads_per_transcript"] == saturation.reads_per_item(pooled)
    # added up, more transcripts are found than in any one sample, and fewer than the samples' numbers summed
    per = [cell(sec, keys[a], "pass")["curves"]["rna_saturation_curve"]["y"]["transcripts"][-1] for a in data["proj"].aliases]
    assert max(per) < s["rna_tx_transcripts_detected"] < sum(per)
    assert sec["params"]["reference"]["kind"] == "transcriptome" and "gene named in each reference name" in sec["params"]["genes_from"]
    assert 0.02 < s["rna_antisense_frac"] < 0.06 and 0.85 < s["rna_tx_aligned_frac"] < 0.95


def test_coverage_along_the_transcript_body_per_sample_and_class(data):
    rd, keys = data["rd"], keys_of(data)
    sec = rd.modules["rna.coverage"]
    classes = ["under_1kb", "1_to_2kb", "2_to_4kb", "4kb_and_longer"]
    assert [c["id"] for c in sec["params"]["length_classes"]] == classes
    reads = 0
    for alias in data["proj"].aliases:
        x = data["tx"][alias]["groups"][(data["proj"].run.run_key, alias)]
        want = np.zeros((4, 100), np.int64)
        n = np.zeros(4, np.int64)
        for t in x["tx"]:
            k = 0 if t["length"] < 1000 else 1 if t["length"] < 2000 else 2 if t["length"] < 4000 else 3
            n[k] += 1
            for b in range(100):
                if 200 * t["start"] <= (2 * b + 1) * t["length"] < 200 * t["end"]:
                    want[k, b] += 1
        for k, cid in enumerate(classes):
            c = cell(sec, keys[alias], cid)
            assert c["scalars"]["reads"] == int(n[k]) and build.hist_full(c["hists"]["rna_body_coverage"])[1:-1].tolist() == want[k].tolist()
        reads += int(n.sum())
    assert sec["summary"]["rna_body_coverage_reads"] == reads


# ----------------------------------------------------------------------------------------- other shapes
def test_thirteen_samples_keep_the_large_members_once_for_the_run(tmp_path, monkeypatch):
    monkeypatch.setenv(link.INPROCESS_ENV, "1")
    names = tuple(f"sample_{i:02d}" for i in range(13))
    d = build_project(tmp_path, samples=names, reads=120, tx=False)
    rd, v = d["rd"], section_validator()
    for name in ("aln", "aln.time", "aln.windows", "dyn.basecalled", "rna.polya"):
        assert rd.modules[name]["detail"] == "coarse" and [e.message for e in v.iter_errors(rd.modules[name])] == [], name
    aln = rd.modules["aln"]
    assert all("hists2d" not in c and "counts" in c for c in aln["cells"])
    assert set(aln["run_level"]["hists2d"]) == {"aln_length_vs_aligned", "aln_identity_vs_aligned", "aln_q_calibration"}
    mapped = sum(c["scalars"]["mapped_reads"] for c in aln["cells"])
    g = aln["run_level"]["hists2d"]["aln_length_vs_aligned"]
    assert sum(g["n"]) + g["outside"] == mapped
    win = rd.modules["aln.windows"]
    assert win["cells"] == [] and win["state"] == "ok"
    total = win["run_level"]["series"]["aln_depth_windows"]["values"]["aligned_bases"]
    assert len(total) == 2744 and sum(total) == sum(c["scalars"]["match_or_mismatch_bases"] + c["scalars"]["deleted_bases"] for c in aln["cells"])
    assert all("curves" not in c and "series" not in c for c in rd.modules["dyn.basecalled"]["cells"])
    assert all("curves" not in c for c in rd.modules["aln.time"]["cells"]) and "aln_identity_over_time_pass" in rd.modules["aln.time"]["run_level"]["hists2d"]
    assert rd.modules["rna.polya"]["cells"][0]["hists"]["rna_polya_hist"]["axis"] == "polya_5nt"
    assert aln["cells"][0]["hists"]["aln_aligned_length_hist"]["axis"] == "read_length_log20"


def test_two_files_of_one_sample_add_up_and_another_reference_is_kept_apart(tmp_path, monkeypatch):
    monkeypatch.setenv(link.INPROCESS_ENV, "1")
    fr = make_run(tmp_path / "p", RunSpec(experiment="demo_experiment", sample="pool_1", seed=21))
    ref = T.run_ref(fr)
    reads = T.make_reads(900, ref, seed=5, alias="alpha", barcode="barcode01")
    d = tmp_path / "p" / "analysis" / "alignments"
    d.mkdir(parents=True)
    whole = T.write_aligned_bam(str(d / "alpha_genome_sorted.bam"), reads, seed=7, secondary_rate=0, supplementary_frac=0)
    # the same reads in two files of one output folder (a pipeline that writes chunks)
    chunks = tmp_path / "p" / "chunks" / "bam_pass" / "alpha"
    chunks.mkdir(parents=True)
    import pysam
    with pysam.AlignmentFile(str(d / "alpha_genome_sorted.bam")) as src, pysam.AlignmentFile(str(chunks / "a_0.bam"), "wb", header=src.header) as o1, \
            pysam.AlignmentFile(str(chunks / "a_1.bam"), "wb", header=src.header) as o2:
        for k, rec in enumerate(src):
            (o1 if k % 2 else o2).write(rec)
    run = ident(fr)
    tier1 = tier1_of(fr)
    tier1["identity"]["protocol_run"] = {"time_origin": ORIGIN}

    def sections(cands, sub):
        linked = api.link(cands, [run], config={}, probes_dir=str(tmp_path / sub / "probes"))
        do_passes(linked, tmp_path / sub / "partials", [run.run_key])
        return api.read_data(run, linked, str(tmp_path / sub / "partials"), tier1=tier1, config={})
    one = sections([cand("align.genome_bam", str(d / "alpha_genome_sorted.bam"), sample="alpha", rule="r")], "one")
    two = sections([cand("align.genome_bam", str(chunks / "a_0.bam"), sample="alpha", rule="chunks"),
                    cand("align.genome_bam", str(chunks / "a_1.bam"), sample="alpha", rule="chunks")], "two")
    assert whole["groups"][(ref.run_key, "alpha")]["mapped"] == one.modules["aln"]["cells"][0]["scalars"]["mapped_reads"]
    for name in ("aln", "aln.time", "aln.windows", "dyn.basecalled", "rna.polya"):
        a, b = one.modules[name], two.modules[name]
        assert a["cells"] == b["cells"] and a.get("run_level") == b.get("run_level") and a["summary"] == b["summary"], name
    # a second sample aligned to another reference: its contigs and its depth are not added to the first one's
    other = T.make_reads(300, ref, seed=6, alias="beta", barcode="barcode02")
    T.write_aligned_bam(str(d / "beta_genome_sorted.bam"), other, seed=8, contigs=(("chrT1", 2_000_000), ("chrT2", 1_500_000), ("chrM", 16_500)))
    both = sections([cand("align.genome_bam", str(d / "alpha_genome_sorted.bam"), sample="alpha", rule="r"),
                     cand("align.genome_bam", str(d / "beta_genome_sorted.bam"), sample="beta", rule="r")], "both")
    assert [n["code"] for n in both.modules["aln"]["notes"]] == ["several_references"]
    assert len(both.modules["aln.windows"]["cells"]) == 1 and len(both.modules["aln"]["cells"]) == 2
    assert "counts" not in both.modules["aln"]["cells"][1] and both.modules["aln"]["cells"][1]["scalars"]["mapped_reads"] > 200


def test_without_a_time_origin_bins_count_from_the_first_read_and_early_reads_fall_in_bin_zero(tmp_path, monkeypatch):
    monkeypatch.setenv(link.INPROCESS_ENV, "1")
    d = build_project(tmp_path / "a", reads=150, origin=None, tx=False)
    clock = d["rd"].modules["dyn.basecalled"]["params"]["clock"]
    assert clock["origin_known"] is False and clock["time_origin"] is None
    s = d["rd"].modules["dyn.basecalled"]["run_level"]["series"]["dyn_output_series_600s"]["values"]["reads_pass"]
    assert s[0] > 0 and sum(s) == d["rd"].modules["yield.basecalled"]["summary"]["yield_basecalled_reads"]
    # a run whose time origin lies ten hours after its first read: nothing is lost, the early reads are in bin 0
    late = build_project(tmp_path / "b", reads=150, origin="2026-03-02T19:04:12.000000000Z", tx=False)
    sec = late["rd"].modules["dyn.basecalled"]
    s2 = sec["run_level"]["series"]["dyn_output_series_600s"]["values"]["reads_pass"]
    assert sum(s2) == sum(s) and s2[0] > 100 and len(s2) < len(s)
    assert sec["params"]["clock"]["bin0_600s"] == "2026-03-02T19:00:00Z"


def test_files_without_the_basecallers_tags_give_no_over_time_section(tmp_path, monkeypatch):
    monkeypatch.setenv(link.INPROCESS_ENV, "1")
    fr = make_run(tmp_path / "p", RunSpec(experiment="demo_experiment", sample="pool_1", seed=21))
    reads = T.make_reads(300, T.run_ref(fr), seed=5, alias="alpha")
    p = tmp_path / "p" / "alpha_genome_sorted.bam"
    T.write_aligned_bam(str(p), reads, dorado_tags=False, wrap=None)
    run = ident(fr)
    linked = api.link([cand("align.genome_bam", str(p), sample="alpha", rule="r")], [run], config={"link": {"rules": []}}, probes_dir=str(tmp_path / "probes"))
    for li in linked:                                                             # no identifier inside: linked by hand for this test
        li.run_keys, li.counted_for = [run.run_key], {run.run_key: ["reads", "genome"]}
        li.link.update(method="path_rule", confidence="high", confirmed=True)
    do_passes(linked, tmp_path / "partials", [run.run_key])
    rd = api.read_data(run, linked, str(tmp_path / "partials"), tier1=tier1_of(fr), config={})
    assert "aln" in rd.modules and "aln.windows" in rd.modules
    assert "dyn.basecalled" not in rd.modules and "rna.polya" not in rd.modules and "aln.time" not in rd.modules
    assert rd.modules["aln"]["cells"][0]["scalars"]["forward_reads"] > 100


def test_whole_number_qualities_give_the_value_itself_as_the_hourly_median(tmp_path, monkeypatch):
    """Files whose read quality is a whole number (integer ``qs`` tags): every value of a 1 Q step is that number, so
    the median of an hour is a whole number, the one a sort of the values gives; never a point inside the step."""
    import dataclasses
    import math
    monkeypatch.setenv(link.INPROCESS_ENV, "1")
    fr = make_run(tmp_path / "project_one", RunSpec(experiment="demo_experiment", sample="pool_1", seed=23))
    ref = T.run_ref(fr)
    reads = [dataclasses.replace(r, q=float(round(r.q))) for r in T.make_reads(900, ref, seed=5, alias="alpha", barcode="barcode01")]
    p = str(tmp_path / "project_one" / "calls.bam")
    T.write_dorado_bam(p, reads)
    run = ident(fr)
    linked = api.link([cand("dorado.bam", p, sample="alpha")], [run], config={}, probes_dir=str(tmp_path / "probes"))
    do_passes(linked, tmp_path / "partials", [run.run_key])
    tier1 = tier1_of(fr)
    tier1["identity"]["protocol_run"] = {"time_origin": ORIGIN}
    rd = api.read_data(run, linked, str(tmp_path / "partials"), tier1=tier1, config={})
    sec = rd.modules["dyn.basecalled"]
    assert sec["params"]["qual_whole_numbers"] is True and "whole number" in sec["params"]["qual_quartile_rule"]
    c, = [x for x in sec["cells"] if x["at"][1] == "pass"]
    cv = c["curves"]["qual_over_time"]
    by_hour: dict[int, list[float]] = {}
    for r in reads:
        by_hour.setdefault(T.clock_bin(r.epoch, 3600) - BIN0_3600, []).append(r.q)
    assert all(float(v).is_integer() for name in ("q25", "median", "q75") for v in cv["y"][name])
    assert cv["y"]["reads"] == [len(by_hour[h]) for h in cv["x"]]
    for k, h in enumerate(cv["x"]):
        vals = sorted(by_hour[h])
        assert cv["y"]["median"][k] == vals[max(0, math.ceil(0.5 * len(vals)) - 1)], (h, cv["y"]["median"][k])
    # length is never treated this way, and float qualities keep the rank inside a step
    assert not all(float(v).is_integer() for v in cell(build_project(tmp_path / "floats", tx=False)["rd"].modules["dyn.basecalled"],
                                                       "s1", "pass")["curves"]["qual_over_time"]["y"]["median"])


def test_table_quantiles_of_whole_numbers():
    rows = np.zeros((2, axes.n_cells("qscore_1")), np.int64)
    rows[0, 1 + 18] = 60                                                        # sixty reads of Q 18, forty of Q 19
    rows[0, 1 + 19] = 40
    rows[1, 1 + 20] = 10
    q = sections_align.table_quantiles(rows, "qscore_1", whole=True)
    assert q.tolist() == [[18.0, 18.0, 19.0], [20.0, 20.0, 20.0]]
    inside = sections_align.table_quantiles(rows, "qscore_1")
    assert 18.0 < inside[0, 1] < 19.0 and inside[0, 2] > 19.0
