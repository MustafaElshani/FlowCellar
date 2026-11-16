"""Attaching files to runs: the ladder of evidence, one input for equal files, one source of record."""
from __future__ import annotations

import os
import shutil

import pytest
from rd_helpers import cand, do_passes, ident, ident_of, tier1_of

from flowcellar.model import partial as pm
from flowcellar.readdata import api, link, record
from flowcellar.testing import readdata as T

A = T.fake_run_ref("TST00001", 1)
B = T.fake_run_ref("TST00002", 2)
C = T.fake_run_ref("TST00003", 3)


def by_name(linked):
    return {os.path.basename(li.candidate.path): li for li in linked}


@pytest.fixture(scope="module")
def tree(tmp_path_factory):
    """An experiment folder with one sample sequenced on two flow cells, and files of every kind of evidence."""
    d = tmp_path_factory.mktemp("link_tree")
    exp = d / "experiment_one"
    (exp / "analysis" / "alignments").mkdir(parents=True)
    ra = T.make_reads(700, A, seed=1, alias="gamma", barcode="barcode01")
    rb = T.make_reads(500, B, seed=2, alias="gamma", barcode="barcode01")
    out = {"dir": d, "exp": str(exp), "ra": ra, "rb": rb}
    out["dorado_a"] = str(exp / "dorado_a.bam")
    out["e_dorado_a"] = T.write_dorado_bam(out["dorado_a"], ra)
    out["genome"] = str(exp / "analysis" / "alignments" / "gamma_genome_sorted.bam")
    out["e_genome"] = T.write_aligned_bam(out["genome"], ra + rb, seed=4)
    out["pychopper"] = str(exp / "analysis" / "pychopper_reports" / "gamma_pychopper_stats.tsv")
    T.write_pychopper_stats(out["pychopper"], pass_reads=1150, qc_fail=50, len_fail=5, primers_found=1100, rescued={2: 5}, unusable=45)
    out["plain"] = str(exp / "analysis" / "alignments" / "delta_genome_sorted.bam")       # no identifier inside
    out["e_plain"] = T.write_aligned_bam(out["plain"], T.make_reads(300, A, seed=7), seed=8, dorado_tags=False, wrap=None)
    out["sheet"] = T.write_sample_sheet(str(exp / "sheet.csv"), [{"experiment_id": "demo_experiment", "kit": T.KIT, "flow_cell_id": "TST00001",
                                                                 "barcode": "barcode01", "sample_id": "gamma", "alias": "gamma"}])
    return out


def test_dorado_bam_is_linked_by_its_header(tree, tmp_path):
    li, = api.link([cand("dorado.bam", tree["dorado_a"], sample="gamma", exp=tree["exp"])], [ident_of(A), ident_of(B)], probes_dir=str(tmp_path))
    assert li.run_keys == [A.run_key] and li.link["method"] == "bam_header_rg" and li.link["confidence"] == "certain"
    assert "PU and DS runid" in li.link["evidence"] and li.scope == "run"
    assert li.counted_for == {A.run_key: ["reads"]} and li.partial_id and li.needs == ["reads"]


def test_aligned_bam_is_linked_by_the_tags_on_its_reads_with_shares(tree, tmp_path):
    li, = api.link([cand("align.genome_bam", tree["genome"], sample="gamma", exp=tree["exp"])], [ident_of(A), ident_of(B)], probes_dir=str(tmp_path))
    assert li.run_keys == sorted([A.run_key, B.run_key]) and li.link["method"] == "read_tag" and li.link["confidence"] == "certain"
    assert li.scope == "shared" and abs(li.shares[A.run_key] - 700 / 1200) < 0.05 and abs(sum(li.shares.values()) - 1) < 1e-5
    assert li.counted_for == {A.run_key: ["reads", "genome"], B.run_key: ["reads", "genome"]}


def test_a_run_the_library_does_not_know_is_still_named(tree, tmp_path):
    li, = api.link([cand("align.genome_bam", tree["genome"], sample="gamma", exp=tree["exp"])], [ident_of(A)], probes_dir=str(tmp_path))
    assert set(li.run_keys) == {A.run_key, B.run_key}                    # reads of a run without a folder are not dropped


def test_sample_sheet_links_by_flow_cell_and_experiment(tree, tmp_path):
    sheet = cand("basecaller.sample_sheet", tree["sheet"], exp=tree["exp"])
    li, = api.link([sheet], [ident_of(A), ident_of(B)], probes_dir=str(tmp_path / "p1"))
    assert li.run_keys == [A.run_key] and li.link["method"] == "sample_sheet" and li.link["confidence"] == "certain"
    # the same flow cell loaded twice: the id alone does not say which loading
    import dataclasses
    a2 = dataclasses.replace(ident_of(A), run_key="TST00001_0000aaaa", run_short="0000aaaa", facts={"experiment_id": "another_experiment"})
    li, = api.link([sheet], [ident_of(A), a2], probes_dir=str(tmp_path / "p2"))
    assert li.run_keys == [A.run_key] and li.link["confidence"] == "certain"          # the experiment id of the sheet decides
    a3 = dataclasses.replace(a2, facts={"experiment_id": "demo_experiment"})
    li, = api.link([sheet], [ident_of(A), a3], probes_dir=str(tmp_path / "p3"))
    assert sorted(li.run_keys) == sorted([A.run_key, a3.run_key]) and li.link["confidence"] == "high"
    a4, a5 = dataclasses.replace(a2, facts={}), dataclasses.replace(ident_of(A), facts={})
    li, = api.link([sheet], [a5, a4], probes_dir=str(tmp_path / "p4"))
    assert len(li.run_keys) == 2 and li.link["confidence"] == "medium" and "does not say which" in li.link["evidence"]


def test_file_without_identifiers_follows_the_sample_name_inside_its_experiment_folder(tree, tmp_path):
    cands = [cand("align.genome_bam", tree["genome"], sample="gamma", exp=tree["exp"]),
             cand("pychopper.stats", tree["pychopper"], sample="gamma", exp=tree["exp"]),
             cand("pychopper.stats", tree["pychopper"].replace("gamma", "gamma"), sample="gamma", exp="/somewhere/else")]
    linked = api.link(cands[:2], [ident_of(A), ident_of(B)], config={"link": {"rules": [{"name": "experiment-folder", "roles": ["pychopper.stats"], "via": "sample", "within": "experiment"}]}},
                      probes_dir=str(tmp_path))
    py = by_name(linked)["gamma_pychopper_stats.tsv"]
    assert py.link["method"] == "sample_name_chain" and py.link["confidence"] == "high" and py.link["rule"] == "experiment-folder"
    assert sorted(py.run_keys) == sorted([A.run_key, B.run_key]) and py.scope == "sample_all_runs"      # one table, two flow cells: attached unsplit
    assert "gamma_genome_sorted.bam" in py.link["evidence"]
    # the same name in another experiment folder is not the same sample
    alone, = api.link([cands[2]], [ident_of(A), ident_of(B)], probes_dir=str(tmp_path))
    assert alone.run_keys == [] and alone.link["method"] is None and "no link rule" in alone.link["evidence"]


def test_sample_listed_in_a_sheet_links_a_table_when_no_read_file_is_left(tree, tmp_path):
    cands = [cand("basecaller.sample_sheet", tree["sheet"], exp=tree["exp"]), cand("pychopper.stats", tree["pychopper"], sample="gamma", exp=tree["exp"])]
    linked = by_name(api.link(cands, [ident_of(A), ident_of(B)], probes_dir=str(tmp_path)))
    py = linked["gamma_pychopper_stats.tsv"]
    assert py.run_keys == [A.run_key] and py.link["method"] == "sample_name_chain" and "sample sheet sheet.csv" in py.link["evidence"]


def test_file_without_identifiers_linked_by_a_path_rule_has_lower_confidence(tree, tmp_path):
    plain = cand("align.genome_bam", tree["plain"], sample="delta", exp=tree["exp"], relpath="experiment_one/analysis/alignments/delta_genome_sorted.bam")
    runs = [ident_of(A, run_dir=os.path.join(tree["exp"], "run_a")), ident_of(B, run_dir="/elsewhere/run_b")]
    unlinked, = api.link([plain], runs, config={}, probes_dir=str(tmp_path / "probes"))
    assert unlinked.run_keys == [] and unlinked.counted_for == {}
    cfg = {"link": {"rules": [{"name": "delta files", "roles": ["align.genome_bam"], "match": "*/analysis/alignments/delta_*"}]}}
    li, = api.link([plain], runs, config=cfg, probes_dir=str(tmp_path / "probes"))
    assert li.run_keys == [A.run_key]                                    # the one run folder inside its experiment folder
    assert li.link["method"] == "path_rule" and li.link["confidence"] == "medium" and li.link["rule"] == "delta files" and not li.link["confirmed"]
    # lower than every link by content
    by_content, = api.link([cand("align.genome_bam", tree["genome"], sample="gamma", exp=tree["exp"])], runs, probes_dir=str(tmp_path / "probes"))
    order = ["certain", "high", "medium", "low"]
    assert order.index(li.link["confidence"]) > order.index(by_content.link["confidence"])
    # its numbers are shown, but nothing is judged on them and they stay out of read fate
    refs = do_passes([li], tmp_path / "partials")
    assert refs[tree["plain"]].complete
    part = pm.read_partial(refs[tree["plain"]].path)
    assert part["routing"] == {"runs": {}, "unrouted": 300, "default_run": None}          # the pass itself attributes nothing
    fr = type("F", (), {"run_key": A.run_key, "expected": {"reads": 5000}})()
    rd = api.read_data(runs[0], [li], str(tmp_path / "partials"), tier1=tier1_of(fr, barcoded=False), config=cfg)
    assert sum(c["scalars"]["reads"] for c in rd.modules["len"]["cells"]) == 300
    assert any(n["code"] == "weak_link" for n in rd.modules["len"]["notes"])
    assert rd.inputs[0]["link"] == {"method": "path_rule", "confidence": "medium", "evidence": li.link["evidence"], "rule": "delta files", "confirmed": False}
    table = rd.modules["funnel"]["cells"][0]["counts"]["funnel_stage_table"]
    assert "genome_mapped" in table["missing"] and table["n"][table["cats"].index("genome_mapped")] == 0
    bundle = {"identity": {"run_key": A.run_key, "assay": {"group": "cdna", "chemistry": "r10.4.1"}, "basecalling": rd.basecalling},
              "modules": rd.modules, "samples": rd.samples, "tiers": {"read_level": rd.read_level}, "provenance": {"inputs": rd.inputs}}
    rules = {r["id"]: r for r in api.read_checks(bundle)}
    assert rules["len_n50"]["outcome"] == "not_evaluated" and "linked by a rule" in rules["len_n50"]["message"]
    assert rules["aln_mapped_frac"]["outcome"] == "not_evaluated"
    # a confirmed rule is trusted
    cfg["link"]["rules"][0]["confirmed"] = True
    li2, = api.link([plain], runs, config=cfg, probes_dir=str(tmp_path / "probes"))
    rd2 = api.read_data(runs[0], [li2], str(tmp_path / "partials"), tier1=tier1_of(fr, barcoded=False), config=cfg)
    assert not any(n["code"] == "weak_link" for n in rd2.modules["len"]["notes"])


def test_file_name_of_the_basecaller_links_with_lower_confidence_when_nothing_inside_does(tmp_path):
    p = tmp_path / f"{A.flow_cell_id}_pass_{A.run8}_{A.acq8}_0.bam"
    T.write_dorado_bam(str(p), T.make_reads(40, A, seed=2), header_rg=False, fn_tag=False, rg_tag=False, bc_tag=False)
    li, = api.link([cand("dorado.bam", str(p))], [ident_of(A)], probes_dir=str(tmp_path / "probes"))
    assert li.run_keys == [A.run_key] and li.link["method"] == "path_rule" and li.link["confidence"] == "medium" and li.link["rule"] == "file name"


def test_duplicate_file_under_two_roots_is_one_input_and_counted_once(tree, tmp_path):
    root_b = tmp_path / "root_b" / "experiment_one" / "analysis" / "alignments"
    root_b.mkdir(parents=True)
    copy = str(root_b / "gamma_genome_sorted.bam")
    shutil.copy(tree["genome"], copy)
    shutil.copy(tree["genome"] + ".bai", copy + ".bai")
    os.utime(copy, (1_700_000_000, 1_700_000_000))                        # another mtime, the same content
    c1 = cand("align.genome_bam", tree["genome"], sample="gamma", exp=tree["exp"], root_id="root_a")
    c2 = cand("align.genome_bam", copy, sample="gamma", exp=str(tmp_path / "root_b" / "experiment_one"), root_id="root_b")
    linked = api.link([c1, c2], [ident_of(A), ident_of(B)], probes_dir=str(tmp_path / "probes"))
    assert len(linked) == 1
    li = linked[0]
    assert li.candidate.path == tree["genome"] and li.alt_paths == [copy] and li.not_counted_reason is None
    refs = do_passes(linked, tmp_path / "partials")
    assert len(refs) == 1                                                 # read once
    fr = type("F", (), {"run_key": A.run_key, "expected": {"reads": 5000}})()
    rd = api.read_data(ident_of(A), linked, str(tmp_path / "partials"), tier1=tier1_of(fr), config={})
    assert sum(c["scalars"]["reads"] for c in rd.modules["len"]["cells"]) == 700          # not 1,400
    assert rd.inputs[0]["alt_paths"] == [copy] and len(rd.inputs) == 1


@pytest.fixture()
def two_sources(tree, tmp_path):
    """One sample aligned twice, into two folders, by the same rule: nothing separates the two files."""
    second = tmp_path / "experiment_one_again" / "analysis" / "alignments"
    second.mkdir(parents=True)
    other = str(second / "gamma_genome_sorted.bam")
    T.write_aligned_bam(other, tree["ra"], seed=99, mapped_frac=0.5)       # other content: another alignment of the same reads
    c1 = cand("align.genome_bam", tree["genome"], sample="gamma", exp=tree["exp"], rule="configured rule align.genome_bam")
    c2 = cand("align.genome_bam", other, sample="gamma", exp=tree["exp"], rule="configured rule align.genome_bam")
    return c1, c2, other


def test_two_candidate_sources_are_never_summed(tree, two_sources, tmp_path):
    c1, c2, other = two_sources
    linked = api.link([c1, c2], [ident_of(A), ident_of(B)], probes_dir=str(tmp_path / "probes"))
    a, b = linked
    # for run A both files claim the reads of the sample: neither is counted, both are named
    assert a.counted_for[A.run_key] == [] and b.counted_for[A.run_key] == []
    msg = b.not_counted_reason
    assert msg and msg.startswith("two candidate sources for") and tree["genome"] in msg and other in msg and "nothing was added up" in msg
    assert a.link["conflicts"][A.run_key]["reads"] == msg and "primary = true" in msg
    # the first file is still the one source for run B, whose reads the second file does not hold
    assert a.counted_for[B.run_key] == ["reads", "genome"] and a.not_counted_reason is None
    refs = do_passes(linked, tmp_path / "partials")
    assert list(refs) == [tree["genome"]]                                 # the file that is counted for nothing is not read at all
    fr = type("F", (), {"run_key": A.run_key, "expected": {"reads": 5000}})()
    rd = api.read_data(ident_of(A), linked, str(tmp_path / "partials"), tier1=tier1_of(fr), config={})
    for mod in ("len", "qual", "aln"):
        sec = rd.modules[mod]
        assert sec["cells"] == [] and sec["state"] == "skipped"            # not a sum, not one of the two picked silently
        note = next(n for n in sec["notes"] if n["code"] == "two_candidate_sources")
        assert note["message"].startswith("Not assessed: two candidate sources") and other in note["message"]
    assert rd.read_level["reads_covered"] is None and rd.read_level["state"] == "partial"
    assert any(w["code"] == "two_candidate_sources" for w in rd.warnings)
    table = rd.modules["funnel"]["cells"][0]["counts"]["funnel_stage_table"]
    assert table["sources"]["genome_mapped"] == "not assessed: two candidate sources" and "genome_mapped" in table["missing"]
    for i in rd.inputs:
        assert i["counted_for"] == [] and i["state"] == "skipped" and i["not_counted_reason"].startswith("two candidate sources")
    bundle = {"identity": {"run_key": A.run_key, "assay": {"group": "cdna", "chemistry": "r10.4.1"}, "basecalling": []}, "modules": rd.modules,
              "samples": rd.samples, "tiers": {"read_level": rd.read_level}, "provenance": {"inputs": rd.inputs}}
    rules = {r["id"]: r for r in api.read_checks(bundle)}
    assert rules["len_n50"]["outcome"] == "not_evaluated" and "two candidate sources" in rules["len_n50"]["message"]
    # run B is untouched by the conflict
    frb = type("F", (), {"run_key": B.run_key, "expected": {"reads": 5000}})()
    rdb = api.read_data(ident_of(B), linked, str(tmp_path / "partials"), tier1=tier1_of(frb), config={})
    assert sum(c["scalars"]["reads"] for c in rdb.modules["len"]["cells"]) == 500


def test_a_primary_rule_or_a_variant_separates_two_sources(tree, two_sources, tmp_path):
    import dataclasses
    c1, c2, other = two_sources
    runs = [ident_of(A), ident_of(B)]
    a, b = api.link([dataclasses.replace(c1, primary=True), c2], runs, probes_dir=str(tmp_path / "probes"))
    assert a.counted_for[A.run_key] == ["reads", "genome"] and b.counted_for[A.run_key] == [] and b.not_counted_reason is None
    assert "primary" in b.link["not_counted"]
    # variants: the configured primary variant is counted, the other one is listed and not read
    a, b = api.link([dataclasses.replace(c1, variant="main"), dataclasses.replace(c2, variant="alt")], runs,
                    config={"variants": {"primary": "main"}}, probes_dir=str(tmp_path / "probes"))
    assert a.counted_for[A.run_key] == ["reads", "genome"] and b.counted_for == {A.run_key: []} and b.not_counted_reason is None
    assert "version 'alt'" in b.link["not_counted"]
    assert list(do_passes([a, b], tmp_path / "partials")) == [tree["genome"]]
    # a file without a variant name is the primary one when a variant file exists beside it
    a, b = api.link([c1, dataclasses.replace(c2, variant="alt")], runs, probes_dir=str(tmp_path / "probes"))
    assert a.counted_for[A.run_key] == ["reads", "genome"] and b.counted_for[A.run_key] == []
    # two named variants and no rule: still two candidate sources
    a, b = api.link([dataclasses.replace(c1, variant="one"), dataclasses.replace(c2, variant="two")], runs, probes_dir=str(tmp_path / "probes"))
    assert a.counted_for[A.run_key] == [] and b.not_counted_reason.startswith("two candidate sources")


def test_an_alignment_without_a_variant_name_chains_to_the_primary_basecalled_file(tmp_path):
    """Basecalled files in two versions (main, alt), each linked by its header; alignment files that kept no
    identifier: ``x_genome_sorted.bam`` (no variant in its name) and ``x_alt_genome_sorted.bam``. With
    ``[variants] primary = "main"`` the file without a variant name is the main one: it chains to the main basecalled
    file and is the source of record of the alignment; the alt one is listed and not counted."""
    exp = tmp_path / "experiment_one"
    (exp / "analysis" / "alignments").mkdir(parents=True)
    reads = T.make_reads(400, A, seed=1, alias="kappa", barcode="barcode01")
    paths = {}
    for v in ("main", "alt"):
        (exp / "calls_basecalled" / v).mkdir(parents=True)
        paths[v] = str(exp / "calls_basecalled" / v / "kappa.bam")
        T.write_dorado_bam(paths[v], reads, min_qscore=7.0 if v == "main" else 8.0)
    plain = str(exp / "analysis" / "alignments" / "kappa_genome_sorted.bam")
    T.write_aligned_bam(plain, reads, seed=4, dorado_tags=False, wrap=None)
    alt = str(exp / "analysis" / "alignments" / "kappa_alt_genome_sorted.bam")
    T.write_aligned_bam(alt, reads, seed=5, dorado_tags=False, wrap=None)
    cands = [cand("dorado.bam", paths["main"], sample="kappa", variant="main", exp=str(exp)), cand("dorado.bam", paths["alt"], sample="kappa", variant="alt", exp=str(exp)),
             cand("align.genome_bam", plain, sample="kappa", exp=str(exp)), cand("align.genome_bam", alt, sample="kappa", variant="alt", exp=str(exp))]
    got = by_name(api.link(cands, [ident_of(A)], config={"variants": {"primary": "main"}}, probes_dir=str(tmp_path / "probes")))
    g, ga = got["kappa_genome_sorted.bam"], got["kappa_alt_genome_sorted.bam"]
    assert g.run_keys == [A.run_key] and g.link["method"] == "sample_name_chain" and "kappa.bam" in g.link["evidence"]
    assert g.counted_for == {A.run_key: ["genome"]}                              # the reads come from the main basecalled file
    assert ga.run_keys == [A.run_key] and ga.counted_for == {A.run_key: []} and "version 'alt'" in ga.link["not_counted"]
    # a section whose source is not the primary version says so; one that reads the primary version says nothing
    fr = type("F", (), {"run_key": A.run_key, "expected": {"reads": 400}})()
    only_alt = api.link([cands[0], cands[1], cands[3]], [ident_of(A)], config={"variants": {"primary": "main"}}, probes_dir=str(tmp_path / "probes"))
    do_passes(only_alt, tmp_path / "partials")
    rd = api.read_data(ident_of(A), only_alt, str(tmp_path / "partials"), tier1=tier1_of(fr), config={"variants": {"primary": "main"}})
    note, = [n for n in rd.modules["aln"]["notes"] if n["code"] == "source_variant"]
    assert "version 'alt'" in note["message"] and "primary version 'main'" in note["message"]
    assert not [n for n in rd.modules["len"]["notes"] if n["code"] == "source_variant"]
    # without a configured primary variant nothing says which basecalled file the plain alignment belongs to
    got = by_name(api.link(cands, [ident_of(A)], probes_dir=str(tmp_path / "probes2")))
    assert got["kappa_genome_sorted.bam"].run_keys == [] and got["kappa_alt_genome_sorted.bam"].run_keys == [A.run_key]


def test_role_preference_basecalled_file_gives_the_reads_alignment_gives_the_alignment(tree, tmp_path):
    cands = [cand("dorado.bam", tree["dorado_a"], sample="gamma", exp=tree["exp"]), cand("align.genome_bam", tree["genome"], sample="gamma", exp=tree["exp"])]
    d, g = api.link(cands, [ident_of(A), ident_of(B)], probes_dir=str(tmp_path))
    assert d.counted_for == {A.run_key: ["reads"]}
    assert g.counted_for == {A.run_key: ["genome"], B.run_key: ["reads", "genome"]}      # run B has no basecalled file: its reads come from the alignment
    assert record.needs_of("align.transcriptome_bam") == ("reads", "transcriptome")


def test_chunks_of_one_output_are_one_set_and_pass_and_fail_folders_belong_together(tmp_path):
    base = tmp_path / "out"
    paths = []
    for sub, seed in (("bam_pass/alpha", 1), ("bam_pass/alpha", 2), ("bam_fail/alpha", 3)):
        (base / sub).mkdir(parents=True, exist_ok=True)
        p = base / sub / f"{A.flow_cell_id}_{sub.split('/')[0][4:]}_alpha_{A.run8}_{A.acq8}_{seed}.bam"
        T.write_dorado_bam(str(p), T.make_reads(100 * seed, A, seed=seed, alias="alpha", barcode="barcode01"), min_qscore=None)
        paths.append(str(p))
    linked = api.link([cand("dorado.bam", p, sample="alpha", rule="built-in rule dorado.bam") for p in paths], [ident_of(A)], probes_dir=str(tmp_path / "probes"))
    assert all(li.counted_for == {A.run_key: ["reads"]} and li.not_counted_reason is None for li in linked)
    do_passes(linked, tmp_path / "partials")
    fr = type("F", (), {"run_key": A.run_key, "expected": {"reads": 5000}})()
    rd = api.read_data(ident_of(A), linked, str(tmp_path / "partials"), tier1=tier1_of(fr), config={})
    cells = {tuple(c["at"]): c["scalars"]["reads"] for c in rd.modules["len"]["cells"]}
    assert cells == {("s1", "pass"): 300, ("s1", "fail"): 300}             # two chunks added up; the fail folder as filter "fail"
    # only the reads given to samples are here: the share that passed cannot be judged on them
    assert rd.modules["len"]["population"]["stage"] == "classified" and rd.modules["yield.basecalled"]["summary"]["yield_pass_frac"] is None
    # with the reads without a barcode, passed and failed, the files hold everything the basecaller wrote
    more = []
    for sub, n, seed in (("bam_pass/unclassified", 50, 8), ("bam_fail/unclassified", 150, 9)):
        (base / sub).mkdir(parents=True)
        p = base / sub / f"{A.flow_cell_id}_{sub.split('/')[0][4:]}_unclassified_{A.run8}_{A.acq8}_0.bam"
        T.write_dorado_bam(str(p), T.make_reads(n, A, seed=seed), min_qscore=None)
        more.append(cand("dorado.bam", str(p), sample="unclassified", rule="built-in rule dorado.bam"))
    linked = api.link([li.candidate for li in linked] + more, [ident_of(A)], probes_dir=str(tmp_path / "probes"))
    do_passes(linked, tmp_path / "partials")
    rd = api.read_data(ident_of(A), linked, str(tmp_path / "partials"), tier1=tier1_of(fr), config={})
    assert rd.modules["len"]["population"]["stage"] == "basecalled_all"
    assert rd.modules["yield.basecalled"]["summary"] == {"yield_basecalled_reads": 800, "yield_basecalled_bases": sum(c["scalars"]["bases"] for c in rd.modules["len"]["cells"]),
                                                         "yield_pass_frac": 350 / 800}


def test_second_link_opens_no_input(tree, tmp_path, monkeypatch):
    cands = [cand("align.genome_bam", tree["genome"], sample="gamma", exp=tree["exp"]), cand("pychopper.stats", tree["pychopper"], sample="gamma", exp=tree["exp"]),
             cand("basecaller.sample_sheet", tree["sheet"], exp=tree["exp"]), cand("dorado.bam", tree["dorado_a"], sample="gamma", exp=tree["exp"])]
    runs = [ident_of(A), ident_of(B)]
    first = api.link(cands, runs, probes_dir=str(tmp_path / "probes"))
    assert os.path.isfile(tmp_path / "probes" / "stat.json.gz")

    def boom(*a, **k):
        raise AssertionError("an input was opened on the second link")
    monkeypatch.setattr(link, "identify_task", boom)
    monkeypatch.setattr(link, "file_fingerprint", boom)
    second = api.link(cands, runs, probes_dir=str(tmp_path / "probes"))
    view = lambda xs: [(li.candidate.path, li.fingerprint, li.run_keys, li.link, li.counted_for, li.partial_id) for li in xs]  # noqa: E731
    assert view(first) == view(second)


def test_unreadable_file_is_reported_and_does_not_stop_the_others(tree, tmp_path):
    bad = tmp_path / "broken.bam"
    bad.write_bytes(b"this is not a BAM file" * 100)
    linked = by_name(api.link([cand("dorado.bam", str(bad)), cand("dorado.bam", tree["dorado_a"], sample="gamma")], [ident_of(A)], probes_dir=str(tmp_path / "probes")))
    assert linked["broken.bam"].run_keys == [] and linked["broken.bam"].link["evidence"].startswith("not readable")
    assert linked["dorado_a.bam"].run_keys == [A.run_key]


@pytest.mark.slow
def test_link_through_the_worker_pool_gives_the_same_links(tree, tmp_path):
    cands = [cand("align.genome_bam", tree["genome"], sample="gamma", exp=tree["exp"]), cand("dorado.bam", tree["dorado_a"], sample="gamma", exp=tree["exp"]),
             cand("pychopper.stats", tree["pychopper"], sample="gamma", exp=tree["exp"])]
    assert not link.inprocess()
    pooled = api.link(cands, [ident_of(A), ident_of(B)], probes_dir=str(tmp_path / "pool"))
    os.environ[link.INPROCESS_ENV] = "1"
    try:
        local = api.link(cands, [ident_of(A), ident_of(B)], probes_dir=str(tmp_path / "local"))
    finally:
        del os.environ[link.INPROCESS_ENV]
    assert [(li.run_keys, li.link["method"], li.counted_for, li.fingerprint) for li in pooled] == [(li.run_keys, li.link["method"], li.counted_for, li.fingerprint) for li in local]


def test_sample_alias_map_joins_two_spellings(tmp_path):
    exp = tmp_path / "exp"
    (exp / "analysis").mkdir(parents=True)
    reads = T.make_reads(200, A, seed=3, alias="Sample-One_A", barcode="barcode02")
    d = str(exp / "d.bam")
    T.write_dorado_bam(d, reads)
    py = str(exp / "analysis" / "SampleOne_A_pychopper_stats.tsv")
    T.write_pychopper_stats(py, pass_reads=190, qc_fail=10, len_fail=0, primers_found=180, rescued={}, unusable=10)
    cands = [cand("dorado.bam", d, exp=str(exp)), cand("pychopper.stats", py, sample="SampleOne_A", exp=str(exp))]
    without = by_name(api.link(cands, [ident_of(A)], config={}, probes_dir=str(tmp_path / "p1")))
    assert without["SampleOne_A_pychopper_stats.tsv"].run_keys == []        # never merged silently
    cfg = {"samples": {"alias": {"Sample-One_A": "SampleOne_A"}}}
    with_map = by_name(api.link(cands, [ident_of(A)], config=cfg, probes_dir=str(tmp_path / "p2")))
    assert with_map["d.bam"].sample == "SampleOne_A"
    assert with_map["SampleOne_A_pychopper_stats.tsv"].run_keys == [A.run_key]


def test_candidates_for_run_finds_files_by_pattern_and_by_content(project):
    fr, proj = project
    run = ident(fr)
    found = api.candidates_for_run(run, [proj.project_dir])
    by_role = {}
    for c in found:
        by_role.setdefault(c.role, []).append(c)
    assert {r: len(v) for r, v in by_role.items()} == {"dorado.bam": 4, "align.genome_bam": 3, "pychopper.stats": 3, "quant.salmon": 3,
                                                       "quant.oarfish": 3, "quant.kallisto": 3, "basecaller.sample_sheet": 1}
    assert sorted(c.sample for c in by_role["dorado.bam"]) == ["alpha_ctrl", "alpha_trt", "beta_ctrl", "unclassified"]
    assert {c.sample for c in by_role["quant.salmon"]} == {c.sample for c in by_role["pychopper.stats"]} == set(proj.aliases)
    assert not any(c.path.startswith(fr.run_dir) for c in found)          # the run folder's own small files belong to run health
    explicit = api.candidates_for_run(run, explicit={"align.genome_bam": [os.path.join(proj.project_dir, "analysis", "alignments")]},
                                      sample_sheet=proj.sample_sheet)
    assert sorted(c.role for c in explicit) == ["align.genome_bam"] * 3 + ["basecaller.sample_sheet"]
    assert all(c.rule.startswith("named on the command line") for c in explicit)


def test_file_named_on_the_command_line_without_identifiers_is_linked_as_stated(tree, tmp_path):
    run = ident_of(A, run_dir=str(tmp_path / "run"))
    cands = api.candidates_for_run(run, [tree["plain"]])
    assert [c.role for c in cands] == ["align.genome_bam"]
    li, = api.link(cands, [run], probes_dir=str(tmp_path / "probes"))
    assert li.run_keys == [A.run_key] and li.link["method"] == "manual" and li.link["confidence"] == "certain"


def test_candidates_in_tree_resolves_patterns_from_their_anchors(project):
    from flowcellar.contracts import InputRule
    fr, proj = project
    top = os.path.dirname(proj.project_dir)
    rules = [InputRule(role="align.genome_bam", glob="**/analysis/alignments/{sample}_genome_sorted.bam", exclude=("beta_*",), source="config"),
             InputRule(role="quant.salmon", glob="**/analysis/counts/salmon/{sample}_salmon/aux_info/meta_info.json", source="config"),
             InputRule(role="quant.kallisto", glob="**/analysis/counts/kallisto/*/{sample}/run_info.json", source="config"),
             InputRule(role="dorado.bam", glob="**/*_basecalled/**/bam_pass/{sample}/*.bam", source="config"),
             InputRule(role="dorado.bam", glob="**/{variant:demo|other}_basecalled/bam_pass/{sample}/*.bam", source="config")]
    found = api.candidates_in_tree(top, rules, root_id="root_a")
    got = {}
    for c in found:
        got.setdefault(c.role, []).append(c.sample)
    assert sorted(got["align.genome_bam"]) == ["alpha_ctrl", "alpha_trt"]                 # the excluded file is left out
    assert sorted(got["quant.salmon"]) == sorted(got["quant.kallisto"]) == sorted(proj.aliases)
    assert sorted(got["dorado.bam"]) == ["alpha_ctrl", "alpha_trt", "beta_ctrl", "unclassified"]   # the first rule that matches a file wins
    assert all(c.root_id == "root_a" and not os.path.isabs(c.relpath) for c in found)
    assert all(c.variant is None for c in found if c.role == "dorado.bam")


def test_run_folder_with_live_basecalling_gives_its_own_files(tmp_path):
    """MinKNOW with basecalling on: bam_pass/ and bam_fail/ inside the run folder, and a summary with read lengths."""
    import dataclasses

    from flowcellar.testing.fixtures import RunSpec, make_run
    fr = make_run(tmp_path, RunSpec(flow_cell_id="TST00004", seed=4))
    ref = T.run_ref(fr)
    reads = T.make_reads(300, ref, seed=1, alias="barcode01", barcode="barcode01")
    for sub, part in (("bam_pass", reads[:200]), ("bam_fail", reads[200:])):
        d = os.path.join(fr.run_dir, sub, "barcode01")
        os.makedirs(d)
        T.write_dorado_bam(os.path.join(d, f"{ref.flow_cell_id}_{sub[4:]}_barcode01_{ref.run8}_{ref.acq8}_0.bam"), part, min_qscore=None)
    summary = os.path.join(fr.run_dir, "sequencing_summary_basecalled.txt")
    T.write_dorado_summary(summary, reads)
    run = dataclasses.replace(ident(fr), files={**fr.files, "minknow.sequencing_summary": summary})
    found = api.candidates_for_run(run)
    assert sorted(c.role for c in found) == ["dorado.bam", "dorado.bam", "dorado.summary"]
    linked = api.link(found, [run], probes_dir=str(tmp_path / "probes"))
    assert all(li.link["method"] == "run_dir" and li.link["confidence"] == "certain" and li.run_keys == [run.run_key] for li in linked)
    by_role = {}
    for li in linked:
        by_role.setdefault(li.candidate.role, []).append(li)
    assert all(li.counted_for == {run.run_key: ["reads"]} for li in by_role["dorado.bam"])
    assert by_role["dorado.summary"][0].counted_for == {run.run_key: []}         # the same reads: the BAM files are the source, the table is not added
    do_passes(linked, tmp_path / "partials")
    rd = api.read_data(run, linked, str(tmp_path / "partials"), tier1=tier1_of(fr), config={})
    assert {tuple(c["at"]): c["scalars"]["reads"] for c in rd.modules["len"]["cells"]} == {("s1", "pass"): 200, ("s1", "fail"): 100}
    assert rd.read_level["origin"] == "minknow_live_basecalling" and rd.basecalling[0]["live"] is True
    # without the BAM files the summary is the source, with pass and fail told per read
    only = api.link([c for c in found if c.role == "dorado.summary"], [run], probes_dir=str(tmp_path / "probes"))
    do_passes(only, tmp_path / "partials")
    rd2 = api.read_data(run, only, str(tmp_path / "partials"), tier1=tier1_of(fr), config={})
    assert sum(c["scalars"]["reads"] for c in rd2.modules["len"]["cells"]) == 300 and {c["at"][1] for c in rd2.modules["len"]["cells"]} == {"pass", "fail"}
    assert rd2.modules["len"]["population"]["source_role"] == "dorado.summary"


def test_the_two_alignments_of_one_sample_are_one_sample_whatever_the_ending():
    from flowcellar.readdata.finder import guess_sample
    for genome, tx in (("alpha_ctrl_genome_sorted.bam", "alpha_ctrl_transcriptome_sorted.bam"), ("alpha_ctrl_gnome_sorted.bam", "alpha_ctrl_txome_sorted.bam"),
                       ("alpha_ctrl.genome.bam", "alpha_ctrl.transcripts.bam"), ("alpha_ctrl_sorted.bam", "alpha_ctrl.bam")):
        assert guess_sample("align.genome_bam", "/x/" + genome) == "alpha_ctrl", genome
        assert guess_sample("align.transcriptome_bam", "/x/" + tx) == "alpha_ctrl", tx

