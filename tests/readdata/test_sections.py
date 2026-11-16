"""The read-level sections of a run: exact totals, the read population on every section, read fate with unknown
stages, the cache of passes, and what the index and the status engine get."""
from __future__ import annotations

import copy
import dataclasses
import os

import jsonschema
import pytest
from rd_helpers import bundle_of, cand, do_passes, ident, ident_of, project_candidates, tier1_of

from flowcellar.aggregates import derive, merge
from flowcellar.contracts import PassOptions
from flowcellar.model import partial as pm
from flowcellar.model import schema
from flowcellar.readdata import api, bam, fate, link
from flowcellar.testing import readdata as T
from flowcellar.testing.fixtures import RunSpec, make_run

READ_LEVEL = ("yield.basecalled", "len", "qual", "aln")


def stage_table(sec, at=None):
    cell = sec["cells"][0] if at is None else next(c for c in sec["cells"] if c["at"] == [at])
    t = cell["counts"]["funnel_stage_table"]
    return {s: (None if s in (t.get("missing") or []) else n) for s, n in zip(t["cats"], t["n"])}, t


def defs(name):
    sch = schema.load_schema()
    ref = "#/$defs/" + name
    return jsonschema.Draft202012Validator({"$ref": ref, "$defs": sch["$defs"]})


def test_per_sample_reads_and_bases_equal_what_the_generator_wrote(project_data):
    d = project_data
    proj, rd = d["proj"], d["rd"]
    keys = {s["name"]: s["key"] for s in rd.samples}
    assert [s["name"] for s in rd.samples] == proj.aliases + ["No barcode found"]
    assert [s["barcode"] for s in rd.samples[:3]] == proj.barcodes and all(s["expected"] and s["source"] == "basecaller_sample_sheet" for s in rd.samples[:3])
    len_cells = {tuple(c["at"]): c for c in rd.modules["len"]["cells"]}
    for alias in proj.aliases:
        want = proj.expected["files"][proj.files["dorado.bam"][alias]]["groups"][(proj.run.run_key, alias)]
        c = len_cells[(keys[alias], "pass")]
        assert c["scalars"] == {"reads": want["reads"], "bases": want["bases"]}
        assert derive.total(c["hists"]["len_hist"], "n") == want["reads"] and derive.total(c["hists"]["len_hist"], "w") == want["bases"]
        q = next(x for x in rd.modules["qual"]["cells"] if x["at"] == [keys[alias], "pass"])
        assert derive.total(q["hists"]["qual_hist"]) == want["reads"] and q["hists"]["qual_hist"]["missing"] == 0
        g = q["hists2d"]["qual_vs_length"]
        assert (g["x"], g["y"]) == ("read_length_log20", "qscore_1") and sum(g["n"]) + g["outside"] == want["reads"] and sum(g["w"]) + g["w_outside"] == want["bases"]
        gw = proj.expected["files"][proj.files["align.genome_bam"][alias]]["groups"][(proj.run.run_key, alias)]
        a = next(x for x in rd.modules["aln"]["cells"] if x["at"] == [keys[alias], "pass"])["scalars"]
        assert (a["reads"], a["mapped_reads"], a["unmapped_reads"], a["secondary_records"], a["supplementary_records"], a["edit_distance"], a["alignment_columns"], a["perfect_reads"]) == \
            (gw["reads"], gw["mapped"], gw["unmapped"], gw["secondary"], gw["supplementary"], gw["nm"], gw["cols"], gw["perfect"])
    un = len_cells[("unclassified", "pass")]["scalars"]
    assert un["reads"] == proj.expected["unclassified"]
    assert sum(c["scalars"]["reads"] for c in rd.modules["len"]["cells"]) == proj.expected["basecalled"] == rd.read_level["reads_covered"]
    assert rd.modules["yield.basecalled"]["summary"]["yield_basecalled_reads"] == proj.expected["basecalled"]


def test_every_read_level_section_names_its_population(project_data):
    rd = project_data["rd"]
    for mod in READ_LEVEL:
        pop = rd.modules[mod]["population"]
        assert pop and pop["stage"] and pop["source_role"] and pop["evidence"], mod
        assert rd.modules[mod]["coverage"]["kind"] == "full"
    reads, aln = rd.modules["len"]["population"], rd.modules["aln"]["population"]
    # the basecalled files with the reads that got no barcode: everything the basecaller kept
    assert (reads["stage"], reads["min_qscore"], reads["trimmed"], reads["both_end_barcodes"], reads["source_role"]) == ("basecalled_pass", 7.0, False, True, "dorado.bam")
    # the alignment of the cDNA step's output: full-length reads, trimmed, told by the read names and the aligner's input
    assert (aln["stage"], aln["trimmed"], aln["source_role"]) == ("full_length_trimmed", True, "align.genome_bam")
    assert any("segment coordinates" in e for e in aln["evidence"]) and any("full-length FASTQ" in e for e in aln["evidence"])
    assert rd.read_level["population"] == reads and rd.read_level["state"] == "complete"
    assert rd.read_level["samples_expected"] == rd.read_level["samples_with_data"] == 3
    assert rd.modules["cdna"].get("population") is None and rd.modules["cdna"]["tier"] == "assay"       # collected: names its file, not a population


def test_read_fate_from_sequenced_to_counted_towards_a_transcript(project_data):
    d = project_data
    proj, rd, fr = d["proj"], d["rd"], d["fr"]
    keys = {s["name"]: s["key"] for s in rd.samples}
    run_counts, table = stage_table(rd.modules["funnel"])
    assert table["cats"] == list(fate.STAGES) and table["role"] == "funnel" and not table.get("missing")
    per = {}
    for alias in proj.aliases:
        e = proj.expected["samples"][alias]
        g = proj.expected["files"][proj.files["dorado.bam"][alias]]["groups"][(proj.run.run_key, alias)]
        once = e["classified"] - (g["split_children"] - g["split_parents"])              # a split read counts once
        counts, t = stage_table(rd.modules["funnel.samples"], keys[alias])
        assert counts == {"classified": once, "full_length": e["full_length"], "genome_mapped": e["genome_mapped"], "quantified": e["quantified"]}
        assert "pychopper" in t["sources"]["full_length"] and t["sources"]["quantified"].startswith("salmon")
        per[alias] = counts
    un = proj.expected["files"][proj.files["dorado.bam"]["unclassified"]]["groups"][(proj.run.run_key, None)]
    assert run_counts == {"sequenced": fr.expected["reads"],
                          "basecalled": sum(p["classified"] for p in per.values()) + un["reads"] - (un["split_children"] - un["split_parents"]),
                          "classified": sum(p["classified"] for p in per.values()), "full_length": sum(p["full_length"] for p in per.values()),
                          "genome_mapped": sum(p["genome_mapped"] for p in per.values()), "quantified": sum(p["quantified"] for p in per.values())}
    sm = rd.modules["funnel"]["summary"]
    assert sm["funnel_usable_stage"] == "quantified" and sm["funnel_usable_frac"] == round(run_counts["quantified"] / fr.expected["reads"], 6)
    assert sm["funnel_reads_sequenced"] == fr.expected["reads"]
    # stages never grow: along the stages each one descends from, and here also in the stored order
    assert fate.never_grow(run_counts) == []
    vals = [run_counts[s] for s in fate.STAGES]
    assert all(a >= b for a, b in zip(vals, vals[1:]))
    for p in per.values():
        assert fate.never_grow(p, fate.SAMPLE_STAGES) == []


def test_a_stage_without_a_source_is_not_known_never_zero(project, tmp_path):
    fr, proj = project
    run = ident(fr)
    cands = project_candidates(proj, roles=("pychopper.stats", "align.genome_bam"))         # no basecalled file, no quantification
    linked = api.link(cands, [run], probes_dir=str(tmp_path / "probes"))
    do_passes(linked, tmp_path / "partials")
    rd = api.read_data(run, linked, str(tmp_path / "partials"), tier1=tier1_of(fr), config={})
    counts, t = stage_table(rd.modules["funnel"])
    assert counts["basecalled"] is None and counts["quantified"] is None
    assert t["missing"] == ["basecalled", "quantified"] and t["n"][1] == 0                  # a number only because the array needs one
    assert t["sources"]["basecalled"].startswith("not known: no basecalled file is linked")
    assert t["sources"]["quantified"].startswith("not known")
    # without the basecalled files, the reads given to each sample come from the pychopper table, and say so
    assert counts["classified"] == sum(s["pychopper"]["input"] for s in proj.expected["samples"].values())
    assert "PassReads + QcFail" in t["sources"]["classified"]
    assert counts["full_length"] == sum(s["full_length"] for s in proj.expected["samples"].values())
    assert rd.modules["funnel"]["state"] == "partial" and rd.modules["funnel"]["partial_reason"] == "inputs_missing"
    sm = rd.modules["funnel"]["summary"]
    assert sm["funnel_usable_frac"] is None and "not known" in rd.modules["funnel"]["params"]["usable_not_assessed_reason"]
    # the reads here are the cDNA step's output: no basecalled yield is claimed from them
    assert "yield.basecalled" not in rd.modules and rd.modules["len"]["population"]["stage"] == "full_length_trimmed"
    rules = {r["id"]: r for r in api.read_checks(bundle_of(tier1_of(fr), rd))}
    assert rules["funnel_usable_frac"]["outcome"] == "not_evaluated" and rules["funnel_usable_frac"]["value"] is None
    assert "yield_pass_frac" not in rules and "funnel_basecall_loss_frac" not in rules
    # no tier 1 at all: "sequenced" itself is not known
    rd2 = api.read_data(run, linked, str(tmp_path / "partials"), tier1=None, config={})
    c2, t2 = stage_table(rd2.modules["funnel"])
    assert c2["sequenced"] is None and "sequenced" in t2["missing"] and rd2.read_level["reads_expected"] is None


def test_a_stage_that_grows_is_stored_as_not_known_with_a_warning(project, tmp_path):
    fr, proj = project
    run = ident(fr)
    alias = proj.aliases[0]
    big = T.write_salmon_meta(str(tmp_path / "counts" / "salmon" / f"{alias}_salmon"), 10_000_000, 9_000_000)      # more than was ever sequenced
    cands = [c for c in project_candidates(proj) if not (c.role == "quant.salmon" and c.sample == alias)]
    cands.append(cand("quant.salmon", big, sample=alias, rule="configured rule quant.salmon", exp=proj.project_dir))
    linked = api.link(cands, [run], probes_dir=str(tmp_path / "probes"))
    do_passes(linked, tmp_path / "partials")
    rd = api.read_data(run, linked, str(tmp_path / "partials"), tier1=tier1_of(fr), config={})
    key = next(s["key"] for s in rd.samples if s["name"] == alias)
    counts, t = stage_table(rd.modules["funnel.samples"], key)
    assert counts["quantified"] is None and "is larger than the reads full-length" in t["sources"]["quantified"] and "9,000,000" in t["sources"]["quantified"]
    assert fate.never_grow(counts, fate.SAMPLE_STAGES) == []
    assert any(w["code"] == "funnel_stage_grew" for w in rd.warnings)
    run_counts, _ = stage_table(rd.modules["funnel"])
    assert run_counts["quantified"] is None and fate.never_grow(run_counts) == []


def test_never_grow_follows_descent_not_position():
    ok = {"sequenced": 100, "basecalled": 90, "classified": 80, "full_length": 70, "genome_mapped": 60, "quantified": 65}
    assert fate.never_grow(ok) == []                         # counted towards a transcript without a genome alignment: not a subset, no violation
    assert fate.never_grow({**ok, "quantified": 71}) == [("quantified", 71, "full_length", 70)]
    assert fate.never_grow({"sequenced": 100, "basecalled": None, "classified": 101}) == [("classified", 101, "sequenced", 100)]
    assert fate.never_grow({"sequenced": None, "full_length": 5, "genome_mapped": 6}) == [("genome_mapped", 6, "full_length", 5)]


def test_inputs_not_read_yet_are_waiting_and_nothing_is_judged_on_them(project, tmp_path):
    fr, proj = project
    run = ident(fr)
    linked = api.link(project_candidates(proj), [run], probes_dir=str(tmp_path / "probes"))
    rd = api.read_data(run, linked, str(tmp_path / "partials"), tier1=tier1_of(fr), config={})
    bams = sorted(list(proj.files["dorado.bam"].values()) + list(proj.files["align.genome_bam"].values()))
    assert rd.waiting == bams and rd.partial_ids == []
    for mod in ("len", "qual", "aln"):
        assert rd.modules[mod]["state"] == "partial" and rd.modules[mod]["partial_reason"] == "waiting" and rd.modules[mod]["cells"] == []
    assert rd.modules["cdna"]["state"] == "ok"                              # the small files are there at once
    assert rd.read_level["state"] == "partial" and rd.read_level["reads_covered"] is None
    b = bundle_of(tier1_of(fr), rd)
    rules = {r["id"]: r for r in api.read_checks(b)}
    assert rules["len_n50"]["outcome"] == "not_evaluated" and "not every read" in rules["len_n50"]["message"]
    cols = api.card_columns(b)
    assert cols["read_data"] == "partial" and cols["t2_files_waiting"] == len(bams) and cols["t2_n50"] is None and cols["t2_passed_gb"] is None
    assert all(r["reads"] is None for r in api.sample_rows(b))
    counts, t = stage_table(rd.modules["funnel"])
    assert counts["basecalled"] is None and counts["genome_mapped"] is None and "not read whole" in t["sources"]["genome_mapped"]
    # a section whose files wait does not claim to be whole: its coverage is the share of bytes read, here none
    for mod in ("len", "qual", "aln"):
        cov = rd.modules[mod]["coverage"]
        assert cov["kind"] == "byte_ranges" and cov["fraction"] == 0.0 and cov["records"] == 0, cov
    assert rd.read_level["coverage"]["kind"] != "full"
    # half read: the share follows the bytes
    do_passes([li for li in linked if li.candidate.role == "align.genome_bam"][:1], tmp_path / "partials")
    rd = api.read_data(run, linked, str(tmp_path / "partials"), tier1=tier1_of(fr), config={})
    cov = rd.modules["aln"]["coverage"]
    assert rd.modules["aln"]["partial_reason"] == "waiting" and 0.0 < cov["fraction"] < 1.0 and cov["kind"] == "byte_ranges"


def test_a_conflict_on_one_sample_does_not_blank_the_other_samples(project, tmp_path):
    """Two candidate genome alignments for one sample: that sample is "not assessed", the others keep their numbers."""
    fr, proj = project
    run = ident(fr)
    first = proj.aliases[0]
    second = tmp_path / "again" / "analysis" / "alignments"
    second.mkdir(parents=True)
    other = str(second / f"{first}_gnome_sorted.bam")
    T.write_aligned_bam(other, T.make_reads(300, T.run_ref(fr), seed=77, alias=first), seed=78, mapped_frac=0.5)
    cands = project_candidates(proj, roles=("pychopper.stats", "align.genome_bam", "quant.salmon"))
    cands.append(cand("align.genome_bam", other, sample=first, rule="configured rule align.genome_bam", exp=proj.project_dir))
    linked = api.link(cands, [run], probes_dir=str(tmp_path / "probes"))
    do_passes(linked, tmp_path / "partials")
    rd = api.read_data(run, linked, str(tmp_path / "partials"), tier1=tier1_of(fr), config={})
    keys = {s["name"]: s["key"] for s in rd.samples}
    for alias in proj.aliases[1:]:
        counts, t = stage_table(rd.modules["funnel.samples"], keys[alias])
        want = proj.expected["files"][proj.files["align.genome_bam"][alias]]["groups"][(run.run_key, alias)]["mapped"]
        assert counts["genome_mapped"] == want, "a sample whose own file is read whole has its number"
        assert "every record read" in t["sources"]["genome_mapped"]
    counts, t = stage_table(rd.modules["funnel.samples"], keys[first])
    assert counts["genome_mapped"] is None and t["sources"]["genome_mapped"] == "not assessed: two candidate sources"
    run_counts, rt = stage_table(rd.modules["funnel"])
    assert run_counts["genome_mapped"] is None, "the flow cell total is not known: never a sum of two of three samples"
    assert rt["sources"]["genome_mapped"] in ("not assessed: two candidate sources",
                                              "not known for the whole flow cell: 2 of 3 samples have a source")
    assert "branches" in rd.modules["funnel"]["params"]


def test_second_pass_reads_nothing(project_data, monkeypatch):
    d = project_data
    opened = []
    monkeypatch.setattr(bam, "scan_bam", lambda *a, **k: opened.append(a) or (_ for _ in ()).throw(AssertionError("a BAM was read again")))
    for li in d["linked"]:
        ref = api.read_pass(li, str(d["lib"] / "partials"), [d["run"].run_key], PassOptions())
        if li.candidate.role in link.PASS_ROLES:
            first = d["refs"][li.candidate.path]
            assert ref.cached and ref.complete and ref.partial_id == first.partial_id == li.partial_id and ref.path == first.path
            assert not first.cached
        else:
            assert ref is None                                              # small collected files need no pass
    assert opened == []
    rd = api.read_data(d["run"], d["linked"], str(d["lib"] / "partials"), tier1=d["tier1"], config={})
    assert rd.modules == d["rd"].modules and rd.partial_ids == d["rd"].partial_ids and len(rd.partial_ids) == 7


def test_partials_are_in_the_planned_format_and_pure(project_data):
    d = project_data
    li = next(x for x in d["linked"] if x.candidate.role == "align.genome_bam")
    p = pm.read_partial(d["refs"][li.candidate.path].path)
    assert p["format"] == "flowcellar.partial" and p["role"] == "align.genome_bam" and pm.is_complete(p)
    assert p["partial_id"] == pm.partial_id(li.fingerprint, "align.genome_bam", link.PASS_MODULES, "") == li.partial_id
    assert p["input"]["fingerprint"] == li.fingerprint and p["modules"] == link.PASS_MODULES and p["config_hash"] == ""
    assert p["coverage"]["kind"] == "full" and p["coverage"]["records"] == p["coverage"]["records_total"] == p["facts"]["records"]
    assert p["routing"]["runs"] == {d["run"].run_key: p["facts"]["primary"]} and p["routing"]["unrouted"] == 0
    assert {f["module"] for f in p["fragments"]} == {"len", "qual", "aln", "time", "polya"}     # a genome: no fragment tx
    assert all(f["dims"] == [] and f["cells"][0]["at"] == [] and f["run_key"] == d["run"].run_key for f in p["fragments"])
    assert p["facts"]["index"]["records"] == p["facts"]["records"] and p["timings"]["workers"] >= 1
    assert pm.fragments_for(p, d["run"].run_key, "len")[0]["cells"][0]["scalars"]["reads"] == p["facts"]["primary"]


def test_sections_inputs_samples_and_rules_validate_against_the_schema(project_data):
    d = project_data
    rd = d["rd"]
    for name, sec in rd.modules.items():
        assert [e.message for e in defs("section").iter_errors(sec)] == [], name
    for i in rd.inputs:
        assert [e.message for e in defs("input").iter_errors(i)] == [], i["id"]
    for s in rd.samples:
        assert [e.message for e in defs("sample").iter_errors(s)] == []
    assert [e.message for e in defs("tiers/properties/read_level").iter_errors(rd.read_level)] == []
    assert [e.message for e in defs("identity/properties/basecalling").iter_errors(rd.basecalling)] == []
    rules = api.read_checks(d["bundle"])
    assert rules and [e.message for e in defs("status/properties/rules").iter_errors(rules)] == []
    assert all(r["tier"] == "read_level" and (r["severity"] == "none") == (r["outcome"] != "triggered") for r in rules)
    for w in rd.warnings:
        assert [e.message for e in defs("note").iter_errors(w)] == []
    # ids are rd1, rd2, ... and every section points at inputs that exist
    ids = [i["id"] for i in rd.inputs]
    assert ids == [f"rd{n}" for n in range(1, len(ids) + 1)]
    assert all(set(sec.get("inputs", [])) <= set(ids) for sec in rd.modules.values())
    assert all(i["needs"] and i["link"]["method"] and i["fingerprint"]["value"] for i in rd.inputs)


def test_provenance_inputs_say_what_each_file_was_counted_for(project_data):
    rd, proj = project_data["rd"], project_data["proj"]
    by_role = {}
    for i in rd.inputs:
        by_role.setdefault(i["role"], []).append(i)
    assert all(i["counted_for"] == ["reads"] and i["link"]["method"] == "bam_header_rg" and i["state"] == "ok" for i in by_role["dorado.bam"])
    assert all(i["counted_for"] == ["genome"] and i["needs"] == ["genome", "reads"] and i["link"]["method"] == "read_tag" for i in by_role["align.genome_bam"])
    assert all(i["counted_for"] == ["collected"] and i["link"]["method"] == "sample_name_chain" and i["link"]["confidence"] == "high" for i in by_role["pychopper.stats"])
    assert by_role["basecaller.sample_sheet"][0]["link"]["method"] == "sample_sheet"
    g = by_role["align.genome_bam"][0]
    assert g["records"] == proj.expected["samples"][proj.aliases[0]]["full_length"] and g["used_by"] == ["aln", "aln.time", "aln.windows", "funnel", "funnel.samples"]
    assert set(by_role["dorado.bam"][0]["used_by"]) == {"len", "qual", "yield.basecalled", "dyn.basecalled", "rna.polya", "funnel", "funnel.samples"}
    bc, = rd.basecalling
    assert (bc["basecaller"], bc["model"], bc["model_tier"], bc["min_qscore"], bc["trim"], bc["barcode_both_ends"], bc["primary"]) == \
        ("dorado", T.MODEL, "sup", 7.0, False, True, True)


def test_summaries_are_read_from_the_cells(project_data):
    rd = project_data["rd"]
    named = [c for c in rd.modules["len"]["cells"] if c["at"][1] == "pass"]
    tot = merge.rollup(named, ["sample", "filter"], [])[0]
    h = tot["hists"]["len_hist"]
    sm = rd.modules["len"]["summary"]
    assert sm["len_n50"] == round(derive.nx(h, 0.5)) and sm["len_median"] == round(derive.median(h)) and sm["len_mean"] == round(derive.mean(h), 1)
    lengths = []
    for f in project_data["proj"].expected["files"].values():
        for g in f["groups"].values():
            if "mapq" in g and not g["mapq"]:
                lengths += g["lengths"]                     # the unaligned files are the reads source
    assert sm["len_short_read_frac"] == round(sum(1 for x in lengths if x < 200) / len(lengths), 6)        # exact: 200 b is an edge for whole lengths
    assert abs(sm["len_mean"] - sum(lengths) / len(lengths)) < 0.06
    aln = rd.modules["aln"]["summary"]
    a = merge.rollup(rd.modules["aln"]["cells"], ["sample", "filter"], [])[0]["scalars"]
    assert aln["aln_mapped_frac"] == round(a["mapped_reads"] / a["reads"], 6) and 0.9 < aln["aln_identity_median"] < 1
    assert aln["aln_supplementary_frac"] == round(a["reads_with_supplementary"] / a["mapped_reads"], 6)
    cd = rd.modules["cdna"]
    sc = [c["scalars"] for c in cd["cells"]]
    assert cd["summary"]["cdna_full_length_frac"] == round(sum(s["primers_found"] + s["rescued_reads"] for s in sc) / sum(s["pass_reads"] for s in sc), 6)
    assert cd["summary"]["cdna_samples"] == 3 and all({"input", "written", "len_fail", "qc_fail", "unusable", "rescue_segments", "umi_detected"} <= set(s) for s in sc)
    assert all(sum(c["counts"]["strand"]["n"]) == c["scalars"]["primers_found"] for c in cd["cells"])


def test_index_columns_and_sample_rows(project_data):
    d = project_data
    cols = api.card_columns(d["bundle"])
    rd, proj = d["rd"], d["proj"]
    assert set(cols) == {"read_data", "t2_population", "t2_samples", "t2_samples_with_data", "t2_reads", "t2_passed_gb", "t2_n50", "t2_q_med",
                         "t2_ident_med", "t2_mapped_frac", "t2_usable", "coverage", "t2_files_waiting", "t2_conflicts", "t2_full_length"}
    assert cols["t2_conflicts"] == 0, "no file of this run is left out by two candidate sources"
    assert cols["read_data"] == "complete" and cols["t2_population"] == "basecalled_pass" and cols["t2_samples"] == cols["t2_samples_with_data"] == 3
    assert cols["t2_reads"] == proj.expected["basecalled"] and cols["t2_files_waiting"] == 0
    assert cols["t2_passed_gb"] == round(sum(c["scalars"]["bases"] for c in rd.modules["len"]["cells"]) / 1e9, 4)
    assert cols["t2_n50"] == rd.modules["len"]["summary"]["len_n50"] and cols["t2_q_med"] == rd.modules["qual"]["summary"]["qual_median"]
    assert cols["t2_ident_med"] == rd.modules["aln"]["summary"]["aln_identity_median"] and cols["t2_mapped_frac"] == rd.modules["aln"]["summary"]["aln_mapped_frac"]
    assert cols["t2_usable"] == rd.modules["funnel"]["summary"]["funnel_usable_frac"] and cols["coverage"] == round(proj.expected["basecalled"] / d["fr"].expected["reads"], 6)
    rows = api.sample_rows(d["bundle"])
    assert [r["name"] for r in rows] == proj.aliases and all(r["run"] == d["run"].run_key and r["key"].startswith("s") for r in rows)
    for r in rows:
        assert r["reads"] == proj.expected["samples"][r["name"]]["classified"] and r["bases"] == proj.expected["samples"][r["name"]]["bases"]
        assert r["n50"] > 500 and 5 < r["q_med"] < 30 and r["barcode"].startswith("barcode") and r["sample_uid"] == r["name"]
    empty = {"tiers": {"read_level": {"state": "absent"}}, "modules": {}, "samples": []}
    assert api.card_columns(empty)["read_data"] == "absent" and api.card_columns(empty)["t2_reads"] is None and api.sample_rows(empty) == []


def test_run_without_linked_files_has_no_read_data(project):
    fr, _ = project
    rd = api.read_data(ident(fr), [], None, tier1=tier1_of(fr))
    assert rd.modules == {} and rd.read_level == {"state": "absent"} and rd.inputs == [] and rd.waiting == []


def test_project_of_a_sample_follows_the_root_configuration(project, tmp_path):
    fr, proj = project
    run = ident(fr)
    top = os.path.dirname(proj.project_dir)
    cands = [dataclasses.replace(c, root_id="root_a", relpath=os.path.relpath(c.path, top)) for c in project_candidates(proj, roles=("pychopper.stats",))]
    cfg = {"roots": [{"id": "root_a", "path": top, "project_depth": 1}]}
    linked = api.link(cands, [run], config=cfg, probes_dir=str(tmp_path / "probes"))
    rd = api.read_data(run, linked, None, tier1=tier1_of(fr), config=cfg)
    assert {s["project"] for s in rd.samples if s["key"] != "unclassified"} == {"project_one"}
    assert rd.samples[0]["sample_uid"] == f"project_one/{proj.aliases[0]}"
    assert [r["project"] for r in api.sample_rows(bundle_of(tier1_of(fr), rd))] == ["project_one"] * 3


@pytest.fixture(scope="module")
def two_cells(tmp_path_factory):
    """One sample sequenced on two flow cells: one alignment file and one pychopper table for both."""
    d = tmp_path_factory.mktemp("two_cells")
    fa = make_run(d / "exp", RunSpec(flow_cell_id="TST00001", sample="pool_a", seed=41))
    fb = make_run(d / "exp", RunSpec(flow_cell_id="TST00002", sample="pool_b", seed=42, start="2026-03-04T09:00:00Z"))
    A, B = T.run_ref(fa), T.run_ref(fb)
    ra, rb = T.make_reads(800, A, seed=1, alias="omega", barcode="barcode04"), T.make_reads(600, B, seed=2, alias="omega", barcode="barcode04")
    exp = str(d / "exp")
    g = os.path.join(exp, "analysis", "alignments", "omega_genome_sorted.bam")
    os.makedirs(os.path.dirname(g))
    eg = T.write_aligned_bam(g, ra + rb, seed=3)
    py = os.path.join(exp, "analysis", "pychopper_reports", "omega_pychopper_stats.tsv")
    T.write_pychopper_stats(py, pass_reads=1500, qc_fail=60, len_fail=10, primers_found=1380, rescued={2: 15}, unusable=105)
    sal = T.write_salmon_meta(os.path.join(exp, "analysis", "counts", "salmon", "omega_salmon"), 1400, 1300)
    cands = [cand("align.genome_bam", g, sample="omega", exp=exp), cand("pychopper.stats", py, sample="omega", exp=exp), cand("quant.salmon", sal, sample="omega", exp=exp)]
    return {"fa": fa, "fb": fb, "A": A, "B": B, "cands": cands, "eg": eg, "dir": d}


def test_sample_on_two_flow_cells_is_split_per_read_and_unsplittable_files_are_marked(two_cells, tmp_path):
    t = two_cells
    runs = [ident(t["fa"]), ident(t["fb"])]
    linked = api.link(t["cands"], runs, probes_dir=str(tmp_path / "probes"))
    refs = do_passes(linked, tmp_path / "partials")
    assert len(refs) == 1
    total = 0
    for run, ref, n in ((runs[0], t["A"], 800), (runs[1], t["B"], 600)):
        fr = t["fa"] if run is runs[0] else t["fb"]
        rd = api.read_data(run, linked, str(tmp_path / "partials"), tier1=tier1_of(fr), config={})
        want = t["eg"]["groups"][(ref.run_key, "omega")]
        cell, = rd.modules["len"]["cells"]
        assert cell["scalars"] == {"reads": n, "bases": want["bases"]} == {"reads": want["reads"], "bases": want["bases"]}
        total += cell["scalars"]["reads"]
        gi = next(i for i in rd.inputs if i["role"] == "align.genome_bam")
        assert gi["scope"] == "shared" and gi["shared_with"] == [runs[1].run_key if run is runs[0] else runs[0].run_key] and gi["records"] == n
        assert any(x["code"] == "shared_file_routed" for x in gi["notes"])
        py = next(i for i in rd.inputs if i["role"] == "pychopper.stats")
        assert py["scope"] == "sample_all_runs" and any(x["code"] == "covers_several_runs" for x in py["notes"])
        # the table covers both flow cells, so the per-flow-cell count of full-length reads comes from the reads of this run
        counts, tb = stage_table(rd.modules["funnel.samples"], rd.samples[0]["key"])
        assert counts["full_length"] == n and "alignment file of full-length reads" in tb["sources"]["full_length"]
        assert counts["genome_mapped"] == want["mapped"] and counts["classified"] is None
        assert counts["quantified"] is None and "merged before salmon ran" in tb["sources"]["quantified"]
        # the flow cell's own table gives the same reason, not "no source file": the file is there, it is one table for both flow cells
        _, whole = stage_table(rd.modules["funnel"])
        assert "merged before salmon ran" in whole["sources"]["quantified"], whole["sources"]
        assert rd.modules["cdna"]["cells"][0]["scalars"]["pass_reads"] == 1500 and any(x["code"] == "covers_several_runs" for x in rd.modules["cdna"]["notes"])
    assert total == t["eg"]["primary"]                                      # the two runs together are the file's primary records, once


def test_one_file_of_two_flow_cells_with_an_alias_per_sheet_goes_to_the_sample_of_each_sheet(tmp_path):
    """One sample on two flow cells in one alignment file named after the first alias; the sample sheet of each flow
    cell names the sample its own way and every read carries the alias of its own sheet. Each run must show the
    samples of its sheet, with the reads under them: never a second sample named after the file."""
    d = tmp_path
    fa = make_run(d / "exp", RunSpec(flow_cell_id="TST00001", sample="pool_a", seed=51))
    fb = make_run(d / "exp", RunSpec(flow_cell_id="TST00002", sample="pool_b", seed=52, start="2026-03-04T09:00:00Z"))
    A, B = T.run_ref(fa), T.run_ref(fb)
    exp = str(d / "exp")
    reads = T.make_reads(700, A, seed=1, alias="omega_LR", barcode="barcode02") + T.make_reads(500, B, seed=2, alias="omega_LR2", barcode="barcode02") \
        + T.make_reads(300, A, seed=3, alias="sigma_LR", barcode="barcode03") + T.make_reads(200, B, seed=4, alias="sigma_LR2", barcode="barcode03")
    by = {"omega_LR": [r for r in reads if r.alias in ("omega_LR", "omega_LR2")], "sigma_LR": [r for r in reads if r.alias in ("sigma_LR", "sigma_LR2")]}
    os.makedirs(os.path.join(exp, "analysis", "alignments"))
    cands, eg = [], {}
    for name, rs in by.items():
        g = os.path.join(exp, "analysis", "alignments", f"{name}_genome_sorted.bam")
        eg[name] = T.write_aligned_bam(g, rs, seed=3)
        cands.append(cand("align.genome_bam", g, sample=name, exp=exp))
        sal = T.write_salmon_meta(os.path.join(exp, "analysis", "counts", "salmon", f"{name}_salmon"), 900 + len(rs), 800)    # two files, two contents
        cands.append(cand("quant.salmon", sal, sample=name, exp=exp))
    for ref, suffix in ((A, ""), (B, "2")):
        rows = [{"experiment_id": ref.experiment, "kit": T.KIT, "flow_cell_id": ref.flow_cell_id, "barcode": bc, "sample_id": al + suffix, "alias": al + suffix}
                for bc, al in (("barcode02", "omega_LR"), ("barcode03", "sigma_LR"))]
        cands.append(cand("basecaller.sample_sheet", T.write_sample_sheet(os.path.join(exp, f"sheet{suffix or 1}.csv"), rows), rule="told by its content", exp=exp))
    runs = [ident(fa), ident(fb)]
    linked = api.link(cands, runs, probes_dir=str(tmp_path / "probes"))
    shared = [li for li in linked if li.candidate.role == "align.genome_bam"]
    assert all(li.scope == "shared" and sorted(li.run_keys) == sorted(r.run_key for r in runs) for li in shared)

    def names(rd):
        return [(s["name"], s["barcode"], s["source"]) for s in rd.samples]

    want_names = {fa.run_key: [("omega_LR", "barcode02", "basecaller_sample_sheet"), ("sigma_LR", "barcode03", "basecaller_sample_sheet")],
                  fb.run_key: [("omega_LR2", "barcode02", "basecaller_sample_sheet"), ("sigma_LR2", "barcode03", "basecaller_sample_sheet")]}
    # before any file is read: the link's probe already saw the label of each run's reads
    for run, fr in ((runs[0], fa), (runs[1], fb)):
        rd = api.read_data(run, linked, str(tmp_path / "partials"), tier1=tier1_of(fr), config={})
        assert names(rd) == want_names[run.run_key], names(rd)
        assert rd.modules["aln"]["state"] == "partial" and rd.modules["aln"]["partial_reason"] == "waiting"
    do_passes(linked, tmp_path / "partials")
    total = 0
    for run, fr, ref, suffix in ((runs[0], fa, A, ""), (runs[1], fb, B, "2")):
        rd = api.read_data(run, linked, str(tmp_path / "partials"), tier1=tier1_of(fr), config={})
        assert names(rd) == want_names[run.run_key], names(rd)
        keys = {s["name"]: s["key"] for s in rd.samples}
        cells = {c["at"][0]: c["scalars"] for c in rd.modules["aln"]["cells"]}
        assert set(cells) == set(keys.values())
        for name in ("omega_LR", "sigma_LR"):
            want = eg[name]["groups"][(ref.run_key, name + suffix)]
            got = cells[keys[name + suffix]]
            assert (got["reads"], got["mapped_reads"]) == (want["reads"], want["mapped"])
            total += got["reads"]
        # a table of the same file name follows the alignment file to the sheet's sample
        sal = [i for i in rd.inputs if i["role"] == "quant.salmon"]
        assert sorted(s for i in sal for s in i["samples"]) == sorted(keys.values())
        assert rd.read_level["samples_expected"] == 2 and rd.read_level["samples_with_data"] == 2
        said = [w for w in rd.warnings if w["code"] == "file_sample_is_sheet_alias"]
        assert len(said) == (2 if suffix else 0) and not any(w["code"] == "sample_not_in_sheet" for w in rd.warnings)
    assert total == sum(e["primary"] for e in eg.values())
    # a bundle put together before this rule held a second, file-named sample for each file: assembling again from that
    # bundle (a scan hands the samples it already has to the next assembly) must not keep them
    t1 = tier1_of(fb)
    t1["samples"] = [{"key": "s1", "barcode": "barcode02", "alias": "omega_LR2", "name": "omega_LR2", "sample_uid": "omega_LR2", "type": "test_sample", "expected": True,
                      "source": "basecaller_sample_sheet", "project": None},
                     {"key": "s2", "barcode": "barcode03", "alias": "sigma_LR2", "name": "sigma_LR2", "sample_uid": "sigma_LR2", "type": "test_sample", "expected": True,
                      "source": "basecaller_sample_sheet", "project": None}] \
        + [{"key": f"s{3 + k}", "barcode": None, "alias": n, "name": n, "sample_uid": n, "type": None, "expected": False, "source": "file_name", "project": None}
           for k, n in enumerate(("omega_LR", "sigma_LR"))]
    rd = api.read_data(runs[1], linked, str(tmp_path / "partials"), tier1=t1, config={})
    assert names(rd) == want_names[fb.run_key], names(rd)


def test_two_file_names_with_one_read_label_are_never_added_up_under_the_sheet_sample(tmp_path):
    """Two alignment files named differently, neither as the sheet names the sample, both holding reads labelled
    with the sheet's alias: only the first stands for the sheet's sample; the second keeps its own name."""
    fa = make_run(tmp_path / "exp", RunSpec(flow_cell_id="TST00001", sample="pool_a", seed=61))
    A = T.run_ref(fa)
    exp = str(tmp_path / "exp")
    os.makedirs(os.path.join(exp, "analysis", "alignments"))
    cands = []
    for k, name in enumerate(("first_copy", "second_copy")):
        g = os.path.join(exp, "analysis", "alignments", f"{name}_genome_sorted.bam")
        T.write_aligned_bam(g, T.make_reads(300 + 50 * k, A, seed=5 + k, alias="omega", barcode="barcode02"), seed=3 + k)
        cands.append(cand("align.genome_bam", g, sample=name, exp=exp))
    rows = [{"experiment_id": A.experiment, "kit": T.KIT, "flow_cell_id": A.flow_cell_id, "barcode": "barcode02", "sample_id": "omega", "alias": "omega"}]
    cands.append(cand("basecaller.sample_sheet", T.write_sample_sheet(os.path.join(exp, "sheet.csv"), rows), rule="told by its content", exp=exp))
    run = ident(fa)
    linked = api.link(cands, [run], probes_dir=str(tmp_path / "probes"))
    do_passes(linked, tmp_path / "partials")
    rd = api.read_data(run, linked, str(tmp_path / "partials"), tier1=tier1_of(fa), config={})
    reads = {s["name"]: next(c["scalars"]["reads"] for c in rd.modules["aln"]["cells"] if c["at"][0] == s["key"]) for s in rd.samples if s["key"] != "unclassified"}
    assert reads == {"omega": 300, "second_copy": 350}, reads


def test_thirteen_samples_are_stored_at_the_coarse_level(tmp_path):
    A = T.fake_run_ref("TST00007", 7)
    reads = []
    for i in range(13):
        reads += T.make_reads(60, A, seed=100 + i, alias=f"smp_{i:02d}", barcode=f"barcode{i + 1:02d}")
    p = tmp_path / "calls.bam"                                               # one file, every barcode: the label is on each read
    exp = T.write_dorado_bam(str(p), reads)
    run = ident_of(A)
    linked = api.link([cand("dorado.bam", str(p))], [run], probes_dir=str(tmp_path / "probes"))
    assert linked[0].sample is None and len(linked[0].link["labels"]) == 13
    do_passes(linked, tmp_path / "partials")
    fr = type("F", (), {"run_key": A.run_key, "expected": {"reads": 5000}})()
    rd = api.read_data(run, linked, str(tmp_path / "partials"), tier1=tier1_of(fr), config={})
    ln, ql = rd.modules["len"], rd.modules["qual"]
    assert ln["detail"] == ql["detail"] == "coarse" and len(ln["cells"]) == 13
    assert {c["hists"]["len_hist"]["axis"] for c in ln["cells"]} == {"read_length_log20"} and {c["hists"]["qual_hist"]["axis"] for c in ql["cells"]} == {"qscore_1"}
    assert all("hists2d" not in c for c in ql["cells"])
    grid = ql["run_level"]["hists2d"]["qual_vs_length"]                       # what the per-sample cells drop is kept once for the run
    assert sum(grid["n"]) + grid["outside"] == 13 * 60
    names = {s["key"]: s["name"] for s in rd.samples}
    for c in ln["cells"]:
        want = exp["groups"][(A.run_key, names[c["at"][0]])]
        assert c["scalars"] == {"reads": want["reads"], "bases": want["bases"]}
        assert derive.total(c["hists"]["len_hist"], "w") == want["bases"] and derive.total(c["hists"]["len_hist"], "n") == 60     # coarsening is exact
    assert all(s["source"] == "bam_read_group" and not s["expected"] for s in rd.samples)
    assert ln["population"]["stage"] == "classified"
    assert [e.message for e in defs("section").iter_errors(ql)] == []


def test_aligned_file_with_mapped_reads_only_does_not_judge_the_mapped_share(tmp_path):
    A = T.fake_run_ref("TST00008", 8)
    p = tmp_path / "mapped_only_genome_sorted.bam"
    T.write_aligned_bam(str(p), T.make_reads(300, A, seed=3), keep_unmapped=False, pychopper_names=False, full_length_input=False)
    run = ident_of(A)
    linked = api.link([cand("align.genome_bam", str(p), sample="solo")], [run], probes_dir=str(tmp_path / "probes"))
    do_passes(linked, tmp_path / "partials")
    fr = type("F", (), {"run_key": A.run_key, "expected": {"reads": 5000}})()
    t1 = tier1_of(fr, group="gdna", barcoded=False)
    rd = api.read_data(run, linked, str(tmp_path / "partials"), tier1=t1, config={})
    assert rd.modules["aln"]["population"]["stage"] == "aligned_primary" == rd.modules["len"]["population"]["stage"]
    assert rd.modules["aln"]["summary"]["aln_mapped_frac"] == 1.0
    rules = {r["id"]: r for r in api.read_checks(bundle_of(t1, rd))}
    assert rules["aln_mapped_frac"]["outcome"] == "not_evaluated" and "mapped reads only" in rules["aln_mapped_frac"]["message"]
    assert rules["aln_identity_hist"]["outcome"] in ("ok", "triggered") and rules["aln_identity_hist"]["population"]["stage"] == "aligned_primary"
    counts, t = stage_table(rd.modules["funnel"])
    assert list(counts) == ["sequenced", "basecalled", "genome_mapped"]       # a genomic run without barcodes has no other stage
    assert rd.modules["funnel"]["summary"]["funnel_usable_stage"] == "genome_mapped" and rd.modules["funnel"]["summary"]["funnel_usable_frac"] == round(counts["genome_mapped"] / 5000, 6)


def test_truncated_input_makes_sections_partial_and_checks_not_assessed(tmp_path):
    A = T.fake_run_ref("TST00006", 6)
    d = tmp_path / "bam_pass" / "alpha"
    d.mkdir(parents=True)
    p = d / f"{A.flow_cell_id}_pass_alpha_{A.run8}_{A.acq8}_0.bam"
    T.write_dorado_bam(str(p), T.make_reads(2500, A, seed=2, alias="alpha", barcode="barcode01"), truncate_frac=0.5)
    run = ident_of(A)
    linked = api.link([cand("dorado.bam", str(p), sample="alpha")], [run], probes_dir=str(tmp_path / "probes"))
    assert linked[0].run_keys == [A.run_key]
    refs = do_passes(linked, tmp_path / "partials")
    ref = refs[str(p)]
    assert ref.state == "truncated" and not ref.complete
    assert api.read_pass(linked[0], str(tmp_path / "partials"), []).cached      # cut short stays cut short while the file is the same
    fr = type("F", (), {"run_key": A.run_key, "expected": {"reads": 5000}})()
    rd = api.read_data(run, linked, str(tmp_path / "partials"), tier1=tier1_of(fr), config={})
    ln = rd.modules["len"]
    assert ln["state"] == "partial" and ln["partial_reason"] == "partly_read" and ln["coverage"]["kind"] == "truncated"
    assert 0 < ln["cells"][0]["scalars"]["reads"] < 2500
    assert rd.inputs[0]["state"] == "truncated" and rd.inputs[0]["coverage"]["kind"] == "truncated" and rd.read_level["state"] == "partial"
    assert rd.read_level["coverage"]["kind"] == "truncated"
    b = bundle_of(tier1_of(fr), rd)
    assert all(r["outcome"] == "not_evaluated" for r in api.read_checks(b) if r["module"] in ("len", "qual", "yield.basecalled"))
    assert api.card_columns(b)["t2_n50"] is None
    counts, t = stage_table(rd.modules["funnel"])
    assert counts["basecalled"] is None and counts["classified"] is None


def test_time_limited_pass_is_cached_with_its_coverage_and_done_again(project, tmp_path):
    import time
    fr, proj = project
    run = ident(fr)
    g = proj.files["align.genome_bam"][proj.aliases[0]]
    li, = api.link([cand("align.genome_bam", g, sample=proj.aliases[0])], [run], probes_dir=str(tmp_path / "probes"))
    ref = api.read_pass(li, str(tmp_path / "partials"), [], PassOptions(deadline=time.monotonic() - 5))
    assert ref.state == "partially_read" and not ref.complete and not ref.cached and os.path.isfile(ref.path)
    again = api.read_pass(li, str(tmp_path / "partials"), [], PassOptions())
    assert again.complete and not again.cached and again.path == ref.path     # the complete pass replaced it
    assert api.read_pass(li, str(tmp_path / "partials"), [], PassOptions()).cached


def test_without_pysam_collected_files_work_and_a_bam_names_the_extra(project, tmp_path, monkeypatch):
    import builtins

    from flowcellar.util.errors import MissingExtraError
    fr, proj = project
    run = ident(fr)
    real_import = builtins.__import__

    def no_pysam(name, *a, **k):
        if name == "pysam":
            raise ImportError("No module named 'pysam'")
        return real_import(name, *a, **k)
    monkeypatch.setattr(builtins, "__import__", no_pysam)
    with pytest.raises(MissingExtraError) as e:
        bam.pysam_module()
    assert "pysam" in str(e.value) and e.value.hint == 'pip install "flowcellar[bam]"'
    cands = project_candidates(proj, roles=("pychopper.stats", "quant.salmon", "dorado.bam"))
    linked = api.link(cands, [run], probes_dir=str(tmp_path / "probes"))
    by_role = {}
    for li in linked:
        by_role.setdefault(li.candidate.role, []).append(li)
    assert all(li.run_keys == [] and "pysam" in li.link["evidence"] for li in by_role["dorado.bam"])          # one clear line, no crash
    assert all(li.run_keys == [run.run_key] for li in by_role["pychopper.stats"] + by_role["quant.salmon"])    # through the sample sheet
    rd = api.read_data(run, linked, None, tier1=tier1_of(fr), config={})
    assert rd.modules["cdna"]["state"] == "ok" and len(rd.modules["cdna"]["cells"]) == 3
    counts, _ = stage_table(rd.modules["funnel"])
    assert counts["full_length"] == sum(s["full_length"] for s in proj.expected["samples"].values()) and counts["quantified"] is not None


def test_stored_copy_of_inputs_survives_json(project_data):
    import json
    d = project_data
    again = json.loads(json.dumps([dataclasses.asdict(li) for li in d["linked"]]))       # the library keeps links as JSON
    from flowcellar.contracts import Candidate, LinkedInput
    back = [LinkedInput(**{**x, "candidate": Candidate(**x["candidate"])}) for x in again]
    rd = api.read_data(d["run"], back, str(d["lib"] / "partials"), tier1=copy.deepcopy(d["tier1"]), config={})
    assert rd.modules == d["rd"].modules and rd.inputs == d["rd"].inputs


def test_stray_barcode_keeps_its_counts_and_loses_its_distributions(tmp_path):
    """A barcode that the sample sheet does not list and that holds under 0.1 % of the reads gets scalars only."""
    A = T.fake_run_ref("TST00005", 5)
    reads = T.make_reads(3000, A, seed=1, alias="kept_one", barcode="barcode01") + T.make_reads(2, A, seed=2, barcode="barcode77")
    p = tmp_path / "calls.bam"
    T.write_dorado_bam(str(p), reads)
    sheet = T.write_sample_sheet(str(tmp_path / "sheet.csv"), [{"experiment_id": "demo_experiment", "kit": T.KIT, "flow_cell_id": A.flow_cell_id,
                                                               "barcode": "barcode01", "sample_id": "kept_one", "alias": "kept_one"}])
    run = ident_of(A)
    linked = api.link([cand("dorado.bam", str(p)), cand("basecaller.sample_sheet", sheet)], [run], probes_dir=str(tmp_path / "probes"))
    do_passes(linked, tmp_path / "partials")
    fr = type("F", (), {"run_key": A.run_key, "expected": {"reads": 5000}})()
    rd = api.read_data(run, linked, str(tmp_path / "partials"), tier1=tier1_of(fr), config={})
    stray = next(s for s in rd.samples if s["name"] == "barcode77")
    assert not stray["expected"] and stray["source"] == "bam_read_group"
    cells = {c["at"][0]: c for c in rd.modules["len"]["cells"]}
    assert cells[stray["key"]] == {"at": [stray["key"], "pass"], "scalars": {"reads": 2, "bases": sum(r.length for r in reads[-2:])}}
    assert "hists" in cells["s1"] and sum(c["scalars"]["reads"] for c in cells.values()) == 3002         # nothing leaves the total
    assert any(w["code"] == "samples_counts_only" for w in rd.warnings) and any(w["code"] == "sample_not_in_sheet" for w in rd.warnings)
    assert [e.message for e in defs("section").iter_errors(rd.modules["qual"])] == []
    rules = {r["id"]: r for r in api.read_checks(bundle_of(tier1_of(fr), rd))}
    assert rules["len_n50"]["outcome"] in ("ok", "triggered") and "kept_one" not in rules["len_n50"]["message"]       # one judged sample
