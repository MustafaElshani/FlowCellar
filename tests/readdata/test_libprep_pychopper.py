"""The whole pychopper statistics table: every category against the true values of the fixture, exactly."""
from __future__ import annotations

import dataclasses

import jsonschema
import pytest
from rd_helpers import do_passes, ident, project_candidates, tier1_of

from flowcellar.aggregates import build, merge
from flowcellar.model import axes, schema
from flowcellar.readdata import api, tables
from flowcellar.readdata import pychopper as P
from flowcellar.testing import libprep as F
from flowcellar.testing import readdata as T
from flowcellar.testing.fixtures import RunSpec, make_run
from flowcellar.util.errors import InputError


def _section_validator():
    sch = schema.load_schema()
    return jsonschema.Draft202012Validator({"$ref": "#/$defs/section", "$defs": sch["$defs"]})


def _entry(key: str, path: str, **kw) -> dict:
    return {"key": key, "name": key.upper(), "input": "rd" + key[1:], "stats": P.read_stats(path), **kw}


def test_every_category_is_read_exactly(tmp_path):
    p = str(tmp_path / "s_pychopper_stats.tsv")
    truth = F.write_pychopper_stats_full(p)
    st = P.read_stats(p)
    assert st["scalars"] == truth["scalars"]
    for k in ("cutoff_q", "segments", "unusable_hits", "rescue_hits", "trimmed", "configs", "autotune"):
        assert st[k] == truth[k], k
    assert st["categories"] == list(P.KNOWN_CATEGORIES) and st["unknown_categories"] == {}
    assert all(v is True for v in st["identities"].values()) and st["failed"] == []
    assert (st["rows_skipped"], st["rows_repeated"]) == (0, 0)
    # the reduction read fate uses and the whole table agree on every number they share
    old = tables.read_pychopper(p)["scalars"]
    assert {k: st["scalars"][k] for k in old} == old


@pytest.mark.parametrize("kw", [{"whole_numbers": True}, {"crlf": True}, {"extra_category": True}])
def test_layout_variants_give_the_same_numbers(tmp_path, kw):
    a, b = str(tmp_path / "a.tsv"), str(tmp_path / "b.tsv")
    F.write_pychopper_stats_full(a)
    F.write_pychopper_stats_full(b, **kw)
    sa, sb = P.read_stats(a), P.read_stats(b)
    for k in ("scalars", "cutoff_q", "segments", "unusable_hits", "rescue_hits", "trimmed", "configs", "autotune", "identities"):
        assert sa[k] == sb[k], k
    if kw.get("extra_category"):
        assert sb["unknown_categories"] == {"FutureCategory": 2} and sb["categories"][-3] == "FutureCategory"


def test_table_without_umi_and_without_tuning(tmp_path):
    p = str(tmp_path / "s.tsv")
    truth = F.write_pychopper_stats_full(p, F.PychopperSpec(umi=None, cutoff_index=None, rescued={}, unusable_hits={0: 10, 3: 5}))
    st = P.read_stats(p)
    assert st["scalars"] == truth["scalars"]
    assert st["scalars"]["umi_detected"] is None and st["scalars"]["umi_detected_final"] is None and st["cutoff_q"] is None
    assert st["autotune"] == [] and st["segments"] == {} and st["scalars"]["rescued_reads"] == 0
    assert st["identities"]["cutoff_in_curve"] is None and st["failed"] == []
    cell = P.cell_of(st, P.config_cats([st]))
    assert "curves" not in cell and cell["hists"]["cdna_rescue_segment_hist"]["n"] == []


def test_a_table_whose_rows_do_not_add_up_is_stored_as_written_and_named(tmp_path):
    p = str(tmp_path / "old_pychopper_stats.tsv")
    T.write_pychopper_stats(p, pass_reads=1000, qc_fail=30, len_fail=9, primers_found=930, rescued={2: 12, 3: 3}, unusable=55)     # the reduced fixture
    st = P.read_stats(p)
    assert st["scalars"]["pass_reads"] == 1000 and st["scalars"]["rescue_segments"] == 33 and st["scalars"]["umi_detected"] is None
    assert st["identities"]["classes_add_up"] is True and st["identities"]["rescue_hits_sum"] is None      # no RescueHitNr rows: not testable
    assert set(st["failed"]) == {"trimmed_sum", "hits_sum"}
    sec = P.build_section([_entry("s1", p)], order=["s1"], expected=["s1"])
    msgs = [n["message"] for n in sec["notes"] if n["code"] == "pychopper_layout"]
    assert len(msgs) == 2 and all("does not hold" in m and m.startswith("S1: ") for m in msgs)


def test_not_a_statistics_table_is_refused(tmp_path):
    p = tmp_path / "x.tsv"
    p.write_text("a\tb\tc\n1\t2\t3\n")
    with pytest.raises(InputError):
        P.read_stats(str(p))
    q = tmp_path / "y.tsv"
    q.write_text("Category\tName\tValue\nStrand\t+\t5.0\n")
    with pytest.raises(InputError):
        P.read_stats(str(q))


def test_unreadable_rows_are_counted_and_repeated_names_added(tmp_path):
    p = str(tmp_path / "s.tsv")
    truth = F.write_pychopper_stats_full(p)
    with open(p, "a") as fh:
        fh.write("Hits\tSSP,-VNP\t5.0\nUnclassHitNr\tmany\t3.0\nStrand\t+\tnot a number\nbroken line\n")
    st = P.read_stats(p)
    assert st["rows_repeated"] == 1 and st["rows_skipped"] == 3
    assert dict(st["configs"])["SSP,-VNP"] == dict(truth["configs"])["SSP,-VNP"] + 5
    assert st["failed"] == ["hits_sum"]                                    # the added 5 reads break exactly that identity


def test_cell_aggregates_are_the_table(tmp_path):
    p = str(tmp_path / "s.tsv")
    truth = F.write_pychopper_stats_full(p)
    st = P.read_stats(p)
    cats = P.config_cats([st])
    cell = P.cell_of(st, cats)
    sc, t = cell["scalars"], truth["scalars"]
    assert {k: sc[k] for k in t} == t
    assert cell["counts"]["cdna_class_table"] == {"cats": ["full_length", "rescued", "unusable"], "n": [t["primers_found"], t["rescued_reads"], t["unusable"]], "role": "categories"}
    assert sum(cell["counts"]["cdna_class_table"]["n"]) == t["pass_reads"]
    assert cell["counts"]["cdna_strand_table"]["n"] == [t["strand_plus"], t["strand_minus"]]
    assert cell["counts"]["cdna_rescue_strand_table"]["n"] == [t["rescue_strand_plus"], t["rescue_strand_minus"]]
    # primer configurations: the twelve most common by name, the rest in one count
    table = cell["counts"]["cdna_primer_config_table"]
    assert table["cats"] == [k for k, _ in truth["configs"][:12]] and table["n"] == [v for _, v in truth["configs"][:12]]
    assert table["other"] == sum(v for _, v in truth["configs"][12:]) == sc["primer_config_other"] and table["other"] > 0
    assert sum(table["n"]) + table["other"] == t["hits_reads"] == t["pass_reads"] - t["unusable_no_hit"]
    # pieces per fused read: bin k of count_1 is the whole number k; the weight is pieces
    seg = cell["hists"]["cdna_rescue_segment_hist"]
    full, w = build.hist_full(seg), build.hist_full(seg, "w")
    assert {k: int(full[1 + k]) for k in truth["segments"]} == truth["segments"] and int(full.sum()) == t["rescued_reads"]
    assert int(w.sum()) == t["rescue_segments"] and (seg["min"], seg["max"]) == (2, 5)
    # hits per unusable read: 25 hits lie beyond the axis and are counted in "over", not lost
    uh = cell["hists"]["cdna_unusable_hits_hist"]
    full = build.hist_full(uh)
    assert uh["over"] == truth["unusable_hits"][25] == 1 and int(full.sum()) == t["unusable"] and int(full[1]) == truth["unusable_hits"][0]
    rh = build.hist_full(cell["hists"]["cdna_rescue_hits_hist"])
    assert int(rh.sum()) == t["rescue_segments"] and {k: int(rh[1 + k]) for k in truth["rescue_hits"]} == truth["rescue_hits"]
    # bases not kept: 2-base bins, 300 and more in "over"; reads and bases are exact totals
    tb = cell["hists"]["cdna_trimmed_bases_hist"]
    full, w = build.hist_full(tb), build.hist_full(tb, "w")
    assert tb["axis"] == "trimmed_bases_2" and int(full.sum()) == t["trimmed_reads"] == t["primers_found"] + t["rescued_reads"]
    assert int(w.sum()) == t["trimmed_bases"] and tb["over"] == truth["trimmed"][300] + truth["trimmed"][1234]
    assert tb["w_over"] == 300 * truth["trimmed"][300] + 1234 * truth["trimmed"][1234] and (tb["min"], tb["max"]) == (0, 1234)
    assert int(full[1 + 36]) == truth["trimmed"][72] + truth["trimmed"][73] and int(full[1 + 149]) == truth["trimmed"][299]
    # the tuning curve, and the cut-off that was chosen: stored with the curve, never among the numbers that add up
    cv = cell["curves"]["cdna_autotune_curve"]
    assert list(zip(cv["x"], cv["y"]["reads"])) == truth["autotune"] and len(cv["x"]) == 30
    assert cv["y"]["chosen"] == [0, 1] + [0] * 28 and truth["cutoff_q"] == cv["x"][1] and max(cv["y"]["reads"]) == cv["y"]["reads"][1]
    assert "cutoff_q" not in sc and all(isinstance(v, int) or v is None for v in sc.values())


def _three(tmp_path):
    specs = {"s1": F.PychopperSpec(seed=1),
             "s2": F.PychopperSpec(primers_found=20000, plus=11000, rescued={2: 90, 4: 2}, unusable_hits={0: 400, 3: 300, 5: 7}, qc_fail=900,
                                   len_fail=120, umi=(19000, 18800), seed=2, configs=F.DEFAULT_CONFIGS[::-1]),
             "s3": F.PychopperSpec(primers_found=5000, plus=2400, rescued={2: 10}, unusable_hits={0: 50, 2: 20, 3: 60}, qc_fail=10, len_fail=4,
                                   umi=(4700, 4650), cutoff_index=0, seed=3)}
    truths, paths = {}, {}
    for k, spec in specs.items():
        paths[k] = str(tmp_path / f"{k}_pychopper_stats.tsv")
        truths[k] = F.write_pychopper_stats_full(paths[k], spec)
    return specs, truths, paths


def test_section_per_sample_and_summed_per_run(tmp_path):
    _, truths, paths = _three(tmp_path)
    sec = P.build_section([_entry(k, paths[k]) for k in ("s3", "s1", "s2")], order=["s1", "s2", "s3"], expected=["s1", "s2", "s3"])
    assert sec["state"] == "ok" and "partial_reason" not in sec and sec["dims"] == ["sample"] and sec["tier"] == "assay"
    assert [c["at"] for c in sec["cells"]] == [["s1"], ["s2"], ["s3"]] and sec["inputs"] == ["rd1", "rd2", "rd3"]
    assert list(_section_validator().iter_errors(sec)) == []
    # one list of named configurations for the run: the twelve most common of the three tables together
    pooled_cfg: dict[str, int] = {}
    for t in truths.values():
        for name, n in t["configs"]:
            pooled_cfg[name] = pooled_cfg.get(name, 0) + n
    top = sorted(pooled_cfg.items(), key=lambda kv: (-kv[1], kv[0]))[:12]
    rl = sec["run_level"]
    assert rl["counts"]["cdna_primer_config_table"]["cats"] == [k for k, _ in top]
    assert rl["counts"]["cdna_primer_config_table"]["n"] == [v for _, v in top]
    assert rl["counts"]["cdna_primer_config_table"]["other"] == sum(pooled_cfg.values()) - sum(v for _, v in top)
    for c in sec["cells"]:
        t = truths[c["at"][0]]
        assert {k: c["scalars"][k] for k in t["scalars"]} == t["scalars"]
        tab = c["counts"]["cdna_primer_config_table"]
        assert tab["cats"] == [k for k, _ in top] and tab["n"] == [dict(t["configs"]).get(k, 0) for k, _ in top]
        assert sum(tab["n"]) + tab["other"] == t["scalars"]["hits_reads"]
    # the run level is the exact sum of the cells: every additive scalar, every class, every histogram bin
    additive = [k for k in truths["s1"]["scalars"] if k != "primer_configs"]
    for k in additive:
        assert rl["scalars"][k] == sum(t["scalars"][k] for t in truths.values()), k
    assert rl["scalars"]["primer_configs"] == max(t["scalars"]["primer_configs"] for t in truths.values())
    assert rl["counts"]["cdna_class_table"]["n"] == [sum(t["scalars"][k] for t in truths.values()) for k in ("primers_found", "rescued_reads", "unusable")]
    for name in ("cdna_rescue_segment_hist", "cdna_unusable_hits_hist", "cdna_rescue_hits_hist", "cdna_trimmed_bases_hist"):
        total = sum(build.hist_full(c["hists"][name]) for c in sec["cells"])
        assert (build.hist_full(rl["hists"][name]) == total).all(), name
    assert "curves" not in rl and all("cdna_autotune_curve" in c["curves"] for c in sec["cells"])
    # the same sum through the generic roll-up (what any other reader of the bundle would do)
    for order in (sec["cells"], sec["cells"][::-1]):
        again = merge.rollup([{k: v for k, v in c.items() if k != "curves"} for c in order], ["sample"], [], sec["scalar_merge"])[0]
        assert again["scalars"] == rl["scalars"] and again["hists"] == rl["hists"]
    chosen = {c["at"][0]: [x for x, f in zip(c["curves"]["cdna_autotune_curve"]["x"], c["curves"]["cdna_autotune_curve"]["y"]["chosen"]) if f] for c in sec["cells"]}
    assert chosen == {k: [t["cutoff_q"]] for k, t in truths.items()} and chosen["s3"] != chosen["s1"]
    # derived numbers, by hand
    s = sec["summary"]
    tot = {k: sum(t["scalars"][k] for t in truths.values()) for k in additive}
    assert s["cdna_pychopper_samples"] == 3
    assert s["cdna_full_length_frac"] == round((tot["primers_found"] + tot["rescued_reads"]) / tot["pass_reads"], 6)
    assert s["cdna_unusable_frac"] == round(tot["unusable"] / tot["pass_reads"], 6)
    assert s["cdna_no_primer_frac"] == round(tot["unusable_no_hit"] / tot["pass_reads"], 6)
    assert s["cdna_multi_primer_frac"] == round(tot["unusable_multi_hit"] / tot["pass_reads"], 6)
    assert s["cdna_strand_plus_frac"] == round(tot["strand_plus"] / tot["primers_found"], 6)
    assert s["cdna_qc_fail_frac"] == round(tot["qc_fail"] / tot["input"], 6) and s["cdna_len_fail_frac"] == round(tot["len_fail"] / tot["oriented"], 6)
    assert s["cdna_umi_detected_frac"] == round(tot["umi_detected"] / tot["oriented"], 6)
    assert s["cdna_umi_final_frac"] == round(tot["umi_detected_final"] / tot["oriented"], 6)
    assert s["cdna_rescue_pieces_per_read"] == round(tot["rescue_segments"] / tot["rescued_reads"], 4)
    every = sorted(b for t in truths.values() for b, n in t["trimmed"].items() for _ in range(n))
    assert s["cdna_trimmed_bases_median"] == every[(len(every) + 1) // 2 - 1]
    assert s["cdna_trimmed_bases_mean"] == round(sum(every) / len(every), 2)
    assert s["cdna_primer_configs"] == len(pooled_cfg) and s["cdna_pychopper_cutoff"] is None      # s3 chose another cut-off: no single value
    assert sec["coverage"] == {"kind": "full", "records": tot["input"], "records_total": None, "fraction": 1.0}


def test_section_says_what_is_missing_or_not_known(tmp_path):
    _, truths, paths = _three(tmp_path)
    no_umi = str(tmp_path / "s4_pychopper_stats.tsv")
    F.write_pychopper_stats_full(no_umi, F.PychopperSpec(umi=None, seed=9))
    sec = P.build_section([_entry("s1", paths["s1"]), _entry("s4", no_umi, several=2, weak=True)], order=["s1", "s2", "s4"],
                          expected=["s1", "s2", "s4"], conflicts=["two candidate sources for the cDNA classification of sample S2: a and b"])
    assert sec["state"] == "partial" and sec["partial_reason"] == "samples_missing"
    assert list(_section_validator().iter_errors(sec)) == []
    codes = [n["code"] for n in sec["notes"]]
    assert codes.count("covers_several_runs") == 1 and "weak_link" in codes and "two_candidate_sources" in codes
    # one sample ran without UMI detection: its UMI counts are not known, and so is the sum; nothing reads as zero
    assert sec["cells"][1]["scalars"]["umi_detected"] is None and sec["run_level"]["scalars"]["umi_detected"] is None
    assert sec["summary"]["cdna_umi_detected_frac"] is None and sec["summary"]["cdna_umi_final_frac"] is None
    assert sec["summary"]["cdna_pychopper_cutoff"] == truths["s1"]["cutoff_q"]        # both chose the same cut-off
    assert P.build_section([], order=[], expected=[]) is None


def test_a_table_without_some_categories_says_what_is_not_known(tmp_path):
    """A table of a pychopper version that wrote no RescueSegmentNr and no UnclassHitNr rows: the numbers read from
    those rows cannot be known. They are zero in the cell, so the section must name them (the page shows the note)."""
    full = str(tmp_path / "s1_pychopper_stats.tsv")
    F.write_pychopper_stats_full(full, F.PychopperSpec(seed=3))
    older = str(tmp_path / "s2_pychopper_stats.tsv")
    with open(full) as fh:
        lines = [ln for ln in fh if not ln.startswith(("RescueSegmentNr\t", "UnclassHitNr\t"))]
    with open(older, "w") as fh:
        fh.writelines(lines)
    st = P.read_stats(older)
    assert set(st["absent"]) == {"RescueSegmentNr", "UnclassHitNr"} and st["scalars"]["rescued_reads"] == 0 and st["scalars"]["rescue_segments"] > 0
    assert P.read_stats(full)["absent"] == {}
    sec = P.build_section([_entry("s1", full), _entry("s2", older)], order=["s1", "s2"], expected=["s1", "s2"])
    assert list(_section_validator().iter_errors(sec)) == []
    notes = [n for n in sec["notes"] if n["code"] == "pychopper_categories_absent"]
    assert len(notes) == 1 and notes[0]["level"] == "warning" and notes[0]["message"].startswith("S2: ")
    assert "RescueSegmentNr" in notes[0]["message"] and "UnclassHitNr" in notes[0]["message"] and "Not known" in notes[0]["message"]
    # a table without any rescued read has no RescueSegmentNr rows and misses nothing
    none = str(tmp_path / "s3_pychopper_stats.tsv")
    F.write_pychopper_stats_full(none, F.PychopperSpec(rescued={}, seed=4))
    assert "RescueSegmentNr" not in P.read_stats(none)["absent"]


def test_trimmed_axis_is_wide_enough_for_two_primers():
    assert axes.EDGES["trimmed_bases_2"][-1] == 300 and axes.n_bins("count_1") == 21


@pytest.fixture(scope="module")
def chain(tmp_path_factory):
    """A fixture project whose pychopper tables are complete ones, linked, read and assembled."""
    root = tmp_path_factory.mktemp("lp_pych")
    fr = make_run(root / "project_one", RunSpec(experiment="demo_experiment", sample="pool_1", seed=41))
    proj = T.make_project(str(root / "project_one"), T.run_ref(fr), seed=43, quant=False)
    truths = {}
    for i, alias in enumerate(proj.aliases):
        old = proj.expected["samples"][alias]["pychopper"]
        spec = F.PychopperSpec(primers_found=old["primers_found"], plus=old["primers_found"] // 2, rescued={2: 4, 3: 1},
                               unusable_hits={0: old["unusable"] // 3, 3: old["unusable"] - old["unusable"] // 3},
                               qc_fail=old["qc_fail"], len_fail=old["len_fail"], umi=(old["primers_found"] // 2, old["primers_found"] // 2 - 3), seed=60 + i)
        truths[alias] = F.write_pychopper_stats_full(proj.files["pychopper.stats"][alias], spec)
        assert truths[alias]["scalars"]["pass_reads"] == old["pass_reads"] and truths[alias]["scalars"]["written"] == old["written"]
    return fr, proj, truths, root


def test_whole_chain_from_files_to_section(chain, monkeypatch):
    from flowcellar.readdata import link
    monkeypatch.setenv(link.INPROCESS_ENV, "1")
    fr, proj, truths, root = chain
    run = ident(fr)
    # the sample sheet and the pychopper tables alone: these tests do not depend on any pass over a BAM
    linked = api.link(project_candidates(proj, roles=("pychopper.stats",)), [run], config={}, probes_dir=str(root / "probes"))
    do_passes(linked, root / "partials", [run.run_key])
    rd = api.read_data(run, linked, str(root / "partials"), tier1=tier1_of(fr), config={})
    sec = rd.modules["cdna.pychopper"]
    assert sec["state"] == "ok" and len(sec["cells"]) == 3 and list(_section_validator().iter_errors(sec)) == []
    name_of = {s["key"]: s["name"] for s in rd.samples}
    for c in sec["cells"]:
        t = truths[name_of[c["at"][0]]]["scalars"]
        assert {k: c["scalars"][k] for k in t} == t
    assert not [n for n in sec["notes"] if n["level"] != "info"]
    # the three classes of read fate still come from the same files and agree
    cd = {c["at"][0]: c["scalars"] for c in rd.modules["cdna"]["cells"]}
    for c in sec["cells"]:
        assert all(cd[c["at"][0]][k] == c["scalars"][k] for k in ("pass_reads", "primers_found", "rescued_reads", "unusable", "written", "input"))
    assert sec["summary"]["cdna_full_length_frac"] == rd.modules["cdna"]["summary"]["cdna_full_length_frac"]
    py = [i for i in rd.inputs if i["role"] == "pychopper.stats"]
    assert len(py) == 3 and all("cdna.pychopper" in i["used_by"] and "cdna" in i["used_by"] for i in py)
    assert sec["inputs"] == [i["id"] for i in py]
    # the index column: the full-length share of the run, the number the section's summary holds
    from rd_helpers import bundle_of
    cols = api.card_columns(bundle_of(tier1_of(fr), rd))
    assert cols["t2_full_length"] == sec["summary"]["cdna_full_length_frac"] and 0.9 < cols["t2_full_length"] < 1
    assert api.card_columns({"tiers": {"read_level": {"state": "absent"}}})["t2_full_length"] is None


def test_reduced_tables_of_round_one_still_assemble(tmp_path, monkeypatch):
    """The basic fixture project writes a reduced table: the section holds it as written and names what fails."""
    from flowcellar.readdata import link
    monkeypatch.setenv(link.INPROCESS_ENV, "1")
    fr = make_run(tmp_path / "p", RunSpec(experiment="demo_experiment", sample="pool_1", seed=44))
    proj = T.make_project(str(tmp_path / "p"), T.run_ref(fr), seed=45, quant=False, dorado=False)
    run = ident(fr)
    cands = [dataclasses.replace(c) for c in project_candidates(proj, roles=("pychopper.stats",))]
    linked = api.link(cands, [run], config={}, probes_dir=str(tmp_path / "probes"))
    do_passes(linked, tmp_path / "partials", [run.run_key])
    rd = api.read_data(run, linked, str(tmp_path / "partials"), tier1=tier1_of(fr), config={})
    sec = rd.modules["cdna.pychopper"]
    assert len(sec["cells"]) == 3 and list(_section_validator().iter_errors(sec)) == []
    # the reduced tables have no RescueHitNr, Hits or Unusable rows: named as not known, beside the identities that fail
    assert {n["code"] for n in sec["notes"] if n["level"] == "warning"} == {"pychopper_layout", "pychopper_categories_absent"}
    for c in sec["cells"]:
        cv = c["curves"]["cdna_autotune_curve"]                       # a cut-off without tuning rows: the one point, no reads counted at it
        assert (cv["x"], cv["y"]) == ([0.8], {"reads": [None], "chosen": [1]}) and c["scalars"]["umi_detected"] is not None
    assert sec["summary"]["cdna_pychopper_cutoff"] == 0.8
