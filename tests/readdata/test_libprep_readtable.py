"""The streamed pass over an annotated read table, and the sections ``rna.biotype`` and ``rna.genes`` built from it:
every stored number against the true counts of the fixture rows, exactly."""
from __future__ import annotations

import collections
import os
import time

import jsonschema
import numpy as np
import pytest
from rd_helpers import cand, do_passes, ident, project_candidates, tier1_of

from flowcellar.aggregates import build, derive, merge
from flowcellar.contracts import LinkedInput, PassOptions
from flowcellar.model import axes, schema
from flowcellar.model import partial as pm
from flowcellar.readdata import api, finder, link, readtable, rnacomp
from flowcellar.testing import libprep as F
from flowcellar.testing import readdata as T
from flowcellar.testing.fixtures import RunSpec, make_run
from flowcellar.util.fingerprint import fingerprint

ROLE = "polya.read_table"


def _validator(name: str = "section"):
    sch = schema.load_schema()
    return jsonschema.Draft202012Validator({"$ref": f"#/$defs/{name}", "$defs": sch["$defs"]})


def _input(path: str, sample: str = "alpha", runs=("TST00001_0a1b2c3d",)) -> LinkedInput:
    c = cand(ROLE, path, sample=sample)
    return LinkedInput(candidate=c, fingerprint=fingerprint(path), run_keys=list(runs), link={"method": "sample_name_chain", "confidence": "high"},
                       sample=sample, needs=[ROLE], counted_for={r: [ROLE] for r in runs})


def _cells(partial: dict) -> dict[str, dict]:
    frag = next(f for f in partial["fragments"] if f["module"] == "rna.biotype")
    assert frag["dims"] == ["read_class"]
    return {c["at"][0]: c for c in frag["cells"]}


def _gene_cell(partial: dict) -> dict:
    return next(f for f in partial["fragments"] if f["module"] == "rna.genes")["cells"][0]


def _check_biotypes(partial: dict, truth: dict) -> None:
    cells = _cells(partial)
    assert set(cells) == set(truth["biotypes"])
    for name, e in truth["biotypes"].items():
        c = cells[name]
        assert c["scalars"] == {k: e[k] for k in ("reads", "tail_estimated", "tail_zero", "tail_not_found", "tail_no_value", "tail_sum")}, name
        h = c["hists"]["rna_biotype_polya_hist"]
        full = build.hist_full(h)
        want = np.zeros(axes.n_cells("polya_1nt"), np.int64)
        for v, n in collections.Counter(e["tails"]).items():
            want[1 + v if v < 1000 else -1] += n                  # 1 nt bins; 1,000 nt and longer in "over"
        assert (full == want).all(), name
        assert h["missing"] == e["reads"] - e["tail_estimated"] and h["axis"] == "polya_1nt"


def _check_genes(partial: dict, truth: dict, names: int = 10) -> None:
    g = _gene_cell(partial)
    sc = g["scalars"]
    assert sc["reads"] == truth["primary"] == sc["reads_with_gene"] and sc["reads_no_gene"] == 0
    for k in ("genes_ge1", "genes_ge10", "genes_ge100", "transcripts_ge1", "reads_top1", "reads_top10", "reads_top100", "reads_mito_genes"):
        assert sc[k] == truth[k], k
    counts = truth["gene_counts"]
    cv = g["curves"]["rna_gene_rank_curve"]
    assert cv["x"][0] == 1 and cv["x"][-1] == len(counts) and len(cv["x"]) <= 64 and cv["x"] == sorted(set(cv["x"]))
    assert cv["y"]["reads_at_rank"] == [counts[x - 1] for x in cv["x"]]
    assert cv["y"]["reads_cum"] == [sum(counts[:x]) for x in cv["x"]] and cv["y"]["reads_cum"][-1] == truth["primary"]
    if names:
        t = g["counts"]["rna_top_gene_table"]
        assert list(zip(t["cats"], t["n"])) == truth["top_genes"][:names] and t["other"] == truth["primary"] - sum(t["n"])
    else:
        assert "counts" not in g


@pytest.mark.parametrize("tail", ["dorado", "na", "mixed"])
def test_pass_counts_every_row_exactly(tmp_path, tail):
    p = str(tmp_path / "alpha_readid_annotated.csv")
    rows = F.make_table_rows(6000, seed=7, tail=tail, not_primary=41, twice=9)
    truth = F.write_annotated_table(p, rows)
    part = readtable.run_pass(_input(p), inprocess=True)
    assert part["state"] == "ok" and part["role"] == ROLE and pm.is_complete(part)
    assert part["coverage"] == {"kind": "full", "records": 6041, "records_total": 6041, "fraction": 1.0}
    f = part["facts"]
    assert (f["records"], f["primary"], f["secondary"], f["supplementary"], f["rows_skipped"]) == (6041, 6000, truth["secondary"], truth["supplementary"], 0)
    assert truth["secondary"] + truth["supplementary"] == 41 and f["tail"]["column"] == "ePolyA_tail_length" and f["tail"]["fractional"] is False
    assert part["routing"] == {"runs": {}, "unrouted": 0, "default_run": None}            # a table names no run: the link decides
    _check_biotypes(part, truth)
    _check_genes(part, truth)
    total = {k: sum(e[k] for e in truth["biotypes"].values()) for k in ("tail_estimated", "tail_zero", "tail_not_found", "tail_no_value")}
    if tail == "dorado":
        assert total["tail_no_value"] == 0 and total["tail_zero"] > 0 and total["tail_not_found"] > 0
    elif tail == "na":
        assert total["tail_zero"] == total["tail_not_found"] == 0 and total["tail_no_value"] > 0
    else:
        assert min(total.values()) > 0
    assert sum(total.values()) == 6000


def test_older_table_without_a_tail_column(tmp_path):
    p = str(tmp_path / "old_readid_annotated.csv")
    rows = F.make_table_rows(1500, seed=8)
    truth = F.write_annotated_table(p, rows, tail_column=False)
    part = readtable.run_pass(_input(p), inprocess=True)
    assert part["state"] == "ok" and part["facts"]["tail"]["column"] is None
    _check_biotypes(part, truth)
    assert all(c["scalars"]["tail_no_value"] == c["scalars"]["reads"] and c["scalars"]["tail_estimated"] == 0 for c in _cells(part).values())
    _check_genes(part, truth)


def test_many_blocks_give_the_same_partial_as_one(tmp_path, monkeypatch):
    p = str(tmp_path / "alpha_readid_annotated.csv")
    truth = F.write_annotated_table(p, F.make_table_rows(9000, seed=9, tail="mixed", not_primary=20))
    one = readtable.run_pass(_input(p), inprocess=True)
    monkeypatch.setattr(readtable, "BLOCK_BYTES", 16 << 10)
    assert os.path.getsize(p) > 40 * (16 << 10)
    many = readtable.run_pass(_input(p), inprocess=True)
    assert many["fragments"] == one["fragments"] and many["coverage"] == one["coverage"] and many["partial_id"] == one["partial_id"]
    _check_biotypes(many, truth)
    _check_genes(many, truth)


def test_columns_in_another_order_and_other_names(tmp_path):
    p = str(tmp_path / "x.csv")
    rows = F.make_table_rows(800, seed=10)
    order = ("gene_name", "ePolyA_tail_length", "FLAG", "read_id", "transcript_type", "gene_id", "transcript_id", "mapping_quality", "length")
    truth = F.write_annotated_table(p, rows, columns=order)
    part = readtable.run_pass(_input(p), inprocess=True)
    _check_biotypes(part, truth)
    _check_genes(part, truth)
    assert api.sniff(p) == ROLE


def test_gene_names_are_stored_only_as_far_as_the_catalogue_says(tmp_path, monkeypatch):
    from flowcellar import catalogue
    assert catalogue.metric("rna_top_gene_share")["names_stored"] == 10 and readtable.names_kept() == 10
    p = str(tmp_path / "alpha_readid_annotated.csv")
    truth = F.write_annotated_table(p, F.make_table_rows(2000, seed=11))
    with_names = readtable.run_pass(_input(p), inprocess=True)
    monkeypatch.setattr(readtable, "names_kept", lambda: 0)
    without = readtable.run_pass(_input(p), inprocess=True)
    _check_genes(without, truth, names=0)
    assert "GENE0" not in str(without["fragments"]) and "MT-TEST" not in str(without["fragments"])
    assert without["partial_id"] != with_names["partial_id"]            # what a pass stores is part of its cache key
    assert without["facts"]["gene_names_kept"] == 0 and with_names["facts"]["gene_names_kept"] == 10
    # no gene id or transcript id is ever stored, with names or without
    assert "FAKEG" not in str(with_names) and "FAKET" not in str(with_names)


def test_rows_that_cannot_be_read(tmp_path):
    p = str(tmp_path / "alpha_readid_annotated.csv")
    rows = F.make_table_rows(1000, seed=12)
    truth = F.write_annotated_table(p, rows)
    with open(p, "a") as fh:
        fh.write("one,row,with,too,few,columns\n")
    part = readtable.run_pass(_input(p), inprocess=True)
    assert part["state"] == "ok" and part["facts"]["rows_skipped"] == 1 and part["coverage"]["records_total"] == 1001
    _check_biotypes(part, truth)
    with open(p, "a") as fh:
        fh.write("short,row\n" * 30)                                     # 3 % of the rows: the table is not trusted
    bad = readtable.run_pass(_input(p), inprocess=True)
    assert bad["state"] == "error" and bad["errors"][0]["code"] == "bad_rows"


def test_table_cut_off_inside_its_last_row(tmp_path):
    p = str(tmp_path / "alpha_readid_annotated.csv")
    F.write_annotated_table(p, F.make_table_rows(500, seed=13), final_newline=False)
    part = readtable.run_pass(_input(p), inprocess=True)
    assert part["state"] == "truncated" and part["coverage"]["kind"] == "truncated" and not pm.is_complete(part)
    assert "line break" in part["errors"][0]["message"]


def test_deadline_stops_the_pass_and_says_how_far_it_got(tmp_path, monkeypatch):
    p = str(tmp_path / "alpha_readid_annotated.csv")
    F.write_annotated_table(p, F.make_table_rows(9000, seed=14))
    monkeypatch.setattr(readtable, "BLOCK_BYTES", 16 << 10)
    part = readtable.run_pass(_input(p), deadline=time.monotonic() - 1.0, inprocess=True)
    assert part["state"] == "partially_read" and part["coverage"]["kind"] == "byte_ranges" and not pm.is_complete(part)
    assert 0 < part["coverage"]["records"] < 9000 and 0 < part["coverage"]["fraction"] < 1


def test_header_only_and_wrong_tables(tmp_path):
    empty = str(tmp_path / "empty_readid_annotated.csv")
    F.write_annotated_table(empty, [])
    part = readtable.run_pass(_input(empty), inprocess=True)
    assert part["state"] == "ok" and part["facts"]["records"] == 0 and _cells(part) == {} and _gene_cell(part)["scalars"]["genes_ge1"] == 0
    other = tmp_path / "other.csv"
    other.write_text("name,label,path\nx,y,z\n")
    bad = readtable.run_pass(_input(str(other)), inprocess=True)
    assert bad["state"] == "error" and "not an annotated read table" in bad["errors"][0]["message"]
    assert api.sniff(str(other)) is None and api.sniff(empty) == ROLE
    assert finder.guess_sample(ROLE, empty) == "empty"


def test_a_column_that_is_no_biotype_is_refused(tmp_path, monkeypatch):
    p = str(tmp_path / "alpha_readid_annotated.csv")
    F.write_annotated_table(p, F.make_table_rows(300, seed=15))
    monkeypatch.setattr(readtable, "MAX_BIOTYPES_SEEN", 5)
    part = readtable.run_pass(_input(p), inprocess=True)
    assert part["state"] == "error" and "not a transcript type" in part["errors"][0]["message"]


def test_very_many_biotypes_are_folded_not_dropped(tmp_path, monkeypatch):
    p = str(tmp_path / "alpha_readid_annotated.csv")
    truth = F.write_annotated_table(p, F.make_table_rows(3000, seed=16))
    monkeypatch.setattr(readtable, "MAX_BIOTYPES_KEPT", 4)
    part = readtable.run_pass(_input(p), inprocess=True)
    cells = _cells(part)
    ranked = sorted(truth["biotypes"], key=lambda b: (-truth["biotypes"][b]["reads"], b))
    assert list(cells) == ranked[:4] + [readtable.FOLDED] and part["facts"]["biotypes_folded"] == len(ranked) - 4
    assert cells[readtable.FOLDED]["scalars"]["reads"] == sum(truth["biotypes"][b]["reads"] for b in ranked[4:])
    assert sum(c["scalars"]["reads"] for c in cells.values()) == 3000


def test_pool_pass_through_a_worker_process_and_the_cache(tmp_path):
    p = str(tmp_path / "alpha_readid_annotated.csv")
    truth = F.write_annotated_table(p, F.make_table_rows(2500, seed=17, not_primary=7))
    inp = _input(p)
    inp.partial_id = pm.partial_id(inp.fingerprint, ROLE, readtable.MODULES, readtable.config_hash())
    ref = api.read_pass(inp, str(tmp_path / "partials"), ["TST00001_0a1b2c3d"], PassOptions())
    assert ref.complete and not ref.cached and ref.partial_id == inp.partial_id and ref.role == ROLE
    part = pm.read_partial(ref.path)
    _check_biotypes(part, truth)
    _check_genes(part, truth)
    again = api.read_pass(inp, str(tmp_path / "partials"), ["TST00001_0a1b2c3d"], PassOptions())
    assert again.cached and again.partial_id == ref.partial_id and again.path == ref.path


# ------------------------------------------------------------------------------------------------ sections
def _entry(key: str, part: dict, **kw) -> dict:
    return {"key": key, "name": key.upper(), "input": "rd" + key[1:], "partial": part, "state": part["state"], "size": part["input"]["size"],
            "several": None, "weak": False, **kw}


@pytest.fixture(scope="module")
def three(tmp_path_factory):
    d = tmp_path_factory.mktemp("lp_tables")
    truths, parts = {}, {}
    for i, key in enumerate(("s1", "s2", "s3")):
        p = str(d / f"{key}_readid_annotated.csv")
        truths[key] = F.write_annotated_table(p, F.make_table_rows(4000 + 1500 * i, seed=20 + i, tail="mixed" if i == 1 else "dorado", not_primary=5 * i))
        parts[key] = readtable.run_pass(_input(p, sample=key), inprocess=True)
    return truths, parts


def test_biotype_section_is_exact_per_sample_and_class(three):
    truths, parts = three
    out = rnacomp.build_sections([_entry(k, parts[k]) for k in ("s2", "s3", "s1")], order=["s1", "s2", "s3"], expected=["s1", "s2", "s3"])
    assert set(out) == {"rna.biotype", "rna.genes"}
    sec = out["rna.biotype"]
    assert sec["state"] == "ok" and sec["dims"] == ["sample", "read_class"] and sec["tier"] == "assay" and sec.get("population") is None
    assert list(_validator().iter_errors(sec)) == []
    # one list of classes for the run: the biotypes with at least 0.1 % of the reads, at most eight, then "other"
    run_total: dict[str, int] = {}
    for t in truths.values():
        for b, e in t["biotypes"].items():
            run_total[b] = run_total.get(b, 0) + e["reads"]
    reads = sum(run_total.values())
    ranked = sorted(run_total, key=lambda b: (-run_total[b], b))
    own = [b for b in ranked if run_total[b] >= 0.001 * reads][:8]
    assert sec["params"]["classes"] == own + ["other"] and len(own) == 8 and len(ranked) == 10
    cells = {tuple(c["at"]): c for c in sec["cells"]}
    for key, t in truths.items():
        for cls in own + ["other"]:
            members = [b for b in t["biotypes"] if (b if b in own else "other") == cls]
            c = cells[(key, cls)]
            for k in ("reads", "tail_estimated", "tail_zero", "tail_not_found", "tail_no_value", "tail_sum"):
                assert c["scalars"][k] == sum(t["biotypes"][b][k] for b in members), (key, cls, k)
            tails = [v for b in members for v in t["biotypes"][b]["tails"]]
            h = c["hists"]["rna_biotype_polya_hist"]
            want = np.zeros(axes.n_cells("polya_5nt"), np.int64)
            for v in tails:
                want[1 + v // 5 if v < 1000 else -1] += 1                 # 5 nt bins
            assert h["axis"] == "polya_5nt" and (build.hist_full(h) == want).all() and h["missing"] == c["scalars"]["reads"] - len(tails)
            tab = c["counts"]["rna_biotype_table"]
            assert dict(zip(tab["cats"], tab["n"])) == {b: t["biotypes"][b]["reads"] for b in members}
            assert dict(zip(tab["cats"], tab["w"])) == {b: t["biotypes"][b]["tail_estimated"] for b in members}
    assert len(cells) == 3 * 9
    # cells never overlap: summing them gives each sample's rows, and the run's
    per_sample = {c["at"][0]: c for c in merge.rollup(sec["cells"], ["sample", "read_class"], ["sample"])}
    for key, t in truths.items():
        assert per_sample[key]["scalars"]["reads"] == t["primary"]
        assert dict(zip(per_sample[key]["counts"]["rna_biotype_table"]["cats"], per_sample[key]["counts"]["rna_biotype_table"]["n"])) == \
            {b: e["reads"] for b, e in t["biotypes"].items()}
    run = merge.rollup(sec["cells"], ["sample", "read_class"], [])[0]
    assert run["scalars"]["reads"] == reads == sec["summary"]["rna_biotype_reads"] == sec["coverage"]["records"]
    # exact quartiles per sample and class, and for the samples together
    table = sec["run_level"]["tables"]["rna_polya_by_biotype"]
    assert [c["id"] for c in table["columns"]] == ["sample", "read_class", "reads", "with_tail", "median_nt", "q25_nt", "q75_nt", "mean_nt"]
    rows = {(r[0], r[1]): r for r in table["rows"]}
    assert len(rows) == len(table["rows"]) == 4 * 10

    def tails_of(keys, cls):
        return [v for k in keys for b, e in truths[k]["biotypes"].items() if cls == "all" or (b if b in own else "other") == cls for v in e["tails"]]
    for sample, keys in (("s1", ["s1"]), ("s2", ["s2"]), ("s3", ["s3"]), ("*", ["s1", "s2", "s3"])):
        for cls in ["all"] + own + ["other"]:
            tails = [min(v, 1000) for v in tails_of(keys, cls)]
            r = rows[(sample, cls)]
            assert r[3] == len(tails) and r[4] == F.lower_median(tails) and r[5] == F.lower_quantile(tails, 0.25) and r[6] == F.lower_quantile(tails, 0.75)
            real = tails_of(keys, cls)
            assert r[7] == round(sum(real) / len(real), 2)
    every = tails_of(["s1", "s2", "s3"], "all")
    s = sec["summary"]
    assert s["rna_biotype_samples"] == 3 and s["rna_biotypes_seen"] == 10
    assert s["rna_biotype_polya_median_nt"] == F.lower_median(every) and s["rna_biotype_polya_frac"] == round(len(every) / reads, 6)
    assert s["rna_protein_coding_frac"] == round(run_total["protein_coding"] / reads, 6)
    assert s["rna_rrna_frac"] == round((run_total["rRNA"] + run_total["Mt_rRNA"]) / reads, 6)
    # the median read from the stored 5 nt histogram lies within one bin of the exact one
    pooled = merge.rollup(sec["cells"], ["sample", "read_class"], [])[0]["hists"]["rna_biotype_polya_hist"]
    assert abs(derive.median(pooled) - F.lower_median(every)) < 5
    assert [n["code"] for n in sec["notes"]] == ["source_is_report", "rows_not_primary", "rows_not_primary"]


def test_gene_section_is_exact_per_sample_and_never_adds_genes_up(three):
    truths, parts = three
    sec = rnacomp.build_sections([_entry(k, parts[k]) for k in ("s1", "s2", "s3")], order=["s1", "s2", "s3"], expected=["s1", "s2", "s3"])["rna.genes"]
    assert sec["state"] == "ok" and sec["dims"] == ["sample"] and list(_validator().iter_errors(sec)) == []
    per = ("genes_ge1", "genes_ge10", "genes_ge100", "transcripts_ge1", "reads_top1", "reads_top10", "reads_top100")
    detected = sec["run_level"]["tables"]["rna_genes_detected"]
    assert [c["id"] for c in detected["columns"]] == ["sample", "reads", *per] and detected["merge"] == "replace"
    assert {r[0]: r[1:] for r in detected["rows"]} == {k: [t["primary"]] + [t[x] for x in per] for k, t in truths.items()}
    top = sec["run_level"]["tables"]["rna_top_gene_table"]
    assert [c["id"] for c in top["columns"]] == ["sample", "rank", "gene", "reads"] and len(top["rows"]) == 30
    for key, t in truths.items():
        assert [(r[2], r[3]) for r in top["rows"] if r[0] == key] == t["top_genes"]
        assert [r[1] for r in top["rows"] if r[0] == key] == list(range(1, 11))
    for c in sec["cells"]:
        t = truths[c["at"][0]]
        assert c["scalars"] == {"reads": t["primary"], "reads_with_gene": t["primary"], "reads_no_gene": 0, "reads_mito_genes": t["reads_mito_genes"]}
        assert c["curves"]["rna_gene_rank_curve"]["y"]["reads_cum"][-1] == t["primary"] and "counts" not in c
    s = sec["summary"]
    assert s["rna_genes_detected_min"] == min(t["genes_ge10"] for t in truths.values()) and s["rna_genes_detected_max"] == max(t["genes_ge10"] for t in truths.values())
    assert s["rna_top_gene_max_frac"] == round(max(t["reads_top10"] / t["primary"] for t in truths.values()), 6)
    assert s["rna_mito_frac"] == round(sum(t["reads_mito_genes"] for t in truths.values()) / sum(t["primary"] for t in truths.values()), 6)
    # the cells hold only what adds up: a roll-up over samples, in any order, gives the run's reads and nothing about genes
    for order in (sec["cells"], sec["cells"][::-1]):
        run = merge.rollup([{k: v for k, v in c.items() if k != "curves"} for c in order], ["sample"], [])[0]["scalars"]
        assert run == {"reads": sum(t["primary"] for t in truths.values()), "reads_with_gene": sum(t["primary"] for t in truths.values()),
                       "reads_no_gene": 0, "reads_mito_genes": sum(t["reads_mito_genes"] for t in truths.values())}


def test_sections_say_waiting_partly_read_and_missing(three):
    truths, parts = three
    e1 = _entry("s1", parts["s1"])
    out = rnacomp.build_sections([e1], order=["s1", "s2"], expected=["s1", "s2"], waiting=[{"input": "rd2", "size": e1["size"] * 3}])
    b = out["rna.biotype"]
    assert b["state"] == "partial" and b["partial_reason"] == "waiting" and b["coverage"]["kind"] == "byte_ranges" and b["coverage"]["fraction"] == 0.25
    assert b["inputs"] == ["rd1", "rd2"] and any(n["code"] == "waiting" for n in b["notes"]) and list(_validator().iter_errors(b)) == []
    only_waiting = rnacomp.build_sections([], order=["s1"], expected=["s1"], waiting=[{"input": "rd1", "size": 10}])
    assert only_waiting["rna.biotype"]["cells"] == [] and only_waiting["rna.biotype"]["partial_reason"] == "waiting"
    assert list(_validator().iter_errors(only_waiting["rna.biotype"])) == [] and list(_validator().iter_errors(only_waiting["rna.genes"])) == []
    missing = rnacomp.build_sections([e1], order=["s1", "s2"], expected=["s1", "s2"])
    assert missing["rna.genes"]["partial_reason"] == "samples_missing" and missing["rna.biotype"]["partial_reason"] == "samples_missing"
    cut = dict(parts["s2"], state="truncated", coverage={"kind": "truncated", "records": 10, "records_total": None, "fraction": 0.5},
               errors=[{"code": "truncated", "message": "cut"}])
    part = rnacomp.build_sections([e1, _entry("s2", cut, several=2, weak=True)], order=["s1", "s2"], expected=["s1", "s2"])["rna.biotype"]
    assert part["partial_reason"] == "partly_read" and part["coverage"] == {"kind": "truncated", "records": part["summary"]["rna_biotype_reads"], "records_total": None, "fraction": 0.5}
    assert {"partly_read", "covers_several_runs", "weak_link"} <= {n["code"] for n in part["notes"]}
    assert rnacomp.build_sections([], order=[], expected=[]) == {}


def test_tables_without_tails_give_no_tail_numbers(tmp_path):
    p = str(tmp_path / "s1_readid_annotated.csv")
    F.write_annotated_table(p, F.make_table_rows(900, seed=30), tail_column=False)
    part = readtable.run_pass(_input(p, sample="s1"), inprocess=True)
    sec = rnacomp.build_sections([_entry("s1", part)], order=["s1"], expected=["s1"])["rna.biotype"]
    assert sec["summary"]["rna_biotype_polya_frac"] is None and sec["summary"]["rna_biotype_polya_median_nt"] is None
    assert sec["summary"]["rna_biotype_reads"] == 900 and any(n["code"] == "no_tail_column" for n in sec["notes"])
    assert all(r[4] is None and r[7] is None and r[3] == 0 for r in sec["run_level"]["tables"]["rna_polya_by_biotype"]["rows"])
    assert list(_validator().iter_errors(sec)) == []


def test_quantile_rule():
    full = np.zeros(axes.n_cells("polya_1nt"), np.int64)
    for v in (5, 5, 7, 9):
        full[1 + v] += 1
    assert rnacomp.quantile(full, 0.5) == 5.0 and rnacomp.quantile(full, 0.75) == 7.0 and rnacomp.quantile(full, 1.0) == 9.0
    full[-1] += 10                                                        # ten tails of 1,000 nt or more
    assert rnacomp.quantile(full, 0.5) == 1000.0 and rnacomp.quantile(np.zeros_like(full), 0.5) is None


# ------------------------------------------------------------------------------------------------ the chain
@pytest.fixture(scope="module")
def chain(tmp_path_factory):
    root = tmp_path_factory.mktemp("lp_chain")
    fr = make_run(root / "project_one", RunSpec(experiment="demo_experiment", sample="pool_1", seed=51))
    proj = T.make_project(str(root / "project_one"), T.run_ref(fr), seed=53, quant=False, reads_per_sample=300)
    truths, paths = {}, {}
    for i, alias in enumerate(proj.aliases):
        paths[alias] = os.path.join(proj.project_dir, "analysis", "dorado_polyA", "annotated_readid", f"{alias}_readid_annotated.csv")
        truths[alias] = F.write_annotated_table(paths[alias], F.make_table_rows(1200 + 300 * i, seed=70 + i, not_primary=3))
    return fr, proj, truths, paths, root


def _chain_candidates(proj, paths, skip=()):
    """The sample sheet, the pychopper tables and the annotated read tables: no BAM, so these tests stand on the tables alone."""
    return project_candidates(proj, roles=("pychopper.stats",)) + [cand(ROLE, p, sample=a, rule=f"configured rule {ROLE}", exp=proj.project_dir) for a, p in paths.items() if a not in skip]


def test_whole_chain_link_pass_sections(chain, monkeypatch):
    monkeypatch.setenv(link.INPROCESS_ENV, "1")
    fr, proj, truths, paths, root = chain
    run = ident(fr)
    linked = api.link(_chain_candidates(proj, paths), [run], config={}, probes_dir=str(root / "probes"))
    tables = [li for li in linked if li.candidate.role == ROLE]
    assert len(tables) == 3
    for li in tables:
        assert li.run_keys == [run.run_key] and li.link["method"] == "sample_name_chain" and li.link["confidence"] == "high"
        assert li.counted_for == {run.run_key: [ROLE]} and li.not_counted_reason is None and li.needs == [ROLE]
        assert li.partial_id == pm.partial_id(li.fingerprint, ROLE, readtable.MODULES, readtable.config_hash())
        assert li.partial_id != pm.partial_id(li.fingerprint, ROLE, link.PASS_MODULES, "")       # not the cache key of the BAM passes
    # before any pass: the sections exist and say they wait
    rd0 = api.read_data(run, linked, str(root / "partials"), tier1=tier1_of(fr), config={})
    assert rd0.modules["rna.biotype"]["partial_reason"] == "waiting" and rd0.modules["rna.biotype"]["cells"] == []
    assert sorted(paths.values()) == [p for p in rd0.waiting if p.endswith("_readid_annotated.csv")]
    refs = do_passes(linked, root / "partials", [run.run_key])
    assert all(refs[p].complete and refs[p].role == ROLE for p in paths.values())
    rd = api.read_data(run, linked, str(root / "partials"), tier1=tier1_of(fr), config={})
    bio, genes = rd.modules["rna.biotype"], rd.modules["rna.genes"]
    assert bio["state"] == "ok" and genes["state"] == "ok" and not [p for p in rd.waiting if p.endswith(".csv")]
    assert list(_validator().iter_errors(bio)) == [] and list(_validator().iter_errors(genes)) == []
    name_of = {s["key"]: s["name"] for s in rd.samples}
    per_sample = {name_of[c["at"][0]]: c for c in merge.rollup(bio["cells"], ["sample", "read_class"], ["sample"])}
    for alias, t in truths.items():
        assert per_sample[alias]["scalars"]["reads"] == t["primary"]
        assert per_sample[alias]["scalars"]["tail_estimated"] == sum(e["tail_estimated"] for e in t["biotypes"].values())
    assert {name_of[r[0]]: r[3] for r in genes["run_level"]["tables"]["rna_genes_detected"]["rows"]} == {a: t["genes_ge10"] for a, t in truths.items()}
    ins = [i for i in rd.inputs if i["role"] == ROLE]
    assert len(ins) == 3 and list(_validator("input").iter_errors(ins[0])) == []
    for i in ins:
        assert i["state"] == "ok" and i["used_by"] == ["rna.biotype", "rna.genes"] and i["counted_for"] == ["collected"]
        assert i["records"] == truths[name_of[i["samples"][0]]]["primary"] and i["coverage"]["kind"] == "full"
    assert bio["inputs"] == [i["id"] for i in ins] and set(rd.partial_ids) >= {refs[p].partial_id for p in paths.values()}
    # a second round of passes reads nothing
    again = do_passes(linked, root / "partials", [run.run_key])
    assert all(again[p].cached for p in paths.values())


def test_chain_with_one_table_missing_and_one_shared_by_two_runs(chain, monkeypatch, tmp_path):
    monkeypatch.setenv(link.INPROCESS_ENV, "1")
    fr, proj, truths, paths, root = chain
    run = ident(fr)
    linked = api.link(_chain_candidates(proj, paths, skip=(proj.aliases[2],)), [run], config={}, probes_dir=str(root / "probes"))
    do_passes(linked, root / "partials", [run.run_key])
    rd = api.read_data(run, linked, str(root / "partials"), tier1=tier1_of(fr), config={})
    assert rd.modules["rna.biotype"]["partial_reason"] == "samples_missing" and rd.modules["rna.biotype"]["summary"]["rna_biotype_samples"] == 2
    # two candidate tables for one sample: neither is counted, both are named, nothing is added up
    twin = str(tmp_path / "elsewhere" / f"{proj.aliases[0]}_readid_annotated.csv")
    F.write_annotated_table(twin, F.make_table_rows(700, seed=99))
    cands = _chain_candidates(proj, paths) + [cand(ROLE, twin, sample=proj.aliases[0], rule="another rule", exp=proj.project_dir)]
    linked = api.link(cands, [run], config={}, probes_dir=str(root / "probes"))
    do_passes(linked, root / "partials", [run.run_key])
    rd = api.read_data(run, linked, str(root / "partials"), tier1=tier1_of(fr), config={})
    bio = rd.modules["rna.biotype"]
    assert bio["summary"]["rna_biotype_samples"] == 2 and bio["state"] == "partial"
    assert any(n["code"] == "two_candidate_sources" and "nothing was added up" in n["message"] for n in bio["notes"])
    assert sum(1 for i in rd.inputs if i["role"] == ROLE and i["not_counted_reason"]) == 2
