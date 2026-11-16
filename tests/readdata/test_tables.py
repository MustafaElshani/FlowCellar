"""Table sources: pychopper statistics, quantification metadata, sample sheets, read summaries, and telling a file
by its first bytes."""
from __future__ import annotations

import gzip
import json
import os

import pytest
from rd_helpers import cand, ident_of

from flowcellar.model import partial as pm
from flowcellar.readdata import api, tables
from flowcellar.readdata.sample_sheet import read_sample_sheet
from flowcellar.testing import readdata as T
from flowcellar.util.errors import InputError

A = T.fake_run_ref("TST00001", 1)
B = T.fake_run_ref("TST00002", 2)


def test_pychopper_statistics(tmp_path):
    p = str(tmp_path / "s_pychopper_stats.tsv")
    exp = T.write_pychopper_stats(p, pass_reads=1000, qc_fail=30, len_fail=9, primers_found=930, rescued={2: 12, 3: 3}, unusable=55, umi=400)
    d = tables.read_pychopper(p)
    sc = d["scalars"]
    assert sc == {"pass_reads": 1000, "len_fail": 9, "qc_fail": 30, "primers_found": 930, "rescue_segments": 33, "unusable": 55,
                  "rescued_reads": 15, "umi_detected": 400, "input": 1030, "written": 930 + 33 - 9}
    assert sc["input"] == exp["input"] and sc["written"] == exp["written"] and d["consistent"]
    assert d["strand"] == {"+": 465, "-": 465} and sum(d["strand"].values()) == sc["primers_found"]
    assert d["primer_configs"][0][1] >= d["primer_configs"][-1][1] and d["cutoff_q"] == 0.8
    assert d["segments_per_rescued_read"] == {"2": 12, "3": 3}


def test_pychopper_table_that_is_not_one_is_refused(tmp_path):
    p = tmp_path / "x.tsv"
    p.write_text("a\tb\tc\n1\t2\t3\n")
    with pytest.raises(InputError):
        tables.read_pychopper(str(p))


def test_quantification_metadata(tmp_path):
    s = T.write_salmon_meta(str(tmp_path / "x_salmon"), 900, 880)
    o = T.write_oarfish_meta(str(tmp_path / "oar" / "x.meta_info.json"), 870, discard_score=4, discard_supp=20)
    k = T.write_kallisto_info(str(tmp_path / "kal" / "x"), 900, 640)
    ds, do, dk = tables.read_quant("quant.salmon", s), tables.read_quant("quant.oarfish", o), tables.read_quant("quant.kallisto", k)
    assert (ds["tool"], ds["processed"], ds["counted"], ds["version"]) == ("salmon", 900, 880, "1.10.3")
    assert (do["tool"], do["processed"], do["counted"]) == ("oarfish", None, 870) and do["discards"]["discard_supp"] == 20
    assert (dk["tool"], dk["processed"], dk["counted"]) == ("kallisto", 900, 640)
    with pytest.raises(InputError):
        tables.read_quant("quant.salmon", k)


def test_sample_sheet_in_two_column_orders(tmp_path):
    rows = [{"experiment_id": "demo_experiment", "kit": T.KIT, "flow_cell_id": "TST00001", "barcode": f"barcode{i:02d}", "sample_id": f"s{i}", "alias": f"alias_{i}"}
            for i in (1, 2, 3)]
    a = T.write_sample_sheet(str(tmp_path / "one.csv"), rows)
    b = T.write_sample_sheet(str(tmp_path / "two.csv"), rows, columns=("alias", "barcode", "kit", "flow_cell_id", "experiment_id", "sample_id"))
    sa, sb = read_sample_sheet(a), read_sample_sheet(b)
    assert sa["rows"] == sb["rows"] and sa["flow_cells"] == ["TST00001"] and sa["experiments"] == ["demo_experiment"]
    assert [r["alias"] for r in sa["rows"]] == ["alias_1", "alias_2", "alias_3"]
    assert api.sniff(a) == api.sniff(b) == "basecaller.sample_sheet"
    other = tmp_path / "config.csv"
    other.write_text("name,label,path\nx,y,z\n")
    assert api.sniff(str(other)) is None
    with pytest.raises(InputError):
        read_sample_sheet(str(other))


def test_sniff_tells_files_by_content(tmp_path):
    reads = T.make_reads(50, A, seed=3, alias="alpha")
    d, g = tmp_path / "anything.bam", tmp_path / "whatever.bam"
    T.write_dorado_bam(str(d), reads)
    T.write_aligned_bam(str(g), reads, index=None)
    tx = tmp_path / "tx.bam"
    T.write_aligned_bam(str(tx), reads, index=None, contigs=(("ENST00000000001.1|ENSG1|x", 5_000_000), ("ENST00000000002.1|ENSG2|y", 3_000_000), ("ENST00000000003.1|ENSG3|z", 20_000)))
    assert api.sniff(str(d)) == "dorado.bam" and api.sniff(str(g)) == "align.genome_bam" and api.sniff(str(tx)) == "align.transcriptome_bam"
    py = str(tmp_path / "p.tsv")
    T.write_pychopper_stats(py, pass_reads=10, qc_fail=1, len_fail=0, primers_found=8, rescued={2: 1}, unusable=1)
    assert api.sniff(py) == "pychopper.stats"
    assert api.sniff(T.write_salmon_meta(str(tmp_path / "s"), 5, 4)) == "quant.salmon"
    assert api.sniff(T.write_oarfish_meta(str(tmp_path / "o.meta_info.json"), 4)) == "quant.oarfish"
    assert api.sniff(T.write_kallisto_info(str(tmp_path / "k"), 5, 3)) == "quant.kallisto"
    sm = str(tmp_path / "sum.txt")
    T.write_dorado_summary(sm, reads)
    assert api.sniff(sm) == "dorado.summary"
    part = tmp_path / "x.fcp.json.gz"
    part.write_bytes(gzip.compress(json.dumps({"format": "flowcellar.partial", "schema_version": "0.1.1"}).encode()))
    assert api.sniff(str(part)) == "flowcellar.partial"
    junk = tmp_path / "junk.bin"
    junk.write_bytes(os.urandom(300))
    assert api.sniff(str(junk)) is None and api.sniff(str(tmp_path / "missing")) is None


def test_default_rules_are_generic():
    rules = api.default_input_rules()
    assert any(r.role == "dorado.bam" and r.glob == "**/bam_pass/{sample}/*.bam" for r in rules)
    assert any(r.role == "basecaller.sample_sheet" and r.sniff for r in rules)
    assert all(r.source == "default" for r in rules)
    text = " ".join((r.glob or "") + (r.sniff or "") for r in rules)
    assert "analysis" not in text and "gnome" not in text           # no laboratory's folder or file names


def test_read_summary_pass_is_exact_and_splits_runs_and_filters(tmp_path):
    reads = T.make_reads(700, A, seed=4, alias="alpha", barcode="barcode01") + T.make_reads(500, B, seed=5, alias="beta", barcode="barcode02")
    p = str(tmp_path / "summary.txt")
    exp = T.write_dorado_summary(p, reads, pass_min_q=12.0)
    c = cand("dorado.summary", p)
    li, = api.link([c], [ident_of(A), ident_of(B)], config={}, probes_dir=str(tmp_path / "probes"))
    assert li.link["method"] == "summary_column" and li.link["confidence"] == "certain" and sorted(li.run_keys) == sorted([A.run_key, B.run_key])
    assert li.scope == "shared"
    ref = api.read_pass(li, str(tmp_path / "partials"), [])
    assert ref.complete and not ref.cached
    part = pm.read_partial(ref.path)
    assert part["routing"]["runs"] == {A.run_key: 700, B.run_key: 500} and part["coverage"]["records"] == 1200
    for (key, label), g in exp["groups"].items():
        fr = next(f for f in part["fragments"] if f["run_key"] == key and f["sample"] == label and f["module"] == "len")
        assert fr["dims"] == ["filter"]
        by = {c["at"][0]: c["scalars"]["reads"] for c in fr["cells"]}
        assert by == {k: g[k] for k in ("pass", "fail") if g[k]}
        assert sum(c["scalars"]["bases"] for c in fr["cells"]) == g["bases"]
    again = api.read_pass(li, str(tmp_path / "partials"), [])
    assert again.cached and again.partial_id == ref.partial_id


def test_summary_without_basecalled_columns_is_an_error_not_a_crash(tmp_path):
    p = tmp_path / "sequencing_summary_x.txt"
    p.write_text("filename_pod5\tread_id\trun_id\tchannel\n" + "a.pod5\tr1\tabc\t1\n" * 5)
    assert api.sniff(str(p)) is None
    c = cand("dorado.summary", str(p))
    li, = api.link([c], [ident_of(A)], config={"link": {"rules": [{"name": "all", "match": "*", "run": A.run_key}]}}, probes_dir=None)
    li.counted_for = {A.run_key: ["reads"]}
    ref = api.read_pass(li, str(tmp_path / "partials"), [])
    assert ref.state == "error" and not ref.complete and ref.path == ""


def test_sniff_never_raises_on_arbitrary_bytes(tmp_path):
    import random
    rnd = random.Random(7)
    for i in range(300):
        p = tmp_path / f"junk_{i}"
        body = bytes(rnd.randrange(256) for _ in range(rnd.randrange(1, 400)))
        if i % 3 == 0:
            body = b"a,b\r,c\n" + body                # the start of a table, then noise
        if i % 7 == 0:
            body = b"\x1f\x8b\x08\x04" + body          # the start of a compressed block, then noise
        p.write_bytes(body)
        assert api.sniff(str(p)) in (None, "basecaller.sample_sheet", "pychopper.stats", "dorado.summary")
    empty = tmp_path / "empty"
    empty.write_bytes(b"")
    assert api.sniff(str(empty)) is None and api.sniff(str(tmp_path)) is None
