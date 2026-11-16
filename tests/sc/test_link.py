"""The evidence ladder: from an analysis to the flow cell runs it came from."""
import pytest

from flowcellar.sc import link
from flowcellar.sc.link import RunFact
from flowcellar.util.errors import InputError, UsageError

RUNS = [RunFact("TST00001_0a1b2c3d", "TST00001", "0a1b2c3d", sample_id="demo_spatial_a", paths=("/data/proj/demo_spatial_a/20260101_1000_X_TST00001_0a1b2c3d",)),
        RunFact("TST00003_aaaa1111", "TST00003", "aaaa1111", sample_id="_Demo-Spatial_B", paths=("/data/elsewhere/run_b",)),
        RunFact("TST00009_bbbb2222", "TST00009", "bbbb2222", sample_id="other", paths=("/data/proj/other/run",)),
        RunFact("TST00009_cccc3333", "TST00009", "cccc3333", sample_id="other", paths=("/data/proj/other/run2",))]


def _fn(fc, prot, acq, reads, share):
    return {"flow_cell_id": fc, "protocol_run_id_prefix": prot, "acquisition_run_id_prefix": acq, "reads": reads, "share": share,
            "st_min": None, "st_max": None}


def test_the_file_name_on_a_read_names_flow_cell_run_and_acquisition():
    m = link.FN_RE.match("TST00001_0a1b2c3d_11112222_77.pod5")
    assert (m["fc"], m["prot"], m["acq"]) == ("TST00001", "0a1b2c3d", "11112222")
    assert link.FN_RE.match("reads.pod5") is None and link.FN_RE.match("TST00001_0a1b2c3d_77.pod5") is None


def test_links_from_read_tags_with_shares_and_two_acquisitions_of_one_run():
    ev = link.Evidence(fn_runs=[_fn("TST00001", "0a1b2c3d", "11112222", 500, 0.5), _fn("TST00002", "4e5f6a7b", "33334444", 400, 0.4),
                                _fn("TST00001", "0a1b2c3d", "99990000", 100, 0.1)], records=1000)
    links = link.intrinsic_links(ev)
    assert [(ln["flow_cell_id"], ln["share"], ln["acquisition_run_id_prefix"]) for ln in links] == [("TST00001", 0.6, "11112222"), ("TST00002", 0.4, "33334444")]
    assert all(ln["method"] == "read_tag" and ln["confidence"] == "high" for ln in links)
    assert "600 of 1,000 records" in links[0]["evidence"]


def test_a_run_that_is_not_in_the_library_stays_named_and_unresolved():
    ev = link.Evidence(fn_runs=[_fn("TST00001", "0a1b2c3d", "11112222", 6, 0.6), _fn("TST00002", "4e5f6a7b", "33334444", 4, 0.4)], records=10)
    links = link.resolve(intrinsic=link.intrinsic_links(ev), runs=RUNS)
    assert [(ln["run_key"], ln["resolved"]) for ln in links] == [("TST00001_0a1b2c3d", True), (None, False)]
    assert links[1]["flow_cell_id"] == "TST00002" and links[1]["protocol_run_id_prefix"] == "4e5f6a7b" and links[1]["share"] == 0.4
    # nothing in the library: every link stays, none resolves
    assert [ln["resolved"] for ln in link.resolve(intrinsic=link.intrinsic_links(ev), runs=[])] == [False, False]


def test_read_tags_win_over_every_weaker_rung():
    ev = link.Evidence(fn_runs=[_fn("TST00002", "4e5f6a7b", "33334444", 4, 1.0)], records=4)
    links = link.resolve(intrinsic=link.intrinsic_links(ev), runs=RUNS, names=["demo_spatial_a"], analysis_path="/data/proj/demo_spatial_a/analysis/x")
    assert [ln["method"] for ln in links] == ["read_tag"]


def test_a_spaceranger_output_takes_the_links_of_the_workflow_analysis_of_the_same_reads():
    sib = [{"run_key": "TST00001_0a1b2c3d", "flow_cell_id": "TST00001", "protocol_run_id_prefix": "0a1b2c3d", "acquisition_run_id_prefix": "11112222",
            "share": 1.0, "method": "read_tag", "confidence": "high", "resolved": True, "evidence": "file names stored on the reads"}]
    links = link.resolve(intrinsic=[], runs=RUNS, sibling_links=sib, names=["other"])
    assert [(ln["method"], ln["confidence"], ln["run_key"], ln["share"]) for ln in links] == [("sibling_analysis", "high", "TST00001_0a1b2c3d", 1.0)]
    assert "the workflow analysis of the same reads" in links[0]["evidence"]


def test_read_id_sketch_gives_runs_and_shares():
    ids = [f"0000{i:012x}" for i in range(40)]
    runs = [RunFact("TST00001_0a1b2c3d", "TST00001", "0a1b2c3d", sketch=("0000", tuple(ids[:30] + ["0000ffffffffffff"]))),
            RunFact("TST00002_4e5f6a7b", "TST00002", "4e5f6a7b", sketch=("0000", tuple(ids[30:40]))),
            RunFact("TST00003_aaaa1111", "TST00003", "aaaa1111", sketch=("0000", tuple(f"0000{i:012x}" for i in range(500, 560))))]
    links = link.resolve(intrinsic=[], runs=runs, sketch={"prefix": "0000", "ids": ids})
    assert [(ln["run_key"], ln["share"], ln["method"]) for ln in links] == [("TST00001_0a1b2c3d", 0.75, "read_id_sketch"), ("TST00002_4e5f6a7b", 0.25, "read_id_sketch")]
    # a sketch taken with five zeros compares with one taken with four on the longer prefix
    five = [f"00000{i:011x}" for i in range(20)]
    runs5 = [RunFact("TST00001_0a1b2c3d", "TST00001", "0a1b2c3d", sketch=("0000", tuple(five + [f"0000a{i:011x}" for i in range(50)])))]
    assert [ln["share"] for ln in link.resolve(intrinsic=[], runs=runs5, sketch={"prefix": "00000", "ids": five})] == [1.0]
    assert link.resolve(intrinsic=[], runs=runs, sketch={"prefix": "0000", "ids": ids[:3]}) == []       # too few ids to say anything


def test_folder_and_sample_name_then_sample_name_alone():
    beside = link.resolve(intrinsic=[], runs=RUNS, names=["Demo_Spatial_A"], analysis_path="/data/proj/demo_spatial_a/analysis/wf/demo_spatial_a")
    assert [(ln["run_key"], ln["method"], ln["confidence"]) for ln in beside] == [("TST00001_0a1b2c3d", "path_rule", "medium")]
    # the same name spelled with other punctuation, and no run folder beside the analysis: by name only, low confidence
    alone = link.resolve(intrinsic=[], runs=RUNS, names=["demo_spatial_b"], analysis_path="/somewhere/else/outs")
    assert [(ln["run_key"], ln["method"], ln["confidence"]) for ln in alone] == [("TST00003_aaaa1111", "sample_name_chain", "low")]
    # two runs carry the name and none lies beside the analysis: no guess
    assert link.resolve(intrinsic=[], runs=RUNS, names=["other"], analysis_path="/somewhere/else/outs") == []
    assert link.resolve(intrinsic=[], runs=RUNS, names=["nobody"], analysis_path="/x") == []
    assert link.normalise_name("_Demo-Spatial_B") == link.normalise_name("demo_spatial_b") == "demospatialb"


def test_what_a_person_states_wins_and_an_ambiguous_flow_cell_is_refused():
    stated = link.parse_run_args(["tst00001", "TST00002@4E5F6A7B"])
    assert stated == [("TST00001", None), ("TST00002", "4e5f6a7b")]
    links = link.resolve(intrinsic=[{"flow_cell_id": "TST00077", "method": "read_tag"}], runs=RUNS, stated=stated)
    assert [(ln["method"], ln["run_key"], ln["resolved"]) for ln in links] == [("manual", "TST00001_0a1b2c3d", True), ("manual", None, False)]
    with pytest.raises(InputError) as e:
        link.resolve(intrinsic=[], runs=RUNS, stated=[("TST00009", None)])
    assert "TST00009@bbbb2222" in e.value.hint and "TST00009@cccc3333" in e.value.hint
    assert link.resolve(intrinsic=[], runs=RUNS, stated=[("TST00009", "cccc3333")])[0]["run_key"] == "TST00009_cccc3333"
    with pytest.raises(UsageError):
        link.parse_run_args(["not a flow cell"])


def test_no_evidence_is_no_link():
    assert link.resolve(intrinsic=[], runs=RUNS) == [] and link.worst_confidence([]) is None
    assert link.worst_confidence([{"confidence": "high"}, {"confidence": "low"}]) == "low"


def test_tally_of_a_bam_without_index_or_tags(tmp_path, wf_v2):
    fake, _ = wf_v2
    import os
    bam = os.path.join(fake.files["sample"], "tagged.bam")
    ev = link.bam_evidence(bam)
    assert ev.fn_runs == [] and ev.records > 0 and ev.records_without_tag == ev.records       # workflow 2.x keeps no tag
    assert link.intrinsic_links(ev) == []
    assert "no index" in link.bam_evidence(bam, index=str(tmp_path / "missing.bai")).note
