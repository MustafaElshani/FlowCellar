"""Quick detail from the small files, full detail from the two passes: every source against what the generator wrote."""
import copy
import json
import os
import time

import numpy as np
import pytest
from engine_support import problems

from flowcellar.aggregates import build, derive
from flowcellar.contracts import RunContext
from flowcellar.engine import api
from flowcellar.model import bundle as B
from flowcellar.model import schema
from flowcellar.testing import minknow_variants as V
from flowcellar.testing.fixtures import RunSpec, make_run


def test_quick_bundle_matches_the_generator(fake_run, quick8):
    exp, b = fake_run.expected, quick8.bundle
    ident = b["identity"]
    assert ident["run_key"] == exp["run_key"] and ident["flow_cell"]["channel_count"] == exp["channel_count"]
    assert (ident["flow_cell"]["family"], ident["device"]["type"], ident["assay"]["group"]) == (exp["family"], exp["device_type"], exp["assay_group"])
    assert ident["protocol_run"]["time_origin"] == exp["time_origin"] and ident["protocol_run"]["planned_duration_s"] == exp["planned_s"]
    assert [a["purpose"] for a in ident["acquisitions"]] == ["calibration", "sequencing"]
    y = b["modules"]["yield"]
    assert y["cells"][0]["scalars"] == {"reads": exp["reads"], "selected_events": exp["minknow_events"], "selected_raw_samples": exp["selected_raw_samples"]}
    assert y["summary"]["yield_reads"] == exp["reads"] and y["summary"]["yield_estimated_bases"] == exp["estimated_bases"]
    assert y["summary"]["yield_reads_at_6h"] == exp["reads_at_h"][6] and y["summary"]["yield_reads_at_12h"] is None
    fc = b["modules"]["fc"]
    scans = fc["cells"][0]["series"]["fc_pore_scan_series"]["values"]
    assert [a + r for a, r in zip(scans["single_pore"], scans["reserved_pore"])] == exp["pores_per_scan"]
    assert fc["summary"]["fc_pores_first_scan"] == exp["pores_first_scan"]
    assert fc["summary"]["fc_channels_with_pore_first_scan"] == exp["channels_with_pore_first_scan"]
    assert fc["summary"]["fc_occupancy_first_hour"] == pytest.approx(exp["occupancy_first_hour"], abs=5e-5)
    assert fc["cells"][0]["scalars"]["total_channel_seconds"] == exp["channel_count"] * exp["minutes"] * 60
    assert b["modules"]["run"]["summary"]["run_stop_cause"] == exp["state"]
    assert b["tiers"]["run_health"] == {"state": "complete", "detail": "quick", "missing_roles": []}
    assert fc["state"] == "partial" and fc["partial_reason"] == "waiting"
    assert [i["role"] for i in b["provenance"]["inputs"]][-2:] == ["minknow.sequencing_summary", "minknow.pore_scan"]
    assert all(i["state"] == "skipped" and i["fingerprint"] is None for i in b["provenance"]["inputs"][-2:])   # nothing large was opened
    assert problems(quick8.composed) == []


def test_quick_bundle_takes_seconds_at_most(fake_run):
    run = api.identify(fake_run.run_dir)
    t0 = time.perf_counter()
    api.quick_bundle(run)
    assert time.perf_counter() - t0 < 2.0


def test_the_temperature_sample_taken_as_the_run_stops_is_not_stored(fake_run, quick8):
    fc = quick8.bundle["modules"]["fc"]
    t = fc["cells"][0]["series"]["fc_temperature_series"]
    assert fake_run.expected["temperature_last_sample_low"]
    assert len(t["values"]["flowcell_c"]) == fake_run.expected["temperature_samples"] - 1
    assert fc["cells"][0]["scalars"]["temperature_samples"] == fake_run.expected["temperature_samples"] - 1
    assert fc["summary"]["fc_temperature_max_dev_c"] == 0.0
    assert fc["params"]["temperature_stop_sample_dropped"] is True
    assert "temperature_stop_sample_dropped" in [n["code"] for n in fc["notes"]]
    assert sorted(t["levels"]) == ["chamber_c", "flowcell_c", "target_max_c", "target_min_c"]


def test_a_low_sample_inside_the_run_is_kept_and_counted(tmp_path):
    f = V.set_temperature(make_run(tmp_path, RunSpec(flow_cell_id="TST00050", seed=50)), -1.2, sample=2)
    b = api.quick_bundle(api.identify(f.run_dir))
    assert b["modules"]["fc"]["summary"]["fc_temperature_max_dev_c"] == pytest.approx(1.1, abs=0.02)


def test_full_detail_matches_the_generator(fake_run, full8):
    exp, b = fake_run.expected, full8.bundle
    assert [r.complete for r in full8.refs] == [True, True]
    dyn = b["modules"]["dyn"]
    assert {c["at"][0]: c["scalars"]["reads"] for c in dyn["cells"]} == {k: v for k, v in exp["reads_by_end_reason"].items() if v}
    assert sum(c["scalars"]["minknow_events"] for c in dyn["cells"]) == exp["minknow_events"]
    assert dyn["params"]["completeness"] == {"reads_in_summary": exp["reads"], "reads_minknow": exp["reads"],
                                             "events_in_summary": exp["minknow_events"], "events_minknow": exp["minknow_events"], "equal": True}
    assert dyn["coverage"] == {"kind": "full", "records": exp["reads"], "records_total": exp["reads"], "fraction": 1.0}
    assert dyn["run_level"]["counts"]["dyn_pore_type"] == {"cats": ["not_set"], "n": [exp["reads"]], "role": "categories"}
    assert dyn["summary"]["dyn_signal_positive_frac"] == pytest.approx(exp["reads_by_end_reason"]["signal_positive"] / exp["reads"], abs=5e-5)
    fc = b["modules"]["fc"]
    assert fc["state"] == "ok" and fc["summary"]["fc_active_channels"] == exp["channels_with_reads"]
    assert fc["cells"][0]["scalars"]["pore_scan_rows"] == exp["pore_scan_rows"]
    grid = fc["cells"][0]["grids"]["fc_channel_time_grid"]
    assert (grid["rows"], grid["cols"], grid["cols_axis"]) == (512, 16, "run_time_1800s")
    assert set(fc["cells"][0]["hists2d"]) == {"fc_open_pore_current", "fc_open_pore_noise"}
    assert b["tiers"]["run_health"] == {"state": "complete", "detail": "full", "missing_roles": []}
    assert b["provenance"]["partials"] == sorted(r.partial_id for r in full8.refs)
    assert problems(full8.composed) == []


def test_numbers_derived_from_the_summary_equal_a_direct_computation(fake_run, full8):
    """Median speed, N50 and per-channel totals, recomputed from the file with numpy, equal the bundle to bin resolution."""
    import pyarrow.csv as pc
    t = pc.read_csv(fake_run.files["minknow.sequencing_summary"], parse_options=pc.ParseOptions(delimiter="\t"))
    ev, dur = t["minknow_events"].to_numpy().astype(float), t["duration"].to_numpy()
    er = np.asarray(t["end_reason"].to_pylist())
    est = ev * 0.6
    sel = (er == "signal_positive") & (dur >= 1.0)
    b = full8.bundle
    assert b["modules"]["dyn"]["summary"]["dyn_speed_est_median"] == pytest.approx(float(np.median(est[sel] / dur[sel])), abs=5.0)
    order = np.sort(est)[::-1]
    n50 = float(order[np.searchsorted(np.cumsum(order), order.sum() / 2)])
    assert b["modules"]["len.est"]["summary"]["len_est_n50"] == pytest.approx(n50, rel=0.06)
    assert b["modules"]["len.est"]["summary"]["len_est_mean"] == pytest.approx(float(est.mean()), abs=0.06)
    cm = b["modules"]["fc"]["cells"][0]["channels"]["fc_channel_map"]["values"]
    ch = t["channel"].to_numpy()
    assert cm["reads"] == [int(x) for x in np.bincount(ch, minlength=513)[1:]]
    assert cm["estimated_bases"] == [round(int(x) * 0.6) for x in np.bincount(ch, weights=ev, minlength=513)[1:]]
    per_reason = {c["at"][0]: c for c in b["modules"]["len.est"]["cells"]}
    h = per_reason["signal_positive"]["hists"]["len_est_hist"]
    assert derive.total(h, "w") == pytest.approx(float(est[er == "signal_positive"].sum()), abs=len(h["w"]) + 2)
    assert h["w_under"] == 0 and h["w_over"] == 0 and "w_name" in h


def test_bundle_bytes_are_deterministic(fake_run, full8):
    ctx = RunContext(deterministic=True)
    a = api.assemble(full8.run, full8.refs, None, ctx)
    b = api.assemble(full8.run, full8.refs, None, ctx)
    assert B.to_gzip_bytes(a) == B.to_gzip_bytes(b)
    assert "timings" not in a["provenance"] and a["provenance"]["created"] == a["provenance"]["updated"]
    c = api.assemble(full8.run, full8.refs)                         # times differ, the content hash does not
    assert c["content_sha256"] == a["content_sha256"] == B.content_sha256(c)
    assert api.quick_bundle(full8.run)["content_sha256"] != a["content_sha256"]


def test_context_is_recorded(fake_run):
    run = api.identify(fake_run.run_dir)
    b = api.quick_bundle(run, RunContext(library_id="lib_0123456789ab", library_name="demo", root_id="root_a", relpath="demo_experiment/demo_sample/x"))
    prov = b["provenance"]
    assert prov["library"] == {"id": "lib_0123456789ab", "name": "demo"} and prov["run_dir"]["root"] == "root_a"
    assert prov["inputs"][0]["root"] == "root_a" and prov["inputs"][0]["relpath"].startswith("demo_experiment/demo_sample/x/")
    assert schema.schema_errors(b) == []


@pytest.mark.parametrize("generation", ["5.3", "5.9"])
def test_older_layouts_in_full_detail(tmp_path, build, generation):
    f = V.make_old_layout(tmp_path, generation=generation, flow_cell_id="TST00051", seed=51)
    built = build(f)
    assert f.expected["summary_columns"] == 15
    assert built.bundle["modules"]["dyn"]["params"]["completeness"]["equal"] is True
    assert built.bundle["tiers"]["run_health"]["detail"] == "full"
    assert built.bundle["modules"]["fc"]["summary"]["fc_pores_first_scan"] == f.expected["pores_first_scan"]
    assert problems(built.composed) == []


def test_report_without_scan_table_takes_pores_from_the_messages(tmp_path, build):
    f = V.make_without_scan_table(tmp_path, flow_cell_id="TST00052", seed=52)
    built = build(f)
    fc = built.bundle["modules"]["fc"]
    assert fc["summary"]["fc_pores_first_scan"] == f.expected["pores_first_scan"]
    codes = [n["code"] for n in fc["notes"]]
    assert "pore_scans_from_messages" in codes and "pore_scan_sources_differ" not in codes
    assert problems(built.composed) == []


def test_promethion_72h(tmp_path, build):
    f = V.make_promethion_72h(tmp_path)
    built = build(f)
    b = built.bundle
    assert b["identity"]["flow_cell"]["family"] == "promethion" and b["identity"]["flow_cell"]["channel_layout"] == "promethion_3000"
    assert b["modules"]["yield"]["summary"]["yield_reads_at_72h"] == f.expected["reads_at_h"][72]
    grid = b["modules"]["fc"]["cells"][0]["grids"]["fc_channel_time_grid"]
    assert (grid["rows"], grid["cols_axis"]) == (3000, "run_time_7200s") and grid["cols"] == 37
    assert built.rules["yield_reads_at_72h"]["role"] == "key" and built.rules["yield_estimated_bases_at_72h"]["role"] == "supporting"
    assert problems(built.composed) == []
    assert len(B.to_gzip_bytes(b)) < 260_000


# ------------------------------------------------------------------------------------------ damaged input
def _cut(path, fraction=0.5):
    data = open(path, "rb").read()
    with open(path, "wb") as fh:
        fh.write(data[: int(len(data) * fraction)])


@pytest.mark.parametrize("role", ["minknow.report_json", "minknow.throughput", "minknow.pore_activity", "minknow.final_summary",
                                  "minknow.sample_sheet"])
def test_a_cut_small_file_costs_a_section_never_the_bundle(tmp_path, role):
    f = make_run(tmp_path, RunSpec(flow_cell_id="TST00053", seed=53))
    _cut(f.files[role])
    run = api.identify(f.run_dir)
    b = api.quick_bundle(run)
    assert schema.schema_errors(b) == []
    states = {k: v["state"] for k, v in b["modules"].items()}
    if role == "minknow.report_json":
        assert states["run"] == "error" and "report_unreadable" in [w["code"] for w in b["provenance"]["warnings"]]
        assert states["yield"] == "ok"                               # the other sections still come from the CSV files
        assert b["modules"]["yield"]["params"]["events_to_base_ratio_source"] == "throughput"
    if role == "minknow.throughput":
        assert states["yield"] == "partial" and b["modules"]["yield"]["partial_reason"] == "partly_read"
    if role == "minknow.pore_activity":
        assert "pore_activity_incomplete" in [n["code"] for n in b["modules"]["fc"]["notes"]]
    sf = api.judge(b, [], in_library=False)
    assert schema.schema_errors(B.compose(b, status=sf["status"], relations=sf["relations"])) == []
    assert sf["status"]["verdict"] in ("pass", "warn", "fail", "incomplete")


@pytest.mark.parametrize("role", ["minknow.report_json", "minknow.throughput", "minknow.pore_activity"])
def test_garbage_in_a_small_file_never_crashes(tmp_path, role):
    f = make_run(tmp_path, RunSpec(flow_cell_id="TST00054", seed=54))
    with open(f.files[role], "wb") as fh:
        fh.write(b"\x00\xff not what MinKNOW wrote \n\n,,,\n" * 50)
    b = api.quick_bundle(api.identify(f.run_dir))
    assert schema.schema_errors(b) == []
    sf = api.judge(b, [], in_library=False)
    assert sf["status"]["verdict"] in ("pass", "warn", "fail", "incomplete")
    rules = {r["id"]: r for r in sf["status"]["rules"]}
    if role == "minknow.report_json":
        assert rules["run_messages"]["outcome"] == "not_evaluated" and rules["run_stop_cause"]["outcome"] == "not_evaluated"


def test_missing_end_of_run_files_give_an_incomplete_verdict_with_its_cause(tmp_path):
    f = make_run(tmp_path, RunSpec(flow_cell_id="TST00055", seed=55))
    for role in ("minknow.report_json", "minknow.throughput", "minknow.pore_activity"):
        os.unlink(f.files[role])
    b = api.quick_bundle(api.identify(f.run_dir))
    assert schema.schema_errors(b) == []
    st = api.judge(b, [], in_library=False)["status"]
    assert (st["verdict"], st["cause"]) == ("incomplete", "inputs_missing")
    assert set(b["tiers"]["run_health"]["missing_roles"]) >= {"minknow.report_json", "minknow.throughput", "minknow.pore_activity"}


def test_two_files_for_one_role(tmp_path):
    f = make_run(tmp_path, RunSpec(flow_cell_id="TST00056", seed=56))
    extra = f.files["minknow.throughput"].replace("throughput_", "throughput_copy_")
    with open(extra, "w") as fh:
        fh.write("Experiment Time (minutes),Reads\n")
    b = api.quick_bundle(api.identify(f.run_dir))
    assert b["modules"]["yield"]["summary"]["yield_reads"] == f.expected["reads"]             # the larger file is read
    assert "several_files_for_role" in [w["code"] for w in b["provenance"]["warnings"]]


def test_sparkline_and_card_columns(full8):
    card = api.card(full8.bundle, full8.status_file)
    c = card["columns"]
    assert c["id"] == full8.run.run_key and c["detail_t1"] == "full" and c["run_state"] == "finished"
    assert len(c["spark_output"]) == 32 and len(c["spark_channels"]) == 32
    assert build.u6_decode(c["spark_output"], build.B64URL).max() == 63
    assert c["reads"] == full8.fixture.expected["reads"] and c["active_channels"] == full8.fixture.expected["channels_with_reads"]
    assert c["verdict"] == full8.status["verdict"] and c["headline"] == full8.status["headline"]
    assert len(c["reasons"]) == len(c["reason_levels"]) == len(c["reason_checks"]) == len(c["reason_tiers"]) <= 3
    assert c["limits"].startswith("cdna@")
    json.dumps(card)                                                 # a card is plain JSON
    assert copy.deepcopy(card["facts"])["run_key"] == c["id"]
