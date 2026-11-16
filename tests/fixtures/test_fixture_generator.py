"""The fixture generator keeps its promises: a small but complete fake run folder, invented names only."""
import csv
import json
import os

import numpy as np
import pyarrow.csv as pc

from flowcellar.model import runfolder, runkey
from flowcellar.testing.fixtures import CHANNEL_STATES, RunSpec, make_library_tree, make_run


def test_folder_is_a_run_folder_with_every_tier_1_role(fake_run):
    names = os.listdir(fake_run.run_dir)
    assert runfolder.classify_listing(os.path.basename(fake_run.run_dir), names) == "run"
    assert runkey.run_key_from_folder_name(fake_run.run_dir) == fake_run.run_key
    for role in ("minknow.final_summary", "minknow.sample_sheet", "minknow.report_json", "minknow.throughput",
                 "minknow.pore_activity", "minknow.pore_scan", "minknow.sequencing_summary"):
        assert os.path.getsize(fake_run.files[role]) > 0, role
        import fnmatch
        rel = os.path.relpath(fake_run.files[role], fake_run.run_dir)
        assert fnmatch.fnmatchcase(rel, runfolder.TIER1_ROLES[role]), role
    assert not any(n.startswith("final_summary") is False and "fastq" in n for n in names)     # basecalling is off


def test_summary_throughput_and_report_agree_exactly(fake_run):
    exp = fake_run.expected
    t = pc.read_csv(fake_run.files["minknow.sequencing_summary"], parse_options=pc.ParseOptions(delimiter="\t"))
    assert t.num_columns == 16 and t.num_rows == exp["reads"] == 4000
    assert t.column_names[3] == "filename_bam" and "sequence_length_template" not in t.column_names
    events = t["minknow_events"].to_numpy()
    assert int(events.sum()) == exp["minknow_events"]
    assert set(t["run_id"].to_pylist()) == {fake_run.protocol_run_id}          # newest layout: the protocol run id
    assert t["channel"].to_numpy().min() >= 1 and t["channel"].to_numpy().max() <= exp["channel_count"]
    assert (t["start_time"].to_numpy() + t["duration"].to_numpy()).max() <= exp["run_hours"] * 3600
    reasons = dict(zip(*np.unique(np.asarray(t["end_reason"].to_pylist()), return_counts=True)))
    assert {k: int(v) for k, v in reasons.items()} == {k: v for k, v in exp["reads_by_end_reason"].items() if v}
    with open(fake_run.files["minknow.throughput"]) as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == exp["minutes"] and [int(r["Experiment Time (minutes)"]) for r in rows] == list(range(exp["minutes"]))
    assert int(rows[-1]["Reads"]) == exp["reads"] and int(rows[-1]["Selected Events"]) == exp["minknow_events"]
    assert int(rows[6 * 60 - 1]["Reads"]) == exp["reads_at_h"][6] and exp["reads_at_h"][24] is None
    assert all(int(a["Reads"]) <= int(b["Reads"]) for a, b in zip(rows, rows[1:]))            # cumulative
    rep = json.load(open(fake_run.files["minknow.report_json"]))
    seq = [a for a in rep["acquisitions"] if a["acquisition_run_info"]["config_summary"]["purpose"] == "SEQUENCING"]
    assert len(seq) == 1 and len(rep["acquisitions"]) == 2                                    # chosen by purpose, not position
    ys = seq[0]["acquisition_run_info"]["yield_summary"]
    assert (int(ys["read_count"]), int(ys["selected_events"])) == (exp["reads"], exp["minknow_events"])
    assert rep["protocol_run_info"]["run_id"] == fake_run.protocol_run_id
    assert rep["protocol_run_info"]["software_versions"]["minknow"]["full"].startswith("6.10")


def test_pore_activity_gives_every_channel_its_full_minute(fake_run):
    exp = fake_run.expected
    t = pc.read_csv(fake_run.files["minknow.pore_activity"])
    assert set(t["Channel State"].to_pylist()) == set(CHANNEL_STATES)
    minute = t["Experiment Time (minutes)"].to_numpy()
    samples = t["State Time (samples)"].to_numpy()
    per_minute = np.bincount(minute, weights=samples)
    assert len(per_minute) == exp["minutes"] and (per_minute == exp["channel_samples_per_minute"]).all()
    assert 0.85 < exp["occupancy_first_hour"] < 0.95


def test_pore_scan_table_matches_the_report(fake_run):
    exp = fake_run.expected
    t = pc.read_csv(fake_run.files["minknow.pore_scan"])
    for col in ("channel", "well", "scan_number", "mux_scan_assessment", "pore_median", "pore_sd",
                "scan_running_voltage_mv", "seconds_since_start_of_run"):
        assert col in t.column_names, col
    assert t.num_rows == exp["pore_scan_rows"] == exp["channel_count"] * 4 * exp["pore_scans"]
    scan = t["scan_number"].to_numpy()
    good = np.isin(np.asarray(t["mux_scan_assessment"].to_pylist()), ["single_pore", "reserved_pore"])
    assert [int(good[scan == k + 1].sum()) for k in range(exp["pore_scans"])] == exp["pores_per_scan"]
    rep = json.load(open(fake_run.files["minknow.report_json"]))
    seq = next(a for a in rep["acquisitions"] if a["acquisition_run_info"]["config_summary"]["purpose"] == "SEQUENCING")
    results = seq["acquisition_run_info"]["bream_info"]["mux_scan_results"]
    assert [r["counts"]["single_pore"] + r["counts"]["reserved_pore"] for r in results] == exp["pores_per_scan"]
    assert exp["pores_per_scan"] == sorted(exp["pores_per_scan"], reverse=True)               # pores only go down
    first = next(m for m in rep["user_messages"] if m["identifier"] == "mux_scan_result")
    assert int(first["extra_data"]["total_pores"]) == exp["pores_first_scan"]
    assert int(first["extra_data"]["num_pores"]) == exp["channels_with_pore_first_scan"]
    temps = seq["temperature"][0]["temperatures"]
    assert len(temps) == exp["temperature_samples"] and temps[-1]["minion"]["heatsink_temperature"] < 34   # the stop sample reads low


def test_equal_specs_give_identical_folders_and_names_are_invented(tmp_path):
    a = make_run(tmp_path / "a", RunSpec(reads=300, hours=2, seed=5))
    b = make_run(tmp_path / "b", RunSpec(reads=300, hours=2, seed=5))
    c = make_run(tmp_path / "c", RunSpec(reads=300, hours=2, seed=6))
    assert a.run_key == b.run_key != c.run_key and a.run_key.startswith("TST00001_")
    for role in ("minknow.sequencing_summary", "minknow.throughput", "minknow.pore_activity", "minknow.pore_scan", "minknow.final_summary"):
        assert open(a.files[role], "rb").read() == open(b.files[role], "rb").read(), role
    fs = dict(line.rstrip("\n").split("=", 1) for line in open(a.files["minknow.final_summary"]))
    assert fs["instrument"] == "fakehost-01" and fs["sample_id"] == "demo_sample" and fs["basecalling_enabled"] == "0"
    assert fs["protocol_run_id"] == a.protocol_run_id and fs["acquisition_run_id"] == a.acquisition_id


def test_promethion_spec_and_the_two_root_tree(tmp_path):
    p = make_run(tmp_path / "p", RunSpec(family="promethion", reads=500, hours=1.5, flow_cell_id="TST00009", kit="SQK-NBD114-24"))
    assert p.expected["channel_count"] == 3000 and p.expected["assay_group"] == "gdna" and p.expected["pore_scans"] == 2
    runs = make_library_tree(tmp_path / "tree")
    assert [r.flow_cell_id for r in runs] == ["TST00001", "TST00001", "TST00002"]
    assert len({r.run_key for r in runs}) == 3
    assert sorted(os.listdir(tmp_path / "tree")) == ["root_a", "root_b"]
    assert runs[1].expected["pores_first_scan"] < runs[0].expected["pores_first_scan"]        # the later loading starts lower
