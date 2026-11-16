"""Run keys, axes, the bundle envelope, the status helpers, partials and the schema."""
import gzip
import json

import numpy as np
import pytest

from flowcellar import SCHEMA_VERSION
from flowcellar.model import axes, bundle, partial, runfolder, runkey, schema, status
from flowcellar.util.errors import ValidationError


def test_run_keys():
    assert runkey.make_run_key("TST00001", "0a1b2c3d-1111-4222-a333-444455556666") == "TST00001_0a1b2c3d"
    assert runkey.split_run_key("TST00001_0a1b2c3d") == ("TST00001", "0a1b2c3d")
    assert runkey.is_run_key("TST00001_0a1b2c3d") and not runkey.is_run_key("demo_sample")
    with pytest.raises(ValueError):
        runkey.make_run_key("not a flow cell", "0a1b2c3d")
    # the status token and the barcode token sit between the flow cell and the run id
    p = runkey.parse_minknow_name("TST00001_skip_0a1b2c3d_4e5f6a7b_146.pod5")
    assert (p.run_key, p.status, p.acquisition_short, p.index) == ("TST00001_0a1b2c3d", "skip", "4e5f6a7b", 146)
    p = runkey.parse_minknow_name("TST00001_pass_barcode07_0a1b2c3d_00000000_0.bam")
    assert (p.run_key, p.status, p.barcode) == ("TST00001_0a1b2c3d", "pass", "barcode07")
    assert runkey.parse_minknow_name("demo_sample_gnome_sorted.bam") is None
    # MinKNOW's report names carry the start of the run before the run id: the date is not a run id
    for name in ("report_TST00001_20260302_0900_0a1b2c3d.json", "sample_sheet_TST00001_20260302_0900_0a1b2c3d.csv"):
        p = runkey.parse_minknow_name(name)
        assert (p.flow_cell_id, p.run_short, p.acquisition_short, p.barcode) == ("TST00001", "0a1b2c3d", None, None), name
    p = runkey.parse_minknow_name("TST00001_20260302_0a1b2c3d_1.bam")        # no time after it: taken as written
    assert p.run_short == "20260302"
    # a device part with an underscore does not matter: the name is read from its end
    assert runkey.run_key_from_folder_name("20260302_0900_P2S_00000-A_TST00001_0a1b2c3d") == "TST00001_0a1b2c3d"
    assert runkey.run_key_from_folder_name("analysis") is None


def test_run_folder_is_recognised_by_its_files():
    assert runfolder.classify_listing("anything", ["final_summary_TST00001_0a1b2c3d_4e5f6a7b.txt", "pod5"]) == "run"
    assert runfolder.classify_listing("20260302_0900_X_TST00001_0a1b2c3d", ["pod5_skip", "other_reports"]) == "run_raw_only"
    assert runfolder.classify_listing("no_key_here", ["pod5"]) is None
    assert runfolder.classify_listing("20260302_0900_X_TST00001_0a1b2c3d", ["bam_pass", "bam_fail"]) == "lookalike"
    assert runfolder.classify_listing("analysis", ["alignments", "counts"]) is None
    # reports of several runs copied into one folder: a collection, not a run folder
    several = ["final_summary_TST00001_0a1b2c3d_4e5f6a7b.txt", "final_summary_TST00002_1b2c3d4e_5f6a7b8c.txt", "throughput_TST00001_0a1b2c3d_4e5f6a7b.csv",
               "report_TST00002_20260302_0900_1b2c3d4e.json"]
    assert runfolder.classify_listing("run_reports", several) is None
    assert runfolder.classify_listing("run_reports", several[:1] + several[2:3]) == "run", "the files of one run are a run folder whatever it is called"
    assert runfolder.is_pruned("pod5_skip") and runfolder.is_pruned(".git") and not runfolder.is_pruned("analysis")
    assert runfolder.is_pruned("notebooks", extra=["notebooks"])


def test_axes_are_nested_and_binning_follows_the_edge_table():
    for name, edges in axes.EDGES.items():
        assert len(edges) - 1 == axes.AXES[name]["n"], name
        assert np.all(np.diff(edges) > 0), name
        for coarse in axes.AXES[name].get("coarsens_to", []):
            assert np.isin(axes.EDGES[coarse], edges).all(), (name, coarse)     # exact re-binning
    assert list(axes.bin_index([5, 10, 999.9, 1000, 1e7, 2e7], "read_length_log40")) == [0, 1, 80, 81, 241, 241]
    assert axes.EDGES["read_length_log40"][80] == 1000.0                       # whole decades are exact edges
    assert list(axes.bin_index([0.0, 0.999, 1.0], "fraction_0p01")) == [1, 100, 100]   # 1.00 falls in the last bin
    assert list(axes.time_bin([0, 599.9, 600, 345600], "run_time_600s")) == [0, 0, 1, 576]
    assert axes.coarsen_factor("read_length_log40", "read_length_log20") == 2
    full = np.arange(axes.n_cells("qscore_0p25"))
    assert axes.coarsen(full, 4).sum() == full.sum() and len(axes.coarsen(full, 4)) == axes.n_cells("qscore_1")
    with pytest.raises(ValueError):
        axes.coarsen_factor("qscore_1", "read_length_log20")


def _bundle():
    ident = {"run_key": "TST00001_0a1b2c3d",
             "flow_cell": {"id": "TST00001", "product_code": "FLO-MIN114", "channel_count": 512},
             "protocol_run": {"id": "0a1b2c3d-1111-4222-a333-444455556666", "short_id": "0a1b2c3d", "started": "2026-03-02T09:00:00Z"}}
    b = bundle.new_bundle(ident)
    b["tiers"]["run_health"] = {"state": "complete", "detail": "quick"}
    b["modules"]["yield"] = {"version": "0.1.0", "tier": "run_health", "state": "ok", "dims": [],
                             "cells": [{"at": [], "scalars": {"reads": 4000},
                                        "series": {"yield_series": {"axis": "run_time_600s", "start": 0, "values": {"reads": [1000, 3000]}}}}],
                             "summary": {"yield_reads": 4000}}
    b["axes"] = axes.used_axes(b["modules"])
    return b


def test_bundle_round_trip_is_deterministic(tmp_path, monkeypatch):
    monkeypatch.setenv("SOURCE_DATE_EPOCH", "1772442000")
    b = _bundle()
    p1, p2 = tmp_path / "a.fcb.json.gz", tmp_path / "b.fcb.json.gz"
    bundle.write_bundle(p1, b)
    bundle.write_bundle(p2, _bundle())
    assert p1.read_bytes() == p2.read_bytes()                       # identical inputs give identical bytes
    assert p1.read_bytes()[4:8] == b"\x00\x00\x00\x00"              # gzip header time 0
    back = bundle.read_bundle(p1)
    assert list(back)[:5] == ["format", "schema_version", "kind", "bundle_id", "content_sha256"]
    assert back["schema_version"] == SCHEMA_VERSION and back["axes"].keys() == {"run_time_600s"}
    assert json.loads(gzip.decompress(bundle.to_gzip_bytes(back))) == back


def test_measured_content_and_judgement_are_stored_apart():
    b = _bundle()
    st = {"verdict": "pass", "rules": [], "by_tier": {"run_health": "pass", "read_level": "not_available"}}
    composed = bundle.compose(b, status=st)
    assert list(composed)[-3:] == ["status", "user", "relations"]
    assert composed["relations"]["flow_cell_use"]["source"] == "unknown"
    measured, parts = bundle.split(composed)
    assert measured == bundle.order_members(b) and set(parts) == {"status", "user", "relations"}
    # the hash covers what was measured: not the verdict, not where the folder is, not when it was read
    h = bundle.content_sha256(b)
    moved = json.loads(json.dumps(composed))
    moved["provenance"]["run_dir"] = {"path": "/somewhere/else"}
    moved["provenance"]["created"] = "2030-01-01T00:00:00Z"
    assert bundle.content_sha256(moved) == h
    moved["modules"]["yield"]["cells"][0]["scalars"]["reads"] += 1
    assert bundle.content_sha256(moved) != h


def test_a_newer_schema_is_refused(tmp_path):
    b = _bundle()
    b["schema_version"] = "0.2.0"
    with pytest.raises(ValidationError) as e:
        bundle.check_version(b)
    assert e.value.exit_code == 7
    b["schema_version"] = "0.1.0"
    bundle.check_version(b)                                          # an older 0.1.x is read


def test_schema_ships_in_the_package_and_validates_a_bundle():
    s = schema.load_schema()
    assert s["$id"].endswith(":" + SCHEMA_VERSION)
    b = bundle.compose(_bundle(), status={"verdict": "pass", "rules": [],
                                          "by_tier": {"run_health": "pass", "read_level": "not_available"}})
    assert schema.schema_errors(b) == []
    b["modules"]["yield"]["state"] = "partial"                      # partial must say why (schema 0.1.1)
    assert any("partial_reason" in e for e in schema.schema_errors(b))
    b["modules"]["yield"]["partial_reason"] = "inputs_missing"
    b["modules"]["yield"]["coverage"] = {"kind": "byte_ranges", "records": 10, "records_total": 40, "fraction": 0.25}
    assert schema.schema_errors(b) == []
    b["user"]["project"] = "x"                                      # one project per run is gone
    assert schema.schema_errors(b)


def test_verdict_roll_up():
    rules = [{"id": "a", "role": "key", "tier": "run_health", "outcome": "triggered", "severity": "warn", "message": "A."},
             {"id": "b", "role": "supporting", "tier": "run_health", "outcome": "triggered", "severity": "fail", "message": "B."},
             {"id": "c", "role": "key", "tier": "read_level", "outcome": "triggered", "severity": "fail", "message": "C."},
             {"id": "d", "role": "key", "tier": "run_health", "outcome": "ok", "severity": "none", "message": "D."}]
    assert status.worst_of_rules(rules, "run_health") == "warn"       # only a key rule can produce Fail
    assert status.worst_of_rules(rules, "read_level") == "fail"
    assert status.combine_verdicts("warn", "fail") == "fail"
    assert status.combine_verdicts("pass", "not_available") == "pass"   # "not available yet" never lowers a verdict
    assert status.combine_verdicts("warn", "unknown") == "warn"
    assert status.combine_verdicts("incomplete", "pass") == "incomplete"
    assert status.combine_verdicts("incomplete", "warn") == "warn"
    assert status.combine_verdicts("not_available", "pass") == "pass"
    r = status.reasons({"rules": rules})
    assert [x["check"] for x in r] == ["a", "b", "c"] and [x["severity"] for x in r] == ["warn", "warn", "fail"]


def test_status_file_round_trip(tmp_path):
    sf = status.new_status_file("TST00001_0a1b2c3d", "ab" * 32, {"verdict": "pass", "rules": []}, bundle.default_relations())
    p = tmp_path / status.status_filename("TST00001_0a1b2c3d")
    status.write_status(p, sf)
    assert status.read_status(p) == sf
    (tmp_path / "x.json").write_text("{}")
    with pytest.raises(ValidationError):
        status.read_status(tmp_path / "x.json")


def test_partial_id_is_a_pure_function_and_partials_round_trip(tmp_path):
    fp = {"method": "quick-v1", "value": "ab" * 32}
    mods = {"dyn": "0.1.0", "len.est": "0.1.0"}
    pid = partial.partial_id(fp, "minknow.sequencing_summary", mods)
    assert pid == partial.partial_id(dict(fp), "minknow.sequencing_summary", dict(reversed(list(mods.items()))))
    assert pid != partial.partial_id(fp, "minknow.sequencing_summary", {**mods, "dyn": "0.1.1"})   # a module bump invalidates
    assert pid != partial.partial_id(fp, "minknow.pore_scan", mods)
    p = partial.new_partial(role="minknow.sequencing_summary", input_path="/x/sequencing_summary.txt", fp=fp, size=10,
                            mtime="2026-03-02T17:00:00Z", modules=mods)
    p["fragments"].append(partial.fragment("TST00001_0a1b2c3d", "dyn", dims=["end_reason"],
                                           cells=[{"at": ["signal_positive"], "scalars": {"reads": 5}}]))
    p["routing"]["runs"]["TST00001_0a1b2c3d"] = 5
    path = partial.write_partial(tmp_path, p)
    assert path == partial.partial_path(tmp_path, pid) == partial.find_partial(tmp_path, pid)
    back = partial.read_partial(path)
    assert partial.is_complete(back) and len(partial.fragments_for(back, "TST00001_0a1b2c3d", "dyn")) == 1
    assert partial.fragments_for(back, "TST00002_ffffffff") == []
    back["coverage"]["kind"] = "byte_ranges"
    assert not partial.is_complete(back)
