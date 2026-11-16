"""Recognising a run folder by its files, identity with its precedence, assay detection, run states."""
import os
import shutil
import time

import pytest

from flowcellar.engine import api
from flowcellar.engine.identity import chemistry_of, detect_assay
from flowcellar.model import schema
from flowcellar.testing import minknow_variants as V
from flowcellar.testing.fixtures import RunSpec, make_run
from flowcellar.util.errors import ExitCode, InputError


def test_identify_reads_small_files_only_and_is_fast(fake_run):
    t0 = time.perf_counter()
    run = api.identify(fake_run.run_dir)
    assert time.perf_counter() - t0 < 1.0
    assert (run.run_key, run.flow_cell_id, run.protocol_run_id) == (fake_run.run_key, fake_run.flow_cell_id, fake_run.protocol_run_id)
    assert run.kind == "run" and run.run_state == "finished"
    assert run.started == fake_run.expected["started"]
    assert set(run.files) >= {"minknow.report_json", "minknow.final_summary", "minknow.sequencing_summary", "minknow.pore_scan"}
    assert run.sizes["minknow.sequencing_summary"] == os.path.getsize(fake_run.files["minknow.sequencing_summary"])
    assert run.facts["device_type"] == fake_run.expected["device_type"] and run.facts["channel_count"] == 512
    assert run.facts["kit"] == fake_run.expected["kit"] and run.facts["product_code"] == fake_run.expected["product_code"]
    assert run.newest_mtime > 0


def test_a_run_is_known_by_its_files_not_by_the_folder_name(tmp_path, fake_run):
    moved = tmp_path / "renamed_by_someone"
    shutil.copytree(fake_run.run_dir, moved)
    run = api.identify(str(moved))
    assert run.run_key == fake_run.run_key and run.kind == "run"


def test_identity_precedence_report_then_final_summary_then_sample_sheet_then_names(tmp_path):
    f = make_run(tmp_path, RunSpec(flow_cell_id="TST00040", seed=40))
    want = (f.run_key, f.protocol_run_id)
    for role in ("minknow.report_json", "minknow.final_summary", "minknow.sample_sheet"):
        run = api.identify(f.run_dir)
        assert (run.run_key, run.protocol_run_id) == want, role
        os.unlink(f.files[role])
    run = api.identify(f.run_dir)                 # only file names are left: the short id, no full id
    assert run.run_key == f.run_key and run.protocol_run_id is None and run.run_state != "finished"


def test_report_head_is_enough_for_identity_even_when_the_report_is_cut(tmp_path):
    f = make_run(tmp_path, RunSpec(flow_cell_id="TST00041", seed=41))
    data = open(f.files["minknow.report_json"], "rb").read()
    with open(f.files["minknow.report_json"], "wb") as fh:
        fh.write(data[: len(data) // 2])
    run = api.identify(f.run_dir)
    assert run.run_key == f.run_key and run.facts["channel_count"] == 512


@pytest.mark.parametrize("generation", ["5.3", "5.5", "5.9", "6.2"])
def test_older_layouts_give_the_same_identity_fields(tmp_path, generation):
    f = V.make_old_layout(tmp_path, generation=generation, flow_cell_id="TST00042", seed=42)
    run = api.identify(f.run_dir)
    assert run.run_key == f.run_key and run.protocol_run_id == f.protocol_run_id
    b = api.quick_bundle(run)
    ident = b["identity"]
    assert ident["software"]["minknow"] == V.GENERATIONS[generation]          # top-level member before MinKNOW 5.7
    assert ident["device"]["position"] == f.expected["device_position"]         # instrument = device, position empty
    assert ident["protocol_run"]["planned_duration_s"] == f.expected["planned_s"]   # from --experiment_time on 5.3
    assert ident["kit"]["sequencing_kit"] == f.expected["kit"]                  # no kit_info: from the tags or the protocol id
    assert [a["purpose"] for a in ident["acquisitions"]] == ["calibration", "sequencing"]
    assert ident["protocol_run"]["time_origin"] == f.expected["time_origin"]
    assert schema.schema_errors(b) == []


def test_raw_only_folder_is_a_stub(tmp_path):
    f = V.make_raw_only(tmp_path, flow_cell_id="TST00043", seed=43)
    run = api.identify(f.run_dir)
    assert (run.kind, run.run_state, run.run_key, run.protocol_run_id) == ("run_raw_only", "raw_only", f.run_key, None)
    b = api.quick_bundle(run)
    assert schema.schema_errors(b) == []
    assert b["modules"] == {} and b["tiers"]["run_health"]["state"] == "absent"
    assert b["identity"]["protocol_run"]["run_state"] == "raw_only" and b["identity"]["protocol_run"]["id"] is None
    assert [w["code"] for w in b["provenance"]["warnings"]] == ["raw_only"]
    assert api.full_pass(run, str(tmp_path / "partials")) == []
    st = api.judge(b, [])["status"]
    assert (st["verdict"], st["cause"], st["by_tier"]) == ("incomplete", "raw_only", {"run_health": "incomplete", "read_level": "not_available"})


def test_lookalike_is_not_a_run(tmp_path):
    look = V.make_lookalike(tmp_path)
    with pytest.raises(InputError) as e:
        api.identify(look)
    assert e.value.exit_code == ExitCode.INPUT and "basecaller output" in e.value.message


def test_a_folder_that_is_nothing_is_an_input_error(tmp_path):
    (tmp_path / "notes.txt").write_text("hello")
    with pytest.raises(InputError):
        api.identify(str(tmp_path))
    with pytest.raises(InputError):
        api.identify(str(tmp_path / "missing"))


def test_running_and_interrupted(tmp_path):
    f = V.make_running(tmp_path / "a", flow_cell_id="TST00044", seed=44)
    run = api.identify(f.run_dir)
    assert run.run_state == "running" and run.run_key == f.run_key
    b = api.quick_bundle(run)
    assert schema.schema_errors(b) == [] and b["tiers"]["run_health"]["state"] == "live"
    assert api.full_pass(run, str(tmp_path / "p")) == []                   # files still change: nothing is cached
    st = api.judge(b, [])["status"]
    assert st["verdict"] == "running" and st["by_tier"]["run_health"] == "running" and st["judged"]["frozen"] is False
    # the same folder, untouched for two days, is an interrupted run
    old = time.time() - 48 * 3600
    for base, dirs, files in os.walk(f.run_dir):
        for n in dirs + files:
            os.utime(os.path.join(base, n), (old, old))
    run = api.identify(f.run_dir)
    assert run.run_state == "interrupted"
    st = api.judge(api.quick_bundle(run), [])["status"]
    assert (st["verdict"], st["cause"]) == ("incomplete", "no_final_summary")


def test_stale_after_hours_comes_from_the_context(tmp_path):
    from flowcellar.contracts import RunContext
    f = V.make_interrupted(tmp_path, with_reports=False, age_hours=3.0, flow_cell_id="TST00045", seed=45)
    run = api.identify(f.run_dir)
    assert run.run_state == "running"                                       # default: changed within 6 h
    b = api.quick_bundle(run, RunContext(stale_after_hours=2.0))
    assert b["identity"]["protocol_run"]["run_state"] == "interrupted"


def test_interrupted_run_with_reports_is_judged_on_what_it_has(tmp_path, partials_dir):
    f = V.make_interrupted(tmp_path, flow_cell_id="TST00046", seed=46, hours=3.0, planned_hours=72.0)
    run = api.identify(f.run_dir)
    b = api.assemble(run, api.full_pass(run, str(partials_dir)))
    assert b["tiers"]["run_health"]["state"] == "partial" and "minknow.final_summary" in b["tiers"]["run_health"]["missing_roles"]
    st = api.judge(b, [])["status"]
    assert st["verdict"] == "incomplete" and st["cause"] == "no_final_summary" and st["rules"]


@pytest.mark.parametrize("kit,product,category,want", [
    ("SQK-PCB114-24", "FLO-PRO114M", ["RNA", "PCR"], ("cdna_pcr", "cdna", "high")),
    ("SQK-PCS111", "FLO-PRO002", [], ("cdna_pcr", "cdna", "medium")),
    ("SQK-DCS109", "FLO-MIN106", ["RNA"], ("cdna_direct", "cdna", "high")),
    ("SQK-RNA004", "FLO-PRO004RA", ["RNA"], ("direct_rna", "drna", "high")),
    ("SQK-NBD114-24", "FLO-MIN114", ["DNA", "PCR-Free"], ("gdna", "gdna", "high")),
    ("SQK-LSK114", "FLO-PRO114M", [], ("gdna", "gdna", "medium")),
    ("SQK-RPB114-24", "FLO-MIN114", ["DNA"], ("amplicon", "gdna", "high")),
    (None, "FLO-MIN114", ["RNA"], ("cdna_pcr", "cdna", "medium")),
    (None, "FLO-MIN114", [], ("unknown", "unknown", "low")),
])
def test_assay_detection(kit, product, category, want):
    a = detect_assay(kit, product, f"sequencing/sequencing_X_DNA:{product}:{kit}:400", category, None)
    assert (a["type"], a["group"], a["confidence"]) == want
    assert a["evidence"]


def test_a_stated_assay_always_wins_and_chemistry_follows_the_product():
    a = detect_assay("SQK-LSK114", "FLO-PRO114M", None, ["DNA"], None, override="cdna_pcr")
    assert (a["type"], a["group"], a["source"], a["confidence"]) == ("cdna_pcr", "cdna", "user", "certain")
    assert chemistry_of("FLO-PRO114M") == "r10.4.1" and chemistry_of("FLO-MIN106") == "r9.4.1" and chemistry_of("FLO-PRO004RA") == "rna004"
    assert detect_assay("SQK-PCB114-24", "FLO-MIN114", None, [], None)["barcoded"] is True


@pytest.mark.parametrize("kind", [k for k in V.KINDS if k not in ("series", "promethion-72h", "reused-flow-cell")])
def test_every_kind_of_fake_folder_is_recognised_for_what_it_is(tmp_path, kind):
    (path,) = V.make_kind(tmp_path, kind)
    if kind == "lookalike":
        with pytest.raises(InputError):
            api.identify(path)
        return
    run = api.identify(path)
    want = {"running": "running", "interrupted": "interrupted", "interrupted-no-reports": "interrupted", "raw-only": "raw_only"}.get(kind, "finished")
    assert run.run_state == want and run.kind == ("run_raw_only" if kind == "raw-only" else "run")
    b = api.quick_bundle(run)
    assert schema.schema_errors(b) == []
    assert api.judge(b, [], in_library=False)["status"]["verdict"] in ("pass", "warn", "fail", "incomplete", "running")


def test_fake_folders_hold_invented_names_only(tmp_path):
    import re
    real = re.compile(r"(?<![A-Za-z0-9])[PF][A-Z]{2}[0-9]{5}(?![0-9])")
    for kind in ("old-5.3", "raw-only", "running"):
        for base, _dirs, files in os.walk(V.make_kind(tmp_path / kind, kind)[0]):
            for n in files:
                assert not real.search(n) and "TST" in n or n.endswith(".pod5")
                if n.endswith((".txt", ".csv", ".json", ".tsv")) and os.path.getsize(os.path.join(base, n)) < 2_000_000:
                    assert not real.search(open(os.path.join(base, n), errors="replace").read())
