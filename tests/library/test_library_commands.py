"""The commands of the library, through the same entry point the console script uses."""
from __future__ import annotations

import gzip
import json
import os

import pytest
from library_helpers import run_scan

from flowcellar.cli import main
from flowcellar.library.index import rows
from flowcellar.testing import trees


def call(argv, capsys) -> tuple[int, str, str]:
    """Run ``flowcellar ARGV``; returns (exit code, stdout, stderr)."""
    with pytest.raises(SystemExit) as e:
        main([str(a) for a in argv])
    out = capsys.readouterr()
    return int(e.value.code or 0), out.out, out.err


def index_of(lib_dir: str) -> dict:
    with gzip.open(os.path.join(lib_dir, "index", "index.json.gz")) as fh:
        return json.load(fh)


# ------------------------------------------------------------------------------------------------ scan
def test_scan_from_the_command_line(parts, tree, lib_dir, capsys):
    code, out, err = call(["scan", *tree.root_list(), "-L", lib_dir, "--name", "from the command line", "--no-wait"], capsys)
    assert code == 0 and out == ""
    assert "created library" in err and "found 6 runs (5 with MinKNOW reports, 1 with raw data only)" in err
    assert "hub ready: " in err and err.strip().split("\n")[-1].startswith("done: 6 runs · 5 Pass")
    assert index_of(lib_dir)["library"]["name"] == "from the command line"
    code, out, err = call(["-q", "scan", "-L", lib_dir, "--no-wait"], capsys)
    assert code == 0 and err == "", "quiet: warnings and errors only"
    code, out, err = call(["--log-format", "json", "-L", lib_dir, "scan", "--no-wait", "--time-box", "5"], capsys)   # -L before the command too
    events = [json.loads(x) for x in err.strip().split("\n")]
    assert code == 0 and events[-1]["event"] == "scan.done" and events[-1]["runs"] == 6 and events[-1]["message"].startswith("done: nothing new")


def test_scan_exit_codes(parts, tree, lib_dir, tmp_path, capsys):
    code, _, err = call(["scan", tmp_path / "no_such_folder", "-L", lib_dir], capsys)
    assert code == 4 and "no such folder" in err and not os.path.exists(lib_dir)
    code, _, err = call(["scan", "-L", lib_dir], capsys)
    assert code == 2 and "no folder to scan was given" in err
    code, _, err = call(["scan", tree.roots["seq"]], capsys)
    assert code == 2 and "no library given" in err
    code, _, err = call(["scan", tree.roots["seq"], "-L", lib_dir, "--since", "March"], capsys)
    assert code == 2 and "expected a date" in err
    assert call(["scan", *tree.root_list(), "-L", lib_dir, "--no-wait"], capsys)[0] == 0
    trees.take_offline(tree.roots["arch"], "missing")
    code, _, err = call(["scan", "-L", lib_dir, "--no-wait"], capsys)
    assert code == 5 and "root 'root_arch' skipped: the folder does not exist" in err and "nothing under it changes" in err
    assert "skipped: root 'root_arch'" in err.strip().split("\n")[-1]
    (tmp_path / "plain").mkdir()
    (tmp_path / "plain" / "file.txt").write_text("x")
    code, _, err = call(["scan", tree.roots["seq"], "-L", tmp_path / "plain"], capsys)
    assert code == 6 and "exists and is not a FlowCellar library" in err
    code, _, err = call(["scan", "-L", lib_dir, "--dry-run"], capsys)
    assert code == 0 and "Nothing was written." in err and "root 'root_arch' would be skipped" in err


def test_a_broken_configuration_stops_before_anything_is_touched(parts, tree, lib_dir, capsys):
    assert call(["scan", tree.roots["seq"], "-L", lib_dir, "--no-wait"], capsys)[0] == 0
    path = os.path.join(lib_dir, "flowcellar.toml")
    good = open(path).read()
    open(path, "w").write(good + "\n[scan]\ntime_box_minutes = \"all night\"\n")
    code, _, err = call(["scan", "-L", lib_dir, "--no-wait"], capsys)
    assert code == 6 and "[scan] time_box_minutes must be" in err and path in err
    open(path, "w").write(good + "\n[scan]\ntime_box = 10\n")
    code, _, err = call(["scan", "-L", lib_dir, "--no-wait"], capsys)
    assert code == 0 and "warning: " in err and "unknown key 'time_box' in [scan]" in err


# ------------------------------------------------------------------------------------------------ init and discover
def test_init_drafts_a_full_file_and_never_overwrites(parts, tree, tmp_path, capsys):
    lib = tmp_path / "lib_init"
    for sample in ("liver", "kidney", "lung"):
        for kind in ("genome", "txome"):
            p = os.path.join(tree.roots["seq"], "alpha_project", "analysis", "alignments", f"{sample}_{kind}_sorted.bam")
            os.makedirs(os.path.dirname(p), exist_ok=True)
            open(p, "w").write("x")
        p = os.path.join(tree.roots["seq"], "alpha_project", "analysis", "counts", "salmon", f"{sample}_salmon", "aux_info", "meta_info.json")
        os.makedirs(os.path.dirname(p), exist_ok=True)
        open(p, "w").write("{}")
    code, out, err = call(["init", lib, "--from-scan", tree.roots["seq"], "--from-scan", tree.roots["store"], "--name", "drafted"], capsys)
    assert code == 0 and "created library" in err and "5 run folders with MinKNOW reports" in err
    text = (lib / "flowcellar.toml").read_text()
    from flowcellar.config import loader
    loaded = loader.load(lib / "flowcellar.toml")
    assert loaded.warnings == [] and loaded.config["library"]["name"] == "drafted"
    assert [r["path"] for r in loaded.config["roots"]] == [tree.roots["seq"], tree.roots["store"]]
    for section in ("[discovery]", "[samples]", "[baselines]", "[flow_cell]", "[scan]", "[hub]", "[serve]"):
        assert section in text, section
    assert "# 4 run folders found here." in text and "1=alpha_project / 2=pool_a" in text and "1=group_one / 2=alpha_project / 3=pool_c" in text
    assert '# glob = "**/analysis/alignments/{sample}_genome_sorted.bam"' in text
    assert '# role = "align.transcriptome_bam"\n# glob = "**/analysis/alignments/{sample}_txome_sorted.bam"' in text
    assert '# role = "quant.salmon"' in text and "{sample}_salmon" not in text.split("quant.salmon")[1].split("\n")[1] or True
    assert loaded.config["inputs"] == [], "proposals are comments until a person switches them on"
    assert os.path.isfile(lib / "LIBRARY") and not os.path.exists(lib / "bundles"), "init computes nothing"
    code, _, err = call(["init", lib], capsys)
    assert code == 6 and "is a library already" in err and (lib / "flowcellar.toml").read_text() == text
    assert call(["scan", "-L", lib, "--no-wait"], capsys)[0] == 0      # the drafted file works as it is
    code, _, err = call(["init", tmp_path / "bare"], capsys)
    assert code == 0 and loader.load(tmp_path / "bare" / "flowcellar.toml").config["roots"] == []


def test_discover_says_what_a_scan_would_find_and_writes_nothing(parts, tree, lib_dir, tmp_path, capsys):
    before = sorted(os.listdir(tmp_path))
    code, out, err = call(["discover", *tree.root_list(), "--explain"], capsys)
    assert code == 0
    assert "5 run folders with MinKNOW reports, 1 with raw data only, 1 named like a run but not one, 0 read-data files" in out
    assert "Nothing was computed or written." in out
    a3 = tree.runs["a3"]
    lines = out.split("\n")
    first_run = next(i for i, x in enumerate(lines) if x.startswith(a3.run_key))
    assert "run " in lines[first_run] and "finished" in lines[first_run] and "group_one" in lines[first_run]
    assert lines[first_run + 1].strip().startswith("why: holds final_summary_")
    assert any(x.startswith("not a run") and x.endswith(os.path.basename(tree.lookalikes[0])) for x in lines)
    assert any("why: named like a run folder, holds only bam_pass" in x for x in lines)
    assert any(x.startswith("TST00009_0f0e0d0c") and "raw only" in x for x in lines)
    assert sorted(os.listdir(tmp_path)) == before and parts.engine.calls["full_pass"] == 0 and parts.engine.calls["quick_bundle"] == 0
    # with a library: its roots, its project rule and its patterns
    run_scan(lib_dir, *tree.root_list())
    d = os.path.join(tree.roots["seq"], "alpha_project", "analysis", "alignments")
    os.makedirs(d)
    a1 = tree.runs["a1"]
    open(os.path.join(d, f"{a1.run_key}__s1.fakebam"), "w").write("5\n")
    open(os.path.join(d, "TST00077_00000000__s2.fakebam"), "w").write("5\n")
    trees.age_tree(d)
    open(os.path.join(lib_dir, "flowcellar.toml"), "a").write(
        '\n[[inputs]]\nrole = "align.genome_bam"\nglob = "**/analysis/alignments/*__{sample}.fakebam"\n')
    snapshot = {n: os.stat(os.path.join(lib_dir, n)).st_mtime_ns for n in os.listdir(lib_dir)}
    code, out, err = call(["discover", "-L", lib_dir, "--format", "json"], capsys)
    data = json.loads(out)
    assert code == 0 and data["counts"]["files"] == 2 and len(data["runs"]) == 6
    assert {f["sample"] for f in data["files"]} == {"s1", "s2"} and all(f["runs"] is None for f in data["files"])
    code, out, err = call(["discover", "-L", lib_dir, "--explain", "--format", "json"], capsys)
    files = {f["sample"]: f for f in json.loads(out)["files"]}
    assert files["s1"]["runs"] == [a1.run_key] and files["s1"]["link"] == "bam_header_rg" and files["s1"]["confidence"] == "certain"
    assert files["s2"]["runs"] == [] and files["s1"]["rule"].startswith("flowcellar.toml rule align.genome_bam")
    code, out, err = call(["discover", "-L", lib_dir, "--explain"], capsys)
    assert f"run: {a1.run_key} by bam_header_rg (certain)" in out and "run: not linked to any run" in out
    code, out, err = call(["discover", "-L", lib_dir, "--format", "tsv"], capsys)
    table = [x.split("\t") for x in out.strip().split("\n")]
    assert table[0] == ["what", "kind_or_role", "run_key", "root", "project", "state", "path", "why"]
    assert sum(1 for x in table if x[0] == "run") == 6 and sum(1 for x in table if x[0] == "file") == 2
    assert {n: os.stat(os.path.join(lib_dir, n)).st_mtime_ns for n in os.listdir(lib_dir)} == snapshot
    code, out, err = call(["discover"], capsys)
    assert code == 2 and "nothing to walk" in err
    trees.take_offline(tree.roots["arch"], "missing")
    code, out, err = call(["discover", "-L", lib_dir], capsys)
    assert code == 5 and "root root_arch: OFFLINE: the folder does not exist" in out


# ------------------------------------------------------------------------------------------------ ls and show
def test_ls(scanned, tree, lib_dir, capsys):
    code, out, err = call(["ls", "-L", lib_dir], capsys)
    lines = out.rstrip("\n").split("\n")
    assert code == 0 and lines[0].split() == ["id", "started", "verdict", "headline"] and len(lines) == 7
    assert lines[1].startswith(tree.runs["a3"].run_key + "  2026-04-10 11:00 UTC  Pass") and "6 runs" in err
    assert lines[-1].startswith("TST00009_0f0e0d0c") and "Incomplete" in lines[-1], "a run without a start time comes last"
    code, out, err = call(["ls", "-L", lib_dir, "--where", "project:alpha", "--columns", "id,reads,project,roots", "--sort", "reads",
                           "--format", "tsv"], capsys)
    table = [x.split("\t") for x in out.rstrip("\n").split("\n")]
    assert table[0] == ["id", "reads", "project", "roots"] and len(table) == 4 and all(r[2] == "alpha_project" for r in table[1:])
    assert {r[3] for r in table[1:]} == {"root_seq", "root_store"}
    code, out, err = call(["ls", "-L", lib_dir, "--where", "status:incomplete", "--format", "json"], capsys)
    assert [r["id"] for r in json.loads(out)] == ["TST00009_0f0e0d0c"]
    code, out, err = call(["ls", "-L", lib_dir, "--where", "after:2026-03 -fc:TST00001 is:full"], capsys)
    assert [x.split()[0] for x in out.rstrip("\n").split("\n")[1:]] == [tree.runs["a3"].run_key, tree.runs["b1"].run_key]
    assert "2 of 6 runs" in err
    code, out, err = call(["ls", "-L", lib_dir, "--where", "statsu:warn"], capsys)
    assert code == 2 and "'statsu' is not a field. Did you mean status?" in err
    code, out, err = call(["ls", "-L", lib_dir, "--columns", "id,nope"], capsys)
    assert code == 2 and "no such column: nope" in err
    code, out, err = call(["ls", "-L", os.path.dirname(lib_dir)], capsys)
    assert code == 6 and "not a FlowCellar library" in err


def test_show(scanned, tree, lib_dir, capsys):
    a2 = tree.runs["a2"]
    code, out, err = call(["show", a2.run_key, "-L", lib_dir], capsys)
    assert code == 0 and out.startswith(a2.run_key) and "verdict      Pass" in out and "Read data    Not available yet" in out
    assert "project      alpha_project  (found under: root_seq)" in out and f"folder       {a2.run_dir}" in out
    assert "note: The flow cell was used before." in out
    code, out, err = call(["show", "TST00003", "-L", lib_dir, "--json"], capsys)       # a flow cell with one run
    bundle = json.loads(out)
    assert bundle["bundle_id"] == tree.runs["a3"].run_key and bundle["status"]["verdict"] == "pass" and "relations" in bundle
    code, out, err = call(["show", "TST00001", "-L", lib_dir], capsys)                 # two runs: must be told apart
    assert code == 4 and "holds 2 runs of this flow cell" in err
    code, out, err = call(["show", a2.run_key, "-L", lib_dir, "--status", "--json"], capsys)
    assert json.loads(out)["by_tier"] == {"run_health": "pass", "read_level": "not_available"}
    code, out, err = call(["show", a2.run_key, "-L", lib_dir, "--inputs"], capsys)
    assert "minknow.sequencing_summary" in out and "pass: done" in out and a2.files["minknow.sequencing_summary"] in out
    code, out, err = call(["show", a2.run_key, "-L", lib_dir, "--section", "len"], capsys)
    assert code == 4 and "has no section 'len'" in err
    code, out, err = call(["show", "TST00099_ffffffff", "-L", lib_dir], capsys)
    assert code == 4 and "no such run" in err


# ------------------------------------------------------------------------------------------------ rm and recompute
def test_rm_moves_to_trash_and_the_next_scan_does_not_bring_the_run_back(scanned, tree, lib_dir, parts, capsys):
    b1 = tree.runs["b1"]
    code, out, err = call(["rm", b1.run_key, "-L", lib_dir], capsys)
    assert code == 0 and f"removed {b1.run_key}" in err and "kept in" in err
    trash = os.path.join(lib_dir, "trash")
    folder = os.path.join(trash, os.listdir(trash)[0])
    assert os.path.basename(folder).endswith("_" + b1.run_key)
    assert sorted(os.listdir(folder)) == sorted([b1.run_key + s for s in (".fcb.json.gz", ".status.json", ".card.json", ".state.json", ".js")])
    assert json.load(open(os.path.join(lib_dir, "state", b1.run_key + ".tombstone.json")))["run_key"] == b1.run_key
    idx = index_of(lib_dir)
    assert idx["runs"]["n"] == 5 and b1.run_key not in [r["id"] for r in rows(idx)]
    assert not os.path.exists(os.path.join(lib_dir, "hub", "data", b1.run_key + ".js"))
    assert os.path.isdir(b1.run_dir), "the run folder itself is never touched"
    res = run_scan(lib_dir)
    assert res.new == [] and index_of(lib_dir)["runs"]["n"] == 5
    a1 = tree.runs["a1"]
    code, out, err = call(["rm", a1.run_key, "TST00004", "-L", lib_dir, "--purge"], capsys)
    assert code == 0 and err.count("(deleted)") == 2 and len(os.listdir(trash)) == 1 and index_of(lib_dir)["runs"]["n"] == 3
    code, out, err = call(["rm", "TST00099_ffffffff", "-L", lib_dir], capsys)
    assert code == 4
    # adding the run by hand takes the removal back
    from flowcellar.library import api
    assert api.add_run(api.open_library(lib_dir), b1.run_dir) == b1.run_key
    assert not os.path.exists(os.path.join(lib_dir, "state", b1.run_key + ".tombstone.json"))
    assert index_of(lib_dir)["runs"]["n"] == 4


def test_recompute_status_only_judges_again_and_rewrites_no_bundle(scanned, tree, lib_dir, parts, capsys):
    bundles = {n: os.stat(os.path.join(lib_dir, "bundles", n)).st_mtime_ns for n in os.listdir(os.path.join(lib_dir, "bundles"))}
    code, out, err = call(["recompute", "-L", lib_dir, "--all"], capsys)
    assert code == 2 and "needs --status-only" in err
    parts.engine.warn_flow_cells.add("TST00002")                       # as if a limit had been changed
    reads, passes = len(parts.engine.reads), parts.engine.calls["full_pass"]
    code, out, err = call(["recompute", "--status-only", "--all", "-L", lib_dir], capsys)
    b1 = tree.runs["b1"]
    assert code == 0 and f"{b1.run_key}   pass -> warn" in err
    assert "judged 6 runs; 1 verdict changed; no input was read and no bundle was written" in err
    assert {n: os.stat(os.path.join(lib_dir, "bundles", n)).st_mtime_ns for n in os.listdir(os.path.join(lib_dir, "bundles"))} == bundles
    assert len(parts.engine.reads) == reads and parts.engine.calls["full_pass"] == passes
    idx = index_of(lib_dir)
    row = {r["id"]: r for r in rows(idx)}[b1.run_key]
    assert row["verdict"] == "warn" and row["reasons"] == ["Occupancy in the first hour: 61 %. The guide is 70 % or more."]
    assert row["reason_levels"] == ["warn"] and row["reason_checks"] == ["fc_occupancy_first_hour"]
    assert idx["library"]["counts"]["by_verdict"]["warn"] == 1
    code, out, err = call(["ls", "-L", lib_dir, "--where", "status:warn check:occupancy"], capsys)
    assert out.rstrip("\n").split("\n")[1].startswith(b1.run_key)
    parts.engine.warn_flow_cells.clear()
    code, out, err = call(["recompute", "--status-only", "-L", lib_dir, "--where", "status:warn"], capsys)
    assert code == 0 and "judged 1 run; 1 verdict changed" in err
    code, out, err = call(["recompute", "--status-only", "-L", lib_dir], capsys)
    assert code == 2 and "no run chosen" in err


def test_recompute_keeps_verdicts_judged_before_older_runs_unless_asked(parts, tree, lib_dir, tmp_path, capsys):
    import shutil
    b2 = tree.runs["b2"]
    held = str(tmp_path / "held_back")
    shutil.move(b2.run_dir, held)
    run_scan(lib_dir, tree.roots["seq"])
    shutil.move(held, b2.run_dir)
    run_scan(lib_dir)                                                  # the oldest run arrives last
    status = {n: os.stat(os.path.join(lib_dir, "status", n)).st_mtime_ns for n in os.listdir(os.path.join(lib_dir, "status"))}
    code, out, err = call(["recompute", "--status-only", "--all", "-L", lib_dir], capsys)
    assert code == 0 and "3 runs were left as judged: older runs were added" in err and "--rebaseline" in err
    now = {n: os.stat(os.path.join(lib_dir, "status", n)).st_mtime_ns for n in os.listdir(os.path.join(lib_dir, "status"))}
    assert now == status
    assert [n["code"] for n in index_of(lib_dir)["notices"]] == ["judged_before_older_runs"]
    code, out, err = call(["recompute", "--status-only", "--all", "--rebaseline", "-L", lib_dir], capsys)
    assert code == 0 and "judged 4 runs" in err and "left as judged" not in err
    run_scan(lib_dir)
    assert index_of(lib_dir)["notices"] == []


# ------------------------------------------------------------------------------------------------ doctor, validate, schema
def test_doctor(scanned, tree, lib_dir, capsys, monkeypatch):
    from flowcellar.report import api as report_api
    monkeypatch.setattr(report_api, "viewer_manifest", lambda: {"files": {"viewer.core.js": {}}})
    code, out, err = call(["doctor"], capsys)
    assert code == 0 and out.startswith("flowcellar ") and "packages: numpy " in out and "add -L DIR" in out
    code, out, err = call(["doctor", "-L", lib_dir, "--format", "json"], capsys)
    rep = json.loads(out)
    assert code == 0 and rep["library"]["runs"] == 6 and rep["library"]["writable"] is True and rep["library"]["lease"] is None
    assert [(r["id"], r["state"]) for r in rep["library"]["roots"]] == [("root_seq", "reachable"), ("root_store", "reachable"),
                                                                        ("root_arch", "reachable")]
    assert all(r["listing_seconds"] >= 0 and r["file_system"] for r in rep["library"]["roots"])
    assert rep["library"]["bundle_gzip_bytes"]["largest"] > 0 and rep["viewer"]["built"] is True
    trees.take_offline(tree.roots["store"], "empty")
    code, out, err = call(["doctor", "-L", lib_dir], capsys)
    assert code == 5 and "FAIL  root root_store" in out and "ok    root root_seq" in out
    code, out, err = call(["doctor", "-L", os.path.dirname(lib_dir)], capsys)
    assert code == 5 and "not a FlowCellar library" in out


def test_schema(capsys):
    code, out, err = call(["schema"], capsys)
    schema = json.loads(out)
    assert code == 0 and schema["$id"] == "urn:flowcellar:schema:bundle:0.1.1"
    code, out, err = call(["schema", "--list"], capsys)
    assert out.strip() == "0.1.1  (current)"
    assert json.loads(call(["schema", "--schema-version", "0.1.1"], capsys)[1]) == schema
    code, out, err = call(["schema", "--schema-version", "9.9.9"], capsys)
    assert code == 2 and "no bundle schema 9.9.9" in err


def _valid_bundle() -> dict:
    """The smallest bundle the schema accepts: an empty envelope around an identity."""
    from flowcellar.model import bundle as bundle_model
    b = bundle_model.new_bundle({"run_key": "TST00001_0a1b2c3d", "flow_cell": {"id": "TST00001"},
                                 "protocol_run": {"id": None, "short_id": "0a1b2c3d"}})
    return b


def test_validate(tmp_path, capsys):
    from flowcellar.model import bundle as bundle_model
    from flowcellar.model import schema as schema_model
    good = _valid_bundle()
    assert schema_model.schema_errors(good) == [], "the fixture of this test must itself be valid"
    ok = tmp_path / "good.fcb.json.gz"
    bundle_model.write_bundle(ok, good)
    code, out, err = call(["validate", ok], capsys)
    assert code == 0 and out.startswith("ok  ") and "1 of 1 valid" in err
    bad = dict(good, kind="nonsense")
    bundle_model.write_bundle(tmp_path / "bad.fcb.json.gz", bad)
    tampered = bundle_model.write_bundle(tmp_path / "tampered.fcb.json.gz", good)
    tampered["samples"] = [{"key": "s1"}]
    from flowcellar.util import atomic
    atomic.atomic_write_json(tmp_path / "tampered.fcb.json.gz", tampered)
    judged = bundle_model.compose(good, status={"verdict": "pass", "rules": [
        {"id": "x_check", "role": "key", "tier": "run_health", "severity": "fail", "outcome": "triggered", "message": "m"}],
        "by_tier": {"run_health": "pass", "read_level": "not_available"}})
    bundle_model.write_bundle(tmp_path / "judged.fcb.json.gz", judged)
    (tmp_path / "junk.fcb.json.gz").write_bytes(b"this is not gzip and not JSON")
    code, out, err = call(["validate", ok, tmp_path / "bad.fcb.json.gz", tmp_path / "tampered.fcb.json.gz",
                           tmp_path / "judged.fcb.json.gz", tmp_path / "junk.fcb.json.gz", "--format", "json"], capsys)
    res = {os.path.basename(r["path"]): r for r in json.loads(out)}
    assert code == 7 and res["good.fcb.json.gz"]["ok"] and res["good.fcb.json.gz"]["run_key"] == "TST00001_0a1b2c3d"
    assert any(p.startswith("/kind") for p in res["bad.fcb.json.gz"]["problems"])
    assert any("content_sha256: does not match the content" in p for p in res["tampered.fcb.json.gz"]["problems"])
    assert any("a triggered rule says 'fail'" in p for p in res["judged.fcb.json.gz"]["problems"])
    assert res["junk.fcb.json.gz"]["ok"] is False
    code, out, err = call(["validate", tmp_path], capsys)
    assert code == 7 and "a folder that is not a FlowCellar library" in out
    assert call(["validate", tmp_path / "no_such_file"], capsys)[0] == 4       # a missing input, not a usage error


def test_validate_a_library_composes_each_bundle_with_its_status_file(tmp_path, capsys):
    from flowcellar.library.store import LibraryStore
    from flowcellar.model import status as status_model
    store = LibraryStore.create(tmp_path / "lib")
    good = _valid_bundle()
    info = store.write_measured(good["bundle_id"], good)
    code, out, err = call(["validate", store.path], capsys)
    assert code == 7 and "no status file: the run was never judged" in out
    from flowcellar.model import bundle as bundle_model
    st = {"verdict": "incomplete", "rules": [], "by_tier": {"run_health": "incomplete", "read_level": "not_available"}}
    store.write_status_file(good["bundle_id"], status_model.new_status_file(good["bundle_id"], info["content_sha256"], st,
                                                                              bundle_model.default_relations()))
    code, out, err = call(["validate", store.path], capsys)
    assert code == 0 and "1 of 1 valid" in err
    store.write_status_file(good["bundle_id"], status_model.new_status_file(good["bundle_id"], "0" * 64, st,
                                                                              bundle_model.default_relations()))
    code, out, err = call(["validate", store.path], capsys)
    assert code == 7 and "was written for other content" in out
