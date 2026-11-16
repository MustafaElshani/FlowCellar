"""The places where two parts of the package meet: each of these was a fault found when the parts were put together."""
import click
import pytest

from flowcellar import cli as cli_mod
from flowcellar.library.projects import NO_PROJECT, merge_projects, project_id
from flowcellar.model import axes, runfolder, schema


def test_an_exit_code_given_with_ctx_exit_reaches_the_shell(monkeypatch):
    """With standalone_mode off click RETURNS the code of ctx.exit(n); the entry point must not turn it into 0."""
    @click.command("probe")
    @click.argument("code", type=int)
    @click.pass_context
    def probe(ctx, code):
        ctx.exit(code)

    monkeypatch.setattr(cli_mod._Group, "get_command", lambda self, ctx, name: probe if name == "probe" else None)
    for code in (0, 3, 5):
        with pytest.raises(SystemExit) as e:
            cli_mod.main(["probe", str(code)])
        assert e.value.code == code


def test_the_folder_above_a_run_is_not_a_raw_only_run():
    run_name = "20260302_0900_X_TST00001_0a1b2c3d"
    names = ["fast5_to_pod5.sh", run_name, "notes.txt"]
    # names alone cannot tell: a script named like raw data beside a folder named like a run
    assert runfolder.classify_listing("some_project", names) == "run_raw_only"
    # with the kinds the rule holds: a raw marker is a folder or a raw file, the key is the folder's own or a raw file's
    assert runfolder.classify_listing("some_project", names, dirs=[run_name]) is None
    assert runfolder.classify_listing(run_name, ["pod5", "other_reports"], dirs=["pod5", "other_reports"]) == "run_raw_only"
    assert runfolder.classify_listing("moved_here", ["pod5"], dirs=["pod5"]) is None, "no run key anywhere"
    assert runfolder.classify_listing("x", ["final_summary_TST00001_0a1b2c3d_4e5f6a7b.txt"], dirs=[]) == "run"
    assert runfolder.classify_listing(run_name, ["bam_pass", "bam_fail"], dirs=["bam_pass", "bam_fail"]) == "lookalike"


def test_a_project_id_is_a_safe_id_whatever_the_folder_is_called():
    assert project_id("Alpha Project (2026)") == "alpha_project_2026"
    assert project_id("_drafts_") == "drafts" and project_id("__") == "project" and project_id("--") != NO_PROJECT
    for name in ("_x", "x_", "9 lives", "ÅÄÖ", "_", "a" * 200 + "_"):
        pid = project_id(name)
        assert pid[0].isalnum() and pid[-1].isalnum() and len(pid) <= 80, (name, pid)
    merged = merge_projects([("_Same", "r1"), ("same_", "r2"), (None, "r3")], ["r2", "r1"])
    assert [(p["id"], p["roots"]) for p in merged] == [("same", ["r2", "r1"])]


def test_every_axis_of_the_registry_is_valid_under_the_schema():
    """The registry may carry ``note`` and ``closed_last``; the schema must allow what it writes."""
    axis_schema = {"$schema": schema.load_schema()["$schema"], "$ref": "#/$defs/axis", "$defs": schema.load_schema()["$defs"]}
    import jsonschema
    validator = jsonschema.Draft202012Validator(axis_schema)
    carried = set()
    for axis_id in axes.AXES:
        definition = axes.axis_def(axis_id)
        problems = [e.message for e in validator.iter_errors(definition)]
        assert problems == [], (axis_id, problems)
        carried |= set(definition)
    assert {"note", "closed_last"} <= carried, "the registry no longer uses them: the test would prove nothing"
