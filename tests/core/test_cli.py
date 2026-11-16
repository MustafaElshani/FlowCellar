"""The command group: every command is registered and answers --help; usage errors and exit codes."""
import subprocess
import sys

import pytest
from click.testing import CliRunner

from flowcellar import __version__
from flowcellar.cli import COMMANDS, cli, main

ALL_COMMANDS = {"run", "scan", "hub", "serve", "discover", "ls", "show", "validate", "schema", "doctor", "init", "recompute", "rm",
                  "analysis"}          # analysis: the single-cell and spatial add-on


def test_version_and_help():
    r = CliRunner().invoke(cli, ["--version"])
    assert r.exit_code == 0 and r.output.strip() == f"flowcellar {__version__}"
    r = CliRunner().invoke(cli, ["--help"])
    assert r.exit_code == 0
    for name in ALL_COMMANDS:
        assert f"\n  {name} " in r.output, name


def test_every_command_is_registered_and_has_help():
    assert set(COMMANDS) == ALL_COMMANDS
    for name in COMMANDS:
        r = CliRunner().invoke(cli, [name, "--help"])
        assert r.exit_code == 0 and "Usage: " in r.output, name


def test_no_command_carries_an_option_of_a_feature_that_is_out_of_scope():
    out = "".join(CliRunner().invoke(cli, [name, "--help"]).output for name in COMMANDS)
    for gone in ("--redact", "--mqc", "--allow-edit", "--partials-only", "--minknow-tmp", "--format py"):
        assert gone not in out, gone


@pytest.mark.parametrize("argv", [["run", "{run_dir}"], ["scan", "-L", "{tmp}/lib"], ["schema"], ["ls", "-L", "{tmp}/lib"]])
def test_a_stub_fails_with_a_clear_message(argv, fake_run, tmp_path, capsys):
    argv = [a.format(run_dir=fake_run.run_dir, tmp=tmp_path) for a in argv]
    with pytest.raises(SystemExit) as e:
        main(argv)
    err = capsys.readouterr().err
    if "not built yet" in err:                         # a command that is not in this build says so
        assert e.value.code == 1 and f"flowcellar {argv[0]}" in err
    else:                                              # afterwards the command must not crash on these arguments
        assert e.value.code in (0, 2, 4, 5, 6)


def test_usage_errors_exit_with_code_2(capsys):
    with pytest.raises(SystemExit) as e:
        main(["no-such-command"])
    assert e.value.code == 2
    with pytest.raises(SystemExit) as e:
        main(["serve", "--no-such-option"])
    assert e.value.code == 2


@pytest.mark.parametrize("argv, code", [
    (["serve", "/no/such/library/folder"], 6), (["hub", "/no/such/library/folder"], 6),
    (["run", "/no/such/run/folder", "-o", "/tmp/never-written"], 4), (["discover", "/no/such/folder"], 4),
    (["validate", "/no/such/file.fcb.json.gz"], 4), (["scan", "/no/such/folder", "-L", "/no/such/library"], 4),
])
def test_a_wrong_path_is_an_input_or_library_error_not_a_usage_error(argv, code, capsys):
    """Exit code 4 for an input that is not there, 6 for a library that is not there; 2 is for usage."""
    with pytest.raises(SystemExit) as e:
        main(argv)
    assert e.value.code == code
    assert "no such" in capsys.readouterr().err


def test_console_script_and_module_entry_point():
    out = subprocess.run([sys.executable, "-m", "flowcellar", "--version"], capture_output=True, text=True, check=True)
    assert out.stdout.strip() == f"flowcellar {__version__}"
