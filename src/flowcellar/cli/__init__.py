"""The command line: one click group, one file per command.

The group, the global options and the error handling are here. Each command lives in its own module
``flowcellar/cli/<name>.py`` that exposes ``cmd`` (a ``click.Command``) and is loaded when it is called.
No command contacts the network or prompts.
"""
from __future__ import annotations

import importlib
import sys
import traceback

import click

from flowcellar import __version__
from flowcellar.util import log as _log
from flowcellar.util.errors import ExitCode, FlowCellarError

#: command -> module under flowcellar.cli
COMMANDS: dict[str, str] = {name: name for name in (
    "run", "scan", "hub", "serve", "init", "discover", "ls", "show", "rm", "recompute", "doctor", "validate", "schema",
    "analysis",        # the single-cell and spatial add-on (flowcellar.sc), loaded by this command only
)}
_ORDER = ("scan", "run", "analysis", "hub", "serve", "init", "discover", "ls", "show", "rm", "recompute", "doctor", "validate", "schema")


class _Group(click.Group):
    def list_commands(self, ctx):
        return list(_ORDER)

    def get_command(self, ctx, name):
        if name not in COMMANDS:
            return None
        return importlib.import_module(f"flowcellar.cli.{COMMANDS[name]}").cmd


@click.group(cls=_Group, context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 110})
@click.option("-L", "--library", "library", metavar="DIR", envvar="FLOWCELLAR_LIBRARY", type=click.Path(file_okay=False),
              help="The library folder (or the environment variable FLOWCELLAR_LIBRARY).")
@click.option("--config", "config", metavar="FILE", type=click.Path(dir_okay=False),
              help="A configuration file to use instead of the library's flowcellar.toml; for 'run' without a library, the file whose limits judge the run.")
@click.option("--jobs", type=click.IntRange(min=0), default=0, show_default=True, metavar="N",
              help="Worker processes for reading large files; 0 chooses by itself.")
@click.option("--threads", type=click.IntRange(min=1), default=4, show_default=True, metavar="N",
              help="Decompression threads for a BAM without an index.")
@click.option("-v", "--verbose", count=True, help="Say more (repeat for more).")
@click.option("-q", "--quiet", is_flag=True, help="Warnings and errors only.")
@click.option("--log-format", type=click.Choice(["text", "json"]), default="text", show_default=True,
              help="Events on stderr as text lines or as JSON lines.")
@click.version_option(__version__, "--version", prog_name="flowcellar", message="%(prog)s %(version)s")
@click.pass_context
def cli(ctx, library, config, jobs, threads, verbose, quiet, log_format):
    """FlowCellar: one QC library for all your Oxford Nanopore flow cells, wherever their runs are stored.

    \b
      flowcellar scan /place/one /place/two -L ~/flowcellar-library     build or update a library
      flowcellar serve ~/flowcellar-library --open                      look at it
      flowcellar run /path/to/run_folder -o qc/                         one run: a bundle and an offline report

    Nothing is copied from, moved in or written to the folders that are scanned.
    """
    ctx.ensure_object(dict)
    ctx.obj.update(library=library, config=config, jobs=jobs, threads=threads,
                   verbosity=-1 if quiet else verbose, log_format=log_format)
    _log.setup(ctx.obj["verbosity"], log_format)


def main(argv: list[str] | None = None) -> None:
    """Entry point of the ``flowcellar`` console script: maps every outcome to an exit code."""
    rc = None
    try:
        # with standalone_mode=False click RETURNS the code a command gave to ctx.exit(n) instead of raising
        rc = cli.main(args=argv, prog_name="flowcellar", standalone_mode=False)
    except click.exceptions.Exit as e:
        sys.exit(e.exit_code)
    except click.exceptions.Abort:
        click.echo("flowcellar: interrupted", err=True)
        sys.exit(ExitCode.INTERRUPTED)
    except click.ClickException as e:           # usage errors
        e.show()
        sys.exit(ExitCode.USAGE)
    except FlowCellarError as e:
        click.echo(f"flowcellar: {e.message}", err=True)
        if e.hint:
            click.echo(f"  {e.hint}", err=True)
        sys.exit(int(e.exit_code))
    except KeyboardInterrupt:
        click.echo("flowcellar: interrupted", err=True)
        sys.exit(ExitCode.INTERRUPTED)
    except BrokenPipeError:
        sys.exit(ExitCode.OK)
    except Exception as e:  # noqa: BLE001  an unexpected internal error
        if "-v" in sys.argv or "--verbose" in sys.argv:
            traceback.print_exc()
        click.echo(f"flowcellar: internal error: {type(e).__name__}: {e}  (run with -v for the traceback)", err=True)
        sys.exit(ExitCode.INTERNAL)
    sys.exit(rc if isinstance(rc, int) and not isinstance(rc, bool) else ExitCode.OK)
