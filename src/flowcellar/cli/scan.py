"""flowcellar scan: build or update a library from any number of roots."""
from __future__ import annotations

import os
import re

import click

from flowcellar.cli._options import library_option, resolve_library
from flowcellar.util.errors import InputError, UsageError


@click.command("scan", short_help="Build or update a library from one or more folders.")
@click.argument("roots", nargs=-1, type=click.Path(file_okay=False))
@library_option(required=True)
@click.option("--name", metavar="TEXT", help="Name of the library, when this scan creates it.")
@click.option("--tier", type=click.Choice(["1", "2", "all"]), default="all", show_default=True,
              help="1: run health only. 2: read data only. all: both, run health first.")
@click.option("--only", multiple=True, metavar="ID", help="Only these run keys or flow cell ids (repeatable).")
@click.option("--since", metavar="DATE", help="Only runs started on or after this date (YYYY-MM-DD).")
@click.option("--time-box", type=click.IntRange(min=0), default=None, metavar="MIN",
              help="Stop starting new passes after this many minutes; 0 = no limit. [default: 240, or [scan] time_box_minutes]")
@click.option("--read", "read", multiple=True, type=click.Choice(["bam", "tables", "analyses"]), metavar="KIND",
              help="Of the read data, read only this kind of large file in this scan (repeatable): bam, tables (annotated read "
                   "tables) or analyses (the per-read tables of single-cell and spatial analyses). [default: all three, in this order]")
@click.option("--quick", is_flag=True, help="Verdicts from the small files only: count no reads and read no read data in this scan.")
@click.option("--wait/--no-wait", default=None, help="Wait for another scan of this library to finish, or exit with code 5. [default: wait when not on a terminal]")
@click.option("--dry-run", is_flag=True, help="Say what would be done and write nothing.")
@click.option("--retry-failed", is_flag=True, help="Try again the inputs whose last pass ended in an error.")
@click.option("--force", is_flag=True, help="Recompute even where nothing changed.")
@click.option("--no-index", is_flag=True, help="Do not rebuild the index.")
@click.option("--no-hub", is_flag=True, help="Do not rewrite the hub.")
@click.pass_context
def cmd(ctx, roots, library, name, tier, only, since, time_box, read, quick, wait, dry_run, retry_failed, force, no_index, no_hub):
    """Find the runs under ROOTS, add them to the library DIR and judge them. Creates the library when DIR does
    not exist and remembers the roots; without ROOTS the remembered ones are scanned again.

    \b
      flowcellar scan /data/minknow /data/projects -L ~/flowcellar-library     the first scan: creates the library
      flowcellar scan -L ~/flowcellar-library                                  later: looks for what is new
      flowcellar scan /another/folder -L ~/flowcellar-library                  adds a folder and scans all of them

    Nothing is copied from, moved in or written to ROOTS. Verdicts from the small files come first; the reads
    are counted afterwards, newest runs first, inside the time box, and a later scan goes on where this one
    stopped. Exit code 5: finished, but a root or a file was skipped.
    """
    from flowcellar.library.scan import ScanOptions, scan
    lib = resolve_library(ctx, library)
    if since and not re.fullmatch(r"\d{4}-\d{2}-\d{2}", since):
        raise UsageError(f"--since {since}: expected a date such as 2026-03-01")
    for r in roots:
        if not os.path.isdir(os.path.expanduser(r)):
            raise InputError(f"{r}: no such folder", hint="a folder to scan must exist when it is first given")
    obj = ctx.obj or {}
    res = scan(lib, ScanOptions(roots=tuple(roots), name=name, tier=tier, only=tuple(only), since=since,
                                time_box_minutes=time_box, read=tuple(read), quick=quick, wait=wait, dry_run=dry_run, retry_failed=retry_failed, force=force,
                                no_index=no_index, no_hub=no_hub, jobs=obj.get("jobs") or 0,
                                threads=obj.get("threads") or 0, config_file=obj.get("config")))
    if res.exit_code:
        raise SystemExit(res.exit_code)   # not ctx.exit: the group would turn it into exit code 0
