"""flowcellar validate: check bundles against the schema and the invariants."""
from __future__ import annotations

import json
import os

import click

from flowcellar.cli._options import format_option, library_option, resolve_library
from flowcellar.util.errors import InputError, UsageError


@click.command("validate", short_help="Check bundle files, or every bundle of a library.")
@click.argument("paths", nargs=-1, type=click.Path())
@library_option()
@format_option(choices=("text", "json"), default="text")
@click.pass_context
def cmd(ctx, paths, library, fmt):
    """PATHS are bundle files (*.fcb.json.gz), reports are not read. A library folder (as a path, or with -L)
    checks all its bundles with their status files composed in. Exit code 7 when anything fails."""
    from flowcellar.library.maintain import validate_paths
    paths = list(paths)
    lib = resolve_library(ctx, library, required=False)
    if lib and lib not in paths and (library or not paths):
        paths.append(lib)                     # -L DIR, or FLOWCELLAR_LIBRARY when nothing else was named
    if not paths:
        raise UsageError("nothing to check", hint="give bundle files or a library folder: flowcellar validate -L DIR")
    for p in paths:
        if not os.path.exists(p):
            raise InputError(f"{p}: no such file or folder")
    results = validate_paths(paths)
    bad = [r for r in results if not r["ok"]]
    if fmt == "json":
        click.echo(json.dumps(results, indent=1))
    else:
        for r in results:
            click.echo(f"{'ok  ' if r['ok'] else 'FAIL'}  {r['path']}")
            for p in r["problems"][:20]:
                click.echo(f"        {p}")
            if len(r["problems"]) > 20:
                click.echo(f"        ... and {len(r['problems']) - 20} more")
        click.echo(f"{len(results) - len(bad)} of {len(results)} valid", err=True)
    if bad:
        raise SystemExit(7)   # not ctx.exit: the group would turn it into exit code 0
