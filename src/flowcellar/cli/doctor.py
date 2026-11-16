"""flowcellar doctor: check an installation and a library."""
from __future__ import annotations

import json

import click

from flowcellar.cli._options import format_option, library_option, resolve_library


@click.command("doctor", short_help="Check the installation, a library and its roots.")
@library_option()
@format_option(choices=("text", "json"), default="text")
@click.pass_context
def cmd(ctx, library, fmt):
    """Versions and optional packages; the built viewer and its hashes; for a library: marker, layout, lease,
    whether each root answers, file-system type, and the time of one directory listing per root."""
    from flowcellar.library.maintain import doctor
    rep = doctor(resolve_library(ctx, library, required=False))
    if fmt == "json":
        click.echo(json.dumps(rep, indent=1))
    else:
        click.echo(f"flowcellar {rep['flowcellar']} · bundle schema {rep['schema']} · Python {rep['python']} · {rep['platform']}")
        click.echo("packages: " + ", ".join(f"{k} {v}" if v else f"{k} (not installed)" for k, v in rep["packages"].items()))
        for c in rep["checks"]:
            mark = "ok  " if c["ok"] else "note" if c["ok"] is None else "FAIL"
            click.echo(f"  {mark}  {c['name']}: {c['message']}")
        if "library" not in rep:
            click.echo("(add -L DIR to check a library and its roots)")
    if any(c["ok"] is False for c in rep["checks"]):
        raise SystemExit(5)   # not ctx.exit: the group would turn it into exit code 0
