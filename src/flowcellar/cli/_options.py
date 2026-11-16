"""Options several commands share."""
from __future__ import annotations

import click

from flowcellar.util.errors import UsageError


def library_option(required: bool = False):
    """``-L DIR`` after the command name (``flowcellar scan ROOT -L DIR``); falls back to the global option and to
    FLOWCELLAR_LIBRARY. Use ``resolve_library(ctx, value)`` to read it."""
    return click.option("-L", "--library", "library", metavar="DIR", type=click.Path(file_okay=False),
                        help="The library folder." + (" Required." if required else ""))


def resolve_library(ctx: click.Context, value: str | None, *, required: bool = True) -> str | None:
    lib = value or (ctx.obj or {}).get("library")
    if required and not lib:
        raise UsageError("no library given", hint="add -L DIR or set FLOWCELLAR_LIBRARY")
    return lib


def format_option(choices=("table", "tsv", "json"), default="table"):
    return click.option("--format", "fmt", type=click.Choice(list(choices)), default=default, show_default=True,
                        help="How to print the result.")
