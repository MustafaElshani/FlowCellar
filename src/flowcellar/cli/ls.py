"""flowcellar ls: list the runs of a library."""
from __future__ import annotations

import json

import click

from flowcellar.cli._options import format_option, library_option, resolve_library
from flowcellar.util.errors import UsageError

DEFAULT_COLUMNS = ("id", "started", "verdict", "headline")


def _text(v) -> str:
    if v is None:
        return ""
    if isinstance(v, bool):
        return "yes" if v else "no"
    if isinstance(v, float):
        return f"{v:.6g}"
    if isinstance(v, list):
        return ",".join(_text(x) for x in v)
    return str(v)


def _tsv(v) -> str:
    s = _text(v).replace("\t", " ").replace("\n", " ")
    return "'" + s if s[:1] in ("=", "+", "-", "@") else s


@click.command("ls", short_help="List the runs of a library.")
@library_option(required=True)
@click.option("--where", metavar="EXPR", help="Filter in the hub's search grammar, for example 'status:warn fc:TST00001'.")
@click.option("--columns", metavar="LIST", help="Comma-separated column ids of the index.")
@click.option("--sort", metavar="COLUMN", help="Column to sort by; a leading '-' reverses. [default: -started]")
@click.option("--analyses", is_flag=True, help="List the single-cell and spatial analyses of the library instead of its runs.")
@format_option()
@click.pass_context
def cmd(ctx, library, where, columns, sort, fmt, analyses):
    """One line per run, newest first: run key, start, verdict, the stored headline.

    \b
      flowcellar ls -L DIR --where 'status:warn,fail after:2026-03'
      flowcellar ls -L DIR --where 'project:alpha reads:>1M' --columns id,started,reads,gb,project --format tsv
    """
    from flowcellar.library import query
    from flowcellar.library.store import LibraryStore
    from flowcellar.model.status import VERDICT_WORD
    from flowcellar.util.timeutil import now_z
    store = LibraryStore.open(resolve_library(ctx, library), config_file=(ctx.obj or {}).get("config"))
    index = store.index()
    if analyses:
        from flowcellar.sc import listing
        listing.list_analyses(index, columns=columns, fmt=fmt)
        return
    rows = query.filter_rows(query.rows_with_samples(index), where, now=now_z())
    known = list(index["runs"]["columns"])
    cols = [c.strip() for c in columns.split(",") if c.strip()] if columns else list(DEFAULT_COLUMNS)
    bad = [c for c in cols if c not in known]
    if bad:
        raise UsageError(f"no such column: {', '.join(bad)}", hint="columns: " + ", ".join(known))
    key = sort or "-started"
    reverse, key = key.startswith("-"), key.lstrip("-")
    if key not in known:
        raise UsageError(f"--sort {key}: no such column", hint="columns: " + ", ".join(known))
    have = [r for r in rows if r.get(key) is not None]
    try:
        have.sort(key=lambda r: r[key], reverse=reverse)
    except TypeError:
        have.sort(key=lambda r: _text(r[key]), reverse=reverse)
    rows = have + [r for r in rows if r.get(key) is None]
    if fmt == "json":
        click.echo(json.dumps([{c: r.get(c) for c in cols} for r in rows], indent=1))
        return
    if fmt == "tsv":
        click.echo("\t".join(cols))
        for r in rows:
            click.echo("\t".join(_tsv(r.get(c)) for c in cols))
        return

    def show(c: str, v) -> str:
        if c == "started" and v:
            return str(v)[:16].replace("T", " ") + " UTC"       # as stored; the pages show the reader's own clock
        if c in ("verdict", "health", "data") and v:
            return VERDICT_WORD.get(v, str(v))
        return _text(v)
    table = [[show(c, r.get(c)) for c in cols] for r in rows]
    widths = [max([len(c)] + [len(t[i]) for t in table]) for i, c in enumerate(cols)]
    widths[-1] = 0
    click.echo("  ".join(c.ljust(w) for c, w in zip(cols, widths)).rstrip())
    for t in table:
        click.echo("  ".join(x.ljust(w) for x, w in zip(t, widths)).rstrip())
    total = index["runs"]["n"]
    click.echo(f"{len(rows)} of {total} run{'s' if total != 1 else ''}" if where else f"{total} run{'s' if total != 1 else ''}", err=True)
