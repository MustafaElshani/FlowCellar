"""flowcellar recompute: judge runs again without reading any input."""
from __future__ import annotations

import click

from flowcellar.cli._options import library_option, resolve_library
from flowcellar.util.errors import UsageError


@click.command("recompute", short_help="Judge runs again with the current limits; reads no input.")
@click.argument("run_ids", nargs=-1)
@library_option(required=True)
@click.option("--status-only", is_flag=True, help="Rewrite status files and cards only (the one mode of this build).")
@click.option("--rebaseline", is_flag=True, help="Also take the history as the library holds it now.")
@click.option("--all", "all_", is_flag=True, help="Every run of the library.")
@click.option("--where", metavar="EXPR", help="Runs matching this search.")
@click.pass_context
def cmd(ctx, run_ids, library, status_only, rebaseline, all_, where):
    """Each run is judged against its history as of its own start; no bundle is rewritten.

    Use it after changing limits in flowcellar.toml. A run that was judged before older runs were added keeps
    its verdict unless --rebaseline is given.
    """
    from flowcellar.library import query
    from flowcellar.library.maintain import recompute_status
    from flowcellar.library.store import LibraryStore
    from flowcellar.library.trash import resolve_run
    from flowcellar.util.timeutil import now_z
    if not status_only:
        raise UsageError("recompute needs --status-only: judging again is the only mode of this version",
                         hint="flowcellar recompute --status-only --all -L DIR")
    store = LibraryStore.open(resolve_library(ctx, library), config_file=(ctx.obj or {}).get("config"))
    if all_:
        keys = store.run_keys()
    elif where:
        keys = [r["id"] for r in query.filter_rows(query.rows_with_samples(store.index()), where, now=now_z())]
    else:
        keys = [resolve_run(store, r) for r in run_ids]
    if not keys:
        raise UsageError("no run chosen", hint="give run keys, --where EXPR or --all")
    out = recompute_status(store, keys, rebaseline=rebaseline)
    for c in out["changed"]:
        click.echo(f"{c['run_key']}   {c['was']} -> {c['now']}", err=True)
    for k, msg in out["failed"].items():
        click.echo(f"{k}   not judged ({msg})", err=True)
    click.echo(f"judged {len(out['judged'])} run{'s' if len(out['judged']) != 1 else ''}; {len(out['changed'])} verdict"
               f"{'s' if len(out['changed']) != 1 else ''} changed; no input was read and no bundle was written", err=True)
    if out["held"]:
        click.echo(f"{len(out['held'])} run{'s were' if len(out['held']) != 1 else ' was'} left as judged: older runs were added "
                   "to the library since. Add --rebaseline to judge them against the library as it is now.", err=True)
    if out["failed"]:
        raise SystemExit(5)   # not ctx.exit: the group would turn it into exit code 0
