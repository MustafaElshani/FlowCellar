"""flowcellar rm: remove runs from a library."""
from __future__ import annotations

import click

from flowcellar.cli._options import library_option, resolve_library


@click.command("rm", short_help="Remove runs from a library (moved to its trash; the data folders are never touched).")
@click.argument("run_ids", nargs=-1, required=True)
@library_option(required=True)
@click.option("--purge", is_flag=True, help="Delete instead of moving to trash/.")
@click.pass_context
def cmd(ctx, run_ids, library, purge):
    """Moves bundle, status, card, state, sidecar and report of each run to trash/ and writes a tombstone, so the
    next scan does not add it again. 'flowcellar run RUN_DIR -L DIR' adds a run back on purpose."""
    from flowcellar.library.store import LibraryStore
    from flowcellar.library.trash import remove_runs
    store = LibraryStore.open(resolve_library(ctx, library), config_file=(ctx.obj or {}).get("config"))
    for r in remove_runs(store, list(run_ids), purge=purge):
        click.echo(f"removed {r['run_key']}" + (f"  (kept in {r['trash']})" if r["trash"] else "  (deleted)"), err=True)
