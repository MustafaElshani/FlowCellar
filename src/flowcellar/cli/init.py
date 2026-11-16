"""flowcellar init: create a library folder and draft its configuration."""
from __future__ import annotations

import os

import click

from flowcellar.util.errors import InputError, LibraryError


@click.command("init", short_help="Create an empty library and draft its flowcellar.toml.")
@click.argument("library", type=click.Path(file_okay=False))
@click.option("--from-scan", "from_scan", multiple=True, metavar="ROOT", type=click.Path(file_okay=False),
              help="Walk this folder and propose roots and file patterns from what is there (repeatable). Computes nothing.")
@click.option("--name", metavar="TEXT", help="Name of the library.")
@click.pass_context
def cmd(ctx, library, from_scan, name):
    """Optional: 'flowcellar scan' creates a library by itself. Use init to shape flowcellar.toml first.

    Writes LIBRARY/flowcellar.toml with every section and its default. With --from-scan the roots are filled in
    and the read-data file patterns that the folders suggest are added as comments to check. Then edit the
    file and run: flowcellar scan -L LIBRARY
    """
    from flowcellar.library.draft import draft_config
    from flowcellar.library.store import LibraryStore, is_library
    for r in from_scan:
        if not os.path.isdir(os.path.expanduser(r)):
            raise InputError(f"--from-scan {r}: no such folder")
    path = os.path.abspath(os.path.expanduser(library))
    if os.path.isdir(path) and is_library(path):
        raise LibraryError(f"{path}: is a library already", hint="edit its flowcellar.toml; init never overwrites one")
    text, facts = draft_config(name or os.path.basename(path.rstrip("/")) or "library", list(from_scan))
    store = LibraryStore.create(path, name=name, config_text=text)
    click.echo(f"created library {store.path}", err=True)
    if from_scan:
        click.echo(f"  looked at {facts['roots']} folder{'s' if facts['roots'] != 1 else ''}: {facts['runs']} run folders with "
                   f"MinKNOW reports, {facts['raw_only']} with raw data only; {facts['proposals']} file pattern"
                   f"{'s' if facts['proposals'] != 1 else ''} proposed as comments", err=True)
        for rid in facts.get("offline") or []:
            click.echo(f"  root '{rid}' did not answer and was not looked at", err=True)
    click.echo(f"  edit {store.layout.config}, then: flowcellar scan -L {store.path}", err=True)
