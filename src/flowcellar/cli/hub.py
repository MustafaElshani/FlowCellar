"""flowcellar hub: (re)build the static hub of a library."""
from __future__ import annotations

import os
import pathlib
import webbrowser

import click


@click.command("hub", short_help="Rebuild the static hub (hub/index.html) of a library.")
@click.argument("library", type=click.Path(file_okay=False))
@click.option("--no-sidecars", is_flag=True, help="Write the page and the index only; runs then open under 'flowcellar serve' only.")
@click.option("--reports", is_flag=True, help="Also write one offline report per run into reports/ (not for analyses: 'flowcellar analysis DIR -o OUT' writes theirs).")
@click.option("--open", "open_", is_flag=True, help="Open the hub in the default browser.")
@click.pass_context
def cmd(ctx, library, no_sidecars, reports, open_):
    """Write LIBRARY/hub/index.html, the index and one data file per run. Scan does this by itself; this command
    only rebuilds the page."""
    from flowcellar.library import api as library_api
    from flowcellar.report import api as report_api

    report_api.viewer_manifest()                      # one clear line when the viewer is not built, before anything is opened
    lib = library_api.open_library(library)
    result = report_api.write_hub(lib, sidecars=not no_sidecars, reports=reports)
    hub_dir = os.path.dirname(result["index_html"])
    in_data = [p for p in result["written"] if os.path.dirname(p) == os.path.join(hub_dir, "data")]
    n_reports = sum(1 for p in result["written"] if p.endswith(".flowcellar.html"))
    n_runs = len(lib.run_keys())
    parts = [f"{n_runs} run{'' if n_runs == 1 else 's'}"]
    if no_sidecars:
        parts.append("no data files: open it with 'flowcellar serve'")
    else:
        parts.append(f"{len(in_data)} data file{'' if len(in_data) == 1 else 's'} written")
    if result["removed"]:
        parts.append(f"{len(result['removed'])} removed")
    if reports:
        parts.append(f"{n_reports} report{'' if n_reports == 1 else 's'} written")
    click.echo(f"hub ready: {result['index_html']}     ({'; '.join(parts)})")
    if reports:
        # a report file is written for runs only; say so when the library also holds analyses, instead of leaving them out silently
        n_an = len((lib.index() or {}).get("analyses") or [])
        if n_an:
            click.echo(f"no report file was written for the {n_an} single-cell and spatial analys{'is' if n_an == 1 else 'es'}: "
                       "'flowcellar analysis ANALYSIS_DIR -o OUT' writes the report of one analysis.")
    if open_:
        webbrowser.open(pathlib.Path(result["index_html"]).resolve().as_uri())
