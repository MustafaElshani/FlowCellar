"""flowcellar discover: say what a scan would find, and why."""
from __future__ import annotations

import json
import os
from collections import Counter

import click

from flowcellar.cli._options import format_option, library_option, resolve_library
from flowcellar.util.errors import FlowCellarError, InputError, NotBuiltError, UsageError


def _cell(s) -> str:
    """Text for a TSV cell: no tabs or line breaks, and nothing a spreadsheet would take for a formula."""
    s = "" if s is None else str(s).replace("\t", " ").replace("\n", " ")
    return "'" + s if s[:1] in ("=", "+", "-", "@") else s


@click.command("discover", short_help="List the runs and files a scan would find; computes nothing.")
@click.argument("roots", nargs=-1, type=click.Path(file_okay=False))
@library_option()
@click.option("--explain", is_flag=True, help="For every file: the run it would be attached to, by which evidence.")
@format_option()
@click.pass_context
def cmd(ctx, roots, library, explain, fmt):
    """Walk ROOTS (or the roots of the library given with -L) and print run folders, raw-only folders, look-alikes
    and read-data files. Nothing is read beyond identity evidence and nothing is written.

    With -L the library's flowcellar.toml gives the file patterns and the project rule of each root; ROOTS then
    limits the walk to those folders.
    """
    from flowcellar import discovery
    from flowcellar.config import loader
    from flowcellar.engine import api as engine_api
    from flowcellar.library.projects import project_name
    from flowcellar.library.store import LibraryStore
    from flowcellar.model.runkey import run_key_from_folder_name
    from flowcellar.readdata import api as readdata_api
    lib = resolve_library(ctx, library, required=False)
    config_file = (ctx.obj or {}).get("config")
    for r in roots:
        if not os.path.isdir(os.path.expanduser(r)):
            raise InputError(f"{r}: no such folder")
    if lib:
        cfg = LibraryStore.open(lib, config_file=config_file).config
    else:
        cfg = loader.load(config_file).config
    if roots:
        want = [loader.normalise_root_path(r) for r in roots]
        known = {r["path"]: r for r in cfg["roots"]}
        new = [p for p in want if p not in known]
        ids = loader.root_ids_for(new, {r["id"] for r in cfg["roots"]})
        extra = {p: {"id": i, "path": p, "live": False, "project_depth": 1, "project_pattern": None} for p, i in zip(new, ids)}
        cfg["roots"] = [known.get(p) or extra[p] for p in want]
    if not cfg["roots"]:
        raise UsageError("nothing to walk", hint="give one or more folders, or a library with -L")
    found = discovery.discover(cfg)
    by_id = {r["id"]: r for r in cfg["roots"]}

    runs, identities, note_engine = [], [], None
    for root_id, f in found.runs:
        row = {"kind": f.kind, "run_key": None, "root": root_id, "project": project_name(by_id[root_id], f.relpath),
               "state": None, "started": None, "path": f.path, "relpath": f.relpath, "why": ", ".join(f.evidence),
               "key_from": "files"}
        try:
            ident = engine_api.identify(f.path)
            row.update(run_key=ident.run_key, state=ident.run_state, started=ident.started, kind=ident.kind)
            identities.append(ident)
        except NotBuiltError as e:
            note_engine = e.message
            row.update(run_key=run_key_from_folder_name(f.path), key_from="folder name")
        except FlowCellarError as e:
            row.update(run_key=run_key_from_folder_name(f.path), key_from="folder name", state=f"not read: {e.message}")
        runs.append(row)
    links = {}
    if explain and found.candidates:
        try:
            for li in readdata_api.link(found.candidates, identities, config=cfg, probes_dir=None):
                links[li.candidate.path] = li
        except FlowCellarError as e:
            click.echo(f"flowcellar: the files could not be linked ({e.message})", err=True)
    files = []
    for c in found.candidates:
        li = links.get(c.path)
        files.append({"role": c.role, "root": c.root_id, "project": project_name(by_id.get(c.root_id), c.relpath, loader.not_project_names(cfg)),
                      "sample": c.sample, "variant": c.variant, "size": c.size, "path": c.path, "rule": c.rule,
                      "runs": li.run_keys if li else None,
                      "link": (li.link or {}).get("method") if li else None,
                      "confidence": (li.link or {}).get("confidence") if li else None,
                      "evidence": (li.link or {}).get("evidence") if li else None,
                      "counted_for": li.counted_for if li else None,
                      "not_counted": li.not_counted_reason if li else None})
    roots_out = [{"id": rid, "path": o.path, "state": o.state, "reason": o.reason, "seconds": round(o.seconds, 3),
                  "listings": o.walk.counts.get("listings") if o.walk else None,
                  "anchors": len(o.walk.anchors) if o.walk else None,
                  "over_budget": o.walk.over_budget if o.walk else []} for rid, o in found.roots.items()]
    looks = [{"root": rid, "path": f.path, "why": "named like a run folder, holds only " + ", ".join(f.evidence)}
             for rid, f in found.lookalikes]
    counts = found.counts()
    if fmt == "json":
        click.echo(json.dumps({"roots": roots_out, "runs": runs, "lookalikes": looks, "files": files, "counts": counts,
                               "seconds": round(found.seconds, 3), "problems": found.problems}, indent=1))
    elif fmt == "tsv":
        click.echo("\t".join(["what", "kind_or_role", "run_key", "root", "project", "state", "path", "why"]))
        for r in runs:
            click.echo("\t".join(_cell(x) for x in ("run", r["kind"], r["run_key"], r["root"], r["project"], r["state"], r["path"], r["why"])))
        for x in looks:
            click.echo("\t".join(_cell(v) for v in ("lookalike", "basecaller_output", "", x["root"], "", "", x["path"], x["why"])))
        for f in files:
            why = f["rule"] or ""
            if f["link"]:
                why += f"; linked by {f['link']} ({f['confidence']}): {f['evidence']}"
            click.echo("\t".join(_cell(v) for v in ("file", f["role"], ",".join(f["runs"] or []), f["root"], f["project"], "", f["path"], why)))
    else:
        for r in roots_out:
            if r["state"] != "reachable":
                click.echo(f"root {r['id']}: {r['state'].upper()}: {r['reason']}  ({r['path']})")
            else:
                click.echo(f"root {r['id']}: {r['path']}  ({r['listings']:,} folders listed in {r['seconds']:.1f} s"
                           + (f", {r['anchors']} anchors" if r["anchors"] else "") + ")")
            for sub in r["over_budget"]:
                click.echo(f"  given up: {sub} (more folders than subtree_listing_budget without a run folder)")
        click.echo("")
        w = max([len(r["run_key"] or "?") for r in runs] + [7])
        for r in sorted(runs, key=lambda r: (r["started"] or "", r["run_key"] or ""), reverse=True):
            kind = "run" if r["kind"] == "run" else "raw only"
            click.echo(f"{(r['run_key'] or '?'):<{w}}  {kind:<8}  {(r['state'] or ''):<11}  {(r['project'] or '-'):<24}  "
                       f"{r['root']}:{r['relpath'] or '.'}")
            if explain:
                click.echo(f"{'':<{w}}  why: holds {r['why']}" + ("" if r["key_from"] == "files" else "; run key read from the folder name"))
        for x in looks:
            click.echo(f"{'not a run':<{w}}  {x['path']}")
            if explain:
                click.echo(f"{'':<{w}}  why: {x['why']}")
        if files:
            click.echo("")
            if explain:
                for f in files:
                    click.echo(f"{f['role']:<24}  {f['path']}")
                    line = f"{'':<24}  why: {f['rule']}"
                    if f["sample"]:
                        line += f"; sample {f['sample']}" + (f", variant {f['variant']}" if f["variant"] else "")
                    click.echo(line)
                    if f["runs"] is not None:
                        if f["runs"]:
                            click.echo(f"{'':<24}  run: {', '.join(f['runs'])} by {f['link']} ({f['confidence']}): {f['evidence'] or ''}")
                        else:
                            click.echo(f"{'':<24}  run: not linked to any run")
                        if f["not_counted"]:
                            click.echo(f"{'':<24}  not counted: {f['not_counted']}")
            else:
                for role, n in sorted(Counter(f["role"] for f in files).items()):
                    size = sum(f["size"] for f in files if f["role"] == role)
                    click.echo(f"{role:<24}  {n:>6,} files  {size / 1e9:>10,.1f} GB")
                click.echo("(--explain lists every file with the rule that found it and the run it belongs to)")
        click.echo("")
        click.echo(f"{counts['runs']} run folders with MinKNOW reports, {counts['raw_only']} with raw data only, "
                   f"{counts['lookalikes']} named like a run but not one, {counts['files']:,} read-data files; "
                   f"{counts['listings']:,} folders listed in {found.seconds:.1f} s. Nothing was computed or written.")
        if counts.get("analyses"):                  # anchor files of the single-cell and spatial add-on; a scan reads them
            click.echo(f"{counts['analyses']} single-cell or spatial analysis folders (a folder that is a copy of another counts once in a scan).")
        for p in found.problems:
            click.echo(f"pattern not used: {p}")
        if note_engine:
            click.echo(f"note: run keys were read from folder names ({note_engine})")
    if any(o.state != "reachable" for o in found.roots.values()):
        raise SystemExit(5)   # not ctx.exit: the group would turn it into exit code 0
