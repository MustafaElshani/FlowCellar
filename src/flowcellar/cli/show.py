"""flowcellar show: inspect one run of a library."""
from __future__ import annotations

import json

import click

from flowcellar.cli._options import library_option, resolve_library
from flowcellar.util.errors import InputError


@click.command("show", short_help="Inspect one run: verdict, checks, inputs or one section.")
@click.argument("run_id")
@library_option(required=True)
@click.option("--inputs", "what", flag_value="inputs", help="The input files with role, state and link.")
@click.option("--status", "what", flag_value="status", help="Every check with its outcome and sentence.")
@click.option("--section", metavar="ID", help="The summary of one section, for example fc or len.")
@click.option("--json", "as_json", is_flag=True, help="Print the composed bundle (or the chosen part) as JSON.")
@click.pass_context
def cmd(ctx, run_id, library, what, section, as_json):
    """RUN_ID is a run key, or a flow cell id when the library holds one run of it; or the id of an analysis
    (an_...), as "flowcellar ls --analyses" lists them."""
    from flowcellar.library.store import LibraryStore
    from flowcellar.library.trash import resolve_run
    from flowcellar.model.status import VERDICT_WORD
    store = LibraryStore.open(resolve_library(ctx, library), config_file=(ctx.obj or {}).get("config"))
    if run_id.startswith("an_"):                    # an analysis of the single-cell and spatial add-on
        from flowcellar.sc import listing
        listing.show(store.path, run_id, what=what, section=section, as_json=as_json)
        return
    key = resolve_run(store, run_id)
    bundle = store.composed(key)
    state = store.read_state(key) or {}
    status = bundle.get("status") or {}
    if section:
        mods = bundle.get("modules") or {}
        if section not in mods:
            raise InputError(f"{key} has no section '{section}'", hint="sections: " + (", ".join(sorted(mods)) or "none yet"))
        part = mods[section]
        if as_json:
            click.echo(json.dumps(part, indent=1))
            return
        click.echo(f"{key}  section {section}  state {part.get('state', 'ok')}"
                   + (f" ({part.get('partial_reason')})" if part.get("partial_reason") else ""))
        for k, v in (part.get("summary") or {}).items():
            click.echo(f"  {k:<40} {v}")
        return
    if what == "inputs":
        rows = [{"role": i.get("role"), "state": i.get("state"), "path": i.get("path"), "size": i.get("size"),
                 "link": (i.get("link") or {}).get("method"), "confidence": (i.get("link") or {}).get("confidence"),
                 "counted_for": i.get("counted_for"), "pass": None}
                for i in (bundle.get("provenance") or {}).get("inputs") or []]
        passes = {i.get("path"): i.get("pass") for i in state.get("inputs") or []}
        seen = {r["path"] for r in rows}
        for r in rows:
            r["pass"] = passes.get(r["path"])
        rows += [{"role": i.get("role"), "state": i.get("state") or "known to the library", "path": i.get("path"), "size": i.get("size"),
                  "link": None, "confidence": None, "counted_for": None, "pass": i.get("pass")}
                 for i in state.get("inputs") or [] if i.get("path") not in seen]
        if as_json:
            click.echo(json.dumps(rows, indent=1))
            return
        for r in rows:
            size = f"{r['size'] / 1e6:,.1f} MB" if isinstance(r["size"], (int, float)) else ""
            link = f"  linked by {r['link']} ({r['confidence']})" if r["link"] else ""
            read = f"  pass: {r['pass']}" if r["pass"] and r["pass"] != "none" else ""
            click.echo(f"{(r['role'] or ''):<28} {(r['state'] or ''):<10} {size:>12}  {r['path']}{link}{read}")
        return
    if what == "status":
        if as_json:
            click.echo(json.dumps(status, indent=1))
            return
        for r in status.get("rules") or []:
            mark = {"triggered": (r.get("severity") or "").upper() or "TRIGGERED", "ok": "ok", "not_evaluated": "not assessed"}.get(r.get("outcome"), "")
            click.echo(f"{mark:<12} {(r.get('tier') or ''):<10} {(r.get('role') or ''):<10} {r.get('id')}")
            if r.get("message"):
                click.echo(f"{'':<12} {r['message']}")
        if not status.get("rules"):
            click.echo("no checks stored for this run")
        return
    if as_json:
        click.echo(json.dumps(bundle, indent=1))
        return
    card = (store.card(key) or {}).get("columns") or {}
    tiers = status.get("by_tier") or {}
    click.echo(key + "  " + " · ".join(str(x) for x in (card.get("device_type"), card.get("product"), card.get("kit"), card.get("assay"),
                                                       f"{card['hours']:g} h" if isinstance(card.get("hours"), (int, float)) else None) if x))
    word = VERDICT_WORD.get(status.get("verdict"), "Not assessed") + (" so far" if status.get("provisional") else "")
    click.echo(f"  verdict      {word}")
    click.echo(f"  Run health   {VERDICT_WORD.get(tiers.get('run_health'), 'Not assessed')}")
    click.echo(f"  Read data    {VERDICT_WORD.get(tiers.get('read_level'), 'Not assessed')}")
    if status.get("headline"):
        click.echo(f"  {status['headline']}")
    if status.get("detail"):
        click.echo(f"  {status['detail']}")
    for m in status.get("modifiers") or []:
        click.echo(f"  note: {m.get('message')}")
    for r in status.get("rules") or []:
        if r.get("outcome") == "triggered" and r.get("severity") in ("warn", "fail"):
            click.echo(f"  {r['severity']:<5} {r.get('message')}")
    rel = bundle.get("relations") or {}
    for p in rel.get("projects") or []:
        click.echo(f"  project      {p.get('name')}  (found under: {', '.join(p.get('roots') or []) or '-'})")
    for p in state.get("paths") or []:
        root = next((r for r in store.config["roots"] if r["id"] == p.get("root")), None)
        where = (root["path"].rstrip("/") + "/" + p["relpath"]) if root and p.get("relpath") else (p.get("path") or "")
        click.echo(f"  folder       {where}" + ("  (no longer there)" if state.get("source_state") == "missing" else ""))
    click.echo(f"  bundle       {store.layout.bundle(key)}")
    click.echo("  more: --status (every check), --inputs (files), --section ID, --json")
