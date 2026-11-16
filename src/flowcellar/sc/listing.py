"""``flowcellar ls --analyses`` and ``flowcellar show <analysis id>``: the analyses of a library on the command line."""
from __future__ import annotations

import json

import click

from flowcellar.model.status import VERDICT_WORD
from flowcellar.sc import ASSAY_WORD, MODULE_ID, PIPELINE_WORD, is_analysis_id
from flowcellar.sc import store as store_mod
from flowcellar.util.errors import InputError, UsageError

DEFAULT_COLUMNS = ("id", "date", "assay", "verdict", "label", "headline")

#: how a link from an analysis to a run was made, in words (the same words as the analysis page)
LINK_WORD = {"read_tag": "matched by the file names stored on the reads", "bam_header_rg": "matched by the read group in the BAM header",
             "sibling_analysis": "taken from the ONT workflow analysis of the same reads", "read_id_sketch": "matched by read names",
             "path_rule": "the run folder lies beside the analysis and carries the same sample name",
             "sample_name_chain": "matched by sample name only", "manual": "stated by a person"}


def _text(v) -> str:
    if v is None:
        return ""
    if isinstance(v, bool):
        return "yes" if v else "no"
    if isinstance(v, float):
        return f"{v:.6g}"
    if isinstance(v, list):
        return ",".join(_text(x.get("run") or x.get("fc") if isinstance(x, dict) else x) for x in v)
    return str(v)


def list_analyses(index: dict, *, columns: str | None, fmt: str) -> None:
    rows = list(index.get("analyses") or [])
    known = sorted({k for r in rows for k in r}) or list(DEFAULT_COLUMNS)
    cols = [c.strip() for c in columns.split(",") if c.strip()] if columns else list(DEFAULT_COLUMNS)
    bad = [c for c in cols if c not in known and c not in DEFAULT_COLUMNS]
    if bad:
        raise UsageError(f"no such column: {', '.join(bad)}", hint="columns: " + ", ".join(known))
    if fmt == "json":
        click.echo(json.dumps([{c: r.get(c) for c in cols} for r in rows], indent=1))
        return
    if fmt == "tsv":
        click.echo("\t".join(cols))
        for r in rows:
            click.echo("\t".join(_text(r.get(c)).replace("\t", " ").replace("\n", " ") for c in cols))
        return

    def show(c: str, v) -> str:
        if c == "date" and v:
            return str(v)[:10]
        if c == "verdict" and v:
            return VERDICT_WORD.get(v, str(v))
        if c == "assay":
            return ASSAY_WORD.get(v, str(v or ""))
        return _text(v)
    table = [[show(c, r.get(c)) for c in cols] for r in rows]
    widths = [max([len(c)] + [len(t[i]) for t in table]) for i, c in enumerate(cols)]
    widths[-1] = 0
    click.echo("  ".join(c.ljust(w) for c, w in zip(cols, widths)).rstrip())
    for t in table:
        click.echo("  ".join(x.ljust(w) for x, w in zip(t, widths)).rstrip())
    click.echo(f"{len(rows)} analys{'is' if len(rows) == 1 else 'es'}", err=True)


def resolve(store: store_mod.AnalysisStore, text: str) -> str:
    """An analysis id, or the label of exactly one analysis of the library."""
    if is_analysis_id(text) and store.has(text):
        return text
    cards = store.cards()
    hits = [aid for aid, c in cards.items() if text in (c["row"].get("label"), c["row"].get("group"), c["row"].get("slide"))]
    if len(hits) == 1:
        return hits[0]
    if len(hits) > 1:
        raise InputError(f"{text}: {len(hits)} analyses carry this name", hint="one of: " + ", ".join(sorted(hits)))
    raise InputError(f"{text}: no such analysis in this library", hint="flowcellar ls -L DIR --analyses lists them")


def show(library_path: str, text: str, *, what: str | None, section: str | None, as_json: bool) -> None:
    store = store_mod.AnalysisStore(library_path)
    aid = resolve(store, text)
    bundle = store.composed(aid)
    status = bundle.get("status") or {}
    ana = bundle["identity"]["analysis"]
    sec = (bundle.get("modules") or {}).get(MODULE_ID) or {}
    if section:
        if section != MODULE_ID:
            raise InputError(f"{aid} has no section '{section}'", hint="sections: " + MODULE_ID)
        if as_json:
            click.echo(json.dumps(sec, indent=1))
            return
        click.echo(f"{aid}  section {MODULE_ID}  state {sec.get('state', 'ok')}" + (f" ({sec.get('partial_reason')})" if sec.get("partial_reason") else ""))
        from_report = sec.get("value_sources") or {}
        for k, v in (sec.get("summary") or {}).items():
            click.echo(f"  {k:<40} {v}" + ("   (from report)" if (from_report.get(k) or {}).get("kind") in ("html_table", "embedded_json") else ""))
        return
    if what == "inputs":
        rows = [{"role": i.get("role"), "state": i.get("state"), "path": i.get("path"), "size": i.get("size"), "records": i.get("records")}
                for i in (bundle.get("provenance") or {}).get("inputs") or []]
        if as_json:
            click.echo(json.dumps(rows, indent=1))
            return
        for r in rows:
            size = f"{r['size'] / 1e6:,.1f} MB" if isinstance(r["size"], (int, float)) else ""
            click.echo(f"{(r['role'] or ''):<32} {(r['state'] or ''):<14} {size:>12}  {r['path']}")
        return
    if what == "status":
        if as_json:
            click.echo(json.dumps(status, indent=1))
            return
        for r in status.get("rules") or []:
            mark = {"triggered": (r.get("severity") or "").upper() or "TRIGGERED", "ok": "ok", "not_evaluated": "not assessed"}.get(r.get("outcome"), "")
            click.echo(f"{mark:<12} {(r.get('role') or ''):<10} {r.get('id')}")
            if r.get("message"):
                click.echo(f"{'':<12} {r['message']}")
        return
    if as_json:
        click.echo(json.dumps(bundle, indent=1))
        return
    click.echo(aid + "  " + " · ".join(str(x) for x in (ana.get("label"), ASSAY_WORD.get(ana.get("assay"), ""),
                                                       (PIPELINE_WORD.get(ana.get("pipeline"), "") + " " + (ana.get("pipeline_version") or "")).strip(),
                                                       ana.get("slide_serial_capture_area")) if x))
    click.echo(f"  verdict      {VERDICT_WORD.get(status.get('verdict'), 'Not assessed')}")
    if status.get("headline"):
        click.echo(f"  {status['headline']}")
    if status.get("detail"):
        click.echo(f"  {status['detail']}")
    for r in status.get("rules") or []:
        if r.get("outcome") == "triggered" and r.get("severity") in ("warn", "fail") and r.get("role") != "info":
            click.echo(f"  {r['severity']:<5} {r.get('message')}")
    links = ana.get("links") or []
    if not links:
        click.echo("  runs         none named")
    for ln in links:
        share = f"{ln['share'] * 100:.1f} % of the reads" if isinstance(ln.get("share"), (int, float)) else "share not known"
        key = ln.get("run_key") or "_".join(x for x in (ln.get("flow_cell_id"), ln.get("protocol_run_id_prefix")) if x)
        click.echo(f"  run          {key}  {share}  ({LINK_WORD.get(ln.get('method'), 'linked')}, {ln.get('confidence')} confidence)"
                   + ("" if ln.get("resolved") else "  run not in library"))
    for s in ana.get("sibling_analyses") or []:
        click.echo(f"  same reads   {s}")
    click.echo(f"  folder       {ana.get('path')}")
    click.echo(f"  bundle       {store.bundle_path(aid)}")
    click.echo("  more: --status (every check), --inputs (files), --section sc, --json")
