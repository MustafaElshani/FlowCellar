"""flowcellar analysis: one single-cell or spatial analysis folder -> a bundle and one offline HTML report.

The add-on (``flowcellar.sc``) is loaded by this command and by a scan that finds an analysis; nothing else
loads it.
"""
from __future__ import annotations

import logging
import os
import time

import click

from flowcellar.cli._options import library_option, resolve_library
from flowcellar.util import log as _log
from flowcellar.util.errors import ExitCode, FlowCellarError, InputError, UsageError

_LEVELS = {"warn": ("warn", "fail"), "fail": ("fail",)}


def _secs(s: float) -> str:
    s = int(round(s))
    return f"{s} s" if s < 60 else f"{s // 60} min {s % 60} s" if s < 3600 else f"{s // 3600} h {s % 3600 // 60} min"


def _size(n: int) -> str:
    return f"{n / 1e6:.1f} MB" if n >= 1e6 else f"{max(n // 1000, 1)} kB"


def _say(event: str, label: str, text: str, **fields) -> None:
    _log.event(f"analysis.{event}", f"  {label:<14} {text}" if label else text, **fields)


def _first_line(bundle: dict) -> str:
    from flowcellar.sc import ASSAY_WORD, PIPELINE_WORD
    ana = bundle["identity"]["analysis"]
    parts = [ASSAY_WORD.get(ana.get("assay"), "assay not known"),
             (PIPELINE_WORD.get(ana.get("pipeline"), ana.get("pipeline") or "") + " " + (ana.get("pipeline_version") or "")).strip()]
    if ana.get("slide_serial_capture_area"):
        parts.append("slide " + ana["slide_serial_capture_area"])
    return f"{ana.get('label') or bundle['bundle_id']}  " + " · ".join(p for p in parts if p)


def _links_line(links: list[dict], in_library: bool) -> str:
    if not links:
        return "no flow cell run could be named (state it with --run FLOWCELL)"
    words = {"read_tag": "file names stored on the reads", "bam_header_rg": "read group in the BAM header",
             "sibling_analysis": "the workflow analysis of the same reads", "read_id_sketch": "read names",
             "path_rule": "folder and sample name", "sample_name_chain": "sample name only", "manual": "stated"}
    parts = []
    for ln in links:
        share = f" {ln['share'] * 100:.1f} %" if isinstance(ln.get("share"), (int, float)) else ""
        key = ln.get("run_key") or "_".join(x for x in (ln.get("flow_cell_id"), ln.get("protocol_run_id_prefix")) if x)
        missing = " (run not in library)" if in_library and not ln.get("resolved") else ""
        parts.append(f"{key}{share}{missing}")
    how = sorted({words.get(ln.get("method"), ln.get("method") or "") for ln in links})
    conf = sorted({ln.get("confidence") or "" for ln in links})
    return ", ".join(parts) + f"   by {', '.join(how)} ({', '.join(conf)} confidence)"


@click.command("analysis", short_help="One single-cell or spatial analysis folder: a bundle and one offline HTML report.")
@click.argument("analysis_dir", type=click.Path(file_okay=False))
@click.option("-o", "--outdir", type=click.Path(file_okay=False), help="Where the bundle and the report are written.")
@click.option("--prefix", metavar="TEXT", help="File name prefix. [default: the analysis id]")
@click.option("--run", "runs", multiple=True, metavar="FLOWCELL[@RUNID8]",
              help="A flow cell run the reads came from, when the files do not say (repeatable).")
@click.option("--sibling", type=click.Path(file_okay=False), metavar="DIR",
              help="The analysis of the same reads by the other pipeline (a Space Ranger outs/ folder, or a wf-single-cell sample folder).")
@click.option("--no-sibling", is_flag=True, help="Do not read the analysis of the same reads by the other pipeline with it.")
@click.option("--label", metavar="TEXT", help="The name to show instead of the pipeline's sample id.")
@click.option("--no-read-pass", is_flag=True, help="Do not read the workflow's per-read table (minutes): read counts come from its report.")
@click.option("--map-um", type=click.Choice(["64", "32"]), default="64", show_default=True, help="Bin size of a Visium HD map, in µm.")
@click.option("--no-tissue-image", is_flag=True, help="Leave the grey tissue picture out of a Visium HD map.")
@click.option("--no-report", is_flag=True, help="Write the bundle only.")
@click.option("--no-bundle", is_flag=True, help="Write the report only.")
@click.option("--fail-on", type=click.Choice(["warn", "fail"]), help="Exit with code 3 when the verdict is at or above this level.")
@library_option()
@click.pass_context
def cmd(ctx, analysis_dir, outdir, prefix, runs, sibling, no_sibling, label, no_read_pass, map_um, no_tissue_image, no_report,
        no_bundle, fail_on, library):
    """Read one analysis folder ANALYSIS_DIR: a Space Ranger outs/ folder (or the folder that holds it), or the
    output folder of wf-single-cell (or one sample folder inside it). Writes <id>.fcb.json.gz and
    <id>.flowcellar.html into --outdir, and with -L adds the analysis to a library.

    \b
      flowcellar analysis spaceranger/outs -o qc/
      flowcellar analysis wf-single-cell/ -o qc/ --run TST00001
      flowcellar analysis spaceranger/outs -L ~/flowcellar-library

    The analysis has a verdict of its own. It never changes the verdict of a flow cell run.
    """
    from flowcellar.model import bundle as B
    from flowcellar.model.runkey import is_safe_id
    from flowcellar.model.status import VERDICT_WORD
    from flowcellar.sc import PIPELINE_SPACERANGER, PIPELINE_WF, detect, judge, link
    from flowcellar.sc import collect as collect_mod
    from flowcellar.sc import report as sc_report
    from flowcellar.sc import store as store_mod

    lib_path = resolve_library(ctx, library, required=False)
    if not os.path.isdir(analysis_dir):
        raise InputError(f"{analysis_dir}: no such folder", hint="ANALYSIS_DIR is a Space Ranger outs/ folder or a wf-single-cell output folder")
    if sibling and not os.path.isdir(sibling):
        raise InputError(f"--sibling {sibling}: no such folder")
    if not outdir and not lib_path:
        raise UsageError("nowhere to write: give -o DIR, or -L DIR to add the analysis to a library")
    if no_report and no_bundle and not lib_path:
        raise UsageError("--no-report and --no-bundle together leave nothing to write")
    if prefix and not is_safe_id(prefix):
        raise UsageError(f"--prefix {prefix!r} cannot be used in a file name", hint="letters, digits, dot, dash and underscore")
    stated = link.parse_run_args(runs)
    obj = ctx.obj or {}
    t_start = time.monotonic()

    claims = detect.claims_at(analysis_dir)
    if not claims:
        raise InputError(f"{analysis_dir}: not an analysis folder this version knows",
                         hint="a Space Ranger outs/ folder holds metrics_summary.csv, web_summary.html and spatial/ or binned_outputs/; "
                              "a wf-single-cell output folder holds wf-single-cell-report.html and one folder per sample")
    if prefix and len(claims) > 1:
        raise UsageError(f"{analysis_dir} holds {len(claims)} samples: --prefix names one file", hint="give one sample folder, or leave --prefix out")
    # the analysis of the same reads by the other pipeline, as far as the files say (or as stated)
    sib_claims: dict[str, detect.Claim] = {}
    if sibling:
        found = detect.claims_at(sibling)
        if len(found) != 1 or found[0].pipeline == claims[0].pipeline:
            raise InputError(f"--sibling {sibling}: not one analysis of the other pipeline",
                             hint="for a wf-single-cell folder give the Space Ranger outs/ folder, and the other way round")
        sib_claims[claims[0].path] = found[0]

    def options(claim, stated_sibling=None) -> collect_mod.Options:
        return collect_mod.Options(read_pass=not no_read_pass, map_um=int(map_um), tissue_image=not no_tissue_image,
                                   workers=int(obj.get("jobs") or 0), label=label,
                                   sibling=(stated_sibling.path if stated_sibling is not None and stated_sibling.pipeline == PIPELINE_SPACERANGER
                                            and claim.pipeline == PIPELINE_WF else None),
                                   progress=lambda m: _log.event("analysis.progress", f"  {'reading':<14} {m}", level=logging.DEBUG))

    lib = None
    refresher = None
    lease = None
    run_facts: list = []
    if lib_path:
        from flowcellar.library.lease import Lease
        from flowcellar.library.store import LibraryStore
        lib = LibraryStore.open(lib_path, config_file=obj.get("config"))
        lease = Lease(lib.layout.lease, wait=True, on_wait=lambda holder: _log.event(
            "lease.wait", "a scan of this library is running" + (f" ({holder.words()})" if holder else "") + "; waiting for it"))
        lease.acquire()
    worst = "pass"
    exit_gate = False
    try:
        if lib is not None:
            refresher = store_mod.Refresher(lib, say=lambda event, message, **kw: _log.event(event, message, level=kw.get("level", 20)))
            run_facts = store_mod.run_facts_of_library(lib)
        for claim in claims:
            t0 = time.monotonic()
            group = [claim]
            stated_sib = sib_claims.get(claim.path)
            if stated_sib is not None:
                group.append(stated_sib)
            elif not no_sibling:
                auto = _auto_sibling(claim, detect)
                if auto is not None:
                    group.append(auto)
            measured: list[dict] = []
            for c in group:
                what = "the analysis of the same reads" if c is not claim else "the analysis"
                _log.event("analysis.reading", f"  {'reading':<14} {what}: {c.path}", level=logging.DEBUG)
                opts = options(c, stated_sib if c is claim else None)
                if refresher is not None:
                    aid = collect_mod.quick_id(c)
                    m = refresher.collect(aid, [c], opts, force=False)
                    if m is None:
                        if c is claim:
                            raise FlowCellarError(f"{c.path}: {refresher.failed.get(aid, 'could not be read')}")
                        continue
                else:
                    try:
                        m = collect_mod.collect(c, opts)
                    except FlowCellarError:
                        if c is claim:
                            raise
                        continue
                measured.append(m)
            main_b = measured[0]
            aid = main_b["bundle_id"]
            if refresher is not None:
                if stated:
                    links = refresher.store.stated_links()
                    links[aid] = [list(x) for x in stated]
                    refresher.store.write_stated_links(links)
                if stated_sib is not None and len(measured) == 2:
                    wf_b = next((m for m in measured if m["identity"]["analysis"]["pipeline"] == PIPELINE_WF), None)
                    sr_b = next((m for m in measured if m["identity"]["analysis"]["pipeline"] == PIPELINE_SPACERANGER), None)
                    if wf_b is not None and sr_b is not None:
                        wf_b["modules"]["sc"]["params"].setdefault("sibling_hint", {})["stated"] = sr_b["identity"]["analysis"]["path"]
                refresher.finish({m["bundle_id"]: m for m in measured}, run_facts, profiles=lib.config.get("profiles"))
                composed = [refresher.store.composed(m["bundle_id"]) for m in measured]
            else:
                composed = _standalone(measured, stated, link, judge, store_mod, PIPELINE_WF)
            b = composed[0]
            ana = b["identity"]["analysis"]
            st = b.get("status") or {}
            _log.event("analysis.identity", _first_line(b), analysis=aid)
            sec = b["modules"]["sc"]
            for n in sec.get("notes") or []:
                if n.get("code") in ("h5py_missing", "read_pass_skipped", "read_pass_partial", "input_unreadable", "report_unreadable",
                                     "link_tags_not_read", "map_needs_positions"):
                    _say("note", "note", n["message"], level=logging.WARNING if n["level"] == "warning" else logging.INFO)
            _say("verdict", "verdict", f"{VERDICT_WORD.get(st.get('verdict'), 'Not assessed')} · {st.get('headline') or ''}", verdict=st.get("verdict"))
            for r in judge.reasons(st, limit=6):
                _say("reason", "", f"  {'':<14} {r['message']}")
            _say("runs", "runs", _links_line(ana.get("links") or [], lib is not None))
            if len(composed) > 1:
                other = composed[1]["identity"]["analysis"]
                ost = composed[1].get("status") or {}
                _say("sibling", "same reads", f"{other.get('pipeline')} {other.get('pipeline_version') or ''}: "
                                               f"{VERDICT_WORD.get(ost.get('verdict'), 'Not assessed')} · {ost.get('headline') or ''}")
            if st.get("verdict") in ("warn", "fail") and (worst == "pass" or st["verdict"] == "fail"):
                worst = st["verdict"]
            name = prefix or aid
            if outdir:
                os.makedirs(outdir, exist_ok=True)
                if not no_bundle:
                    for cb in composed:
                        path = os.path.join(outdir, B.bundle_filename(name if cb is b else cb["bundle_id"]))
                        written = B.write_bundle(path, cb, with_hash=False)
                        if cb is b:
                            _say("bundle", "bundle", f"{path}  ({_size(os.path.getsize(path))})", path=path)
                        del written
                if not no_report:
                    path = os.path.join(outdir, B.report_filename(name))
                    sc_report.write_report(composed, path)
                    _say("report", "report", f"{path}  ({_size(os.path.getsize(path))})", path=path)
            _say("done", "done", f"in {_secs(time.monotonic() - t0)}")
        if lib is not None:
            from flowcellar.library.maintain import rebuild_index_and_hub
            rebuild_index_and_hub(lib)
            _say("library", "library", f"{lib.path}  ({len(refresher.rows())} analys{'is' if len(refresher.rows()) == 1 else 'es'})")
    finally:
        if lease is not None:
            lease.release()
    del t_start
    if fail_on and worst in _LEVELS[fail_on]:
        exit_gate = True
    if exit_gate:
        ctx.exit(int(ExitCode.GATE))


def _auto_sibling(claim, detect):
    """The analysis of the same reads by the other pipeline, when the files of this one say where it is: the
    Space Ranger folder named in the workflow's parameters or lying inside its output folder; for a Space Ranger
    folder, the workflow output folder it lies in."""
    from flowcellar.sc import PIPELINE_WF, wfsc
    if claim.pipeline == PIPELINE_WF:
        below = detect.spaceranger_outs_below(claim.path) or detect.spaceranger_outs_below(claim.top or claim.path)
        if len(below) == 1:
            return detect.spaceranger_claim(below[0])
        if "params" in claim.files:
            try:
                p = wfsc.read_json(claim.files["params"])
            except (FlowCellarError, OSError):
                return None
            sbam, out_dir = p.get("spaceranger_bam"), os.path.normpath(str(p.get("out_dir") or "")).lstrip("./")
            top = (claim.top or "").rstrip("/")
            if isinstance(sbam, str) and sbam and out_dir and top.endswith(out_dir):
                launch = top[: len(top) - len(out_dir)].rstrip("/") or "/"
                cand = os.path.dirname(os.path.normpath(os.path.join(launch, sbam))) if not os.path.isabs(sbam) else os.path.dirname(sbam)
                return detect.spaceranger_claim(cand)
        return None
    # a workflow output folder beside the Space Ranger folder whose parameters name this outs/ folder
    holder = os.path.dirname(os.path.dirname(claim.path))
    try:
        beside = sorted(e.path for e in os.scandir(holder) if e.is_dir(follow_symlinks=False) and not e.name.startswith("."))
    except OSError:
        beside = []
    for d in beside[:40]:
        for c in detect.wf_claims(d):
            back = _auto_sibling(c, detect)
            if back is not None and os.path.normpath(back.path) == os.path.normpath(claim.path):
                return c
    # a Space Ranger folder inside a workflow output folder (the older layout)
    d = os.path.dirname(claim.path)
    for _ in range(4):
        d = os.path.dirname(d)
        if not d or d == os.path.dirname(d):
            break
        found = [c for c in detect.wf_claims(d)]
        if len(found) == 1:
            return found[0]
        if found:
            inside = [c for c in found if claim.path.startswith(c.path.rstrip("/") + "/")]
            return inside[0] if len(inside) == 1 else None
    return None


def _standalone(measured: list[dict], stated, link, judge, store_mod, pipeline_wf) -> list[dict]:
    """Without a library: the links the analysis carries by itself (or the stated ones), siblings named to each
    other, and a verdict that does not look anything up."""
    by_pipe = {m["identity"]["analysis"]["pipeline"]: m for m in measured}
    wf = by_pipe.get(pipeline_wf)
    out = []
    for m in sorted(measured, key=lambda x: x["identity"]["analysis"]["pipeline"] != pipeline_wf):
        ana = m["identity"]["analysis"]
        others = [o for o in measured if o is not m]
        if stated:
            ana["links"] = link.join(link.manual_links(stated, []), [])
        elif not ana.get("links") and wf is not None and m is not wf and wf["identity"]["analysis"].get("links"):
            ana["links"] = link.resolve(intrinsic=[], runs=[], sibling_links=wf["identity"]["analysis"]["links"])
        ana["sibling_analyses"] = [o["bundle_id"] for o in others]
        if m is wf and others:
            ana["group_key"] = others[0]["identity"]["analysis"].get("group_key") or ana.get("group_key")
        from flowcellar.model import bundle as B
        m["content_sha256"] = B.content_sha256(m)
        sib = {"barcodes_unique_molecules": ((others[0]["modules"]["sc"].get("summary") or {}).get("barcodes_unique_molecules"))} if others else None
        out.append(store_mod.compose(m, judge.judge(m, in_library=False, sibling=sib)))
    return sorted(out, key=lambda b: [m["bundle_id"] for m in measured].index(b["bundle_id"]))
