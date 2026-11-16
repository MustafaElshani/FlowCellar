"""From the small files and the cached partials of one run to its measured bundle.

``assemble`` reads the small files of the run folder (quick detail), merges the fragments of the partials it is
given (full detail when the summary partial is complete), adds the read-data sections, and returns the bundle
without ``status``, ``user`` and ``relations``. It opens no large input.
"""
from __future__ import annotations

import os
import time
from typing import Sequence

from flowcellar.contracts import PartialRef, ReadData, RunContext, RunIdentity
from flowcellar.engine import identity as ident
from flowcellar.model import axes as axes_mod
from flowcellar.model import bundle as B
from flowcellar.model import partial as P
from flowcellar.modules import dyn as mod_dyn
from flowcellar.modules import fc as mod_fc
from flowcellar.modules import len_est as mod_len
from flowcellar.modules import run as mod_run
from flowcellar.modules import yield_ as mod_yield
from flowcellar.modules.common import Inputs, note
from flowcellar.sources import minknow as mk
from flowcellar.sources import porescan, seqsummary
from flowcellar.util.errors import InputError
from flowcellar.util.timeutil import from_epoch, now_z

SECTION_ORDER = ("run", "fc", "yield", "dyn", "len.est", "yield.basecalled", "len", "qual", "aln", "cdna", "funnel", "funnel.samples")
_SMALL_ROLES = ("minknow.final_summary", "minknow.sample_sheet", "minknow.barcode_alignment", "minknow.report_json",
                "minknow.throughput", "minknow.pore_activity")
_TIER1_REQUIRED = ("minknow.final_summary", "minknow.report_json", "minknow.throughput", "minknow.pore_activity",
                   "minknow.sequencing_summary")


def _axis_as_the_schema_allows(axis: dict) -> dict:
    """An axis definition without the members the shipped bundle schema does not know for its type.

    The axis registry may carry ``note`` and ``closed_last``, which schema 0.1.1 does not list yet;
    a bundle must validate against the schema it names. Nothing is dropped once the schema lists them.
    """
    from flowcellar.model import schema as _schema
    for variant in (_schema.load_schema().get("$defs", {}).get("axis", {}).get("oneOf") or []):
        props = variant.get("properties") or {}
        if (props.get("type") or {}).get("const") == axis.get("type"):
            if variant.get("additionalProperties") is False:
                return {k: v for k, v in axis.items() if k in props}
            break
    return dict(axis)


def _timestamps(ctx: RunContext, inputs: Inputs) -> str:
    """``created`` and ``updated``: now; pinned with SOURCE_DATE_EPOCH or ``deterministic`` (then the time of the
    newest input, so that equal inputs give equal bytes)."""
    if os.environ.get("SOURCE_DATE_EPOCH", "").isdigit() or not ctx.deterministic:
        return now_z()
    times = [i["mtime"] for i in inputs.items if i.get("mtime")]
    return max(times) if times else from_epoch(0)


def _provenance(ctx: RunContext, run_dir: str, layout: str, inputs: Inputs, warnings: list[dict], partial_ids: list[str],
                timings: dict) -> dict:
    stamp = _timestamps(ctx, inputs)
    prov = {"generator": B.generator(), "created": stamp, "updated": stamp, "revision": 1,
            "library": {"id": ctx.library_id, "name": ctx.library_name or ctx.library_id} if ctx.library_id else None,
            "run_dir": {"path": run_dir, "root": ctx.root_id, "relpath": ctx.relpath, "layout": layout, "alternates": []},
            "inputs": inputs.items, "warnings": warnings, "profile": ctx.profile if ctx.profile in ("full", "standard", "compact") else "standard",
            "partials": sorted(set(partial_ids))}
    if not ctx.deterministic and not os.environ.get("SOURCE_DATE_EPOCH", "").isdigit():
        prov["timings"] = {k: round(v, 3) for k, v in timings.items()}
    return prov


def stub_bundle(run: RunIdentity, ctx: RunContext) -> dict:
    """The bundle of a folder with raw data only: identity from names, no section."""
    listing = mk.list_run_folder(run.run_dir)
    identity, warnings = ident.raw_only_identity(listing)
    b = B.new_bundle(identity)
    inputs = Inputs(ctx.root_id, ctx.relpath, run.run_dir)
    b["provenance"] = _provenance(ctx, run.run_dir, "raw_only", inputs, warnings, [], {})
    b["tiers"] = {"run_health": {"state": "absent", "detail": None, "missing_roles": list(_TIER1_REQUIRED)},
                  "read_level": {"state": "absent", "origin": None, "samples_expected": None, "samples_with_data": 0,
                                 "reads_expected": None, "reads_covered": 0, "missing_roles": [], "population": None}}
    return B.order_members(b)


#: the kinds of file MinKNOW writes, in words (a note never shows an internal name)
_FILE_KIND = {"minknow.final_summary": "final summary", "minknow.sample_sheet": "sample sheet", "minknow.report_json": "report",
              "minknow.throughput": "throughput table", "minknow.pore_activity": "pore activity table",
              "minknow.sequencing_summary": "sequencing summary", "minknow.pore_scan": "pore-scan table",
              "minknow.barcode_alignment": "barcode table", "minknow.output_hash": "list of output files",
              "minknow.report_html": "report page", "minknow.report_md": "report text"}


def _load_partials(refs: Sequence[PartialRef], run_key: str, warnings: list[dict], unopened: set[str] | None = None) -> dict[str, dict]:
    """role -> partial, for the two tier 1 roles; a partial that cannot be read is left out with a warning. A ref
    without a partial behind it (an input that could not be opened) is named in ``unopened``."""
    out: dict[str, dict] = {}
    for ref in refs:
        if ref.role not in (seqsummary.ROLE, porescan.ROLE):
            continue
        if not ref.partial_id or not ref.path:
            if unopened is not None and ref.state == "error":
                unopened.add(ref.role)
            continue
        try:
            p = P.read_partial(ref.path)
        except Exception as e:  # noqa: BLE001
            warnings.append(note("warning", "partial_unreadable", f"A cached result could not be read and was left out ({type(e).__name__})."))
            continue
        if run_key not in (p.get("routing") or {}).get("runs", {run_key: 0}):
            continue
        out[ref.role] = p
    return out


def assemble(run: RunIdentity, partials: Sequence[PartialRef] = (), read_data: ReadData | None = None,
             ctx: RunContext = RunContext()) -> dict:
    t_start = time.perf_counter()
    if run.kind == "run_raw_only":
        return stub_bundle(run, ctx)
    listing = mk.list_run_folder(run.run_dir)
    if listing.kind != "run":
        raise InputError(f"{run.run_dir}: no longer a run folder")
    files = listing.files
    inputs = Inputs(ctx.root_id, ctx.relpath, listing.run_dir)
    fs = mk.read_final_summary(files["minknow.final_summary"]) if "minknow.final_summary" in files else {}
    sheet = mk.read_sample_sheet(files["minknow.sample_sheet"]) if "minknow.sample_sheet" in files else []
    rep, report_problem = (mk.read_report(files["minknow.report_json"]) if "minknow.report_json" in files else ({}, None))
    identity, warnings = ident.build_identity(listing, rep, fs, sheet, assay_override=ctx.assay_override,
                                              stale_after_hours=ctx.stale_after_hours)
    run_state = identity["protocol_run"]["run_state"]
    run_key = identity["run_key"]
    in_fs = inputs.add("minknow.final_summary", files["minknow.final_summary"], ["run"]) if "minknow.final_summary" in files else None
    in_ss_sheet = inputs.add("minknow.sample_sheet", files["minknow.sample_sheet"], ["run"]) if "minknow.sample_sheet" in files else None
    in_ba = inputs.add("minknow.barcode_alignment", files["minknow.barcode_alignment"], ["run"]) if "minknow.barcode_alignment" in files else None
    in_rj = None
    if "minknow.report_json" in files:
        in_rj = inputs.add("minknow.report_json", files["minknow.report_json"], ["run", "yield", "fc"], state="error" if report_problem else "ok")
        if report_problem:
            warnings.append(note("error", "report_unreadable", "MinKNOW's report could not be used: " + report_problem + ".", input=in_rj))
    tp = mk.read_throughput(files.get("minknow.throughput"))
    in_tp = None
    if "minknow.throughput" in files:
        in_tp = inputs.add("minknow.throughput", files["minknow.throughput"], ["yield"],
                           state={"ok": "ok", "truncated": "truncated"}.get(tp.state, "error"), records=tp.rows)
    pa = mk.read_pore_activity(files.get("minknow.pore_activity"))
    in_pa = None
    if "minknow.pore_activity" in files:
        in_pa = inputs.add("minknow.pore_activity", files["minknow.pore_activity"], ["fc"],
                           state={"ok": "ok", "truncated": "truncated"}.get(pa.state, "error"), records=pa.rows)
    for role, others in sorted(listing.extra.items()):
        kind = _FILE_KIND.get(role) or role.split(".")[-1].replace("_", " ")
        warnings.append(note("info", "several_files_for_role",
                             f"The run folder holds {len(others) + 1} files of the kind \"{kind}\"; the largest is read "
                             f"({os.path.basename(files[role])})."))
    t_small = time.perf_counter() - t_start

    # ---------------------------------------------------------------- quick sections
    seq, _n_seq = mk.sequencing_acquisition(rep)
    ari = (seq or {}).get("acquisition_run_info") or {}
    cfg = ari.get("config_summary") if isinstance(ari.get("config_summary"), dict) else {}
    channels = identity["flow_cell"]["channel_count"]
    ratio, ratio_source = mod_yield.ratio_of(cfg, tp)
    run_hours = (int(tp.minute[-1]) / 60.0) if tp.usable else None
    has_ss, has_ps = "minknow.sequencing_summary" in files, "minknow.pore_scan" in files
    unopened: set[str] = set()
    loaded = _load_partials(partials, run_key, warnings, unopened)
    sec_fc, d_fc = mod_fc.build_quick(pa, rep, seq, channels=channels, run_hours=run_hours, in_rj=in_rj, in_pa=in_pa,
                                      waiting=has_ss or has_ps)
    run_minutes, duration_source = (int(tp.minute[-1]), "throughput") if tp.usable else (None, None)
    if run_hours is None and d_fc.get("pa_hours") is not None:
        run_hours, run_minutes, duration_source = d_fc["pa_hours"], pa.last_minute, "pore_activity"
    sec_yield, d_y = mod_yield.build_section(tp, ari, ratio=ratio, ratio_source=ratio_source, in_rj=in_rj if seq else None, in_tp=in_tp,
                                             p_first=d_fc.get("p_first"), act_h=d_fc.get("act_h") or {})
    sec_run, _d_run = mod_run.build_section(rep, ari, fs, run_hours=run_hours, planned_s=identity["protocol_run"].get("planned_duration_s"),
                                            report_problem=report_problem, inputs=[i for i in (in_fs, in_rj, in_ss_sheet, in_ba) if i],
                                            in_rj=in_rj, run_minutes=run_minutes, duration_source=duration_source)
    modules: dict[str, dict] = {}
    for name, sec in (("run", sec_run), ("fc", sec_fc), ("yield", sec_yield)):
        if sec is not None:
            modules[name] = sec

    # ---------------------------------------------------------------- the large inputs and their partials
    partial_ids: list[str] = []
    timings = {"small_files_s": t_small}
    detail, rh_state = "quick", "complete"
    ps_part, ss_part = loaded.get(porescan.ROLE), loaded.get(seqsummary.ROLE)
    in_ss = in_ps = None

    def large_input(role: str, path: str, part: dict | None, used: list[str]) -> str:
        if part is None and role in unopened:
            return inputs.add(role, path, [], state="error", fingerprint=None,
                              notes=[note("error", "not_opened", "Present, but it could not be opened, so it was not read.")])
        if part is None:
            return inputs.add(role, path, [], state="skipped", fingerprint=None,
                              notes=[note("info", "not_read", "Present, not read yet: the read count has not run for this file.")])
        cov = dict(part.get("coverage") or {})
        state = {"ok": "ok", "partially_read": "partially_read", "truncated": "truncated"}.get(part["state"], "error")
        extra = {"records": cov.get("records"), "coverage": cov}
        if part.get("errors"):
            extra["notes"] = [note("error", e.get("code", "read_error"), e.get("message", "")) for e in part["errors"]]
        if (part.get("facts") or {}).get("bad_rows"):
            extra["records_skipped"] = int(part["facts"]["bad_rows"])
        return inputs.add(role, path, used, state=state, fingerprint=part["input"].get("fingerprint"), **extra)

    def checked(part: dict | None, with_ratio: bool) -> dict | None:
        """The partial when it was made with the parameters the report gives now, else None (and a warning)."""
        if part is None:
            return None
        f = part.get("facts") or {}
        if f.get("channel_count") in (None, channels) and (not with_ratio or f.get("events_to_base_ratio") in (None, ratio)):
            return part
        warnings.append(note("warning", "partial_parameters_differ",
                             "A cached read count was made with another channel count or another number of bases per event "
                             "than the report gives now; it was left out and will be made again."))
        return None
    if has_ss:
        ss_part = checked(ss_part, True)
        in_ss = large_input("minknow.sequencing_summary", files["minknow.sequencing_summary"], ss_part,
                            ["dyn", "len.est", "fc"] if ss_part is not None and ss_part.get("fragments") else [])
    else:
        ss_part = None
    if has_ps:
        ps_part = checked(ps_part, False)
        in_ps = large_input("minknow.pore_scan", files["minknow.pore_scan"], ps_part,
                            ["fc"] if ps_part is not None and ps_part.get("fragments") else [])
    else:
        ps_part = None
    summary_complete = False
    if ss_part is not None and ss_part.get("fragments"):
        partial_ids.append(ss_part["partial_id"])
        cov = dict(ss_part["coverage"])
        state = ss_part["state"]
        expected_reads, expected_events = d_y.get("reads"), d_y.get("selected_events")
        got_reads, got_events = ss_part["facts"]["reads"], ss_part["facts"]["minknow_events"]
        check = {"reads_in_summary": got_reads, "reads_minknow": expected_reads, "events_in_summary": got_events,
                 "events_minknow": expected_events,
                 "equal": (expected_reads is not None and got_reads == expected_reads and got_events == expected_events) if state == "ok" else None}
        extra_notes: list[dict] = []
        if state == "ok" and expected_reads is not None and not check["equal"]:
            # the completeness test: a whole pass reproduces MinKNOW's own read count and event sum exactly
            if got_reads < expected_reads:
                state = "truncated"
                cov = {"kind": "truncated", "records": got_reads, "records_total": expected_reads,
                       "fraction": round(got_reads / expected_reads, 6) if expected_reads else None}
                extra_notes.append(note("warning", "summary_incomplete",
                                        f"The sequencing summary holds {got_reads:,} reads, MinKNOW counted {expected_reads:,}: "
                                        "the file is incomplete, and the numbers below are for the reads it holds.", input=in_ss))
            else:
                extra_notes.append(note("warning", "summary_differs_from_minknow",
                                        f"The sequencing summary holds {got_reads:,} reads and {got_events:,} events, "
                                        f"MinKNOW's throughput file {expected_reads:,} and {expected_events:,}.", input=in_ss))
        elif state == "ok":
            cov["records_total"] = got_reads
        elif expected_reads:
            # a part of the file: the share of MinKNOW's own read count says more than the share of bytes
            cov["records_total"] = expected_reads
            cov["fraction"] = round(min(got_reads / expected_reads, 1.0), 6)
        summary_complete = state == "ok"
        frags = {f["module"]: f for f in P.fragments_for(ss_part, run_key)}
        speed_min = (modules.get("fc", {}).get("params") or {}).get("target_speed_min_bps")
        hours = run_hours if run_hours is not None else 0.0
        if "dyn" in frags:
            modules["dyn"] = mod_dyn.build_section(frags["dyn"], state=state, coverage=cov, ratio=ratio, hours=hours, speed_min=speed_min,
                                                   in_ss=in_ss, completeness=check, extra_notes=extra_notes)
        if "len.est" in frags:
            modules["len.est"] = mod_len.build_section(frags["len.est"], state=state, coverage=cov, ratio=ratio, in_ss=in_ss,
                                                       extra_notes=extra_notes, exact=(ss_part.get("facts") or {}).get("exact"))
            short = modules["len.est"]["summary"].pop("dyn_short_read_est_frac", None)
            if "dyn" in modules:
                modules["dyn"]["summary"]["dyn_short_read_est_frac"] = short
        inp = inputs.get(in_ss)
        inp["coverage"] = cov
        if state != ss_part["state"]:
            inp["state"] = "truncated"
        if "fc" in modules:
            ps_frag = next(iter(P.fragments_for(ps_part, run_key, "fc")), None) if ps_part is not None else None
            if ps_frag is not None:
                partial_ids.append(ps_part["partial_id"])
            mod_fc.add_full(modules["fc"], summary_frag=frags.get("fc"), summary_complete=summary_complete, scan_frag=ps_frag, ratio=ratio,
                            channels=channels, layout=identity["flow_cell"].get("channel_layout"), hours=hours, in_ss=in_ss, in_ps=in_ps,
                            scan_missing=not has_ps)
        detail = "full" if summary_complete else "quick"
        rh_state = "complete" if summary_complete else "partial"
        timings["summary_pass_s"] = float((ss_part.get("timings") or {}).get("seconds") or 0.0)
        if ps_part is not None:
            timings["pore_scan_s"] = float((ps_part.get("timings") or {}).get("seconds") or 0.0)
        if not summary_complete:
            frac = cov.get("fraction")
            warnings.append(note("warning", "summary_partly_read",
                                 "The sequencing summary was not read whole"
                                 + (f" ({frac * 100:.1f} % of it)" if isinstance(frac, (int, float)) else "")
                                 + "; the sections on read dynamics and estimated length say so.", input=in_ss, value=frac))
    elif "fc" in modules and ps_part is not None and ps_part.get("fragments"):
        ps_frag = next(iter(P.fragments_for(ps_part, run_key, "fc")), None)
        if ps_frag is not None:
            partial_ids.append(ps_part["partial_id"])
            mod_fc.add_full(modules["fc"], summary_frag=None, summary_complete=False, scan_frag=ps_frag, ratio=ratio, channels=channels,
                            layout=identity["flow_cell"].get("channel_layout"), hours=run_hours or 0.0, in_ss=in_ss, in_ps=in_ps,
                            scan_missing=False)
            if ps_part is not None:
                timings["pore_scan_s"] = float((ps_part.get("timings") or {}).get("seconds") or 0.0)
    if ss_part is not None and not ss_part.get("fragments") and ss_part.get("errors"):
        # the pass ran and failed: an error section, never a crash and never silence
        msg = "; ".join(e.get("message", "") for e in ss_part["errors"])
        for name, version in (("dyn", mod_dyn.VERSION), ("len.est", mod_len.VERSION)):
            modules[name] = {"version": version, "tier": "run_health", "state": "error", "inputs": [in_ss] if in_ss else [],
                             "dims": ["end_reason"], "cells": [],
                             "notes": [note("error", "summary_unreadable", "The sequencing summary could not be read: " + msg, input=in_ss)]}
        partial_ids.append(ss_part["partial_id"])
        rh_state = "partial"                     # the read table is there and could not be read: not a complete picture
    if unopened:
        # the state stays what the small files give (quick detail): the read count is still owed, and says so
        for role in sorted(unopened):
            warnings.append(note("warning", "input_not_opened",
                                 f"{os.path.basename(files.get(role, role))} could not be opened; what it holds is not in this report."))

    # ---------------------------------------------------------------- tiers
    missing = [r for r in _TIER1_REQUIRED if r not in files]
    if "minknow.pore_scan" not in files:
        missing.append("minknow.pore_scan")
    if run_state == "running":
        rh_state, detail = "live", ("quick" if modules else None)
    elif not modules:
        rh_state, detail = "absent", None
    elif run_state == "interrupted" or any(modules[m]["state"] == "error" for m in ("run", "fc", "yield") if m in modules) \
            or any(m not in modules for m in ("run", "fc", "yield")):
        rh_state = "partial"
    tiers = {"run_health": {"state": rh_state, "detail": detail, "missing_roles": missing},
             "read_level": {"state": "absent", "origin": None, "samples_expected": None, "samples_with_data": 0,
                            "reads_expected": d_y.get("reads"), "reads_covered": 0, "missing_roles": [], "population": None}}
    samples: list[dict] = []

    # ---------------------------------------------------------------- read data
    if read_data is not None:
        for sid, sec in (read_data.modules or {}).items():
            modules[sid] = sec
        samples = list(read_data.samples or [])
        if read_data.read_level:
            rl = dict(read_data.read_level)
            rl.setdefault("reads_expected", d_y.get("reads"))
            tiers["read_level"] = rl
        inputs.items.extend(read_data.inputs or [])
        if read_data.basecalling:
            identity["basecalling"] = list(read_data.basecalling)
        partial_ids.extend(read_data.partial_ids or [])
        warnings.extend(read_data.warnings or [])
    ordered = {k: modules[k] for k in SECTION_ORDER if k in modules}
    ordered.update((k, v) for k, v in modules.items() if k not in ordered)
    timings["assemble_s"] = time.perf_counter() - t_start
    layout = ident.run_dir_layout(listing.run_dir, identity["protocol_run"].get("experiment_id"), identity["protocol_run"].get("sample_id"))
    b = B.new_bundle(identity)
    b["provenance"] = _provenance(ctx, listing.run_dir, layout, inputs, warnings, partial_ids, timings)
    b["tiers"] = tiers
    b["samples"] = samples
    b["modules"] = ordered
    used = axes_mod.used_axes(ordered)
    if any((c.get("channels") or (sec.get("run_level") or {}).get("channels")) for sec in ordered.values() for c in sec.get("cells", [])):
        used["channel"] = axes_mod.axis_def("channel")
    b["axes"] = {k: _axis_as_the_schema_allows(v) for k, v in sorted(used.items())}
    b["content_sha256"] = B.content_sha256(b)
    return B.order_members(b)
