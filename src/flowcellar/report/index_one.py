"""The one-run index a report embeds (mode A).

Two sources:

* ``slice_index``: one row cut out of a library's index (``flowcellar hub --reports``);
* ``index_for_bundle``: built from a composed bundle alone. It first asks LIBRARY
  (``library.api.build_index([card_for(bundle)])``), which owns the index; until that is built, or when it
  fails, ``plain_index`` writes the columns that can be read straight from the bundle. The viewer draws a run
  from its bundle and needs the row only for the page title and for the line it shows while the bundle decodes.
"""
from __future__ import annotations

from flowcellar import SCHEMA_VERSION, __version__
from flowcellar.model import status as status_model
from flowcellar.util.errors import FlowCellarError
from flowcellar.util.timeutil import now_z

FORMAT = "flowcellar.index"
_VERDICT_KEYS = ("pass", "warn", "fail", "incomplete", "running", "unknown")


def _summaries(bundle: dict) -> dict:
    out: dict = {}
    for m in (bundle.get("modules") or {}).values():
        out.update(m.get("summary") or {})
    return out


def plain_columns(bundle: dict) -> dict:
    """The index row of one run, as far as the bundle alone tells it."""
    ident, st, tiers = bundle.get("identity") or {}, bundle.get("status") or {}, bundle.get("tiers") or {}
    fc, pr, dev = ident.get("flow_cell") or {}, ident.get("protocol_run") or {}, ident.get("device") or {}
    kit, assay, rel = ident.get("kit") or {}, ident.get("assay") or {}, bundle.get("relations") or {}
    use, by, s = rel.get("flow_cell_use") or {}, st.get("by_tier") or {}, _summaries(bundle)
    rh, rl = tiers.get("run_health") or {}, tiers.get("read_level") or {}
    rs = status_model.reasons(st)
    projects = rel.get("projects") or []
    gb = s.get("yield_estimated_bases")
    return {
        "id": bundle.get("bundle_id") or ident.get("run_key"), "fc": fc.get("id"),
        "use": use.get("index") or 1, "use_of": use.get("of") or 1, "use_source": use.get("source") or "unknown",
        "restart": bool(use.get("restart_of")), "restart_of": use.get("restart_of"),
        "started": pr.get("started"), "hours": s.get("run_duration_h"),
        "planned_h": (pr.get("planned_duration_s") or 0) / 3600 or None,
        "run_state": pr.get("run_state") or "finished",
        "device_type": dev.get("type"), "position": dev.get("position"), "product": fc.get("product_code"),
        "kit": kit.get("sequencing_kit"), "assay": assay.get("type"), "assay_group": assay.get("group") or "unknown",
        "chemistry": assay.get("chemistry"), "minknow": (ident.get("software") or {}).get("minknow"),
        "experiment": pr.get("experiment_id"), "run_name": pr.get("sample_id"), "channels": fc.get("channel_count"),
        "group": " | ".join(str(x) for x in (dev.get("type"), fc.get("product_code"), kit.get("sequencing_kit"))),
        "folder": None, "roots": sorted({r for p in projects for r in (p.get("roots") or [])}),
        "project": projects[0]["name"] if projects else None, "project_ids": [p["id"] for p in projects],
        "source_state": "ok",
        "verdict": st.get("verdict") or "unknown", "health": by.get("run_health") or "not_available",
        "data": by.get("read_level") or "not_available", "provisional": bool(st.get("provisional")),
        "pending": len(st.get("pending") or []), "checks_done": len(st.get("rules") or []), "cause": st.get("cause"),
        "modifiers": [m.get("code") for m in st.get("modifiers") or []],
        "headline": st.get("headline") or "", "detail": st.get("detail") or "",
        "reasons": [r["message"] for r in rs], "reason_levels": [r["severity"] for r in rs],
        "reason_checks": [r["check"] for r in rs], "reason_tiers": [r["tier"] for r in rs],
        "reads": s.get("yield_reads"), "gb": gb / 1e9 if gb is not None else None,
        "pores": s.get("fc_pores_first_scan"), "pores_last": s.get("fc_pores_last_scan"),
        "occ_first": s.get("fc_occupancy_first_hour"),
        "detail_t1": rh.get("detail") if rh.get("state") not in (None, "absent") else None,
        "read_data": rl.get("state") or "absent", "t2_population": (rl.get("population") or {}).get("stage"),
        "t2_samples": rl.get("samples_expected"), "t2_samples_with_data": rl.get("samples_with_data"),
        "t2_reads": rl.get("reads_covered"),
        "t2_passed_gb": s["yield_basecalled_bases"] / 1e9 if s.get("yield_basecalled_bases") is not None else None,
        "t2_n50": s.get("len_n50"), "t2_q_med": s.get("qual_median"), "t2_usable": s.get("funnel_usable_frac"),
        "t2_files_waiting": 0, "spark_output": "", "spark_output_max_gb_per_3h": 0,
    }


def plain_index(bundle: dict) -> dict:
    """A one-run index in the shape of the library index, from the bundle alone."""
    cols = plain_columns(bundle)
    prov_lib = (bundle.get("provenance") or {}).get("library") or {}
    samples = [{"sample_uid": s.get("sample_uid"), "name": s.get("name") or s.get("alias") or s.get("key"), "alias": s.get("alias"),
                "barcode": s.get("barcode"), "project": s.get("project"), "run": cols["id"], "key": s.get("key"),
                "reads": None, "bases": None, "n50": None, "q_med": None}
               for s in bundle.get("samples") or [] if s.get("key") != "unclassified"]
    verdicts = {k: 0 for k in _VERDICT_KEYS}
    verdicts[cols["verdict"] if cols["verdict"] in verdicts else "unknown"] = 1
    started = cols["started"]
    return {
        "format": FORMAT, "schema": SCHEMA_VERSION, "written": now_z(),
        "meta": {"library": prov_lib.get("name"), "runs": 1, "generator": {"name": "flowcellar", "version": __version__}},
        "library": {"id": prov_lib.get("id"), "name": prov_lib.get("name"), "path": None, "created": None, "scanned": None,
                    "counts": {"runs": 1, "flow_cells": 1, "projects": len(cols["project_ids"]), "samples": len(samples),
                               "by_verdict": verdicts, "by_detail": {}, "with_read_data": int(cols["read_data"] != "absent")},
                    "estimated_bases": None, "reads": cols["reads"], "first_run": started, "last_run": started},
        "runs": {"n": 1, "columns": {k: [v] for k, v in cols.items()}},
        "samples": samples, "projects": [], "roots": [], "activity": {"months": [], "runs": [], "estimated_bases": [], "reads": []},
        "baselines": {"groups": []}, "facets": {}, "notices": [],
        "progress": {"state": "idle"},
    }


def index_for_bundle(bundle: dict) -> dict:
    """The index of a report made from one composed bundle: LIBRARY's when it can make one, else the plain one."""
    try:
        from flowcellar.library import api as library_api

        status_file = None
        if bundle.get("status") is not None:
            status_file = status_model.new_status_file(bundle.get("bundle_id"), bundle.get("content_sha256"),
                                                       bundle["status"], bundle.get("relations") or {})
        index = library_api.build_index([library_api.card_for(bundle, status_file)])
        if isinstance(index, dict) and (index.get("runs") or {}).get("n") == 1:
            return index
    except FlowCellarError:
        pass                                    # not built yet, or it cannot make a card of this bundle
    return plain_index(bundle)


def _value(col, i):
    if isinstance(col, dict) and "dict" in col:
        return col["dict"][col["codes"][i]]
    return col[i]


def slice_index(index: dict, run_key: str) -> dict | None:
    """One row of a library index as an index of its own: what a report of that run embeds. Names of things that
    are not in the file (other runs, other projects, the folders that were scanned) are left out."""
    cols = (index.get("runs") or {}).get("columns") or {}
    n = (index.get("runs") or {}).get("n") or 0
    ids = cols.get("id")
    i = next((k for k in range(n) if _value(ids, k) == run_key), None) if ids is not None else None
    if i is None:
        return None
    row = {k: _value(c, i) for k, c in cols.items()}
    out = {k: v for k, v in index.items() if k not in ("runs", "samples", "projects", "roots", "activity", "baselines", "notices", "progress", "facets", "library")}
    lib = index.get("library") or {}
    out["library"] = {"id": lib.get("id"), "name": lib.get("name"), "path": None, "created": lib.get("created"), "scanned": lib.get("scanned")}
    out["runs"] = {"n": 1, "columns": {k: [v] for k, v in row.items()}}
    out["samples"] = [s for s in index.get("samples") or [] if s.get("run") == run_key]
    out["projects"] = [p for p in index.get("projects") or [] if p.get("id") in (row.get("project_ids") or [])]
    out["roots"] = []
    out["activity"] = {"months": [], "runs": [], "estimated_bases": [], "reads": []}
    out["baselines"] = {"groups": [g for g in (index.get("baselines") or {}).get("groups") or [] if g.get("group") == row.get("group")]}
    out["facets"] = {}
    out["notices"] = []
    out["progress"] = {"state": "idle"}
    return out
