"""Cards: the index row of one run, kept as one small file per run (``cards/<run_key>.card.json``).

A card is ``{"run_key", "columns", "samples", "facts", "projects"}``. ``columns`` holds the engine's columns, the
read-data columns and the library's own; ``facts`` is what the status engine needs
from this run when it is history for a later one; ``projects`` are the ``relations.projects`` entries, kept
so that the index can name every project of a run without opening its status file.
"""
from __future__ import annotations

import os

from flowcellar.contracts import Card
from flowcellar.engine import api as engine_api
from flowcellar.model import bundle as bundle_model
from flowcellar.readdata import api as readdata_api
from flowcellar.util import atomic

LIBRARY_COLUMNS = ("folder", "roots", "project", "project_ids", "source_state", "bundle", "bundle_bytes",
                   "bundle_gzip_bytes")


def library_columns(bundle: dict, status_file: dict | None, extra: dict | None = None) -> dict:
    rel = (status_file or {}).get("relations") or {}
    paths = rel.get("paths") or []
    projects = rel.get("projects") or []
    run_dir = ((bundle.get("provenance") or {}).get("run_dir") or {})
    where = (paths[0].get("relpath") if paths else None) or run_dir.get("relpath") or run_dir.get("path") or ""
    roots: list[str] = []
    for p in paths:
        if p.get("root") and p["root"] not in roots:
            roots.append(p["root"])
    cols = {
        "folder": os.path.basename(where.rstrip("/")) or None,
        "roots": roots,
        "project": projects[0]["name"] if projects else None,
        "project_ids": [p["id"] for p in projects],
        "source_state": "ok",
        "bundle": bundle_model.bundle_filename(bundle["bundle_id"]),
    }
    extra = dict(extra or {})
    if "bundle_bytes" not in extra or "bundle_gzip_bytes" not in extra:
        measured, _ = bundle_model.split(bundle)
        raw = atomic.dumps(bundle_model.order_members(measured)).encode("utf-8")
        extra.setdefault("bundle_bytes", len(raw))
        extra.setdefault("bundle_gzip_bytes", len(atomic.gzip_bytes(raw)))
    cols.update(extra)
    return cols


def card_for(bundle: dict, status_file: dict | None = None, *, extra: dict | None = None) -> Card:
    """The complete index row of one run: ENGINE's columns, READ DATA's columns and the library's own."""
    base = engine_api.card(bundle, status_file)
    columns = dict(base.get("columns") or {})
    columns.update(readdata_api.card_columns(bundle) or {})
    columns.update(library_columns(bundle, status_file, extra))
    columns.setdefault("id", bundle["bundle_id"])
    return {"run_key": bundle["bundle_id"], "columns": columns, "samples": list(readdata_api.sample_rows(bundle) or []),
            "facts": base.get("facts") or {},
            "projects": list(((status_file or {}).get("relations") or {}).get("projects") or [])}


def judge_settings(config: dict) -> dict:
    """What ``engine.api.judge`` takes from ``flowcellar.toml``: limit overrides per assay group, the four
    ``[baselines]`` settings (handed over under the key ``_baselines``, which the status engine reads), the
    dated comparison groups and the flow cell re-use settings."""
    b = config.get("baselines") or {}
    profile = dict(config.get("profiles") or {})
    profile["_baselines"] = {k: b[k] for k in ("window", "min_runs", "trusted_runs", "spread_guard") if b.get(k) is not None}
    return {"profile": profile, "baseline_groups": tuple(b.get("groups") or ()), "flow_cell": config.get("flow_cell")}
