"""The index: one columnar payload built from the cards, which is all Home, Runs and Projects ever load
.

``build_index`` is pure: it needs cards and nothing else, so a one-run report builds its index the same way a
library does. Numbers in ``library``, ``projects``, ``roots`` and ``activity`` are counted here from the run
columns; nothing in them is typed in by hand.
"""
from __future__ import annotations

from collections import Counter
from typing import Sequence

from flowcellar import SCHEMA_VERSION, __version__
from flowcellar import catalogue as catalogue_mod
from flowcellar.contracts import Card
from flowcellar.engine import api as engine_api
from flowcellar.library.projects import NO_PROJECT, project_id
from flowcellar.util.errors import NotBuiltError
from flowcellar.util.timeutil import now_z

INDEX_FORMAT = "flowcellar.index"
DICT_MAX = 64
VERDICT_KEYS = ("pass", "warn", "fail", "incomplete", "running", "unknown")

# the order of the columns in the payload (section 2.8); anything else a card carries follows in name order
COLUMN_ORDER = (
    "id", "fc", "use", "use_of", "use_source", "restart", "restart_of", "started", "hours", "planned_h", "run_state",
    "device_type", "position", "product", "kit", "assay", "assay_group", "chemistry", "minknow", "experiment", "run_name",
    "channels", "group", "folder", "roots", "project", "project_ids", "source_state", "verdict", "health", "data",
    "provisional", "pending", "checks_done", "cause", "modifiers", "headline", "detail", "reasons", "reason_levels",
    "reason_checks", "reason_tiers", "outcome_metric", "outcome_h", "outcome_value", "outcome_level", "outcome_n",
    "outcome_trust", "mark_flow_cell", "mark_loading", "mark_library", "reads", "gb",
    "reads_6h", "reads_12h", "reads_24h", "reads_48h", "reads_72h", "gb_6h", "gb_12h", "gb_24h", "gb_48h", "gb_72h",
    "pores", "pores_last", "active_peak", "pore_half_h", "duty_half_h", "occ_first", "occ_24h", "occ_run", "recovering",
    "kb_per_ch_h_48h", "mb_per_pore_48h", "speed", "n50_est", "normal_end_frac", "active_channels", "detail_t1",
    "spark_output", "spark_output_max_gb_per_3h", "spark_channels", "spark_channels_max", "history_hash", "limits",
    "read_data", "t2_population", "t2_samples", "t2_samples_with_data", "t2_reads", "t2_passed_gb", "t2_n50", "t2_q_med",
    "t2_ident_med", "t2_mapped_frac", "t2_usable", "t2_full_length", "coverage", "t2_files_waiting", "t2_conflicts", "bundle", "bundle_bytes",
    "bundle_gzip_bytes",
)
# identifiers and free text: never dictionary-encoded, whatever their number of distinct values
PLAIN_TEXT = frozenset(("id", "fc", "started", "experiment", "run_name", "folder", "headline", "detail", "history_hash",
                        "spark_output", "spark_channels", "bundle", "restart_of"))
LIST_COLUMNS = frozenset(("roots", "project_ids", "modifiers", "reasons", "reason_levels", "reason_checks", "reason_tiers"))
DEFAULTS = {"roots": [], "project_ids": [], "modifiers": [], "reasons": [], "reason_levels": [], "reason_checks": [],
            "reason_tiers": [], "source_state": "ok", "read_data": "absent", "t2_files_waiting": 0, "t2_conflicts": 0}

META_LAYOUT = "columnar: one array per field; low-cardinality strings as {dict, codes}"
META_SPARKLINES = ("32 points on a fixed 0 to 96 h axis (3 h per point), 6 bits per point, base64url alphabet "
                   "A-Z a-z 0-9 - _, scaled to the row's own maximum")
META_MARKS = "-1 low, 0 typical, 1 high, null unknown (no history of 5 comparable runs before the run)"


def _encode(name: str, values: list):
    """A string column with at most 64 distinct values becomes ``{"dict", "codes"}``; null is an entry of dict."""
    if name in PLAIN_TEXT or name in LIST_COLUMNS or not values:
        return values
    if not all(v is None or isinstance(v, str) for v in values) or all(v is None for v in values):
        return values
    words: dict = {}
    codes = [words.setdefault(v, len(words)) for v in values]
    if len(words) > DICT_MAX:
        return values
    return {"dict": list(words), "codes": codes}


def column(index: dict, name: str) -> list:
    """One run column of an index as a plain list, whichever way it is stored."""
    col = index["runs"]["columns"].get(name)
    if col is None:
        return [None] * index["runs"]["n"]
    if isinstance(col, dict):
        d = col["dict"]
        return [d[c] for c in col["codes"]]
    return col


def rows(index: dict, names: Sequence[str] | None = None) -> list[dict]:
    """The runs of an index as row objects (all columns, or the named ones)."""
    names = list(names) if names is not None else list(index["runs"]["columns"])
    cols = {n: column(index, n) for n in names}
    return [{n: cols[n][i] for n in names} for i in range(index["runs"]["n"])]


def _tally(verdicts) -> dict:
    c = Counter(v if v in VERDICT_KEYS else "unknown" for v in verdicts)
    return {k: c.get(k, 0) for k in VERDICT_KEYS}


def _bases(row: dict) -> float:
    gb = row.get("gb")
    return float(gb) * 1e9 if isinstance(gb, (int, float)) else 0.0


def _months(first: str, last: str) -> list[str]:
    y, m = int(first[:4]), int(first[5:7])
    ly, lm = int(last[:4]), int(last[5:7])
    out = []
    while (y, m) <= (ly, lm):
        out.append(f"{y:04d}-{m:02d}")
        y, m = (y + 1, 1) if m == 12 else (y, m + 1)
    return out


def _detail(row: dict) -> str:
    """full, quick or raw_only: how much of a run the library has."""
    if row.get("run_state") == "raw_only" or row.get("cause") == "raw_only" or row.get("detail_t1") in (None, "none"):
        return "raw_only"
    return "full" if row.get("detail_t1") == "full" else "quick"


def idle_progress() -> dict:
    return {"state": "idle", "started": None, "host": None, "time_box_minutes": None,
            "tier1": {"runs_waiting": 0, "bytes_waiting": 0, "bytes_done": 0, "mb_per_s": None, "eta_minutes": None},
            "read_data": {"files_total": 0, "files_read": 0, "files_waiting": 0, "bytes_waiting": 0, "eta_minutes": None},
            "message": None}


def build_index(cards: Sequence[Card], *, library: dict | None = None, roots: Sequence[dict] = (),
                samples: Sequence[dict] = (), notices: Sequence[dict] = (), progress: dict | None = None,
                baselines: dict | None = None, analyses: Sequence[dict] = ()) -> dict:
    """The index payload from cards. Pure: needs no library folder. ``analyses`` are the rows of the single-cell
    and spatial analyses (``flowcellar.sc.store``); each of their links is marked resolved when its run is a row
    of this index, so a link resolves by itself once its run is in the library."""
    cards = sorted(cards, key=lambda c: (str(c["columns"].get("started") or ""), str(c["columns"].get("id") or "")),
                   reverse=True)                                   # newest run first; runs without a start last
    n = len(cards)
    # a library without runs still has every column (as an empty array): a reader finds the names it asks for
    names = [c for c in COLUMN_ORDER if not cards or any(c in card["columns"] for card in cards) or c in DEFAULTS]
    names += sorted({k for card in cards for k in card["columns"]} - set(names))
    table = [{k: card["columns"].get(k, DEFAULTS.get(k)) for k in names} for card in cards]
    # loadings of a flow cell are library knowledge, not a judgement: a status file is frozen at the run's start
    # and cannot know of loadings that came later, so these four columns are counted again at every build
    try:
        uses = engine_api.flow_cell_columns(cards)
    except NotBuiltError:
        uses = {}
    for row in table:
        u = uses.get(row.get("id"))
        if u and row.get("use_source") == "library":
            row.update({k: u[k] for k in ("use", "use_of", "restart", "restart_of") if k in row})
    for row, card in zip(table, cards):
        for k in LIST_COLUMNS:
            if row.get(k) is None:
                row[k] = []
        row["_projects"] = card.get("projects") or [{"id": pid, "name": row.get("project") if i == 0 else pid, "roots": row["roots"]}
                                                    for i, pid in enumerate(row.get("project_ids") or [])]
    columns = {k: _encode(k, [row[k] for row in table]) for k in names}

    sample_rows = list(samples) or [dict(s) for card in cards for s in (card.get("samples") or [])]
    run_projects = {row.get("id"): row["_projects"] for row in table}

    # ---- projects ------------------------------------------------------------------------------------------
    by_project: dict[str, dict] = {}
    for row in table:
        entries = row["_projects"] or [{"id": NO_PROJECT, "name": None, "roots": row["roots"]}]
        for p in entries:
            g = by_project.setdefault(p["id"], {"id": p["id"], "name": p.get("name"), "rows": [], "roots": [], "samples": set()})
            g["rows"].append(row)
            for r in (p.get("roots") or []):
                if r not in g["roots"]:
                    g["roots"].append(r)
    for s in sample_rows:
        own = project_id(s["project"]) if s.get("project") else None
        ids = [own] if own in by_project else [p["id"] for p in run_projects.get(s.get("run"), [])][:1] or [NO_PROJECT]
        for pid in ids:
            if pid in by_project:
                by_project[pid]["samples"].add(s.get("sample_uid") or s.get("name") or s.get("key"))
    root_order = [r["id"] for r in roots]
    projects = []
    for g in by_project.values():
        starts = sorted(str(r["started"]) for r in g["rows"] if r.get("started"))
        projects.append({
            "id": g["id"], "name": g["name"], "runs": len(g["rows"]),
            "flow_cells": len({r.get("fc") for r in g["rows"] if r.get("fc")}), "samples": len(g["samples"]),
            "first": starts[0][:10] if starts else None, "last": starts[-1][:10] if starts else None,
            "verdicts": _tally(r.get("verdict") for r in g["rows"]),
            "estimated_bases": int(round(sum(_bases(r) for r in g["rows"]))),
            "reads": int(sum(r.get("reads") or 0 for r in g["rows"])),
            "roots": sorted(g["roots"], key=lambda r: root_order.index(r) if r in root_order else len(root_order)),
        })
    projects.sort(key=lambda p: (p["id"] == NO_PROJECT, _neg(p["last"]), str(p["name"] or "")))

    # ---- roots ---------------------------------------------------------------------------------------------
    roots_out = []
    for r in roots:
        mine = [row for row in table if r["id"] in row["roots"]]
        det = Counter(_detail(row) for row in mine)
        roots_out.append({
            "id": r["id"], "path": r.get("path"), "live": bool(r.get("live", False)),
            "project_depth": r.get("project_depth"), "project_pattern": r.get("project_pattern"),
            "state": r.get("state", "reachable"), "reason": r.get("reason"), "offline_since": r.get("offline_since"),
            "last_scanned": r.get("last_scanned"), "runs": len(mine),
            "detail": {"full": det.get("full", 0), "quick": det.get("quick", 0), "raw_only": det.get("raw_only", 0)},
            "read_data_files": dict(r.get("read_data_files") or {"found": 0, "read": 0}),
        })
    offline = {r["id"] for r in roots_out if r["state"] != "reachable"}
    if offline and "source_state" in columns:
        states = [("offline" if row["roots"] and all(x in offline for x in row["roots"]) else row.get("source_state"))
                  for row in table]
        for row, s in zip(table, states):
            row["source_state"] = s
        columns["source_state"] = _encode("source_state", states)

    # ---- activity ------------------------------------------------------------------------------------------
    started = sorted(str(r["started"]) for r in table if r.get("started"))
    activity = {"months": [], "runs": [], "estimated_bases": [], "reads": []}
    if started:
        months = _months(started[0], started[-1])
        runs_m, bases_m, reads_m = Counter(), Counter(), Counter()
        for r in table:
            if r.get("started"):
                m = str(r["started"])[:7]
                runs_m[m] += 1
                bases_m[m] += _bases(r)
                reads_m[m] += r.get("reads") or 0
        activity = {"months": months, "runs": [runs_m.get(m, 0) for m in months],
                    "estimated_bases": [int(round(bases_m.get(m, 0))) for m in months],
                    "reads": [int(reads_m.get(m, 0)) for m in months]}

    # ---- library -------------------------------------------------------------------------------------------
    lib = dict(library or {})
    det = Counter(_detail(row) for row in table)
    lib_out = {
        "id": lib.get("id"), "name": lib.get("name"), "path": lib.get("path"), "created": lib.get("created"),
        "scanned": lib.get("scanned"),
        "counts": {"runs": n, "flow_cells": len({r.get("fc") for r in table if r.get("fc")}),
                   "projects": sum(1 for p in projects if p["id"] != NO_PROJECT),
                   "samples": len({s.get("sample_uid") or s.get("name") or (s.get("run"), s.get("key")) for s in sample_rows}),
                   "by_verdict": _tally(r.get("verdict") for r in table),
                   "by_detail": {"full": det.get("full", 0), "quick": det.get("quick", 0), "none": det.get("raw_only", 0)},
                   "with_read_data": sum(1 for r in table if r.get("read_data") in ("partial", "complete"))},
        "estimated_bases": int(round(sum(_bases(r) for r in table))),
        "reads": int(sum(r.get("reads") or 0 for r in table)),
        "first_run": started[0] if started else None, "last_run": started[-1] if started else None,
    }

    if baselines is None:
        try:
            baselines = engine_api.baselines_today(cards)
        except NotBuiltError:
            baselines = {"groups": []}

    def values(name: str) -> list:
        return sorted({r[name] for r in table if r.get(name) is not None}, key=str)

    facets = {k: values(k) for k in ("verdict", "health", "data", "assay_group", "kit", "product", "device_type")}
    facets["project"] = sorted({p["name"] for p in projects if p["name"]}, key=str)
    facets["detail_t1"] = values("detail_t1")
    facets["root"] = sorted({x for r in table for x in r["roots"]})

    limits = next((r.get("limits") for r in table if r.get("limits")), None)
    return {
        "format": INDEX_FORMAT, "schema": SCHEMA_VERSION, "written": now_z(),
        "meta": {"library": lib.get("name"), "runs": n, "generator": {"name": "flowcellar", "version": __version__},
                 "catalogue": catalogue_mod.version(),
                 "limits": {"id": limits, "version": catalogue_mod.version(), "hash": limits},
                 "layout": META_LAYOUT, "sparklines": META_SPARKLINES, "marks": META_MARKS},
        "library": lib_out,
        "runs": {"n": n, "columns": columns},
        "samples": sample_rows,
        "projects": projects,
        "roots": roots_out,
        "activity": activity,
        "baselines": baselines,
        "facets": facets,
        "notices": list(notices),
        "progress": progress or idle_progress(),
        "analyses": _join_analyses(analyses, {row.get("id") for row in table}),
    }


def _join_analyses(analyses: Sequence[dict], run_ids: set) -> list[dict]:
    """The analysis rows with every link joined to the runs of this index."""
    out = []
    for a in analyses or ():
        a = dict(a)
        a["links"] = [dict(ln, resolved=ln.get("run") in run_ids) for ln in a.get("links") or []]
        out.append(a)
    return out


def analysis_rows(store) -> list[dict]:
    """The analysis rows a library holds on disk. The add-on is loaded only when the library has any."""
    import os
    try:
        if not any(n.startswith("an_") for n in os.listdir(store.layout.cards)):
            return []
    except OSError:
        return []
    from flowcellar.sc import store as sc_store
    return sc_store.index_rows(store.path)


def _neg(date: str | None):
    """Sort key that puts the latest date first and missing dates last."""
    return tuple(-ord(ch) for ch in date) if date else (1,)


def roots_rows(store) -> list[dict]:
    """The configured roots of a library with what the last scan recorded about each."""
    recorded = (store.read_aux(store.layout.roots_state, {}) or {}).get("roots", {})
    out = []
    for r in store.config["roots"]:
        s = recorded.get(r["id"], {})
        out.append({**r, "state": s.get("state", "reachable"), "reason": s.get("reason"),
                    "offline_since": s.get("offline_since"), "last_scanned": s.get("last_scanned"),
                    "read_data_files": s.get("read_data_files") or {"found": 0, "read": 0}})
    return out


def index_of_library(store, *, notices: Sequence[dict] | None = None, progress: dict | None = None,
                     scanned: str | None = None) -> dict:
    """The index of a library from what is on disk: its cards, its roots, the notices of the last scan."""
    cards = list(store.cards().values())
    saved = store.read_aux(store.layout.notices, {}) or {}
    return build_index(
        cards,
        library={"id": store.id, "name": store.name, "path": store.path, "created": store.created,
                 "scanned": scanned or saved.get("scanned")},
        roots=roots_rows(store),
        notices=list(notices) if notices is not None else list(saved.get("notices") or []),
        progress=progress or saved.get("progress"), analyses=analysis_rows(store))
