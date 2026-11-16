"""Which loading of its flow cell a run is.

A flow cell is used in loadings. The use index of a run is 1 + the number of earlier loadings of that flow cell.
An earlier run counts as a loading when it sequenced for 2 h or more, or produced 2 % or more of the fixed 24 h
warning level of its device and assay in estimated bases; a shorter earlier run is a restart of the next loading.
"""
from __future__ import annotations

from flowcellar.engine.limits import loading_floor_bases

DEFAULTS = {"restart_max_hours": 2.0, "restart_max_output_frac": 0.02}


def _settings(cfg: dict | None) -> tuple[float, float]:
    cfg = cfg or {}
    return (float(cfg.get("restart_max_hours", DEFAULTS["restart_max_hours"])),
            float(cfg.get("restart_max_output_frac", DEFAULTS["restart_max_output_frac"])))


def is_restart(f: dict, cfg: dict | None = None) -> bool:
    """A start too short and too small to be a loading."""
    max_h, frac = _settings(cfg)
    hours, bases = f.get("run_hours") or 0.0, f.get("est_bases") or 0
    return bool(hours < max_h and bases < loading_floor_bases(f.get("family"), f.get("agroup"), frac))


def annotate(facts: list[dict], cfg: dict | None = None) -> dict[str, dict]:
    """run_key -> use of its flow cell, for a set of runs in any order. Each entry: ``index``, ``of``, ``restart``,
    ``restart_of``, ``prior_runs``, ``prior_hours``, ``prior_estimated_bases``, ``runs`` (every known run of the
    flow cell in order of start)."""
    by_fc: dict[str, list[dict]] = {}
    for f in facts:
        if f.get("flow_cell_id"):
            by_fc.setdefault(f["flow_cell_id"], []).append(f)
    out: dict[str, dict] = {}
    for rs in by_fc.values():
        rs.sort(key=lambda r: (r.get("started") or "", r.get("run_key") or ""))
        loadings = 0
        uses = []
        for i, r in enumerate(rs):
            # a folder with raw data only tells nothing about its length: it is neither a loading nor a restart
            stub = not r.get("has_tier1", True)
            restart = (not stub) and is_restart(r, cfg)
            u = {"index": loadings + 1, "restart": restart, "prior_runs": i,
                 "prior_hours": round(sum((p.get("run_hours") or 0.0) for p in rs[:i]), 2),
                 "prior_estimated_bases": int(sum((p.get("est_bases") or 0) for p in rs[:i]))}
            uses.append(u)
            if not stub and not restart:
                loadings += 1
        n_load = max(u["index"] for u in uses)
        listing = [{"run_key": r["run_key"], "started": _z(r.get("started")), "hours": round(r["run_hours"], 2) if r.get("run_hours") is not None else None,
                    "estimated_bases": int(r["est_bases"]) if r.get("est_bases") is not None else None,
                    "pores_first_scan": r.get("mux0_pores"), "loading": u["index"], "restart": u["restart"]} for r, u in zip(rs, uses)]
        for k, (r, u) in enumerate(zip(rs, uses)):
            nxt = next((x for x, ux in zip(rs[k + 1:], uses[k + 1:]) if not ux["restart"] and ux["index"] == u["index"]), None)
            u["restart_of"] = nxt["run_key"] if (u["restart"] and nxt) else None
            u["of"] = n_load
            u["runs"] = listing
            out[r["run_key"]] = u
    return out


def _z(ts: str | None) -> str | None:
    return (ts[:19] + "Z") if ts and len(ts) >= 19 else None


def use_of(f: dict, history: list[dict], cfg: dict | None = None, declared_use: int | None = None) -> dict:
    """The use of one run's flow cell given the other runs the library knows, as ``relations.flow_cell_use`` plus
    ``restart``. Runs that started later are listed and never change the index."""
    same = [h for h in history if h.get("flow_cell_id") == f.get("flow_cell_id") and h.get("run_key") != f.get("run_key")]
    u = annotate(same + [f], cfg)[f["run_key"]] if f.get("flow_cell_id") else \
        {"index": 1, "of": None, "restart": False, "restart_of": None, "prior_runs": 0, "prior_hours": 0.0, "prior_estimated_bases": 0, "runs": []}
    source = "library" if u["prior_runs"] else "unknown"
    if declared_use is not None:
        u = dict(u, index=int(declared_use), of=max(u.get("of") or 1, int(declared_use)))
        source = "declared"
    u["source"] = source
    return u


def relation(u: dict) -> dict:
    """The stored form (``relations.flow_cell_use``)."""
    return {"index": u["index"], "of": u.get("of"), "source": u["source"], "prior_runs": u.get("prior_runs"),
            "prior_hours": u.get("prior_hours"), "prior_estimated_bases": u.get("prior_estimated_bases"),
            "restart_of": u.get("restart_of"), "runs": u.get("runs") or []}
