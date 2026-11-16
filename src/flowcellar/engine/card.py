"""The index row of one run and the facts the status engine reads.

``facts_of`` turns a measured bundle into the flat numbers the checks are made on; the same dict travels in the
card of the run (``card["facts"]``), so that a later run can be compared with this one without its bundle.
"""
from __future__ import annotations

import re

import numpy as np

from flowcellar.aggregates import build
from flowcellar.model import status as S

HORIZONS = (6, 12, 24, 48, 72)
SPARK_POINTS = 32
SPARK_BINS = 18                  # 600 s bins per sparkline point: 3 h


def _summary(b: dict, section: str) -> dict:
    return ((b.get("modules") or {}).get(section) or {}).get("summary") or {}


def _scalars(b: dict, section: str) -> dict:
    cells = ((b.get("modules") or {}).get(section) or {}).get("cells") or []
    return (cells[0].get("scalars") or {}) if cells else {}


def _series(b: dict, section: str, name: str) -> dict:
    cells = ((b.get("modules") or {}).get(section) or {}).get("cells") or []
    return ((cells[0].get("series") or {}).get(name) or {}).get("values") or {} if cells else {}


def minknow_parts(v: str | None) -> tuple[str | None, str | None]:
    """(major.minor.patch, major.minor) of a MinKNOW version text such as ``5.3.0-rc6-p2solo``."""
    m = re.match(r"^(\d+)\.(\d+)(?:\.(\d+))?", v or "")
    if not m:
        return None, None
    return ".".join(x for x in m.groups() if x is not None), f"{m.group(1)}.{m.group(2)}"


def _teaches(rh: dict, mods: dict, run_state: str | None) -> bool:
    """True when the numbers of the small files are whole, so that the run may serve as lab history.

    Complete run health, in quick or full detail; and also a run whose only gap is its read table (cut off on
    disk, or not read to the end yet): output over time, pores and occupancy do not come from that table, and the
    numbers that do (speed, estimated length) are kept out of the history on their own (``speed_med`` below).
    A run whose folder lacks a report, or that was interrupted, does not teach."""
    if rh.get("state") == "complete":
        return True
    if rh.get("state") != "partial" or run_state in ("interrupted", "running"):
        return False
    if any(m not in mods or mods[m].get("state") == "error" for m in ("run", "fc", "yield")):
        return False
    dyn = mods.get("dyn") or {}
    if dyn.get("state") == "error":
        return True                              # the read table could not be read at all: the same gap, only larger
    return dyn.get("state") == "partial" and dyn.get("partial_reason") == "partly_read"


def facts_of(bundle: dict) -> dict:
    """Everything the checks and the history need from one measured bundle."""
    ident = bundle.get("identity") or {}
    fc, pr = ident.get("flow_cell") or {}, ident.get("protocol_run") or {}
    dev, kit, assay = ident.get("device") or {}, ident.get("kit") or {}, ident.get("assay") or {}
    mods = bundle.get("modules") or {}
    tiers = bundle.get("tiers") or {}
    rh = tiers.get("run_health") or {}
    ys, fs, rs = _summary(bundle, "yield"), _summary(bundle, "fc"), _summary(bundle, "run")
    run_sc, fc_params = _scalars(bundle, "run"), (mods.get("fc") or {}).get("params") or {}
    started = pr.get("started")
    mk_full, mk_mm = minknow_parts((ident.get("software") or {}).get("minknow"))
    has_t1 = any(k in mods for k in ("run", "fc", "yield"))
    hours = rs.get("run_duration_h")
    minutes = ((mods.get("run") or {}).get("params") or {}).get("duration_minutes")
    if minutes is not None:
        hours = minutes / 60.0
    p_first = fs.get("fc_pores_first_scan")
    f: dict = {
        "run_key": ident.get("run_key"), "flow_cell_id": fc.get("id"), "started": started, "day": (started or "")[:10],
        "run_state": pr.get("run_state") or "unknown", "has_tier1": has_t1,
        "t1_complete": _teaches(rh, mods, pr.get("run_state")), "detail": rh.get("detail"),
        "device_type": dev.get("type"), "product_code": fc.get("product_code"), "kit": kit.get("sequencing_kit"),
        "family": fc.get("family"), "agroup": assay.get("group") or "unknown", "assay_type": assay.get("type") or "unknown",
        "chemistry": assay.get("chemistry"), "mk": mk_mm, "minknow": mk_full, "channel_count": fc.get("channel_count"),
        "experiment_id": pr.get("experiment_id"), "sample_id": pr.get("sample_id"), "position": dev.get("position"),
        "run_hours": hours, "planned_s": pr.get("planned_duration_s"),
        "reads": ys.get("yield_reads"), "est_bases": ys.get("yield_estimated_bases"),
        "mux0_pores": p_first, "mux_last_pores": fs.get("fc_pores_last_scan"),
        "msg0_num_pores": fs.get("fc_channels_with_pore_first_scan"), "chan_pore_frac": fs.get("fc_channels_with_pore_first_scan_frac"),
        "occ_30_90": fs.get("fc_occupancy_first_hour"), "occ_median_run": fs.get("fc_occupancy_median"),
        "occ_median_24h": fs.get("fc_occupancy_median_24h"), "rec_of_working": fs.get("fc_recovering_frac"),
        "rec_last_third_of_working": fs.get("fc_recovering_last_third_frac"), "adapter_frac": fs.get("fc_adapter_frac"),
        "seq_duty": fs.get("fc_sequencing_duty"), "temp_max_dev_c": fs.get("fc_temperature_max_dev_c"),
        "speed_min": fc_params.get("target_speed_min_bps"), "speed_max": fc_params.get("target_speed_max_bps"),
        "has_run_section": "run" in mods and mods["run"].get("state") in ("ok",),
        "state": rs.get("run_stop_cause"), "stop_reason": rs.get("run_stop_reason"),
        "msg_warn": run_sc.get("messages_warning"), "msg_error": run_sc.get("messages_error"),
        "msg_warn_disk": run_sc.get("messages_warning_computer"), "fallback_files": run_sc.get("files_in_fallback"),
        "write_failed": run_sc.get("raw_bytes_failed"), "write_backlog": run_sc.get("raw_bytes_backlog"),
        "pore_survival_24h": fs.get("fc_pore_survival_24h_frac"), "pore_half_life_h": fs.get("fc_pore_half_life_h"),
        "duty_half_life_h": fs.get("fc_duty_half_life_h"), "active_peak": fs.get("fc_active_channels_peak"),
        "kb_per_active_ch_hour_24h": ys.get("yield_kb_per_active_channel_hour_24h"),
        "kb_per_active_ch_hour_48h": ys.get("yield_kb_per_active_channel_hour_48h"),
        "summary_present": "minknow.sequencing_summary" not in (rh.get("missing_roles") or []),
        "exclude": False, "baseline_group": None,
    }
    for h in HORIZONS:
        f[f"reads_{h}h"] = ys.get(f"yield_reads_at_{h}h")
        f[f"bases_{h}h"] = ys.get(f"yield_estimated_bases_at_{h}h")
        f[f"mbpp_{h}h"] = (f[f"bases_{h}h"] / p_first / 1e6) if f[f"bases_{h}h"] and p_first else None
    dyn = mods.get("dyn")
    if dyn is not None:
        ds, ls = dyn.get("summary") or {}, _summary(bundle, "len.est")
        complete = dyn.get("state") == "ok"
        cov = (dyn.get("coverage") or {}).get("fraction")
        f["full"] = {"complete": complete, "error": dyn.get("state") == "error", "coverage_frac": 1.0 if complete else cov,
                     "speed_med": ds.get("dyn_speed_est_median"), "speed_drop": ds.get("dyn_speed_drop_frac"),
                     "sp_frac": ds.get("dyn_signal_positive_frac"), "unblock_frac": ds.get("dyn_unblock_frac"),
                     "est_n50": ls.get("len_est_n50"), "active_channels": fs.get("fc_active_channels") if complete else None,
                     "active_ch_frac": fs.get("fc_active_channels_frac") if complete else None, "reads_read": ds.get("dyn_reads")}
        f["speed_med"] = ds.get("dyn_speed_est_median") if complete else None
    else:
        f["full"] = None
        f["speed_med"] = None
    # sparklines: 32 points on a fixed 0 to 96 h axis (3 h per point)
    yv = _series(bundle, "yield", "yield_series").get("estimated_bases") or []
    sp_y = np.zeros(SPARK_POINTS)
    for i in range(min(len(yv), SPARK_POINTS * SPARK_BINS)):
        sp_y[i // SPARK_BINS] += yv[i]
    cs = _series(bundle, "fc", "fc_channel_state_series")
    act = np.zeros(0)
    if cs:
        n = len(next(iter(cs.values())))
        act = sum((np.asarray(cs.get(s, [0] * n), np.float64) for s in ("strand", "adapter", "pore")), np.zeros(n)) / 600.0
    sp_c, sp_n = np.zeros(SPARK_POINTS), np.zeros(SPARK_POINTS)
    for i in range(min(len(act), SPARK_POINTS * SPARK_BINS)):
        sp_c[i // SPARK_BINS] += act[i]
        sp_n[i // SPARK_BINS] += 1
    sp_c = np.where(sp_n > 0, sp_c / np.maximum(sp_n, 1), 0)
    f["spark_output_gb_per_3h"] = [round(float(v) / 1e9, 3) for v in sp_y]
    f["spark_working_channels"] = [round(float(v), 1) for v in sp_c]
    return f


def _spark(vals: list[float], hours: float | None) -> tuple[str, float]:
    n = min(SPARK_POINTS, int(np.ceil(hours / 3.0))) if hours else 0
    v = list(vals[:n]) + [0.0] * (SPARK_POINTS - n)
    text, mx = build.sparkline(v, points=SPARK_POINTS)
    return text, round(mx, 3)


def columns_of(f: dict, bundle: dict, status_file: dict | None) -> dict:
    """The ENGINE columns of the index row."""
    st = (status_file or {}).get("status") or {}
    rel = ((status_file or {}).get("relations") or {}).get("flow_cell_use") or {}
    rules = st.get("rules") or []
    reasons = S.reasons(st) if st else []
    outcome = st.get("_outcome") or outcome_of(st)
    ex = {e.get("aspect"): e for e in st.get("explanations") or [] if e.get("aspect") in ("flow_cell", "loading", "library")}

    def mark(aspect):
        return None if aspect not in ex else {"low": -1, "typical": 0, "high": 1}.get(ex[aspect].get("direction"))
    full = f.get("full") or {}
    hours = f.get("run_hours")
    sy, sy_max = _spark(f.get("spark_output_gb_per_3h") or [], hours)
    sc, sc_max = _spark(f.get("spark_working_channels") or [], hours)
    has_t1 = f.get("has_tier1")
    by_tier = st.get("by_tier") or {}
    restart = bool(rel.get("restart_of")) or any(r.get("run_key") == f.get("run_key") and r.get("restart") for r in rel.get("runs") or [])
    cols = {
        "id": f.get("run_key"), "fc": f.get("flow_cell_id"), "use": rel.get("index", 1), "use_of": rel.get("of"),
        "use_source": rel.get("source", "unknown"), "restart": restart, "restart_of": rel.get("restart_of"),
        "started": f.get("started"), "hours": round(hours, 2) if hours is not None and has_t1 else None,
        "planned_h": round(f["planned_s"] / 3600, 1) if f.get("planned_s") else None,
        "run_state": f.get("run_state") if f.get("run_state") in ("finished", "running", "interrupted", "raw_only") else "finished",
        "device_type": f.get("device_type"), "position": f.get("position"), "product": f.get("product_code"), "kit": f.get("kit"),
        "assay": f.get("assay_type"), "assay_group": f.get("agroup"), "chemistry": f.get("chemistry"), "minknow": f.get("minknow"),
        "experiment": f.get("experiment_id"), "run_name": f.get("sample_id"), "channels": f.get("channel_count"),
        "group": group_text(f),
        "verdict": st.get("verdict", "unknown"), "health": by_tier.get("run_health", "unknown"),
        "data": by_tier.get("read_level", "not_available"), "provisional": bool(st.get("provisional")),
        "pending": len(st.get("pending") or []), "checks_done": sum(1 for x in rules if x.get("outcome") != "not_evaluated"),
        "cause": st.get("cause"), "modifiers": [m["code"] for m in st.get("modifiers") or []],
        "headline": st.get("headline", ""), "detail": st.get("detail") or "",
        "reasons": [r["message"] for r in reasons], "reason_levels": [r["severity"] for r in reasons],
        "reason_checks": [r["check"] for r in reasons], "reason_tiers": [r["tier"] for r in reasons],
        "outcome_metric": outcome.get("metric"), "outcome_h": outcome.get("horizon_h"), "outcome_value": outcome.get("value"),
        "outcome_level": outcome.get("level"), "outcome_n": outcome.get("n"), "outcome_trust": outcome.get("trust"),
        "mark_flow_cell": mark("flow_cell"), "mark_loading": mark("loading"), "mark_library": mark("library"),
        "reads": f.get("reads"), "gb": round(f["est_bases"] / 1e9, 3) if f.get("est_bases") is not None else None,
        **{f"reads_{h}h": f.get(f"reads_{h}h") for h in HORIZONS},
        **{f"gb_{h}h": (round(f[f"bases_{h}h"] / 1e9, 3) if f.get(f"bases_{h}h") is not None else None) for h in HORIZONS},
        "pores": f.get("mux0_pores"), "pores_last": f.get("mux_last_pores"), "active_peak": f.get("active_peak"),
        "pore_half_h": f.get("pore_half_life_h"), "duty_half_h": f.get("duty_half_life_h"),
        "occ_first": f.get("occ_30_90"), "occ_24h": f.get("occ_median_24h"), "occ_run": f.get("occ_median_run"),
        "recovering": f.get("rec_of_working"),
        "kb_per_ch_h_48h": f.get("kb_per_active_ch_hour_48h"),
        "mb_per_pore_48h": round(f["mbpp_48h"], 2) if f.get("mbpp_48h") else None,
        "speed": full.get("speed_med"), "n50_est": full.get("est_n50"), "normal_end_frac": full.get("sp_frac"),
        "active_channels": full.get("active_channels"),
        "detail_t1": f.get("detail") or "none",
        "spark_output": sy, "spark_output_max_gb_per_3h": sy_max, "spark_channels": sc, "spark_channels_max": sc_max,
        "history_hash": (st.get("judged") or {}).get("history_hash"),
        "limits": limits_text(st.get("profile") or {}),
    }
    return cols


def card(bundle: dict, status_file: dict | None = None) -> dict:
    """ENGINE ``card``: ``{"columns": {...}, "facts": {...}}`` of one run."""
    f = facts_of(bundle)
    eng = (status_file or {}).get("engine") or {}
    if eng.get("group"):
        f["group"] = eng["group"]
    cols = columns_of(f, bundle, status_file)
    if eng.get("restart"):
        cols["restart"] = True
    return {"columns": cols, "facts": f}


def group_text(f: dict) -> str | None:
    from flowcellar.engine.baselines import group_key
    return f.get("group") or group_key(f, f.get("baseline_group"))


def limits_text(profile: dict) -> str:
    if not profile:
        return ""
    h = profile.get("sha256") or ""
    return f"{profile.get('id', 'unknown')}@{profile.get('version', '0')}" + (f":{h[:12]}" if h else "")


def outcome_of(st: dict) -> dict:
    """The judged outcome of a stored status: the key output check (or, on a re-used flow cell, output per starting
    pore), with the level that applied."""
    for r in st.get("rules") or []:
        if r.get("role") != "key" or r.get("outcome") == "not_evaluated":
            continue
        m = re.match(r"^yield_(reads|estimated_bases)_at_(\d+)h$", r.get("id") or "") or \
            re.match(r"^yield_(per_pore_mb)_(\d+)h$", r.get("id") or "")
        if not m:
            continue
        th = r.get("threshold") or {}
        bl = th.get("baseline") or {}
        level = bl.get("warn_level") if th.get("applied") == "baseline" else th.get("warn")
        return {"metric": r["id"], "horizon_h": int(m.group(2)), "value": r.get("value"), "level": level,
                "n": bl.get("n", th.get("peers") or 0), "trust": bl.get("trust", "none"),
                "unit": "reads" if m.group(1) == "reads" else "bases"}
    return {}
