"""History: comparison groups, the robust score, the trust steps and the spread guard.

A run is placed against the last 40 eligible runs of its comparison group that started before it: first use of
the flow cell, at least 24 h long, tier 1 complete, not on the same day. ``z = (x - median) / max(1.4826 x MAD,
floor x |median|)``. Fewer than 5 such runs: no history. 5 to 7: provisional (can warn, never fail). 8 or more:
trusted. A window whose spread (1.4826 x MAD / median) exceeds 0.45 is provisional whatever its size.
"""
from __future__ import annotations

import math
import statistics

MADK = 1.4826
DEFAULTS = {"window": 40, "min_runs": 5, "trusted_runs": 8, "spread_guard": 0.45}


def settings(cfg: dict | None) -> dict:
    out = dict(DEFAULTS)
    for k in DEFAULTS:
        if isinstance((cfg or {}).get(k), (int, float)):
            out[k] = (cfg or {})[k]
    return out


def sig(x, n: int = 4):
    """Round to ``n`` significant digits."""
    if x is None or x == 0 or isinstance(x, bool):
        return x
    if not math.isfinite(x):
        return None
    return round(x, n - 1 - int(math.floor(math.log10(abs(x)))))


def group_label(f: dict, groups, declared: str | None = None) -> str | None:
    """The optional label that splits a comparison group: stated for the run, else the first matching dated rule of
    ``[[baselines.groups]]`` (``match``: kit, product, experiment, started_after, started_before)."""
    if declared:
        return str(declared)
    if f.get("baseline_group"):
        return str(f["baseline_group"])
    started = f.get("started") or ""
    for g_ in groups or ():
        m = g_.get("match") if isinstance(g_, dict) else None
        if not isinstance(m, dict) or not g_.get("label"):
            continue
        if "kit" in m and m["kit"] != f.get("kit"):
            continue
        if "product" in m and m["product"] != f.get("product_code"):
            continue
        if "experiment" in m and m["experiment"] != f.get("experiment_id"):
            continue
        if "started_after" in m and not (started and started[:len(str(m["started_after"]))] >= str(m["started_after"])):
            continue
        if "started_before" in m and not (started and started[:len(str(m["started_before"]))] < str(m["started_before"])):
            continue
        return str(g_["label"])
    return None


def group_key(f: dict, label: str | None = None) -> str | None:
    if not f.get("device_type") and not f.get("product_code") and not f.get("kit"):
        return None
    base = f"{f.get('device_type')} | {f.get('product_code')} | {f.get('kit')}"
    return f"{base} | {label}" if label else base


def robust(x, ref: list, floor: float, floor_abs: float = 0.0) -> dict | None:
    ref = [v for v in ref if v is not None and not (isinstance(v, float) and math.isnan(v))]
    if not ref or x is None:
        return None
    med = statistics.median(ref)
    mad = statistics.median([abs(v - med) for v in ref]) * MADK
    sc = max(mad, floor * abs(med), floor_abs)
    if sc <= 0:
        return None
    return {"z": (x - med) / sc, "n": len(ref), "median": med, "scale": sc, "rel": (mad / abs(med)) if med else float("nan"),
            "rank_below": sum(1 for v in ref if v < x)}


class History:
    """The runs a run may be compared with. ``runs``: facts of the other runs, each with ``group`` and
    ``use_index`` set; they are filtered to those that started before the run in question."""

    def __init__(self, runs: list[dict], cfg: dict | None = None):
        self.runs = sorted(runs, key=lambda r: (r.get("started") or "", r.get("run_key") or ""))
        self.cfg = settings(cfg)

    def peers(self, r: dict, col: str, extra=None) -> list[dict]:
        out = []
        started, day = r.get("started") or "", r.get("day")
        for p in self.runs:
            ps = p.get("started") or ""
            if p.get("run_key") == r.get("run_key") or not ps or not started or not (ps < started) or p.get("day") == day:
                continue
            if p.get("group") is None or p.get("group") != r.get("group") or p.get(col) is None:
                continue
            if p.get("use_index") != 1 or (p.get("run_hours") or 0) < 24 or p.get("exclude") or not p.get("t1_complete", True):
                continue
            if extra and not extra(p):
                continue
            out.append(p)
        return out[-int(self.cfg["window"]):]

    def base(self, r: dict, col: str, floor: float, extra=None, floor_abs: float = 0.0):
        """(placement or None, trust, number of peers)."""
        ps = self.peers(r, col, extra)
        n = len(ps)
        trust = "trusted" if n >= self.cfg["trusted_runs"] else "provisional" if n >= self.cfg["min_runs"] else "none"
        b = robust(r.get(col), [p[col] for p in ps], floor, floor_abs) if trust != "none" else None
        if b:
            b["guard"] = bool(b["rel"] > self.cfg["spread_guard"])
            if b["guard"] and trust == "trusted":
                trust = "provisional"
            b["trust"] = trust
            b["peers"] = [(p["run_key"], p[col]) for p in ps]
        return b, trust, n


def block(r: dict, b: dict, warn_z: float | None = None, fail_z: float | None = None, horizon: int | None = None,
          label: str | None = None) -> dict:
    """The stored ``baseline`` block of a rule or an explanation."""
    med, sc = b["median"], b["scale"]
    return {"group": r.get("group") or "", "group_label": label, "as_of": (r.get("started") or "")[:19] + "Z", "trust": b["trust"],
            "n": b["n"], "median": sig(med, 6), "scale": sig(sc, 6), "z": round(b["z"], 2), "rank_below": b["rank_below"],
            "band": [sig(med - 1.5 * sc, 6), sig(med + 1.5 * sc, 6)], "warn_z": warn_z, "fail_z": fail_z, "horizon_h": horizon,
            "spread": round(b["rel"], 3) if b["rel"] == b["rel"] else None, "spread_guard": b["guard"],
            "warn_level": sig(med + warn_z * sc, 6) if warn_z is not None else None,
            "fail_level": sig(med + fail_z * sc, 6) if fail_z is not None else None}
