"""Limits of the tier 1 checks: the catalogue's defaults, and what a library's configuration changes.

Order of precedence: ``[profiles.<assay group>.rules.<metric id>]`` of the configuration, then the catalogue.
"""
from __future__ import annotations

from functools import lru_cache

from flowcellar import catalogue
from flowcellar.util.fingerprint import hash_obj

HORIZONS = (72, 48, 24, 12, 6)
#: floors of the robust scale, as a share of the median
FLOOR = {"yield": 0.15, "pores": 0.08, "occ": 0.04, "speed": 0.04, "half": 0.20, "per_pore": 0.12, "kb": 0.12, "other": 0.10}
BAD_STOPS = ("STOPPED_NO_DISK_SPACE", "STOPPED_DEVICE_ERROR", "STOPPED_BAD_TEMPERATURE", "STOPPED_INTERNAL_ERROR",
             "STOPPED_DEVICE_STOPPED_ACQUISITION", "STOPPED_SHUTDOWN")
SPEED_TOLERANCE = 0.05
SPEED_FAIL_BELOW = 0.15
#: checks that wait for the read count while tier 1 is in quick detail
PENDING = ("dyn_speed_est_median", "dyn_speed_drop_frac", "dyn_signal_positive_frac", "dyn_unblock_frac", "fc_active_channels",
           "len_est_n50")
FIXED_CHECKS = ("fc_pores_first_scan", "fc_channels_with_pore_first_scan", "fc_occupancy_first_hour", "fc_occupancy_median",
                "fc_recovering_frac", "fc_recovering_last_third_frac", "fc_temperature_max_dev_c", "fc_active_channels",
                "dyn_speed_drop_frac", "dyn_signal_positive_frac", "dyn_unblock_frac", "len_est_n50")
RULE_CHECKS = ("run_stop_cause", "run_messages", "run_output_files", "dyn_speed_est_median")


def _when_matches(when: dict, family: str | None, agroup: str | None) -> bool:
    for k, v in (when or {}).items():
        have = {"family": family, "assay": agroup}.get(k)
        if have != v:
            return False
    return True


def _catalogue_rule(metric_id: str, family: str | None, agroup: str | None) -> dict | None:
    m = catalogue.metric(metric_id) or {}
    th = m.get("thresholds") or {}
    for r in th.get("rules") or []:
        if _when_matches(r.get("when") or {}, family, agroup):
            return {"kind": th.get("kind"), "warn": r.get("warn"), "fail": r.get("fail")}
    return None


def role_of(metric_id: str, override: dict | None = None) -> str:
    role = (override or {}).get("role") or (catalogue.metric(metric_id) or {}).get("verdict_role") or "info"
    return role if role in ("key", "supporting", "info") else "info"


class Limits:
    """The limits one run is judged with: catalogue defaults for its family and assay group, with the overrides of
    the configuration's ``[profiles.<assay group>]`` table."""

    def __init__(self, family: str | None, agroup: str | None, profile: dict | None = None):
        self.family, self.agroup = family, agroup
        prof = (profile or {}).get(agroup or "", {}) if isinstance(profile, dict) else {}
        self.overrides: dict[str, dict] = {k: dict(v) for k, v in (prof.get("rules") or {}).items() if isinstance(v, dict)}
        self.outcome_unit: str | None = prof.get("outcome_unit") if prof.get("outcome_unit") in ("reads", "bases") else None
        self.custom = bool(self.overrides) or self.outcome_unit is not None

    def enabled(self, metric_id: str) -> bool:
        return self.overrides.get(metric_id, {}).get("enabled", True) is not False

    def role(self, metric_id: str) -> str:
        return role_of(metric_id, self.overrides.get(metric_id))

    def fixed(self, metric_id: str) -> dict | None:
        """``{"kind", "warn", "fail", "configured"}`` of a check with fixed limits, or None when it has none here."""
        base = _catalogue_rule(metric_id, self.family, self.agroup)
        ov = self.overrides.get(metric_id) or {}
        if base is None and not ({"warn", "fail"} & set(ov)):
            return None
        kind = (base or {}).get("kind") or (catalogue.metric(metric_id) or {}).get("thresholds", {}).get("kind") or "min"
        out = {"kind": kind, "warn": (base or {}).get("warn"), "fail": (base or {}).get("fail"), "configured": False}
        for k in ("warn", "fail"):
            if k in ov:
                out[k] = ov[k]
                out["configured"] = True
        return out

    def output_units(self) -> list[tuple[str, str]]:
        """[(unit, role)]: the judged outcome first. Reads for cDNA and direct RNA with estimated bases beside
        them (supporting); estimated bases for everything else."""
        unit = self.outcome_unit or ("reads" if self.agroup in ("cdna", "drna") else "bases")
        if unit == "reads":
            return [("reads", "key"), ("bases", "supporting")]
        return [("bases", "key")]

    def output_levels(self, unit: str, horizon: int) -> tuple[float, float | None] | None:
        """(warn, fail) of the fixed fall-back at a horizon: the 72 h levels times the measured shape of the assay."""
        mid = f"yield_{'reads' if unit == 'reads' else 'estimated_bases'}_at_72h"
        ov = self.overrides.get(f"yield_{'reads' if unit == 'reads' else 'estimated_bases'}_at_{horizon}h") or {}
        if "warn" in ov:
            return float(ov["warn"]), (float(ov["fail"]) if ov.get("fail") is not None else None)
        base = self.fixed(mid)
        th = (catalogue.metric(mid) or {}).get("thresholds") or {}
        factor = ((th.get("horizon_factors") or {}).get(self.agroup or "") or {}).get(str(horizon))
        if base is None or base.get("warn") is None or factor is None:
            return None
        return float(base["warn"]) * factor, (float(base["fail"]) * factor if base.get("fail") is not None else None)

    def warranty(self) -> float | None:
        f = self.fixed("fc_pores_first_scan")
        return f["warn"] if f else None

    def profile(self) -> dict:
        """``status.profile``: which limits judged the run, with a hash that changes when any of them does."""
        table = {m: self.fixed(m) for m in FIXED_CHECKS}
        table.update({f"{u}_{h}": self.output_levels(u, h) for u in ("reads", "bases") for h in HORIZONS})
        table["roles"] = {m: self.role(m) for m in FIXED_CHECKS + RULE_CHECKS}
        table["disabled"] = sorted(m for m in self.overrides if not self.enabled(m))
        table["outcome_unit"] = self.outcome_unit
        # limits the configuration sets for checks judged elsewhere (read data): part of the hash only when there are any
        known = set(FIXED_CHECKS) | set(RULE_CHECKS) | {f"yield_{u}_at_{h}h" for u in ("reads", "estimated_bases") for h in HORIZONS}
        other = {m: {k: v for k, v in sorted(o.items())} for m, o in sorted(self.overrides.items()) if m not in known}
        if other:
            table["other_rules"] = other
        return {"id": self.agroup or "unknown", "version": catalogue.version(), "sha256": hash_obj(table),
                "source": "library_config" if self.custom else "builtin"}


@lru_cache(maxsize=64)
def horizon_factor(agroup: str, horizon: int) -> float | None:
    th = (catalogue.metric("yield_estimated_bases_at_72h") or {}).get("thresholds") or {}
    return ((th.get("horizon_factors") or {}).get(agroup) or {}).get(str(horizon))


def loading_floor_bases(family: str | None, agroup: str | None, frac: float = 0.02) -> float:
    """Estimated bases below which an earlier run was a restart and not a loading: 2 % of the fixed 24 h
    warning level of its device and assay (0.54 Gb for cDNA on PromethION)."""
    lim = Limits(family or "promethion", agroup if agroup in ("cdna", "drna", "gdna") else "cdna")
    lv = lim.output_levels("bases", 24)
    if lv is None:
        lv = Limits("promethion", "cdna").output_levels("bases", 24) or (27e9, None)
    return lv[0] * frac
