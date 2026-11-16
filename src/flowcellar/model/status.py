"""Verdicts and the status file.

The status file ``status/<run_key>.status.json`` is ``{"format": "flowcellar.status", "schema_version",
"bundle_id", "content_sha256", "status": {...}, "relations": {...}}``: how one bundle was judged, and its
library relations. It is rewritten when the bundle changed or on ``recompute --status-only``, never because a
later run arrived.
"""
from __future__ import annotations

import os

from flowcellar import SCHEMA_VERSION
from flowcellar.util import atomic
from flowcellar.util.errors import ValidationError

STATUS_FORMAT = "flowcellar.status"
STATUS_SUFFIX = ".status.json"

VERDICTS = ("pass", "warn", "fail", "running", "incomplete", "unknown")
TIER_VALUES = VERDICTS + ("not_available",)
VERDICT_WORD = {"pass": "Pass", "warn": "Warning", "fail": "Fail", "running": "Running", "incomplete": "Incomplete",
                "unknown": "Not assessed", "not_available": "Not available yet"}
SEVERITIES = ("none", "info", "warn", "fail")
_RANK = {"pass": 0, "warn": 1, "fail": 2}

CAUSES = ("stopped_early", "output_not_judged", "no_final_summary", "raw_only", "inputs_missing")
MODIFIERS = ("reused_flow_cell", "flow_cell_use_unknown", "baseline_provisional", "baseline_spread_guard", "no_history",
             "judged_before_older_runs", "custom_profile", "quick_detail", "reads_partly_counted")


def status_filename(bundle_id: str) -> str:
    return bundle_id + STATUS_SUFFIX


def rule_level(rule: dict) -> str:
    """pass, warn or fail: what one stored rule contributes. Only a key rule can contribute fail."""
    if rule.get("outcome") != "triggered":
        return "pass"
    sev = rule.get("severity")
    if sev == "fail":
        return "fail" if rule.get("role") == "key" else "warn"
    return "warn" if sev == "warn" else "pass"


def worst_of_rules(rules, tier: str | None = None) -> str:
    """The worst outcome among key and supporting rules of one part (or of all): pass, warn or fail."""
    v = "pass"
    for r in rules:
        if tier is not None and r.get("tier") != tier:
            continue
        if r.get("role") == "info":
            continue
        lv = rule_level(r)
        if _RANK[lv] > _RANK[v]:
            v = lv
    return v


def combine_verdicts(run_health: str, read_level: str) -> str:
    """The overall verdict from its two parts: the worse of the two, leaving out a part without data.

    not_available and unknown never lower a verdict. An incomplete Run health stays incomplete unless Read data
    warns or fails.
    """
    quiet = ("not_available", "unknown", "pass")
    if run_health == "not_available":
        return "unknown" if read_level == "not_available" else read_level
    if read_level in quiet:
        return run_health
    if run_health in ("incomplete", "running", "unknown"):
        return read_level if read_level in ("warn", "fail") else run_health
    return read_level if _RANK.get(read_level, 0) > _RANK.get(run_health, 0) else run_health


def reasons(status: dict, limit: int = 3) -> list[dict]:
    """The triggered checks behind a verdict, most important first, as the index stores them:
    ``[{"check": id, "severity": "warn" | "fail", "tier": ..., "message": sentence}]``."""
    trig = [r for r in status.get("rules", []) if r.get("outcome") == "triggered" and r.get("severity") in ("warn", "fail")]
    trig.sort(key=lambda r: (r.get("tier") != "run_health", r.get("role") != "key", r.get("severity") != "fail"))
    return [{"check": r["id"], "severity": rule_level(r), "tier": r.get("tier"), "message": r.get("message", "")}
            for r in trig[:limit]]


def new_status_file(bundle_id: str, content_sha256: str | None, status: dict, relations: dict) -> dict:
    return {"format": STATUS_FORMAT, "schema_version": SCHEMA_VERSION, "bundle_id": bundle_id,
            "content_sha256": content_sha256, "status": status, "relations": relations}


def write_status(path: str | os.PathLike, status_file: dict) -> None:
    atomic.atomic_write_json(path, status_file, indent=1)


def read_status(path: str | os.PathLike) -> dict:
    d = atomic.read_json(path)
    if not isinstance(d, dict) or d.get("format") != STATUS_FORMAT:
        raise ValidationError(f"{path}: not a FlowCellar status file")
    return d
