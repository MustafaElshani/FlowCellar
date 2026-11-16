"""The metric catalogue: ids, sections, limits, display copy, the panel registry and the tile order.

The source is ``catalogue/metrics_nanopore.yaml`` in the repository; ``tools/build_catalogue.py`` writes the JSON
beside this file, and only that JSON is read at run time (no YAML parser is a run-time dependency).
"""
from __future__ import annotations

import json
from functools import lru_cache
from importlib import resources

CATALOGUE_FILE = "metrics_nanopore.json"


@lru_cache(maxsize=1)
def load() -> dict:
    """The whole catalogue as one dict (keys: metrics, derived, panels, tiles, limit_text, axes, modules, ...)."""
    return json.loads(resources.files(__name__).joinpath(CATALOGUE_FILE).read_text("utf-8"))


@lru_cache(maxsize=1)
def metrics() -> dict[str, dict]:
    """Every metric and derived summary id -> its entry."""
    cat = load()
    out = {m["id"]: m for m in cat["metrics"]}
    out.update({d["id"]: d for d in cat.get("derived", [])})
    return out


def metric(metric_id: str) -> dict | None:
    return metrics().get(metric_id)


def display(metric_id: str) -> dict:
    """The display copy of a metric (label, short_label, reading, message, advice, panel); empty when it has none."""
    m = metrics().get(metric_id) or {}
    return m.get("display") or {}


def label(metric_id: str) -> str:
    return display(metric_id).get("label") or metric_id


def checks(tier: int | None = None) -> list[dict]:
    """The key and supporting metrics (the ones a status rule exists for), optionally of one tier."""
    return [m for m in load()["metrics"]
            if m.get("verdict_role") in ("key", "supporting") and (tier is None or m.get("tier") == tier)]


def limit_text() -> dict[str, str]:
    return load()["limit_text"]


def panels() -> list[dict]:
    return load()["panels"]


def tiles() -> dict[str, list[str]]:
    return load()["tiles"]


def version() -> str:
    return str(load().get("catalogue_version", "0"))
