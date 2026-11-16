"""The catalogue of the add-on: metric ids, labels, default limits and the sentences of its checks.

The source is ``catalogue/metrics_singlecell_spatial.yaml`` in the repository; ``tools/build_sc_catalogue.py``
writes the JSON beside this file, and only that JSON is read at run time.
"""
from __future__ import annotations

import json
from functools import lru_cache
from importlib import resources

CATALOGUE_FILE = "catalogue.json"


@lru_cache(maxsize=1)
def load() -> dict:
    return json.loads(resources.files(__package__).joinpath(CATALOGUE_FILE).read_text("utf-8"))


@lru_cache(maxsize=1)
def metrics() -> dict[str, dict]:
    return {m["id"]: m for m in load()["metrics"]}


def metric(metric_id: str) -> dict | None:
    return metrics().get(metric_id)


def label(metric_id: str) -> str:
    m = metrics().get(metric_id) or {}
    return (m.get("display") or {}).get("label") or metric_id


def checks() -> list[dict]:
    """The checks in the order they are judged and listed."""
    return list(load().get("checks") or [])


def limits(metric_id: str) -> dict:
    """The default limits of one metric (``warn``, ``fail``, ``warn_low``, ``warn_high``, ``direction``)."""
    th = (metrics().get(metric_id) or {}).get("thresholds")
    return dict(th) if isinstance(th, dict) else {}


def unit_words(assay: str | None) -> dict:
    """``{"name", "plural"}`` of what a barcode stands for in one assay."""
    return dict((load().get("units") or {}).get(assay or "", {"name": "unit", "plural": "units"}))


def version() -> str:
    return str(load().get("catalogue_version", "0"))
