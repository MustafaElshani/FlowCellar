"""Partials: the cached result of one pass over one input file.

A partial is ``<partials dir>/<id[:2]>/<id>.fcp.json.gz``. It holds cube fragments in RAW units (samples, events;
the conversion happens at assembly), the routing summary, and enough about its input to describe itself.
Its id is a pure function of the input's fingerprint, the role, the module versions and the cell-shaping
configuration, so two processes computing the same partial write the same file.

Layout of the file::

    {"format": "flowcellar.partial", "schema_version": "0.1.1", "partial_id": "<32 hex>",
     "role": "minknow.sequencing_summary",
     "input": {"path", "size", "mtime", "fingerprint": {"method", "value"}},
     "generator": {"name", "version"}, "modules": {"dyn": "0.1.0", ...}, "config_hash": "<16 hex>",
     "state": "ok" | "partially_read" | "truncated" | "error",
     "coverage": {"kind", "records", "records_total", "fraction"},
     "sample": null | alias, "variant": null | label, "population": null | {...},
     "routing": {"runs": {"<run_key>": records}, "unrouted": 0, "default_run": null | run_key},
     "fragments": [{"run_key", "sample", "module", "dims", "cells", "run_level", "params", "facts"}],
     "facts": {...},            # free-form, per role: what assembly and linking need besides the cells
     "errors": [{"code", "message"}], "timings": {"seconds", "mb_per_s", "workers"}, "created": "..."}
"""
from __future__ import annotations

import os

from flowcellar import SCHEMA_VERSION, __version__
from flowcellar.util import atomic
from flowcellar.util.errors import ValidationError
from flowcellar.util.fingerprint import hash_obj
from flowcellar.util.timeutil import now_z

PARTIAL_FORMAT = "flowcellar.partial"
PARTIAL_SUFFIX = ".fcp.json.gz"
STATES = ("ok", "partially_read", "truncated", "error")


def partial_id(fp: dict, role: str, modules: dict[str, str], config_hash: str = "") -> str:
    """32 hex characters: H(input fingerprint, role, module versions, cell-shaping configuration)."""
    return hash_obj({"fingerprint": [fp["method"], fp["value"]], "role": role,
                     "modules": sorted(modules.items()), "config": config_hash}, 32)


def partial_path(partials_dir: str | os.PathLike, pid: str) -> str:
    return os.path.join(os.fspath(partials_dir), pid[:2], pid + PARTIAL_SUFFIX)


def new_partial(*, role: str, input_path: str, fp: dict, size: int, mtime: str | None, modules: dict[str, str],
                config_hash: str = "", sample: str | None = None, variant: str | None = None) -> dict:
    return {
        "format": PARTIAL_FORMAT, "schema_version": SCHEMA_VERSION,
        "partial_id": partial_id(fp, role, modules, config_hash), "role": role,
        "input": {"path": input_path, "size": size, "mtime": mtime, "fingerprint": fp},
        "generator": {"name": "flowcellar", "version": __version__}, "modules": dict(modules), "config_hash": config_hash,
        "state": "ok", "coverage": {"kind": "full", "records": None, "records_total": None, "fraction": 1.0},
        "sample": sample, "variant": variant, "population": None,
        "routing": {"runs": {}, "unrouted": 0, "default_run": None},
        "fragments": [], "facts": {}, "errors": [], "timings": {}, "created": now_z(),
    }


def fragment(run_key: str | None, module: str, *, sample: str | None = None, dims: list[str] | None = None,
             cells: list[dict] | None = None, run_level: dict | None = None, params: dict | None = None,
             facts: dict | None = None) -> dict:
    """One piece of one section addressed to one run (and one sample alias, when the input is per sample)."""
    return {"run_key": run_key, "sample": sample, "module": module, "dims": dims or [], "cells": cells or [],
            "run_level": run_level, "params": params or {}, "facts": facts or {}}


def write_partial(partials_dir: str | os.PathLike, partial: dict) -> str:
    """Write under ``<partials_dir>/<hh>/``; returns the path. A partial of a growing or running input is never
    cached: the caller must not call this for one (state partially_read from a time limit is cached, with its
    coverage, and replaced by the next full pass because ``is_complete`` is false)."""
    path = partial_path(partials_dir, partial["partial_id"])
    atomic.atomic_write_json(path, partial)
    return path


def read_partial(path: str | os.PathLike) -> dict:
    d = atomic.read_json(path)
    if not isinstance(d, dict) or d.get("format") != PARTIAL_FORMAT:
        raise ValidationError(f"{path}: not a FlowCellar partial")
    return d


def find_partial(partials_dir: str | os.PathLike, pid: str) -> str | None:
    p = partial_path(partials_dir, pid)
    return p if os.path.isfile(p) else None


def is_complete(partial: dict) -> bool:
    return partial.get("state") == "ok" and (partial.get("coverage") or {}).get("kind") == "full"


def fragments_for(partial: dict, run_key: str, module: str | None = None) -> list[dict]:
    return [f for f in partial.get("fragments", [])
            if f.get("run_key") == run_key and (module is None or f.get("module") == module)]
