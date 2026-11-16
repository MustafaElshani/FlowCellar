"""Bundles: the envelope, the split into measured content and judgement, reading and writing.

A bundle file inside a library (``bundles/<run_key>.fcb.json.gz``) holds measured content only. ``status``,
``user`` and ``relations`` are the *composed* members: they live in ``status/<run_key>.status.json`` (and, in a
later build, annotations) and are put back into every bundle that leaves the library.
"""
from __future__ import annotations

import copy
import os
import platform
import sys

from flowcellar import SCHEMA_VERSION, __version__
from flowcellar.util import atomic, fingerprint
from flowcellar.util.errors import ValidationError
from flowcellar.util.timeutil import now_z

FORMAT = "flowcellar.bundle"
BUNDLE_SUFFIX = ".fcb.json.gz"
REPORT_SUFFIX = ".flowcellar.html"
COMPOSED = ("status", "user", "relations")
MEMBER_ORDER = ("format", "schema_version", "kind", "bundle_id", "content_sha256", "identity", "provenance", "tiers",
                "samples", "axes", "modules", "status", "user", "relations", "extensions")
# left out of content_sha256: library knowledge, paths and times (two computations of the same inputs compare equal)
_PROVENANCE_OUTSIDE_HASH = ("created", "updated", "revision", "timings", "library", "run_dir")
_INPUT_OUTSIDE_HASH = ("path", "root", "relpath", "alt_paths", "mtime")


def bundle_filename(bundle_id: str) -> str:
    return bundle_id + BUNDLE_SUFFIX


def report_filename(prefix: str) -> str:
    return prefix + REPORT_SUFFIX


def generator() -> dict:
    return {"name": "flowcellar", "version": __version__,
            "python": "%d.%d.%d" % sys.version_info[:3], "platform": platform.system().lower()}


def new_bundle(identity: dict, *, kind: str = "run") -> dict:
    """An empty, valid envelope around ``identity`` (which must hold ``run_key`` for a run)."""
    bundle_id = identity["run_key"] if kind == "run" else identity["analysis"]["id"]
    now = now_z()
    return {
        "format": FORMAT, "schema_version": SCHEMA_VERSION, "kind": kind, "bundle_id": bundle_id,
        "identity": identity,
        "provenance": {"generator": generator(), "created": now, "updated": now, "inputs": [], "partials": [], "warnings": []},
        "tiers": {"run_health": {"state": "absent"}, "read_level": {"state": "absent"}},
        "samples": [], "axes": {}, "modules": {},
    }


def order_members(bundle: dict) -> dict:
    """The same bundle with its top-level members in the fixed order (``format`` and ``schema_version`` first)."""
    out = {k: bundle[k] for k in MEMBER_ORDER if k in bundle}
    out.update((k, v) for k, v in bundle.items() if k not in out)
    return out


def default_user() -> dict:
    return {"tags": [], "notes": "", "display_name": None, "fields": {}, "verdict_override": None, "updated": None,
            "exclude_from_baselines": False, "baseline_group": None, "flow_cell_use": None, "flow_cell_check_pores": None}


def default_relations() -> dict:
    return {"analyses": [], "superseded_by": None,
            "flow_cell_use": {"index": 1, "of": None, "source": "unknown", "prior_runs": None, "prior_hours": None,
                              "prior_estimated_bases": None, "restart_of": None},
            "projects": [], "paths": []}


def split(bundle: dict) -> tuple[dict, dict]:
    """(measured content, the composed members that were present). Neither shares structure with the input's top level."""
    measured = {k: v for k, v in bundle.items() if k not in COMPOSED}
    composed = {k: bundle[k] for k in COMPOSED if k in bundle}
    return order_members(measured), composed


def compose(measured: dict, *, status: dict | None = None, user: dict | None = None,
            relations: dict | None = None) -> dict:
    """The bundle as it leaves a library: measured content plus status, user and relations."""
    out = {k: v for k, v in measured.items() if k not in COMPOSED}
    if status is not None:
        out["status"] = status
    out["user"] = user if user is not None else default_user()
    out["relations"] = relations if relations is not None else default_relations()
    return order_members(out)


def content_sha256(bundle: dict) -> str:
    """Hash of the measured content: no status, user, relations, extensions, paths or times."""
    b = {k: v for k, v in bundle.items() if k not in COMPOSED and k not in ("content_sha256", "extensions")}
    prov = copy.deepcopy(b.get("provenance") or {})
    for k in _PROVENANCE_OUTSIDE_HASH:
        prov.pop(k, None)
    (prov.get("generator") or {}).pop("python", None)
    (prov.get("generator") or {}).pop("platform", None)
    for i in prov.get("inputs") or []:
        for k in _INPUT_OUTSIDE_HASH:
            i.pop(k, None)
    b["provenance"] = prov
    return fingerprint.sha256_hex(fingerprint.canonical_json(b))


def check_version(bundle: dict, where: str = "bundle") -> None:
    """Refuse a file that is not a FlowCellar bundle or was written by a newer schema (exit code 7)."""
    if bundle.get("format") != FORMAT:
        raise ValidationError(f"{where}: not a FlowCellar bundle (format is {bundle.get('format')!r})")
    v = str(bundle.get("schema_version", ""))
    mine = tuple(int(x) for x in SCHEMA_VERSION.split("."))
    try:
        theirs = tuple(int(x) for x in v.split("."))
    except ValueError:
        raise ValidationError(f"{where}: unreadable schema_version {v!r}") from None
    if theirs[:2] > mine[:2]:
        raise ValidationError(f"{where}: written with schema {v}; this FlowCellar reads up to {SCHEMA_VERSION}",
                              hint="install a newer FlowCellar")


def write_bundle(path: str | os.PathLike, bundle: dict, *, with_hash: bool = True) -> dict:
    """Write ``bundle`` as gzip JSON (level 9, header time 0, fixed member order). Returns what was written."""
    b = dict(bundle)
    if with_hash:
        b["content_sha256"] = content_sha256(b)
    b = order_members(b)
    atomic.atomic_write_json(path, b)
    return b


def to_gzip_bytes(bundle: dict) -> bytes:
    """The bytes of a ``.fcb.json.gz``: what a report embeds and a hub sidecar carries (base64 of these)."""
    return atomic.gzip_bytes(atomic.dumps(order_members(bundle)).encode("utf-8"))


def read_bundle(path: str | os.PathLike) -> dict:
    b = atomic.read_json(path)
    if not isinstance(b, dict):
        raise ValidationError(f"{path}: not a FlowCellar bundle")
    check_version(b, os.fspath(path))
    return b
