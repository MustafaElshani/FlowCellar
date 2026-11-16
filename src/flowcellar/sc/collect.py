"""One analysis folder -> the measured content of its bundle (``kind: "analysis"``).

``collect(claim, options)`` reads the files of one claim (``flowcellar.sc.detect``) and returns the bundle
without its verdict. It needs no library: what it holds is what the files of the analysis say, with the links
the analysis carries by itself. A library adds the links that need its runs, and the verdict
(``flowcellar.sc.store``, ``flowcellar.sc.judge``).

The section ``modules.sc`` is one cell of the generic aggregates of the bundle schema:

    scalars   counts as read (reads, molecules, units)
    hists     per-unit distributions and reads per molecule, on the three fixed axes
    counts    the read funnel, the adapter layouts
    curves    the barcode rank plot
    tables    bin sizes, segmentation, registration, parameters, the slowest workflow steps
    grids     the spatial map
    summary   every number a tile, a check or the index shows (fractions, medians, values as the pipeline gave them)
    notes     the pipeline's own alarms, word for word, and what could not be read
    value_sources   the values that were taken from a report's HTML and not from a data file
"""
from __future__ import annotations

import hashlib
import os
import time
from dataclasses import dataclass
from typing import Callable

import numpy as np

from flowcellar.model import bundle as bundle_model
from flowcellar.sc import (
    MODULE_ID,
    MODULE_VERSION,
    PIPELINE_SPACERANGER,
    SECTION_CAP_BYTES,
    SECTION_SCHEMA_VERSION,
    axes,
)
from flowcellar.sc import bundle_id as make_bundle_id
from flowcellar.sc.detect import Claim
from flowcellar.util import atomic
from flowcellar.util.errors import ResourceLimitError
from flowcellar.util.text import clean_text
from flowcellar.util.timeutil import from_epoch

UNIT = {"singlecell": {"name": "cell", "size_um": None}, "visium": {"name": "spot", "size_um": 55.0},
        "visium_hd": {"name": "bin", "size_um": 8.0}}


@dataclass
class Options:
    read_pass: bool = True                  # the streaming pass over the workflow's per-read table
    dedup: bool = True                      # reads per molecule from that pass (the saturation curve of the workflow)
    map_um: int = 64                        # bin size of a Visium HD map: 64 or 32
    tissue_image: bool = True               # the grey tissue thumbnail of a Visium HD map
    link_tags: bool = True                  # tally the read tags of the workflow's BAM
    workers: int = 0
    deadline: float | None = None           # time.time() value after which the read pass stops
    label: str | None = None
    sibling: str | None = None              # --sibling DIR: the Space Ranger outs/ of the same reads
    progress: Callable[[str], None] | None = None

    def shaping(self) -> dict:
        """The settings that change what a collection holds: part of its cache key."""
        return {"read_pass": self.read_pass, "dedup": self.dedup, "map_um": self.map_um, "tissue_image": self.tissue_image,
                "link_tags": self.link_tags, "label": self.label, "sibling": self.sibling}


class Builder:
    """Gathers the parts of one section while the files are read."""

    def __init__(self, claim: Claim, opts: Options):
        self.claim, self.opts = claim, opts
        self.scalars: dict = {}
        self.summary: dict = {}
        self.hists: dict = {}
        self.counts: dict = {}
        self.curves: dict = {}
        self.tables: dict = {}
        self.grids: dict = {}
        self.notes: list[dict] = []
        self.sources: dict = {}
        self.inputs: list[dict] = []
        self.params: dict = {}
        self.timings: dict = {}
        self.identity: dict = {}
        self.links: list[dict] = []
        self._input_ids: dict[str, str] = {}

    # ---- inputs -------------------------------------------------------------------------------------------
    def input(self, role: str, path: str, *, state: str = "ok", records: int | None = None, note: str | None = None,
              fingerprint: bool = False) -> str:
        """Record one file that was read (or could not be) and return its input id."""
        if path in self._input_ids:
            for i in self.inputs:
                if i["path"] == path:
                    if state != "ok":
                        i["state"] = state
                    if records is not None:
                        i["records"] = int(records)
                    if note:
                        i.setdefault("notes", []).append({"level": "info", "code": "input_note", "message": clean_text(note)})
            return self._input_ids[path]
        iid = f"an{len(self.inputs) + 1}"
        self._input_ids[path] = iid
        entry: dict = {"id": iid, "role": role, "path": path, "state": state, "used_by": [MODULE_ID]}
        try:
            st = os.stat(path)
            entry["size"], entry["mtime"] = int(st.st_size), from_epoch(st.st_mtime)
            if fingerprint and st.st_size <= (8 << 20):
                with open(path, "rb") as fh:
                    entry["fingerprint"] = {"method": "sha256", "value": hashlib.sha256(fh.read()).hexdigest()}
        except OSError:
            entry["state"] = "missing"
        if records is not None:
            entry["records"] = int(records)
        if note:
            entry["notes"] = [{"level": "info", "code": "input_note", "message": clean_text(note)}]
        self.inputs.append(entry)
        return iid

    # ---- values -------------------------------------------------------------------------------------------
    def put(self, key: str, value, *, kind: str = "data", input_id: str | None = None) -> None:
        """One summary value. ``kind`` other than data marks it as taken from a report's HTML."""
        if value is None:
            return
        if isinstance(value, (np.integer,)):
            value = int(value)
        elif isinstance(value, (np.floating,)):
            value = float(value)
        self.summary[key] = value
        if kind != "data":
            self.sources[key] = {"input": input_id, "kind": kind}
        else:
            self.sources.pop(key, None)

    def scalar(self, key: str, value) -> None:
        if value is not None:
            self.scalars[key] = int(value) if float(value).is_integer() else float(value)

    def stats(self, prefix: str, values) -> dict:
        s = axes.summary_stats(values)
        for k in ("median", "mean", "p10", "p90"):
            if s[k] is not None:
                self.summary[f"{prefix}_{k}"] = float(s[k])
        return s

    def note(self, level: str, code: str, message: str, **extra) -> None:
        n = {"level": level, "code": code, "message": clean_text(message)}
        n.update({k: v for k, v in extra.items() if v is not None})
        self.notes.append(n)

    def table(self, key: str, columns: list[tuple[str, str]], rows: list[list], *, unit: dict | None = None) -> None:
        cols = [{"id": c, "type": t, **({"unit": unit[c]} if unit and c in unit else {})} for c, t in columns]
        self.tables[key] = {"columns": cols, "merge": "replace", "rows": rows}

    def say(self, message: str) -> None:
        if self.opts.progress:
            self.opts.progress(message)

    def timed(self, name: str):
        return _Timer(self, name)


class _Timer:
    def __init__(self, b: Builder, name: str):
        self.b, self.name = b, name

    def __enter__(self):
        self.t0 = time.monotonic()
        return self

    def __exit__(self, *exc):
        self.b.timings[self.name] = round(self.b.timings.get(self.name, 0.0) + time.monotonic() - self.t0, 3)
        return False


# ---------------------------------------------------------------------------------------------------------
# identity of the analysis itself
# ---------------------------------------------------------------------------------------------------------

def _read_bytes(path: str, limit: int = 4 << 20) -> bytes:
    with open(path, "rb") as fh:
        return fh.read(limit)


def analysis_hash(claim: Claim, sample_id: str | None = None) -> str:
    """SHA-1 over the pipeline's name, the sample id and the bytes of the primary metrics file, so that copies
    of one output folder are one analysis and a re-run is a new one. For Space Ranger the file is
    ``metrics_summary.csv`` (its first column is the sample id); for wf-single-cell it is ``config_stats.json``
    or ``params.json``, with the size line of the gene matrix."""
    h = hashlib.sha1()
    h.update(claim.pipeline.encode() + b"\0")
    if claim.pipeline == PIPELINE_SPACERANGER:
        h.update((sample_id or "").encode("utf-8") + b"\0")
        h.update(_read_bytes(claim.files["metrics_summary"]))
    else:
        from flowcellar.sc import wfsc
        h.update((claim.alias or "").encode("utf-8") + b"\0")
        for role in ("config_stats", "params"):
            if role in claim.files:
                h.update(_read_bytes(claim.files[role]))
                break
        for role in ("gene_matrix_8um", "gene_matrix", "gene_matrix_2um"):
            if role in claim.files:
                h.update(("%d %d %d" % wfsc.mtx_size_line(claim.files[role])).encode())
                break
    return h.hexdigest()


def quick_id(claim: Claim) -> str:
    """The bundle id of a claim, from its primary metrics file alone (two small reads)."""
    sample = None
    if claim.pipeline == PIPELINE_SPACERANGER:
        from flowcellar.sc import spaceranger
        sample = str(spaceranger.metrics_summary(claim.files["metrics_summary"]).get("Sample ID") or "")
    return make_bundle_id(claim.pipeline, analysis_hash(claim, sample))


def signature(claim: Claim) -> list:
    """Name, size and modification time of every file of a claim: "did this analysis change" without opening it."""
    out = []
    for role, path in sorted(claim.files.items()):
        try:
            st = os.stat(os.path.join(path, "matrix.mtx.gz") if os.path.isdir(path) else path)
            out.append([role, int(st.st_size), int(st.st_mtime_ns)])
        except OSError:
            out.append([role, None, None])
    return out


# ---------------------------------------------------------------------------------------------------------
# the bundle
# ---------------------------------------------------------------------------------------------------------

def collect(claim: Claim, opts: Options | None = None) -> dict:
    """Read one analysis and return its measured bundle."""
    opts = opts or Options()
    b = Builder(claim, opts)
    t0 = time.monotonic()
    if claim.pipeline == PIPELINE_SPACERANGER:
        from flowcellar.sc import collect_sr
        collect_sr.collect(b)
    else:
        from flowcellar.sc import collect_wf
        collect_wf.collect(b)
    b.timings["total"] = round(time.monotonic() - t0, 3)
    return _bundle(b)


def _section(b: Builder) -> dict:
    cell: dict = {"at": []}
    for name in ("scalars", "hists", "counts", "curves", "tables", "grids"):
        part = getattr(b, name)
        if part:
            cell[name] = part
    missing = [n for n in b.notes if n["code"] in ("h5py_missing", "read_pass_skipped", "read_pass_partial", "input_unreadable")]
    section: dict = {"version": MODULE_VERSION, "tier": "assay", "state": "partial" if missing else "ok", "dims": [],
                     "cells": [cell], "summary": dict(sorted(b.summary.items())),
                     "params": {"section_schema_version": SECTION_SCHEMA_VERSION, **b.params},
                     "inputs": [i["id"] for i in b.inputs]}
    if missing:
        section["partial_reason"] = "partly_read" if any(n["code"] == "read_pass_partial" for n in missing) else "inputs_missing"
    if b.sources:
        section["value_sources"] = dict(sorted(b.sources.items()))
    if b.notes:
        section["notes"] = b.notes
    return section


def _bundle(b: Builder) -> dict:
    claim = b.claim
    ident = b.identity
    a_hash = analysis_hash(claim, ident.get("sample_id"))
    assay = ident.get("assay")
    analysis = {
        "id": a_hash, "pipeline": claim.pipeline, "pipeline_version": ident.get("pipeline_version"),
        "label": clean_text(b.opts.label or ident.get("label") or claim.alias or claim.variant, 256),
        "path": claim.path, "sibling_analyses": [], "links": b.links,
        "group_key": clean_text(ident.get("slide") or b.opts.label or ident.get("label") or claim.alias or claim.variant, 256),
        "variant": clean_text(claim.variant, 256), "assay": assay,
        "unit": dict(UNIT[assay]) if assay in UNIT else None,
        "slide_serial_capture_area": ident.get("slide"), "reference": ident.get("reference"),
    }
    bundle = bundle_model.new_bundle({"analysis": {**analysis, "id": a_hash}}, kind="analysis")
    bundle["bundle_id"] = make_bundle_id(claim.pipeline, a_hash)
    section = _section(b)
    size = len(atomic.dumps(section))
    section["params"]["section_bytes"] = size
    if size > SECTION_CAP_BYTES:
        raise ResourceLimitError(f"{claim.path}: the analysis section would be {size / 1000:.0f} kB; the limit is "
                                 f"{SECTION_CAP_BYTES // 1000} kB", hint="use the default map size (--map-um 64)")
    bundle["modules"] = {MODULE_ID: section}
    bundle["axes"] = axes.axis_defs(h["axis"] for h in b.hists.values())
    prov = bundle["provenance"]
    prov["inputs"] = b.inputs
    prov["timings"] = {k: float(v) for k, v in b.timings.items()}
    prov["run_dir"] = {"path": claim.path, "layout": "unknown"}
    return bundle
