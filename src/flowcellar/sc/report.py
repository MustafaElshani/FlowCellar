"""Pages and data files of analyses: the one-file report of ``flowcellar analysis`` and the hub's data files.

Both use the same built viewer and the same page parts as a run (``flowcellar.report.html``): the data is
base64 of gzip bytes in script elements the browser does not run, the title and the text summary are escaped,
and the page asks the network for nothing.

    <prefix>.flowcellar.html        one analysis (and the analysis of the same reads by the other pipeline, when
                                    it was read with it): the index of the file and the bundles, inside the file
    hub/data/<analysis id>.js       one analysis of a library, as one call of ``__FCB_SIDECAR__``
"""
from __future__ import annotations

import hashlib
import json
import os

from flowcellar import SCHEMA_VERSION, __version__
from flowcellar.contracts import Served
from flowcellar.model import bundle as bundle_model
from flowcellar.model import status as status_model
from flowcellar.report import html as page
from flowcellar.sc import ASSAY_WORD, MODULE_ID, PIPELINE_WORD, is_analysis_id
from flowcellar.sc import store as store_mod
from flowcellar.util import atomic
from flowcellar.util.errors import FlowCellarError, InputError
from flowcellar.util.timeutil import now_z

JS = "text/javascript; charset=utf-8"
SIDECAR_SUFFIX = ".js"


# ---------------------------------------------------------------------------------------------------------
# the index of a file that holds analyses and no run
# ---------------------------------------------------------------------------------------------------------

def index_for(bundles: list[dict]) -> dict:
    """An index in the shape the viewer reads, with no run and one row per analysis."""
    rows = []
    for b in bundles:
        sf = None
        if b.get("status") is not None:
            sf = status_model.new_status_file(b["bundle_id"], b.get("content_sha256"), b["status"], {})
        rows.append(store_mod.card_for(b, sf)["row"])
    return {
        "format": "flowcellar.index", "schema": SCHEMA_VERSION, "written": now_z(),
        "meta": {"library": None, "runs": 0, "generator": {"name": "flowcellar", "version": __version__}},
        "library": {"id": None, "name": None, "path": None, "created": None, "scanned": None,
                    "counts": {"runs": 0, "flow_cells": 0, "projects": 0, "samples": 0, "by_verdict": {}, "by_detail": {},
                               "with_read_data": 0}, "estimated_bases": None, "reads": None, "first_run": None, "last_run": None},
        "runs": {"n": 0, "columns": {}}, "samples": [], "projects": [], "roots": [],
        "activity": {"months": [], "runs": [], "estimated_bases": [], "reads": []}, "baselines": {"groups": []}, "facets": {},
        "notices": [], "progress": {"state": "idle"}, "analyses": rows,
    }


# ---------------------------------------------------------------------------------------------------------
# the page
# ---------------------------------------------------------------------------------------------------------

def _boot(analysis_id: str) -> str:
    if not is_analysis_id(analysis_id):
        raise InputError(f"not an analysis id: {analysis_id!r}")
    text = json.dumps({"mode": "report", "generator": __version__, "schema": SCHEMA_VERSION, "run": None, "analysis": analysis_id},
                      ensure_ascii=True, separators=(",", ":"))
    return text.replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026")


def _pct(x) -> str:
    return f"{float(x) * 100:.1f} %"


def key_numbers(bundle: dict) -> list[tuple[str, str]]:
    s = ((bundle.get("modules") or {}).get(MODULE_ID) or {}).get("summary") or {}
    unit = {"bin": "8 µm bin", "spot": "spot", "cell": "cell"}.get(s.get("unit_name"), "unit")
    rows: list[tuple[str, str]] = []
    if s.get("spatial_units_under_tissue") is not None:
        rows.append((f"{unit.capitalize()}s under tissue", f"{int(s['spatial_units_under_tissue']):,}"))
    elif s.get("barcodes_cells") is not None:
        rows.append((f"{unit.capitalize()}s", f"{int(s['barcodes_cells']):,}"))
    for key, label, fmt in (("spatial_umis_per_unit", f"UMIs per {unit}, {s.get('unit_stat', 'median')}", "{:,.1f}"),
                            ("barcodes_umis_median", f"UMIs per {unit}, median", "{:,.0f}"),
                            ("spatial_genes_per_unit", f"Genes per {unit}, {s.get('unit_stat', 'median')}", "{:,.1f}"),
                            ("barcodes_genes_median", f"Genes per {unit}, median", "{:,.0f}")):
        if s.get(key) is not None and not any(r[0].split(",")[0] == label.split(",")[0] for r in rows):
            rows.append((label, fmt.format(s[key])))
    for key, label in (("spatial_fraction_reads_under_tissue", "Reads under tissue"), ("sr_saturation", "Saturation, Space Ranger"),
                       ("saturation_value", "Saturation, ONT workflow"), ("reads_full_length_fraction", "Full-length reads"),
                       ("reads_gene_assigned_fraction", "Reads assigned to a gene")):
        if s.get(key) is not None:
            rows.append((label, _pct(s[key])))
    return rows[:8]


def report_title(bundle: dict) -> str:
    ana = (bundle.get("identity") or {}).get("analysis") or {}
    word = status_model.VERDICT_WORD.get((bundle.get("status") or {}).get("verdict"), "")
    return " · ".join(x for x in (ana.get("label") or bundle.get("bundle_id"), word, "FlowCellar analysis report") if x)


def report_summary(bundles: list[dict]) -> str:
    """The text summary of the page: what a mail preview and a reader without scripts see."""
    esc = page.esc
    out = ['<div class="fc-page fc-nojs">', '<p class="fc-caps">FlowCellar analysis report</p>']
    for i, b in enumerate(bundles):
        ana, st = (b.get("identity") or {}).get("analysis") or {}, b.get("status") or {}
        if i == 0:
            out.append(f'<h1 class="fc-runhead__id">{esc(ana.get("label") or b.get("bundle_id"))}</h1>')
        what = " · ".join(str(x) for x in (ASSAY_WORD.get(ana.get("assay"), ""), PIPELINE_WORD.get(ana.get("pipeline"), ana.get("pipeline")),
                                           ana.get("pipeline_version"), ana.get("slide_serial_capture_area")) if x)
        out.append(f'<p class="fc-runhead__line">{esc(what)}</p>')
        out.append('<section class="fc-verdict"><div class="fc-verdict__head">' + page._badge(st.get("verdict") or "unknown", large=True)
                   + f'<div><p class="fc-verdict__headline">{esc(st.get("headline"))}</p>'
                   + (f'<p class="fc-verdict__detail">{esc(st.get("detail"))}</p>' if st.get("detail") else "") + "</div></div>")
        trig = [r for r in st.get("rules") or [] if r.get("outcome") == "triggered" and r.get("severity") in ("warn", "fail") and r.get("role") != "info"]
        if trig:
            out.append('<ul class="fc-reasons">' + "".join(f"<li><span></span><span>{esc(r.get('message'))}</span></li>" for r in trig[:6]) + "</ul>")
        out.append("</section>")
        nums = key_numbers(b)
        if nums:
            out.append('<div class="fc-keynumbers"><div class="fc-tiles">' + "".join(
                f'<div class="fc-tile"><span class="fc-tile__label"><span>{esc(k)}</span></span><span class="fc-tile__value">{esc(v)}</span></div>'
                for k, v in nums) + "</div></div>")
    out.append('<p class="fc-nojs__note fc-secondary">This is the text summary of the report. The map, the charts and every other number '
               "appear when the file is opened in a browser with scripts switched on (Firefox 113, Chrome or Edge 80, Safari 16.4, or newer). "
               f"Self-contained: works offline, sends nothing. Made by FlowCellar {esc(__version__)}.</p>")
    out.append("</div>")
    return "".join(out)


def report_page(bundles: list[dict]) -> tuple[str, int]:
    """The one-file report of one analysis (``bundles[0]``) and, when given, of its sibling. Each bundle is
    composed (it carries ``status``). Returns the page text and the size of what it embeds."""
    if not bundles:
        raise InputError("no analysis to write a report for")
    for b in bundles:
        if b.get("kind") != "analysis" or not is_analysis_id(b.get("bundle_id")):
            raise InputError(f"cannot write an analysis report: {b.get('bundle_id')!r} is not an analysis")
        if "status" not in b:
            raise InputError(f"cannot write a report for {b['bundle_id']}: the bundle carries no verdict")
    v = page.load_viewer()
    gz_index = atomic.gzip_bytes(atomic.dumps(index_for(bundles)).encode("utf-8"))
    payload = page.data_element(page.INDEX_TYPE, "fci", gz_index)
    total = len(gz_index)
    for b in bundles:
        gz = bundle_model.to_gzip_bytes(b)
        total += len(gz)
        payload += page.data_element(page.BUNDLE_TYPE, "fcb:" + b["bundle_id"], gz, sha256=True)
    text = page.fill(v, mode="report", title=report_title(bundles[0]), summary=report_summary(bundles), boot=_boot(bundles[0]["bundle_id"]),
                     payload=payload)
    return text, total


def write_report(bundles: list[dict], out_path: str) -> str:
    text, _ = report_page(bundles)
    atomic.atomic_write_bytes(os.fspath(out_path), text.encode("utf-8"))
    return os.fspath(out_path)


# ---------------------------------------------------------------------------------------------------------
# the hub's data files
# ---------------------------------------------------------------------------------------------------------

def sidecar(analysis_id: str, gz_bundle: bytes) -> bytes:
    """``hub/data/<analysis id>.js``: the same call a run's data file makes."""
    if not is_analysis_id(analysis_id):
        raise InputError(f"not an analysis id: {analysis_id!r}")
    return ('__FCB_SIDECAR__("' + analysis_id + '","' + page.b64(gz_bundle) + '","' + hashlib.sha256(gz_bundle).hexdigest() + '");\n').encode("ascii")


def sidecar_path(hub_dir: str, analysis_id: str) -> str:
    return os.path.join(hub_dir, "data", analysis_id + SIDECAR_SUFFIX)


def _mtime_ns(path: str) -> int | None:
    try:
        return os.stat(path).st_mtime_ns
    except OSError:
        return None


def _fresh(store: store_mod.AnalysisStore, hub_dir: str, aid: str) -> bool:
    mine = _mtime_ns(sidecar_path(hub_dir, aid))
    times = [t for t in (_mtime_ns(store.bundle_path(aid)), _mtime_ns(store.status_path(aid))) if t is not None]
    return mine is not None and bool(times) and mine > max(times)


def write_sidecars(library, hub_dir: str, *, sidecars: bool = True, only=None) -> dict:
    """Write the data file of every analysis whose bundle or verdict is newer than it (or that ``only`` names),
    and remove the data files of analyses the library no longer holds. ``{"written": [...], "removed": [...]}``."""
    store = store_mod.AnalysisStore(library.path)
    ids = store.ids()
    known = set(ids)
    data_dir = os.path.join(hub_dir, "data")
    written, removed = [], []
    if os.path.isdir(data_dir):
        for name in sorted(os.listdir(data_dir)):
            key = name[: -len(SIDECAR_SUFFIX)] if name.endswith(SIDECAR_SUFFIX) else ""
            if is_analysis_id(key) and (key not in known or not sidecars):
                os.unlink(os.path.join(data_dir, name))
                removed.append(os.path.join(data_dir, name))
    if sidecars:
        wanted = set(only or ())
        for aid in ids:
            if aid not in wanted and _fresh(store, hub_dir, aid):
                continue
            try:
                data = sidecar(aid, store.composed_gzip(aid))
            except (FlowCellarError, OSError, ValueError):
                continue
            atomic.atomic_write_bytes(sidecar_path(hub_dir, aid), data)
            written.append(sidecar_path(hub_dir, aid))
    return {"written": written, "removed": removed}


def serve_sidecar(library, hub_dir: str, analysis_id: str, headers: dict) -> Served | None:
    """What ``flowcellar serve`` answers for ``/data/<analysis id>.js``; None for an id the library does not hold."""
    if not is_analysis_id(analysis_id):
        return None
    store = store_mod.AnalysisStore(library.path)
    if not store.has(analysis_id):
        return None
    if _fresh(store, hub_dir, analysis_id):
        try:
            return Served(atomic.read_bytes(sidecar_path(hub_dir, analysis_id), retries=1), JS, headers)
        except OSError:
            pass
    try:
        return Served(sidecar(analysis_id, store.composed_gzip(analysis_id)), JS, headers)
    except (FlowCellarError, OSError, ValueError):
        return None
