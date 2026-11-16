"""Pages made from the built viewer.

A page is ``viewer_dist/template.html`` with eight placeholders filled, each exactly once and in one pass, so
that no text put into the page can be read as a placeholder. No template engine and no Node at run time.

What may reach a page from a run folder (names typed at the sequencer, sample names, paths) arrives in three
ways only, and each is closed:

* the title and the text summary: every string goes through ``html.escape``;
* the boot block: JSON with ``<``, ``>`` and ``&`` written as escapes, inside a script element that is not run;
* the data: base64 of gzip bytes (alphabet asserted) inside script elements of a non-JavaScript type, which the
  viewer decodes and draws as text nodes (it never builds markup from data).
"""
from __future__ import annotations

import base64
import hashlib
import html
import json
import re
from functools import lru_cache
from importlib import resources

from flowcellar import SCHEMA_VERSION, __version__
from flowcellar.model import status as status_model
from flowcellar.util import atomic
from flowcellar.util.errors import InputError

PLACEHOLDERS = ("%%CSP%%", "%%TITLE%%", "%%FAVICON%%", "%%CSS%%", "%%SUMMARY%%", "%%BOOT%%", "%%PAYLOAD%%", "%%JS%%")
VIEWER_FILES = ("viewer.core.js", "viewer.css", "favicon.svg", "template.html")
#: a run key as it may appear in a file name, an element id and an address (the viewer checks the same pattern)
RUN_KEY_RE = re.compile(r"^[A-Z]{2,3}[0-9]{4,6}_[0-9a-f]{8,12}$")
_B64_RE = re.compile(r"^[A-Za-z0-9+/]*={0,2}$")
_PLACEHOLDER_RE = re.compile(r"(%%[A-Z]+%%)")
INDEX_TYPE = "application/vnd.flowcellar.index+gzip"
BUNDLE_TYPE = "application/vnd.flowcellar.bundle+gzip"
NOT_BUILT = ("the viewer is not built: src/flowcellar/viewer_dist/ is missing or incomplete",
             "developers: run 'node viewer/build.mjs' in the repository; everybody else: install a released wheel")


# ---------------------------------------------------------------------------------------------------------
# the built viewer
# ---------------------------------------------------------------------------------------------------------

def _dist():
    return resources.files("flowcellar").joinpath("viewer_dist")


def viewer_manifest() -> dict:
    """``viewer_dist/manifest.json``. InputError with one clear line when the viewer was not built."""
    try:
        m = json.loads(_dist().joinpath("manifest.json").read_text("utf-8"))
    except (FileNotFoundError, NotADirectoryError, json.JSONDecodeError):
        raise InputError(NOT_BUILT[0], hint=NOT_BUILT[1]) from None
    if not isinstance(m, dict) or m.get("format") != "flowcellar.viewer" or not all(f in m.get("files", {}) for f in VIEWER_FILES):
        raise InputError(NOT_BUILT[0], hint=NOT_BUILT[1])
    return m


class Viewer:
    """The four files of the built viewer, read once and checked against the manifest."""

    def __init__(self, manifest: dict, files: dict[str, str]):
        self.manifest = manifest
        self.js = files["viewer.core.js"]
        self.css = files["viewer.css"]
        self.template = files["template.html"]
        self.favicon_uri = "data:image/svg+xml;base64," + base64.b64encode(files["favicon.svg"].encode("utf-8")).decode("ascii")
        self.js_csp = csp_hash(self.js)
        self.css_csp = csp_hash(self.css)
        #: changes whenever any file of the viewer changes
        self.build_id = hashlib.sha256("".join(manifest["files"][f]["sha256"] for f in VIEWER_FILES).encode()).hexdigest()[:16]


def csp_hash(text: str) -> str:
    """The Content-Security-Policy source of an inline script or style: ``'sha256-<base64>'`` of its exact text."""
    return "'sha256-" + base64.b64encode(hashlib.sha256(text.encode("utf-8")).digest()).decode("ascii") + "'"


@lru_cache(maxsize=2)
def _load(manifest_text: str) -> Viewer:
    manifest = json.loads(manifest_text)
    files = {}
    for name in VIEWER_FILES:
        try:
            text = _dist().joinpath(name).read_bytes().decode("utf-8")
        except FileNotFoundError:
            raise InputError(NOT_BUILT[0], hint=NOT_BUILT[1]) from None
        if hashlib.sha256(text.encode("utf-8")).hexdigest() != manifest["files"][name]["sha256"]:
            raise InputError(f"the built viewer is damaged: {name} does not match viewer_dist/manifest.json", hint=NOT_BUILT[1])
        files[name] = text
    found = _PLACEHOLDER_RE.findall(files["template.html"])
    if sorted(found) != sorted(PLACEHOLDERS):
        raise InputError("the built viewer is damaged: template.html must hold each placeholder exactly once", hint=NOT_BUILT[1])
    if re.search(r"</script", files["viewer.core.js"], re.I) or re.search(r"</style", files["viewer.css"], re.I):
        raise InputError("the built viewer is damaged: the script or the stylesheet would end its own element", hint=NOT_BUILT[1])
    return Viewer(manifest, files)


def load_viewer() -> Viewer:
    viewer_manifest()                                                # the clear error when nothing is built
    return _load(_dist().joinpath("manifest.json").read_text("utf-8"))


# ---------------------------------------------------------------------------------------------------------
# the parts of a page
# ---------------------------------------------------------------------------------------------------------

def policy(viewer: Viewer, mode: str) -> str:
    """The Content-Security-Policy of a page. A hub may load its sibling data files ('self'); a report loads nothing."""
    script = viewer.js_csp + (" 'self'" if mode == "hub" else "")
    return (f"default-src 'none'; script-src {script}; style-src {viewer.css_csp}; img-src data: blob:; "
            "font-src 'none'; connect-src 'none'; base-uri 'none'; form-action 'none'")


def boot_json(mode: str, run: str | None = None) -> str:
    """The boot block. It sits inside a script element, so ``<``, ``>``, ``&`` and the two line separators that
    JSON allows raw are written as escapes."""
    if run is not None and not RUN_KEY_RE.match(run):
        raise InputError(f"not a run key: {run!r}")
    text = json.dumps({"mode": mode, "generator": __version__, "schema": SCHEMA_VERSION, "run": run},
                      ensure_ascii=True, separators=(",", ":"))
    return text.replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026")


def b64(data: bytes) -> str:
    text = base64.b64encode(data).decode("ascii")
    assert _B64_RE.match(text), "payload alphabet"
    return text


def data_element(mime: str, element_id: str, gz: bytes, *, sha256: bool = False) -> str:
    """One embedded payload: base64 of gzip bytes in a script element the browser does not run."""
    if not re.fullmatch(r"[A-Za-z0-9:_]+", element_id):
        raise InputError(f"not usable as an element id: {element_id!r}")
    payload = b64(gz)
    sha = f' data-sha256="{hashlib.sha256(gz).hexdigest()}"' if sha256 else ""
    return f'<script type="{mime}" id="{element_id}"{sha}>{payload}</script>\n'


def fill(viewer: Viewer, *, mode: str, title: str, summary: str, boot: str, payload: str = "") -> str:
    """The page text. ``title`` is escaped here; ``summary`` must already be markup built with ``esc``."""
    if "<" in boot:
        raise InputError("the boot block must not contain '<'")
    values = {"%%CSP%%": html.escape(policy(viewer, mode), quote=True), "%%TITLE%%": html.escape(title, quote=True),
              "%%FAVICON%%": viewer.favicon_uri, "%%CSS%%": viewer.css, "%%SUMMARY%%": summary, "%%BOOT%%": boot,
              "%%PAYLOAD%%": payload, "%%JS%%": viewer.js}
    # one pass over the template: a placeholder name that turns up inside a value is never looked at again
    parts = _PLACEHOLDER_RE.split(viewer.template)
    seen = [p for p in parts if p in values]
    if sorted(seen) != sorted(PLACEHOLDERS):
        raise InputError("template.html must hold each placeholder exactly once")
    return "".join(values.get(p, p) for p in parts)


# ---------------------------------------------------------------------------------------------------------
# the text summary: what a mail preview, a reader without scripts and an old browser show
# ---------------------------------------------------------------------------------------------------------

def esc(text) -> str:
    return html.escape("" if text is None else str(text), quote=True)


def _fmt_count(x) -> str:
    if x is None:
        return "–"
    if abs(x) >= 1e6:
        v = x / 1e6
        return (f"{v:.0f}" if abs(v) >= 100 else f"{v:.1f}") + " M"
    return f"{int(round(x)):,}"


def _fmt_bases(x) -> str:
    if x is None:
        return "–"
    for div, unit in ((1e9, "Gb"), (1e6, "Mb"), (1e3, "kb")):
        if abs(x) >= div:
            v = x / div
            return (f"{v:.0f}" if abs(v) >= 100 else f"{v:.1f}") + " " + unit
    return f"{int(round(x))} b"


def _fmt_hours(h) -> str:
    if h is None:
        return "–"
    if h < 1:
        return f"{max(1, round(h * 60))} min"
    return f"{h:.0f} h" if abs(h - round(h)) < 0.05 else f"{h:.1f} h"


def _badge(verdict: str, *, large: bool = False, so_far: bool = False, word: str | None = None) -> str:
    cls = {"pass": "pass", "warn": "warn", "fail": "fail", "incomplete": "incomplete", "running": "running"}.get(verdict, "unknown")
    text = word or status_model.VERDICT_WORD.get(verdict, verdict)
    if so_far and verdict in ("pass", "warn"):
        text += " so far"
    return f'<span class="fc-badge fc-badge--{cls}{" fc-badge--lg" if large else ""}">{esc(text)}</span>'


def _summaries(bundle: dict) -> dict:
    out: dict = {}
    for m in (bundle.get("modules") or {}).values():
        out.update(m.get("summary") or {})
    return out


def key_numbers(bundle: dict) -> list[tuple[str, str]]:
    """Label and value of the numbers a reader looks for first, from the summaries the bundle stores."""
    s = _summaries(bundle)
    rows: list[tuple[str, str]] = []
    if s.get("yield_reads") is not None:
        rows.append(("Reads", _fmt_count(s["yield_reads"])))
    if s.get("yield_estimated_bases") is not None:
        rows.append(("Estimated bases, from signal", _fmt_bases(s["yield_estimated_bases"])))
    if s.get("run_duration_h") is not None:
        rows.append(("Run time", _fmt_hours(s["run_duration_h"])))
    if s.get("fc_pores_first_scan") is not None:
        rows.append(("Pores at the first scan", _fmt_count(s["fc_pores_first_scan"])))
    if s.get("fc_occupancy_first_hour") is not None and (s.get("run_duration_h") or 0) >= 1.5:
        rows.append(("Occupancy, first hour", f"{s['fc_occupancy_first_hour'] * 100:.0f} %"))
    if s.get("yield_basecalled_bases") is not None:
        rows.append(("Passed bases", _fmt_bases(s["yield_basecalled_bases"])))
    if s.get("len_n50") is not None:
        rows.append(("Read N50", _fmt_bases(s["len_n50"])))
    if s.get("qual_median") is not None:
        rows.append(("Median read quality", f"Q {s['qual_median']:.1f}"))
    return rows


def report_title(bundle: dict) -> str:
    st = bundle.get("status") or {}
    fc = ((bundle.get("identity") or {}).get("flow_cell") or {}).get("id") or bundle.get("bundle_id") or "run"
    word = status_model.VERDICT_WORD.get(st.get("verdict"), "")
    if st.get("provisional") and st.get("verdict") in ("pass", "warn"):
        word += " so far"
    return " · ".join(x for x in (fc, word, "FlowCellar run report") if x)


def report_summary(bundle: dict) -> str:
    """Verdict, sentences, the two parts and the key numbers as escaped text with the badge and tile classes."""
    ident, st = bundle.get("identity") or {}, bundle.get("status") or {}
    fc = (ident.get("flow_cell") or {}).get("id") or bundle.get("bundle_id")
    name = (ident.get("protocol_run") or {}).get("sample_id")
    started = ((ident.get("protocol_run") or {}).get("started") or "")[:16].replace("T", " ")
    by = st.get("by_tier") or {}
    rules = st.get("rules") or []
    out = ['<div class="fc-page fc-nojs">', '<p class="fc-caps">FlowCellar run report</p>',
           f'<h1 class="fc-runhead__id fc-mono">{esc(fc)}</h1>']
    # MinKNOW's files give UTC; the start of a folder with raw data only is the time in its name: the sequencer's own clock
    raw_only = (ident.get("protocol_run") or {}).get("run_state") == "raw_only"
    when = (f"started {started} (from the folder name, the sequencer's clock)" if raw_only else f"started {started} UTC") if started else None
    line = " · ".join(esc(x) for x in (name, when) if x)
    if line:
        out.append(f'<p class="fc-runhead__line">{line}</p>')
    out.append('<section class="fc-verdict"><div class="fc-verdict__head">'
               + _badge(st.get("verdict") or "unknown", large=True, so_far=bool(st.get("provisional")))
               + f'<div><p class="fc-verdict__headline">{esc(st.get("headline"))}</p>'
               + (f'<p class="fc-verdict__detail">{esc(st.get("detail"))}</p>' if st.get("detail") else "") + "</div></div>")
    for m in st.get("modifiers") or []:
        out.append(f'<p class="fc-verdict__note">{esc(m.get("message"))}</p>')
    out.append('<div class="fc-verdict__parts">')
    for tier, title in (("run_health", "Run health"), ("read_level", "Read data")):
        v = by.get(tier) or "not_available"
        mine = [r for r in rules if r.get("tier") == tier]
        trig = [r for r in mine if r.get("outcome") == "triggered" and r.get("severity") in ("warn", "fail")]
        ok = sum(1 for r in mine if r.get("outcome") == "ok")
        skipped = sum(1 for r in mine if r.get("outcome") == "not_evaluated")
        if not mine:
            count = ""
        elif not trig and not skipped:
            count = f"{ok} of {len(mine)} checks passed"
        else:
            count = " · ".join(x for x in (f"{ok} passed", f"{len(trig)} triggered" if trig else "", f"{skipped} not assessed" if skipped else "") if x)
        pending = len(st.get("pending") or []) if tier == "run_health" else 0
        if pending:
            count += f" · {pending} wait for the read count"
        out.append('<div class="fc-verdict__part"><div class="fc-verdict__parthead">'
                   + _badge(v, so_far=tier == "run_health" and bool(st.get("provisional")))
                   + f"<b>{title}</b>" + (f'<span class="fc-secondary">{esc(count)}</span>' if count else "") + "</div>")
        if trig:
            out.append('<ul class="fc-reasons">' + "".join(f"<li><span></span><span>{esc(r.get('message'))}</span></li>" for r in trig[:6]) + "</ul>")
        out.append("</div>")
    out.append("</div></section>")
    nums = key_numbers(bundle)
    if nums:
        out.append('<div class="fc-keynumbers"><div class="fc-tiles">' + "".join(
            f'<div class="fc-tile"><span class="fc-tile__label"><span>{esc(k)}</span></span><span class="fc-tile__value">{esc(v)}</span></div>'
            for k, v in nums) + "</div></div>")
    out.append('<p class="fc-nojs__note fc-secondary">This is the text summary of the report. The charts and every other number appear '
               "when the file is opened in a browser with scripts switched on (Firefox 113, Chrome or Edge 80, Safari 16.4, or newer). "
               f"Self-contained: works offline, sends nothing. Made by FlowCellar {esc(__version__)}.</p>")
    out.append("</div>")
    return "".join(out)


HUB_TITLE = "FlowCellar library"
HUB_SUMMARY = ('<div class="fc-page fc-nojs"><p class="fc-caps">FlowCellar</p><h1 class="fc-runhead__id">Library</h1>'
               '<p class="fc-nojs__note">This page shows the runs of a FlowCellar library. It needs a browser with scripts switched on '
               "(Firefox 113, Chrome or Edge 80, Safari 16.4, or newer); it reads only the files beside it and sends nothing.</p></div>")


# ---------------------------------------------------------------------------------------------------------
# whole pages
# ---------------------------------------------------------------------------------------------------------

def hub_page() -> str:
    """``hub/index.html``: the viewer and no data. It changes only when the viewer changes."""
    v = load_viewer()
    return fill(v, mode="hub", title=HUB_TITLE, summary=HUB_SUMMARY, boot=boot_json("hub"))


def hub_policy() -> str:
    return policy(load_viewer(), "hub")


def report_page(bundle: dict, index: dict) -> tuple[str, int]:
    """The one-run report and the size of what it embeds. ``bundle`` is composed; ``index`` holds its one row."""
    from flowcellar.model import bundle as bundle_model

    run_key = bundle.get("bundle_id") or (bundle.get("identity") or {}).get("run_key")
    if not run_key or not RUN_KEY_RE.match(run_key):
        raise InputError(f"cannot write a report: {run_key!r} is not a run key")
    if "status" not in bundle:
        raise InputError(f"cannot write a report for {run_key}: the bundle carries no verdict (compose it with its status first)")
    v = load_viewer()
    gz_index = atomic.gzip_bytes(atomic.dumps(index).encode("utf-8"))
    gz_bundle = bundle_model.to_gzip_bytes(bundle)
    payload = data_element(INDEX_TYPE, "fci", gz_index) + data_element(BUNDLE_TYPE, "fcb:" + run_key, gz_bundle, sha256=True)
    text = fill(v, mode="report", title=report_title(bundle), summary=report_summary(bundle),
                boot=boot_json("report", run_key), payload=payload)
    return text, len(gz_index) + len(gz_bundle)


def index_sidecar(gz_index: bytes) -> bytes:
    """``hub/index.js``."""
    return ('__FC_INDEX__("' + b64(gz_index) + '");\n').encode("ascii")


def run_sidecar(run_key: str, gz_bundle: bytes) -> bytes:
    """``hub/data/<run_key>.js``."""
    if not RUN_KEY_RE.match(run_key or ""):
        raise InputError(f"not a run key: {run_key!r}")
    return ('__FCB_SIDECAR__("' + run_key + '","' + b64(gz_bundle) + '","' + hashlib.sha256(gz_bundle).hexdigest() + '");\n').encode("ascii")
