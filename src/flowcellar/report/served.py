"""What ``flowcellar serve`` answers (mode C).

Four addresses and nothing else. No path that came from a request is ever joined to a folder: the page is the
one fixed file ``hub/index.html``, the index is made from the library, and a run's data file is found from a run
key that matched ``RUN_KEY_RE`` and is one of the library's own.
"""
from __future__ import annotations

import html as html_lib
import os
import re

from flowcellar.contracts import Served
from flowcellar.report import html as page
from flowcellar.report import hub as hub_files
from flowcellar.util import atomic
from flowcellar.util.errors import FlowCellarError

HTML = "text/html; charset=utf-8"
JS = "text/javascript; charset=utf-8"
_META_CSP = re.compile(rb'<meta http-equiv="Content-Security-Policy" content="([^"]*)">')


def headers(policy: str | None = None) -> dict[str, str]:
    """On every answer. ``policy`` is the one of the page being served; by default that of the installed viewer."""
    return {"Content-Security-Policy": policy or page.hub_policy(), "X-Content-Type-Options": "nosniff",
            "Referrer-Policy": "no-referrer", "Cross-Origin-Opener-Policy": "same-origin", "Cache-Control": "no-store"}


def page_policy(body: bytes) -> str | None:
    """The policy a hub page states in its own meta element. The header must say the same: a browser applies
    both, and a page written by an older viewer carries the hash of that viewer's script."""
    m = _META_CSP.search(body[:4096])
    return html_lib.unescape(m.group(1).decode("utf-8")) if m else None


def serve_file(library, url_path: str) -> Served | None:
    path = (url_path or "").split("#", 1)[0].split("?", 1)[0]      # the viewer asks for index.js?t=<counter>
    hub_dir = library.hub_dir()
    if path in ("/", "/index.html"):
        try:
            body = atomic.read_bytes(os.path.join(hub_dir, "index.html"), retries=1)
        except OSError:
            return None                                             # no hub was written: serve says how to write one
        return Served(body, HTML, headers(page_policy(body)))
    if path == "/index.js":
        try:
            return Served(page.index_sidecar(library.index_gzip()), JS, headers())
        except (FlowCellarError, OSError):
            return None
    if path.startswith("/data/") and path.endswith(".js"):
        run_key = path[len("/data/"): -len(".js")]
        if run_key.startswith("an_"):               # an analysis of the single-cell and spatial add-on
            from flowcellar.sc import report as sc_report
            return sc_report.serve_sidecar(library, hub_dir, run_key, headers())
        if not page.RUN_KEY_RE.match(run_key) or run_key not in set(library.run_keys()):
            return None
        if hub_files.is_fresh(library, hub_dir, run_key):
            try:
                return Served(atomic.read_bytes(hub_files.sidecar_path(hub_dir, run_key), retries=1), JS, headers())
            except OSError:
                pass
        try:
            return Served(page.run_sidecar(run_key, library.composed_gzip(run_key)), JS, headers())
        except (FlowCellarError, OSError, KeyError):
            return None
    return None
