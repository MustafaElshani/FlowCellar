"""What the other parts call in VIEWER's Python. Signatures are the contract."""
from __future__ import annotations

import logging
import os
from typing import Sequence

from flowcellar.contracts import Served
from flowcellar.util import atomic
from flowcellar.util import log as _log

WARN_REPORT_BYTES = 12 << 20        # a report above 12 MB opens slowly in WebKit


def write_report(bundle: dict, out_path: str, *, index: dict | None = None) -> str:
    """Write one offline HTML report for one run. ``bundle`` is composed (it carries ``status``). ``index`` is a
    one-run index; when None it is built with ``flowcellar.library.api.build_index``. Returns ``out_path``."""
    from flowcellar.report import html as page
    from flowcellar.report import index_one

    out_path = os.fspath(out_path)
    text, embedded = page.report_page(bundle, index if index is not None else index_one.index_for_bundle(bundle))
    data = text.encode("utf-8")
    if len(data) > WARN_REPORT_BYTES:
        _log.event("report.large", f"{os.path.basename(out_path)} is {len(data) / 1e6:.0f} MB "
                   f"({embedded / 1e6:.0f} MB of it data): above 12 MB a report opens slowly in Safari",
                   level=logging.WARNING, path=out_path, bytes=len(data))
    atomic.atomic_write_bytes(out_path, data)
    return out_path


def write_hub(library, *, sidecars: bool = True, reports: bool = False, only: Sequence[str] | None = None) -> dict:
    """Write the static hub of a library: ``hub/index.html`` (when the viewer changed), ``hub/index.js`` (always) and
    ``hub/data/<run_key>.js`` for the run keys in ``only`` (None = every run whose sidecar is missing or older than
    its bundle or status file). Returns ``{"index_html": path, "written": [...], "removed": [...]}``."""
    from flowcellar.report import hub

    return hub.write_hub(library, sidecars=sidecars, reports=reports, only=only)


def serve_file(library, url_path: str) -> Served | None:
    """What ``flowcellar serve`` answers for one GET path: ``/``, ``/index.html``, ``/index.js``, ``/data/<run_key>.js``.
    None for anything else (the server answers 404). Carries the Content-Security-Policy header."""
    from flowcellar.report import served

    return served.serve_file(library, url_path)


def viewer_manifest() -> dict:
    """``viewer_dist/manifest.json``: file names, sizes and SHA-256 hashes of the built viewer. Raises InputError with
    one clear line when the viewer has not been built."""
    from flowcellar.report import html as page

    return page.viewer_manifest()
