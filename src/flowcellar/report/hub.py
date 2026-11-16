"""The static hub of a library (mode B).

    hub/index.html            the viewer and no data; rewritten only when the viewer changed
    hub/index.js              the index, as one call of __FC_INDEX__; rewritten every time
    hub/data/<run_key>.js     one run, as one call of __FCB_SIDECAR__; rewritten only when its bundle or its
                              status file is newer than it, or when the caller names the run

Every write is a whole-file replace, so a browser that reads while a scan writes never sees half a file.
"""
from __future__ import annotations

import os

from flowcellar.model.bundle import bundle_filename, report_filename
from flowcellar.model.status import status_filename
from flowcellar.report import html as page
from flowcellar.report import index_one
from flowcellar.util import atomic

SIDECAR_SUFFIX = ".js"


def sidecar_path(hub_dir: str, run_key: str) -> str:
    return os.path.join(hub_dir, "data", run_key + SIDECAR_SUFFIX)


def _mtime_ns(path: str) -> int | None:
    try:
        return os.stat(path).st_mtime_ns
    except OSError:
        return None


def source_mtime_ns(library, run_key: str) -> int | None:
    """When the run's bundle or its status file last changed; None when neither can be found (then nothing can
    be trusted to be up to date)."""
    times = [_mtime_ns(os.path.join(library.path, "bundles", bundle_filename(run_key))),
             _mtime_ns(os.path.join(library.path, "status", status_filename(run_key)))]
    times = [t for t in times if t is not None]
    return max(times) if times else None


def is_fresh(library, hub_dir: str, run_key: str) -> bool:
    """True when the sidecar exists and is newer than the run's bundle and status file."""
    mine, source = _mtime_ns(sidecar_path(hub_dir, run_key)), source_mtime_ns(library, run_key)
    return mine is not None and source is not None and mine > source


def _read(path: str) -> bytes | None:
    try:
        with open(path, "rb") as fh:
            return fh.read()
    except OSError:
        return None


def _head(path: str, n: int = 2048) -> bytes:
    try:
        with open(path, "rb") as fh:
            return fh.read(n)
    except OSError:
        return b""


def _write_if_changed(path: str, data: bytes) -> bool:
    if _read(path) == data:
        return False
    atomic.atomic_write_bytes(path, data)
    return True


def _has_analyses(library_path: str, data_dir: str) -> bool:
    for folder in (os.path.join(library_path, "bundles"), data_dir):
        try:
            if any(n.startswith("an_") for n in os.listdir(folder)):
                return True
        except OSError:
            pass
    return False


def write_hub(library, *, sidecars: bool = True, reports: bool = False, only=None) -> dict:
    hub_dir = library.hub_dir()
    data_dir = os.path.join(hub_dir, "data")
    os.makedirs(hub_dir, exist_ok=True)
    keys = [k for k in library.run_keys() if page.RUN_KEY_RE.match(k)]
    known = set(keys)
    written: list[str] = []
    removed: list[str] = []

    # the page: only when the viewer (or this package's text around it) changed
    index_html = os.path.join(hub_dir, "index.html")
    if _write_if_changed(index_html, page.hub_page().encode("utf-8")):
        written.append(index_html)

    # the index: always
    index_js = os.path.join(hub_dir, "index.js")
    atomic.atomic_write_bytes(index_js, page.index_sidecar(library.index_gzip()))
    written.append(index_js)

    # data files of runs that are no longer in the library
    if os.path.isdir(data_dir):
        for name in sorted(os.listdir(data_dir)):
            if not name.endswith(SIDECAR_SUFFIX) or name.startswith("."):
                continue
            key = name[: -len(SIDECAR_SUFFIX)]
            if key.startswith("an_"):
                continue                            # the data file of an analysis: looked after below
            stale = key not in known or (not sidecars and not is_fresh(library, hub_dir, key))
            if stale:
                os.unlink(os.path.join(data_dir, name))
                removed.append(os.path.join(data_dir, name))

    todo: list[str] = []
    if sidecars:
        wanted = set(only) & known if only is not None else None
        for key in keys:
            if (wanted is not None and key in wanted) or not is_fresh(library, hub_dir, key):
                todo.append(key)
        for key in todo:
            path = sidecar_path(hub_dir, key)
            data = page.run_sidecar(key, library.composed_gzip(key))
            if _read(path) == data and (wanted is None or key not in wanted):
                os.utime(path)                      # same content: only say that it is up to date
                continue
            atomic.atomic_write_bytes(path, data)
            written.append(path)

    # single-cell and spatial analyses keep their data files beside those of the runs; the add-on writes them and
    # is loaded only when the library holds an analysis or the folder still has such a file
    if _has_analyses(library.path, data_dir):
        from flowcellar.sc import report as sc_report
        out = sc_report.write_sidecars(library, hub_dir, sidecars=sidecars, only=only)
        written.extend(out["written"])
        removed.extend(out["removed"])

    if reports:
        index = library.index()
        rep_dir = os.path.join(library.path, "reports")
        named = set(only or ())
        script = page.load_viewer().js_csp.encode("ascii").replace(b"'", b"&#x27;")       # as the policy element spells it
        for key in keys:
            path = os.path.join(rep_dir, report_filename(key))
            source, mine = source_mtime_ns(library, key), _mtime_ns(path)
            # a report is kept when it is newer than the run's data, was made with this viewer and was not asked for
            if mine is not None and source is not None and mine > source and key not in named and key not in todo and script in _head(path):
                continue
            bundle = library.composed(key)
            one = index_one.slice_index(index, key) or index_one.plain_index(bundle)
            text, _ = page.report_page(bundle, one)
            atomic.atomic_write_bytes(path, text.encode("utf-8"))
            written.append(path)
        if os.path.isdir(rep_dir):
            for name in sorted(os.listdir(rep_dir)):
                if name.endswith(report_filename("")) and name[: -len(report_filename(""))] not in known:
                    os.unlink(os.path.join(rep_dir, name))
                    removed.append(os.path.join(rep_dir, name))

    return {"index_html": index_html, "written": written, "removed": removed}
