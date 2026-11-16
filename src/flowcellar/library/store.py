"""An open library: reading for everybody, writing for LIBRARY code only.

Every write is a whole-file replace (``util.atomic``). Bundles hold measured content only; the status file, the
card and the scan state of a run are separate small files, so judging a run again never rewrites its bundle.
"""
from __future__ import annotations

import os
import secrets

from flowcellar import LAYOUT_VERSION, __version__
from flowcellar.config import loader
from flowcellar.library.layout import (
    CARD_SUFFIX,
    MARKER,
    MARKER_FORMAT,
    STATE_SUFFIX,
    TOMBSTONE_SUFFIX,
    Layout,
    check_key,
)
from flowcellar.model import bundle as bundle_model
from flowcellar.model import status as status_model
from flowcellar.util import atomic
from flowcellar.util.errors import LibraryError, ResourceLimitError
from flowcellar.util.timeutil import now_z

BUNDLE_WARN_BYTES = 1 << 20        # a bundle above 1 MB gzip is reported
BUNDLE_MAX_BYTES = 5 << 20         # and above 5 MB it is refused


def is_library(path: str | os.PathLike) -> bool:
    return os.path.isfile(os.path.join(os.fspath(path), MARKER))


class LibraryStore:
    """Implements ``flowcellar.library.api.Library``."""

    def __init__(self, path: str, marker: dict, loaded: loader.Loaded):
        self.layout = Layout(path)
        self.path = self.layout.path
        self.marker = marker
        self.id: str = marker["library_id"]
        self.created: str | None = marker.get("created")
        self.loaded = loaded
        self.config: dict = loaded.config
        self.name: str = self.config["library"].get("name") or os.path.basename(self.path.rstrip("/")) or "library"

    # ------------------------------------------------------------------------------------------ opening
    @classmethod
    def open(cls, path: str | os.PathLike, *, config_file: str | None = None) -> LibraryStore:
        path = os.path.abspath(os.fspath(path))
        marker_path = os.path.join(path, MARKER)
        if not os.path.isdir(path):
            raise LibraryError(f"{path}: no such library folder", hint="flowcellar scan ROOT... -L DIR creates one")
        if not os.path.isfile(marker_path):
            raise LibraryError(f"{path}: not a FlowCellar library (it has no {MARKER} file)",
                               hint="flowcellar scan ROOT... -L DIR creates a library in a new or empty folder")
        try:
            marker = atomic.read_json(marker_path)
        except (OSError, ValueError) as e:
            raise LibraryError(f"{marker_path}: cannot be read ({e})") from None
        if not isinstance(marker, dict) or marker.get("format") != MARKER_FORMAT or "library_id" not in marker:
            raise LibraryError(f"{marker_path}: not the marker of a FlowCellar library")
        if int(marker.get("layout_version", 0)) > LAYOUT_VERSION:
            raise LibraryError(f"{path}: made by a newer FlowCellar (layout {marker.get('layout_version')}; this version "
                               f"reads layout {LAYOUT_VERSION})", hint="install a newer FlowCellar")
        loaded = loader.load(config_file or os.path.join(path, loader.CONFIG_NAME))
        return cls(path, marker, loaded)

    @classmethod
    def create(cls, path: str | os.PathLike, *, name: str | None = None, roots: list[dict] | None = None,
               config_text: str | None = None, config_file: str | None = None) -> LibraryStore:
        """Make a library in a new or empty folder: the marker and ``flowcellar.toml``."""
        path = os.path.abspath(os.fspath(path))
        if os.path.exists(path):
            if not os.path.isdir(path):
                raise LibraryError(f"{path}: exists and is not a folder")
            if is_library(path):
                raise LibraryError(f"{path}: is a library already")
            if [n for n in os.listdir(path) if n != loader.CONFIG_NAME]:
                raise LibraryError(f"{path}: exists and is not a library",
                                   hint="choose a new or empty folder; FlowCellar never turns a folder with other files into a library")
        lay = Layout(path)
        try:
            os.makedirs(path, exist_ok=True)
            if not os.path.isfile(lay.config):
                text = config_text if config_text is not None else loader.render_new(
                    name or os.path.basename(path.rstrip("/")) or "library", roots or [])
                atomic.atomic_write_text(lay.config, text)
            marker = {"format": MARKER_FORMAT, "layout_version": LAYOUT_VERSION,
                      "library_id": "lib_" + secrets.token_hex(6), "created": now_z(), "generator": __version__}
            atomic.atomic_write_json(lay.marker, marker, indent=1)
        except OSError as e:
            raise LibraryError(f"{path}: cannot create the library ({e.strerror or e})") from None
        return cls.open(path, config_file=config_file)

    def reload_config(self, config_file: str | None = None) -> None:
        self.loaded = loader.load(config_file or self.layout.config)
        self.config = self.loaded.config
        self.name = self.config["library"].get("name") or os.path.basename(self.path.rstrip("/")) or "library"

    # ------------------------------------------------------------------------------------------ reading
    def run_keys(self) -> list[str]:
        return self.layout.keys_in(self.layout.bundles, bundle_model.BUNDLE_SUFFIX)

    def has_run(self, run_key: str) -> bool:
        return os.path.isfile(self.layout.bundle(run_key))

    def read_measured(self, run_key: str) -> dict:
        path = self.layout.bundle(run_key)
        if not os.path.isfile(path):
            raise LibraryError(f"{run_key}: not in this library")
        return bundle_model.read_bundle(path)

    def read_status_file(self, run_key: str) -> dict | None:
        path = self.layout.status_file(run_key)
        return status_model.read_status(path) if os.path.isfile(path) else None

    def composed(self, run_key: str) -> dict:
        measured = self.read_measured(run_key)
        sf = self.read_status_file(run_key) or {}
        return bundle_model.compose(measured, status=sf.get("status"), relations=sf.get("relations"))

    def composed_gzip(self, run_key: str) -> bytes:
        return bundle_model.to_gzip_bytes(self.composed(run_key))

    def card(self, run_key: str) -> dict | None:
        path = self.layout.card(run_key)
        return atomic.read_json(path) if os.path.isfile(path) else None

    def cards(self) -> dict[str, dict]:
        out = {}
        for k in self.layout.keys_in(self.layout.cards, CARD_SUFFIX):
            c = self.card(k)
            if isinstance(c, dict) and "columns" in c:
                out[k] = c
        return out

    def index(self) -> dict:
        if os.path.isfile(self.layout.index):
            return atomic.read_json(self.layout.index)
        from flowcellar.library.index import index_of_library
        return index_of_library(self)

    def index_gzip(self) -> bytes:
        if os.path.isfile(self.layout.index):
            return atomic.read_bytes(self.layout.index)
        return atomic.gzip_bytes(atomic.dumps(self.index()).encode("utf-8"))

    def hub_dir(self) -> str:
        return self.layout.hub

    def read_state(self, run_key: str) -> dict | None:
        path = self.layout.run_state(run_key)
        return atomic.read_json(path) if os.path.isfile(path) else None

    def states(self) -> dict[str, dict]:
        out = {}
        for k in self.layout.keys_in(self.layout.state, STATE_SUFFIX):
            s = self.read_state(k)
            if isinstance(s, dict):
                out[k] = s
        return out

    def tombstones(self) -> set[str]:
        return set(self.layout.keys_in(self.layout.state, TOMBSTONE_SUFFIX))

    def read_aux(self, path: str, default):
        """One of the small bookkeeping files under state/ (roots, notices, links); ``default`` when absent or broken."""
        try:
            return atomic.read_json(path) if os.path.isfile(path) else default
        except (OSError, ValueError):
            return default

    # ------------------------------------------------------------------------------------------ writing
    def _write(self, fn, *args, **kwargs):
        try:
            return fn(*args, **kwargs)
        except PermissionError as e:
            raise LibraryError(f"{self.path}: not writable ({e.filename or e})") from None
        except OSError as e:
            if e.errno in (28, 30, 122):          # no space, read-only file system, quota
                raise LibraryError(f"{self.path}: cannot write ({e.strerror})") from None
            raise

    def write_measured(self, run_key: str, bundle: dict) -> dict:
        """Write a measured bundle (no status, user or relations). Returns ``{"content_sha256", "bytes",
        "gzip_bytes", "warn"}``. Refuses a bundle above 5 MB gzip."""
        check_key(run_key)
        measured, _ = bundle_model.split(bundle)
        path = self.layout.bundle(run_key)
        prov = measured.setdefault("provenance", {})
        if os.path.isfile(path):
            try:
                old = bundle_model.read_bundle(path).get("provenance") or {}
                prov["created"] = old.get("created") or prov.get("created")
                prov["revision"] = int(old.get("revision") or 1) + 1
            except Exception:  # noqa: BLE001  a damaged old file is simply replaced
                prov.setdefault("revision", 1)
        else:
            prov.setdefault("revision", 1)
        measured["content_sha256"] = bundle_model.content_sha256(measured)
        raw = atomic.dumps(bundle_model.order_members(measured)).encode("utf-8")
        gz = atomic.gzip_bytes(raw)
        if len(gz) > BUNDLE_MAX_BYTES:
            raise ResourceLimitError(f"{run_key}: the bundle would be {len(gz) / 1e6:.1f} MB compressed; the limit is 5 MB")
        self._write(atomic.atomic_write_bytes, path, gz)
        return {"content_sha256": measured["content_sha256"], "bytes": len(raw), "gzip_bytes": len(gz),
                "warn": len(gz) > BUNDLE_WARN_BYTES}

    def write_status_file(self, run_key: str, status_file: dict) -> None:
        self._write(status_model.write_status, self.layout.status_file(run_key), status_file)

    def write_card(self, run_key: str, card: dict) -> None:
        self._write(atomic.atomic_write_json, self.layout.card(run_key), card)

    def write_state(self, run_key: str, state: dict) -> None:
        self._write(atomic.atomic_write_json, self.layout.run_state(run_key), state, indent=1)

    def write_aux(self, path: str, obj) -> None:
        self._write(atomic.atomic_write_json, path, obj, indent=1)

    def write_index(self, index: dict) -> None:
        self._write(atomic.atomic_write_json, self.layout.index, index)
