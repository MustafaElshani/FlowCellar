"""A read-only stand-in for ``flowcellar.library.api.Library`` over a folder in the layout of a library.

The viewer's Python talks to a library only through that protocol, so its tests need nothing of the scan: this
class reads bundles/, status/ and index/index.json.gz and composes on request, as the real library does.
"""
from __future__ import annotations

import os

from flowcellar.model import bundle as bundle_model
from flowcellar.model import status as status_model
from flowcellar.util import atomic


class FakeLibrary:
    def __init__(self, path: str):
        self.path = os.path.abspath(path)
        marker = atomic.read_json(os.path.join(self.path, "LIBRARY"))
        self.id = marker["library_id"]
        self.name = os.path.basename(self.path)
        self.config: dict = {}

    def run_keys(self) -> list[str]:
        d = os.path.join(self.path, "bundles")
        return sorted(f[: -len(bundle_model.BUNDLE_SUFFIX)] for f in os.listdir(d) if f.endswith(bundle_model.BUNDLE_SUFFIX))

    def read_measured(self, run_key: str) -> dict:
        return bundle_model.read_bundle(os.path.join(self.path, "bundles", bundle_model.bundle_filename(run_key)))

    def read_status_file(self, run_key: str) -> dict | None:
        p = os.path.join(self.path, "status", status_model.status_filename(run_key))
        return status_model.read_status(p) if os.path.exists(p) else None

    def composed(self, run_key: str) -> dict:
        sf = self.read_status_file(run_key) or {}
        return bundle_model.compose(self.read_measured(run_key), status=sf.get("status"), relations=sf.get("relations"))

    def composed_gzip(self, run_key: str) -> bytes:
        return bundle_model.to_gzip_bytes(self.composed(run_key))

    def card(self, run_key: str):
        return None

    def index(self) -> dict:
        return atomic.read_json(os.path.join(self.path, "index", "index.json.gz"))

    def index_gzip(self) -> bytes:
        return atomic.read_bytes(os.path.join(self.path, "index", "index.json.gz"))

    def hub_dir(self) -> str:
        return os.path.join(self.path, "hub")
