"""Where everything lives inside a library folder. Paths only: no I/O."""
from __future__ import annotations

import os

from flowcellar.config.loader import CONFIG_NAME
from flowcellar.model.bundle import BUNDLE_SUFFIX, REPORT_SUFFIX
from flowcellar.model.runkey import RUN_KEY_RE
from flowcellar.model.status import STATUS_SUFFIX
from flowcellar.util.errors import InputError

MARKER = "LIBRARY"
MARKER_FORMAT = "flowcellar.library"
CARD_SUFFIX = ".card.json"
STATE_SUFFIX = ".state.json"
TOMBSTONE_SUFFIX = ".tombstone.json"


def check_key(run_key: str) -> str:
    """A run key is the only thing ever joined into a library path: refuse anything else."""
    if not isinstance(run_key, str) or not RUN_KEY_RE.match(run_key):
        raise InputError(f"not a run key: {run_key!r}", hint="a run key looks like TST00001_0a1b2c3d")
    return run_key


class Layout:
    def __init__(self, path: str | os.PathLike):
        self.path = os.path.abspath(os.fspath(path))
        j = lambda *p: os.path.join(self.path, *p)  # noqa: E731
        self.marker = j(MARKER)
        self.config = j(CONFIG_NAME)
        self.bundles, self.status, self.cards, self.state = j("bundles"), j("status"), j("cards"), j("state")
        self.partials, self.probes = j("partials"), j("probes")
        self.index_dir, self.index = j("index"), j("index", "index.json.gz")
        self.hub, self.hub_data, self.reports = j("hub"), j("hub", "data"), j("reports")
        self.claims, self.lease = j("claims"), j("claims", "scan.lease")
        self.trash, self.logs = j("trash"), j("logs")
        self.roots_state, self.notices = j("state", "roots.json"), j("state", "notices.json")
        self.links, self.sniffed = j("state", "links.json"), j("state", "sniffed.json")

    def bundle(self, k: str) -> str:
        return os.path.join(self.bundles, check_key(k) + BUNDLE_SUFFIX)

    def status_file(self, k: str) -> str:
        return os.path.join(self.status, check_key(k) + STATUS_SUFFIX)

    def card(self, k: str) -> str:
        return os.path.join(self.cards, check_key(k) + CARD_SUFFIX)

    def run_state(self, k: str) -> str:
        return os.path.join(self.state, check_key(k) + STATE_SUFFIX)

    def tombstone(self, k: str) -> str:
        return os.path.join(self.state, check_key(k) + TOMBSTONE_SUFFIX)

    def sidecar(self, k: str) -> str:
        return os.path.join(self.hub_data, check_key(k) + ".js")

    def report(self, k: str) -> str:
        return os.path.join(self.reports, check_key(k) + REPORT_SUFFIX)

    def keys_in(self, folder: str, suffix: str) -> list[str]:
        """Run keys that have a file with ``suffix`` in ``folder`` (names that are not run keys are ignored)."""
        try:
            names = os.listdir(folder)
        except FileNotFoundError:
            return []
        out = [n[: -len(suffix)] for n in names if n.endswith(suffix) and not n.startswith(".")]
        return sorted(k for k in out if RUN_KEY_RE.match(k))
