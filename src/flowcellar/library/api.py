"""What the other parts call in LIBRARY. Signatures are the contract."""
from __future__ import annotations

import os
from typing import Protocol, Sequence

from flowcellar.contracts import Card


class Library(Protocol):
    """An open library. Read access for VIEWER's Python and for commands; only LIBRARY code writes."""
    path: str               # absolute path of the library folder
    id: str                 # library_id of the LIBRARY marker
    name: str
    config: dict            # the parsed flowcellar.toml with defaults filled in

    def run_keys(self) -> list[str]: ...
    def read_measured(self, run_key: str) -> dict: ...
    def read_status_file(self, run_key: str) -> dict | None: ...
    def composed(self, run_key: str) -> dict: ...           # measured + status + user + relations
    def composed_gzip(self, run_key: str) -> bytes: ...     # the bytes a sidecar or a report carries
    def card(self, run_key: str) -> Card | None: ...
    def index(self) -> dict: ...                            # index/index.json.gz, parsed
    def index_gzip(self) -> bytes: ...                      # its bytes
    def hub_dir(self) -> str: ...


def open_library(path: str, *, create: bool = False, name: str | None = None, roots: Sequence[str] = ()) -> Library:
    """Open a library folder. ``create=True`` makes it when it does not exist (LIBRARY marker, flowcellar.toml with
    ``roots``). Raises LibraryError (exit code 6) for a folder that is not a library or has a newer layout."""
    from flowcellar.config import loader
    from flowcellar.library.store import LibraryStore, is_library
    path = os.path.abspath(os.path.expanduser(os.fspath(path)))
    if create and not (os.path.isdir(path) and is_library(path)):
        paths = [loader.normalise_root_path(p) for p in roots]
        ids = loader.root_ids_for(paths, set())
        entries = [{"id": i, "path": p, "live": False, "project_depth": 1, "project_pattern": None} for i, p in zip(ids, paths)]
        return LibraryStore.create(path, name=name, roots=entries)
    return LibraryStore.open(path)


def card_for(bundle: dict, status_file: dict | None = None, *, extra: dict | None = None) -> Card:
    """The complete index row of one run: ENGINE's columns, READ DATA's columns and the library's own."""
    from flowcellar.library.cards import card_for as _card_for
    return _card_for(bundle, status_file, extra=extra)


def build_index(cards: Sequence[Card], *, library: dict | None = None, roots: Sequence[dict] = (),
                samples: Sequence[dict] = (), notices: Sequence[dict] = (), progress: dict | None = None,
                baselines: dict | None = None) -> dict:
    """The index payload from cards. Pure: needs no library folder, so a one-run report
    builds its index with ``build_index([card_for(bundle, status_file)])``."""
    from flowcellar.library.index import build_index as _build_index
    return _build_index(cards, library=library, roots=roots, samples=samples, notices=notices, progress=progress,
                        baselines=baselines)


def add_run(library: Library, run_dir: str, *, quick: bool = False) -> str:
    """``flowcellar run RUN_DIR -L DIR``: add or refresh one run under the scan lease; returns its run key."""
    from flowcellar.library.scan import Scan, ScanOptions
    return Scan(library.path, ScanOptions()).add_single(run_dir, quick=quick)
