"""``flowcellar rm``: take runs out of a library.

Everything the library holds for the run (bundle, status file, card, scan state, hub sidecar, report) moves to
``trash/<time>_<run_key>/`` and a tombstone is written, so that the next scan does not add the run again. The
folders that were scanned are never touched. ``flowcellar run RUN_DIR -L DIR`` brings a run back on purpose.
"""
from __future__ import annotations

import os
import shutil
import time

from flowcellar.library.lease import Lease
from flowcellar.library.store import LibraryStore
from flowcellar.util import atomic
from flowcellar.util.errors import InputError
from flowcellar.util.timeutil import now_z


def resolve_run(store: LibraryStore, run_id: str) -> str:
    """A run key, or a flow cell id when the library holds exactly one run of it."""
    keys = store.run_keys()
    if run_id in keys:
        return run_id
    mine = [k for k in keys if k.split("_", 1)[0] == run_id]
    if len(mine) == 1:
        return mine[0]
    if len(mine) > 1:
        raise InputError(f"{run_id}: the library holds {len(mine)} runs of this flow cell: " + ", ".join(mine),
                         hint="give the run key")
    near = [k for k in keys if k.startswith(run_id)]
    raise InputError(f"{run_id}: no such run in this library" + (f" (did you mean {near[0]}?)" if len(near) == 1 else ""),
                     hint="flowcellar ls -L DIR lists the runs")


def remove_runs(store: LibraryStore, run_ids: list[str], *, purge: bool = False, wait: bool = True) -> list[dict]:
    """Remove the runs and rebuild index and hub. Returns ``[{"run_key", "trash" | None, "files"}]``."""
    from flowcellar.library.maintain import rebuild_index_and_hub
    keys = []
    for rid in run_ids:
        k = resolve_run(store, rid)
        if k not in keys:
            keys.append(k)
    lay = store.layout
    out = []
    with Lease(lay.lease, wait=wait):
        stamp = time.strftime("%Y%m%dT%H%M%SZ", time.gmtime())
        for k in keys:
            files = [p for p in (lay.bundle(k), lay.status_file(k), lay.card(k), lay.run_state(k), lay.sidecar(k), lay.report(k))
                     if os.path.isfile(p)]
            dest = None
            if purge:
                for p in files:
                    os.unlink(p)
            else:
                dest = os.path.join(lay.trash, f"{stamp}_{k}")
                os.makedirs(dest, exist_ok=True)
                for p in files:
                    shutil.move(p, os.path.join(dest, os.path.basename(p)))
            atomic.atomic_write_json(lay.tombstone(k), {"run_key": k, "removed": now_z()}, indent=1)
            out.append({"run_key": k, "trash": dest, "files": len(files)})
        rebuild_index_and_hub(store)
    return out
