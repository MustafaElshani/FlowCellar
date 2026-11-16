"""Names and small helpers the tier 1 section builders share."""
from __future__ import annotations

import os

from flowcellar.util import fingerprint as _fp
from flowcellar.util.timeutil import from_epoch

#: versions of the section builders; a bump changes the id of the partials they read
VERSIONS = {"run": "0.1.0", "fc": "0.1.0", "yield": "0.1.0", "dyn": "0.1.0", "len.est": "0.1.0"}

END_REASONS = ("signal_positive", "signal_negative", "unblock_mux_change", "mux_change", "data_service_unblock_mux_change",
               "analysis_config_change", "device_data_error", "api_request", "paused", "unknown", "other")
CHANNEL_STATES = ("strand", "adapter", "pore", "unavailable", "unblocking", "no_pore", "zero", "saturated", "multiple",
                  "unknown_positive", "unknown_negative", "disabled", "locked", "membrane", "unclassified",
                  "unclassified_following_reset", "pending_manual_reset", "pending_mux_change")
SCAN_CATS = ("single_pore", "reserved_pore", "unavailable", "saturated", "zero", "multiple", "other")
GOOD_SCAN = ("single_pore", "reserved_pore")
HORIZONS = (6, 12, 24, 48, 72)
TIME_AXIS_DAYS = 21                      # run-time accumulators of a pass cover three weeks; later reads go to the last bin


def grid_bin_s(channels: int | None) -> int:
    """Column width of the channel x time grid per flow cell class."""
    c = channels or 0
    return 7200 if c >= 2000 else 1800 if c >= 400 else 600


class Inputs:
    """``provenance.inputs`` of one bundle, with the stable local ids in1, in2, ..."""

    def __init__(self, root_id: str | None = None, run_relpath: str | None = None, run_dir: str | None = None):
        self.items: list[dict] = []
        self.root_id, self.run_relpath, self.run_dir = root_id, run_relpath, run_dir

    def add(self, role: str, path: str, used_by: list[str], *, state: str = "ok", fingerprint: dict | bool | None = True,
            size: int | None = None, mtime: float | None = None, **extra) -> str:
        iid = f"in{len(self.items) + 1}"
        rel = None
        if self.root_id and self.run_dir and path.startswith(self.run_dir.rstrip("/") + "/"):
            rel = os.path.join(self.run_relpath or "", path[len(self.run_dir.rstrip("/")) + 1:])
        d: dict = {"id": iid, "role": role, "path": path, "root": self.root_id if rel is not None else None, "relpath": rel}
        try:
            st = os.stat(path)
            d["size"] = st.st_size if size is None else size
            d["mtime"] = from_epoch(st.st_mtime if mtime is None else mtime)
            if fingerprint is True:
                d["fingerprint"] = _fp.fingerprint(path)
            else:
                d["fingerprint"] = fingerprint or None
        except OSError:
            d.update(size=size, mtime=from_epoch(mtime) if mtime else None, fingerprint=fingerprint if isinstance(fingerprint, dict) else None)
            state = "missing"
        d["state"] = state
        d["scope"] = "run"
        d["link"] = {"method": "run_dir", "confidence": "certain"}
        d["used_by"] = list(used_by)
        d.update(extra)
        self.items.append(d)
        return iid

    def get(self, iid: str) -> dict:
        return next(i for i in self.items if i["id"] == iid)


def note(level: str, code: str, message: str, **extra) -> dict:
    d = {"level": level, "code": code, "message": message}
    d.update({k: v for k, v in extra.items() if v is not None})
    return d


def r4(x):
    return None if x is None else round(float(x), 4)
