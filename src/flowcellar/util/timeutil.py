"""Timestamps: everything stored is RFC 3339 UTC with a Z suffix."""
from __future__ import annotations

import datetime as _dt
import os
import re

_ISO = re.compile(r"^(\d{4}-\d{2}-\d{2})[T ](\d{2}:\d{2}:\d{2})(\.\d+)?(Z|[+-]\d{2}:?\d{2})?$")


def iso_z(ts: str | None, digits: int | None = None) -> str | None:
    """Any ISO timestamp of the MinKNOW files -> RFC 3339 UTC with a Z suffix; None when it cannot be read."""
    if not ts:
        return None
    m = _ISO.match(ts.strip())
    if not m:
        return None
    frac, tz = m.group(3) or "", m.group(4) or "Z"
    dt = _dt.datetime.strptime(m.group(1) + "T" + m.group(2), "%Y-%m-%dT%H:%M:%S")
    if tz not in ("Z", "+00:00", "+0000"):
        sign = 1 if tz[0] == "+" else -1
        dt = dt - sign * _dt.timedelta(hours=int(tz[1:3]), minutes=int(tz[-2:]))
    frac = frac[:10]
    if digits is not None:
        frac = frac[: digits + 1] if digits else ""
    return dt.strftime("%Y-%m-%dT%H:%M:%S") + frac + "Z"


def parse_ts(ts: str | None) -> float | None:
    """Seconds since the epoch of an ISO timestamp (fraction kept); None when it cannot be read."""
    z = iso_z(ts)
    if not z:
        return None
    m = re.match(r"^(.{19})(\.\d+)?Z$", z)
    dt = _dt.datetime.strptime(m.group(1), "%Y-%m-%dT%H:%M:%S").replace(tzinfo=_dt.UTC)
    return dt.timestamp() + (float(m.group(2)) if m.group(2) else 0.0)


def from_epoch(seconds: float) -> str:
    return _dt.datetime.fromtimestamp(seconds, _dt.UTC).strftime("%Y-%m-%dT%H:%M:%SZ")


def now_z() -> str:
    """Now, or SOURCE_DATE_EPOCH when it is set (deterministic output)."""
    sde = os.environ.get("SOURCE_DATE_EPOCH")
    if sde and sde.isdigit():
        return from_epoch(int(sde))
    return _dt.datetime.now(_dt.UTC).strftime("%Y-%m-%dT%H:%M:%SZ")


def mtime_z(path: str | os.PathLike) -> str:
    return from_epoch(os.path.getmtime(path))
