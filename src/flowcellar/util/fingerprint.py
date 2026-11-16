"""Fingerprints of input files and hashes of JSON values."""
from __future__ import annotations

import hashlib
import json
import os

QUICK_ABOVE = 64 << 20      # files above 64 MiB get the quick-v1 fingerprint
_MIB = 1 << 20


def fingerprint(path: str | os.PathLike) -> dict:
    """``{"method": "sha256" | "quick-v1", "value": hex}``.

    sha256 of the whole file up to 64 MiB; above that quick-v1 = SHA-256 of the size (decimal text), the
    first MiB and the last MiB. The file is opened read-only and never memory-mapped.
    """
    size = os.path.getsize(path)
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        if size <= QUICK_ABOVE:
            for blk in iter(lambda: fh.read(_MIB), b""):
                h.update(blk)
            return {"method": "sha256", "value": h.hexdigest()}
        h.update(str(size).encode())
        h.update(fh.read(_MIB))
        fh.seek(size - _MIB)
        h.update(fh.read(_MIB))
        return {"method": "quick-v1", "value": h.hexdigest()}


def stat_key(path: str | os.PathLike) -> tuple[int, int]:
    """(size, mtime in nanoseconds): the cheap test of "did this file change" before a fingerprint is taken."""
    st = os.stat(path)
    return st.st_size, st.st_mtime_ns


def canonical_json(obj) -> bytes:
    """Sorted keys, no spaces, UTF-8: the form every hash is taken over."""
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False).encode("utf-8")


def sha256_hex(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def hash_obj(obj, length: int = 64) -> str:
    """SHA-256 of the canonical JSON of ``obj``, cut to ``length`` hex characters."""
    return sha256_hex(canonical_json(obj))[:length]
