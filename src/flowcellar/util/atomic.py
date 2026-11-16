"""Writing files so that a reader never sees half of one, and reading them back.

Every write is a whole-file replace: temporary name in the same directory, fsync, os.replace.
"""
from __future__ import annotations

import errno
import gzip
import io
import json
import os
import tempfile
import time

MAX_DECOMPRESSED = 256 << 20          # readers cap the decompressed size at 256 MB


class _Encoder(json.JSONEncoder):
    def default(self, o):
        try:
            import numpy as np
        except ImportError:      # pragma: no cover
            np = None
        if np is not None:
            if isinstance(o, np.integer):
                return int(o)
            if isinstance(o, np.floating):
                return float(o)
            if isinstance(o, np.ndarray):
                return o.tolist()
            if isinstance(o, np.bool_):
                return bool(o)
        return super().default(o)


def dumps(obj, *, indent: int | None = None) -> str:
    """JSON text in the member order given, shortest round-trip numbers, no NaN."""
    return json.dumps(obj, cls=_Encoder, ensure_ascii=False, allow_nan=False, indent=indent,
                      separators=(",", ": ") if indent else (",", ":"))


def gzip_bytes(data: bytes) -> bytes:
    """gzip level 9 with header time 0 and no file name: identical input gives identical bytes."""
    buf = io.BytesIO()
    with gzip.GzipFile(filename="", mode="wb", fileobj=buf, compresslevel=9, mtime=0) as gz:
        gz.write(data)
    return buf.getvalue()


def atomic_write_bytes(path: str | os.PathLike, data: bytes) -> None:
    path = os.fspath(path)
    d = os.path.dirname(path) or "."
    os.makedirs(d, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix="." + os.path.basename(path) + ".", suffix=".tmp", dir=d)
    try:
        with os.fdopen(fd, "wb") as fh:
            fh.write(data)
            fh.flush()
            os.fsync(fh.fileno())
        # mkstemp makes 0600; give the file the mode a plain open() would have given it
        umask = os.umask(0)
        os.umask(umask)
        os.chmod(tmp, 0o666 & ~umask)
        os.replace(tmp, path)
    except BaseException:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise


def atomic_write_text(path: str | os.PathLike, text: str) -> None:
    atomic_write_bytes(path, text.encode("utf-8"))


def atomic_write_json(path: str | os.PathLike, obj, *, indent: int | None = None) -> None:
    """Plain JSON, or gzip JSON when the name ends in .gz."""
    data = (dumps(obj, indent=indent) + ("\n" if indent else "")).encode("utf-8")
    atomic_write_bytes(path, gzip_bytes(data) if os.fspath(path).endswith(".gz") else data)


def read_bytes(path: str | os.PathLike, *, retries: int = 3, pause: float = 0.5) -> bytes:
    """Read a library file whole. ESTALE and ENOENT are retried: on NFS a file held open may have been replaced."""
    last: OSError | None = None
    for attempt in range(retries):
        try:
            with open(path, "rb") as fh:
                return fh.read()
        except OSError as e:
            if e.errno not in (errno.ESTALE, errno.ENOENT) or attempt == retries - 1:
                raise
            last = e
            time.sleep(pause)
    raise last  # pragma: no cover


def gunzip_capped(data: bytes, cap: int = MAX_DECOMPRESSED) -> bytes:
    with gzip.GzipFile(fileobj=io.BytesIO(data)) as gz:
        out = gz.read(cap + 1)
    if len(out) > cap:
        from flowcellar.util.errors import ResourceLimitError
        raise ResourceLimitError(f"decompressed size is above {cap >> 20} MB")
    return out


def read_json(path: str | os.PathLike, *, retries: int = 3):
    """Read plain or gzip JSON (told apart by the gzip magic bytes, not by the name)."""
    data = read_bytes(path, retries=retries, pause=0.2)
    if data[:2] == b"\x1f\x8b":
        data = gunzip_capped(data)
    return json.loads(data)
