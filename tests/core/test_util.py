"""Atomic writing, fingerprints, times, exit codes."""
import os

import pytest

from flowcellar.util import atomic, fingerprint, timeutil
from flowcellar.util.errors import ExitCode, NotBuiltError, ResourceLimitError, not_built


def test_atomic_write_replaces_whole_files_and_leaves_no_temporary_file(tmp_path):
    p = tmp_path / "deep" / "card.json"
    atomic.atomic_write_json(p, {"a": 1})
    atomic.atomic_write_json(p, {"a": 2, "b": [1.5, None]})
    assert atomic.read_json(p) == {"a": 2, "b": [1.5, None]}
    assert os.listdir(p.parent) == ["card.json"]
    gz = tmp_path / "index.json.gz"
    atomic.atomic_write_json(gz, {"n": 3})
    assert gz.read_bytes()[:2] == b"\x1f\x8b" and atomic.read_json(gz) == {"n": 3}
    with pytest.raises(ValueError):                    # NaN never reaches a file
        atomic.dumps({"x": float("nan")})


def test_a_failed_write_leaves_the_old_file(tmp_path):
    p = tmp_path / "status.json"
    atomic.atomic_write_json(p, {"ok": True})
    with pytest.raises(TypeError):
        atomic.atomic_write_json(p, {"bad": object()})
    assert atomic.read_json(p) == {"ok": True} and os.listdir(tmp_path) == ["status.json"]


def test_decompressed_size_is_capped():
    data = atomic.gzip_bytes(b"0" * 4096)
    assert atomic.gunzip_capped(data) == b"0" * 4096
    with pytest.raises(ResourceLimitError):
        atomic.gunzip_capped(data, cap=1024)


def test_fingerprints(tmp_path):
    small = tmp_path / "small.txt"
    small.write_bytes(b"abc")
    fp = fingerprint.fingerprint(small)
    assert fp == {"method": "sha256", "value": "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"}
    big = tmp_path / "big.bin"
    with open(big, "wb") as fh:                        # sparse: 65 MiB without writing 65 MiB
        fh.seek((65 << 20) - 1)
        fh.write(b"\x01")
    assert fingerprint.fingerprint(big)["method"] == "quick-v1"
    assert fingerprint.stat_key(small)[0] == 3
    assert fingerprint.hash_obj({"b": 1, "a": [1, 2]}) == fingerprint.hash_obj({"a": [1, 2], "b": 1})
    assert len(fingerprint.hash_obj({"a": 1}, 16)) == 16


def test_timestamps_are_utc_with_z(monkeypatch):
    assert timeutil.iso_z("2026-03-02T10:00:00.123456+01:00") == "2026-03-02T09:00:00.123456Z"
    assert timeutil.iso_z("2026-03-02T09:00:00.5Z", digits=0) == "2026-03-02T09:00:00Z"
    assert timeutil.iso_z("yesterday") is None and timeutil.parse_ts(None) is None
    assert timeutil.parse_ts("2026-03-02T09:00:00.25Z") - timeutil.parse_ts("2026-03-02T09:00:00Z") == pytest.approx(0.25)
    monkeypatch.setenv("SOURCE_DATE_EPOCH", "0")
    assert timeutil.now_z() == "1970-01-01T00:00:00Z"


def test_exit_codes():
    assert (ExitCode.GATE, ExitCode.INPUT, ExitCode.PARTIAL, ExitCode.LIBRARY, ExitCode.VALIDATION, ExitCode.INTERRUPTED) == (3, 4, 5, 6, 7, 130)
    e = not_built("flowcellar scan")
    assert isinstance(e, NotBuiltError) and "not built yet" in e.message and "flowcellar scan" in e.message
