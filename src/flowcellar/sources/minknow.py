"""Readers of the small files of a MinKNOW run folder (tier 1, quick detail).

Final summary, sample sheet, barcode table, report JSON, throughput and pore activity. Every reader is header or
key driven, so the file layouts of every MinKNOW generation are read by the same code: nothing here depends on a
column position or on a version number. A reader never raises on damaged input: it returns what it could read and
says in ``state`` and ``problems`` what went wrong, so that one cut file costs one section and never the bundle.
"""
from __future__ import annotations

import csv
import fnmatch
import io
import json
import os
import re
from dataclasses import dataclass, field

import numpy as np

from flowcellar.model import runfolder
from flowcellar.model.runkey import parse_minknow_name, run_key_from_folder_name
from flowcellar.util.errors import InputError

REPORT_CAP = 256 << 20            # a report above this size is not parsed
HEAD_BYTES = 256 << 10            # how much of a report the identity probe reads
SMALL_CAP = 64 << 20


def g(d, *keys, default=None):
    """Nested ``dict.get`` that survives a missing level or a level of another type."""
    for k in keys:
        if not isinstance(d, dict) or k not in d:
            return default
        d = d[k]
    return d


def to_int(v, default=0):
    try:
        return int(v)
    except (TypeError, ValueError):
        try:
            return int(float(v))
        except (TypeError, ValueError, OverflowError):
            return default


def to_float(v, default=None):
    try:
        f = float(v)
    except (TypeError, ValueError):
        return default
    return f if np.isfinite(f) else default


# ------------------------------------------------------------------------------------------- the folder
@dataclass
class Listing:
    """What one look at a run folder shows: no file is opened."""
    run_dir: str
    name: str
    kind: str | None                                    # run, run_raw_only, lookalike, None
    files: dict[str, str] = field(default_factory=dict)         # tier 1 role -> absolute path (the chosen file)
    sizes: dict[str, int] = field(default_factory=dict)
    mtimes: dict[str, float] = field(default_factory=dict)
    extra: dict[str, list[str]] = field(default_factory=dict)   # role -> other files that matched the pattern
    names: list[str] = field(default_factory=list)
    subdirs: list[str] = field(default_factory=list)
    newest_mtime: float = 0.0
    oldest_mtime: float = 0.0
    raw_name_keys: list[tuple[str, str]] = field(default_factory=list)   # (flow cell, run short id) read from raw file names


def list_run_folder(run_dir: str) -> Listing:
    """List a folder once (and ``other_reports`` inside it) and give every tier 1 role its file.

    When several files match one role the largest non-empty one is taken and the others are remembered in
    ``extra``. Raises InputError when the folder cannot be listed.
    """
    run_dir = os.path.abspath(run_dir)
    name = os.path.basename(run_dir.rstrip("/"))
    try:
        entries = list(os.scandir(run_dir))
    except OSError as e:
        raise InputError(f"{run_dir}: cannot be read ({e.strerror or e})") from None
    out = Listing(run_dir=run_dir, name=name, kind=None)
    stats: dict[str, os.stat_result] = {}
    for e in entries:
        if e.name.startswith("."):
            continue
        try:
            st = e.stat(follow_symlinks=False)
            is_dir = e.is_dir(follow_symlinks=False)
        except OSError:
            continue
        out.names.append(e.name)
        if is_dir:
            out.subdirs.append(e.name)
        else:
            stats[e.name] = st
        out.newest_mtime = max(out.newest_mtime, st.st_mtime)
        out.oldest_mtime = st.st_mtime if not out.oldest_mtime else min(out.oldest_mtime, st.st_mtime)
    out.names.sort()
    out.subdirs.sort()
    out.kind = runfolder.classify_listing(name, out.names, out.subdirs)
    cands: dict[str, list[tuple[int, str, float]]] = {}
    for role, pattern in runfolder.TIER1_ROLES.items():
        if "/" in pattern:
            continue
        for n, st in stats.items():
            if fnmatch.fnmatchcase(n, pattern):
                cands.setdefault(role, []).append((st.st_size, os.path.join(run_dir, n), st.st_mtime))
    if "other_reports" in out.subdirs:
        sub = os.path.join(run_dir, "other_reports")
        try:
            for e in os.scandir(sub):
                if e.name.startswith(".") or not e.is_file(follow_symlinks=False):
                    continue
                st = e.stat(follow_symlinks=False)
                out.newest_mtime = max(out.newest_mtime, st.st_mtime)
                for role, pattern in runfolder.TIER1_ROLES.items():
                    if "/" in pattern and fnmatch.fnmatchcase("other_reports/" + e.name, pattern):
                        cands.setdefault(role, []).append((st.st_size, e.path, st.st_mtime))
                p = parse_minknow_name(e.name)
                if p:
                    out.raw_name_keys.append((p.flow_cell_id, p.run_short))
        except OSError:
            pass
    for role, cs in cands.items():
        cs.sort(key=lambda c: (-c[0], c[1]))
        nonempty = [c for c in cs if c[0] > 0] or cs
        size, path, mtime = nonempty[0]
        out.files[role], out.sizes[role], out.mtimes[role] = path, size, mtime
        if len(cs) > 1:
            out.extra[role] = [c[1] for c in cs if c[1] != path]
    if out.kind == "run_raw_only":
        for d in out.subdirs:
            if not (d.startswith("pod5") or d.startswith("fast5")):
                continue
            try:
                with os.scandir(os.path.join(run_dir, d)) as it:
                    for i, e in enumerate(it):
                        p = parse_minknow_name(e.name)
                        if p:
                            out.raw_name_keys.append((p.flow_cell_id, p.run_short))
                        if i >= 40:
                            break
            except OSError:
                pass
    return out


def key_from_names(listing: Listing) -> tuple[str, str] | None:
    """(flow cell id, run short id) from the names of the identity and raw files, then from the folder name."""
    for path in listing.files.values():
        p = parse_minknow_name(os.path.basename(path))
        if p:
            return p.flow_cell_id, p.run_short
    if listing.raw_name_keys:
        return listing.raw_name_keys[0]
    k = run_key_from_folder_name(listing.name)
    if k:
        fc, short = k.split("_", 1)
        return fc, short
    return None


def acquisition_from_names(listing: Listing) -> str | None:
    for role in ("minknow.final_summary", "minknow.sequencing_summary", "minknow.throughput", "minknow.pore_activity"):
        path = listing.files.get(role)
        p = parse_minknow_name(os.path.basename(path)) if path else None
        if p and p.acquisition_short:
            return p.acquisition_short
    return None


# ------------------------------------------------------------------------------------- small text files
def _read_text(path: str, cap: int = SMALL_CAP) -> str:
    with open(path, "rb") as fh:
        data = fh.read(cap + 1)
    if len(data) > cap:
        raise ValueError(f"larger than {cap >> 20} MB")
    return data.decode("utf-8", errors="replace")


def read_final_summary(path: str) -> dict[str, str]:
    """``key=value`` lines. An unreadable file gives an empty dict."""
    out: dict[str, str] = {}
    try:
        for line in _read_text(path, 4 << 20).splitlines():
            if "=" in line:
                k, v = line.split("=", 1)
                out[k.strip()] = v.strip()
    except (OSError, ValueError):
        pass
    return out


def read_csv_rows(path: str, delimiter: str = ",", cap: int = 16 << 20) -> list[dict[str, str]]:
    """A small delimited table as dicts by header name; empty on any problem."""
    try:
        text = _read_text(path, cap)
    except (OSError, ValueError):
        return []
    try:
        return [dict(r) for r in csv.DictReader(io.StringIO(text, newline=""), delimiter=delimiter)]
    except csv.Error:
        return []


def read_sample_sheet(path: str) -> list[dict[str, str]]:
    return read_csv_rows(path, ",")


def read_barcode_alignment(path: str) -> list[dict[str, str]]:
    return read_csv_rows(path, "\t")


# ------------------------------------------------------------------------------------------ report JSON
def read_report_head(path: str) -> dict:
    """``protocol_run_info``, ``software_versions`` and ``host`` of a report without loading the rest.

    MinKNOW writes these before the acquisitions, and one MinKNOW generation writes 60 to 70 MB of voltage
    samples after them. Reads 256 kB; when that does not hold the whole ``protocol_run_info`` object the file is
    read whole if it is under 16 MB. Returns ``{}`` when nothing could be read.
    """
    try:
        size = os.path.getsize(path)
        with open(path, "rb") as fh:
            head = fh.read(HEAD_BYTES).decode("utf-8", errors="replace")
    except OSError:
        return {}
    out: dict = {}
    dec = json.JSONDecoder()
    for key in ("protocol_run_info", "software_versions", "host"):
        m = re.search(r'"%s"\s*:\s*' % key, head)
        if not m:
            continue
        try:
            out[key], _ = dec.raw_decode(head, m.end())
        except ValueError:
            pass
    if "protocol_run_info" not in out and size <= (16 << 20):
        rep, _ = read_report(path)
        out = {k: rep[k] for k in ("protocol_run_info", "software_versions", "host") if isinstance(rep, dict) and k in rep}
    return out


def read_report(path: str) -> tuple[dict, str | None]:
    """(the report, None) or ({}, why it could not be read)."""
    try:
        size = os.path.getsize(path)
        if size > REPORT_CAP:
            return {}, f"the report file is larger than {REPORT_CAP >> 20} MB and was not read"
        with open(path, "rb") as fh:
            data = fh.read()
    except OSError as e:
        return {}, f"the report file could not be read ({e.strerror or e})"
    try:
        rep = json.loads(data)
    except ValueError:
        return {}, "the report file is not complete JSON (cut off or damaged)"
    if not isinstance(rep, dict):
        return {}, "the report file does not hold a report"
    return rep, None


def software_versions(rep: dict) -> dict:
    """Inside ``protocol_run_info`` from MinKNOW 5.7, at the top level before."""
    sw = g(rep, "protocol_run_info", "software_versions") or rep.get("software_versions") or {}
    return sw if isinstance(sw, dict) else {}


def minknow_version(sw: dict) -> str | None:
    mk = sw.get("minknow")
    if isinstance(mk, dict):
        if mk.get("full"):
            return str(mk["full"])
        if mk.get("major") is not None:
            return ".".join(str(mk.get(k, 0)) for k in ("major", "minor", "patch"))
        return None
    return str(mk) if mk else None


def sequencing_acquisition(rep: dict) -> tuple[dict | None, int]:
    """(the sequencing acquisition, how many acquisitions say they are one). Chosen by purpose, never by
    position; among several the one with the most reads; without any purpose the one with the most reads."""
    acqs = [a for a in (rep.get("acquisitions") or []) if isinstance(a, dict)]
    if not acqs:
        return None, 0

    def reads(a):
        return to_int(g(a, "acquisition_run_info", "yield_summary", "read_count"), 0)
    seq = [a for a in acqs if g(a, "acquisition_run_info", "config_summary", "purpose") == "SEQUENCING"]
    return (max(seq, key=reads) if seq else max(acqs, key=reads)), len(seq)


# --------------------------------------------------------------------------------------------- throughput
@dataclass
class Throughput:
    """Cumulative output per minute, as MinKNOW wrote it."""
    state: str = "missing"                  # ok, truncated, error, missing
    problems: list[str] = field(default_factory=list)
    rows: int = 0
    minute: np.ndarray = field(default_factory=lambda: np.zeros(0, np.int64))
    reads: np.ndarray = field(default_factory=lambda: np.zeros(0, np.int64))
    events: np.ndarray = field(default_factory=lambda: np.zeros(0, np.int64))
    est_bases: np.ndarray = field(default_factory=lambda: np.zeros(0, np.int64))
    raw_samples: np.ndarray = field(default_factory=lambda: np.zeros(0, np.int64))
    basecalled_pass: np.ndarray = field(default_factory=lambda: np.zeros(0, np.int64))
    basecalled_fail: np.ndarray = field(default_factory=lambda: np.zeros(0, np.int64))
    basecalled_bases: np.ndarray = field(default_factory=lambda: np.zeros(0, np.int64))

    @property
    def usable(self) -> bool:
        return self.rows > 0


_TP_COLS = {"minute": "Experiment Time (minutes)", "reads": "Reads", "events": "Selected Events",
            "est_bases": "Estimated Bases", "raw_samples": "Selected Raw Samples",
            "basecalled_pass": "Basecalled Reads Passed", "basecalled_fail": "Basecalled Reads Failed",
            "basecalled_bases": "Basecalled Bases"}


def _complete_lines(text: str) -> tuple[list[str], bool]:
    """(lines, cut): the last line is dropped when the text does not end in a line break."""
    cut = bool(text) and not text.endswith("\n")
    lines = text.splitlines()
    if cut and lines:
        lines = lines[:-1]
    return lines, cut


def read_throughput(path: str | None) -> Throughput:
    t = Throughput()
    if not path:
        return t
    try:
        lines, cut = _complete_lines(_read_text(path))
    except (OSError, ValueError) as e:
        t.state = "error"
        t.problems.append(f"could not be read ({e})")
        return t
    if not lines:
        t.state = "error"
        t.problems.append("the file is empty")
        return t
    header = [h.strip() for h in lines[0].split(",")]
    idx = {k: header.index(v) for k, v in _TP_COLS.items() if v in header}
    need = ("minute", "reads", "events")
    if any(k not in idx for k in need):
        t.state = "error"
        t.problems.append("columns missing: " + ", ".join(_TP_COLS[k] for k in need if k not in idx))
        return t
    cols: dict[str, list[int]] = {k: [] for k in idx}
    bad = 0
    ncol = len(header)
    for line in lines[1:]:
        parts = line.split(",")
        if len(parts) != ncol:
            bad += 1
            continue
        try:
            vals = {k: int(parts[i]) for k, i in idx.items()}
        except ValueError:
            try:
                vals = {k: int(float(parts[i])) for k, i in idx.items()}
            except ValueError:
                bad += 1
                continue
        for k, v in vals.items():
            cols[k].append(v)
    for k in _TP_COLS:
        setattr(t, k, np.asarray(cols.get(k, []), np.int64) if k in cols else np.zeros(0, np.int64))
    t.rows = len(t.minute)
    if t.rows and (np.diff(t.minute) < 0).any():              # rows out of order: sort by minute
        order = np.argsort(t.minute, kind="stable")
        for k in idx:
            setattr(t, k, getattr(t, k)[order])
    if bad:
        t.problems.append(f"{bad} rows could not be read and were left out")
    if cut:
        t.problems.append("the last line is cut off and was left out")
    t.state = "ok" if t.rows and not cut and not bad else "truncated" if t.rows else "error"
    return t


# ------------------------------------------------------------------------------------------ pore activity
@dataclass
class PoreActivity:
    """Samples per channel state and minute."""
    state: str = "missing"
    problems: list[str] = field(default_factory=list)
    rows: int = 0
    last_minute: int = -1
    states: dict[str, np.ndarray] = field(default_factory=dict)     # state -> samples per minute (length last_minute + 1)

    @property
    def usable(self) -> bool:
        return self.rows > 0 and self.last_minute >= 0


def read_pore_activity(path: str | None) -> PoreActivity:
    pa = PoreActivity()
    if not path:
        return pa
    try:
        lines, cut = _complete_lines(_read_text(path))
    except (OSError, ValueError) as e:
        pa.state = "error"
        pa.problems.append(f"could not be read ({e})")
        return pa
    if not lines:
        pa.state = "error"
        pa.problems.append("the file is empty")
        return pa
    header = [h.strip() for h in lines[0].split(",")]
    try:
        i_s, i_m, i_v = (header.index(c) for c in ("Channel State", "Experiment Time (minutes)", "State Time (samples)"))
    except ValueError:
        pa.state = "error"
        pa.problems.append("columns missing: Channel State, Experiment Time (minutes) or State Time (samples)")
        return pa
    names: list[str] = []
    mins: list[int] = []
    vals: list[int] = []
    bad = 0
    for line in lines[1:]:
        parts = line.split(",")
        try:
            m, v = int(parts[i_m]), int(parts[i_v])
            s = parts[i_s].strip()
        except (ValueError, IndexError):
            bad += 1
            continue
        if m < 0 or v < 0 or not s:
            bad += 1
            continue
        names.append(s)
        mins.append(m)
        vals.append(v)
    pa.rows = len(mins)
    if pa.rows:
        m_arr, v_arr = np.asarray(mins, np.int64), np.asarray(vals, np.int64)
        s_arr = np.asarray(names)
        pa.last_minute = int(m_arr.max())
        for s in sorted(set(names)):
            sel = s_arr == s
            pa.states[s] = np.bincount(m_arr[sel], weights=v_arr[sel], minlength=pa.last_minute + 1).astype(np.int64)
    if bad:
        pa.problems.append(f"{bad} rows could not be read and were left out")
    if cut:
        pa.problems.append("the last line is cut off and was left out")
    pa.state = "ok" if pa.rows and not cut and not bad else "truncated" if pa.rows else "error"
    return pa
