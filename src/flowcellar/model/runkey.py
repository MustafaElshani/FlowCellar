"""Run keys and the names MinKNOW and Dorado write.

``run_key = {flow_cell_id}_{protocol_run_id[:8]}``. File names inside a run and in Dorado output follow
``{flow cell}[_{pass|fail|skip}][_{barcode or alias}]_{protocol run id[:8]}_{acquisition id[:8]}_{n}``:
the reader takes the first field as the flow cell and the first 8-hex field after it as the run id, and
never relies on position.
"""
from __future__ import annotations

import re
from dataclasses import dataclass

FLOW_CELL_RE = re.compile(r"^[A-Z]{2,3}[0-9]{4,6}$")
HEX8_RE = re.compile(r"^[0-9a-f]{8}$")
RUN_KEY_RE = re.compile(r"^([A-Z]{2,3}[0-9]{4,6})_([0-9a-f]{8}(?:[0-9a-f]{4})?)$")
SAFE_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$")
_FOLDER_TAIL_RE = re.compile(r"(?:^|_)([A-Z]{2,3}[0-9]{4,6})_([0-9a-f]{8})$")
_STATUS = ("pass", "fail", "skip")


def make_run_key(flow_cell_id: str, protocol_run_id: str, length: int = 8) -> str:
    """Flow cell id plus the first ``length`` hex characters of the protocol run id (8, or 12 after a collision)."""
    short = protocol_run_id.replace("-", "")[:length].lower()
    if not FLOW_CELL_RE.match(flow_cell_id or "") or not re.fullmatch(r"[0-9a-f]{%d}" % length, short):
        raise ValueError(f"cannot make a run key from {flow_cell_id!r} and {protocol_run_id!r}")
    return f"{flow_cell_id}_{short}"


def is_run_key(s: str) -> bool:
    return bool(RUN_KEY_RE.match(s or ""))


def split_run_key(run_key: str) -> tuple[str, str]:
    """``("TST00001", "0a1b2c3d")``; ValueError when the text is not a run key."""
    m = RUN_KEY_RE.match(run_key or "")
    if not m:
        raise ValueError(f"not a run key: {run_key!r}")
    return m.group(1), m.group(2)


def is_safe_id(s: str) -> bool:
    """True when the text may be used as a file name and inside an address."""
    return bool(SAFE_ID_RE.match(s or "")) and ".." not in s


def _is_date(text: str) -> bool:
    """``20260302``: eight digits that name a day of this century."""
    return bool(re.fullmatch(r"20\d{2}(0[1-9]|1[0-2])(0[1-9]|[12]\d|3[01])", text))


@dataclass(frozen=True)
class ParsedName:
    flow_cell_id: str
    run_short: str                  # first 8 hex of the protocol run id
    acquisition_short: str | None   # first 8 hex of the acquisition id
    status: str | None              # pass, fail, skip
    barcode: str | None             # barcode or alias token, when there is one
    index: int | None               # trailing file number

    @property
    def run_key(self) -> str:
        return f"{self.flow_cell_id}_{self.run_short}"


def parse_minknow_name(name: str) -> ParsedName | None:
    """Read a POD5, BAM or FASTQ name, an ``fn`` tag or an ``RG`` value. None when it holds no run key.

    ``TST00001_skip_0a1b2c3d_4e5f6a7b_146.pod5`` -> flow cell TST00001, run 0a1b2c3d, status skip.
    """
    base = name.rsplit("/", 1)[-1]
    base = re.sub(r"\.(pod5|fast5|bam|fastq|fq|fastq\.gz|fq\.gz|txt|csv|tsv|json)$", "", base)
    fields = base.split("_")
    for i, f in enumerate(fields):
        if not FLOW_CELL_RE.match(f):
            continue
        rest = fields[i + 1:]
        # MinKNOW's own reports carry the start of the run between the flow cell and the run id
        # (report_TST00001_20260302_0900_0a1b2c3d.json): eight digits that are a date, followed by a time, are not
        # a run id, although every digit is also a hex digit
        start = 2 if len(rest) >= 3 and _is_date(rest[0]) and re.fullmatch(r"[0-2]\d[0-5]\d", rest[1]) and HEX8_RE.match(rest[2]) else 0
        j = next((k for k, x in enumerate(rest) if k >= start and HEX8_RE.match(x)), None)
        if j is None:
            continue
        between = rest[start:j]
        status = next((x for x in between if x in _STATUS), None)
        other = [x for x in between if x not in _STATUS]
        after = rest[j + 1:]
        acq = after[0] if after and HEX8_RE.match(after[0]) else None
        tail = after[1:] if acq else after
        idx = int(tail[0]) if tail and tail[0].isdigit() else None
        return ParsedName(f, rest[j], acq, status, "_".join(other) or None, idx)
    return None


def run_key_from_folder_name(name: str) -> str | None:
    """The ``{flow cell}_{8 hex}`` a MinKNOW run folder name ends in. The date, time and device parts are never parsed."""
    m = _FOLDER_TAIL_RE.search(name.rstrip("/").rsplit("/", 1)[-1])
    return f"{m.group(1)}_{m.group(2)}" if m else None
