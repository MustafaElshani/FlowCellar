"""Basecaller sample sheets, and telling what a small file is from its first bytes.

A basecaller sample sheet is a CSV whose header holds ``barcode``, ``alias`` and ``flow_cell_id`` or ``position_id``,
in any order (Dorado accepts any order). It says which sample sits on which barcode of which flow cell.
"""
from __future__ import annotations

import csv
import gzip
import io
import json
import os
import struct
import zlib

from flowcellar.util.errors import InputError
from flowcellar.util.text import clean_text

SNIFF_BYTES = 64 << 10
SHEET_NEEDS = ("barcode", "alias")
SHEET_ONE_OF = ("flow_cell_id", "position_id")
MAX_SHEET_ROWS = 4096


def is_sheet_header(columns) -> bool:
    cols = {c.strip().lower() for c in columns}
    return all(c in cols for c in SHEET_NEEDS) and any(c in cols for c in SHEET_ONE_OF)


def read_sample_sheet(path: str) -> dict:
    """``{"columns", "rows": [{flow_cell_id, position_id, experiment_id, kit, barcode, alias, sample_id, type}],
    "flow_cells": [...], "experiments": [...]}``. Values are cleaned text; an empty cell is None."""
    if os.path.getsize(path) > (4 << 20):
        raise InputError(f"{path}: too large for a sample sheet")
    with open(path, "rb") as fh:
        text = fh.read().decode("utf-8-sig", errors="replace")
    rd = csv.reader(io.StringIO(text))
    try:
        header = [c.strip().lower() for c in next(rd)]
    except StopIteration:
        raise InputError(f"{path}: empty file") from None
    if not is_sheet_header(header):
        raise InputError(f"{path}: not a basecaller sample sheet (needs the columns barcode, alias and flow_cell_id or position_id)")
    keep = ("flow_cell_id", "position_id", "experiment_id", "kit", "barcode", "alias", "sample_id", "type", "protocol_run_id")
    rows = []
    for raw in rd:
        if not any(c.strip() for c in raw):
            continue
        d = {h: (clean_text(v.strip(), 256) or None) for h, v in zip(header, raw)}
        rows.append({k: d.get(k) for k in keep})
        if len(rows) >= MAX_SHEET_ROWS:
            break
    return {"columns": header, "rows": rows,
            "flow_cells": sorted({r["flow_cell_id"] for r in rows if r["flow_cell_id"]}),
            "experiments": sorted({r["experiment_id"] for r in rows if r["experiment_id"]}),
            "run_ids": sorted({r["protocol_run_id"] for r in rows if r["protocol_run_id"]})}


def _first_bgzf_block(head: bytes) -> bytes | None:
    """The uncompressed bytes of the first BGZF block, when ``head`` starts with one."""
    if len(head) < 18 or head[:4] != b"\x1f\x8b\x08\x04":
        return None
    try:
        return zlib.decompressobj(31).decompress(head, 1 << 20)
    except zlib.error:
        return None


def sniff(path: str) -> str | None:
    """The role of a file told by at most its first 64 kB, or None.

    Recognised: a basecaller sample sheet (by its columns), a partial (by its ``format``), pychopper statistics,
    salmon, oarfish and kallisto metadata, a read summary with basecalled columns, and BAM files (``dorado.bam``
    when the header names no reference; an aligned BAM otherwise: ``align.transcriptome_bam`` when the references
    look like transcripts, else ``align.genome_bam``).
    """
    try:
        with open(path, "rb") as fh:
            head = fh.read(SNIFF_BYTES)
    except OSError:
        return None
    if not head:
        return None
    if head[:2] == b"\x1f\x8b":
        block = _first_bgzf_block(head)
        if block is not None and block[:4] == b"BAM\x01":
            (l_text,) = struct.unpack_from("<i", block, 4)
            text = block[8: 8 + max(l_text, 0)].decode("utf-8", errors="replace")
            names = [ln.split("SN:", 1)[1].split("\t", 1)[0] for ln in text.splitlines() if ln.startswith("@SQ") and "SN:" in ln]
            if not names:
                # a header longer than the first block may still hold references: only "no @SQ in a complete header" is unaligned
                return "dorado.bam" if l_text <= len(block) - 8 else "align.genome_bam"
            transcript_like = sum(1 for n in names[:200] if "|" in n or n[:4] in ("ENST", "ENSM", "NM_0", "NR_0", "XM_0", "XR_0") or n[:5] == "MSTRG")
            if transcript_like * 2 > min(len(names), 200) or l_text > (4 << 20):
                return "align.transcriptome_bam"
            return "align.genome_bam"
        try:
            text = gzip.GzipFile(fileobj=io.BytesIO(head)).read(4096).decode("utf-8", errors="replace")
        except (OSError, EOFError, zlib.error):
            try:
                text = zlib.decompressobj(31).decompress(head, 4096).decode("utf-8", errors="replace")
            except zlib.error:
                return None
        return "flowcellar.partial" if '"format":"flowcellar.partial"' in text.replace(" ", "") else None
    if b"\x00" in head[:4096]:
        return None                         # binary: none of the text files this function knows
    try:
        text = head.decode("utf-8-sig")
    except UnicodeDecodeError:
        text = head.decode("utf-8", errors="replace")
        if text[:4096].count("\ufffd") > 8:
            return None
    first = text.split("\n", 1)[0].rstrip("\r")
    tabs = [c.strip() for c in first.split("\t")]
    if tabs[:3] == ["Category", "Name", "Value"]:
        return "pychopper.stats"
    if "read_id" in tabs and "sequence_length_template" in tabs and "mean_qscore_template" in tabs:
        return "dorado.summary"
    if "," in first:
        try:
            header = next(csv.reader([first]), [])
        except csv.Error:
            header = []
        if is_sheet_header(header):
            return "basecaller.sample_sheet"
        from flowcellar.readdata.readtable import is_read_table_header
        if is_read_table_header([h.strip() for h in header]):      # one row per read: transcript, gene, biotype, tail
            return "polya.read_table"
    stripped = text.lstrip()
    if stripped.startswith("{"):
        try:
            d = json.loads(text)
        except ValueError:
            compact = text.replace(" ", "").replace("\n", "").replace("\t", "")
            if '"format":"flowcellar.partial"' in compact:
                return "flowcellar.partial"
            if '"salmon_version"' in compact and '"num_processed"' in compact:
                return "quant.salmon"
            if '"discard_table"' in compact and '"valid_best_aln"' in compact:
                return "quant.oarfish"
            if '"n_pseudoaligned"' in compact:
                return "quant.kallisto"
            return None
        if isinstance(d, dict):
            if d.get("format") == "flowcellar.partial":
                return "flowcellar.partial"
            if "salmon_version" in d and "num_processed" in d:
                return "quant.salmon"
            if isinstance(d.get("discard_table"), dict) and "valid_best_aln" in d["discard_table"]:
                return "quant.oarfish"
            if "n_pseudoaligned" in d and "n_processed" in d:
                return "quant.kallisto"
    return None
