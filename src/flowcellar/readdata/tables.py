"""Table sources: the small result files of a pipeline, read whole, and read summaries, streamed.

* pychopper statistics (``Category / Name / Value``): the cDNA classification of one sample;
* quantification metadata of salmon (``aux_info/meta_info.json``), oarfish (``*.meta_info.json``) and kallisto
  (``run_info.json``): reads processed and reads counted towards a transcript, without any pass over a BAM;
* a ``dorado summary`` table, or a basecalled sequencing summary outside a run folder: the same per-read content as
  an unaligned BAM (length, mean quality, run, sample, pass or fail), streamed once.

What the pychopper numbers mean was checked on real files: ``Primers_found + rescued reads + Unusable = PassReads``,
and the reads pychopper writes (and every later step reads) are ``Primers_found + Rescue - LenFail``, where
``Rescue`` counts the pieces cut out of fused reads and ``LenFail`` the pieces that were too short. ``QcFail`` reads
were dropped before classification, so the reads given to pychopper are ``PassReads + QcFail``.
"""
from __future__ import annotations

import json
import os
import time

import numpy as np

from flowcellar.model import runkey
from flowcellar.readdata._pool import beat
from flowcellar.readdata.acc import ReadAcc
from flowcellar.util.errors import InputError, ResourceLimitError
from flowcellar.util.text import clean_text

MAX_COLLECTED = 256 << 20           # a collected file is read whole, up to 256 MB
PRIMER_CONFIGS_KEPT = 8


def _read_small(path: str) -> bytes:
    size = os.path.getsize(path)
    if size > MAX_COLLECTED:
        raise ResourceLimitError(f"{path}: {size >> 20} MB is above the limit of {MAX_COLLECTED >> 20} MB for a file read whole")
    with open(path, "rb") as fh:
        return fh.read()


def _int(v) -> int:
    return int(round(float(v or 0)))


# --------------------------------------------------------------------------------------- pychopper
def read_pychopper(path: str) -> dict:
    """One pychopper statistics table -> integer scalars, strand and primer-arrangement counts.

    ``scalars``: ``pass_reads``, ``qc_fail``, ``len_fail``, ``primers_found``, ``rescue_segments`` (pieces),
    ``rescued_reads`` (fused reads that were cut apart), ``unusable``, ``umi_detected`` (when the table has it),
    ``input`` = reads given to pychopper, ``written`` = reads it wrote. ``consistent`` is False when the classes do
    not add up to ``PassReads`` (a layout this reader does not know).
    """
    cat: dict[str, dict[str, float]] = {}
    text = _read_small(path).decode("utf-8", errors="replace")
    lines = text.splitlines()
    if not lines or [c.strip() for c in lines[0].split("\t")][:3] != ["Category", "Name", "Value"]:
        raise InputError(f"{path}: not a pychopper statistics table (no Category / Name / Value header)")
    skipped = 0
    for line in lines[1:]:
        p = line.rstrip("\r\n").split("\t")
        if len(p) != 3:
            skipped += bool(line.strip())
            continue
        try:
            cat.setdefault(p[0], {})[p[1]] = float(p[2])
        except ValueError:
            skipped += 1
    rs, cl = cat.get("ReadStats", {}), cat.get("Classification", {})
    if "PassReads" not in rs or "Primers_found" not in cl:
        raise InputError(f"{path}: pychopper statistics without ReadStats PassReads or Classification Primers_found")
    seg = cat.get("RescueSegmentNr", {})
    sc = {"pass_reads": _int(rs.get("PassReads")), "len_fail": _int(rs.get("LenFail")), "qc_fail": _int(rs.get("QcFail")),
          "primers_found": _int(cl.get("Primers_found")), "rescue_segments": _int(cl.get("Rescue")),
          "unusable": _int(cl.get("Unusable")), "rescued_reads": _int(sum(seg.values()))}
    if "Umi_detected" in rs:
        sc["umi_detected"] = _int(rs["Umi_detected"])
    sc["input"] = sc["pass_reads"] + sc["qc_fail"]
    sc["written"] = max(sc["primers_found"] + sc["rescue_segments"] - sc["len_fail"], 0)
    consistent = sc["primers_found"] + sc["rescued_reads"] + sc["unusable"] == sc["pass_reads"]
    hits = sorted(cat.get("Hits", {}).items(), key=lambda kv: (-kv[1], kv[0]))
    strand = cat.get("Strand", {})
    return {"scalars": sc, "consistent": consistent, "rows_skipped": skipped,
            "strand": {clean_text(k, 16): _int(v) for k, v in sorted(strand.items())},
            "primer_configs": [(clean_text(k, 120), _int(v)) for k, v in hits[:PRIMER_CONFIGS_KEPT]],
            "primer_configs_other": _int(sum(v for _, v in hits[PRIMER_CONFIGS_KEPT:])),
            "segments_per_rescued_read": {clean_text(k, 8): _int(v) for k, v in sorted(seg.items())},
            "cutoff_q": cat.get("Parameter", {}).get("Cutoff(q)")}


# ---------------------------------------------------------------------------------- quantification
def _json(path: str) -> dict:
    try:
        d = json.loads(_read_small(path))
    except (ValueError, UnicodeDecodeError) as e:
        raise InputError(f"{path}: not readable as JSON ({e})") from None
    if not isinstance(d, dict):
        raise InputError(f"{path}: not a JSON object")
    return d


def read_salmon_meta(path: str) -> dict:
    d = _json(path)
    if "num_processed" not in d or "num_mapped" not in d:
        raise InputError(f"{path}: not a salmon meta_info.json (no num_processed / num_mapped)")
    return {"tool": "salmon", "version": clean_text(d.get("salmon_version") or "", 40) or None, "processed": int(d["num_processed"]),
            "counted": int(d["num_mapped"]), "mode": clean_text(d.get("mapping_type") or "", 40) or None,
            "counted_field": "num_mapped", "processed_field": "num_processed"}


def read_oarfish_meta(path: str) -> dict:
    d = _json(path)
    t = d.get("discard_table")
    if not isinstance(t, dict) or "valid_best_aln" not in t:
        raise InputError(f"{path}: not an oarfish meta_info.json (no discard_table.valid_best_aln)")
    discards = {clean_text(k, 40): int(v) for k, v in t.items() if k.startswith("discard_") and isinstance(v, (int, float))}
    return {"tool": "oarfish", "version": clean_text(d.get("oarfish_version") or "", 40) or None, "processed": None,
            "counted": int(t["valid_best_aln"]), "mode": clean_text(d.get("alignment_source") or "", 40) or None,
            "counted_field": "discard_table.valid_best_aln", "processed_field": None, "discards": discards}


def read_kallisto_info(path: str) -> dict:
    d = _json(path)
    if "n_processed" not in d or "n_pseudoaligned" not in d:
        raise InputError(f"{path}: not a kallisto run_info.json (no n_processed / n_pseudoaligned)")
    return {"tool": "kallisto", "version": clean_text(d.get("kallisto_version") or "", 40) or None, "processed": int(d["n_processed"]),
            "counted": int(d["n_pseudoaligned"]), "mode": "pseudoalignment", "counted_field": "n_pseudoaligned",
            "processed_field": "n_processed", "unique": int(d.get("n_unique") or 0)}


QUANT_READERS = {"quant.salmon": read_salmon_meta, "quant.oarfish": read_oarfish_meta, "quant.kallisto": read_kallisto_info}
QUANT_PREFERENCE = ("quant.salmon", "quant.oarfish", "quant.kallisto")      # first that exists gives "counted towards a transcript"


def read_quant(role: str, path: str) -> dict:
    return QUANT_READERS[role](path)


# ----------------------------------------------------------------------------------- read summaries
_FILE_COLUMNS = ("filename_pod5", "filename", "filename_fast5", "filename_bam")
_LABEL_COLUMNS = ("alias", "barcode_arrangement", "barcode")
SUMMARY_NEEDS = ("sequence_length_template", "mean_qscore_template")


def summary_columns(path: str) -> list[str]:
    with open(path, "rb") as fh:
        first = fh.readline(1 << 16)
    return [c.strip() for c in first.decode("utf-8", errors="replace").rstrip("\r\n").split("\t")]


def _run_of_file(v):
    p = runkey.parse_minknow_name(v) if v else None
    return p.run_key if p else None


def _run_of_id(v):
    v = (v or "").replace("-", "").lower()
    return "?_" + v[:8] if len(v) >= 8 and all(ch in "0123456789abcdef" for ch in v[:8]) else None


def _passes(v):
    return None if v is None else ("pass" if str(v).strip().lower() in ("true", "1", "t", "yes") else "fail")


def _codes(batch, name, mapper):
    """(code per row, values): the distinct values of one text column mapped once each; the last value stands for null."""
    import pyarrow.compute as pc
    n = batch.num_rows
    if name is None:
        return np.zeros(n, np.int64), [None]
    enc = pc.dictionary_encode(batch.column(name))
    values = [mapper(v) for v in enc.dictionary.to_pylist()] + [None]
    idx = enc.indices.to_numpy(zero_copy_only=False)
    if idx.dtype.kind == "f":                       # nulls come back as NaN
        idx = np.where(np.isnan(idx), len(values) - 1, idx)
    return idx.astype(np.int64), values


def summary_task(args: tuple) -> dict:
    """One streaming pass over a read summary. ``args`` = (path, stop_wall). Runs in a worker process.

    Returns groups keyed by (run key, sample label, "pass" | "fail" | None): the filter is known per read only when
    the table has ``passes_filtering``.
    """
    import pyarrow as pa
    import pyarrow.csv as pacsv

    from flowcellar.readdata.bam import norm_label
    path, stop_wall = args
    t0 = time.time()
    cols = summary_columns(path)
    missing = [c for c in SUMMARY_NEEDS if c not in cols]
    if missing:
        raise InputError(f"{path}: a read summary without the columns {', '.join(missing)} (basecalling was off)")
    fcol = next((c for c in _FILE_COLUMNS if c in cols), None)
    lcol = next((c for c in _LABEL_COLUMNS if c in cols), None)
    want = [c for c in (fcol, "run_id" if "run_id" in cols else None, lcol, "passes_filtering" if "passes_filtering" in cols else None,
                        *SUMMARY_NEEDS) if c]
    types = {c: pa.string() for c in want}
    types["sequence_length_template"] = pa.float64()
    types["mean_qscore_template"] = pa.float64()
    size = os.path.getsize(path)
    rows = bad = 0
    groups: dict[tuple, ReadAcc] = {}
    finished = False
    skipped = [0]

    def wrong_columns(row) -> str:              # a row with the wrong number of columns is skipped and counted
        skipped[0] += 1
        return "skip"

    try:
        reader = pacsv.open_csv(path, read_options=pacsv.ReadOptions(block_size=8 << 20),
                                parse_options=pacsv.ParseOptions(delimiter="\t", invalid_row_handler=wrong_columns),
                                convert_options=pacsv.ConvertOptions(include_columns=want, column_types=types, strings_can_be_null=True,
                                                                     null_values=["", "nan", "NaN"]))
        for batch in reader:
            n = batch.num_rows
            if not n:
                continue
            rows += n
            length = batch.column("sequence_length_template").to_numpy(zero_copy_only=False)
            q = batch.column("mean_qscore_template").to_numpy(zero_copy_only=False).astype(np.float64)
            ok = np.isfinite(length)
            bad += int((~ok).sum())
            length = np.where(ok, length, 0.0)

            ki, kv = _codes(batch, fcol, _run_of_file)
            if "run_id" in want:                         # the run id column backs up a file name that names no run
                ri, rv = _codes(batch, "run_id", _run_of_id)
            else:
                ri, rv = np.zeros(n, np.int64), [None]
            li, lv = _codes(batch, lcol, norm_label)
            if "passes_filtering" in want:
                pi, pv = _codes(batch, "passes_filtering", _passes)
            else:
                pi, pv = np.zeros(n, np.int64), [None]
            combo = ((ki * len(rv) + ri) * len(lv) + li) * len(pv) + pi
            zeros = np.zeros(n)
            for c in np.unique(combo[ok]):
                sel = ok & (combo == c)
                c = int(c)
                p_i = c % len(pv); c //= len(pv)
                l_i = c % len(lv); c //= len(lv)
                r_i = c % len(rv); k_i = c // len(rv)
                gid = (kv[k_i] or rv[r_i], lv[l_i], pv[p_i])
                acc = groups.get(gid)
                if acc is None:
                    acc = groups[gid] = ReadAcc()
                z = zeros[: int(sel.sum())]
                acc.add(length[sel], q[sel], z, z - 1, z - 1, z, z)
            beat(n)
            if stop_wall is not None and time.time() >= stop_wall:
                break
        else:
            finished = True
    except (pa.ArrowInvalid, OSError) as e:
        return {"groups": {g: a.to_state() for g, a in groups.items()}, "rows": rows, "rows_bad": bad, "rows_skipped": skipped[0],
                "finished": False, "error": str(e)[:300], "size": size, "seconds": time.time() - t0,
                "columns": {"file": fcol, "label": lcol, "filter": "passes_filtering" in want}}
    return {"groups": {g: a.to_state() for g, a in groups.items()}, "rows": rows, "rows_bad": bad, "rows_skipped": skipped[0],
            "finished": finished, "error": None, "size": size, "seconds": time.time() - t0,
            "columns": {"file": fcol, "label": lcol, "filter": "passes_filtering" in want}}
