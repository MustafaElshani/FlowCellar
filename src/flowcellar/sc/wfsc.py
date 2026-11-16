"""Readers of a wf-single-cell output folder, for both file layouts (workflow 2.x and 3.x).

    params.json                  3.x: kit, parameters, the path of the sibling Space Ranger BAM (no workflow version)
    <alias>/config_stats.json    input reads, full-length reads, adapter layouts, strands (2.x; 3.x per documentation)
    <alias>/*feature_bc_matrix*  MatrixMarket triplets, streamed: molecules and genes per unit
    <alias>/*mito*per*cell*.tsv  mitochondrial share per unit
    <alias>/whitelist.tsv        barcodes kept by cell calling
    execution/trace.txt          tasks, failures, time and memory of the workflow
    wf-single-cell-report.html   last resort, for a named short list: the workflow version, and in 3.x the input
                                 reads, the full-length reads, the adapter layouts and the reported saturation

The per-read table is read by ``flowcellar.sc.readpass``.
"""
from __future__ import annotations

import base64
import csv
import gzip
import json
import os
import re
from html.parser import HTMLParser

import numpy as np

from flowcellar.sc.spaceranger import short_path
from flowcellar.util.errors import InputError

ADAPTER_CLASSES = ("full_length", "single_adapter1", "single_adapter2", "double_adapter1", "double_adapter2",
                   "no_adapters", "other")
#: the parameters that change what the numbers mean, in the order they are shown
PARAMETER_KEYS = ("kit", "full_length_only", "min_read_qual", "barcode_max_ed", "barcode_min_ed_diff", "barcode_min_quality",
                  "gene_assigns_minqv", "matrix_min_genes", "matrix_min_cells", "matrix_max_mito", "mito_prefix",
                  "expected_cells", "estimate_cell_count", "ref_genome_dir")
_PATH_KEYS = frozenset(("ref_genome_dir", "bam", "fastq", "spaceranger_bam", "adapter_stats", "out_dir"))


def read_json(path: str) -> dict:
    try:
        with open(path, encoding="utf-8") as fh:
            d = json.load(fh)
    except ValueError:
        raise InputError(f"{path}: not readable as JSON") from None
    if not isinstance(d, dict):
        raise InputError(f"{path}: a JSON object was expected")
    return d


def assay_of_kit(kit: str | None) -> str | None:
    """``visium_hd:v1`` -> visium_hd, ``visium:v1`` -> visium, ``3prime:v3`` and the other droplet kits -> singlecell."""
    if not kit:
        return None
    name = str(kit).split(":", 1)[0].strip().lower()
    if name == "visium_hd":
        return "visium_hd"
    if name == "visium":
        return "visium"
    return "singlecell"


def parameters(params: dict) -> list[tuple[str, str]]:
    """The parameter table of a bundle: the listed keys that are present, paths cut to two components."""
    out = []
    for k in PARAMETER_KEYS:
        if k in params and params[k] is not None:
            v = params[k]
            text = short_path(v) if k in _PATH_KEYS else (("true" if v else "false") if isinstance(v, bool) else str(v))
            out.append((k, text))
    return out


# ---------------------------------------------------------------------------------------------------------
# count matrices (MatrixMarket), streamed
# ---------------------------------------------------------------------------------------------------------

def mtx_size_line(matrix_dir: str) -> tuple[int, int, int]:
    """(features, barcodes, entries): the size line of ``matrix.mtx.gz``."""
    with gzip.open(os.path.join(matrix_dir, "matrix.mtx.gz"), "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("%"):
                continue
            parts = line.split()
            if len(parts) < 3:
                break
            return int(parts[0]), int(parts[1]), int(parts[2])
    raise InputError(f"{matrix_dir}: matrix.mtx.gz has no size line")


def _lines(path: str) -> list[str]:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        return [ln.rstrip("\n") for ln in fh]


def mtx_per_barcode(matrix_dir: str, *, mito_prefix: str | None = "MT-") -> dict:
    """Molecules, genes and mitochondrial molecules per barcode from the triplets of one matrix folder. The
    matrix is never built: the triplets are read in blocks and summed with ``bincount``."""
    import pyarrow as pa
    import pyarrow.csv as pacsv

    n_feat, n_bc, n_entries = mtx_size_line(matrix_dir)
    barcodes = _lines(os.path.join(matrix_dir, "barcodes.tsv.gz"))
    feats = [ln.split("\t") for ln in _lines(os.path.join(matrix_dir, "features.tsv.gz"))]
    names = np.array([f[1] if len(f) > 1 else f[0] for f in feats]) if feats else np.zeros(0, dtype="U1")
    is_mt = np.char.startswith(np.char.upper(names), mito_prefix.upper()) if (mito_prefix and len(names)) else np.zeros(len(names), bool)
    umis = np.zeros(n_bc, dtype=np.float64)
    genes = np.zeros(n_bc, dtype=np.int64)
    mt = np.zeros(n_bc, dtype=np.float64)
    total = 0.0
    entries = 0
    if n_entries:
        # the lines before the first triplet: comments and the size line
        with gzip.open(os.path.join(matrix_dir, "matrix.mtx.gz"), "rt", encoding="utf-8", errors="replace") as fh:
            skip = 0
            for line in fh:
                skip += 1
                if not line.startswith("%"):
                    break
        reader = pacsv.open_csv(
            os.path.join(matrix_dir, "matrix.mtx.gz"),
            read_options=pacsv.ReadOptions(skip_rows=skip, column_names=["f", "b", "v"], block_size=32 << 20),
            parse_options=pacsv.ParseOptions(delimiter=" "),
            convert_options=pacsv.ConvertOptions(column_types={"f": pa.int32(), "b": pa.int32(), "v": pa.float64()}))
        for batch in reader:
            f = batch.column(0).to_numpy(zero_copy_only=False) - 1
            b = batch.column(1).to_numpy(zero_copy_only=False) - 1
            v = batch.column(2).to_numpy(zero_copy_only=False)
            entries += len(f)
            umis += np.bincount(b, weights=v, minlength=n_bc)
            genes += np.bincount(b, minlength=n_bc)
            if is_mt.any():
                mt += np.bincount(b, weights=np.where(is_mt[f], v, 0.0), minlength=n_bc)
            total += float(v.sum())
    return {"barcodes": barcodes, "umis": np.rint(umis).astype(np.int64), "genes": genes, "mito_umis": np.rint(mt).astype(np.int64),
            "n_features": n_feat, "n_barcodes": n_bc, "entries": entries, "entries_stated": n_entries,
            "total_umis": int(round(total)), "n_mito_features": int(is_mt.sum())}


# ---------------------------------------------------------------------------------------------------------
# small tables
# ---------------------------------------------------------------------------------------------------------

def mito_fractions(path: str) -> np.ndarray:
    """Mitochondrial share per unit, 0..1, from the workflow's two-column table (``CB``, ``mito_pct`` in percent)."""
    import pyarrow as pa
    import pyarrow.csv as pacsv

    try:
        t = pacsv.read_csv(path, parse_options=pacsv.ParseOptions(delimiter="\t"),
                           convert_options=pacsv.ConvertOptions(include_columns=["mito_pct"], column_types={"mito_pct": pa.float64()}))
    except (pa.ArrowInvalid, KeyError):
        raise InputError(f"{path}: no column mito_pct") from None
    v = t.column(0).to_numpy(zero_copy_only=False)
    return v[~np.isnan(v)] / 100.0


def count_lines(path: str) -> int:
    n = 0
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            n += chunk.count(b"\n")
    return n


_DUR = re.compile(r"([\d.]+)\s*(ms|s|m|h|d)")
_MEM = {"B": 1, "KB": 1 << 10, "MB": 1 << 20, "GB": 1 << 30, "TB": 1 << 40}


def _seconds(text: str) -> float:
    unit = {"ms": 0.001, "s": 1, "m": 60, "h": 3600, "d": 86400}
    return sum(float(n) * unit[u] for n, u in _DUR.findall(text or ""))


def _bytes(text: str) -> float | None:
    parts = (text or "").split()
    if len(parts) != 2 or parts[1].upper() not in _MEM:
        return None
    try:
        return float(parts[0]) * _MEM[parts[1].upper()]
    except ValueError:
        return None


def _epoch(text: str) -> float | None:
    from datetime import datetime
    try:
        return datetime.strptime(text.strip()[:19], "%Y-%m-%d %H:%M:%S").timestamp()
    except ValueError:
        return None


def trace(path: str) -> dict:
    """The Nextflow trace of the workflow: tasks, failed tasks, wall time, task hours, peak memory, and the five
    processes that took the longest."""
    tasks = failed = cached = 0
    seconds = 0.0
    peak = 0.0
    first = last = None
    by_proc: dict[str, list[float]] = {}
    with open(path, encoding="utf-8", errors="replace", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            tasks += 1
            status = (r.get("status") or "").strip()
            if status == "CACHED":
                cached += 1
            if status not in ("COMPLETED", "CACHED") or (r.get("exit") or "0").strip() not in ("0", "-"):
                failed += 1
            rt = _seconds(r.get("realtime") or "")
            seconds += rt
            mem = _bytes(r.get("peak_rss") or "")
            if mem:
                peak = max(peak, mem)
            t0 = _epoch(r.get("submit") or "")
            if t0 is not None:
                first = t0 if first is None else min(first, t0)
                end = t0 + _seconds(r.get("duration") or "")
                last = end if last is None else max(last, end)
            name = (r.get("name") or "").split(" (")[0].split(":")[-1]
            e = by_proc.setdefault(name, [0, 0.0, 0.0])
            e[0] += 1
            e[1] += rt
            e[2] = max(e[2], mem or 0.0)
    top = sorted(by_proc.items(), key=lambda kv: -kv[1][1])[:5]
    return {"tasks": tasks, "failed": failed, "cached": cached, "task_hours": seconds / 3600.0,
            "wall_hours": ((last - first) / 3600.0) if (first is not None and last is not None) else None,
            "peak_rss_bytes": int(peak) if peak else None,
            "top": [{"process": k, "tasks": int(v[0]), "hours": v[1] / 3600.0, "peak_rss_bytes": int(v[2]) if v[2] else None} for k, v in top]}


# ---------------------------------------------------------------------------------------------------------
# the workflow's HTML report: tables, the version in the footer, the embedded figure data
# ---------------------------------------------------------------------------------------------------------

class _Tables(HTMLParser):
    """Every table of the page with the heading above it. Script and style blocks are never looked into."""

    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.heading = ""
        self.tables: list[tuple[str, list[list[str]]]] = []
        self._cap: str | None = None
        self._buf: list[str] = []
        self._rows: list[list[str]] | None = None
        self._row: list[str] | None = None
        self._skip = 0
        self.footer: list[str] = []
        self._in_footer = False

    def handle_starttag(self, tag, attrs):
        if tag in ("script", "style"):
            self._skip += 1
        elif tag == "footer":
            self._in_footer = True
        elif tag in ("h1", "h2", "h3", "h4", "h5"):
            self._cap, self._buf = tag, []
        elif tag == "table":
            self._rows = []
        elif tag == "tr" and self._rows is not None:
            self._row = []
        elif tag in ("td", "th") and self._row is not None:
            self._cap, self._buf = tag, []

    def handle_data(self, data):
        if self._skip:
            return
        if self._cap:
            self._buf.append(data)
        if self._in_footer:
            self.footer.append(data)

    def handle_endtag(self, tag):
        if tag in ("script", "style"):
            self._skip = max(0, self._skip - 1)
        elif tag == "footer":
            self._in_footer = False
        if self._cap == tag:
            text = " ".join("".join(self._buf).split())
            if tag.startswith("h"):
                self.heading = text
            elif self._row is not None:
                self._row.append(text)
            self._cap = None
        if tag == "tr" and self._rows is not None and self._row is not None:
            self._rows.append(self._row)
            self._row = None
        if tag == "table" and self._rows is not None:
            self.tables.append((self.heading, self._rows))
            self._rows = None


_VERSION = re.compile(r"wf-single-cell\s+(?:nextflow\s+)?workflow\s*\(\s*v?([0-9][0-9A-Za-z.+_-]*)\s*\)")
_DOCS_JSON = re.compile(r"const docs_json = '(.*?)';\s*\n", re.S)


def _int(text: str) -> int | None:
    t = (text or "").replace(",", "").strip()
    return int(t) if re.fullmatch(r"-?\d+", t) else None


def _bokeh_array(v):
    if isinstance(v, dict) and v.get("type") == "ndarray":
        arr = v.get("array")
        if isinstance(arr, dict) and arr.get("type") == "bytes":
            order = "<" if v.get("order", "little") == "little" else ">"
            return np.frombuffer(base64.b64decode(arr["data"]), dtype=np.dtype(v["dtype"]).newbyteorder(order)).tolist()
        if isinstance(arr, list):
            return arr
    return v if isinstance(v, list) else None


def _find(o, name: str, out: list) -> None:
    if isinstance(o, dict):
        if o.get("type") == "object" and o.get("name") == name:
            out.append(o)
        for v in o.values():
            _find(v, name, out)
    elif isinstance(o, list):
        for v in o:
            _find(v, name, out)


def _figures(text: str) -> list[dict]:
    """The figures of the embedded Bokeh document: ``[{"title", "axes", "columns": {name: list}}]``."""
    m = _DOCS_JSON.search(text)
    if not m:
        return []
    try:
        doc = json.loads(m.group(1).replace("\\\\", "\\").replace("\\'", "'"))
    except ValueError:
        return []
    figs: list = []
    _find(doc, "Figure", figs)
    out = []
    for f in figs:
        a = f.get("attributes") or {}
        titles: list = []
        _find(a.get("title", {}), "Title", titles)
        axes: list = []
        for nm in ("LinearAxis", "LogAxis", "CategoricalAxis"):
            _find(a, nm, axes)
        sources: list = []
        _find(a, "ColumnDataSource", sources)
        cols: dict = {}
        for s in sources[:2]:
            data = (s.get("attributes") or {}).get("data") or {}
            for k, v in (data.get("entries") or []) if isinstance(data, dict) else []:
                arr = _bokeh_array(v)
                if arr is not None and k not in cols:
                    cols[k] = arr
        out.append({"title": (titles[0].get("attributes") or {}).get("text") if titles else None,
                    "axes": [(x.get("attributes") or {}).get("axis_label") for x in axes], "columns": cols})
    return out


def report(path: str, alias: str | None = None) -> dict:
    """What is taken from the workflow's HTML report. Every member may be None: an unknown layout gives less,
    never an error.

    ``{"version", "input_reads", "full_length", "valid_barcode", "gene_assigned", "transcript_assigned",
    "barcodes_with_data", "mapped", "unmapped", "supplementary", "adapter_pct", "saturation_reported",
    "parameters"}``; for a 2.x report ``input_reads`` and ``cells`` come from its sample table.
    """
    with open(path, encoding="utf-8", errors="replace") as fh:
        text = fh.read()
    p = _Tables()
    try:
        p.feed(text)
    except Exception:  # noqa: BLE001  a page the parser cannot follow gives what was read so far
        pass
    out: dict = {k: None for k in ("version", "input_reads", "full_length", "valid_barcode", "gene_assigned",
                                   "transcript_assigned", "barcodes_with_data", "mapped", "unmapped", "supplementary",
                                   "adapter_pct", "saturation_reported", "parameters", "cells")}
    m = _VERSION.search(" ".join(" ".join(p.footer).split()))
    if m:
        out["version"] = m.group(1)
    for heading, rows in p.tables:
        h = heading.lower()
        flat = {r[0].lower(): r[1] for r in rows if len(r) == 2}
        if h.startswith("experiment summary"):
            out["input_reads"] = _int(flat.get("input reads", ""))
            hit = next((v for k, v in flat.items() if k.endswith("with data") or k == "cells"), None)
            out["barcodes_with_data"] = _int(hit or "")
        elif h.startswith("alignment / feature summary"):
            out["mapped"], out["unmapped"] = _int(flat.get("mapped", "")), _int(flat.get("unmapped", ""))
            out["supplementary"] = _int(flat.get("supplementary", ""))
        elif h.startswith("read assignment summary") and rows and len(rows) >= 2:
            head = [c.lower() for c in rows[0]]
            counts = next((r for r in rows[1:] if r and r[0].lower() == "count"), None)
            if counts:
                for key, label in (("full_length", "full length"), ("valid_barcode", "valid barcode"),
                                   ("gene_assigned", "gene assigned"), ("transcript_assigned", "transcript assigned")):
                    if label in head and head.index(label) < len(counts):
                        out[key] = _int(counts[head.index(label)])
        elif h.startswith("single cell sample summary") and len(rows) >= 2:
            head = [c.lower() for c in rows[0]]
            row = next((r for r in rows[1:] if alias is None or (r and r[0] == alias)), None)
            if row:
                if "reads" in head:
                    out["input_reads"] = _int(row[head.index("reads")])
                if "cells" in head:
                    out["cells"] = _int(row[head.index("cells")])
        elif h.startswith("alignment summary") and len(rows) >= 2:
            head = [c.lower() for c in rows[0]]
            row = next((r for r in rows[1:] if alias is None or (r and r[0] == alias)), None)
            if row:
                for key, label in (("mapped", "primary"), ("unmapped", "unmapped"), ("supplementary", "supplementary")):
                    if label in head:
                        out[key] = _int(row[head.index(label)])
        elif h.startswith("workflow parameters"):
            out["parameters"] = {r[0]: r[1] for r in rows if len(r) == 2 and r[0] != "Key"}
    for fig in _figures(text):
        cols = fig["columns"]
        if "groups" in cols and len(cols) >= 2:
            name = alias if alias in cols else next((k for k in cols if k != "groups"), None)
            if name and isinstance(cols[name], list) and len(cols[name]) == len(cols["groups"]):
                out["adapter_pct"] = {str(g): float(v) for g, v in zip(cols["groups"], cols[name])}
        elif fig["title"] == "Sequencing saturation" and cols.get("x") and cols.get("y"):
            i = int(np.argmax(cols["x"]))
            out["saturation_reported"] = {"reads": int(cols["x"][i]), "value": float(cols["y"][i])}
    return out
