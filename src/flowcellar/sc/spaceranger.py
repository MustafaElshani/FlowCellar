"""Readers of a Space Ranger ``outs/`` folder. Data files first; the report's HTML for a short, named list only.

    metrics_summary.csv      every headline number, looked up by column NAME (the set changes with product and version)
    web_summary.html         the JSON object ``const data = {...}``: alarms, diagnostics, the sample table
    molecule_info.h5         reads per molecule (saturation), reads per barcode, the read funnel, chemistry, slide
    feature_slice.h5         Visium HD: every layer of the 2 um grid (molecules, reads, tissue mask, image, genes)
    *_feature_bc_matrix.h5   molecules and genes per barcode
    spatial/tissue_positions.csv   Visium: tissue mask and grid position per spot
    possorted_genome_bam.bam(.bai) mapped reads per contig, from the header and the index only

HDF5 files need ``h5py`` (the ``sc`` extra); every function that opens one says so in one line when it is
missing. Nothing is loaded whole unless it is small: an HDF5 file is read dataset by dataset.
"""
from __future__ import annotations

import csv
import gzip
import json
import os
import re
import struct

import numpy as np

from flowcellar.util.errors import InputError, MissingExtraError

H5_HINT = "maps, saturation and per-unit distributions of Space Ranger output"


def h5py_or_none():
    try:
        import h5py  # noqa: PLC0415  the sc extra: loaded only when an HDF5 file is opened
    except ImportError:
        return None
    return h5py


def _h5():
    h5py = h5py_or_none()
    if h5py is None:
        raise MissingExtraError("h5py", "sc", "Reading " + H5_HINT)
    return h5py


def _number(text: str):
    t = text.strip()
    if t == "":
        return None
    try:
        v = float(t.replace(",", "")) if not t.endswith("%") else float(t[:-1].replace(",", "")) / 100.0
    except ValueError:
        return t
    return int(v) if v.is_integer() and "." not in t and "e" not in t.lower() and not t.endswith("%") else v


# ---------------------------------------------------------------------------------------------------------
# metrics_summary.csv
# ---------------------------------------------------------------------------------------------------------

def metrics_summary(path: str) -> dict:
    """Column name -> value of the one data row. Fractions are 0..1 as written; counts are integers."""
    with open(path, encoding="utf-8-sig", newline="") as fh:
        rows = list(csv.reader(fh))
    if len(rows) < 2 or not rows[0]:
        raise InputError(f"{path}: not a Space Ranger metrics summary (one header row and one row of values expected)")
    return {name.strip(): _number(value) for name, value in zip(rows[0], rows[1])}


def col(ms: dict, *names: str):
    """The first of ``names`` that the summary has (a column changed its name between versions)."""
    for n in names:
        if ms.get(n) is not None:
            return ms[n]
    return None


# ---------------------------------------------------------------------------------------------------------
# web_summary.html: the embedded JSON object, and nothing else of the page
# ---------------------------------------------------------------------------------------------------------
_DATA_START = re.compile(r"const data\s*=\s*")
_TAGS = re.compile(r"<[^>]+>")


def _tables(o, out: dict, path: str = "") -> None:
    """Every ``{"rows": [[label, value], ...]}`` table of the object, by its path."""
    if isinstance(o, dict):
        rows = o.get("rows")
        if isinstance(rows, list) and rows and all(isinstance(r, list) for r in rows):
            out[path] = rows
        for k, v in o.items():
            if k != "_resources":
                _tables(v, out, f"{path}.{k}" if path else k)
    elif isinstance(o, list):
        for i, v in enumerate(o):
            if isinstance(v, (dict, list)):
                _tables(v, out, f"{path}[{i}]")


def _find_key(o, key: str):
    """The first value stored under ``key`` anywhere in the object (the layout moved between versions)."""
    if isinstance(o, dict):
        if key in o and o[key] is not None:
            return o[key]
        for k, v in o.items():
            if k != "_resources":
                hit = _find_key(v, key)
                if hit is not None:
                    return hit
    elif isinstance(o, list):
        for v in o:
            hit = _find_key(v, key)
            if hit is not None:
                return hit
    return None


def web_summary(path: str) -> dict:
    """``{"alarms", "diagnostics", "sample", "command_line"}`` from the JSON object inside the page.

    ``alarms`` are the pipeline's own messages, word for word: ``[{"level", "title", "value", "message"}]``.
    ``sample`` is the sample table (Sample ID, Chemistry, Slide Serial Number, Transcriptome, Pipeline Version,
    alignment, CytAssist run). Raises InputError when the page holds no such object.
    """
    with open(path, encoding="utf-8", errors="replace") as fh:
        text = fh.read()
    m = _DATA_START.search(text)
    if not m:
        raise InputError(f"{path}: the page holds no data object in the layout this version knows")
    try:
        data = json.JSONDecoder().raw_decode(text, m.end())[0]
    except ValueError:
        raise InputError(f"{path}: the data object of the page cannot be read") from None
    del text
    alarms_obj = _find_key(data, "alarms")
    raw_alarms = alarms_obj.get("alarms") if isinstance(alarms_obj, dict) else alarms_obj
    alarms = []
    for a in raw_alarms or []:
        if not isinstance(a, dict):
            continue
        alarms.append({"level": str(a.get("level") or "").upper(), "title": str(a.get("title") or ""),
                       "value": a.get("formatted_value") if a.get("formatted_value") is not None else a.get("raw_value"),
                       "message": _TAGS.sub("", str(a.get("message") or "")).strip()})
    diagnostics = _find_key(data, "diagnostics")
    tables: dict = {}
    _tables(data, tables)
    sample: dict = {}
    for rows in tables.values():
        if any(isinstance(r, list) and len(r) == 2 and r[0] == "Sample ID" for r in rows):
            sample = {str(r[0]): r[1] for r in rows if isinstance(r, list) and len(r) == 2}
            break
    cmd = None
    cl = _find_key(data, "command_line")
    if cl is not None:
        m2 = re.search(r"(spaceranger\s+count[^<]*)", json.dumps(cl))
        cmd = m2.group(1).replace("\\n", " ").strip() if m2 else None
    return {"alarms": alarms, "diagnostics": diagnostics if isinstance(diagnostics, dict) else {}, "sample": sample,
            "command_line": cmd, "tables": {k: v for k, v in tables.items() if len(v) <= 40}}


# ---------------------------------------------------------------------------------------------------------
# molecule_info.h5
# ---------------------------------------------------------------------------------------------------------

def molecule_info(path: str, *, per_barcode: bool = True) -> dict:
    """``{"metrics", "molecules", "reads", "count_hist", "reads_per_barcode", "molecules_per_barcode", "barcodes"}``.

    ``count_hist[c]`` = molecules seen ``c`` times. ``per_barcode=False`` reads ``count`` only (Visium HD: the
    per-barcode arrays would be 9 million names long and are not needed, the grid layers give them).
    """
    h5py = _h5()
    with h5py.File(path, "r") as f:
        raw = f["metrics_json"][()]
        metrics = json.loads(raw.decode("utf-8") if isinstance(raw, bytes) else raw)
        count = f["count"][:]
        out = {"metrics": metrics, "molecules": int(len(count)), "reads": int(count.sum(dtype=np.int64)),
               "count_hist": np.bincount(count).astype(np.int64)}
        if per_barcode:
            idx = f["barcode_idx"][:].astype(np.int64)
            n = int(f["barcodes"].shape[0])
            out["reads_per_barcode"] = np.bincount(idx, weights=count, minlength=n).astype(np.int64)
            out["molecules_per_barcode"] = np.bincount(idx, minlength=n).astype(np.int64)
            out["barcodes"] = [b.decode("ascii", "replace") for b in f["barcodes"][:]]
    return out


def molecule_info_metrics(path: str) -> dict:
    """Only the small JSON of the file: version, chemistry, slide, the read funnel."""
    h5py = _h5()
    with h5py.File(path, "r") as f:
        raw = f["metrics_json"][()]
    return json.loads(raw.decode("utf-8") if isinstance(raw, bytes) else raw)


# ---------------------------------------------------------------------------------------------------------
# feature-barcode matrices (HDF5, CSC)
# ---------------------------------------------------------------------------------------------------------

def matrix_h5(path: str, *, mito_prefix: str | None = "MT-", want_data: bool = True) -> dict:
    """Per-barcode molecules, genes and mitochondrial molecules of one matrix, without building the matrix.

    Genes per barcode = ``diff(indptr)``; molecules = sums of ``data`` between the pointers; the mitochondrial
    part needs ``indices``. ``want_data=False`` reads ``indptr`` and the names only (genes per barcode).
    """
    h5py = _h5()
    with h5py.File(path, "r") as f:
        m = f["matrix"]
        indptr = m["indptr"][:].astype(np.int64)
        barcodes = [b.decode("ascii", "replace") for b in m["barcodes"][:]]
        names = np.array([b.decode("utf-8", "replace") for b in m["features/name"][:]])
        genes = np.diff(indptr)
        out = {"barcodes": barcodes, "genes": genes, "n_features": int(len(names)),
               "software_version": _attr(f, "software_version"), "chemistry": _attr(f, "chemistry_description")}
        if want_data:
            data = m["data"][:].astype(np.int64)
            csum = np.concatenate([[0], np.cumsum(data)])
            out["umis"] = csum[indptr[1:]] - csum[indptr[:-1]]
            out["total_umis"] = int(data.sum())
            if mito_prefix:
                is_mt = np.char.startswith(np.char.upper(names), mito_prefix.upper())
                indices = m["indices"][:]
                out["n_mito_features"] = int(is_mt.sum())
                mt = np.where(is_mt[indices], data, 0)
                msum = np.concatenate([[0], np.cumsum(mt)])
                out["mito_umis"] = msum[indptr[1:]] - msum[indptr[:-1]]
                seen = np.zeros(len(names), dtype=bool)
                seen[indices] = True
                out["features_seen"] = seen
                out["indices"], out["data"] = indices, data
    return out


def _attr(f, name: str) -> str | None:
    v = f.attrs.get(name)
    if v is None:
        return None
    if isinstance(v, bytes):
        return v.decode("utf-8", "replace")
    if isinstance(v, np.ndarray):
        return str(v.tolist())
    return str(v)


# ---------------------------------------------------------------------------------------------------------
# spatial/tissue_positions.csv (Visium)
# ---------------------------------------------------------------------------------------------------------

def tissue_positions(path: str) -> dict[str, tuple[int, int, int]]:
    """barcode -> (in_tissue, array_row, array_col). The older file of the same content has no header row."""
    out: dict[str, tuple[int, int, int]] = {}
    with open(path, encoding="utf-8", newline="") as fh:
        for r in csv.reader(fh):
            if len(r) < 4 or r[0] == "barcode":
                continue
            out[r[0]] = (int(r[1]), int(r[2]), int(r[3]))
    return out


# ---------------------------------------------------------------------------------------------------------
# feature_slice.h5 (Visium HD)
# ---------------------------------------------------------------------------------------------------------

class FeatureSlice:
    """The layers of ``feature_slice.h5`` on the 2 um grid. Each layer is a triple (row, col, data)."""

    def __init__(self, path: str):
        self.h5py = _h5()
        self.f = self.h5py.File(path, "r")
        meta = self.f.attrs.get("metadata_json")
        self.meta = json.loads(meta.decode("utf-8") if isinstance(meta, bytes) else meta) if meta is not None else {}
        self.nrows = int(self.meta.get("nrows") or 0)
        self.ncols = int(self.meta.get("ncols") or 0)
        self.pitch_um = float(self.meta.get("spot_pitch") or 2.0)
        self.software_version = _attr(self.f, "software_version")

    def close(self) -> None:
        self.f.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def has(self, name: str) -> bool:
        return name in self.f

    def layer(self, name: str):
        g = self.f[name]
        return g["row"][:].astype(np.int64), g["col"][:].astype(np.int64), g["data"][:]

    def rebin(self, name: str, factor: int, *, binary: bool = False) -> np.ndarray:
        """One layer summed onto a grid ``factor`` times coarser (float64). ``binary``: count non-zero entries."""
        r, c, v = self.layer(name)
        return rebin_sum(r, c, (v > 0) if binary else v, self.nrows, self.ncols, factor)

    def mito_layers(self, prefix: str = "MT-"):
        """(row, col, data) of every slice of a feature whose name starts with ``prefix``, and how many there are."""
        names = [b.decode("utf-8", "replace") for b in self.f["features/name"][:]]
        idx = [i for i, n in enumerate(names) if n.upper().startswith(prefix.upper())]
        rs, cs, vs = [], [], []
        slices = self.f["feature_slices"]
        for i in idx:
            k = str(i)
            if k in slices:
                g = slices[k]
                rs.append(g["row"][:].astype(np.int64))
                cs.append(g["col"][:].astype(np.int64))
                vs.append(g["data"][:].astype(np.int64))
        if not rs:
            z = np.zeros(0, np.int64)
            return z, z, z, len(idx)
        return np.concatenate(rs), np.concatenate(cs), np.concatenate(vs), len(idx)


def rebin_sum(row, col, val, nrows: int, ncols: int, factor: int) -> np.ndarray:
    """Triples on a fine grid -> sums on a grid ``factor`` times coarser."""
    rr, cc = -(-nrows // factor), -(-ncols // factor)
    flat = (np.asarray(row) // factor) * cc + (np.asarray(col) // factor)
    return np.bincount(flat, weights=np.asarray(val, dtype=np.float64), minlength=rr * cc).reshape(rr, cc)


_HD_BARCODE = re.compile(r"^s_(\d{3})um_(\d{5})_(\d{5})")
_POW10 = np.array([10000, 1000, 100, 10, 1], dtype=np.int64)


def hd_rows_cols(barcodes) -> tuple[np.ndarray, np.ndarray, int]:
    """Row and column of Visium HD barcodes (``s_008um_00012_00345-1``), and their bin size in um. The names
    have a fixed width, so the digits are read as bytes: no string is split."""
    if not len(barcodes):
        return np.zeros(0, np.int64), np.zeros(0, np.int64), 0
    m = _HD_BARCODE.match(barcodes[0])
    if not m:
        raise InputError(f"not a Visium HD barcode: {barcodes[0]!r}")
    um = int(m.group(1))
    raw = np.array(barcodes, dtype="S19")
    digits = raw.view(np.uint8).reshape(len(raw), 19).astype(np.int64) - 48
    rows, cols = digits[:, 8:13] @ _POW10, digits[:, 14:19] @ _POW10
    if (digits[:, 8:13] < 0).any() or (digits[:, 8:13] > 9).any() or (digits[:, 14:19] < 0).any() or (digits[:, 14:19] > 9).any():
        raise InputError("Visium HD barcodes of an unexpected form")
    return rows, cols, um


# ---------------------------------------------------------------------------------------------------------
# BAM: the header and the index only
# ---------------------------------------------------------------------------------------------------------

def bam_header(path: str, *, max_bytes: int = 64 << 20) -> dict:
    """``{"text", "references": [(name, length)]}`` of a BAM, read from its first compressed blocks."""
    with gzip.open(path, "rb") as fh:
        head = fh.read(8)
        if len(head) < 8 or head[:4] != b"BAM\x01":
            raise InputError(f"{path}: not a BAM file")
        l_text = struct.unpack("<i", head[4:8])[0]
        if l_text < 0 or l_text > max_bytes:
            raise InputError(f"{path}: the header is larger than {max_bytes >> 20} MB")
        text = fh.read(l_text).decode("utf-8", "replace").rstrip("\x00")
        raw = fh.read(4)
        if len(raw) < 4:
            raise InputError(f"{path}: the header is cut off")
        refs = []
        for _ in range(struct.unpack("<i", raw)[0]):
            l_name = struct.unpack("<i", fh.read(4))[0]
            name = fh.read(l_name)[:-1].decode("ascii", "replace")
            refs.append((name, struct.unpack("<i", fh.read(4))[0]))
    return {"text": text, "references": refs}


def bai_counts(path: str) -> tuple[list[tuple[int, int]], int]:
    """(mapped, unmapped) records per reference from the pseudo-bins of a ``.bai`` index, and the records without
    a position. A reference without the pseudo-bin counts (0, 0)."""
    with open(path, "rb") as fh:
        data = fh.read()
    if data[:4] != b"BAI\x01":
        raise InputError(f"{path}: not a BAM index")
    view = memoryview(data)
    n_ref = struct.unpack_from("<i", view, 4)[0]
    pos = 8
    out = []
    for _ in range(n_ref):
        n_bin = struct.unpack_from("<i", view, pos)[0]
        pos += 4
        counts = (0, 0)
        for _b in range(n_bin):
            bin_id, n_chunk = struct.unpack_from("<Ii", view, pos)
            pos += 8
            if bin_id == 37450 and n_chunk == 2:
                counts = struct.unpack_from("<QQ", view, pos + 16)
            pos += 16 * n_chunk
        n_intv = struct.unpack_from("<i", view, pos)[0]
        pos += 4 + 8 * n_intv
        out.append((int(counts[0]), int(counts[1])))
    no_coor = struct.unpack_from("<Q", view, pos)[0] if pos + 8 <= len(data) else 0
    return out, int(no_coor)


_MITO_CONTIGS = ("chrM", "MT", "chrMT", "M", "mitochondrion_genome")


def bam_mapping(bam: str, index: str | None = None) -> dict:
    """Mapped and unmapped records and the share on the mitochondrial genome, from the header and the index.

    ``{"mapped", "unmapped", "chrm_mapped", "chrm_fraction", "header_text"}``; the counts are of records (a read
    with a supplementary alignment has two)."""
    index = index or bam + ".bai"
    head = bam_header(bam)
    counts, no_coor = bai_counts(index)
    names = [n for n, _ in head["references"]]
    mapped = sum(m for m, _ in counts)
    unmapped = sum(u for _, u in counts) + no_coor
    chrm = sum(m for (m, _), n in zip(counts, names) if n in _MITO_CONTIGS)
    has_m = any(n in _MITO_CONTIGS for n in names)
    return {"mapped": int(mapped), "unmapped": int(unmapped), "chrm_mapped": int(chrm) if has_m else None,
            "chrm_fraction": (chrm / mapped) if (has_m and mapped) else None, "header_text": head["text"],
            "references": head["references"]}


_SLIDE_CO = re.compile(r"slide_serial_capture_area:([A-Za-z0-9-]+)")


def slide_from_header(header_text: str) -> str | None:
    m = _SLIDE_CO.search(header_text or "")
    return m.group(1) if m else None


# ---------------------------------------------------------------------------------------------------------
# the pipestance, when it was kept beside outs/
# ---------------------------------------------------------------------------------------------------------

def pipestance_metrics(path: str, ms: dict) -> dict | None:
    """The JSON of the pipeline's reporter stage, but only when it describes the same output: such a folder can
    be older than ``outs/`` (a re-run that kept the old stage). Its read count and its unit count must agree
    with ``metrics_summary.csv``."""
    try:
        with open(path, encoding="utf-8") as fh:
            d = json.load(fh)
    except (OSError, ValueError):
        return None
    if not isinstance(d, dict):
        return None
    reads = d.get("total_read_pairs") or d.get("total_reads")
    if reads is not None and ms.get("Number of Reads") is not None and int(reads) != int(ms["Number of Reads"]):
        return None
    units_csv = col(ms, "Number of Spots Under Tissue", "Number of Bins Under Tissue 8 µm")
    units = d.get("filtered_bcs")
    if units is None and isinstance(d.get("square_008um"), dict):
        units = d["square_008um"].get("filtered_bcs")
    if units is not None and units_csv is not None and int(units) != int(units_csv):
        return None
    return d


def short_path(p: str) -> str:
    """A file-system path cut to its last two components: what a parameter table shows."""
    parts = [x for x in str(p).replace("\\", "/").split("/") if x not in ("", ".", "..")]
    return "/".join(parts[-2:]) if len(parts) > 2 else "/".join(parts)


def command_parameters(cmdline: str) -> list[tuple[str, str]]:
    """``spaceranger count --id X --slide Y ...`` -> [(option, value)], with paths cut to two components."""
    toks = cmdline.replace("=", " ").split()
    out = []
    i = 0
    while i < len(toks):
        t = toks[i]
        if not t.startswith("--"):
            i += 1
            continue
        if i + 1 < len(toks) and not toks[i + 1].startswith("--"):
            val, i = toks[i + 1], i + 2
        else:
            val, i = "true", i + 1
        out.append((t[2:], short_path(val) if "/" in val else val))
    return out


def read_small(path: str, limit: int = 1 << 20) -> str:
    with open(path, encoding="utf-8", errors="replace") as fh:
        return fh.read(limit)


def file_facts(path: str) -> dict:
    st = os.stat(path)
    return {"size": int(st.st_size), "mtime_ns": int(st.st_mtime_ns)}
