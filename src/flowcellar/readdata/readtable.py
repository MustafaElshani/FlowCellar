"""One streamed pass over an annotated read table (role ``polya.read_table``) -> one partial.

An annotated read table is the CSV an annotation step writes after the transcriptome alignment: one row per read
that was given to a transcript, with the transcript's gene and biotype and, in newer tables, the poly(A) tail
length the basecaller estimated::

    read_id,transcript_id,gene_id,gene_name,length,transcript_type,mapping_quality,FLAG,ePolyA_tail_length

Such a table is 0.4 to 5 GB. It is read once, in blocks of 16 MB that end at a line break; each block is decoded by
``pyarrow.csv`` (only the columns that are counted, in 1 MB pieces side by side), nothing but the counters below is
kept, and the pass stops at a deadline when it is given one, knowing how many bytes it read.

What is counted, and of what (checked on real tables of four years):

* **rows** = primary rows of the table. A row whose ``FLAG`` has the secondary (256) or supplementary (2048) bit is
  counted apart and used for nothing. A read that a cDNA pre-processor cut into pieces has one row per piece under
  the same ``read_id`` (about 1 row in 15,000 on the tables looked at), so "rows" is reads given to a transcript,
  counting such a read once per piece.
* **per biotype** (``transcript_type``, as written): rows; rows with a tail estimate; the tail histogram.
* **tail estimate**: a value above 0. ``0`` (the basecaller found the primer but no tail), a negative value (it
  found no primer; Dorado writes ``-1``) and an empty cell or ``NA`` (no value: older tables) are three classes of
  "no estimate" and are never binned. A table without a tail column has every row in "no value".
* **per gene** (``gene_id``, or ``gene_name`` when there is no id column): rows. Kept are the number of genes with
  at least 1, 10 and 100 rows, the rows of the top 1, 10 and 100 genes, the rows of genes named ``MT-...``, a
  rank curve on a fixed set of ranks (ten per decade), and the names of the top genes when the catalogue says how many may be stored
  (``rna_top_gene_share.names_stored``). No other gene name leaves the pass.
"""
from __future__ import annotations

import os
import time

import numpy as np

from flowcellar.aggregates import build
from flowcellar.contracts import LinkedInput
from flowcellar.model import axes
from flowcellar.model import partial as partial_model
from flowcellar.readdata._pool import beat, run_tasks
from flowcellar.util.errors import InputError
from flowcellar.util.text import clean_text
from flowcellar.util.timeutil import from_epoch

VERSION = "0.1.0"
ROLE = "polya.read_table"
#: versions of the section builders that read this pass: a bump changes the partial id of the tables only
MODULES = {"rna.biotype": VERSION, "rna.genes": VERSION}
BLOCK_BYTES = 16 << 20
PARSE_BYTES = 1 << 20
TAIL_AXIS = "polya_1nt"
MAX_BIOTYPES_KEPT = 64               # more distinct values than this are folded into FOLDED, most common kept
MAX_BIOTYPES_SEEN = 5000             # more than this and the column is not a biotype column
MAX_GENES_SEEN = 2_000_000
RANK_GRID = tuple(sorted({int(round(10 ** (k / 10))) for k in range(0, 61)}))      # 1, 2, 3, 4, 5, 6, 8, 10, 13, ... 1,000,000
NOT_GIVEN = "(not given)"
FOLDED = "(other values)"
BAD_ROWS_FRACTION = 0.01
SECONDARY, SUPPLEMENTARY = 0x100, 0x800

BIOTYPE_COLUMNS = ("transcript_type", "transcript_biotype", "biotype")
TAIL_COLUMNS = ("ePolyA_tail_length", "polya_tail_length", "poly_tail_length", "polyA_tail_length", "tail_length")
FLAG_COLUMNS = ("FLAG", "flag")
MISSING_WORDS = frozenset(("", "na", "nan", "none", "null", "."))
TAIL_CODING = {"estimated": "value above 0", "zero": "value 0: the primer was found, no tail length",
               "not_found": "negative value: the primer was not found", "no_value": "empty cell or NA, or the table has no tail column"}


def names_kept() -> int:
    """How many gene names of the top of the list may be stored: the catalogue decides (0 when it says nothing)."""
    from flowcellar import catalogue
    try:
        n = int((catalogue.metric("rna_top_gene_share") or {}).get("names_stored") or 0)
    except (TypeError, ValueError):
        n = 0
    return max(0, min(n, 25))


def config_hash() -> str:
    """The setting that shapes what a pass stores (part of the partial id)."""
    return f"gene_names_{names_kept()}"


def header_columns(path: str) -> list[str]:
    with open(path, "rb") as fh:
        first = fh.readline(1 << 16)
    return [c.strip().strip('"') for c in first.decode("utf-8-sig", errors="replace").rstrip("\r\n").split(",")]


def pick_columns(cols: list[str]) -> dict:
    """Which column plays which part; ``biotype`` is None when the table is not an annotated read table."""
    def first(names):
        return next((c for c in names if c in cols), None)
    return {"biotype": first(BIOTYPE_COLUMNS), "gene": first(("gene_id", "gene_name")), "gene_name": first(("gene_name",)),
            "transcript": first(("transcript_id",)), "flag": first(FLAG_COLUMNS), "tail": first(TAIL_COLUMNS)}


def is_read_table_header(cols: list[str]) -> bool:
    c = pick_columns(cols)
    return "read_id" in cols and c["biotype"] is not None and (c["gene"] is not None or c["transcript"] is not None)


def _tail_value(text) -> float:
    """The tail length a cell holds; NaN for "no value"."""
    s = "" if text is None else str(text).strip()
    if s.lower() in MISSING_WORDS:
        return float("nan")
    try:
        v = float(s)
    except ValueError:
        return float("nan")
    return v if v == v and abs(v) != float("inf") else float("nan")


def _flag_value(text) -> int:
    try:
        return int(float(str(text).strip()))
    except (TypeError, ValueError):
        return 0


def _encoded(col):
    """(codes per row as int64, the distinct values): a column of a record batch, dictionary-encoded once."""
    import pyarrow as pa
    import pyarrow.compute as pc
    if isinstance(col, pa.ChunkedArray):
        col = col.combine_chunks()
    arr = col if pa.types.is_dictionary(col.type) else pc.dictionary_encode(col)
    values = arr.dictionary.to_pylist()
    idx = arr.indices.to_numpy(zero_copy_only=False)
    if idx.dtype.kind == "f":                       # nulls come back as NaN: they get a code of their own
        idx = np.where(np.isnan(idx), len(values), idx)
        values = values + [None]
    return idx.astype(np.int64), values


def _blocks(fh, block_bytes: int):
    """(whole rows as bytes, True when the last row has no line break) for the rest of an open file."""
    carry = b""
    while True:
        chunk = fh.read(block_bytes)
        if not chunk:
            if carry:
                yield carry, True
            return
        data = carry + chunk if carry else chunk
        cut = data.rfind(b"\n")
        if cut < 0:
            carry = data
            continue
        yield data[: cut + 1], False
        carry = data[cut + 1:]


def table_task(args: tuple) -> dict:
    """The pass. ``args`` = (path, stop_wall, gene names kept). Runs in a worker process and returns plain data."""
    import pyarrow as pa
    import pyarrow.csv as pacsv

    path, stop_wall, n_names = args
    t0 = time.time()
    size = os.path.getsize(path)
    cols = header_columns(path)
    if not is_read_table_header(cols):
        raise InputError(f"{path}: not an annotated read table (needs the columns read_id, a transcript type and a gene or transcript id)")
    C = pick_columns(cols)
    want = [c for c in dict.fromkeys((C["biotype"], C["gene"], C["gene_name"], C["transcript"], C["flag"], C["tail"])) if c]
    # every column that is counted is a column of repeated words: decoded straight into a dictionary
    types = {c: pa.dictionary(pa.int32(), pa.string()) for c in want}
    n_cells = axes.n_cells(TAIL_AXIS)

    bt_ids: dict[str, int] = {}
    bt_hist: list[np.ndarray] = []                # per biotype: internal histogram array of the estimated tails
    bt_count: list[np.ndarray] = []               # per biotype: [rows, zero, not_found, no_value]
    bt_sum: list[float] = []                      # per biotype: sum of the estimated tails
    genes: dict[str, list] = {}                   # gene key -> [rows, name]
    transcripts: set[str] = set()
    rows = primary = secondary = supplementary = no_gene = 0
    tail_min = tail_max = None
    fractional = False
    skipped = [0]
    finished = False
    error = None
    read_bytes = 0

    def wrong_columns(row) -> str:                # a row with the wrong number of columns is skipped and counted
        skipped[0] += 1
        return "skip"

    ends_with_newline = True
    try:
        with open(path, "rb") as fh:
            fh.readline(1 << 16)                  # the header, read above
            read_bytes = fh.tell()
            parse = pacsv.ParseOptions(delimiter=",", invalid_row_handler=wrong_columns)
            convert = pacsv.ConvertOptions(include_columns=want, column_types=types, strings_can_be_null=False)
            for block, unterminated in _blocks(fh, BLOCK_BYTES):
                ends_with_newline = not unterminated
                table = pacsv.read_csv(pa.BufferReader(block), read_options=pacsv.ReadOptions(column_names=cols, block_size=PARSE_BYTES),
                                       parse_options=parse, convert_options=convert)        # 1 MB pieces, decoded side by side
                read_bytes += len(block)
                batch = table
                n = batch.num_rows
                if not n:
                    continue
                rows += n
                # ---- which rows are primary
                if C["flag"]:
                    fi, fv = _encoded(batch.column(C["flag"]))
                    flags = np.array([_flag_value(v) for v in fv], np.int64)[fi]
                    sec = (flags & SECONDARY) != 0
                    sup = ((flags & SUPPLEMENTARY) != 0) & ~sec
                    keep = ~(sec | sup)
                    secondary += int(sec.sum())
                    supplementary += int(sup.sum())
                else:
                    keep = np.ones(n, bool)
                n_keep = int(keep.sum())
                primary += n_keep
                # ---- biotype and tail
                bi, bv = _encoded(batch.column(C["biotype"]))
                local = np.empty(len(bv), np.int64)
                for k, v in enumerate(bv):
                    name = clean_text(v or "", 60).strip() or NOT_GIVEN
                    g = bt_ids.get(name)
                    if g is None:
                        if len(bt_ids) >= MAX_BIOTYPES_SEEN:
                            raise InputError(f"{path}: the column {C['biotype']} holds more than {MAX_BIOTYPES_SEEN:,} distinct values; it is not a transcript type")
                        g = bt_ids[name] = len(bt_ids)
                        bt_hist.append(np.zeros(n_cells, np.int64))
                        bt_count.append(np.zeros(4, np.int64))
                        bt_sum.append(0.0)
                    local[k] = g
                if C["tail"]:
                    ti, tv = _encoded(batch.column(C["tail"]))
                    tail = np.array([_tail_value(v) for v in tv], np.float64)[ti]
                else:
                    tail = np.full(n, np.nan)
                nb = len(bv)
                has = ~np.isnan(tail)
                est = keep & has & (tail > 0)
                zero = keep & has & (tail == 0)
                neg = keep & has & (tail < 0)
                none = keep & ~has
                rows_bt = np.bincount(bi[keep], minlength=nb)
                zero_bt = np.bincount(bi[zero], minlength=nb)
                neg_bt = np.bincount(bi[neg], minlength=nb)
                none_bt = np.bincount(bi[none], minlength=nb)
                if est.any():
                    te = tail[est]
                    fractional = fractional or bool((te != np.floor(te)).any())
                    lo, hi = float(te.min()), float(te.max())
                    tail_min = lo if tail_min is None else min(tail_min, lo)
                    tail_max = hi if tail_max is None else max(tail_max, hi)
                    grid = np.bincount(bi[est] * n_cells + axes.bin_index(te, TAIL_AXIS), minlength=nb * n_cells).reshape(nb, n_cells)
                    sums = np.bincount(bi[est], weights=te, minlength=nb)
                else:
                    grid, sums = None, None
                for k in np.nonzero(rows_bt)[0]:
                    g = int(local[k])
                    bt_count[g] += (rows_bt[k], zero_bt[k], neg_bt[k], none_bt[k])
                    if grid is not None:
                        bt_hist[g] += grid[k]
                        bt_sum[g] += float(sums[k])
                # ---- genes and transcripts
                if C["gene"]:
                    gi, gv = _encoded(batch.column(C["gene"]))
                    per = np.bincount(gi[keep], minlength=len(gv))
                    present = np.nonzero(per)[0]
                    if C["gene_name"] and C["gene_name"] != C["gene"]:
                        ni, nv = _encoded(batch.column(C["gene_name"]))
                        codes_present, first_row = np.unique(gi, return_index=True)       # the first row of each gene
                        first_of = dict(zip(codes_present.tolist(), first_row.tolist()))
                    else:
                        ni = nv = first_of = None
                    for k in present.tolist():
                        key = gv[k]
                        if key is None or not str(key).strip():
                            no_gene += int(per[k])
                            continue
                        acc = genes.get(key)
                        if acc is None:
                            if len(genes) >= MAX_GENES_SEEN:
                                raise InputError(f"{path}: more than {MAX_GENES_SEEN:,} distinct values in the column {C['gene']}; it is not a gene column")
                            name = nv[int(ni[first_of[k]])] if first_of is not None else key
                            genes[key] = [int(per[k]), name]
                        else:
                            acc[0] += int(per[k])
                else:
                    no_gene += n_keep
                if C["transcript"]:
                    xi, xv = _encoded(batch.column(C["transcript"]))
                    per = np.bincount(xi[keep], minlength=len(xv))
                    transcripts.update(xv[k] for k in np.nonzero(per)[0].tolist() if xv[k])
                beat(n)
                if stop_wall is not None and time.time() >= stop_wall and read_bytes < size:
                    break
            else:
                finished = True
    except (pa.ArrowInvalid, OSError) as e:
        error = str(e)[:300]
        ends_with_newline = True

    # ---- genes: counts, the rank curve, the top of the list
    ranked = sorted(genes.values(), key=lambda a: (-a[0], str(a[1])))
    counts = np.array([a[0] for a in ranked], np.int64)
    cum = np.cumsum(counts)
    n_genes = len(counts)

    def top(k: int) -> int:
        return int(cum[min(k, n_genes) - 1]) if n_genes else 0
    if n_genes:
        # one fixed set of ranks for every table (ten per decade), so the curves of two samples share their points
        xs = sorted({r for r in RANK_GRID if r <= n_genes} | {n_genes})
        curve = {"x": xs, "y": {"reads_cum": [int(cum[x - 1]) for x in xs], "reads_at_rank": [int(counts[x - 1]) for x in xs]},
                 "x_unit": "gene rank, most reads first", "y_unit": "reads", "x_scale": "log10"}
    else:
        curve = None
    mito = sum(a[0] for a in ranked if str(a[1] or "")[:3].upper() == "MT-")
    named = [(clean_text(a[1] or "", 40) or "(no name)", int(a[0])) for a in ranked[:n_names]]
    gene_scalars = {"reads": primary, "reads_with_gene": int(cum[-1]) if n_genes else 0, "reads_no_gene": no_gene,
                    "genes_ge1": n_genes, "genes_ge10": int((counts >= 10).sum()), "genes_ge100": int((counts >= 100).sum()),
                    "transcripts_ge1": len(transcripts) if C["transcript"] else None,
                    "reads_top1": top(1), "reads_top10": top(10), "reads_top100": top(100), "reads_mito_genes": int(mito)}

    # ---- biotypes: most common first, the rare tail of a very long list folded
    order = sorted(bt_ids, key=lambda b: (-int(bt_count[bt_ids[b]][0]), b))
    kept, folded = order[:MAX_BIOTYPES_KEPT], order[MAX_BIOTYPES_KEPT:]
    biotypes = []
    for name in kept:
        g = bt_ids[name]
        biotypes.append({"name": name, "counts": bt_count[g].tolist(), "hist": bt_hist[g], "sum": bt_sum[g]})
    if folded:
        hist, cnt, s = np.zeros(n_cells, np.int64), np.zeros(4, np.int64), 0.0
        for name in folded:
            g = bt_ids[name]
            hist += bt_hist[g]
            cnt += bt_count[g]
            s += bt_sum[g]
        biotypes.append({"name": FOLDED, "counts": cnt.tolist(), "hist": hist, "sum": s})
    return {"rows": rows, "primary": primary, "secondary": secondary, "supplementary": supplementary, "rows_skipped": skipped[0],
            "biotypes": biotypes, "biotypes_seen": len(bt_ids), "biotypes_folded": len(folded),
            "gene_scalars": gene_scalars, "curve": curve, "top_genes": named, "names_kept": n_names,
            "tail_min": tail_min, "tail_max": tail_max, "tail_fractional": fractional,
            "columns": C, "all_columns": [clean_text(c, 60) for c in cols[:40]],
            "finished": finished, "error": error, "ends_with_newline": ends_with_newline, "size": size, "read_bytes": read_bytes,
            "seconds": time.time() - t0}


def _biotype_cell(b: dict) -> dict:
    rows, zero, neg, none = (int(v) for v in b["counts"])
    hist = b["hist"]
    estimated = int(np.asarray(hist).sum())
    s = b["sum"]
    h = build.hist(TAIL_AXIS, hist, missing=rows - estimated)
    return {"at": [b["name"]],
            "scalars": {"reads": rows, "tail_estimated": estimated, "tail_zero": zero, "tail_not_found": neg, "tail_no_value": none,
                        "tail_sum": int(round(s)) if float(s).is_integer() else round(float(s), 3)},
            "hists": {"rna_biotype_polya_hist": h}}


def run_pass(inp: LinkedInput, *, deadline: float | None = None, stall_timeout_s: float = 120.0, inprocess: bool = False) -> dict:
    """The partial of one annotated read table (not written). ``state`` is ok, truncated, partially_read or error."""
    c = inp.candidate
    t0 = time.time()
    p = partial_model.new_partial(role=c.role, input_path=c.path, fp=inp.fingerprint, size=c.size,
                                  mtime=from_epoch(c.mtime_ns / 1e9) if c.mtime_ns else None, modules=dict(MODULES),
                                  config_hash=config_hash(), sample=inp.sample, variant=inp.variant)
    stop_wall = None if deadline is None else time.time() + max(deadline - time.monotonic(), 0.0)
    try:
        out = run_tasks(table_task, [(c.path, stop_wall, names_kept())], workers=1, stall_timeout_s=stall_timeout_s, inprocess=inprocess)
    except Exception as e:      # noqa: BLE001  a table this reader cannot use is an input error, not a crash of the scan
        p["state"] = "error"
        p["errors"] = [{"code": "unreadable", "message": str(e)[:300]}]
        return p
    if out["stalled"] or not out["results"]:
        p["state"] = "error"
        p["errors"] = [{"code": "stalled", "message": f"no progress for {stall_timeout_s:.0f} s; the pass was abandoned"}]
        return p
    res = out["results"][0]
    cells = [_biotype_cell(b) for b in res["biotypes"]]
    p["fragments"].append(partial_model.fragment(None, "rna.biotype", dims=["read_class"], cells=cells))
    gene_cell: dict = {"at": [], "scalars": res["gene_scalars"]}
    if res["curve"]:
        gene_cell["curves"] = {"rna_gene_rank_curve": res["curve"]}
    if res["top_genes"]:
        t = build.counts([k for k, _ in res["top_genes"]], [v for _, v in res["top_genes"]], role="categories")
        t["other"] = int(res["gene_scalars"]["reads_with_gene"] - sum(v for _, v in res["top_genes"]))
        gene_cell["counts"] = {"rna_top_gene_table": t}
    p["fragments"].append(partial_model.fragment(None, "rna.genes", dims=[], cells=[gene_cell]))
    total = res["rows"] + res["rows_skipped"]
    fraction = round(min(res["read_bytes"] / res["size"], 1.0), 6) if res["size"] else 1.0
    if res["error"]:
        p["state"] = "truncated"
        p["errors"] = [{"code": "truncated", "message": res["error"]}]
        p["coverage"] = {"kind": "truncated", "records": res["rows"], "records_total": None, "fraction": fraction}
    elif not res["finished"]:
        p["state"] = "partially_read"
        p["coverage"] = {"kind": "byte_ranges", "records": res["rows"], "records_total": None, "fraction": fraction}
    elif not res["ends_with_newline"]:
        p["state"] = "truncated"
        p["errors"] = [{"code": "truncated", "message": "the table does not end with a line break: its last row may be cut off"}]
        p["coverage"] = {"kind": "truncated", "records": res["rows"], "records_total": None, "fraction": None}
    else:
        p["coverage"] = {"kind": "full", "records": res["rows"], "records_total": total, "fraction": 1.0}
        if total and res["rows_skipped"] / total > BAD_ROWS_FRACTION:
            p["state"] = "error"
            p["errors"] = [{"code": "bad_rows", "message": f"{res['rows_skipped']:,} of {total:,} rows could not be read"}]
    seconds = time.time() - t0
    p["timings"] = {"seconds": round(seconds, 3), "mb_per_s": round(res["read_bytes"] / 1e6 / max(seconds, 1e-6), 1), "workers": 1}
    p["facts"] = {"kind": "read_table", "records": res["rows"], "primary": res["primary"], "secondary": res["secondary"],
                  "supplementary": res["supplementary"], "rows_skipped": res["rows_skipped"],
                  "columns": res["columns"], "all_columns": res["all_columns"],
                  "tail": {"column": res["columns"]["tail"], "coding": dict(TAIL_CODING), "min": res["tail_min"], "max": res["tail_max"],
                           "fractional": res["tail_fractional"]},
                  "biotypes_seen": res["biotypes_seen"], "biotypes_folded": res["biotypes_folded"], "gene_names_kept": res["names_kept"]}
    return p
