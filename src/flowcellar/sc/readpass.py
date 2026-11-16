"""One streaming pass over the per-read table of wf-single-cell (``read_tags.tsv`` in workflow 2.x,
``<alias>.read_summary.tsv`` in 3.x; the same twelve columns, one row per read that got a barcode).

The file is 12 to 30 GB. It is cut into byte ranges; each worker process aligns its range to line boundaries
and parses it block by block, so no process ever holds more than one block of text. What comes back is small:

* the read-assignment counts: rows, rows with a barcode, with a usable UMI, assigned to a gene, to a transcript;
* reads per unit (per barcode; for Visium HD per 8 um bin, the 2 um square name turned into two integers);
* reads per molecule: how many (barcode, gene, UMI) were seen once, twice, ... (the saturation curve follows
  from it). A molecule is a 64-bit key, so the pass keeps eight bytes per molecule, not its three strings;
* the read-id sketch: the read ids whose text starts with a run of zeros, a fixed sub-sample that ties a
  file that kept only read names to the run it came from.

Everything is exact: every row is read. A pass that is stopped by a deadline says so (``complete: False``) and
its counts are then counts of the part that was read, never scaled up.
"""
from __future__ import annotations

import contextlib
import hashlib
import io
import multiprocessing as mp
import os
import queue
import time

import numpy as np

COLUMNS = ("read_id", "corrected_barcode", "corrected_umi", "gene", "transcript")
BLOCK_BYTES = 96 << 20
RANGE_BYTES = 512 << 20
MISSING = "-"
_HD_PREFIX = b"s_0"
_BASE2 = np.full(256, 4, dtype=np.uint64)
for _i, _ch in enumerate(b"ACGT"):
    _BASE2[_ch] = _i


def sketch_prefix(file_bytes: int) -> str:
    """Four zeros (1 read in 65,536) for a table of up to about 20 million rows, five (1 in 1,048,576) above:
    a few hundred ids either way. Ids that start with five zeros also start with four, so the two compare."""
    return "0000" if file_bytes < (3 << 30) else "00000"


def ranges_of(size: int, range_bytes: int = RANGE_BYTES) -> list[tuple[int, int]]:
    n = max(1, -(-size // range_bytes))
    step = -(-size // n)
    return [(i * step, min(size, (i + 1) * step)) for i in range(n)]


# ---------------------------------------------------------------------------------------------------------
# keys
# ---------------------------------------------------------------------------------------------------------

def _mix(x: np.ndarray) -> np.ndarray:
    """A 64-bit finaliser (splitmix64): spreads structured integers over the whole range."""
    x = x.astype(np.uint64, copy=True)
    x ^= x >> np.uint64(30)
    x *= np.uint64(0xBF58476D1CE4E5B9)
    x ^= x >> np.uint64(27)
    x *= np.uint64(0x94D049BB133111EB)
    x ^= x >> np.uint64(31)
    return x


def _hash_strings(values: list) -> np.ndarray:
    """64 bits of a keyed-less BLAKE2 per string: equal strings give equal numbers in every process."""
    out = np.empty(len(values), dtype=np.uint64)
    for i, s in enumerate(values):
        b = s if isinstance(s, bytes) else str(s).encode("utf-8")
        out[i] = int.from_bytes(hashlib.blake2b(b, digest_size=8).digest(), "little")
    return out


def _codes(arr):
    """(indices as int64, the distinct values as bytes) of a string array."""
    import pyarrow.compute as pc

    d = pc.dictionary_encode(arr.combine_chunks() if hasattr(arr, "combine_chunks") else arr)
    return d.indices.to_numpy(zero_copy_only=False).astype(np.int64), [v.as_py().encode("utf-8") for v in d.dictionary]


def hd_ids(names: list) -> tuple[np.ndarray, int] | None:
    """Distinct Visium HD barcode names (``s_002um_RRRRR_CCCCC-1``) -> (row * 65536 + col, bin size in um)."""
    if not names or not names[0].startswith(_HD_PREFIX) or len(names[0]) < 19:
        return None
    try:
        um = int(names[0][2:5])
        raw = np.array(names, dtype="S19").view(np.uint8).reshape(len(names), 19).astype(np.int64) - 48
    except ValueError:
        return None
    p10 = np.array([10000, 1000, 100, 10, 1], dtype=np.int64)
    r, c = raw[:, 8:13], raw[:, 14:19]
    if (r < 0).any() or (r > 9).any() or (c < 0).any() or (c > 9).any():
        return None
    return (r @ p10) * 65536 + (c @ p10), um


def _umi_codes(arr) -> np.ndarray:
    """One 64-bit number per UMI string. UMIs of one length are packed two bits per base straight from the text
    buffer; anything else is hashed value by value."""
    n = len(arr)
    if n == 0:
        return np.zeros(0, dtype=np.uint64)
    import pyarrow as pa

    a = arr.combine_chunks() if hasattr(arr, "combine_chunks") else arr
    offs_t = np.int64 if pa.types.is_large_string(a.type) else np.int32
    bufs = a.buffers()
    offs = np.frombuffer(bufs[1], dtype=offs_t, count=n + 1, offset=a.offset * np.dtype(offs_t).itemsize).astype(np.int64)
    lens = np.diff(offs)
    width = int(lens[0])
    if 0 < width <= 28 and (lens == width).all():
        data = np.frombuffer(bufs[2], dtype=np.uint8, count=n * width, offset=int(offs[0])).reshape(n, width)
        two = _BASE2[data]
        if not (two == 4).any():
            out = np.full(n, width, dtype=np.uint64)          # the length leads, so UMIs of two lengths never meet
            for k in range(width):
                out = (out << np.uint64(2)) | two[:, k]
            return out
    idx, values = _codes(a)
    return _hash_strings(values)[idx]


def _merge(keys: list, counts: list) -> tuple[np.ndarray, np.ndarray]:
    """Sum the counts of equal keys over several (keys, counts) pairs."""
    keys = [k for k in keys if len(k)]
    counts = [c for c in counts if len(c)]
    if not keys:
        return np.zeros(0, dtype=np.uint64), np.zeros(0, dtype=np.int64)
    k = np.concatenate(keys)
    c = np.concatenate(counts)
    order = np.argsort(k, kind="stable")
    k, c = k[order], c[order]
    first = np.concatenate([[True], k[1:] != k[:-1]])
    starts = np.nonzero(first)[0]
    return k[starts], np.add.reduceat(c, starts)


# ---------------------------------------------------------------------------------------------------------
# one byte range, in a worker process
# ---------------------------------------------------------------------------------------------------------

def _empty() -> dict:
    return {"rows": 0, "valid_barcode": 0, "valid_umi": 0, "gene_assigned": 0, "transcript_assigned": 0,
            "molecule_reads": 0, "bytes": 0, "hd_um": None, "unit_ids": np.zeros(0, np.int64), "unit_counts": np.zeros(0, np.int64),
            "unit_names": {}, "mol_keys": np.zeros(0, np.uint64), "mol_counts": np.zeros(0, np.int64), "sketch": []}


def _block(data: bytes, header: list[str], prefix: str, dedup: bool, acc: dict) -> None:
    import pyarrow as pa
    import pyarrow.compute as pc
    import pyarrow.csv as pacsv

    tbl = pacsv.read_csv(
        io.BytesIO(data),
        read_options=pacsv.ReadOptions(column_names=header, use_threads=False, block_size=32 << 20),
        parse_options=pacsv.ParseOptions(delimiter="\t", quote_char=False),
        convert_options=pacsv.ConvertOptions(include_columns=list(COLUMNS), strings_can_be_null=False,
                                             column_types={c: pa.string() for c in COLUMNS}))
    n = tbl.num_rows
    if n == 0:
        return
    tbl = tbl.combine_chunks()
    bc, umi, gene, tr, rid = (tbl.column(c) for c in ("corrected_barcode", "corrected_umi", "gene", "transcript", "read_id"))
    has_bc = pc.not_equal(bc, MISSING)
    has_umi = pc.and_(has_bc, pc.not_equal(umi, MISSING))
    has_gene = pc.not_equal(gene, MISSING)
    acc["rows"] += n
    acc["valid_barcode"] += int(pc.sum(has_bc).as_py() or 0)
    acc["valid_umi"] += int(pc.sum(has_umi).as_py() or 0)
    acc["gene_assigned"] += int(pc.sum(has_gene).as_py() or 0)
    acc["transcript_assigned"] += int(pc.sum(pc.not_equal(tr, MISSING)).as_py() or 0)

    # reads per unit: the distinct barcode names of the block, each with its number of rows
    vc = pc.value_counts(pc.filter(bc, has_bc))
    if hasattr(vc, "combine_chunks"):
        vc = vc.combine_chunks()
    names = [v.as_py().encode("utf-8") for v in vc.field("values")]
    cnt = vc.field("counts").to_numpy(zero_copy_only=False).astype(np.int64)
    hd = hd_ids(names)
    if hd is not None:
        ids, um = hd
        acc["hd_um"] = um
        acc["_unit_ids"].append(ids)
        acc["_unit_counts"].append(cnt)
    else:
        d = acc["unit_names"]
        for name, c in zip(names, cnt.tolist()):
            d[name] = d.get(name, 0) + c

    # reads per molecule: one 64-bit key per (barcode, gene, UMI)
    if dedup:
        m = pc.and_(has_umi, has_gene)
        nm = int(pc.sum(m).as_py() or 0)
        acc["molecule_reads"] += nm
        if nm:
            bi, bvals = _codes(pc.filter(bc, m))
            hb = hd_ids(bvals)
            bkey = (hb[0].astype(np.uint64) if hb is not None else _hash_strings(bvals))[bi]
            gi, gvals = _codes(pc.filter(gene, m))
            gkey = _hash_strings(gvals)[gi]
            ukey = _umi_codes(pc.filter(umi, m))
            key = _mix(_mix(bkey) ^ (gkey * np.uint64(0x9E3779B97F4A7C15)) ^ _mix(ukey + np.uint64(0x632BE59BD9B4E019)))
            k, c = np.unique(key, return_counts=True)
            acc["_mol_keys"].append(k)
            acc["_mol_counts"].append(c.astype(np.int64))
            if len(acc["_mol_keys"]) >= 8:                 # keep the worker's memory at a few blocks of keys
                k, c = _merge(acc["_mol_keys"], acc["_mol_counts"])
                acc["_mol_keys"], acc["_mol_counts"] = [k], [c]

    hit = pc.filter(rid, pc.starts_with(rid, prefix))
    if len(hit):
        acc["sketch"].extend(str(v.as_py()).replace("-", "")[:16] for v in hit)


def scan_range(args: tuple) -> dict:
    """Parse the lines that START inside ``[start, end)``. Module-level so that every start method can find it."""
    path, start, end, header, prefix, dedup, block_bytes, deadline = args
    acc = _empty()
    acc.update(_unit_ids=[], _unit_counts=[], _mol_keys=[], _mol_counts=[], stopped=False)
    with open(path, "rb", buffering=1 << 20) as fh:
        if start > 0:
            fh.seek(start - 1)
            fh.readline()                               # the line that started before the range belongs to it
        else:
            fh.readline()                               # the header line
        pos = fh.tell()
        carry = b""
        while pos < end:
            if deadline is not None and time.time() > deadline:
                acc["stopped"] = True
                break
            buf = fh.read(min(block_bytes, end - pos))
            if not buf:
                break
            pos += len(buf)
            if pos >= end and not buf.endswith(b"\n"):
                buf += fh.readline()                    # finish the line that crosses the end of the range
            data = carry + buf
            cut = data.rfind(b"\n")
            carry, data = data[cut + 1:], data[:cut + 1]
            acc["bytes"] += len(data)
            if data:
                _block(data, header, prefix, dedup, acc)
        if carry and not acc["stopped"]:                # a last line without a line end
            acc["bytes"] += len(carry)
            _block(carry + b"\n", header, prefix, dedup, acc)
    acc["unit_ids"], acc["unit_counts"] = _merge([a.astype(np.uint64) for a in acc.pop("_unit_ids")], acc.pop("_unit_counts"))
    acc["mol_keys"], acc["mol_counts"] = _merge(acc.pop("_mol_keys"), acc.pop("_mol_counts"))
    return acc


# ---------------------------------------------------------------------------------------------------------
# the whole file
# ---------------------------------------------------------------------------------------------------------

def _context():
    return mp.get_context("forkserver" if "forkserver" in mp.get_all_start_methods() else "spawn")


@contextlib.contextmanager
def _without_main():
    """Start worker processes that do not import the caller's main module. A worker here needs nothing of it (its
    task is a function of this module), and importing it again would run a script without a main guard and fails
    for a program read from standard input."""
    from multiprocessing import spawn

    original = spawn.get_preparation_data

    def prepared(name):
        d = original(name)
        d.pop("init_main_from_name", None)
        d.pop("init_main_from_path", None)
        return d
    spawn.get_preparation_data = prepared
    try:
        yield
    finally:
        spawn.get_preparation_data = original


def _worker(jobs_q, out_q) -> None:
    while True:
        item = jobs_q.get()
        if item is None:
            return
        i, args = item
        try:
            out_q.put((i, scan_range(args), None))
        except BaseException as e:  # noqa: BLE001  the parent raises it with the range it belongs to
            out_q.put((i, None, f"{type(e).__name__}: {e}"))


def _in_workers(jobs: list, workers: int, progress) -> list[dict]:
    """``scan_range`` for every job, in ``workers`` processes. All reading of the table happens there."""
    from flowcellar.util.errors import InputError

    ctx = _context()
    jobs_q, out_q = ctx.Queue(), ctx.Queue()
    for item in enumerate(jobs):
        jobs_q.put(item)
    procs = [ctx.Process(target=_worker, args=(jobs_q, out_q), daemon=True) for _ in range(workers)]
    for _ in procs:
        jobs_q.put(None)
    with _without_main():
        for p in procs:
            p.start()
    results: list = [None] * len(jobs)
    got = 0
    try:
        while got < len(jobs):
            try:
                i, res, err = out_q.get(timeout=1.0)
            except queue.Empty:
                if not any(p.is_alive() for p in procs) and out_q.empty():
                    raise InputError("the pass over the per-read table stopped: a worker process ended early") from None
                continue
            if err:
                raise InputError(f"the per-read table could not be read ({err})")
            results[i] = res
            got += 1
            if progress:
                progress(got / len(jobs))
    finally:
        for p in procs:
            if p.is_alive() and got < len(jobs):
                p.terminate()
        for p in procs:
            p.join(timeout=5)
    return results


def read_pass(path: str, *, workers: int = 0, dedup: bool = True, deadline: float | None = None, progress=None,
              range_bytes: int = RANGE_BYTES, block_bytes: int = BLOCK_BYTES, inprocess: bool = False) -> dict:
    """Read the whole per-read table. ``deadline`` is a ``time.time()`` value; ``progress(fraction)`` is called
    as ranges finish. Returns plain numbers and numpy arrays (see the module text)."""
    size = os.path.getsize(path)
    with open(path, "rb") as fh:
        header = fh.readline().decode("utf-8", "replace").rstrip("\r\n").split("\t")
    missing = [c for c in COLUMNS if c not in header]
    if missing:
        from flowcellar.util.errors import InputError
        raise InputError(f"{path}: not the per-read table of wf-single-cell (no column {', '.join(missing)})")
    prefix = sketch_prefix(size)
    ranges = ranges_of(size, range_bytes)
    jobs = [(path, a, b, header, prefix, dedup, block_bytes, deadline) for a, b in ranges]
    workers = workers or min(24, os.cpu_count() or 1)
    results: list[dict] = []
    t0 = time.time()
    n_workers = max(1, min(workers, len(jobs)))
    if inprocess:
        for j in jobs:
            results.append(scan_range(j))
            if progress:
                progress(len(results) / len(jobs))
    else:
        results = _in_workers(jobs, n_workers, progress)
    out = {k: sum(r[k] for r in results) for k in ("rows", "valid_barcode", "valid_umi", "gene_assigned", "transcript_assigned",
                                                    "molecule_reads", "bytes")}
    out["complete"] = not any(r.get("stopped") for r in results) and out["bytes"] + _header_bytes(path) == size
    out["file_bytes"] = size
    out["seconds"] = round(time.time() - t0, 2)
    out["workers"] = 1 if inprocess else n_workers
    hd_um = next((r["hd_um"] for r in results if r["hd_um"]), None)
    out["hd_um"] = hd_um
    ids, counts = _merge([r["unit_ids"] for r in results], [r["unit_counts"] for r in results])
    names: dict = {}
    for r in results:
        for k, v in r["unit_names"].items():
            names[k] = names.get(k, 0) + v
    out["unit_ids"], out["unit_id_counts"] = ids.astype(np.int64), counts
    out["unit_names"] = names
    keys, mcounts = _merge([r["mol_keys"] for r in results], [r["mol_counts"] for r in results])
    out["molecules"] = int(len(keys)) if dedup else None
    out["molecule_hist"] = np.bincount(mcounts).astype(np.int64) if (dedup and len(mcounts)) else (np.zeros(1, np.int64) if dedup else None)
    out["sketch"] = {"prefix": prefix, "ids": sorted(set(x for r in results for x in r["sketch"])), "rows_seen": out["rows"]}
    return out


def _header_bytes(path: str) -> int:
    with open(path, "rb") as fh:
        return len(fh.readline())


def reads_per_unit(result: dict, *, unit_um: int = 8) -> np.ndarray:
    """Reads of every unit that has at least one, from the result of a pass. Visium HD squares are summed into
    bins of ``unit_um``; any other barcode is its own unit."""
    counts = [np.asarray(list(result["unit_names"].values()), dtype=np.int64)] if result["unit_names"] else []
    ids = result["unit_ids"]
    if len(ids):
        um = int(result.get("hd_um") or unit_um)
        f = max(1, unit_um // um)
        rows, cols = ids // 65536, ids % 65536
        key = (rows // f) * 65536 + (cols // f)
        order = np.argsort(key, kind="stable")
        key, c = key[order], result["unit_id_counts"][order]
        starts = np.nonzero(np.concatenate([[True], key[1:] != key[:-1]]))[0]
        counts.append(np.add.reduceat(c, starts))
    return np.concatenate(counts) if counts else np.zeros(0, dtype=np.int64)


def hd_grid(result: dict, *, unit_um: int = 8) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(row, col, reads) of every ``unit_um`` bin with a read, for the map of an analysis that has no Space
    Ranger sibling."""
    ids = result["unit_ids"]
    if not len(ids):
        z = np.zeros(0, np.int64)
        return z, z, z
    um = int(result.get("hd_um") or unit_um)
    f = max(1, unit_um // um)
    key = ((ids // 65536) // f) * 65536 + ((ids % 65536) // f)
    order = np.argsort(key, kind="stable")
    key, c = key[order], result["unit_id_counts"][order]
    starts = np.nonzero(np.concatenate([[True], key[1:] != key[:-1]]))[0]
    return key[starts] // 65536, key[starts] % 65536, np.add.reduceat(c, starts)
