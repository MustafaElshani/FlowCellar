"""Reading a BAM index (.bai or .csi) without htslib: where records start, and how many there are.

The index is the only cheap source of two things the read pass needs:

* **cut points**: every chunk start and every linear-index entry is the virtual offset of the first byte of a
  record, so the file can be split into tasks of about equal compressed size. A task is a half-open range of
  virtual offsets; a record belongs to the task that holds its first byte, so no record is read twice and none
  is missed, whatever the number of tasks;
* **totals**: mapped and unmapped records per reference (the pseudo-bin 37450) and the records without a
  coordinate. A complete pass must reproduce their sum.

A virtual offset is ``compressed_offset << 16 | offset_inside_the_block``.
"""
from __future__ import annotations

import gzip
import os
import struct
from dataclasses import dataclass, field

import numpy as np

PSEUDO_BIN = 37450
MAX_INDEX_BYTES = 512 << 20


@dataclass
class BamIndex:
    path: str
    kind: str                                   # "bai" or "csi"
    n_ref: int
    offsets: np.ndarray                         # sorted unique virtual offsets of record starts (uint64)
    mapped: int | None                          # records with a coordinate that are mapped; None when not stored
    unmapped_placed: int | None                 # unmapped records placed on a reference
    no_coordinate: int | None                   # records without a coordinate (the unmapped tail)
    notes: list[str] = field(default_factory=list)

    @property
    def records(self) -> int | None:
        """Every record of the file according to the index (primary, secondary and supplementary)."""
        if self.mapped is None or self.no_coordinate is None:
            return None
        return self.mapped + (self.unmapped_placed or 0) + self.no_coordinate


def find_index(bam_path: str) -> str | None:
    """The index beside a BAM: ``x.bam.bai``, ``x.bai``, ``x.bam.csi`` or ``x.csi``."""
    stem = bam_path[:-4] if bam_path.endswith(".bam") else bam_path
    for cand in (bam_path + ".bai", stem + ".bai", bam_path + ".csi", stem + ".csi"):
        if os.path.isfile(cand):
            return cand
    return None


def index_is_stale(bam_path: str, index_path: str) -> bool:
    """An index older than its BAM cannot be trusted to describe it."""
    try:
        return os.stat(index_path).st_mtime_ns < os.stat(bam_path).st_mtime_ns
    except OSError:
        return True


def _read_all(path: str) -> bytes:
    if os.path.getsize(path) > MAX_INDEX_BYTES:
        raise ValueError(f"{path}: index larger than {MAX_INDEX_BYTES >> 20} MB")
    with open(path, "rb") as fh:
        data = fh.read()
    if data[:2] == b"\x1f\x8b":                 # a .csi is BGZF-compressed, which is a valid multi-member gzip
        data = gzip.decompress(data)
    return data


def read_index(index_path: str) -> BamIndex:
    """Parse a .bai or .csi. Raises ValueError when the file is not an index or is cut short."""
    buf = _read_all(index_path)
    magic = buf[:4]
    if magic == b"BAI\x01":
        kind, pos, pseudo = "bai", 4, PSEUDO_BIN
    elif magic == b"CSI\x01":
        kind = "csi"
        _min_shift, depth, l_aux = struct.unpack_from("<iii", buf, 4)
        pos = 16 + l_aux
        pseudo = ((1 << (depth * 3 + 3)) - 1) // 7 + 1      # the bin above the last real one
    else:
        raise ValueError(f"{index_path}: not a BAM index")
    try:
        (n_ref,) = struct.unpack_from("<i", buf, pos)
        pos += 4
        parts: list[np.ndarray] = []
        mapped = unmapped = 0
        have_counts = False
        for _ in range(n_ref):
            (n_bin,) = struct.unpack_from("<i", buf, pos)
            pos += 4
            for _ in range(n_bin):
                if kind == "bai":
                    bin_id, n_chunk = struct.unpack_from("<Ii", buf, pos)
                    pos += 8
                else:
                    bin_id, loffset, n_chunk = struct.unpack_from("<IQi", buf, pos)
                    pos += 16
                    if loffset:
                        parts.append(np.array([loffset], dtype=np.uint64))
                chunks = np.frombuffer(buf, dtype="<u8", count=2 * n_chunk, offset=pos)
                pos += 16 * n_chunk
                # the pseudo-bin holds (first offset, last offset) and (mapped, unmapped) of the reference
                if bin_id == pseudo and n_chunk == 2:
                    parts.append(chunks[:2])
                    mapped += int(chunks[2])
                    unmapped += int(chunks[3])
                    have_counts = True
                elif n_chunk:
                    parts.append(chunks[0::2])
            if kind == "bai":
                (n_intv,) = struct.unpack_from("<i", buf, pos)
                pos += 4
                parts.append(np.frombuffer(buf, dtype="<u8", count=n_intv, offset=pos))
                pos += 8 * n_intv
        no_coor = struct.unpack_from("<Q", buf, pos)[0] if pos + 8 <= len(buf) else None
    except (struct.error, ValueError) as e:
        raise ValueError(f"{index_path}: index is cut short or damaged ({e})") from None
    offsets = np.unique(np.concatenate(parts)) if parts else np.zeros(0, np.uint64)
    offsets = offsets[offsets > 0]
    return BamIndex(index_path, kind, n_ref, offsets, mapped if have_counts else None,
                    unmapped if have_counts else None, int(no_coor) if no_coor is not None else None)


def cut_tasks(first_record: int, offsets: np.ndarray, file_size: int, target_bytes: int,
              max_tasks: int = 4096) -> list[tuple[int, int | None]]:
    """Half-open ranges of virtual offsets ``[(start, end), ..., (start, None)]`` that cover the file from its first
    record to its end, each of about ``target_bytes`` compressed bytes. Cuts fall only on record starts known to the
    index, so a region the index cannot split (one deep locus, the unmapped tail) stays one task.
    """
    target = max(int(target_bytes), 1)
    if max_tasks and file_size / target > max_tasks:
        target = int(file_size / max_tasks) + 1
    offs = offsets[offsets > first_record]
    cuts: list[int] = []
    last_c = first_record >> 16
    if len(offs):
        comp = (offs >> np.uint64(16)).astype(np.int64)
        while True:         # the first candidate at least one target beyond the previous cut
            j = int(np.searchsorted(comp, last_c + target, side="left"))
            if j >= len(offs):
                break
            cuts.append(int(offs[j]))
            last_c = int(comp[j])
    bounds = [first_record] + cuts
    return [(bounds[k], bounds[k + 1] if k + 1 < len(bounds) else None) for k in range(len(bounds))]
