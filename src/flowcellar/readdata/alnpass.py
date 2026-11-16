"""From the merged accumulators of a BAM pass to the fragments of its partial: the alignment, over-time and poly(A)
parts (``alnacc``). Called by ``passes._bam_pass`` for every group (one run key and one sample label).

What a group adds to the partial:

* members of its ``aln`` cell (strand, CIGAR totals, mismatches, aligned length, soft-clipped share, the three grids,
  the contig table and the depth windows of a genome) and, in the fragment's ``facts``, the reference the file was
  aligned to and aligned length and identity by clock hour;
* a fragment ``tx`` for a transcriptome: coverage along the transcript body, and in ``facts`` the read count of every
  transcript (and gene) that has one, so that reads per transcript and the saturation curve can be worked out exactly
  for one sample and for several added up;
* a fragment ``time``: reads and bases with a start time, the measured speed, and in ``facts`` the clock tables;
* a fragment ``polya``: the poly(A) tail lengths of the ``pt`` tag.
"""
from __future__ import annotations

import numpy as np

from flowcellar.model import partial as partial_model
from flowcellar.readdata.alnacc import sparse_counts
from flowcellar.readdata.reference import Layout


class _Described:
    """The reference layout as the plan described it: enough to write the cells (names, kind, number of windows)."""

    def __init__(self, d: dict):
        self.kind = d.get("kind")
        self.names = [c["name"] for c in d.get("contigs") or []]
        self.n_windows = int(d.get("windows") or 0)


def reference_summary(d: dict | None) -> dict | None:
    """What the facts of a partial keep about the reference (the contig list is in the ``aln`` fragments)."""
    if not d:
        return None
    return {k: d.get(k) for k in ("kind", "evidence", "id", "sequences", "total_length", "window_bases", "windows", "mitochondrial")
            if k in d}


def extend(p: dict, cells: dict, merged, gid: tuple, plan: dict) -> None:
    """Add the members of ``alnacc`` to ``cells`` (module -> cell) of one group and append its extra fragments.
    The fragments of ``cells`` themselves are appended by the caller afterwards."""
    key, label = gid
    ref = plan.get("reference") or {}
    lay: Layout | _Described | None = _Described(ref) if ref.get("kind") else None
    acc = merged.groups[gid]
    aln = merged.aln.get(gid)
    if aln is not None and "aln" in cells:
        extra = aln.cell(lay, unmapped=int(acc.unmapped))
        for kind, members in extra.items():
            cells["aln"].setdefault(kind, {}).update(members)
        if lay is not None and lay.kind == "transcriptome":
            tx = aln.tx_cell()
            if tx is not None:
                facts: dict = {"reference": {k: ref.get(k) for k in ("kind", "id", "sequences", "total_length", "length_classes")},
                               "transcripts": None, "genes": None}
                if aln.ref_n is not None and not aln.ref_counts_dropped:
                    facts["transcripts"] = sparse_counts(aln.ref_n)
                    genes = plan.get("genes")
                    if genes is not None:
                        gene_of, n_genes = genes
                        per_gene = np.bincount(gene_of, weights=aln.ref_n, minlength=n_genes).astype(np.int64)
                        facts["genes"] = sparse_counts(per_gene)
                p["fragments"].append(partial_model.fragment(key, "tx", sample=label, dims=[], cells=[{"at": [], **tx}], facts=facts))
    tags = merged.tagacc.get(gid)
    if tags is not None:
        tc = tags.time_cell()
        if tc is not None:
            p["fragments"].append(partial_model.fragment(key, "time", sample=label, dims=[], cells=[{"at": [], **tc}],
                                                         facts={"time": tags.time_json()}))
        pc = tags.polya_cell()
        if pc is not None:
            p["fragments"].append(partial_model.fragment(key, "polya", sample=label, dims=[], cells=[{"at": [], **pc}]))


def finish(p: dict, merged, plan: dict) -> None:
    """After the fragments of every group were appended: put the reference and the clock tables of the aligned
    reads into the ``facts`` of the ``aln`` fragments."""
    ref = plan.get("reference") or {}
    for f in p["fragments"]:
        if f.get("module") != "aln":
            continue
        aln = merged.aln.get((f.get("run_key"), f.get("sample")))
        f["facts"] = {"reference": ref if ref.get("kind") else None, "time": aln.time_json() if aln is not None else None}
