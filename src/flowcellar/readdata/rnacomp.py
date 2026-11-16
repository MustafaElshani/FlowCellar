"""RNA composition of a run, from the partials of its annotated read tables (``readtable``).

Two sections, both collected from a pipeline's own result files (they name their files, not a read population):

``rna.biotype``  dims ``sample`` x ``read_class``. A read class is a transcript biotype as the table writes it; the
                 most common biotypes of the run are classes of their own, the rest form the class ``other``. The
                 list of classes is one list for the whole run, so cells add up class by class and sample by sample.
                 A cell holds the reads of the class, how many carry a tail estimate (and the three ways of having
                 none), the histogram of the estimated tail lengths, and the reads of every biotype it is made of.
                 ``run_level.tables.rna_polya_by_biotype`` holds exact quantiles per sample and class and for the
                 samples together (sample ``*``).
``rna.genes``    dims ``sample``. Reads per gene summarised. The cells hold what adds up over samples (reads, reads
                 that name a gene, reads of mitochondrial genes) and each sample's rank curve. What does not add up,
                 because the genes of two samples overlap, is kept as one row per sample in ``run_level.tables``:
                 ``rna_genes_detected`` (genes with at least 1, 10 and 100 reads, transcripts with a read, the reads
                 of the top 1, 10 and 100 genes) and ``rna_top_gene_table`` (the short list of top genes by name,
                 when the catalogue allows names).
"""
from __future__ import annotations

import numpy as np

from flowcellar.aggregates import build
from flowcellar.model import axes
from flowcellar.readdata import readtable
from flowcellar.util.text import clean_text

VERSION = readtable.VERSION
SECTION_BIOTYPE = "rna.biotype"
SECTION_GENES = "rna.genes"
MAX_CLASSES = 8                         # biotypes that are a class of their own
MIN_CLASS_SHARE = 0.001                 # ... when they hold at least this share of the run's reads
OTHER = "other"
ALL = "all"
POOLED = "*"
PARTIAL_AXIS = readtable.TAIL_AXIS      # 1 nt
BUNDLE_AXIS = "polya_5nt"
RRNA_BIOTYPES = ("rRNA", "rRNA_pseudogene", "Mt_rRNA")
BIOTYPES_LISTED = 40                    # biotypes of the class "other" stored by name per cell; the rest is one count
GENE_SUMS = ("reads", "reads_with_gene", "reads_no_gene", "reads_mito_genes")                 # these add up over samples
GENE_PER_SAMPLE = ("genes_ge1", "genes_ge10", "genes_ge100", "transcripts_ge1", "reads_top1", "reads_top10", "reads_top100")   # these do not
TAIL_KEYS = ("reads", "tail_estimated", "tail_zero", "tail_not_found", "tail_no_value", "tail_sum")


def _fragment(partial: dict, module: str) -> dict | None:
    return next((f for f in partial.get("fragments", []) if f.get("module") == module), None)


def _full(cell: dict) -> np.ndarray:
    h = (cell.get("hists") or {}).get("rna_biotype_polya_hist")
    return build.hist_full(h) if h else np.zeros(axes.n_cells(PARTIAL_AXIS), np.int64)


def quantile(full: np.ndarray, q: float) -> float | None:
    """A quantile of whole-number tail lengths from the 1 nt histogram: the value of the bin that holds the rank
    ``ceil(q x n)`` (for the median of an even count: the lower of the two middle values). Exact for whole numbers
    under 1,000 nt; a rank that falls among the longer tails gives 1000."""
    n = int(full.sum())
    if not n:
        return None
    rank = max(1, int(np.ceil(q * n)))
    k = int(np.searchsorted(np.cumsum(full), rank, side="left"))
    edges = axes.EDGES[PARTIAL_AXIS]
    return float(edges[min(max(k - 1, 0), len(edges) - 1)])


def classes_of(totals: dict[str, int]) -> tuple[list[str], dict[str, str]]:
    """(the classes of the run, biotype -> class): the biotypes with the most reads, at most ``MAX_CLASSES`` and each
    at least ``MIN_CLASS_SHARE`` of the reads, then ``other`` when anything is left."""
    total = sum(totals.values())
    ranked = sorted(totals, key=lambda b: (-totals[b], b))
    own = [b for b in ranked if b != OTHER and total and totals[b] >= MIN_CLASS_SHARE * total][:MAX_CLASSES]
    cls = {b: (b if b in own else OTHER) for b in ranked}
    return own + ([OTHER] if any(v == OTHER for v in cls.values()) else []), cls


def _class_cell(members: list[dict], names_order: list[str]) -> tuple[dict, np.ndarray]:
    """The cell of one sample and class from the partial cells of its biotypes (without ``at``), and the 1 nt array."""
    sc = {k: 0 for k in TAIL_KEYS}
    full = np.zeros(axes.n_cells(PARTIAL_AXIS), np.int64)
    by_name: dict[str, tuple[int, int]] = {}
    for c in members:
        s = c.get("scalars") or {}
        for k in TAIL_KEYS:
            sc[k] += s.get(k) or 0
        full += _full(c)
        by_name[c["at"][0]] = (int(s.get("reads") or 0), int(s.get("tail_estimated") or 0))
    names = [b for b in names_order if b in by_name]
    listed, rest = names[:BIOTYPES_LISTED], names[BIOTYPES_LISTED:]
    table = build.counts(listed, [by_name[b][0] for b in listed], w=[by_name[b][1] for b in listed], role="categories")
    table["w_name"] = "reads with a tail estimate"
    if rest:
        table["other"] = sum(by_name[b][0] for b in rest)
    coarse = axes.coarsen(full, axes.coarsen_factor(PARTIAL_AXIS, BUNDLE_AXIS))
    hist = build.hist(BUNDLE_AXIS, coarse, missing=int(sc["reads"] - sc["tail_estimated"]))
    cell = {"scalars": sc, "counts": {"rna_biotype_table": table}, "hists": {"rna_biotype_polya_hist": hist}}
    return cell, full


def _stats_row(sample: str, cls: str, reads: int, full: np.ndarray, tail_sum) -> list:
    n = int(full.sum())
    return [sample, cls, int(reads), n, quantile(full, 0.5), quantile(full, 0.25), quantile(full, 0.75),
            round(float(tail_sum) / n, 2) if n else None]


def _state(entries: list[dict], waiting: list[dict], expected: list[str], conflicts: list[str], records: int) -> tuple[str, str | None, dict]:
    cov = {"kind": "full", "records": int(records), "records_total": None, "fraction": 1.0}
    part = [e for e in entries if e["state"] in ("partially_read", "truncated")]
    if part:
        fr = [(e["partial"].get("coverage") or {}).get("fraction") for e in part]
        fr = [f for f in fr if f is not None]
        cov.update(kind="truncated" if any(e["state"] == "truncated" for e in part) else "byte_ranges", fraction=min(fr) if fr else None)
        return "partial", "partly_read", cov
    if waiting:
        total = sum(int(e["size"] or 0) for e in entries) + sum(int(w["size"] or 0) for w in waiting)
        read = sum(int(e["size"] or 0) for e in entries)
        cov.update(kind="byte_ranges", fraction=round(read / total, 6) if total else None)
        return "partial", "waiting", cov
    have = {e["key"] for e in entries}
    if conflicts or any(k not in have for k in expected):
        return "partial", "samples_missing", cov
    return "ok", None, cov


def _common_notes(entries: list[dict], waiting: list[dict], conflicts: list[str]) -> list[dict]:
    notes = [{"level": "info", "code": "source_is_report",
              "message": "Counted from each sample's annotated read table: one row per read that the pipeline gave to a transcript."}]
    for e in entries:
        f = e["partial"].get("facts") or {}
        if e.get("several"):
            notes.append({"level": "info", "code": "covers_several_runs", "input": e["input"], "value": e["several"],
                          "message": clean_text(f"{e['name']}: one read table for {e['several']} flow cells; their reads were merged before the table was made.", 600)})
        if e.get("weak"):
            notes.append({"level": "warning", "code": "weak_link", "input": e["input"],
                          "message": clean_text(f"{e['name']}: the table is linked by a rule, not by what it contains.", 600)})
        if e["state"] in ("partially_read", "truncated"):
            why = "; ".join(x.get("message", "") for x in e["partial"].get("errors") or []) or "the pass stopped before the end of the file"
            notes.append({"level": "warning", "code": "partly_read", "input": e["input"],
                          "message": clean_text(f"{e['name']}: the table was not read to its end ({why}); its numbers are a lower bound.", 600)})
        if f.get("rows_skipped"):
            notes.append({"level": "warning", "code": "rows_skipped", "input": e["input"], "value": f["rows_skipped"],
                          "message": clean_text(f"{e['name']}: {f['rows_skipped']:,} rows with the wrong number of columns were left out.", 600)})
        extra = (f.get("secondary") or 0) + (f.get("supplementary") or 0)
        if extra:
            notes.append({"level": "info", "code": "rows_not_primary", "input": e["input"], "value": extra,
                          "message": clean_text(f"{e['name']}: {extra:,} rows are secondary or supplementary alignments; they are counted apart and used for nothing.", 600)})
    if waiting:
        notes.append({"level": "info", "code": "waiting", "value": len(waiting),
                      "message": f"{len(waiting)} annotated read table{'s' if len(waiting) != 1 else ''} of this run "
                                 f"{'have' if len(waiting) != 1 else 'has'} not been read yet."})
    for m in conflicts:
        notes.append({"level": "warning", "code": "two_candidate_sources", "message": clean_text("Not assessed: " + m, 1200)})
    return notes


def build_sections(entries: list[dict], *, order: list[str], expected: list[str], waiting: list[dict] = (), conflicts: list[str] = ()) -> dict[str, dict]:
    """``{"rna.biotype": section, "rna.genes": section}`` (empty when the run has no annotated read table).

    ``entries``: one per sample whose table was read, ``{"key", "name", "input", "partial", "state", "size",
    "several", "weak"}``. ``waiting``: the tables not read yet (``{"input", "size"}``). ``order``: the sample keys of
    the bundle in order; ``expected``: the keys of the samples that should have a table."""
    waiting, conflicts = list(waiting), list(conflicts)
    if not entries and not waiting:
        return {}
    pos = {k: i for i, k in enumerate(order)}
    entries = sorted(entries, key=lambda e: (pos.get(e["key"], 1 << 30), e["key"]))
    inputs = [e["input"] for e in entries] + [w["input"] for w in waiting]
    notes = _common_notes(entries, waiting, conflicts)

    # ---------------------------------------------------------------------------------- rna.biotype
    per_sample: list[tuple[dict, list[dict]]] = []
    totals: dict[str, int] = {}
    for e in entries:
        frag = _fragment(e["partial"], "rna.biotype")
        cells = (frag or {}).get("cells") or []
        per_sample.append((e, cells))
        for c in cells:
            totals[c["at"][0]] = totals.get(c["at"][0], 0) + int((c.get("scalars") or {}).get("reads") or 0)
    classes, cls_of = classes_of(totals)
    names_order = sorted(totals, key=lambda b: (-totals[b], b))
    cells_out: list[dict] = []
    stats_rows: list[list] = []
    pooled_full = {c: np.zeros(axes.n_cells(PARTIAL_AXIS), np.int64) for c in classes + [ALL]}
    pooled_sc = {c: {k: 0 for k in TAIL_KEYS} for c in classes + [ALL]}
    for e, cells in per_sample:
        all_full = np.zeros(axes.n_cells(PARTIAL_AXIS), np.int64)
        all_sc = {k: 0 for k in TAIL_KEYS}
        rows_here = []
        for cls in classes:
            members = [c for c in cells if cls_of.get(c["at"][0]) == cls]
            if not members:
                continue
            cell, full = _class_cell(members, names_order)
            if not cell["scalars"]["reads"]:
                continue
            cells_out.append({"at": [e["key"], cls], **cell})
            rows_here.append(_stats_row(e["key"], cls, cell["scalars"]["reads"], full, cell["scalars"]["tail_sum"]))
            all_full += full
            pooled_full[cls] += full
            for k in TAIL_KEYS:
                all_sc[k] += cell["scalars"][k]
                pooled_sc[cls][k] += cell["scalars"][k]
        stats_rows.append(_stats_row(e["key"], ALL, all_sc["reads"], all_full, all_sc["tail_sum"]))
        stats_rows.extend(rows_here)
        pooled_full[ALL] += all_full
        for k in TAIL_KEYS:
            pooled_sc[ALL][k] += all_sc[k]
    pooled_rows = [_stats_row(POOLED, c, pooled_sc[c]["reads"], pooled_full[c], pooled_sc[c]["tail_sum"]) for c in [ALL] + classes
                   if pooled_sc[c]["reads"] or c == ALL]
    tot = pooled_sc[ALL]
    reads = int(tot["reads"])
    has_tail = any(((e["partial"].get("facts") or {}).get("tail") or {}).get("column") for e in entries)
    fractional = any(((e["partial"].get("facts") or {}).get("tail") or {}).get("fractional") for e in entries)
    b_notes = list(notes)
    no_tail = [e for e in entries if not ((e["partial"].get("facts") or {}).get("tail") or {}).get("column")]
    if no_tail:
        b_notes.append({"level": "info", "code": "no_tail_column", "value": len(no_tail),
                        "message": clean_text(f"{len(no_tail)} of the {len(entries)} tables hold no tail length column ({no_tail[0]['name']}"
                                              f"{', ...' if len(no_tail) > 1 else ''}): their reads are counted as having no value.", 600)})
    folded = sum(int((e["partial"].get("facts") or {}).get("biotypes_folded") or 0) for e in entries)
    if folded:
        b_notes.append({"level": "info", "code": "biotypes_folded", "value": folded,
                        "message": f"A table names more than {readtable.MAX_BIOTYPES_KEPT} transcript types; the rarest are counted together as '{readtable.FOLDED}'."})
    state, reason, cov = _state(entries, waiting, expected, conflicts, reads)
    rrna = sum(v for b, v in totals.items() if b in RRNA_BIOTYPES)
    table = {"columns": [{"id": "sample", "type": "str"}, {"id": "read_class", "type": "str"}, {"id": "reads", "type": "int"},
                         {"id": "with_tail", "type": "int"}, {"id": "median_nt", "type": "float", "unit": "nt"},
                         {"id": "q25_nt", "type": "float", "unit": "nt"}, {"id": "q75_nt", "type": "float", "unit": "nt"},
                         {"id": "mean_nt", "type": "float", "unit": "nt"}],
             "key": ["sample", "read_class"], "merge": "replace", "rows": pooled_rows + stats_rows}
    biotype = {
        "version": VERSION, "tier": "assay", "state": state if cells_out or state != "ok" else "skipped", "inputs": inputs,
        "dims": ["sample", "read_class"], "cells": cells_out, "run_level": {"tables": {"rna_polya_by_biotype": table}},
        "summary": {
            "rna_biotype_samples": len(entries), "rna_biotype_reads": reads, "rna_biotypes_seen": len(totals),
            "rna_protein_coding_frac": round(totals.get("protein_coding", 0) / reads, 6) if reads else None,
            "rna_rrna_frac": round(rrna / reads, 6) if reads else None,
            "rna_biotype_polya_frac": round(tot["tail_estimated"] / reads, 6) if reads and has_tail else None,
            "rna_biotype_polya_median_nt": quantile(pooled_full[ALL], 0.5) if has_tail else None,
        },
        "params": {
            "classes": classes, "class_rule": f"a biotype is a class of its own when it holds at least {MIN_CLASS_SHARE:.1%} of the run's reads "
                                              f"(at most {MAX_CLASSES}); the others form the class '{OTHER}'",
            "counts_of": {
                "reads": "primary rows of the annotated read table: reads given to a transcript; a read that was cut into pieces counts once per piece",
                "tail_estimated": "reads with a tail length above 0",
                "tail_zero": "reads with the value 0: the basecaller found the primer but no tail length",
                "tail_not_found": "reads with a negative value: the basecaller did not find the primer",
                "tail_no_value": "reads with an empty cell or NA, or from a table without a tail column",
                "tail_sum": "sum of the estimated tail lengths, in nt (the mean is tail_sum / tail_estimated)",
                "rna_biotype_table": "reads of each transcript biotype in the class, as the table writes it; the weight is reads with a tail estimate",
                "rna_biotype_polya_hist": "reads with a tail estimate by tail length, in 5 nt steps; missing = reads without one",
                "rna_polya_by_biotype": "per sample and class, and for the samples together (sample *): reads, reads with a tail estimate, and the quartiles and mean of the tail length",
            },
            "rrna_biotypes": list(RRNA_BIOTYPES),
            "quantile_rule": "quartiles in run_level are read from a 1 nt histogram: the value at rank ceil(q x n)"
                             + ("; tail lengths with a fraction were cut to whole nt first" if fractional else "; exact for whole-number tail lengths under 1,000 nt"),
            "tail_coding": dict(readtable.TAIL_CODING),
        },
        "notes": b_notes, "detail": "standard", "coverage": cov,
    }
    if biotype["state"] == "partial":
        biotype["partial_reason"] = reason

    # ---------------------------------------------------------------------------------- rna.genes
    g_cells: list[dict] = []
    detected_rows: list[list] = []
    top_rows: list[list] = []
    g_reads = g_mito = g_with = 0
    ge10: list[int] = []
    top10: list[float] = []
    names_kept = 0
    for e in entries:
        frag = _fragment(e["partial"], "rna.genes")
        src = ((frag or {}).get("cells") or [None])[0]
        sc = (src or {}).get("scalars") or {}
        if not sc.get("reads"):
            continue
        cell: dict = {"at": [e["key"]], "scalars": {k: int(sc.get(k) or 0) for k in GENE_SUMS}}
        if src.get("curves"):
            cell["curves"] = src["curves"]
        g_cells.append(cell)
        detected_rows.append([e["key"], int(sc["reads"])] + [sc.get(k) for k in GENE_PER_SAMPLE])
        t = (src.get("counts") or {}).get("rna_top_gene_table")
        if t:
            top_rows.extend([e["key"], rank + 1, name, int(n)] for rank, (name, n) in enumerate(zip(t["cats"], t["n"])))
        g_reads += sc.get("reads") or 0
        g_with += sc.get("reads_with_gene") or 0
        g_mito += sc.get("reads_mito_genes") or 0
        ge10.append(int(sc.get("genes_ge10") or 0))
        if sc.get("reads_with_gene"):
            top10.append(sc.get("reads_top10", 0) / sc["reads_with_gene"])
        names_kept = max(names_kept, int((e["partial"].get("facts") or {}).get("gene_names_kept") or 0))
    g_tables = {"rna_genes_detected": {"columns": [{"id": "sample", "type": "str"}, {"id": "reads", "type": "int"}] + [{"id": k, "type": "int"} for k in GENE_PER_SAMPLE],
                                       "key": ["sample"], "merge": "replace", "rows": detected_rows}}
    if top_rows:
        g_tables["rna_top_gene_table"] = {"columns": [{"id": "sample", "type": "str"}, {"id": "rank", "type": "int"}, {"id": "gene", "type": "str"},
                                                      {"id": "reads", "type": "int"}],
                                          "key": ["sample", "rank"], "merge": "replace", "rows": top_rows}
    g_state, g_reason, g_cov = _state(entries, waiting, expected, conflicts, g_reads)
    genes = {
        "version": VERSION, "tier": "assay", "state": g_state if g_cells or g_state != "ok" else "skipped", "inputs": inputs,
        "dims": ["sample"], "cells": g_cells, "run_level": {"tables": g_tables},
        "summary": {
            "rna_genes_samples": len(g_cells),
            "rna_genes_detected_min": min(ge10) if ge10 else None, "rna_genes_detected_max": max(ge10) if ge10 else None,
            "rna_top_gene_max_frac": round(max(top10), 6) if top10 else None,
            "rna_mito_frac": round(g_mito / g_with, 6) if g_with else None,
        },
        "params": {
            "counts_of": {
                "reads": "primary rows of the annotated read table",
                "reads_with_gene": "of those, the rows that name a gene",
                "reads_no_gene": "of those, the rows that name no gene",
                "reads_mito_genes": "reads of genes whose name starts with MT-",
                "rna_gene_rank_curve": "the genes of one sample in order of their reads: reads of the gene at each rank, and reads of all genes up to it",
                "rna_genes_detected": "one row per sample: genes with at least 1, 10 and 100 reads, transcripts with at least one read, "
                                      "and the reads of the 1, 10 and 100 genes with the most reads; numbers of one sample, never added up",
                "rna_top_gene_table": "one row per sample and rank: the genes with the most reads in that sample, by name",
            },
            "gene_key": "gene_id, or gene_name in a table without an id column",
            "gene_names_stored": names_kept,
            "detected_rule": "rna_genes_detected_min and _max are the smallest and largest number of genes with at least 10 reads among the samples",
        },
        "notes": list(notes), "detail": "standard", "coverage": g_cov,
    }
    if genes["state"] == "partial":
        genes["partial_reason"] = g_reason
    return {SECTION_BIOTYPE: biotype, SECTION_GENES: genes}
