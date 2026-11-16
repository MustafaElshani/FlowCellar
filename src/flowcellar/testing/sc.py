"""Fake analysis folders of the single-cell and spatial add-on, small and with known true values.

Every layout the add-on reads, written from invented data (no real sample, slide or flow cell):

    make_spaceranger_visium   an ``outs/`` folder of a 55 um slide: summary, web summary (older page layout),
                              molecules, raw matrix, spot positions; optionally the pipestance beside it and a BAM
    make_spaceranger_hd       an ``outs/`` folder of a Visium HD slide: summary with the three bin sizes, web summary
                              (newer page layout), molecules, the 2 um layers, the 8 um matrix
    make_wf_v2                a wf-single-cell 2.x output folder: report, trace, and one sample folder with
                              ``read_tags.tsv``, ``config_stats.json``, matrices, mito table, whitelist
    make_wf_v3                a wf-single-cell 3.x output folder: report, ``params.json``, trace, and one sample
                              folder with ``<alias>.read_summary.tsv`` and matrices; Visium HD or droplet; the
                              tagged BAM carries the basecaller's ``fn`` tag of two runs

Each returns a ``FakeAnalysis`` whose ``expected`` holds what a correct reader must find. HDF5 files need
``h5py`` and BAM files ``pysam``; a maker that needs one raises ImportError without it, so tests can skip.
"""
from __future__ import annotations

import csv
import gzip
import json
import os
import random
from collections import Counter
from dataclasses import dataclass, field

import numpy as np

REPORT_NAME = "wf-single-cell-report.html"


@dataclass
class FakeAnalysis:
    path: str                     # the folder to hand to ``flowcellar analysis`` or ``detect.claims_at``
    pipeline: str
    layout: str
    expected: dict = field(default_factory=dict)
    files: dict = field(default_factory=dict)


def _w(path: str, text: str) -> str:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write(text)
    return path


def _gz(path: str, text: str) -> str:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with gzip.GzipFile(path, "wb", mtime=0) as fh:
        fh.write(text.encode("utf-8"))
    return path


def _barcode(rng: random.Random, n: int = 16) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


def gene_names(n: int = 24) -> list[str]:
    """Gene names with three mitochondrial ones (one in lower case: the prefix is matched without regard to case)."""
    return [f"GENE{i}" for i in range(n - 3)] + ["MT-ND1", "MT-CO1", "mt-atp6"]


# ---------------------------------------------------------------------------------------------------------
# Space Ranger: the web summary page and the summary table
# ---------------------------------------------------------------------------------------------------------

def _web_summary(path: str, *, layout: str, sample: dict, alarms: list[dict], diagnostics: dict, command: str | None = None) -> str:
    rows = [[k, v] for k, v in sample.items()]
    if layout == "old":
        data = {"summary": {"alarms": {"alarms": alarms}, "diagnostics": diagnostics,
                            "summary_tab": {"pipeline_info_table": {"header": ["Sample"], "rows": rows}}},
                "_resources": {"000": "data:image/png;base64,AAAA"}}
    else:
        data = {"sample": {"id": sample.get("Sample ID"), "pipeline": "Space Ranger count"},
                "tabs": {"tab_data": [{"run_summary": {"card": {"inner": {"rows": rows}, "title": {"title": "Sample Metadata"}}},
                                       "command_line": {"card": {"data": [["", [f"<code><pre>{command}</pre></code>"]]]}} if command else None}]},
                "diagnostics": diagnostics, "alarms": {"alarms": alarms}, "_resources": {}}
    html = ("<!doctype html><html><head><title>web summary</title></head><body><div id='root'></div>\n<script>\n"
            "const data = " + json.dumps(data) + ";\n</script></body></html>\n")
    return _w(path, html)


def _metrics_csv(path: str, values: dict) -> str:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(list(values))
        w.writerow([values[k] for k in values])
    return path


def _matrix_h5(path: str, names: list[str], barcodes: list[str], triples: dict, *, version: str, chemistry: str) -> None:
    """A CSC feature-barcode matrix from ``{(feature index, barcode index): count}``."""
    import h5py

    os.makedirs(os.path.dirname(path), exist_ok=True)
    by_col: dict[int, list[tuple[int, int]]] = {}
    for (f, b), v in triples.items():
        by_col.setdefault(b, []).append((f, v))
    indptr, indices, data = [0], [], []
    for b in range(len(barcodes)):
        for f, v in sorted(by_col.get(b, [])):
            indices.append(f)
            data.append(v)
        indptr.append(len(indices))
    with h5py.File(path, "w") as h:
        h.attrs["software_version"] = version
        h.attrs["chemistry_description"] = chemistry
        m = h.create_group("matrix")
        m.create_dataset("data", data=np.array(data, dtype=np.int32))
        m.create_dataset("indices", data=np.array(indices, dtype=np.int64))
        m.create_dataset("indptr", data=np.array(indptr, dtype=np.int64))
        m.create_dataset("barcodes", data=np.array([b.encode() for b in barcodes], dtype="S24"))
        m.create_dataset("shape", data=np.array([len(names), len(barcodes)], dtype=np.int32))
        ft = m.create_group("features")
        ft.create_dataset("name", data=np.array([n.encode() for n in names], dtype="S32"))
        ft.create_dataset("id", data=np.array([f"ID{i:05d}".encode() for i in range(len(names))], dtype="S32"))
        ft.create_dataset("feature_type", data=np.array([b"Gene Expression"] * len(names), dtype="S32"))
        ft.create_dataset("genome", data=np.array([b"demo"] * len(names), dtype="S32"))


def _matrix_h5_csc(path: str, names: list[str], barcodes: list[str], indptr, indices, data, *, version: str, chemistry: str) -> None:
    """The same file as ``_matrix_h5`` from the three arrays of a CSC matrix (for matrices too large for a dict)."""
    import h5py

    os.makedirs(os.path.dirname(path), exist_ok=True)
    with h5py.File(path, "w") as h:
        h.attrs["software_version"] = version
        h.attrs["chemistry_description"] = chemistry
        m = h.create_group("matrix")
        big = len(indices) > 100_000
        m.create_dataset("data", data=np.asarray(data, dtype=np.int32), compression="gzip" if big else None)
        m.create_dataset("indices", data=np.asarray(indices, dtype=np.int64), compression="gzip" if big else None)
        m.create_dataset("indptr", data=np.asarray(indptr, dtype=np.int64))
        m.create_dataset("barcodes", data=np.array([b.encode() for b in barcodes], dtype="S24"))
        m.create_dataset("shape", data=np.array([len(names), len(barcodes)], dtype=np.int32))
        ft = m.create_group("features")
        ft.create_dataset("name", data=np.array([n.encode() for n in names], dtype="S32"))
        ft.create_dataset("id", data=np.array([f"ID{i:05d}".encode() for i in range(len(names))], dtype="S32"))
        ft.create_dataset("feature_type", data=np.array([b"Gene Expression"] * len(names), dtype="S32"))
        ft.create_dataset("genome", data=np.array([b"demo"] * len(names), dtype="S32"))


def _molecule_info(path: str, counts, barcode_idx, barcodes: list[str], metrics: dict) -> None:
    import h5py

    os.makedirs(os.path.dirname(path), exist_ok=True)
    with h5py.File(path, "w") as h:
        h.create_dataset("count", data=np.asarray(counts, dtype=np.uint32))
        h.create_dataset("barcode_idx", data=np.asarray(barcode_idx, dtype=np.uint64))
        h.create_dataset("barcodes", data=np.array([b.encode() for b in barcodes], dtype="S43"))
        h.create_dataset("metrics_json", data=json.dumps(metrics))


def _bam(path: str, records: list[dict], references: list[tuple[str, int]], *, header_extra: dict | None = None) -> None:
    """A coordinate-sorted BAM with its ``.bai``. ``records`` = [{"name", "ref", "pos", "tags": {...}, "unmapped"?}]."""
    import pysam

    os.makedirs(os.path.dirname(path), exist_ok=True)
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": n, "LN": ln} for n, ln in references]}
    header.update(header_extra or {})
    names = [n for n, _ in references]
    recs = sorted([r for r in records if not r.get("unmapped")], key=lambda r: (names.index(r["ref"]), r["pos"]))
    recs += [r for r in records if r.get("unmapped")]
    save = pysam.set_verbosity(0)
    try:
        with pysam.AlignmentFile(path, "wb", header=header) as out:
            for r in recs:
                a = pysam.AlignedSegment(out.header)
                a.query_name = r["name"]
                a.query_sequence = "ACGTACGTAC"
                a.query_qualities = pysam.qualitystring_to_array("IIIIIIIIII")
                if r.get("unmapped"):
                    a.is_unmapped = True
                else:
                    a.reference_id = names.index(r["ref"])
                    a.reference_start = r["pos"]
                    a.mapping_quality = 60
                    a.cigartuples = [(0, 10)]
                    a.flag = r.get("flag", 0)
                for k, v in (r.get("tags") or {}).items():
                    a.set_tag(k, v)
                out.write(a)
        pysam.index(path)
    finally:
        pysam.set_verbosity(save)


# ---------------------------------------------------------------------------------------------------------
# Space Ranger, Visium (55 um spots)
# ---------------------------------------------------------------------------------------------------------

def make_spaceranger_visium(out_dir: str, *, seed: int = 1, sample: str = "demo_spatial_a", slide: str = "V00T00-001-A1",
                            antisense: float = 0.04, reads_under_tissue_boost: float = 6.0, alarms: list[dict] | None = None,
                            pipestance: bool = True, stale_pipestance: bool = False, bam: bool = False, folder: str = "spaceranger_out",
                            barcodes: list[str] | None = None) -> FakeAnalysis:
    """``<out_dir>/<folder>/outs`` of a slide with 8 rows of 10 spots; 20 spots lie under tissue."""
    rng = random.Random(seed)
    top = os.path.join(out_dir, folder)
    outs = os.path.join(top, "outs")
    n_rows, n_cols = 8, 10
    names = gene_names()
    spots = [(r, 2 * c + (r % 2)) for r in range(n_rows) for c in range(n_cols)]
    bcs = list(barcodes) if barcodes else [_barcode(rng) for _ in spots]
    under = [1 if (2 <= r <= 5 and 5 <= ac <= 14) else 0 for r, ac in spots]
    # molecules: (spot, feature, reads); two spots have none and are left out of the matrix
    empty = {3, 77}
    mol_spot, mol_feat, mol_count = [], [], []
    for i in range(len(spots)):
        if i in empty:
            continue
        n_mol = rng.randint(30, 60) if under[i] else rng.randint(2, 9)
        for _ in range(n_mol):
            mol_spot.append(i)
            mol_feat.append(rng.randrange(len(names)) if rng.random() > 0.15 else len(names) - 1 - rng.randrange(3))
            mol_count.append(1 + min(rng.randrange(6), rng.randrange(6)) + (int(reads_under_tissue_boost) if under[i] and rng.random() < 0.3 else 0))
    in_matrix = [i for i in range(len(spots)) if i not in empty]
    col_of = {s: k for k, s in enumerate(in_matrix)}
    triples: Counter = Counter()
    for s, f in zip(mol_spot, mol_feat):
        triples[(f, col_of[s])] += 1
    _matrix_h5(os.path.join(outs, "raw_feature_bc_matrix.h5"), names, [bcs[i] + "-1" for i in in_matrix], triples,
               version="spaceranger-3.0.1", chemistry="Visium V1 Slide")
    feature_reads = int(sum(mol_count))
    reads_under = int(sum(c for s, c in zip(mol_spot, mol_count) if under[s]))
    raw_reads = feature_reads * 3
    metrics = {"cellranger_version": "spaceranger-3.0.1", "chemistry_description": "Visium V1 Slide", "chemistry_name": "SPATIAL3Pv1",
               "slide_serial_capture_area": slide, "analysis_parameters": {"include_introns": False},
               "libraries": {"0": {"raw_read_pairs": raw_reads, "feature_read_pairs": feature_reads, "usable_read_pairs": reads_under}}}
    _molecule_info(os.path.join(outs, "molecule_info.h5"), mol_count, [col_of[s] for s in mol_spot], [bcs[i] for i in in_matrix], metrics)
    lines = ["barcode,in_tissue,array_row,array_col,pxl_row_in_fullres,pxl_col_in_fullres"]
    for i, (r, ac) in enumerate(spots):
        lines.append(f"{bcs[i]}-1,{under[i]},{r},{ac},{1000 + 100 * r},{1000 + 50 * ac}")
    _w(os.path.join(outs, "spatial", "tissue_positions.csv"), "\n".join(lines) + "\n")
    _w(os.path.join(outs, "spatial", "scalefactors_json.json"), json.dumps({"spot_diameter_fullres": 30.0, "tissue_hires_scalef": 0.05,
                                                                              "tissue_lowres_scalef": 0.015, "fiducial_diameter_fullres": 48.0}))
    umis = Counter(mol_spot)
    genes = {i: len({f for s, f in zip(mol_spot, mol_feat) if s == i}) for i in range(len(spots))}
    reads = Counter()
    for s, c in zip(mol_spot, mol_count):
        reads[s] += c
    t = [i for i in range(len(spots)) if under[i]]
    n_under = len(t)
    med_umis = float(np.median([umis.get(i, 0) for i in t]))
    med_genes = float(np.median([genes.get(i, 0) for i in t]))
    genes_under = len({f for s, f in zip(mol_spot, mol_feat) if under[s]})
    saturation = 1.0 - len(mol_count) / feature_reads
    is_mt = [n.upper().startswith("MT-") for n in names]
    mito_by_spot = Counter(s for s, f in zip(mol_spot, mol_feat) if is_mt[f])
    mito_frac = [mito_by_spot.get(i, 0) / umis[i] for i in t if umis.get(i, 0) > 0]
    values = {
        "Sample ID": sample, "Number of Spots Under Tissue": n_under, "Number of Reads": raw_reads,
        "Mean Reads per Spot": raw_reads / n_under, "Mean Reads Under Tissue per Spot": reads_under / n_under,
        "Fraction of Spots Under Tissue": n_under / len(spots), "Median Genes per Spot": med_genes,
        "Median UMI Counts per Spot": med_umis, "Valid Barcodes": 0.9471, "Valid UMIs": 0.9999,
        "Sequencing Saturation": saturation, "Q30 Bases in Barcode": 0.85, "Q30 Bases in RNA Read": 0.8, "Q30 Bases in UMI": 0.83,
        "Reads Mapped to Genome": 0.92, "Reads Mapped Confidently to Genome": 0.87,
        "Reads Mapped Confidently to Intergenic Regions": 0.12, "Reads Mapped Confidently to Intronic Regions": 0.04,
        "Reads Mapped Confidently to Exonic Regions": 0.71, "Reads Mapped Confidently to Transcriptome": 0.41,
        "Reads Mapped Antisense to Gene": antisense, "Fraction Reads in Spots Under Tissue": reads_under / feature_reads,
        "Total Genes Detected": genes_under,
    }
    _metrics_csv(os.path.join(outs, "metrics_summary.csv"), values)
    if alarms is None:
        alarms = [{"level": "WARN", "title": "High Fraction of Reads Mapped Antisense to Genes", "formatted_value": f"{antisense * 100:.1f}%",
                   "raw_value": antisense, "message": "Ideal < 10%. This can indicate use of an unsupported chemistry type.", "id": "x"}] \
            if antisense > 0.1 else []
    sample_table = {"Sample ID": sample, "Sample Description": "", "Chemistry": "Visium V1 Slide - 3'", "Slide Serial Number": slide,
                    "Transcriptome": "demo-2024-A", "Pipeline Version": "spaceranger-3.0.1", "Image Reorientation": "On",
                    "Loupe Manual Alignment": "None"}
    diag = {"tso_frac": 0.0003, "low_support_umi_reads_frac": 0.003}
    _web_summary(os.path.join(outs, "web_summary.html"), layout="old", sample=sample_table, alarms=alarms, diagnostics=diag)
    files = {"outs": outs}
    if pipestance:
        _w(os.path.join(top, "_versions"), json.dumps({"martian": "v4.0.12", "pipelines": "spaceranger-3.0.1"}))
        _w(os.path.join(top, "_cmdline"), f"spaceranger count --id {sample} --transcriptome ../../refs/refdata-demo-2024-A --create-bam true "
                                          f"--slide {slide[:-3]} --area {slide[-2:]} --image ../images/demo.tif --output-dir ./{folder} --localcores 8")
        pj = {"total_read_pairs": raw_reads + (7 if stale_pipestance else 0), "filtered_bcs": n_under, "corrected_bc_frac": 0.0466,
              "tso_frac": 0.00031, "low_support_umi_reads_frac": 0.00296, "tissue_registration_type": "automatic"}
        _w(os.path.join(top, "SPATIAL_RNA_COUNTER_CS", "SPATIAL_RNA_COUNTER", "SPATIAL_REPORTER", "fork0", "join-u0000000000", "files",
                        "metrics_summary_json.json"), json.dumps(pj))
    exp_bam = {}
    if bam:
        refs = [("chr1", 100000), ("chrM", 16000)]
        recs = [{"name": f"r{i}_0", "ref": "chr1", "pos": 10 + 13 * i} for i in range(300)]
        recs += [{"name": f"m{i}_0", "ref": "chrM", "pos": 5 + 7 * i} for i in range(100)]
        recs += [{"name": f"u{i}_0", "unmapped": True} for i in range(50)]
        _bam(os.path.join(outs, "possorted_genome_bam.bam"), recs, refs,
             header_extra={"CO": [f"slide_serial_capture_area:{slide}"], "RG": [{"ID": f"{sample}:0:1:unknow_flowcell:0", "SM": sample, "PL": "ILLUMINA"}]})
        exp_bam = {"reads_bam_mapped": 400, "reads_bam_unmapped": 50, "reads_chrm_fraction": 100 / 400}
    expected = {
        "label": sample, "slide": slide, "assay": "visium", "pipeline_version": "3.0.1", "reference": "demo-2024-A", "spots": len(spots),
        "summary": {
            "spatial_units_under_tissue": n_under, "sr_reads": raw_reads, "spatial_genes_per_unit": med_genes, "spatial_umis_per_unit": med_umis,
            "spatial_fraction_under_tissue": n_under / len(spots), "spatial_fraction_reads_under_tissue": reads_under / feature_reads,
            "spatial_genes_detected": genes_under, "sr_saturation": saturation, "sr_antisense": antisense, "sr_valid_barcodes": 0.9471,
            "sr_valid_umis": 0.9999, "sr_mapped_transcriptome": 0.41, "sr_reads_feature": feature_reads, "sr_reads_usable": reads_under,
            "sr_usable_reads_fraction": reads_under / raw_reads, "barcodes_unique_molecules": len(mol_count),
            "saturation_molecule_reads": feature_reads, "barcodes_umis_median": med_umis, "barcodes_genes_median": med_genes,
            "barcodes_reads_median": float(np.median([reads.get(i, 0) for i in t])),
            "barcodes_detected": len(reads), "mito_fraction_median": float(np.median(mito_frac)),
            "mito_umi_fraction": sum(mito_by_spot.values()) / len(mol_spot),
            "spatial_fraction_umis_under_tissue": sum(umis[i] for i in t) / len(mol_spot),
            "sr_corrected_bc_frac": 0.0466 if (pipestance and not stale_pipestance) else None, **exp_bam},
        "count_hist": np.bincount(mol_count).tolist(), "alarms": [a["title"] for a in alarms],
        "umis_by_spot": [umis.get(i, 0) for i in range(len(spots))], "genes_by_spot": [genes.get(i, 0) for i in range(len(spots))],
        "under": under, "spots_rc": spots, "barcodes": bcs, "stale_pipestance": stale_pipestance,
    }
    return FakeAnalysis(outs, "spaceranger", "spaceranger", expected, files)


# ---------------------------------------------------------------------------------------------------------
# Space Ranger, Visium HD
# ---------------------------------------------------------------------------------------------------------

def make_spaceranger_hd(out_dir: str, *, seed: int = 2, sample: str = "demo_hd_day2", slide: str = "H1-TESTSLD-A1", n: int = 100,
                        registration_alarm: bool = True, folder: str = "spaceranger", image: bool = True, tissue=None, density=None,
                        grey=None, n_genes: int = 30, gene_share: tuple[float, float] | None = None, summary_values: dict | None = None,
                        reads_per_molecule: float | None = None) -> FakeAnalysis:
    """``<out_dir>/<folder>/outs`` of a slide of ``n`` x ``n`` squares of 2 um (25 x 25 bins of 8 um for n = 100).
    Tissue is a disc; molecules, reads and mitochondrial molecules are written square by square.

    For a slide that looks like one: ``tissue`` (an ``n`` x ``n`` mask) and ``density`` (expected molecules per square
    under tissue, a number or an ``n`` x ``n`` array) replace the disc and its even cover, ``grey`` the picture,
    ``n_genes`` the number of genes, ``gene_share`` (low, high) the genes of an 8 um bin as a share of its molecules,
    ``reads_per_molecule`` the mean of the reads a molecule was seen in (drawn per molecule, in place of a fixed
    pattern), and ``summary_values`` values of the summary table that this function otherwise writes as fixed numbers."""
    import h5py

    rng = np.random.default_rng(seed)
    outs = os.path.join(out_dir, folder, "outs")
    os.makedirs(outs, exist_ok=True)
    names = gene_names(n_genes)
    rr, cc = np.mgrid[0:n, 0:n]
    if tissue is None:
        tissue = ((rr - n * 0.45) ** 2 + (cc - n * 0.5) ** 2) < (n * 0.33) ** 2
    if density is None:
        umis = np.where(tissue, rng.poisson(3.0, (n, n)), (rng.random((n, n)) < 0.05).astype(np.int64))
    else:
        umis = np.where(tissue, rng.poisson(density, (n, n)), (rng.random((n, n)) < 0.05 * min(1.0, float(np.mean(density)))).astype(np.int64))
    reads = umis * 4 + np.where(umis > 0, rng.integers(0, 6, (n, n)), 0) + (rng.random((n, n)) < 0.02)
    mito = np.minimum(umis, rng.binomial(1, 0.08, (n, n)) * rng.integers(1, 3, (n, n)))
    if grey is None:
        grey = (40 + 2 * rr + np.where(tissue, 60, 0)).clip(0, 255).astype(np.uint8)
    f8 = 4
    n8 = -(-n // f8)

    def to8(a):
        flat = (rr // f8) * n8 + (cc // f8)
        return np.bincount(flat.ravel(), weights=a.ravel().astype(np.float64), minlength=n8 * n8).reshape(n8, n8)
    tissue8 = to8(tissue) > 0
    umis8, reads8, mito8 = to8(umis), to8(reads), to8(mito)
    r8, c8 = np.nonzero(tissue8)
    if gene_share is None:
        genes8 = np.minimum(umis8[r8, c8], rng.integers(1, len(names), len(r8))).astype(np.int64)
    else:
        genes8 = np.clip(np.rint(umis8[r8, c8] * rng.uniform(gene_share[0], gene_share[1], len(r8))), 0, len(names) - 3).astype(np.int64)

    def coo(group, mask, values, dtype):
        r, c = np.nonzero(mask)
        group.create_dataset("row", data=r.astype(np.uint32), compression="gzip")
        group.create_dataset("col", data=c.astype(np.uint32), compression="gzip")
        group.create_dataset("data", data=values[r, c].astype(dtype), compression="gzip")
    with h5py.File(os.path.join(outs, "feature_slice.h5"), "w") as h:
        h.attrs["filetype"] = "feature_slice"
        h.attrs["software_version"] = "4.1.0"
        h.attrs["metadata_json"] = json.dumps({"sample_id": sample, "slide_name": "visium_hd_v1", "nrows": n, "ncols": n, "spot_pitch": 2.0})
        coo(h.create_group("umis/total"), umis > 0, umis, np.uint32)
        coo(h.create_group("reads/sequenced"), reads > 0, reads, np.uint32)
        coo(h.create_group("masks/filtered"), tissue, tissue.astype(np.uint8), np.uint8)
        g = h.create_group("masks/square_008um")
        g.create_dataset("row", data=r8.astype(np.int32))
        g.create_dataset("col", data=c8.astype(np.int32))
        g.create_dataset("data", data=np.ones(len(r8), dtype=np.uint8))
        if image:
            coo(h.create_group("images/microscope"), np.ones((n, n), bool), grey, np.uint8)
        ft = h.create_group("features")
        ft.create_dataset("name", data=np.array([x.encode() for x in names], dtype="S256"))
        sl = h.create_group("feature_slices")
        # the mitochondrial molecules of a square are split over the three mitochondrial features
        parts = [mito // 2, mito - mito // 2, np.zeros_like(mito)]
        for k, part in zip(range(len(names) - 3, len(names)), parts):
            if part.any():
                coo(sl.create_group(str(k)), part > 0, part, np.uint32)
        coo(sl.create_group("0"), umis > 2, umis, np.uint32)
    total_umis = int(umis.sum())
    # molecules: one count per molecule, their reads add up to a little less than the reads on the slide
    counts = np.repeat(np.array([1, 2, 3, 4, 7, 70, 130]), [total_umis // 3, total_umis // 4, total_umis // 6, total_umis // 8, 5, 3, 2])
    counts = np.concatenate([counts, np.ones(total_umis - len(counts), dtype=np.int64)]) if len(counts) < total_umis else counts[:total_umis]
    if reads_per_molecule:
        counts = rng.geometric(1.0 / reads_per_molecule, total_umis)
    feature_reads = int(counts.sum())
    raw_reads = int(reads.sum()) + 9000
    usable = int(feature_reads * 0.95)
    metrics = {"cellranger_version": "4.1.0", "chemistry_description": "Visium HD v1 - 3'", "chemistry_name": "SPATIAL-HD-v1-3P",
               "slide_serial_capture_area": slide, "analysis_parameters": {"include_introns": True},
               "libraries": {"0": {"raw_read_pairs": raw_reads, "feature_read_pairs": feature_reads, "usable_read_pairs": usable}}}
    _molecule_info(os.path.join(outs, "molecule_info.h5"), counts, np.zeros(len(counts), dtype=np.uint64), ["s_002um_00000_00000"], metrics)
    bcs8 = [f"s_008um_{r:05d}_{c:05d}-1" for r, c in zip(r8, c8)]
    # every bin holds one molecule of each of its first genes: the matrix as three arrays, column by column
    indptr = np.concatenate([[0], np.cumsum(genes8)])
    indices = np.arange(int(indptr[-1])) - np.repeat(indptr[:-1], genes8)
    if gene_share is not None:                               # not always the first genes: a stretch that starts anywhere
        first = (rng.random(len(genes8)) * (len(names) - 3 - genes8)).astype(np.int64)
        indices = indices + np.repeat(first, genes8)
    _matrix_h5_csc(os.path.join(outs, "binned_outputs", "square_008um", "filtered_feature_bc_matrix.h5"), names, bcs8, indptr, indices,
                   np.ones(len(indices), dtype=np.int32), version="4.1.0", chemistry="Visium HD v1 - 3'")
    _w(os.path.join(outs, "binned_outputs", "square_008um", "spatial", "scalefactors_json.json"),
       json.dumps({"bin_size_um": 8.0, "microns_per_pixel": 0.22, "spot_diameter_fullres": 36.3}))
    squares = int(tissue.sum())
    bins8 = int(tissue8.sum())
    umis_under = int(umis[tissue].sum())
    saturation = 1.0 - total_umis / feature_reads
    values = {
        "Sample ID": sample, "Number of Reads": raw_reads, "Valid Barcodes": 0.8791, "Valid UMI Sequences": 0.99997,
        "Sequencing Saturation": saturation, "Q30 Bases in Barcode": 0.8375, "Q30 Bases in RNA Read": 0.8308, "Q30 Bases in UMI": 0.8229,
        "Reads Mapped to Genome": 0.9562, "Reads Mapped Confidently to Genome": 0.9364,
        "Reads Mapped Confidently to Intergenic Regions": 0.0283, "Reads Mapped Confidently to Intronic Regions": 0.089,
        "Reads Mapped Confidently to Exonic Regions": 0.819, "Reads Mapped Confidently to Transcriptome": 0.8875,
        "Reads Mapped Antisense to Gene": 0.0134, "Fraction Reads in Squares Under Tissue": 0.9456, "Total Genes Detected": 27,
        "Number of Squares Under Tissue 2 µm": squares, "Mean Reads Under Tissue per Square 2 µm": raw_reads / squares,
        "Fraction of Squares Under Tissue 2 µm": squares / (n * n), "Mean Genes Under Tissue per Square 2 µm": 2.5,
        "Mean UMIs Under Tissue per Square 2 µm": umis_under / squares, "Total Genes Detected Under Tissue 2 µm": 27,
        "Number of Bins Under Tissue 8 µm": bins8, "Mean Reads Under Tissue per Bin 8 µm": raw_reads / bins8,
        "Fraction of Bins Under Tissue 8 µm": bins8 / (n8 * n8), "Mean Genes Under Tissue per Bin 8 µm": float(genes8.mean()),
        "Mean UMIs Under Tissue per Bin 8 µm": float(umis8[tissue8].mean()), "Total Genes Detected Under Tissue 8 µm": 27,
        "UMIs per sq mm of Tissue": 786906.6, "Reads per sq mm of Tissue": 4432198.3,
        "Number of Bins Under Tissue 16 µm": bins8 // 4 + 3, "Mean Reads Under Tissue per Bin 16 µm": raw_reads / (bins8 // 4 + 3),
        "Fraction of Bins Under Tissue 16 µm": 0.77, "Mean Genes Under Tissue per Bin 16 µm": 12.5,
        "Mean UMIs Under Tissue per Bin 16 µm": 60.2, "Total Genes Detected Under Tissue 16 µm": 27,
        "Number of Cells": 412, "Reads in Cells": 0.6, "UMIs in Cells": 0.6429, "Mean Reads per Cell": 1271.2, "Median Genes per Cell": 15.0,
        "Median UMIs per Cell": 18.0, "Median Cell Area (μm²)": 164.0, "Median Nucleus Area (μm²)": 164.0,
        "Microscope Image Pixel Size": 0.22,
    }
    values.update(summary_values or {})
    _metrics_csv(os.path.join(outs, "metrics_summary.csv"), values)
    alarms = [{"level": "WARN", "title": "UMI Registration Warning", "formatted_value": "True",
               "message": "The automatic UMI-based registration could not determine a confident offset."}] if registration_alarm else []
    sample_table = {"Sample ID": sample, "Sample Description": "", "Chemistry": "Visium HD H1 Slide - 3'", "Slide Serial Number": slide,
                    "Transcriptome": "demo-2024-A", "Pipeline Version": "4.1.0", "Image Reorientation": "On",
                    "Loupe Manual Alignment": "demo-fiducials-image-registration.json", "CytAssist Run Name": "2026-01-01_10-00-00"}
    diag = {"corrected_bc_frac": 0.493, "tso_frac": 0.00008, "low_support_umi_reads_frac": 0.0253,
            "hd_layout_offset": {"x_offset": 0.0, "y_offset": 0.0, "gaussian_fit_x_offset": 14.96, "gaussian_fit_y_offset": 16.08,
                                 "has_gaussian_fit_failed": False, "is_gaussian_fit_outside_range": bool(registration_alarm)}}
    _web_summary(os.path.join(outs, "web_summary.html"), layout="new", sample=sample_table, alarms=alarms, diagnostics=diag,
                 command=f"spaceranger count --id={sample} --transcriptome=/refs/refdata-demo-2024-A --slide={slide[:-3]} --area={slide[-2:]} --create-bam=true")
    has8 = tissue8 & (umis8 > 0)
    expected = {
        "label": sample, "slide": slide, "assay": "visium_hd", "pipeline_version": "4.1.0", "reference": "demo-2024-A",
        "summary": {
            "spatial_units_under_tissue": bins8, "spatial_squares_under_tissue_2um": squares, "sr_reads": raw_reads,
            "spatial_umis_per_unit": float(umis8[tissue8].mean()), "spatial_genes_per_unit": float(genes8.mean()),
            "sr_saturation": saturation, "sr_valid_barcodes": 0.8791, "sr_valid_umis": 0.99997, "sr_antisense": 0.0134,
            "spatial_umis_in_cells": 0.6429, "spatial_cells": 412, "sr_reads_feature": feature_reads, "sr_reads_usable": usable,
            "sr_usable_reads_fraction": usable / raw_reads, "barcodes_unique_molecules": total_umis,
            "barcodes_umis_mean": float(umis8[tissue8].mean()), "barcodes_umis_median": float(np.median(umis8[tissue8])),
            "barcodes_genes_mean": float(genes8.mean()), "barcodes_genes_median": float(np.median(genes8)),
            "barcodes_reads_median": float(np.median(reads8[tissue8])), "barcodes_detected": int((reads8 > 0).sum()),
            "mito_umi_fraction": float(mito.sum() / total_umis), "mito_fraction_median": float(np.median((mito8 / np.maximum(umis8, 1))[has8])),
            "spatial_fraction_umis_under_tissue": umis_under / total_umis, "sr_corrected_bc_frac": 0.493,
            "spatial_registration_flag": bool(registration_alarm)},
        "scalars": {"total_umis": total_umis, "spatial_umis_under_tissue": umis_under, "mito_umis": int(mito.sum())},
        "count_hist": np.bincount(counts).tolist(), "alarms": [a["title"] for a in alarms], "n": n, "n8": n8,
        "umis": umis, "reads": reads, "tissue": tissue, "mito": mito, "grey": grey, "umis8": umis8, "tissue8": tissue8,
        "genes8": (r8, c8, genes8), "bin_rows": 3,
    }
    return FakeAnalysis(outs, "spaceranger", "spaceranger", expected, {"outs": outs})


# ---------------------------------------------------------------------------------------------------------
# wf-single-cell
# ---------------------------------------------------------------------------------------------------------
HEADER = ("read_id\tuncorrected_barcode\tcorrected_barcode\tquality_barcode\tuncorrected_umi\tcorrected_umi\tquality_umi\t"
          "gene\ttranscript\tstart\tend\tchr\n")


def _read_id(rng: random.Random, zeros: int = 0) -> str:
    hexes = "".join(rng.choice("0123456789abcdef") for _ in range(32))
    hexes = "0" * zeros + hexes[zeros:]
    if zeros < 4 and hexes.startswith("0000"):
        hexes = "1" + hexes[1:]
    return f"{hexes[:8]}-{hexes[8:12]}-{hexes[12:16]}-{hexes[16:20]}-{hexes[20:]}"


def _reads(rng: random.Random, barcodes: list[str], names: list[str], n: int, umi_len: int, weights: list[float] | None = None) -> list[tuple]:
    """(read id, barcode, umi, gene, transcript, chromosome). A tenth has no gene, a few no UMI; one in forty read
    ids starts with zeros (the sketch)."""
    rows = []
    for i in range(n):
        bc = rng.choices(barcodes, weights=weights)[0] if weights else rng.choice(barcodes)
        umi = "".join(rng.choice("ACGT") for _ in range(3)) + "A" * (umi_len - 3)
        if rng.random() < 0.03:
            umi = "-"
        gene = rng.choice(names) if rng.random() > 0.1 else "-"
        tx = f"TX_{gene}_1" if gene != "-" and rng.random() < 0.6 else "-"
        rid = _read_id(rng, zeros=6 if i % 40 == 0 else 0) + "_0"
        rows.append((rid, bc, umi, gene, tx, f"chr{1 + (names.index(gene) % 3) if gene != '-' else 1}"))
    rows.sort(key=lambda r: r[5])                           # the workflow writes the table chromosome by chromosome
    return rows


def _write_table(path: str, rows: list[tuple]) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write(HEADER)
        for rid, bc, umi, gene, tx, chrom in rows:
            fh.write(f"{rid}\t{bc.split('-')[0][:16] if not bc.startswith('s_') else 'ACGT'}\t{bc}\tIIII\t{umi}\t{umi}\tIIII\t{gene}\t{tx}\t100\t900\t{chrom}\n")


def _mtx(folder: str, names: list[str], barcodes: list[str], triples: dict) -> None:
    lines = ["%%MatrixMarket matrix coordinate integer general",
             '%metadata_json: {"software_version": "ont-single-cell","format_version": 2}',
             f"{len(names)} {len(barcodes)} {len(triples)}"]
    lines += [f"{f + 1} {b + 1} {v}" for (f, b), v in sorted(triples.items())]
    _gz(os.path.join(folder, "matrix.mtx.gz"), "\n".join(lines) + "\n")
    _gz(os.path.join(folder, "barcodes.tsv.gz"), "".join(b + "\n" for b in barcodes))
    _gz(os.path.join(folder, "features.tsv.gz"), "".join(f"unknown_{i:05d}\t{n}\tGene Expression\n" for i, n in enumerate(names)))


def _trace(path: str, failed: int = 0) -> None:
    rows = ["task_id\thash\tnative_id\tname\tstatus\texit\tsubmit\tduration\trealtime\t%cpu\tpeak_rss\tpeak_vmem\trchar\twchar"]
    procs = [("pipeline:preprocess:call_paftools", "1m 59s", "1m 56s", "52.1 MB"), ("pipeline:process_bams:assign_features (3)", "30m", "29m 30s", "5.5 GB"),
             ("pipeline:process_bams:assign_features (4)", "31m", "30m 30s", "4.9 GB"), ("pipeline:process_bams:tag_bam", "1h 2m", "1h", "22.2 GB"),
             ("pipeline:makeReport", "40s", "35s", "1.1 GB")]
    for i, (name, dur, rt, mem) in enumerate(procs):
        bad = i < failed
        rows.append(f"{i + 1}\tab/{i:06d}\t{1000 + i}\t{name}\t{'FAILED' if bad else 'COMPLETED'}\t{1 if bad else 0}\t2026-01-0{1 + i // 3} 10:{10 + i}:00.000"
                    f"\t{dur}\t{rt}\t100.0%\t{mem}\t{mem}\t1 GB\t1 GB")
    _w(path, "\n".join(rows) + "\n")


def _report_html(path: str, *, version: str, v3: bool, alias: str, numbers: dict, params: dict, adapter_pct: dict | None = None,
                 saturation: tuple[int, float] | None = None) -> None:
    def table(heading: str, rows: list[list]) -> str:
        body = "".join("<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>" for r in rows)
        return f"<section><h2>{heading}</h2><table><tbody>{body}</tbody></table></section>\n"
    html = "<!doctype html><html><head><style>td{padding:2px}</style><script>var x = '<table><tr><td>not data</td></tr></table>';</script></head><body><main>\n"
    html += "<h1>Workflow Single Cell report</h1>\n"
    if v3:
        html += table("Experiment summary", [["Input reads", f"{numbers['input']:,}"], ["8 µm bins with data", f"{numbers['units']:,}"],
                                             ["Reads per 8 µm bin (mean)", "263"]])
        html += table("Alignment / feature summary", [["Pass reads", f"{numbers['input']:,}"], ["Mapped", f"{numbers['mapped']:,}"],
                                                      ["Unmapped", f"{numbers['unmapped']:,}"], ["Supplementary", "12"]])
        html += table("Read assignment summary", [["", "Full length", "Valid barcode", "Gene assigned", "Transcript assigned"],
                                                  ["count", f"{numbers['full_length']:,}", f"{numbers['valid_barcode']:,}",
                                                   f"{numbers['gene_assigned']:,}", f"{numbers['transcript_assigned']:,}"],
                                                  ["% full length reads", "100.00%", "88.69%", "81.80%", "43.57%"]])
    else:
        html += table("Single cell sample summary", [["sample ID", "reads", "cells", "genes", "transcripts"],
                                                     [alias, numbers["input"], numbers["units"], 20, 30]])
        html += table("Alignment summary", [["sample", "reads_aligned", "primary", "secondary", "supplementary", "unmapped"],
                                            [alias, numbers["mapped"] + numbers["unmapped"], numbers["mapped"], 5, 3, numbers["unmapped"]]])
    html += table("Workflow parameters", [["Key", "Value"]] + [[k, v] for k, v in params.items()])
    if v3 and (adapter_pct or saturation):
        def cds(cols: dict) -> dict:
            return {"type": "object", "name": "ColumnDataSource", "attributes": {"data": {"type": "map", "entries": [[k, v] for k, v in cols.items()]}}}
        figs = []
        if adapter_pct:
            figs.append({"type": "object", "name": "Figure", "attributes": {
                "below": [{"type": "object", "name": "CategoricalAxis", "attributes": {"axis_label": "Primer configuration"}}],
                "renderers": [cds({"groups": list(adapter_pct), alias: list(adapter_pct.values())})]}})
        if saturation:
            figs.append({"type": "object", "name": "Figure", "attributes": {
                "title": {"type": "object", "name": "Title", "attributes": {"text": "Sequencing saturation"}},
                "renderers": [cds({"x": [saturation[0], saturation[0] // 2], "y": [saturation[1], saturation[1] * 0.9]})]}})
        doc = json.dumps({"doc-0000": {"roots": figs}}).replace("\\", "\\\\").replace("'", "\\'")
        html += "<script type=\"text/javascript\">(function() { const docs_json = '" + doc + "';\n })();</script>\n"
    html += (f"</main><footer><div><h4>About this report</h4><p>This report was produced using the <code>epi2me-labs/wf-single-cell</code> "
             f"nextflow workflow ({version}). </p></div></footer></body></html>\n")
    _w(path, html)


def _wf_truth(rows: list[tuple], names: list[str], bin_of=lambda bc: bc) -> dict:
    """The true numbers of a per-read table, counted the slow and obvious way."""
    mol = Counter((bc, gene, umi) for _, bc, umi, gene, _, _ in rows if umi != "-" and gene != "-")
    per_unit = Counter(bin_of(bc) for _, bc, _, _, _, _ in rows)
    triples: Counter = Counter()
    for (bc, gene, _umi), _n in mol.items():
        triples[(names.index(gene), bin_of(bc))] += 1
    return {"rows": len(rows), "valid_umi": sum(1 for r in rows if r[2] != "-"), "gene_assigned": sum(1 for r in rows if r[3] != "-"),
            "transcript_assigned": sum(1 for r in rows if r[4] != "-"), "molecules": len(mol), "molecule_reads": sum(mol.values()),
            "molecule_hist": np.bincount(list(mol.values())).tolist(), "per_unit": per_unit, "triples": triples}


def _wf_summary(truth: dict, config: dict, units_with_umis: dict, genes_per_unit: dict) -> dict:
    full = config["full_length"]
    u = np.array(sorted(units_with_umis.values()), dtype=np.float64)
    g = np.array(sorted(genes_per_unit.values()), dtype=np.float64)
    per = np.array(sorted(truth["per_unit"].values()), dtype=np.float64)
    return {"reads_input": config["reads"], "reads_full_length": full, "reads_full_length_fraction": min(1.0, full / config["reads"]),
            "reads_valid_barcode": truth["rows"], "reads_valid_barcode_fraction": truth["rows"] / full,
            "reads_valid_umi_fraction": truth["valid_umi"] / truth["rows"], "reads_gene_assigned": truth["gene_assigned"],
            "reads_gene_assigned_fraction": truth["gene_assigned"] / full, "reads_transcript_assigned": truth["transcript_assigned"],
            "reads_transcript_assigned_fraction": truth["transcript_assigned"] / full, "barcodes_detected": len(truth["per_unit"]),
            "barcodes_unique_molecules": truth["molecules"], "saturation_value": 1.0 - truth["molecules"] / truth["gene_assigned"],
            "barcodes_umis_median": float(np.median(u)), "barcodes_genes_median": float(np.median(g)),
            "barcodes_reads_median": float(np.median(per)), "barcodes_reads_mean": float(per.mean())}


def make_wf_v2(out_dir: str, *, seed: int = 3, alias: str = "demo_spatial_a", folder: str = "wf-out", n_reads: int = 3000,
               barcodes: list[str] | None = None, failed_tasks: int = 0, bam: bool = False, plus: float = 0.54) -> FakeAnalysis:
    """A workflow 2.x output folder for a 55 um slide: ``<out_dir>/<folder>/{report, execution/, <alias>/}``."""
    rng = random.Random(seed)
    top = os.path.join(out_dir, folder)
    sdir = os.path.join(top, alias)
    names = gene_names()
    bcs = list(barcodes) if barcodes else [_barcode(rng) for _ in range(40)]
    rows = _reads(rng, bcs, names, n_reads, 12)
    _write_table(os.path.join(sdir, "read_tags.tsv"), rows)
    truth = _wf_truth(rows, names)
    order = sorted({b for (_, b) in truth["triples"]})
    col = {b: i for i, b in enumerate(order)}
    triples = {(f, col[b]): v for (f, b), v in truth["triples"].items()}
    _mtx(os.path.join(sdir, "gene_raw_feature_bc_matrix"), names, [b + "-1" for b in order], triples)
    _mtx(os.path.join(sdir, "transcript_raw_feature_bc_matrix"), [f"TX{i}" for i in range(37)], [b + "-1" for b in order], {(0, 0): 1})
    _mtx(os.path.join(sdir, "gene_processed_feature_bc_matrix"), names[:5], [b + "-1" for b in order[:5]], {(0, 0): 1})
    umis = Counter()
    genes = Counter()
    mito = Counter()
    for (f, b), v in truth["triples"].items():
        umis[b] += v
        genes[b] += 1
        if names[f].upper().startswith("MT-"):
            mito[b] += v
    _w(os.path.join(sdir, "gene.expression.mito-per-cell.tsv"), "CB\tmito_pct\n" + "".join(f"{b}\t{100.0 * mito[b] / umis[b]}\n" for b in order))
    _w(os.path.join(sdir, "whitelist.tsv"), "".join(b + "\n" for b in order))
    reads_in = int(n_reads / 0.71)
    full = int(reads_in * 0.917)
    stranded = int(reads_in * 0.98)
    config = {"reads": reads_in, "full_length": full, "stranded": stranded, "plus": int(stranded * plus), "minus": stranded - int(stranded * plus),
              "full_len": full, "single_adapter2": int(reads_in * 0.035), "single_adapter1": int(reads_in * 0.028), "other": int(reads_in * 0.014),
              "no_adapters": int(reads_in * 0.006), "double_adapter1": 7, "double_adapter2": 3}
    _w(os.path.join(sdir, "config_stats.json"), json.dumps(config, indent=4))
    _trace(os.path.join(top, "execution", "trace.txt"), failed=failed_tasks)
    params = {"kit": "visium:v1", "full_length_only": "True", "matrix_min_genes": "200", "matrix_max_mito": "20", "mito_prefix": "MT-",
              "expected_cells": "500", "ref_genome_dir": "../refs/refdata-gex-demo-2024-A", "out_dir": f"./{folder}"}
    _report_html(os.path.join(top, REPORT_NAME), version="2.2.0", v3=False, alias=alias,
                 numbers={"input": reads_in, "units": len(order), "mapped": int(reads_in * 0.85), "unmapped": reads_in - int(reads_in * 0.85)}, params=params)
    _w(os.path.join(top, ".ipynb_checkpoints", REPORT_NAME), "<html>a checkpoint copy</html>")
    exp_bam = {}
    if bam:
        refs = [("chr1", 100000), ("chrM", 16000)]
        recs = [{"name": f"r{i}_0", "ref": "chr1", "pos": 10 + 13 * i} for i in range(240)]
        recs += [{"name": f"m{i}_0", "ref": "chrM", "pos": 5 + 7 * i} for i in range(60)]
        _bam(os.path.join(sdir, "tagged.bam"), recs, refs)
        exp_bam = {"reads_chrm_fraction": 60 / 300}
    mfrac = np.array([mito[b] / umis[b] for b in order])
    summary = _wf_summary(truth, config, dict(umis), dict(genes))
    summary.update({"barcodes_cells": len(order), "barcodes_genes_detected": len(names), "barcodes_transcripts_detected": 37,
                    "reads_plus_fraction": config["plus"] / (config["plus"] + config["minus"]), "mito_fraction_median": float(np.median(mfrac)),
                    "mito_umi_fraction": sum(mito.values()) / truth["molecules"], "workflow_tasks": 5, "workflow_failed_tasks": failed_tasks,
                    "reads_report_unmapped": reads_in - int(reads_in * 0.85), **exp_bam})
    expected = {"label": alias, "assay": "visium", "pipeline_version": "2.2.0", "reference": "demo-2024-A", "summary": summary,
                "molecule_hist": truth["molecule_hist"], "adapter": {k: config[k] for k in ("full_length", "single_adapter1", "single_adapter2",
                                                                                              "double_adapter1", "double_adapter2", "no_adapters", "other")},
                "sketch_ids": sorted({r[0].replace("-", "")[:16] for r in rows if r[0].startswith("0000")}), "barcodes": bcs,
                "per_unit": dict(truth["per_unit"])}
    return FakeAnalysis(top, "wf-single-cell", "wf_v2", expected, {"sample": sdir, "top": top})


def make_wf_v3(out_dir: str, *, seed: int = 4, alias: str = "SAMPLE_S1_L001", folder: str = "wf-single-cell", kit: str = "visium_hd:v1",
               n_reads: int = 4000, runs=(("TST00001", "0a1b2c3d", "11112222", 3), ("TST00002", "4e5f6a7b", "33334444", 2)),
               spaceranger_bam: str | None = "./demo_analysis/spaceranger/outs/possorted_genome_bam.bam", out_rel: str | None = None,
               bam: bool = True, config_stats: bool = False, expected_cells: int | None = None, grid: int = 100) -> FakeAnalysis:
    """A workflow 3.x output folder: ``<out_dir>/<folder>/{report, params.json, execution/, <alias>/}``.

    ``kit="visium_hd:v1"`` writes 2 um barcodes and ``_2um`` / ``_8um`` matrices; any other kit writes droplet
    barcodes and plain matrix names. ``runs`` gives the flow cell, protocol run, acquisition and weight of the
    runs whose ``fn`` tag the BAM records carry (3 : 2 = 60 % and 40 %)."""
    rng = random.Random(seed)
    top = os.path.join(out_dir, folder)
    sdir = os.path.join(top, alias)
    names = gene_names(30)
    hd = kit.startswith("visium_hd")
    if hd:
        squares = [(r, c) for r in range(10, grid - 10, 3) for c in range(8, grid - 8, 3)]
        bcs = [f"s_002um_{r:05d}_{c:05d}-1" for r, c in squares]
        bin_of = lambda bc: f"s_008um_{int(bc[8:13]) // 4:05d}_{int(bc[14:19]) // 4:05d}"   # noqa: E731
        rows = _reads(rng, bcs, names, n_reads, 9)
    else:
        cells = [_barcode(rng) + "-1" for _ in range(60)]
        weights = [40.0] * 12 + [1.0] * 48                  # twelve real cells, the rest background
        bcs = cells
        bin_of = lambda bc: bc                               # noqa: E731
        rows = _reads(rng, bcs, names, n_reads, 12, weights)
    _write_table(os.path.join(sdir, f"{alias}.read_summary.tsv"), rows)
    truth = _wf_truth(rows, names, bin_of)
    order = sorted({b for (_, b) in truth["triples"]})
    col = {b: i for i, b in enumerate(order)}
    triples = {(f, col[b]): v for (f, b), v in truth["triples"].items()}
    suffix8, suffix2 = ("_8um", "_2um") if hd else ("", None)
    _mtx(os.path.join(sdir, f"{alias}.gene_raw_feature_bc_matrix{suffix8}"), names, [b if b.endswith("-1") else b + "-1" for b in order], triples)
    _mtx(os.path.join(sdir, f"{alias}.transcript_raw_feature_bc_matrix{suffix8}"), [f"TX{i}" for i in range(52)], [order[0] + "-1"], {(0, 0): 1})
    if suffix2:
        _mtx(os.path.join(sdir, f"{alias}.gene_raw_feature_bc_matrix{suffix2}"), names, bcs[:3], {(0, 0): 1})
    umis, genes, mito = Counter(), Counter(), Counter()
    for (f, b), v in truth["triples"].items():
        umis[b] += v
        genes[b] += 1
        if names[f].upper().startswith("MT-"):
            mito[b] += v
    _w(os.path.join(sdir, f"{alias}.gene_expression_mito_per_cell.tsv"), "CB\tmito_pct\n" + "".join(f"{b}\t{100.0 * mito[b] / umis[b]}\n" for b in order))
    reads_in = int(n_reads / 0.88)
    full = int(reads_in * 0.987)
    config = {"reads": reads_in, "full_length": full}
    out_rel = out_rel or f"./demo_analysis/{folder}"
    params = {"bam": "./demo_analysis/prep/SAMPLE_S1_L001.bam", "spaceranger_bam": spaceranger_bam, "adapter_stats": "./demo_analysis/prep/configs.json",
              "out_dir": out_rel, "kit": kit, "threads": 8, "full_length_only": True, "min_read_qual": 10, "barcode_max_ed": 2,
              "barcode_min_ed_diff": 2, "gene_assigns_minqv": 30, "matrix_min_genes": 1, "matrix_min_cells": 1, "matrix_max_mito": 100,
              "expected_cells": expected_cells, "estimate_cell_count": True, "mito_prefix": "MT-",
              "ref_genome_dir": "/refs/refdata-gex-demo-2024-A", "wf": {"container_sha": "sha0000"}}
    _w(os.path.join(top, "params.json"), json.dumps(params, indent=4))
    _trace(os.path.join(top, "execution", "trace.txt"))
    adapter_pct = {"full_length": 98.7, "single_adapter1": 0.4, "double_adapter1": 0.1, "single_adapter2": 0.5, "double_adapter2": 0.1,
                   "no_adapters": 0.1, "other": 0.1}
    reported_sat = 1.0 - truth["molecules"] / truth["rows"]
    if config_stats:
        cs = {"reads": reads_in, "full_length": full, "stranded": reads_in, "plus": reads_in // 2, "minus": reads_in - reads_in // 2,
              "single_adapter1": 11, "single_adapter2": 12, "double_adapter1": 1, "double_adapter2": 2, "no_adapters": 3, "other": 4}
        _w(os.path.join(sdir, f"{alias}.config_stats.json"), json.dumps(cs))
    _report_html(os.path.join(top, REPORT_NAME), version="3.3.4", v3=True, alias=alias,
                 numbers={"input": reads_in, "units": len(truth["per_unit"]), "mapped": reads_in - 9, "unmapped": 9, "full_length": full,
                          "valid_barcode": truth["rows"], "gene_assigned": truth["gene_assigned"], "transcript_assigned": truth["transcript_assigned"]},
                 params={k: str(v) for k, v in params.items() if k != "wf"}, adapter_pct=adapter_pct, saturation=(truth["rows"], reported_sat))
    link_expected = []
    if bam:
        refs = [("chr1", 400000), ("chr2", 300000), ("chrM", 16000)]
        pattern = [i for i, r in enumerate(runs) for _ in range(r[3])]
        recs = []
        for i in range(1500):
            fc, prot, acq, _ = runs[pattern[i % len(pattern)]]
            ref = "chr1" if i % 3 else "chr2"
            recs.append({"name": f"q{i}_0", "ref": ref, "pos": 20 + (i * 197) % (250000 if ref == "chr1" else 200000),
                         "tags": {"fn": f"{fc}_{prot}_{acq}_{i % 9}.pod5", "st": f"2026-01-0{1 + i % 2}T10:00:{i % 60:02d}.000+00:00", "CB": "X"}})
        recs += [{"name": f"mm{i}_0", "ref": "chrM", "pos": 3 + 5 * i, "tags": {"fn": f"{runs[0][0]}_{runs[0][1]}_{runs[0][2]}_1.pod5"}} for i in range(60)]
        _bam(os.path.join(sdir, f"{alias}.tagged.bam"), recs, refs)
        total = sum(r[3] for r in runs)
        link_expected = [{"flow_cell_id": r[0], "protocol_run_id_prefix": r[1], "acquisition_run_id_prefix": r[2], "share": r[3] / total} for r in runs]
    summary = _wf_summary(truth, config, dict(umis), dict(genes))
    mfrac = np.array([mito[b] / umis[b] for b in order])
    summary.update({"barcodes_cells": len(order), "barcodes_genes_detected": len(names), "barcodes_transcripts_detected": 52,
                    "mito_fraction_median": float(np.median(mfrac)), "mito_umi_fraction": sum(mito.values()) / truth["molecules"],
                    "saturation_reported": reported_sat, "reads_report_unmapped": 9,
                    "barcodes_reads_per_cell_mean": reads_in / (len(truth["per_unit"]) if hd else len(order))})
    if bam:
        summary["reads_chrm_fraction"] = 60 / 1560
    expected = {"label": alias, "assay": "visium_hd" if hd else "singlecell", "pipeline_version": "3.3.4", "reference": "demo-2024-A",
                "summary": summary, "molecule_hist": truth["molecule_hist"], "links": link_expected, "from_report": not config_stats,
                "adapter_pct": adapter_pct, "per_unit": dict(truth["per_unit"]),
                "sketch_ids": sorted({r[0].replace("-", "")[:16] for r in rows if r[0].startswith("0000")}),
                "umis_by_unit": dict(umis), "genes_by_unit": dict(genes)}
    return FakeAnalysis(top, "wf-single-cell", "wf_v3", expected, {"sample": sdir, "top": top})


def make_sibling_pair(out_dir: str, *, name: str = "demo_analysis", seed: int = 5) -> tuple[FakeAnalysis, FakeAnalysis]:
    """A Visium HD sample analysed by both pipelines, side by side: ``<out_dir>/<name>/{spaceranger/outs, wf-single-cell}``.
    The workflow's ``params.json`` names the Space Ranger BAM by the path the pipeline was started with."""
    base = os.path.join(out_dir, name)
    sr = make_spaceranger_hd(base, seed=seed)
    wf = make_wf_v3(base, seed=seed + 1, spaceranger_bam=f"./{name}/spaceranger/outs/possorted_genome_bam.bam", out_rel=f"./{name}/wf-single-cell",
                    n_reads=14000)
    return wf, sr


# ---------------------------------------------------------------------------------------------------------
# a tree for library tests: runs and the analyses made from them
# ---------------------------------------------------------------------------------------------------------

@dataclass
class SpatialTree:
    base: str
    roots: dict = field(default_factory=dict)            # short name -> folder
    runs: dict = field(default_factory=dict)             # label -> FixtureRun
    analyses: dict = field(default_factory=dict)         # label -> FakeAnalysis
    late_root: str = ""                                  # a folder outside the roots that holds the second flow cell's run


def make_spatial_tree(out_dir: str, *, seed: int = 21) -> SpatialTree:
    """Two roots with one run and four analyses, and a run kept outside them:

    ``root_seq/proj_spatial/demo_spatial_a/``      a MinKNOW run folder (flow cell TST00001), and under ``results/`` a
                                                    workflow 2.x output folder of the same sample with a Space Ranger
                                                    output inside its sample folder (the older lay-out)
    ``root_store/proj_hd/demo_analysis/``           a Visium HD sample analysed by both pipelines side by side; the reads of
                                                    the workflow's BAM name the run above and a run on flow cell TST00077
    ``root_store/proj_hd/copy/spaceranger/outs``    a copy of that Space Ranger output: the same analysis, a second path
    ``late/``                                       the run on TST00077, to be given to a later scan
    """
    import shutil

    from flowcellar.testing.fixtures import RunSpec, make_run
    from flowcellar.testing.trees import age_tree

    base = os.fspath(out_dir)
    t = SpatialTree(base=base, roots={"seq": os.path.join(base, "root_seq"), "store": os.path.join(base, "root_store")},
                    late_root=os.path.join(base, "late"))
    run_a = make_run(t.roots["seq"], RunSpec(flow_cell_id="TST00001", experiment="proj_spatial", sample="demo_spatial_a",
                                             start="2026-03-02T09:00:00Z", hours=6, reads=1500, seed=seed))
    run_b = make_run(t.late_root, RunSpec(flow_cell_id="TST00077", experiment="proj_hd", sample="demo_hd_day2",
                                          start="2026-03-04T09:00:00Z", hours=6, reads=1500, seed=seed + 1))
    t.runs = {"a": run_a, "b": run_b}
    sample_dir = os.path.dirname(run_a.run_dir)
    sr_inner_dir = os.path.join(sample_dir, "results", "wf-out", "demo_spatial_a")
    sr_v = make_spaceranger_visium(sr_inner_dir, seed=seed + 2, sample="demo_spatial_a", antisense=0.2)
    wf2 = make_wf_v2(os.path.join(sample_dir, "results"), seed=seed + 3, alias="demo_spatial_a", barcodes=sr_v.expected["barcodes"][:60])
    hd_base = os.path.join(t.roots["store"], "proj_hd", "demo_analysis")
    sr_hd = make_spaceranger_hd(hd_base, seed=seed + 4)
    short = lambda r: r.run_key.split("_")[1]             # noqa: E731
    wf3 = make_wf_v3(hd_base, seed=seed + 5, spaceranger_bam="./demo_analysis/spaceranger/outs/possorted_genome_bam.bam",
                     out_rel="./demo_analysis/wf-single-cell", n_reads=14000,
                     runs=(("TST00001", short(run_a), "11112222", 3), ("TST00077", short(run_b), "33334444", 2)))
    shutil.copytree(os.path.dirname(sr_hd.path), os.path.join(t.roots["store"], "proj_hd", "copy", "spaceranger"))
    t.analyses = {"sr_visium": sr_v, "wf_v2": wf2, "sr_hd": sr_hd, "wf_v3": wf3}
    age_tree(base)
    return t
