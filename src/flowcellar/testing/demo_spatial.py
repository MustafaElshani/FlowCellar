"""The single-cell and spatial analysis folders of the FlowCellar example (``flowcellar.testing.demo``).

Two invented samples, each analysed with Space Ranger and with wf-single-cell, side by side::

    <analysis>/visium_a/spaceranger/outs          a Visium slide of 55 um spots (one flow cell)
    <analysis>/visium_a/wf-single-cell            the workflow's output for the same sample
    <analysis>/visium_hd_b/spaceranger/outs       a Visium HD slide
    <analysis>/visium_hd_b/wf-single-cell         the workflow's output; the reads in its BAM name two flow cells

The Space Ranger folders hold a slide with a shape: sections of tissue with denser and thinner parts, written in the
layouts ``flowcellar.testing.sc`` knows, at a depth that follows the reads of the flow cells (so ``--scale`` scales them
too). The workflow folders are written from the same slide: molecules fall where the tissue is and are read a few
times each, so the two pipelines of a sample agree as they do on real data. Their per-read tables hold one row per
read (some 12 GB for the Visium HD sample at full size).
"""
from __future__ import annotations

import json
import os

import numpy as np

from flowcellar.testing import sc

WF_MAX_READS = 130_000_000


def _field(rng, shape: tuple[int, int], coarse: int) -> np.ndarray:
    """A smooth random surface: values drawn on a coarse grid, joined by straight lines along both axes."""
    g = rng.normal(0.0, 1.0, (coarse, coarse))
    ys, xs = np.linspace(0, coarse - 1, shape[0]), np.linspace(0, coarse - 1, shape[1])
    rows = np.stack([np.interp(xs, np.arange(coarse), g[i]) for i in range(coarse)])
    return np.stack([np.interp(ys, np.arange(coarse), rows[:, j]) for j in range(shape[1])], axis=1)


def tissue_and_density(rng, y: np.ndarray, x: np.ndarray, sections) -> tuple[np.ndarray, np.ndarray]:
    """(under tissue, relative density of molecules) at the points ``y``, ``x`` of the unit square: a few rounded
    sections with an uneven edge, denser towards their middle, with thinner patches inside."""
    shape = y.shape
    edge, inner = _field(rng, shape, 9), _field(rng, shape, 22)
    best = np.full(shape, -1.0)
    for cy, cx, ry, rx in sections:
        d = 1.0 - np.sqrt(((y - cy) / ry) ** 2 + ((x - cx) / rx) ** 2)
        best = np.maximum(best, d)
    tissue = (best + 0.16 * edge + 0.05 * inner) > 0.0
    lumen = (inner > 1.2) & (best > 0.12)                           # holes inside the tissue
    tissue &= ~lumen
    density = np.clip(0.55 + 0.9 * np.clip(best, 0, 1) + 0.22 * inner, 0.15, 2.2)
    return tissue, density / density[tissue].mean() if tissue.any() else density


def make_visium(out_dir: str, *, seed: int, sample: str, slide: str, reads: int, folder: str = "spaceranger") -> sc.FakeAnalysis:
    """``<out_dir>/<folder>/outs`` of a Visium slide with all its 4,992 spots, written from arrays: summary, web summary,
    molecules, raw matrix, spot positions. ``reads`` is the number of reads the slide was sequenced with."""
    rng = np.random.default_rng(seed)
    top = os.path.join(out_dir, folder)
    outs = os.path.join(top, "outs")
    n_rows, n_cols, n_genes = 78, 64, 6000
    names = sc.gene_names(n_genes + 3)
    r, c = np.mgrid[0:n_rows, 0:n_cols]
    r, c = r.ravel(), c.ravel()
    acol = 2 * c + (r % 2)
    under, dens = tissue_and_density(rng, (r / (n_rows - 1)).reshape(n_rows, n_cols), (acol / (2 * n_cols - 1)).reshape(n_rows, n_cols),
                                     ((0.36, 0.34, 0.27, 0.25), (0.66, 0.68, 0.24, 0.22)))
    under, dens = under.ravel(), dens.ravel()
    letters = np.frombuffer(b"ACGT", np.uint8)[rng.integers(0, 4, (len(r), 16), dtype=np.uint8)]
    bcs = [b.decode("ascii") for b in letters.view("S16").ravel()]
    feature_reads_goal = reads * 0.72
    saturation_goal = 0.72
    mol_goal = feature_reads_goal * (1.0 - saturation_goal)
    per_spot = mol_goal * 0.93 / max(int(under.sum()), 1)
    expected = np.where(under, per_spot * dens * rng.lognormal(0.0, 0.22, len(r)), per_spot * 0.035 * rng.lognormal(0.0, 0.6, len(r)))
    n_mol = rng.poisson(expected)
    weights = rng.lognormal(0.0, 1.9, n_genes + 3)
    weights[-3:] = weights[:-3].sum() * np.array([0.022, 0.019, 0.012])           # the three mitochondrial genes
    counts = rng.multinomial(n_mol, weights / weights.sum()).astype(np.int32)       # spots x genes
    in_matrix = np.nonzero(n_mol > 0)[0]
    sub = counts[in_matrix]
    spot_i, gene_i = np.nonzero(sub)
    indptr = np.concatenate([[0], np.cumsum(np.bincount(spot_i, minlength=len(in_matrix)))])
    sc._matrix_h5_csc(os.path.join(outs, "raw_feature_bc_matrix.h5"), names, [bcs[i] + "-1" for i in in_matrix], indptr, gene_i, sub[spot_i, gene_i],
                      version="spaceranger-3.0.1", chemistry="Visium V1 Slide")
    mol_spot = np.repeat(np.arange(len(in_matrix)), n_mol[in_matrix])
    mol_count = rng.geometric(1.0 - saturation_goal, len(mol_spot))
    feature_reads = int(mol_count.sum())
    reads_under = int(mol_count[under[in_matrix][mol_spot]].sum())
    raw_reads = int(max(reads, feature_reads / 0.74))
    metrics = {"cellranger_version": "spaceranger-3.0.1", "chemistry_description": "Visium V1 Slide", "chemistry_name": "SPATIAL3Pv1",
               "slide_serial_capture_area": slide, "analysis_parameters": {"include_introns": False},
               "libraries": {"0": {"raw_read_pairs": raw_reads, "feature_read_pairs": feature_reads, "usable_read_pairs": reads_under}}}
    sc._molecule_info(os.path.join(outs, "molecule_info.h5"), mol_count, mol_spot, [bcs[i] for i in in_matrix], metrics)
    lines = ["barcode,in_tissue,array_row,array_col,pxl_row_in_fullres,pxl_col_in_fullres"]
    lines += [f"{bcs[i]}-1,{int(under[i])},{r[i]},{acol[i]},{1000 + 100 * r[i]},{1000 + 50 * acol[i]}" for i in range(len(r))]
    sc._w(os.path.join(outs, "spatial", "tissue_positions.csv"), "\n".join(lines) + "\n")
    sc._w(os.path.join(outs, "spatial", "scalefactors_json.json"), json.dumps({"spot_diameter_fullres": 89.4, "tissue_hires_scalef": 0.15,
                                                                                 "tissue_lowres_scalef": 0.045, "fiducial_diameter_fullres": 144.5}))
    genes_spot = (counts > 0).sum(axis=1)
    n_under = int(under.sum())
    saturation = 1.0 - len(mol_count) / max(feature_reads, 1)
    values = {
        "Sample ID": sample, "Number of Spots Under Tissue": n_under, "Number of Reads": raw_reads,
        "Mean Reads per Spot": raw_reads / max(n_under, 1), "Mean Reads Under Tissue per Spot": reads_under / max(n_under, 1),
        "Fraction of Spots Under Tissue": n_under / len(r), "Median Genes per Spot": float(np.median(genes_spot[under])),
        "Median UMI Counts per Spot": float(np.median(n_mol[under])), "Valid Barcodes": 0.9612, "Valid UMIs": 0.9997,
        "Sequencing Saturation": saturation, "Q30 Bases in Barcode": 0.871, "Q30 Bases in RNA Read": 0.842, "Q30 Bases in UMI": 0.866,
        "Reads Mapped to Genome": 0.951, "Reads Mapped Confidently to Genome": 0.912,
        "Reads Mapped Confidently to Intergenic Regions": 0.041, "Reads Mapped Confidently to Intronic Regions": 0.062,
        "Reads Mapped Confidently to Exonic Regions": 0.809, "Reads Mapped Confidently to Transcriptome": feature_reads / raw_reads,
        "Reads Mapped Antisense to Gene": 0.018, "Fraction Reads in Spots Under Tissue": reads_under / max(feature_reads, 1),
        "Total Genes Detected": int((counts[under].sum(axis=0) > 0).sum()),
    }
    sc._metrics_csv(os.path.join(outs, "metrics_summary.csv"), values)
    sample_table = {"Sample ID": sample, "Sample Description": "", "Chemistry": "Visium V1 Slide - 3'", "Slide Serial Number": slide,
                    "Transcriptome": "demo-2024-A", "Pipeline Version": "spaceranger-3.0.1", "Image Reorientation": "On",
                    "Loupe Manual Alignment": "None"}
    sc._web_summary(os.path.join(outs, "web_summary.html"), layout="old", sample=sample_table, alarms=[],
                    diagnostics={"tso_frac": 0.0003, "low_support_umi_reads_frac": 0.003})
    sc._w(os.path.join(top, "_versions"), json.dumps({"martian": "v4.0.12", "pipelines": "spaceranger-3.0.1"}))
    sc._w(os.path.join(top, "_cmdline"), f"spaceranger count --id {sample} --transcriptome ../../refs/refdata-demo-2024-A --create-bam false "
                                         f"--slide {slide[:-3]} --area {slide[-2:]} --image ../images/{sample}.tif --output-dir ./{folder} --localcores 16")
    return sc.FakeAnalysis(outs, "spaceranger", "spaceranger", {"barcodes": bcs, "under": under, "reads": raw_reads, "weight": expected},
                           {"outs": outs})


def hd_layout(seed: int, n: int) -> tuple[np.ndarray, np.ndarray]:
    """(under tissue, relative density) of the example's Visium HD slide on its grid of ``n`` x ``n`` squares: the same
    for the Space Ranger folder and the workflow folder of the sample."""
    rng = np.random.default_rng(seed)
    yy, xx = np.mgrid[0:n, 0:n] / float(n - 1)
    return tissue_and_density(rng, yy, xx, HD_SECTIONS)


def make_hd(out_dir: str, *, seed: int, sample: str, slide: str, reads: int, n: int) -> sc.FakeAnalysis:
    """A Visium HD ``outs`` folder of ``n`` x ``n`` squares of 2 um, at a depth that gives about ``reads`` reads."""
    rng = np.random.default_rng(seed + 1000)
    tissue, dens = hd_layout(seed, n)
    squares = int(tissue.sum())
    per_square = min(reads / 4.6 / max(squares, 1), 12.0)                       # molecules per square: a read in about 4.6 gives a new one
    grey = np.clip(236.0 - np.where(tissue, 70.0 + 60.0 * dens, 0.0) + rng.normal(0.0, 6.0, (n, n)), 0, 255).astype(np.uint8)
    bins8 = max(squares // 16, 1)
    umis8 = per_square * 16
    area_mm2 = squares * 4e-6
    cells = int(area_mm2 * 1e6 / 190.0)
    total_umis = per_square * squares
    overrides = {
        "Total Genes Detected": 17_412, "Total Genes Detected Under Tissue 2 µm": 17_208, "Total Genes Detected Under Tissue 8 µm": 17_208,
        "Total Genes Detected Under Tissue 16 µm": 17_210, "Mean Genes Under Tissue per Square 2 µm": per_square * 0.94,
        "Number of Bins Under Tissue 16 µm": bins8 // 4, "Mean Reads Under Tissue per Bin 16 µm": reads / max(bins8 // 4, 1),
        "Fraction of Bins Under Tissue 16 µm": squares / float(n * n), "Mean Genes Under Tissue per Bin 16 µm": umis8 * 4 * 0.79,
        "Mean UMIs Under Tissue per Bin 16 µm": umis8 * 4, "Number of Cells": cells, "Reads in Cells": 0.684, "UMIs in Cells": 0.731,
        "Mean Reads per Cell": reads * 0.684 / max(cells, 1), "Median Genes per Cell": round(total_umis * 0.731 / max(cells, 1) * 0.74),
        "Median UMIs per Cell": round(total_umis * 0.731 / max(cells, 1) * 0.86), "Median Cell Area (μm²)": 171.0, "Median Nucleus Area (μm²)": 58.0,
        "UMIs per sq mm of Tissue": total_umis / max(area_mm2, 1e-9), "Reads per sq mm of Tissue": reads / max(area_mm2, 1e-9),
    }
    return sc.make_spaceranger_hd(out_dir, seed=seed, sample=sample, slide=slide, n=n, registration_alarm=False, tissue=tissue,
                                  density=per_square * dens, grey=grey, n_genes=17_415, gene_share=(0.80, 0.93), summary_values=overrides,
                                  reads_per_molecule=3.4)


def hd_grid(scale: float) -> int:
    """Squares of 2 um along one side of the example's Visium HD slide: the whole 6.5 mm from a tenth of full size on,
    a smaller capture area below (so that a small example stays small)."""
    return 3250 if scale >= 0.1 else int(np.clip(3250 * np.sqrt(scale * 10), 1000, 3250))


def _digits(values: np.ndarray, width: int) -> np.ndarray:
    """Whole numbers as zero-padded text, one byte per digit (rows x width)."""
    out = np.empty((len(values), width), np.uint8)
    for k in range(width):
        out[:, width - 1 - k] = (values // 10 ** k) % 10 + 48
    return out


def _mtx_arrays(folder: str, names: list[str], barcodes: list[str], feature: np.ndarray, column: np.ndarray, value: np.ndarray) -> None:
    """A matrix folder as the workflow writes it (``sc._mtx``), from three arrays sorted by feature and column."""
    import pyarrow as pa
    import pyarrow.csv as pacsv
    os.makedirs(folder, exist_ok=True)
    head = ("%%MatrixMarket matrix coordinate integer general\n" '%metadata_json: {"software_version": "ont-single-cell","format_version": 2}\n'
            f"{len(names)} {len(barcodes)} {len(value)}\n")
    with pa.CompressedOutputStream(pa.OSFile(os.path.join(folder, "matrix.mtx.gz"), "wb"), "gzip") as out:
        out.write(head.encode("ascii"))
        if len(value):
            pacsv.write_csv(pa.table({"f": feature + 1, "b": column + 1, "v": value}), out,
                            write_options=pacsv.WriteOptions(delimiter=" ", include_header=False, quoting_style="none"))
    sc._gz(os.path.join(folder, "barcodes.tsv.gz"), "".join(b + "\n" for b in barcodes))
    sc._gz(os.path.join(folder, "features.tsv.gz"), "".join(f"unknown_{i:05d}\t{n}\tGene Expression\n" for i, n in enumerate(names)))


def make_wf(out_dir: str, *, seed: int, alias: str, layout: str, reads: int, place_weight: np.ndarray, place_names, runs=(),
            folder: str = "wf-single-cell", n_genes: int = 17_412, grid: int | None = None) -> sc.FakeAnalysis:
    """A wf-single-cell output folder written from arrays, so that it can hold millions of reads: the per-read table,
    the gene matrix, the report and, for ``layout="v3"``, ``params.json`` and a tagged BAM whose reads name ``runs``.

    ``layout="v2"`` is a 55 um slide (``place_names`` are the spot barcodes), ``"v3"`` a Visium HD slide (``grid`` squares
    of 2 um along one side; a place is a square). ``place_weight`` says where molecules fall. Molecules are drawn first
    and each is then read a few times, as in a library that was amplified: the table holds more reads than molecules."""
    import pyarrow as pa
    import pyarrow.compute as pc
    import pyarrow.csv as pacsv
    rng = np.random.default_rng(seed)
    hd = layout == "v3"
    top = os.path.join(out_dir, folder)
    sdir = os.path.join(top, alias)
    os.makedirs(sdir, exist_ok=True)
    names = sc.gene_names(n_genes + 3)
    n_names = len(names)
    rows_goal = int(reads * 0.86)                                    # reads that are full length and got a place on the slide
    per_molecule = 3.4
    n_mol = max(int(rows_goal / per_molecule), 10)
    cdf = np.cumsum(place_weight.ravel().astype(np.float64))
    place = np.searchsorted(cdf, rng.random(n_mol) * cdf[-1]).astype(np.int64)
    gw = rng.lognormal(0.0, 1.9, n_names)
    gw[-3:] = gw[:-3].sum() * np.array([0.022, 0.019, 0.012])
    gene = np.searchsorted(np.cumsum(gw / gw.sum()), rng.random(n_mol)).clip(0, n_names - 1).astype(np.int64)
    gene[rng.random(n_mol) < 0.09] = -1                              # no gene could be given
    umi_len = 9 if hd else 12
    umi = rng.integers(0, 4 ** 9, n_mol, dtype=np.int64)
    times = rng.geometric(1.0 / per_molecule, n_mol)
    order = rng.permutation(int(times.sum()))
    mol_of_read = np.repeat(np.arange(n_mol), times)[order]
    n = len(mol_of_read)
    no_umi = rng.random(n) < 0.02
    chrom_of_gene = np.array([f"chr{1 + g % 22}" for g in range(n_names)] + ["chr1"])
    r_gene = gene[mol_of_read]
    by_chrom = np.argsort((np.where(r_gene >= 0, r_gene, 0) % 22), kind="stable")     # the workflow writes chromosome by chromosome
    mol_of_read, no_umi, r_gene = mol_of_read[by_chrom], no_umi[by_chrom], r_gene[by_chrom]
    table = os.path.join(sdir, f"{alias}.read_summary.tsv" if hd else "read_tags.tsv")
    gene_names_arr = pa.array(names + ["-"])
    chrom_arr = pa.array(chrom_of_gene)
    acgt = np.frombuffer(b"ACGT", np.uint8)
    hexc = np.frombuffer(b"0123456789abcdef", np.uint8)
    spot_names = None if hd else pa.array(list(place_names))
    with open(table, "wb") as sink:
        sink.write(sc.HEADER.encode("ascii"))
        with pacsv.CSVWriter(sink, pa.schema([(c, pa.string()) for c in sc.HEADER.strip().split("\t")]),
                             write_options=pacsv.WriteOptions(delimiter="\t", include_header=False, quoting_style="none")) as w:
            for lo in range(0, n, 1 << 20):
                hi = min(n, lo + (1 << 20))
                m = hi - lo
                mol = mol_of_read[lo:hi]
                h = hexc[rng.integers(0, 16, (m, 32), dtype=np.uint8)]
                rid = np.empty((m, 38), np.uint8)
                rid[:, [8, 13, 18, 23]] = 45
                rid[:, 0:8], rid[:, 9:13], rid[:, 14:18], rid[:, 19:23], rid[:, 24:36] = h[:, 0:8], h[:, 8:12], h[:, 12:16], h[:, 16:20], h[:, 20:32]
                rid[:, 36], rid[:, 37] = 95, 48                       # "_0"
                if hd:
                    bc = np.empty((m, 21), np.uint8)
                    bc[:, 0:8] = np.frombuffer(b"s_002um_", np.uint8)
                    bc[:, 8:13], bc[:, 13], bc[:, 14:19] = _digits(place[mol] // grid, 5), 95, _digits(place[mol] % grid, 5)
                    bc[:, 19:21] = np.frombuffer(b"-1", np.uint8)
                    corrected = pa.array(bc.view("S21").ravel()).cast(pa.string())
                    raw_bc = pa.repeat(pa.scalar("ACGT"), m)
                else:
                    corrected = spot_names.take(pa.array(place[mol]))
                    raw_bc = corrected
                u = np.empty((m, umi_len), np.uint8)
                u[:] = 65
                code = umi[mol]
                for k in range(9):
                    u[:, k] = acgt[(code >> (2 * k)) & 3]
                umis = pc.if_else(pa.array(no_umi[lo:hi]), pa.scalar("-"), pa.array(u.view(f"S{umi_len}").ravel()).cast(pa.string()))
                g = r_gene[lo:hi]
                g_names = gene_names_arr.take(pa.array(np.where(g >= 0, g, n_names)))
                has_tx = pa.array((g >= 0) & (rng.random(m) < 0.55))
                tx = pc.if_else(has_tx, pc.binary_join_element_wise(pa.scalar("TX_"), g_names, pa.scalar("_1"), ""), pa.scalar("-"))
                four = pa.repeat(pa.scalar("IIII"), m)
                w.write_table(pa.table([pa.array(rid.view("S38").ravel()).cast(pa.string()), raw_bc, corrected, four, umis, umis, four, g_names, tx,
                                        pa.repeat(pa.scalar("100"), m), pa.repeat(pa.scalar("900"), m),
                                        chrom_arr.take(pa.array(np.where(g >= 0, g, n_names)))], names=sc.HEADER.strip().split("\t")))
    # ---- what the workflow counted: molecules with a gene and a UMI, per unit (an 8 um bin, or a spot)
    seen = np.zeros(n_mol, bool)
    seen[mol_of_read[~no_umi]] = True
    counted = seen & (gene >= 0)
    if hd:
        g8 = -(-grid // 4)
        unit = (place // grid // 4) * g8 + (place % grid) // 4
    else:
        unit = place
    units_with_reads = np.unique(unit[np.unique(mol_of_read)])
    key = unit[counted] * n_names + gene[counted]
    pairs, molecules = np.unique(key, return_counts=True)
    p_unit, p_gene = pairs // n_names, pairs % n_names
    units = np.unique(p_unit)
    col = np.searchsorted(units, p_unit)
    if hd:
        unit_names = [f"s_008um_{u // g8:05d}_{u % g8:05d}" for u in units]
    else:
        unit_names = [place_names[u] for u in units]
    by_fb = np.lexsort((col, p_gene))
    suffix8 = "_8um" if hd else ""
    prefix = f"{alias}." if hd else ""
    _mtx_arrays(os.path.join(sdir, f"{prefix}gene_raw_feature_bc_matrix{suffix8}"), names, [b + "-1" for b in unit_names], p_gene[by_fb], col[by_fb],
                molecules[by_fb])
    sc._mtx(os.path.join(sdir, f"{prefix}transcript_raw_feature_bc_matrix{suffix8}"), [f"TX{i}" for i in range(52)], [unit_names[0] + "-1"], {(0, 0): 1})
    umis_unit = np.bincount(col, weights=molecules, minlength=len(units))
    mito_unit = np.bincount(col[p_gene >= n_names - 3], weights=molecules[p_gene >= n_names - 3], minlength=len(units))
    mito_name = f"{alias}.gene_expression_mito_per_cell.tsv" if hd else "gene.expression.mito-per-cell.tsv"
    sc._w(os.path.join(sdir, mito_name), "CB\tmito_pct\n" + "".join(f"{b}\t{100.0 * mito_unit[i] / umis_unit[i]}\n" for i, b in enumerate(unit_names)))
    reads_in = int(n / 0.86)
    full = int(reads_in * 0.987)
    gene_assigned = int((r_gene >= 0).sum())
    sc._trace(os.path.join(top, "execution", "trace.txt"))
    expected = {"label": alias, "rows": n, "molecules": int(counted.sum()), "units": len(units), "reads_in": reads_in}
    if hd:
        sc._mtx(os.path.join(sdir, f"{alias}.gene_raw_feature_bc_matrix_2um"), names, ["s_002um_00000_00000-1"], {(0, 0): 1})
        params = {"bam": f"./{alias}/prep/{alias}.bam", "spaceranger_bam": f"./{alias}/spaceranger/outs/possorted_genome_bam.bam",
                  "adapter_stats": f"./{alias}/prep/configs.json", "out_dir": f"./{alias}/{folder}", "kit": "visium_hd:v1", "threads": 16,
                  "full_length_only": True, "min_read_qual": 10, "barcode_max_ed": 2, "barcode_min_ed_diff": 2, "gene_assigns_minqv": 30,
                  "matrix_min_genes": 1, "matrix_min_cells": 1, "matrix_max_mito": 100, "expected_cells": None, "estimate_cell_count": True,
                  "mito_prefix": "MT-", "ref_genome_dir": "/refs/refdata-gex-demo-2024-A", "wf": {"container_sha": "sha0000"}}
        sc._w(os.path.join(top, "params.json"), json.dumps(params, indent=4))
        adapter_pct = {"full_length": 98.7, "single_adapter1": 0.4, "double_adapter1": 0.1, "single_adapter2": 0.5, "double_adapter2": 0.1,
                       "no_adapters": 0.1, "other": 0.1}
        sc._report_html(os.path.join(top, sc.REPORT_NAME), version="3.3.4", v3=True, alias=alias,
                        numbers={"input": reads_in, "units": len(units_with_reads), "mapped": reads_in - int(reads_in * 0.004),
                                 "unmapped": int(reads_in * 0.004), "full_length": full, "valid_barcode": n, "gene_assigned": gene_assigned,
                                 "transcript_assigned": int(gene_assigned * 0.55)},
                        params={k: str(v) for k, v in params.items() if k != "wf"}, adapter_pct=adapter_pct,
                        saturation=(n, 1.0 - int(counted.sum()) / max(gene_assigned, 1)))
        # the tagged BAM: the file names of the raw data on its reads say which flow cells the reads came from
        refs = [("chr1", 400000), ("chr2", 300000), ("chrM", 16000)]
        pattern = [i for i, r in enumerate(runs) for _ in range(r[3])]
        rng.shuffle(pattern)
        recs = []
        for i in range(3000):
            fc, prot, acq, _ = runs[pattern[i % len(pattern)]]
            ref = "chr1" if i % 3 else "chr2"
            recs.append({"name": f"q{i}_0", "ref": ref, "pos": 20 + (i * 197) % (250000 if ref == "chr1" else 200000),
                         "tags": {"fn": f"{fc}_{prot}_{acq}_{i % 9}.pod5", "st": f"2026-10-2{7 + i % 2}T10:00:{i % 60:02d}.000+00:00", "CB": "X"}})
        recs += [{"name": f"mm{i}_0", "ref": "chrM", "pos": 3 + 5 * i, "tags": {"fn": f"{runs[0][0]}_{runs[0][1]}_{runs[0][2]}_1.pod5"}} for i in range(110)]
        sc._bam(os.path.join(sdir, f"{alias}.tagged.bam"), recs, refs)
    else:
        sc._mtx(os.path.join(sdir, "gene_processed_feature_bc_matrix"), names[:5], [b + "-1" for b in unit_names[:5]], {(0, 0): 1})
        sc._w(os.path.join(sdir, "whitelist.tsv"), "".join(b + "\n" for b in unit_names))
        stranded = int(reads_in * 0.98)
        config = {"reads": reads_in, "full_length": full, "stranded": stranded, "plus": int(stranded * 0.54), "minus": stranded - int(stranded * 0.54),
                  "full_len": full, "single_adapter2": int(reads_in * 0.005), "single_adapter1": int(reads_in * 0.004), "other": int(reads_in * 0.002),
                  "no_adapters": int(reads_in * 0.001), "double_adapter1": int(reads_in * 0.0006), "double_adapter2": int(reads_in * 0.0004)}
        sc._w(os.path.join(sdir, "config_stats.json"), json.dumps(config, indent=4))
        params = {"kit": "visium:v1", "full_length_only": "True", "matrix_min_genes": "200", "matrix_max_mito": "20", "mito_prefix": "MT-",
                  "expected_cells": "500", "ref_genome_dir": "../refs/refdata-gex-demo-2024-A", "out_dir": f"./{folder}"}
        sc._report_html(os.path.join(top, sc.REPORT_NAME), version="2.2.0", v3=False, alias=alias,
                        numbers={"input": reads_in, "units": len(units), "mapped": int(reads_in * 0.97), "unmapped": reads_in - int(reads_in * 0.97)},
                        params=params)
    return sc.FakeAnalysis(top, "wf-single-cell", "wf_v3" if hd else "wf_v2", expected, {"sample": sdir, "top": top})


HD_SECTIONS = ((0.25, 0.30, 0.19, 0.24), (0.33, 0.74, 0.25, 0.17), (0.72, 0.24, 0.20, 0.17), (0.71, 0.66, 0.17, 0.26))


def write_analysis(base: str, sample: str, parts: list[tuple], ids_of, *, scale: float = 0.01, seed: int = 31) -> tuple[str, str]:
    """One sample of the example under ``base``: ``visium_a`` or ``visium_hd_b``. ``parts`` lists the runs its reads
    came from as (run, share in per cent). Returns (Space Ranger outs, workflow folder)."""
    reads = int(sum(r.reads for r, _ in parts) * scale * 0.88)                  # what was left for the pipelines after basecalling
    top = os.path.join(base, sample)
    if sample == "visium_a":
        # one flow cell: the files name no run, the sample name links the analysis to it
        sr = make_visium(top, seed=seed, sample=sample, slide="V00T00-101-A1", reads=reads)
        wf = make_wf(top, seed=seed + 1, alias=sample, layout="v2", reads=min(reads, WF_MAX_READS), place_weight=sr.expected["weight"],
                     place_names=sr.expected["barcodes"], n_genes=6000)
        return sr.path, wf.path
    # two flow cells: the reads of the workflow's BAM say which, and how much of each
    n = hd_grid(scale)
    sr = make_hd(top, seed=seed + 2, sample=sample, slide="H1-DEMOSLD-D1", reads=reads, n=n)
    tissue, dens = hd_layout(seed + 2, n)
    wf = make_wf(top, seed=seed + 3, alias=sample, layout="v3", reads=min(reads, WF_MAX_READS), grid=n,
                 place_weight=np.where(tissue, dens, 0.012), place_names=None,
                 runs=tuple((r.flow_cell_id, ids_of(r)[0][:8], ids_of(r)[1][:8], w) for r, w in parts))
    return sr.path, wf.path
