"""What is read from a Space Ranger ``outs/`` folder, in the order: the summary files first (cheap, exact, always),
then the HDF5 files (per-unit distributions, saturation, the map), then the BAM index.

Every number of the web summary's tables is in ``metrics_summary.csv`` and is stored as the pipeline gave it.
What FlowCellar computes itself is checked against those numbers where both exist (``params.cross_checks``).
Three things are in no data file of ``outs/`` and are taken from the JSON object inside ``web_summary.html``:
the pipeline's own alarms, the diagnostics, and the sample table (chemistry, slide, alignment file).
"""
from __future__ import annotations

import numpy as np

from flowcellar.sc import axes, maps
from flowcellar.sc import spaceranger as sr
from flowcellar.sc.collect import Builder
from flowcellar.util.errors import FlowCellarError, InputError, MissingExtraError

_MICRO = str.maketrans({"μ": "µ"})          # the summary spells "µm" with two different letters


def _get(ms: dict, *names: str):
    for n in names:
        v = ms.get(n.translate(_MICRO))
        if v is not None:
            return v
    return None


CSV_SUMMARY = (
    ("sr_reads", ("Number of Reads",)),
    ("sr_valid_barcodes", ("Valid Barcodes",)),
    ("sr_valid_umis", ("Valid UMIs", "Valid UMI Sequences")),
    ("sr_saturation", ("Sequencing Saturation",)),
    ("sr_q30_barcode", ("Q30 Bases in Barcode",)),
    ("sr_q30_umi", ("Q30 Bases in UMI",)),
    ("sr_q30_rna_read", ("Q30 Bases in RNA Read",)),
    ("sr_mapped_genome", ("Reads Mapped to Genome",)),
    ("sr_mapped_confidently_genome", ("Reads Mapped Confidently to Genome",)),
    ("sr_mapped_intergenic", ("Reads Mapped Confidently to Intergenic Regions",)),
    ("sr_mapped_intronic", ("Reads Mapped Confidently to Intronic Regions",)),
    ("sr_mapped_exonic", ("Reads Mapped Confidently to Exonic Regions",)),
    ("sr_mapped_transcriptome", ("Reads Mapped Confidently to Transcriptome",)),
    ("sr_antisense", ("Reads Mapped Antisense to Gene",)),
    ("spatial_fraction_reads_under_tissue", ("Fraction Reads in Spots Under Tissue", "Fraction Reads in Squares Under Tissue")),
    ("spatial_genes_detected", ("Total Genes Detected",)),
    ("spatial_units_under_tissue", ("Number of Spots Under Tissue", "Number of Bins Under Tissue 8 µm")),
    ("spatial_fraction_under_tissue", ("Fraction of Spots Under Tissue", "Fraction of Bins Under Tissue 8 µm")),
    ("spatial_reads_per_unit", ("Mean Reads Under Tissue per Spot", "Mean Reads Under Tissue per Bin 8 µm")),
    ("spatial_reads_per_unit_all", ("Mean Reads per Spot",)),
    ("spatial_genes_per_unit", ("Median Genes per Spot", "Mean Genes Under Tissue per Bin 8 µm")),
    ("spatial_umis_per_unit", ("Median UMI Counts per Spot", "Mean UMIs Under Tissue per Bin 8 µm")),
    ("spatial_umis_per_mm2", ("UMIs per sq mm of Tissue",)),
    ("spatial_reads_per_mm2", ("Reads per sq mm of Tissue",)),
    ("spatial_squares_under_tissue_2um", ("Number of Squares Under Tissue 2 µm",)),
    ("spatial_cells", ("Number of Cells",)),
    ("spatial_reads_in_cells", ("Reads in Cells",)),
    ("spatial_umis_in_cells", ("UMIs in Cells",)),
    ("spatial_image_um_per_pixel", ("Microscope Image Pixel Size",)),
)
BIN_LEVELS = ((2, "Squares", "Square 2 µm", "2 µm"), (8, "Bins", "Bin 8 µm", "8 µm"), (16, "Bins", "Bin 16 µm", "16 µm"))
SEGMENTATION = (("cells", "Number of Cells"), ("reads_in_cells", "Reads in Cells"), ("umis_in_cells", "UMIs in Cells"),
                ("mean_reads_per_cell", "Mean Reads per Cell"), ("median_genes_per_cell", "Median Genes per Cell"),
                ("median_umis_per_cell", "Median UMIs per Cell"), ("median_cell_area_um2", "Median Cell Area (µm²)"),
                ("median_nucleus_area_um2", "Median Nucleus Area (µm²)"))
ASSAY_OF_CHEMISTRY = (("SPATIAL-HD", "visium_hd"), ("SPATIAL", "visium"))


def collect(b: Builder) -> None:
    claim, f = b.claim, b.claim.files
    ms_raw = sr.metrics_summary(f["metrics_summary"])
    ms = {k.translate(_MICRO): v for k, v in ms_raw.items()}
    b.input("spaceranger.metrics_summary", f["metrics_summary"], records=1, fingerprint=True)
    hd = claim.hd
    ident = b.identity
    ident.update(sample_id=str(ms.get("Sample ID") or ""), label=str(ms.get("Sample ID") or "") or None,
                 assay="visium_hd" if hd else "visium")
    b.params.update(pipeline="spaceranger", layout=claim.layout, hd=hd, summary_columns=len(ms_raw), sample_names=[ident["sample_id"]],
                    cross_checks={})
    for key, names in CSV_SUMMARY:
        b.put(key, _get(ms, *names))
    b.put("unit_name", "bin" if hd else "spot")
    b.put("unit_um", 8.0 if hd else 55.0)
    b.put("unit_stat", "mean" if hd else "median")          # what the pipeline reports per unit: a mean for bins, a median for spots
    for k in ("sr_reads", "spatial_units_under_tissue", "spatial_genes_detected", "spatial_squares_under_tissue_2um", "spatial_cells"):
        b.scalar(k, b.summary.get(k))
    if hd:
        _bin_levels(b, ms)
    _web_summary(b, ms)
    _pipestance(b, ms)
    h5 = sr.h5py_or_none() is not None
    b.params["h5py"] = h5
    if not h5:
        b.note("info", "h5py_missing", str(MissingExtraError("h5py", "sc", "The spatial map, the saturation curve and the per-unit "
                                                              "distributions of Space Ranger output").message)
               + ' Install it with: pip install "flowcellar[sc]"')
        for role in ("molecule_info", "feature_slice", "raw_matrix", "filtered_matrix_008um"):
            if role in f:
                b.input("spaceranger." + ("matrix" if "matrix" in role else role), f[role], state="skipped")
    else:
        _molecule_info(b, hd)
        try:
            if hd:
                _hd_units(b, ms)
            else:
                _visium_units(b, ms)
        except (FlowCellarError, OSError, KeyError, ValueError) as e:
            b.note("warning", "input_unreadable", f"The per-unit numbers and the map could not be read: {getattr(e, 'message', None) or e}")
    _bam(b, ms)
    _parameters(b)
    if ident.get("slide") is None and b.summary.get("slide_serial_capture_area"):
        ident["slide"] = b.summary["slide_serial_capture_area"]


# ---------------------------------------------------------------------------------------------------------
def _bin_levels(b: Builder, ms: dict) -> None:
    rows = []
    for um, plural, per, size in BIN_LEVELS:
        n = _get(ms, f"Number of {plural} Under Tissue {size}")
        if n is None:
            continue
        rows.append([um, n, _get(ms, f"Fraction of {plural} Under Tissue {size}"), _get(ms, f"Mean Reads Under Tissue per {per}"),
                     _get(ms, f"Mean UMIs Under Tissue per {per}"), _get(ms, f"Mean Genes Under Tissue per {per}"),
                     _get(ms, f"Total Genes Detected Under Tissue {size}")])
    if rows:
        b.table("spatial_bin_levels", [("bin_um", "int"), ("units_under_tissue", "int"), ("fraction_under_tissue", "float"),
                                       ("mean_reads", "float"), ("mean_umis", "float"), ("mean_genes", "float"),
                                       ("genes_detected", "int")], rows)
    seg = [[key, float(_get(ms, name))] for key, name in SEGMENTATION if isinstance(_get(ms, name), (int, float))]
    if seg:
        b.table("spatial_segmentation", [("metric", "str"), ("value", "float")], seg)


def _web_summary(b: Builder, ms: dict) -> None:
    path = b.claim.files.get("web_summary")
    if not path:
        return
    try:
        with b.timed("web_summary"):
            ws = sr.web_summary(path)
    except (InputError, OSError) as e:
        b.input("spaceranger.web_summary", path, state="error", note=getattr(e, "message", None) or str(e))
        b.note("warning", "report_unreadable", "The pipeline's alarms and its sample table could not be read from web_summary.html "
                                               "(a layout this version does not know). They are not available; every other number is.")
        b.put("spatial_alarms_known", False)
        return
    wid = b.input("spaceranger.web_summary", path)
    warn = err = 0
    for a in ws["alarms"]:
        level = "error" if a["level"] == "ERROR" else "warning" if a["level"] in ("WARN", "WARNING") else "info"
        warn += level == "warning"
        err += level == "error"
        b.note(level, "pipeline_alarm", a["message"] or a["title"], title=a["title"][:256],
               value=a["value"] if isinstance(a["value"], (int, float, str)) else None, source="spaceranger", input=wid)
    b.put("spatial_alarms_known", True)
    b.put("spatial_alarms_warn", int(warn), kind="embedded_json", input_id=wid)
    b.put("spatial_alarms_error", int(err), kind="embedded_json", input_id=wid)
    s = ws["sample"]
    ident = b.identity
    version = str(s.get("Pipeline Version") or "").replace("spaceranger-", "") or None
    ident.setdefault("pipeline_version", version)
    if s.get("Slide Serial Number"):
        ident["slide"] = str(s["Slide Serial Number"])
        b.put("slide_serial_capture_area", ident["slide"], kind="embedded_json", input_id=wid)
    if s.get("Transcriptome"):
        ident["reference"] = str(s["Transcriptome"])
    if s.get("Chemistry"):
        b.put("chemistry", str(s["Chemistry"]), kind="embedded_json", input_id=wid)
    rows = [[k, str(v)] for k, v in s.items() if k not in ("Sample ID", "Sample Description") and v not in (None, "")]
    if rows:
        b.table("spatial_sample", [("field", "str"), ("value", "str")], rows)
        b.sources["spatial_sample"] = {"input": wid, "kind": "embedded_json"}
    # image registration: what the sample table and the pipeline's own offset check say
    reg = []
    manual = s.get("Loupe Manual Alignment")
    if manual is not None:
        is_manual = str(manual).strip().lower() not in ("", "none")
        reg.append(["alignment", "manual alignment file" if is_manual else "automatic"])
        if is_manual:
            reg.append(["alignment_file", str(manual)])
    if s.get("Image Reorientation") is not None:
        reg.append(["image_reorientation", str(s["Image Reorientation"])])
    d = ws["diagnostics"] or {}
    off = d.get("hd_layout_offset") if isinstance(d.get("hd_layout_offset"), dict) else None
    if off:
        for key in ("x_offset", "y_offset", "gaussian_fit_x_offset", "gaussian_fit_y_offset", "has_gaussian_fit_failed",
                    "is_gaussian_fit_outside_range", "is_custom_offset"):
            if off.get(key) is not None:
                reg.append([key, str(off[key]).lower() if isinstance(off[key], bool) else f"{off[key]:g}"])
        b.put("spatial_registration_flag", bool(off.get("has_gaussian_fit_failed") or off.get("is_gaussian_fit_outside_range")),
              kind="embedded_json", input_id=wid)
    if reg:
        b.table("spatial_registration", [("field", "str"), ("value", "str")], reg)
        b.sources["spatial_registration"] = {"input": wid, "kind": "embedded_json"}
    for key, name in (("sr_corrected_bc_frac", "corrected_bc_frac"), ("sr_tso_frac", "tso_frac"),
                      ("sr_low_support_umi_reads_frac", "low_support_umi_reads_frac")):
        if isinstance(d.get(name), (int, float)):
            b.put(key, float(d[name]), kind="embedded_json", input_id=wid)
    b.params["command_line_in_report"] = bool(ws.get("command_line"))
    b._command_line = ws.get("command_line")


def _pipestance(b: Builder, ms: dict) -> None:
    """Diagnostics from the pipeline's own JSON when the pipestance was kept and describes this output."""
    path = b.claim.files.get("pipestance_metrics")
    if not path:
        return
    d = sr.pipestance_metrics(path, ms)
    if d is None:
        b.input("spaceranger.pipestance_metrics", path, state="skipped",
                note="older than outs/: its read count or unit count differs from metrics_summary.csv, so it is not used")
        return
    pid = b.input("spaceranger.pipestance_metrics", path)
    for key, name in (("sr_corrected_bc_frac", "corrected_bc_frac"), ("sr_tso_frac", "tso_frac"),
                      ("sr_low_support_umi_reads_frac", "low_support_umi_reads_frac")):
        if isinstance(d.get(name), (int, float)):
            b.put(key, float(d[name]), kind="data", input_id=pid)
    if d.get("tissue_registration_type"):
        t = b.tables.get("spatial_registration")
        row = ["registration_type", str(d["tissue_registration_type"])]
        if t:
            t["rows"].append(row)
        else:
            b.table("spatial_registration", [("field", "str"), ("value", "str")], [row])


def _molecule_info(b: Builder, hd: bool) -> None:
    path = b.claim.files.get("molecule_info")
    if not path:
        return
    with b.timed("molecule_info"):
        mi = sr.molecule_info(path, per_barcode=not hd)
    b.input("spaceranger.molecule_info", path, records=mi["molecules"])
    m = mi["metrics"]
    ident = b.identity
    name = str(m.get("chemistry_name") or "")
    for prefix, assay in ASSAY_OF_CHEMISTRY:
        if name.upper().startswith(prefix):
            ident["assay"] = assay
            break
    if m.get("chemistry_description"):
        b.put("chemistry", str(m["chemistry_description"]))
    b.put("chemistry_name", name or None)
    if m.get("slide_serial_capture_area"):
        ident["slide"] = str(m["slide_serial_capture_area"])
        b.put("slide_serial_capture_area", ident["slide"])
    if m.get("cellranger_version"):
        ident["pipeline_version"] = str(m["cellranger_version"]).replace("spaceranger-", "")
    params = m.get("analysis_parameters") or {}
    if "include_introns" in params:
        b.params["include_introns"] = bool(params["include_introns"])
    libs = [v for v in (m.get("libraries") or {}).values() if isinstance(v, dict)]
    if libs:
        raw, feat, usable = (sum(int(v.get(k) or 0) for v in libs) for k in ("raw_read_pairs", "feature_read_pairs", "usable_read_pairs"))
        b.counts["sr_read_funnel"] = {"cats": ["raw", "feature", "usable"], "n": [raw, feat, usable], "role": "funnel"}
        b.put("sr_reads_feature", feat)
        b.put("sr_reads_usable", usable)
        if raw:
            b.put("sr_usable_reads_fraction", usable / raw)
    b.scalar("barcodes_unique_molecules", mi["molecules"])
    b.put("barcodes_unique_molecules", mi["molecules"])
    b.scalar("saturation_molecule_reads", mi["reads"])
    b.put("saturation_molecule_reads", mi["reads"])
    b.hists["saturation_reads_per_molecule"] = axes.reads_per_molecule_hist(mi["count_hist"])
    if mi["reads"]:
        b.params["cross_checks"]["saturation_from_molecules"] = 1.0 - mi["molecules"] / mi["reads"]
    b._molecule_info = mi


def _visium_units(b: Builder, ms: dict) -> None:
    f = b.claim.files
    if "raw_matrix" not in f or "tissue_positions" not in f:
        return
    with b.timed("matrix"):
        mx = sr.matrix_h5(f["raw_matrix"])
    b.input("spaceranger.matrix", f["raw_matrix"], records=len(mx["barcodes"]))
    pos = sr.tissue_positions(f["tissue_positions"])
    b.input("spaceranger.tissue_positions", f["tissue_positions"], records=len(pos))
    names = list(pos)
    in_t = np.array([pos[n][0] for n in names], dtype=np.int64)
    rows = np.array([pos[n][1] for n in names], dtype=np.int64)
    cols = np.array([pos[n][2] for n in names], dtype=np.int64)
    where = {n: i for i, n in enumerate(names)}
    idx = np.array([where.get(n, -1) for n in mx["barcodes"]], dtype=np.int64)
    ok = idx >= 0

    def spread(values) -> np.ndarray:
        out = np.zeros(len(names), dtype=np.float64)
        out[idx[ok]] = np.asarray(values)[ok]
        return out
    umis, genes, mito = spread(mx["umis"]), spread(mx["genes"]), spread(mx["mito_umis"])
    t = in_t == 1
    b.scalar("spatial_spots", len(names))
    b.put("spatial_spots", len(names))
    b.scalar("total_umis", mx["total_umis"])
    b.hists["barcodes_umis_hist"] = axes.hist(axes.LOG_COUNTS, umis[t])
    b.hists["barcodes_genes_hist"] = axes.hist(axes.LOG_COUNTS, genes[t])
    b.stats("barcodes_umis", umis[t])
    b.stats("barcodes_genes", genes[t])
    has = t & (umis > 0)
    frac = np.divide(mito, umis, out=np.zeros_like(umis), where=umis > 0)
    b.hists["mito_fraction_hist"] = axes.hist(axes.FRACTION_50, frac[has])
    if has.any():
        b.put("mito_fraction_median", float(np.median(frac[has])))
        b.put("mito_fraction_p90", float(np.percentile(frac[has], 90)))
    if umis.sum() > 0:
        b.put("mito_umi_fraction", float(mito.sum() / umis.sum()))
        b.put("spatial_fraction_umis_under_tissue", float(umis[t].sum() / umis.sum()))
    b.params["mito_features"] = mx.get("n_mito_features")
    rk, rv = axes.rank_curve(umis[t])
    if rk:
        b.curves["barcodes_knee"] = {"x": rk, "y": {"under_tissue": rv}, "x_unit": "rank", "y_unit": "umis", "x_scale": "log10"}
    rk, rv = axes.rank_curve(umis[~t])
    if rk:
        b.curves["barcodes_knee_off_tissue"] = {"x": rk, "y": {"off_tissue": rv}, "x_unit": "rank", "y_unit": "umis", "x_scale": "log10"}
    seen = mx.get("features_seen")
    checks = b.params["cross_checks"]
    checks.update(units_under_tissue=int(t.sum()), umis_per_unit_median=float(np.median(umis[t])) if t.any() else None,
                  genes_per_unit_median=float(np.median(genes[t])) if t.any() else None)
    if seen is not None and "indices" in mx:
        under = np.zeros(len(mx["barcodes"]), dtype=bool)
        under[ok] = t[idx[ok]]
        per_entry = np.repeat(under, mx["genes"])
        checks["genes_detected_under_tissue"] = int(len(np.unique(mx["indices"][per_entry])))
    # reads per spot, and the saturation of each spot, from the molecule table
    mi = getattr(b, "_molecule_info", None)
    reads = sat = None
    if mi is not None and "reads_per_barcode" in mi:
        stem = {n.split("-")[0]: i for i, n in enumerate(names)}
        ridx = np.array([stem.get(n.split("-")[0], -1) for n in mi["barcodes"]], dtype=np.int64)
        rok = ridx >= 0
        reads = np.zeros(len(names), dtype=np.float64)
        mols = np.zeros(len(names), dtype=np.float64)
        reads[ridx[rok]] = mi["reads_per_barcode"][rok]
        mols[ridx[rok]] = mi["molecules_per_barcode"][rok]
        b.hists["barcodes_reads_hist"] = axes.hist(axes.LOG_COUNTS, reads[t])
        b.stats("barcodes_reads", reads[t])
        b.scalar("barcodes_detected", int((reads > 0).sum()))
        b.put("barcodes_detected", int((reads > 0).sum()))
        sat = 1.0 - np.divide(mols, reads, out=np.zeros_like(reads), where=reads > 0)
        if reads.sum() > 0:
            checks["fraction_reads_under_tissue"] = float(reads[t].sum() / reads.sum())
    b.grids["spatial_map"] = maps.visium_map(rows, cols, in_t, umis=umis, genes=genes, mito=frac, saturation=sat,
                                             sat_valid=reads if reads is not None else None)
    b.params["map"] = {"kind": "visium", "spots": int(len(names))}


def _hd_units(b: Builder, ms: dict) -> None:
    f = b.claim.files
    if "feature_slice" not in f:
        return
    genes8 = None
    genes = None
    mpath = f.get("filtered_matrix_008um")
    if mpath:
        with b.timed("matrix"):
            mx = sr.matrix_h5(mpath, want_data=False)
        b.input("spaceranger.matrix", mpath, records=len(mx["barcodes"]))
        r8, c8, _ = sr.hd_rows_cols(mx["barcodes"])
        genes = mx["genes"]
        genes8 = (r8, c8, genes)
    with b.timed("feature_slice"):
        with sr.FeatureSlice(f["feature_slice"]) as fs:
            grid, ex = maps.hd_map(fs, map_um=b.opts.map_um, genes8=genes8, tissue_image=b.opts.tissue_image)
            if b.identity.get("pipeline_version") is None and fs.software_version:
                b.identity["pipeline_version"] = fs.software_version
    b.input("spaceranger.feature_slice", f["feature_slice"])
    b.grids["spatial_map"] = grid
    b.params["map"] = {"kind": "visium_hd", "bin_um": grid["layout"]["bin_um"], "tissue_image": "image" in grid["channels"]}
    under = ex["tissue8"]
    U8, R8, M8 = ex["umis8"], ex["reads8"], ex["mito8"]
    b.scalar("total_umis", ex["umis_total"])
    b.scalar("spatial_umis_under_tissue", ex["umis_under_tissue"])
    b.scalar("mito_umis", ex["mito_umis"])
    b.params["mito_features"] = ex["mito_features"]
    if ex["umis_total"]:
        b.put("spatial_fraction_umis_under_tissue", ex["umis_under_tissue"] / ex["umis_total"])
        b.put("mito_umi_fraction", ex["mito_umis"] / ex["umis_total"])
    u = U8[under]
    b.hists["barcodes_umis_hist"] = axes.hist(axes.LOG_COUNTS, u)
    b.stats("barcodes_umis", u)
    if R8.sum() > 0:
        r = R8[under]
        b.hists["barcodes_reads_hist"] = axes.hist(axes.LOG_COUNTS, r)
        b.stats("barcodes_reads", r)
        b.scalar("barcodes_detected", int((R8 > 0).sum()))
        b.put("barcodes_detected", int((R8 > 0).sum()))
    if genes is not None:
        b.hists["barcodes_genes_hist"] = axes.hist(axes.LOG_COUNTS, genes)
        b.stats("barcodes_genes", genes)
    has = under & (U8 > 0)
    frac = np.divide(M8, U8, out=np.zeros_like(U8), where=U8 > 0)
    b.hists["mito_fraction_hist"] = axes.hist(axes.FRACTION_50, frac[has])
    if has.any():
        b.put("mito_fraction_median", float(np.median(frac[has])))
        b.put("mito_fraction_p90", float(np.percentile(frac[has], 90)))
    rk, rv = axes.rank_curve(u)
    if rk:
        b.curves["barcodes_knee"] = {"x": rk, "y": {"under_tissue": rv}, "x_unit": "rank", "y_unit": "umis", "x_scale": "log10"}
    rk, rv = axes.rank_curve(U8[~under])
    if rk:
        b.curves["barcodes_knee_off_tissue"] = {"x": rk, "y": {"off_tissue": rv}, "x_unit": "rank", "y_unit": "umis", "x_scale": "log10"}
    n_bins = int(under.sum())
    b.params["cross_checks"].update(
        squares_under_tissue=ex["squares_under_tissue"], units_under_tissue=n_bins,
        umis_per_unit_mean=float(u.mean()) if n_bins else None,
        genes_per_unit_mean=float(genes.mean()) if genes is not None and len(genes) else None,
        umis_per_square_mean=(ex["umis_under_tissue"] / ex["squares_under_tissue"]) if ex["squares_under_tissue"] else None,
        valid_barcode_reads=ex["reads_total"],
        fraction_valid_reads_under_tissue=(ex["reads_under_tissue"] / ex["reads_total"]) if ex["reads_total"] else None)


def _bam(b: Builder, ms: dict) -> None:
    f = b.claim.files
    if "bam" not in f or "bam_index" not in f:
        return
    try:
        with b.timed("bam_index"):
            m = sr.bam_mapping(f["bam"], f["bam_index"])
    except (InputError, OSError, ValueError, EOFError) as e:
        b.input("spaceranger.bam_index", f["bam_index"], state="error", note=getattr(e, "message", None) or str(e))
        return
    b.input("spaceranger.bam_index", f["bam_index"], records=m["mapped"] + m["unmapped"])
    b.scalar("reads_bam_mapped", m["mapped"])
    b.scalar("reads_bam_unmapped", m["unmapped"])
    if m["chrm_fraction"] is not None:
        b.scalar("reads_chrm_mapped", m["chrm_mapped"])
        b.put("reads_chrm_fraction", m["chrm_fraction"])
    slide = sr.slide_from_header(m["header_text"])
    if slide and not b.identity.get("slide"):
        b.identity["slide"] = slide
    total = m["mapped"] + m["unmapped"]
    b.params["cross_checks"].update(bam_records=total, bam_mapped_fraction=(m["mapped"] / total) if total else None)


def _parameters(b: Builder) -> None:
    f = b.claim.files
    text, kind, iid = None, "data", None
    if "cmdline" in f:
        text = sr.read_small(f["cmdline"]).strip()
        iid = b.input("spaceranger.cmdline", f["cmdline"])
    elif getattr(b, "_command_line", None):
        text, kind = b._command_line, "embedded_json"
        iid = b._input_ids.get(f.get("web_summary", ""))
    if not text:
        return
    rows = [[k, v] for k, v in sr.command_parameters(text) if k not in ("localcores", "localmem", "uiport", "output-dir", "jobmode")]
    if rows:
        b.table("workflow_parameters", [("key", "str"), ("value", "str")], rows)
        if kind != "data":
            b.sources["workflow_parameters"] = {"input": iid, "kind": kind}
