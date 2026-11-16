"""What is read from one sample folder of a wf-single-cell output, for both file layouts.

Order: the small files (parameters, adapter statistics, the workflow trace), the count matrices (streamed), the
mitochondrial table, the BAM header and index, a bounded look at the read tags of the BAM (which runs the reads
came from), and last the one step that takes minutes: the pass over the per-read table.

Four numbers of workflow 3.x exist only in its HTML report when the adapter statistics file was not kept: the
input reads, the full-length reads, the adapter layouts and the unmapped count. They are taken from the report
and marked; when the report cannot be read they are "not available" and everything else is still reported.
"""
from __future__ import annotations

import os
import time

import numpy as np

from flowcellar.sc import axes, detect, link, maps, readpass, wfsc
from flowcellar.sc import spaceranger as sr
from flowcellar.sc.collect import Builder
from flowcellar.util.errors import FlowCellarError, InputError

FUNNEL = ("input", "full_length", "valid_barcode", "gene_assigned", "transcript_assigned")


def collect(b: Builder) -> None:
    claim = b.claim
    ident = b.identity
    ident.update(sample_id=claim.alias or "", label=claim.alias)
    b.params.update(pipeline="wf-single-cell", layout=claim.layout, hd=claim.hd, alias=claim.alias, sample_names=[claim.alias or ""],
                    cross_checks={})
    params = _params(b)
    rep = _report(b)
    kit = params.get("kit") or ((rep or {}).get("parameters") or {}).get("kit")
    assay = wfsc.assay_of_kit(kit) or ("visium_hd" if claim.hd else None)
    ident["assay"] = assay
    b.put("kit", str(kit) if kit else None)
    unit = {"visium_hd": ("bin", 8.0), "visium": ("spot", 55.0)}.get(assay, ("cell", None))
    b.put("unit_name", unit[0])
    b.put("unit_um", unit[1])
    b.put("unit_stat", "median")
    mito_prefix = str(params.get("mito_prefix") or ((rep or {}).get("parameters") or {}).get("mito_prefix") or "MT-")
    b.params["mito_prefix"] = mito_prefix
    _config_stats(b, rep)
    mx = _matrices(b, mito_prefix)
    _mito(b, mx)
    _whitelist(b, mx)
    _trace(b)
    _bam(b)
    rp = _read_pass(b, rep, assay, params)
    _derived(b, assay)
    _map(b, mx, rp, assay)
    _sibling_hint(b, params)


# ---------------------------------------------------------------------------------------------------------
def _params(b: Builder) -> dict:
    path = b.claim.files.get("params")
    if not path:
        return {}
    try:
        p = wfsc.read_json(path)
    except (InputError, OSError) as e:
        b.input("wfsc.params", path, state="error", note=getattr(e, "message", None) or str(e))
        return {}
    b.input("wfsc.params", path, fingerprint=True)
    rows = [[k, v] for k, v in wfsc.parameters(p)]
    if rows:
        b.table("workflow_parameters", [("key", "str"), ("value", "str")], rows)
    if p.get("ref_genome_dir"):
        b.identity["reference"] = os.path.basename(str(p["ref_genome_dir"]).rstrip("/")).replace("refdata-gex-", "") or None
    if p.get("expected_cells") is not None:
        b.put("expected_cells", p["expected_cells"])
    return p


def _report(b: Builder) -> dict | None:
    path = b.claim.files.get("report")
    if not path:
        return None
    try:
        with b.timed("report_html"):
            rep = wfsc.report(path, b.claim.alias)
    except OSError as e:
        b.input("wfsc.report_html", path, state="error", note=str(e))
        return None
    rid = b.input("wfsc.report_html", path)
    b._report_id = rid
    if rep.get("version"):
        b.identity["pipeline_version"] = rep["version"]
        b.sources["pipeline_version"] = {"input": rid, "kind": "html_table"}
    else:
        b.note("info", "version_unknown", "The workflow version is written only in the footer of its HTML report, and it was not found there.")
    if "params" not in b.claim.files and rep.get("parameters"):
        par = rep["parameters"]
        rows = [[k, sr.short_path(par[k]) if k == "ref_genome_dir" else str(par[k])] for k in wfsc.PARAMETER_KEYS
                if par.get(k) not in (None, "None", "")]
        if rows:
            b.table("workflow_parameters", [("key", "str"), ("value", "str")], rows)
            b.sources["workflow_parameters"] = {"input": rid, "kind": "html_table"}
        if par.get("ref_genome_dir"):
            b.identity["reference"] = os.path.basename(str(par["ref_genome_dir"]).rstrip("/")).replace("refdata-gex-", "") or None
        if str(par.get("expected_cells") or "").isdigit():
            b.put("expected_cells", int(par["expected_cells"]), kind="html_table", input_id=rid)
    for key, name in (("reads_report_mapped", "mapped"), ("reads_report_unmapped", "unmapped"),
                      ("reads_report_supplementary", "supplementary")):
        if rep.get(name) is not None:
            b.put(key, rep[name], kind="html_table", input_id=rid)
    if rep.get("mapped") is not None and rep.get("unmapped") is not None and rep["mapped"] + rep["unmapped"] > 0:
        b.put("reads_report_unmapped_fraction", rep["unmapped"] / (rep["mapped"] + rep["unmapped"]), kind="html_table", input_id=rid)
    if rep.get("saturation_reported"):
        b.put("saturation_reported", rep["saturation_reported"]["value"], kind="embedded_json", input_id=rid)
        b.put("saturation_reported_reads", rep["saturation_reported"]["reads"], kind="embedded_json", input_id=rid)
    return rep


def _config_stats(b: Builder, rep: dict | None) -> None:
    path = b.claim.files.get("config_stats")
    rid = getattr(b, "_report_id", None)
    if path:
        try:
            cs = wfsc.read_json(path)
        except (InputError, OSError) as e:
            b.input("wfsc.config_stats", path, state="error", note=getattr(e, "message", None) or str(e))
            cs = None
        if cs is not None:
            b.input("wfsc.config_stats", path, fingerprint=True)
            b.scalar("reads_input", cs.get("reads"))
            b.put("reads_input", cs.get("reads"))
            full = cs.get("full_length", cs.get("full_len"))
            b.scalar("reads_full_length", full)
            b.put("reads_full_length", full)
            cats = [c for c in wfsc.ADAPTER_CLASSES if isinstance(cs.get(c), (int, float))]
            if cats:
                b.counts["reads_adapter_config"] = {"cats": cats, "n": [int(cs[c]) for c in cats], "role": "categories"}
            plus, minus = cs.get("plus"), cs.get("minus")
            if isinstance(plus, (int, float)) and isinstance(minus, (int, float)) and plus + minus > 0:
                b.scalar("reads_strand_plus", plus)
                b.scalar("reads_strand_minus", minus)
                b.put("reads_plus_fraction", plus / (plus + minus))
            return
    if not rep:
        return
    # workflow 3.x without its adapter statistics file: the report is the only place that has these numbers
    if rep.get("input_reads") is not None:
        b.scalar("reads_input", rep["input_reads"])
        b.put("reads_input", rep["input_reads"], kind="html_table", input_id=rid)
    if rep.get("full_length") is not None:
        b.scalar("reads_full_length", rep["full_length"])
        b.put("reads_full_length", rep["full_length"], kind="html_table", input_id=rid)
    pct = rep.get("adapter_pct")
    if pct and rep.get("input_reads"):
        cats = [c for c in wfsc.ADAPTER_CLASSES if c in pct]
        b.counts["reads_adapter_config"] = {"cats": cats, "n": [int(round(pct[c] / 100.0 * rep["input_reads"])) for c in cats],
                                            "role": "categories",
                                            "sources": {c: "rebuilt from the percentage shown in the workflow's report" for c in cats}}
        b.sources["reads_adapter_config"] = {"input": rid, "kind": "embedded_json"}
        b.params["adapter_config_from_percentages"] = True


def _matrices(b: Builder, mito_prefix: str) -> dict | None:
    f = b.claim.files
    role = next((r for r in ("gene_matrix_8um", "gene_matrix", "gene_matrix_2um") if r in f), None)
    mx = None
    if role:
        try:
            with b.timed("gene_matrix"):
                mx = wfsc.mtx_per_barcode(f[role], mito_prefix=mito_prefix)
        except (FlowCellarError, OSError, ValueError, EOFError) as e:
            b.input("wfsc.gene_matrix", os.path.join(f[role], "matrix.mtx.gz"), state="error", note=getattr(e, "message", None) or str(e))
            b.note("warning", "input_unreadable", f"The gene count matrix could not be read: {getattr(e, 'message', None) or e}")
            mx = None
    if mx is not None:
        b.input("wfsc.gene_matrix", os.path.join(f[role], "matrix.mtx.gz"), records=mx["entries"],
                state="ok" if mx["entries"] == mx["entries_stated"] else "truncated")
        b.params["gene_matrix"] = role
        has = mx["umis"] > 0
        umis, genes = mx["umis"][has], mx["genes"][has]
        b.hists["barcodes_umis_hist"] = axes.hist(axes.LOG_COUNTS, umis)
        b.hists["barcodes_genes_hist"] = axes.hist(axes.LOG_COUNTS, genes)
        b.stats("barcodes_umis", umis)
        b.stats("barcodes_genes", genes)
        b.scalar("barcodes_cells", mx["n_barcodes"])
        b.put("barcodes_cells", mx["n_barcodes"])
        b.scalar("barcodes_genes_detected", mx["n_features"])
        b.put("barcodes_genes_detected", mx["n_features"])
        b.scalar("barcodes_unique_molecules", mx["total_umis"])
        b.put("barcodes_unique_molecules", mx["total_umis"])
        b.scalar("total_umis", mx["total_umis"])
        b.params["mito_features"] = mx["n_mito_features"]
        if mx["total_umis"] > 0 and mx["n_mito_features"]:
            b.scalar("mito_umis", int(mx["mito_umis"].sum()))
            b.put("mito_umi_fraction", float(mx["mito_umis"].sum() / mx["total_umis"]))
    trole = next((r for r in ("transcript_matrix_8um", "transcript_matrix", "transcript_matrix_2um") if r in f), None)
    if trole:
        try:
            n_tx, _, _ = wfsc.mtx_size_line(f[trole])
            b.input("wfsc.transcript_matrix", os.path.join(f[trole], "matrix.mtx.gz"), note="the size line only")
            b.scalar("barcodes_transcripts_detected", n_tx)
            b.put("barcodes_transcripts_detected", n_tx)
        except (FlowCellarError, OSError, ValueError, EOFError):
            pass
    return mx


def _mito(b: Builder, mx: dict | None) -> None:
    path = b.claim.files.get("mito")
    frac = None
    if path:
        try:
            frac = wfsc.mito_fractions(path)
            b.input("wfsc.mito_table", path, records=len(frac))
        except (FlowCellarError, OSError, ValueError) as e:
            b.input("wfsc.mito_table", path, state="error", note=getattr(e, "message", None) or str(e))
    if frac is None and mx is not None and mx["n_mito_features"]:
        has = mx["umis"] > 0
        frac = mx["mito_umis"][has] / mx["umis"][has]
        b.params["mito_fraction_source"] = "gene matrix"
    if frac is not None and len(frac):
        b.hists["mito_fraction_hist"] = axes.hist(axes.FRACTION_50, np.clip(frac, 0, 1))
        b.put("mito_fraction_median", float(np.median(frac)))
        b.put("mito_fraction_p90", float(np.percentile(frac, 90)))


def _whitelist(b: Builder, mx: dict | None) -> None:
    path = b.claim.files.get("whitelist")
    if not path:
        return
    try:
        n = wfsc.count_lines(path)
    except OSError:
        return
    b.input("wfsc.whitelist", path, records=n)
    b.scalar("barcodes_cells", n)
    b.put("barcodes_cells", n)
    if mx is not None:
        b.params["cross_checks"]["matrix_barcodes"] = mx["n_barcodes"]


def _trace(b: Builder) -> None:
    path = b.claim.files.get("trace")
    if not path:
        return
    try:
        t = wfsc.trace(path)
    except (OSError, ValueError, KeyError) as e:
        b.input("wfsc.trace", path, state="error", note=str(e))
        return
    b.input("wfsc.trace", path, records=t["tasks"])
    b.put("workflow_tasks", t["tasks"])
    b.put("workflow_failed_tasks", t["failed"])
    b.put("workflow_cached_tasks", t["cached"])
    b.put("workflow_task_hours", round(t["task_hours"], 3))
    b.put("workflow_wall_hours", round(t["wall_hours"], 3) if t["wall_hours"] is not None else None)
    b.put("workflow_peak_rss_bytes", t["peak_rss_bytes"])
    rows = [[r["process"], r["tasks"], round(r["hours"], 3), r["peak_rss_bytes"]] for r in t["top"]]
    if rows:
        b.table("workflow_execution", [("process", "str"), ("tasks", "int"), ("hours", "float"), ("peak_rss_bytes", "int")], rows)


def _bam(b: Builder) -> None:
    f = b.claim.files
    if "bam" not in f:
        return
    if "bam_index" in f:
        try:
            with b.timed("bam_index"):
                m = sr.bam_mapping(f["bam"], f["bam_index"])
            b.input("wfsc.tagged_bam", f["bam"], note="the header and the index only; at most 100,000 records for the link to runs")
            b.scalar("reads_bam_mapped", m["mapped"])
            if m["chrm_fraction"] is not None:
                b.scalar("reads_chrm_mapped", m["chrm_mapped"])
                b.put("reads_chrm_fraction", m["chrm_fraction"])
        except (InputError, OSError, ValueError, EOFError) as e:
            b.input("wfsc.tagged_bam", f["bam"], state="error", note=getattr(e, "message", None) or str(e))
    if b.opts.link_tags and "bam_index" in f:
        with b.timed("link_tags"):
            ev = link.bam_evidence(f["bam"], index=f["bam_index"])
        b.links = link.intrinsic_links(ev)
        b.params["link_evidence"] = {"records": ev.records, "records_without_tag": ev.records_without_tag, "note": ev.note,
                                     "runs": [{k: r[k] for k in ("flow_cell_id", "protocol_run_id_prefix", "acquisition_run_id_prefix",
                                                                 "reads", "st_min", "st_max")} for r in ev.fn_runs]}
        if ev.note:
            b.note("info", "link_tags_not_read", ev.note)


def _read_pass(b: Builder, rep: dict | None, assay: str | None, params: dict) -> dict | None:
    path = b.claim.files.get("read_summary")
    rid = getattr(b, "_report_id", None)

    def from_report(reason_code: str, reason: str) -> None:
        b.note("info", reason_code, reason)
        if not rep:
            return
        for key, name in (("reads_valid_barcode", "valid_barcode"), ("reads_gene_assigned", "gene_assigned"),
                          ("reads_transcript_assigned", "transcript_assigned"), ("barcodes_detected", "barcodes_with_data")):
            if rep.get(name) is not None:
                b.scalar(key, rep[name])
                b.put(key, rep[name], kind="html_table", input_id=rid)

    if not path:
        from_report("read_pass_skipped", "The per-read table of the workflow is not in the folder: read counts per stage come from "
                                         "the workflow's report where it has them.")
        return None
    if not b.opts.read_pass:
        b.input("wfsc.read_summary", path, state="skipped")
        from_report("read_pass_skipped", "The per-read table was not read (--no-read-pass): read counts per stage come from the "
                                         "workflow's report where it has them; reads per unit and the saturation curve are not available.")
        return None
    size = os.path.getsize(path)
    b.say(f"reading the per-read table ({size / 1e9:.1f} GB)")
    t0 = time.monotonic()
    try:
        rp = readpass.read_pass(path, workers=b.opts.workers, dedup=b.opts.dedup, deadline=b.opts.deadline,
                                progress=(lambda fr: b.say(f"per-read table {fr * 100:.0f} %")) if b.opts.progress else None)
    except (FlowCellarError, OSError, ValueError) as e:
        b.input("wfsc.read_summary", path, state="error", note=getattr(e, "message", None) or str(e))
        from_report("read_pass_skipped", f"The per-read table could not be read ({getattr(e, 'message', None) or e}): read counts per "
                                         "stage come from the workflow's report where it has them.")
        return None
    b.timings["read_pass"] = round(time.monotonic() - t0, 2)
    b.params["read_pass"] = {"complete": rp["complete"], "bytes": rp["bytes"], "file_bytes": rp["file_bytes"], "seconds": rp["seconds"],
                             "workers": rp["workers"], "rows": rp["rows"], "molecule_reads": rp["molecule_reads"],
                             "molecules": rp["molecules"]}
    if not rp["complete"]:
        b.input("wfsc.read_summary", path, state="partially_read", records=rp["rows"])
        from_report("read_pass_partial", f"The per-read table was read to {rp['bytes'] / max(rp['file_bytes'], 1) * 100:.0f} % before the time "
                                         "limit: read counts per stage come from the workflow's report where it has them. Scan again to finish.")
        return None
    b.input("wfsc.read_summary", path, records=rp["rows"])
    for key, name in (("reads_valid_barcode", "valid_barcode"), ("reads_valid_umi", "valid_umi"), ("reads_gene_assigned", "gene_assigned"),
                      ("reads_transcript_assigned", "transcript_assigned")):
        b.scalar(key, rp[name])
        b.put(key, rp[name])
    if rp["valid_barcode"]:
        b.put("reads_valid_umi_fraction", rp["valid_umi"] / rp["valid_barcode"])
    per = readpass.reads_per_unit(rp, unit_um=8)
    b.scalar("barcodes_detected", len(per))
    b.put("barcodes_detected", int(len(per)))
    if len(per):
        b.hists["barcodes_reads_hist"] = axes.hist(axes.LOG_COUNTS, per)
        b.stats("barcodes_reads", per)
        rk, rv = axes.rank_curve(per)
        b.curves["barcodes_knee"] = {"x": rk, "y": {"all": rv}, "x_unit": "rank", "y_unit": "reads", "x_scale": "log10"}
    if rp["molecule_hist"] is not None and rp["molecules"]:
        b.hists["saturation_reads_per_molecule"] = axes.reads_per_molecule_hist(rp["molecule_hist"])
        b.scalar("saturation_molecule_reads", rp["molecule_reads"])
        b.put("saturation_molecule_reads", rp["molecule_reads"])
        b.params["cross_checks"]["molecules_from_read_table"] = rp["molecules"]
    b.params["read_id_sketch"] = rp["sketch"]
    # droplet data: the read count above which a barcode is a cell, by the rule of the workflow's documentation
    exp = params.get("expected_cells") or b.summary.get("expected_cells")
    if assay == "singlecell" and isinstance(exp, (int, float)) and exp > 0 and len(per):
        top = np.sort(per)[::-1][: int(exp)]
        b.put("barcodes_cell_threshold", float(np.percentile(top, 95) / 20.0))
    if rep:
        b.params["cross_checks"].update({k: rep.get(k) for k in ("valid_barcode", "gene_assigned", "transcript_assigned", "barcodes_with_data")
                                         if rep.get(k) is not None})
    return rp


def _derived(b: Builder, assay: str | None) -> None:
    s = b.summary
    inp, full = s.get("reads_input"), s.get("reads_full_length")
    kinds = {k: b.sources.get(k) for k in ("reads_input", "reads_full_length", "reads_valid_barcode", "reads_gene_assigned",
                                           "reads_transcript_assigned")}

    def derived(key: str, value, *parents: str) -> None:
        src = next((kinds[p] for p in parents if kinds.get(p)), None)
        b.put(key, value, kind=src["kind"] if src else "data", input_id=src["input"] if src else None)
    if inp and full is not None:
        derived("reads_full_length_fraction", min(1.0, full / inp), "reads_input", "reads_full_length")
    if full:
        for key, part in (("reads_valid_barcode_fraction", "reads_valid_barcode"), ("reads_gene_assigned_fraction", "reads_gene_assigned"),
                          ("reads_transcript_assigned_fraction", "reads_transcript_assigned")):
            if s.get(part) is not None:
                derived(key, s[part] / full, "reads_full_length", part)
    mol, gene = s.get("barcodes_unique_molecules"), s.get("reads_gene_assigned")
    if mol is not None and gene:
        derived("saturation_value", 1.0 - mol / gene, "reads_gene_assigned")
    units = s.get("barcodes_cells") if assay == "singlecell" else (s.get("barcodes_detected") or s.get("barcodes_cells"))
    if inp and units:
        derived("barcodes_reads_per_cell_mean", inp / units, "reads_input")
    cats = [c for c in FUNNEL if s.get("reads_" + c) is not None]
    missing = [c for c in FUNNEL if s.get("reads_" + c) is None]
    n = [int(s.get("reads_" + c) or 0) for c in FUNNEL]
    if cats:
        stages, counts = list(FUNNEL), n
        if mol is not None:
            stages, counts = stages + ["unique_molecules"], counts + [int(mol)]
        fun: dict = {"cats": stages, "n": counts, "role": "funnel"}
        if missing:
            fun["missing"] = missing
        words = {"data": "data file", "html_table": "workflow report", "embedded_json": "workflow report"}
        fun["sources"] = {c: words[(b.sources.get("reads_" + c) or {}).get("kind", "data")] for c in cats}
        b.counts["reads_funnel"] = fun


def _map(b: Builder, mx: dict | None, rp: dict | None, assay: str | None) -> None:
    if mx is None or not mx["barcodes"]:
        return
    try:
        if b.claim.hd and mx["barcodes"][0].startswith("s_"):
            r8, c8, um = sr.hd_rows_cols(mx["barcodes"])
            reads = readpass.hd_grid(rp, unit_um=um) if rp is not None else None
            b.grids["spatial_map"] = maps.bins_map(r8, c8, umis=mx["umis"], genes=mx["genes"], reads=reads, bin_um=um, map_um=b.opts.map_um)
            b.params["map"] = {"kind": "visium_hd_workflow", "bin_um": b.grids["spatial_map"]["layout"]["bin_um"], "tissue_image": False}
        elif assay == "visium":
            _visium_map_from_sibling(b, mx, rp)
    except (FlowCellarError, OSError, ValueError, KeyError) as e:
        b.note("info", "map_not_made", f"The spatial map could not be made: {getattr(e, 'message', None) or e}")


def _visium_map_from_sibling(b: Builder, mx: dict, rp: dict | None) -> None:
    """Standard Visium: the workflow does not know where a spot lies. The position list of the Space Ranger
    output of the same reads does; without it there is no map."""
    outs = b.opts.sibling or next(iter(detect.spaceranger_outs_below(b.claim.path) or detect.spaceranger_outs_below(b.claim.top or b.claim.path)), None)
    pos_path = None
    if outs:
        c = detect.spaceranger_claim(outs)
        pos_path = c.files.get("tissue_positions") if c else None
    if not pos_path:
        b.note("info", "map_needs_positions", "No spatial map: the workflow's output does not say where a spot lies on the slide. "
                                              "It is drawn when the Space Ranger output of the same reads is found (or named with --sibling).")
        return
    pos = sr.tissue_positions(pos_path)
    b.input("spaceranger.tissue_positions", pos_path, records=len(pos), note="of the Space Ranger output of the same reads: spot positions for the map")
    names = list(pos)
    stem = {n.split("-")[0]: i for i, n in enumerate(names)}
    idx = np.array([stem.get(n.split("-")[0], -1) for n in mx["barcodes"]], dtype=np.int64)
    ok = idx >= 0
    if not ok.any():
        return

    def spread(values) -> np.ndarray:
        out = np.zeros(len(names), dtype=np.float64)
        out[idx[ok]] = np.asarray(values)[ok]
        return out
    umis, genes, mito = spread(mx["umis"]), spread(mx["genes"]), spread(mx["mito_umis"])
    reads = None
    if rp is not None and rp["unit_names"]:
        reads = np.zeros(len(names), dtype=np.float64)
        for name, n in rp["unit_names"].items():
            i = stem.get(name.decode("ascii", "replace").split("-")[0], -1)
            if i >= 0:
                reads[i] = n
    b.grids["spatial_map"] = maps.visium_map(
        np.array([pos[n][1] for n in names]), np.array([pos[n][2] for n in names]), np.array([pos[n][0] for n in names]),
        umis=umis, genes=genes, mito=np.divide(mito, umis, out=np.zeros_like(umis), where=umis > 0), reads=reads)
    b.params["map"] = {"kind": "visium_workflow", "spots": int(len(names)), "positions_from": "sibling"}


def _sibling_hint(b: Builder, params: dict) -> None:
    """Where the Space Ranger output of the same reads should be, as far as the files of this analysis say."""
    hint: dict = {}
    sbam = params.get("spaceranger_bam")
    if isinstance(sbam, str) and sbam:
        hint["spaceranger_bam"] = sbam
        out_dir = str(params.get("out_dir") or "")
        top = (b.claim.top or "").rstrip("/")
        norm = os.path.normpath(out_dir) if out_dir else ""
        launch = None
        if norm and norm not in (".", "") and top.endswith(norm.lstrip("./")):
            launch = top[: len(top) - len(norm.lstrip("./"))].rstrip("/") or "/"
        if launch and not os.path.isabs(sbam):
            cand = os.path.dirname(os.path.normpath(os.path.join(launch, sbam)))
            hint["expected_outs"] = cand
            if detect.is_spaceranger_outs(cand):
                hint["outs"] = cand
        elif os.path.isabs(sbam) and detect.is_spaceranger_outs(os.path.dirname(sbam)):
            hint["outs"] = os.path.dirname(sbam)
        hint["tail"] = "/".join(os.path.dirname(os.path.normpath(sbam)).split("/")[-3:])
    below = detect.spaceranger_outs_below(b.claim.path) + [p for p in detect.spaceranger_outs_below(b.claim.top or b.claim.path)]
    if below:
        hint["outs_below"] = sorted(set(below))
    if b.opts.sibling:
        hint["stated"] = os.path.abspath(b.opts.sibling)
    if hint:
        b.params["sibling_hint"] = hint
