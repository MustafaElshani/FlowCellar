// Single-cell and spatial analyses: plain functions (no markup), so that they can be checked on their own
// (verify/analysis_units.mjs). An analysis is a bundle of its own kind: one pipeline output for one sample, with
// the flow cell runs it came from. Everything here reads the section "sc" of such a bundle or a row of the index.
import cat from '../gen/sc_catalogue.js';
import * as F from '../core/format.js';

export const ANALYSIS_ID = /^an_[a-z]{2,16}_[0-9a-f]{12}$/;
export const isAnalysisId = id => typeof id === 'string' && ANALYSIS_ID.test(id);
export const ASSAY_WORD = { visium: 'Visium', visium_hd: 'Visium HD', singlecell: 'Single cell' };
export const ASSAY_BADGE = { visium: 'VIS', visium_hd: 'HD', singlecell: 'SC' };
export const PIPELINE_WORD = { spaceranger: 'Space Ranger', 'wf-single-cell': 'wf-single-cell' };
export const PIPELINE_TAB = { spaceranger: 'Space Ranger', 'wf-single-cell': 'ONT workflow' };
// how a link from an analysis to a run was made, in words
export const LINK_WORD = { read_tag: 'matched by the file names stored on the reads', bam_header_rg: 'matched by the read group in the BAM header',
  sibling_analysis: 'taken from the ONT workflow analysis of the same reads', read_id_sketch: 'matched by read names', path_rule: 'the run folder lies beside the analysis and carries the same sample name',
  sample_name_chain: 'matched by sample name only', manual: 'stated by a person' };
export const hrefAnalysis = id => '#/analysis/' + encodeURIComponent(id);
// what a file of an analysis is, in words (the last part of its role)
export const FILE_WORD = { bam_index: 'BAM index', cmdline: 'command line', config_stats: 'adapter statistics', feature_slice: 'feature slice', gene_matrix: 'gene matrix', matrix: 'count matrix',
  metrics_summary: 'metrics summary', mito_table: 'mitochondrial share per unit', molecule_info: 'molecule information', params: 'workflow settings', pipestance_metrics: 'the pipeline’s own metrics',
  read_summary: 'per-read table', report_html: 'report', tagged_bam: 'tagged BAM', tissue_positions: 'tissue positions', trace: 'task trace', transcript_matrix: 'transcript matrix', web_summary: 'web summary',
  whitelist: 'barcode list', scalefactors: 'scale factors', tissue_image: 'tissue image', barcodes: 'barcodes', features: 'features', filtered_matrix: 'filtered matrix', raw_matrix: 'raw matrix' };
export const fileWord = role => { const k = String(role || '').split('.').pop(); return FILE_WORD[k] || k.replace(/_/g, ' '); };
// a number that has no catalogue entry of its own (it feeds a check or a table), in words
const VALUE_WORD = { barcodes_reads_per_cell_mean: 'mean reads per cell', expected_cells: 'expected cells', pipeline_version: 'pipeline version', reads_adapter_config: 'adapter configuration',
  reads_full_length: 'full-length reads', reads_full_length_fraction: 'full-length share', reads_gene_assigned_fraction: 'share assigned to a gene', reads_input: 'input reads',
  reads_report_mapped: 'mapped reads', reads_report_supplementary: 'supplementary alignments', reads_report_unmapped: 'unmapped reads', reads_report_unmapped_fraction: 'unmapped share',
  reads_transcript_assigned_fraction: 'share resolved to one isoform', reads_valid_barcode_fraction: 'valid barcode rate', reads_valid_barcode: 'reads matched to a barcode',
  reads_gene_assigned: 'reads assigned to a gene', reads_transcript_assigned: 'reads resolved to one isoform', saturation_reported: 'saturation as the pipeline states it',
  saturation_reported_reads: 'reads behind the stated saturation', saturation_value: 'saturation', spatial_alarms_error: 'pipeline errors', spatial_alarms_warn: 'pipeline warnings',
  spatial_registration: 'image alignment', spatial_registration_flag: 'image alignment check', spatial_sample: 'sample description', sr_corrected_bc_frac: 'share of corrected barcodes',
  sr_low_support_umi_reads_frac: 'share of reads with a weakly supported UMI', sr_tso_frac: 'share of reads with the template-switch oligo', workflow_parameters: 'workflow settings',
  barcodes_cells: 'cells called', barcodes_detected: 'barcodes with data', barcodes_unique_molecules: 'unique molecules' };
export const valueWord = id => (cat.metrics[id] || {}).label || VALUE_WORD[id] || String(id).replace(/^(sr|sc|sp)_/, '').replace(/_frac(tion)?$/, ' share').replace(/_/g, ' ');

export const metricLabel = id => (cat.metrics[id] || {}).label || id;
// what a barcode stands for: { name: 'spot', plural: 'spots' }
export function unitWords(assay) { return (cat.units || {})[assay] || { name: 'unit', plural: 'units' }; }
export const cap = s => s ? s.charAt(0).toUpperCase() + s.slice(1) : s;

// ---- the section of a bundle ----
export const section = b => (b && b.modules && b.modules.sc) || null;
export const cellOf = b => { const s = section(b); return (s && s.cells && s.cells[0]) || {}; };
export const summaryOf = b => { const s = section(b); return (s && s.summary) || {}; };
export const analysisOf = b => (b && b.identity && b.identity.analysis) || {};
// true when a value was taken from a report's HTML and not from a data file
export function fromReport(b, key) { const s = section(b), v = s && s.value_sources && s.value_sources[key]; return !!v && (v.kind === 'html_table' || v.kind === 'embedded_json'); }
export const alarmsOf = b => ((section(b) || {}).notes || []).filter(n => n.code === 'pipeline_alarm');
export const notesOf = (b, codes) => ((section(b) || {}).notes || []).filter(n => codes.includes(n.code));
export const tableOf = (b, name) => { const t = (cellOf(b).tables || {})[name]; return t ? { columns: t.columns.map(c => c.id), rows: t.rows } : null; };
export const tableMap = (b, name) => { const t = tableOf(b, name), o = {}; if (t) for (const r of t.rows) o[r[0]] = r[1]; return o; };

// ---- axes of the add-on: 8 bins per decade from 1, 50 bins of 0.02, reads per molecule ----
export function edgesOfAxis(ax) {
  if (!ax) return null;
  if (ax.type === 'edges') return Float64Array.from(ax.edges);
  const n = ax.n, e = new Float64Array(n + 1);
  if (ax.type === 'log10') for (let i = 0; i <= n; i++) e[i] = Math.pow(10, ax.min + i / ax.bins_per_decade);
  else for (let i = 0; i <= n; i++) e[i] = ax.start + i * ax.step;
  return e;
}
// a stored histogram on its axis: { edges, n, w, under, over }
export function histOf(b, name) {
  const h = (cellOf(b).hists || {})[name]; if (!h || !b.axes || !b.axes[h.axis]) return null;
  const edges = edgesOfAxis(b.axes[h.axis]), bins = edges.length - 1, n = new Float64Array(bins), w = h.w ? new Float64Array(bins) : null;
  (h.n || []).forEach((v, i) => { if (h.start + i < bins) n[h.start + i] = v; });
  if (w) h.w.forEach((v, i) => { if (h.start + i < bins) w[h.start + i] = v; });
  return { edges, n, w, under: h.under || 0, over: h.over || 0, min: h.min, max: h.max, axis: h.axis };
}

// ---- saturation: the exact expectation under random subsampling of the reads ----
// h: a reads-per-molecule histogram (n molecules and w reads per bin). p: the share of the reads that is kept.
export function moleculesAt(h, p) {
  let unique = 0;
  for (let i = 0; i < h.n.length; i++) { const n = h.n[i]; if (!(n > 0)) continue; const c = h.w[i] / n; unique += n * (p >= 1 ? 1 : 1 - Math.pow(1 - p, c)); }
  return unique;
}
export function saturationAt(h, p) { let reads = 0; for (const v of h.w) reads += v; if (!(reads > 0) || !(p > 0)) return null; return 1 - moleculesAt(h, p) / (p * reads); }
// the curve from a tenth of the reads to twice the reads: { reads[], saturation[], molecules[], total, at1, gain2 }
export function saturationCurve(h, steps = 40) {
  let total = 0, mol = 0; for (let i = 0; i < h.n.length; i++) { total += h.w[i]; mol += h.n[i]; }
  if (!(total > 0)) return null;
  const reads = [], sat = [], molecules = [];
  for (let k = 1; k <= steps; k++) { const p = k / steps; reads.push(p * total); sat.push(saturationAt(h, p)); molecules.push(moleculesAt(h, p)); }
  // twice the reads cannot be read off a subsample: the molecules it would add are predicted from how many were seen
  // once, twice, three times ... (the alternating sum of Good and Toulmin, over the counts that are stored one by one)
  let gain = 0; for (let c = 1; c <= Math.min(64, h.n.length); c++) gain += (c % 2 ? 1 : -1) * (h.n[c - 1] || 0);
  return { reads, saturation: sat, molecules, total, molecules1: mol, at1: 1 - mol / total, gain2: mol > 0 ? Math.max(0, gain) / mol : null };
}

// ---- the spatial map: one byte per cell and channel ----
const NODATA = 255;
export function decodeBytes(b64) { const bin = atob(b64), out = new Uint8Array(bin.length); for (let i = 0; i < bin.length; i++) out[i] = bin.charCodeAt(i); return out; }
// the value a level stands for; NaN where the channel says "no data"
export function levelValue(ch, q) {
  if (ch.nodata && q === NODATA) return NaN;
  if (ch.transform === 'log1p') return Math.expm1(q * ch.scale);
  if (ch.transform === 'sqrt') return (q * ch.scale) * (q * ch.scale) + (ch.offset || 0);
  return (ch.offset || 0) + q * ch.scale;
}
// { q: Uint8Array, value(i), lo, hi }: lo and hi are the levels the colour ramp runs between (the 99.5th percentile
// of the cells with data, so that one bright cell does not darken the whole picture)
export function decodeChannel(ch) {
  const q = decodeBytes(ch.data), hist = new Uint32Array(256); let n = 0, top = 0;
  for (let i = 0; i < q.length; i++) { if (ch.nodata && q[i] === NODATA) continue; hist[q[i]]++; n++; if (q[i] > top) top = q[i]; }
  let hi = 0, acc = 0; const want = n * 0.995;
  for (let k = 0; k < 256; k++) { acc += hist[k]; if (acc >= want && hist[k]) { hi = k; break; } }
  if (ch.transform === 'linear' && ch.unit === 'fraction') hi = Math.round((1 - (ch.offset || 0)) / ch.scale);
  if (ch.unit === 'grey') hi = 255;
  // top: the highest level any cell holds; where it lies above hi, the last colour stands for "hi or more"
  return { q, n, hi: Math.max(1, hi), top, value: i => levelValue(ch, q[i]), nodata: !!ch.nodata };
}
export const CHANNEL_WORD = { umis: 'UMIs', genes: 'Genes', reads: 'Reads', mito: 'Mitochondrial share', saturation: 'Saturation', tissue: 'Tissue', image: 'Tissue image' };
export const CHANNEL_ORDER = ['umis', 'genes', 'reads', 'mito', 'saturation', 'tissue', 'image'];
// how one cell's value is written: channel id -> text. A count is read off one byte per cell, so it is a whole number
// that may be a little off: "about" from ten on (the first levels are exact). A cell without a value says of what.
const about = v => (v >= 10 ? 'about ' : '') + F.int(Math.round(v));
export function channelText(id, v, unit) {
  if (!(v === v)) return (id === 'umis' ? 'UMIs' : (CHANNEL_WORD[id] || id).toLowerCase()) + ': no data';
  if (id === 'umis') return about(v) + ' UMIs';
  if (id === 'reads') return about(v) + ' reads';
  if (id === 'genes') return (unit === 'bin' && v < 10 ? F.sig(v) : about(v)) + (Math.round(v) === 1 && !(unit === 'bin' && v < 10) ? ' gene' : ' genes') + (unit === 'bin' ? ' per 8 µm bin' : '');
  if (id === 'mito') return F.pct1(v) + ' mitochondrial';
  if (id === 'saturation') return F.pct(v) + ' saturation';
  if (id === 'tissue') return unit === 'spot' ? (v >= 0.5 ? 'under tissue' : 'off tissue') : F.pct(v) + ' under tissue';
  if (id === 'image') return 'grey level ' + Math.round(v);
  return String(v);
}
// the two ends of a colour scale, in words
export function channelRange(id, ch, dec) {
  const hi = levelValue(ch, dec.hi), more = dec.top > dec.hi ? ' or more' : '';       // the scale ends below the largest cell: its last colour is "or more"
  if (id === 'mito' || id === 'saturation' || id === 'tissue') return { lo: '0', hi: F.pctNum(Math.min(1, hi)), unit: '%' + (hi < 1 ? more : '') };
  if (id === 'image') return { lo: 'dark', hi: 'bright', unit: '' };
  return { lo: '0', hi: hi >= 100 ? F.int(hi) : F.sig(hi), unit: (id === 'umis' ? 'UMIs' : (CHANNEL_WORD[id] || id).toLowerCase()) + more + (ch.transform === 'log1p' ? ', log scale' : '') };
}

// ---- numbers of a row or a summary, written ----
export const pct = v => F.isNum(v) ? F.pct(v) : F.DASH;
export function unitsOf(row) { const v = row.spatial_units_under_tissue != null ? row.spatial_units_under_tissue : row.barcodes_cells != null ? row.barcodes_cells : row.barcodes_detected; return v; }
export function unitsLabel(row) { const u = unitWords(row.assay); return row.spatial_units_under_tissue != null ? u.plural + ' under tissue' : row.assay === 'singlecell' ? 'cells' : u.plural + ' with data'; }
// genes and UMIs per unit as the pipeline states them: a median for spots and cells, a mean for Visium HD bins
export function perUnit(row, what) {
  const sp = row['spatial_' + what + '_per_unit'], med = row['barcodes_' + what + '_median'];
  if (sp != null) return { v: sp, stat: row.unit_stat || 'median' };
  if (med != null) return { v: med, stat: 'median' };
  return null;
}
export function saturationOf(row) {
  if (row.pipeline === 'spaceranger' && row.sr_saturation != null) return { v: row.sr_saturation, from: 'Space Ranger' };
  if (row.saturation_value != null) return { v: row.saturation_value, from: 'ONT workflow' };
  if (row.sr_saturation != null) return { v: row.sr_saturation, from: 'Space Ranger' };
  return null;
}
export const num1 = v => !F.isNum(v) ? F.DASH : v >= 100 ? F.int(v) : (Math.round(v * 10) / 10).toFixed(1);

// ---- the analyses of an index, grouped ----
// rows that show the same sample (the ONT workflow analysis and the Space Ranger output of the same reads) are one
// entry: the row that was asked for first, then its siblings
export function family(rows, byId, id) {
  const first = byId.get(id); if (!first) return [];
  const out = [first];
  for (const s of first.siblings || []) { const r = byId.get(s); if (r && !out.includes(r)) out.push(r); }
  return out;
}
// one entry per sample for the table: [{ key, rows: [main, sibling...] }]; Space Ranger leads a pair (it names the sample)
export function entries(rows) {
  const byId = new Map(rows.map(r => [r.id, r])), seen = new Set(), out = [];
  for (const r of rows) {
    if (seen.has(r.id)) continue;
    const fam = family(rows, byId, r.id).sort((a, b) => (a.pipeline === 'spaceranger' ? 0 : 1) - (b.pipeline === 'spaceranger' ? 0 : 1));
    for (const x of fam) seen.add(x.id);
    out.push({ key: fam[0].id, rows: fam });
  }
  return out;
}
// how many different samples the rows are about: analyses of one slide area (or, without one, of one name) are one
export const sampleCount = rows => new Set(rows.map(r => r.group || r.label || r.id)).size;
// the labels that more than one entry carries: such entries also name their folder, to be told apart
export function repeatedLabels(list) { const n = new Map(); for (const e of list) { const k = e.rows[0].label || ''; n.set(k, (n.get(k) || 0) + 1); } return new Set([...n].filter(x => x[1] > 1).map(x => x[0])); }
const RANK = { fail: 3, warn: 2, pass: 1, unknown: 0 };
export const worstVerdict = rows => rows.reduce((w, r) => (RANK[r.verdict] || 0) > (RANK[w] || 0) ? r.verdict : w, 'unknown');
// the flow cell runs of an entry, each once, with the largest share stated for it
export function linksOf(rows) {
  const by = new Map();
  for (const r of rows) for (const ln of r.links || []) { const k = ln.run || ln.fc; const e = by.get(k); if (!e || (ln.share || 0) > (e.share || 0)) by.set(k, { ...ln, resolved: !!ln.resolved || !!(e && e.resolved) }); }
  return [...by.values()].sort((a, b) => (b.share || 0) - (a.share || 0));
}
// the share of an analysis's reads that came from one run, as far as it is known: counted on the records that were
// probed, so it is written to the whole per cent with "about" ("about 55 %"), never to a decimal
export const aboutShare = f => !F.isNum(f) ? '' : f >= 0.995 ? 'about 100 %' : f < 0.005 ? 'under 1 %' : 'about ' + Math.round(f * 100) + ' %';
