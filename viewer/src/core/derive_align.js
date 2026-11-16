// Pure functions from a bundle to what the alignment, over-time and transcript charts draw (sections aln, aln.time,
// aln.windows, dyn.basecalled, rna.polya, rna.tx, rna.coverage). Aggregates are summed, never re-binned.
//
// Per sample and added up: up to twelve samples every cell carries its own grids, depth windows and curves; with
// more, they are stored once for the passed reads of all samples under run_level, and a chart of one sample says
// that it cannot be drawn. Quartiles cannot be added: for one sample they are read from its stored curve, for all
// samples they are worked out here from the table under run_level (the same rule: by rank inside a bin).
import { mod, readCells, axisEdges, expand, addInto, quantile, sum } from './derive.js';

const pick = (b, id, keys, filters) => { const want = new Set(keys), f = new Set(filters); return readCells(b, id).filter(c => want.has(c.sample) && f.has(c.filter)).map(c => c.cell); };
// the selection is "everything the run holds" when no sample is chosen: then what is stored once for the run may be used
export const wholeRun = c => !c.sKeys.length;

const ALN = { reads: 'reads', mapped: 'mapped_reads', unmapped: 'unmapped_reads', secondary: 'secondary_records', supplementary: 'supplementary_records', withSupp: 'reads_with_supplementary',
  fwd: 'forward_reads', rev: 'reverse_reads', mm: 'match_or_mismatch_bases', ins: 'inserted_bases', del: 'deleted_bases', skip: 'skipped_bases', soft: 'soft_clipped_bases',
  hard: 'hard_clipped_bases', mism: 'mismatched_bases', mmNm: 'match_or_mismatch_bases_with_nm', columns: 'alignment_columns', edits: 'edit_distance', perfect: 'perfect_reads' };
// counts of the chosen cells of an alignment section, added up; cigar: the files were read with the CIGAR totals
export function alnTotals(b, keys, filters, id = 'aln') {
  const out = { cells: 0, cigar: false }; for (const k of Object.keys(ALN)) out[k] = 0;
  for (const c of pick(b, id, keys, filters)) { const sc = c.scalars || {}; out.cells++; if (sc.inserted_bases != null) out.cigar = true; for (const [k, f] of Object.entries(ALN)) out[k] += sc[f] || 0; }
  return out;
}
// a 1-D histogram of the chosen cells, added up: { edges, n, w, log, missing, over }
export function hist1(b, id, name, keys, filters) {
  const cells = pick(b, id, keys, filters).map(c => c.hists && c.hists[name]).filter(Boolean); if (!cells.length) return null;
  const ax = b.axes[cells[0].axis]; if (!ax) return null;
  const out = { edges: axisEdges(ax), n: new Float64Array(ax.n), w: null, log: ax.type === 'log10', missing: 0, over: 0, under: 0, axis: cells[0].axis };
  for (const h of cells) { if (h.axis !== cells[0].axis) continue; addInto(out.n, expand(h, ax.n)); if (h.w) { if (!out.w) out.w = new Float64Array(ax.n); addInto(out.w, expand(h, ax.n, 'w')); } out.missing += h.missing || 0; out.over += h.over || 0; out.under += h.under || 0; }
  return out;
}
function addGrid(g, h) { g.outside += h.outside || 0; g.missing += h.missing || 0; for (let i = 0; i < h.i.length; i++) if (h.i[i] < g.nx && h.j[i] < g.ny) g.grid[h.i[i]][h.j[i]] += h.n[i]; }
function newGrid(b, h, nx) { const xa = b.axes[h.x], ya = b.axes[h.y]; if (!ya || (!xa && nx == null)) return null; const n = nx != null ? nx : xa.n;
  return { xe: nx != null ? null : axisEdges(xa), ye: axisEdges(ya), xlog: xa && xa.type === 'log10', ylog: ya.type === 'log10', nx: n, ny: ya.n, grid: Array.from({ length: n }, () => new Float64Array(ya.n)), outside: 0, missing: 0 }; }
// a 2-D histogram of the chosen cells, added up; when the cells do not carry it (more than twelve samples) and the whole
// run is asked for, the one stored for the run. { xe, ye, grid[ix][iy], nx, ny, outside, missing, run } or null;
// { unavailable: true } when it exists for the run only and a sample was chosen
export function grid2(b, id, name, keys, filters, whole) {
  const cells = pick(b, id, keys, filters).map(c => c.hists2d && c.hists2d[name]).filter(Boolean), m = mod(b, id);
  const rl = m && m.run_level && m.run_level.hists2d && m.run_level.hists2d[name];
  if (!cells.length) { if (!rl) return null; if (!whole) return { unavailable: true }; const g = newGrid(b, rl); if (!g) return null; addGrid(g, rl); g.run = true; return g; }
  const g = newGrid(b, cells[0]); if (!g) return null; for (const h of cells) addGrid(g, h); return g;
}
// median and quartiles of the y values in every x bin of a grid: [{ i, n, q25, med, q75 }] for bins with at least `min` values
// (50, or 10 when the whole grid holds fewer than 20,000 reads)
export function gridQuantiles(g, min) {
  if (min == null) { let t = 0; for (let i = 0; i < g.nx; i++) t += sum(g.grid[i]); min = t >= 20000 ? 50 : 10; }
  const out = []; for (let i = 0; i < g.nx; i++) { const n = sum(g.grid[i]); if (n < min) continue; out.push({ i, n, q25: quantile(g.ye, g.grid[i], 0.25, g.ylog), med: quantile(g.ye, g.grid[i], 0.5, g.ylog), q75: quantile(g.ye, g.grid[i], 0.75, g.ylog) }); }
  return out;
}

// ---- along the genome ----
export const genomeRef = b => { const m = mod(b, 'aln.windows'); return (m && m.params && m.params.reference) || null; };
// reads and aligned bases per contig of the chosen cells: { rows: [{ name, reads, bases, length }], reads, bases, mito }
export function contigs(b, keys, filters) {
  const tabs = pick(b, 'aln', keys, filters).map(c => c.counts && c.counts.aln_contig_table).filter(Boolean); if (!tabs.length) return null;
  const ref = genomeRef(b) || {}, len = {}; for (const c of ref.contigs || []) len[c.name] = c.length;
  const by = new Map(); for (const t of tabs) t.cats.forEach((name, i) => { const o = by.get(name) || { name, reads: 0, bases: 0, length: name === 'other contigs' ? ((ref.other || {}).length || null) : (len[name] || null) }; o.reads += t.n[i]; o.bases += t.w ? t.w[i] : 0; by.set(name, o); });
  const rows = [...by.values()], other = rows.find(r => r.name === 'other contigs');
  const amr = (mod(b, 'aln').params || {}).reference || {};
  return { rows: rows.filter(r => r !== other).concat(other && (other.reads || other.length) ? [other] : []), reads: rows.reduce((a, r) => a + r.reads, 0), bases: rows.reduce((a, r) => a + r.bases, 0), mito: amr.mitochondrial || ref.mitochondrial || null, total: ref.total_length || amr.total_length || null };
}
// aligned bases and depth per window of the chosen cells: { ref, bases, depth, len, run } (run: the one stored for the run)
export function depthWindows(b, keys, filters, whole) {
  const m = mod(b, 'aln.windows'), ref = genomeRef(b); if (!m || !ref || !ref.windows) return null;
  const cells = pick(b, 'aln.windows', keys, filters).map(c => c.series && c.series.aln_depth_windows).filter(Boolean), rl = m.run_level && m.run_level.series && m.run_level.series.aln_depth_windows;
  if (!cells.length && !(rl && whole)) return rl ? { unavailable: true } : null;
  const n = ref.windows, bases = new Float64Array(n); for (const s of (cells.length ? cells : [rl])) { const v = s.values.aligned_bases; for (let i = 0; i < n && i < v.length; i++) bases[i] += v[i]; }
  const len = new Float64Array(n), depth = new Float64Array(n), W = ref.window_bases;
  for (const c of ref.contigs) for (let k = 0; k < c.windows; k++) { const i = c.first_window + k; len[i] = k === c.windows - 1 ? c.length - k * W : W; depth[i] = len[i] > 0 ? bases[i] / len[i] : 0; }
  return { ref, bases, depth, len, run: !cells.length };
}

// the lower edge of the bin that holds the rank q: the value itself where every value of a bin is its lower edge
export function stepValue(edges, counts, q) { let total = 0; for (const c of counts) total += c; if (!(total > 0)) return null; let acc = 0; const want = q * total; for (let i = 0; i < counts.length; i++) { const c = counts[i]; if (acc + c >= want && c > 0) return edges[i]; acc += c; } return edges[edges.length - 1]; }
// ---- over run time ----
export const clockOf = (b, id) => { const m = mod(b, id); return (m && m.params && m.params.clock) || null; };
// the middle of hour bin k, in hours since the run's time origin
const hourMid = (clock, k) => (clock.offset_3600s_s + (k + 0.5) * 3600) / 3600;
// Median and quartiles of one quantity per hour, for the selection.
//   one sample chosen: its stored curve; several: one line per sample; none: all samples, from the table of the run.
// -> { lines: [{ key, t: [hours], n, q25, med, q75 }], run: bool, unit, unavailable }
export function overTime(b, id, name, c, min = 30) {
  const m = mod(b, id), clock = clockOf(b, id); if (!m || !clock) return null;
  const fromCurve = (key, cv) => { const t = [], n = [], q25 = [], med = [], q75 = []; cv.x.forEach((k, i) => { if (cv.y.reads[i] < min) return; t.push(hourMid(clock, k)); n.push(cv.y.reads[i]); q25.push(cv.y.q25[i]); med.push(cv.y.median[i]); q75.push(cv.y.q75[i]); }); return { key, t, n, q25, med, q75 }; };
  const filt = c.filters.includes('pass') ? 'pass' : c.filters[0];
  if (c.sKeys.length) {
    const lines = []; for (const key of c.sKeys) { const cell = readCells(b, id).find(q => q.sample === key && q.filter === filt); const cv = cell && cell.cell.curves && cell.cell.curves[name]; if (cv) lines.push(fromCurve(key, cv)); }
    if (!lines.length) return { lines: [], unavailable: !!(m.run_level && m.run_level.hists2d), unit: null };
    return { lines, run: false, filter: filt };
  }
  const hs = (m.run_level && m.run_level.hists2d) || {}, parts = c.filters.map(f => hs[name + '_' + f]).filter(Boolean); if (!parts.length) return null;
  const ya = b.axes[parts[0].y]; if (!ya) return null; const ye = axisEdges(ya), log = ya.type === 'log10';
  let nx = 0; for (const h of parts) for (const i of h.i) if (i + 1 > nx) nx = i + 1;
  const grid = Array.from({ length: nx }, () => new Float64Array(ya.n)); for (const h of parts) for (let k = 0; k < h.i.length; k++) if (h.j[k] < ya.n) grid[h.i[k]][h.j[k]] += h.n[k];
  const line = { key: null, t: [], n: [], q25: [], med: [], q75: [] };
  // read qualities that are whole numbers (params.qual_whole_numbers): every read of a 1 Q step has that value, so a
  // quartile is the value of the step that holds its rank, as the engine works it out for one sample
  const whole = name === 'qual_over_time' && m.params && m.params.qual_whole_numbers, q = (counts, p) => whole ? stepValue(ye, counts, p) : quantile(ye, counts, p, log);
  for (let i = 0; i < nx; i++) { const n = sum(grid[i]); if (n < min) continue; line.t.push(hourMid(clock, i)); line.n.push(n); line.q25.push(q(grid[i], 0.25)); line.med.push(q(grid[i], 0.5)); line.q75.push(q(grid[i], 0.75)); }
  return { lines: [line], run: true };
}
// Reads and bases by start time, ten minutes by ten minutes, of all samples: { t (hours, end of each step), th (middle),
// pass: { reads, bases }, fail, hasFail }; cumulative sums are taken by the caller
export function outputFromReads(b) {
  const m = mod(b, 'dyn.basecalled'), clock = clockOf(b, 'dyn.basecalled'), s = m && m.run_level && m.run_level.series && m.run_level.series.dyn_output_series_600s; if (!s || !clock) return null;
  const v = s.values, n = Math.max(...Object.values(v).map(a => a.length)), z = k => { const a = new Float64Array(n), src = v[k] || []; for (let i = 0; i < src.length; i++) a[i] = src[i]; return a; };
  const t = [], th = []; for (let i = 0; i < n; i++) { t.push((clock.offset_600s_s + (i + 1) * 600) / 3600); th.push((clock.offset_600s_s + (i + 0.5) * 600) / 3600); }
  return { t, th, pass: { reads: z('reads_pass'), bases: z('bases_pass') }, fail: { reads: z('reads_fail'), bases: z('bases_fail') }, hasFail: !!v.reads_fail, hasPass: !!v.reads_pass };
}
// per sample, hour by hour (the cells of dyn.basecalled): { t: [hours, end of each hour], reads, bases } added over the chosen cells
export function outputOfSamples(b, keys, filters) {
  const clock = clockOf(b, 'dyn.basecalled'), ss = pick(b, 'dyn.basecalled', keys, filters).map(c => c.series && c.series.dyn_output_series).filter(Boolean); if (!ss.length || !clock) return null;
  const n = Math.max(...ss.map(s => s.values.reads.length)), reads = new Float64Array(n), bases = new Float64Array(n);
  for (const s of ss) { s.values.reads.forEach((x, i) => { reads[i] += x; }); s.values.bases.forEach((x, i) => { bases[i] += x; }); }
  const t = []; for (let i = 0; i < n; i++) t.push(Math.max(0, (clock.offset_3600s_s + (i + 1) * 3600) / 3600));
  return { t, reads, bases };
}
export function speedTotals(b, keys, filters) {
  const o = { reads: 0, timed: 0, withDur: 0, bases: 0, ms: 0 }; for (const c of pick(b, 'dyn.basecalled', keys, filters)) { const sc = c.scalars || {}; o.reads += sc.reads || 0; o.timed += sc.reads_with_start_time || 0; o.withDur += sc.reads_with_duration || 0; o.bases += sc.bases_with_duration || 0; o.ms += sc.duration_ms || 0; }
  return o;
}

// ---- poly(A) and transcripts ----
export function polya(b, keys, filters) {
  const cells = pick(b, 'rna.polya', keys, filters); if (!cells.length) return null;
  const o = { reads: 0, tagged: 0, found: 0, zero: 0, none: 0, hist: hist1(b, 'rna.polya', 'rna_polya_hist', keys, filters) };
  for (const c of cells) { const sc = c.scalars || {}; o.reads += sc.reads || 0; o.tagged += sc.reads_with_tag || 0; o.found += sc.tail_found_reads || 0; o.zero += sc.tail_zero_reads || 0; o.none += sc.tail_not_found_reads || 0; }
  return o;
}
export const TX_CLASSES = [['under_1kb', 'under 1 kb'], ['1_to_2kb', '1 to 2 kb'], ['2_to_4kb', '2 to 4 kb'], ['4kb_and_longer', '4 kb and longer']];
// coverage along the transcript body of the chosen samples: [{ id, label, reads, cover: Float64Array(100) share of reads }]
export function bodyCoverage(b, keys) {
  const m = mod(b, 'rna.coverage'); if (!m) return null; const want = new Set(keys), si = m.dims.indexOf('sample'), ci = m.dims.indexOf('read_class');
  return TX_CLASSES.map(([id, label]) => { const n = new Float64Array(100); let reads = 0;
    for (const c of m.cells) { if (!want.has(c.at[si]) || c.at[ci] !== id) continue; reads += (c.scalars || {}).reads || 0; const h = c.hists && c.hists.rna_body_coverage; if (h) addInto(n, expand(h, 100)); }
    return { id, label, reads, cover: Array.from(n, v => reads ? v / reads : null) }; });
}
// the curves of the alignment to the transcripts for the selection: one sample -> its own; none chosen -> the run's
export function txCurves(b, c, name) {
  const m = mod(b, 'rna.tx'); if (!m) return null; const filt = c.filters.includes('pass') ? 'pass' : c.filters[0];
  if (c.sKeys.length) { const out = []; for (const key of c.sKeys) { const cell = readCells(b, 'rna.tx').find(q => q.sample === key && q.filter === filt), cv = cell && cell.cell.curves && cell.cell.curves[name]; if (cv) out.push({ key, cv }); } return out; }
  const cv = m.run_level && m.run_level.curves && m.run_level.curves[name]; return cv ? [{ key: null, cv }] : [];
}
