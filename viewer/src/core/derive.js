// Pure functions from a bundle to what the panels draw. Aggregates are summed, never re-binned.
const memo = new WeakMap();
function once(b, key, fn) { let m = memo.get(b); if (!m) { m = new Map(); memo.set(b, m); } if (!m.has(key)) m.set(key, fn()); return m.get(key); }
const fin = v => typeof v === 'number' && isFinite(v);

export const mod = (b, id) => (b && b.modules && b.modules[id]) || null;
export const cell0 = (b, id) => { const m = mod(b, id); return m && m.cells && m.cells[0] || null; };
export const summary = (b, id) => { const m = mod(b, id); return (m && m.summary) || {}; };

// ---- axes ----
export function axisEdges(ax) {
  const n = ax.n, e = new Float64Array(n + 1);
  if (ax.type === 'log10') for (let i = 0; i <= n; i++) e[i] = Math.pow(10, ax.min + i / ax.bins_per_decade);
  else for (let i = 0; i <= n; i++) e[i] = ax.start + i * ax.step;
  return e;
}
export const edgesOf = (b, name) => once(b, 'edges:' + name, () => axisEdges(b.axes[name]));
// a stored histogram (start + values) -> a full-length array on its axis
export function expand(h, n, key = 'n') { const out = new Float64Array(n), v = h && h[key]; if (!v) return out; const s = h.start || 0; for (let i = 0; i < v.length; i++) if (s + i < n) out[s + i] = v[i]; return out; }
export function addInto(a, b) { for (let i = 0; i < a.length; i++) a[i] += b[i]; return a; }
// quantile by rank interpolation inside a bin (the specification's rule)
export function quantile(edges, counts, q, log) {
  let total = 0; for (const c of counts) total += c; if (!(total > 0)) return null;
  let acc = 0; const want = q * total;
  for (let i = 0; i < counts.length; i++) { const c = counts[i]; if (acc + c >= want && c > 0) { const f = (want - acc) / c; return log ? Math.pow(10, Math.log10(edges[i]) + f * (Math.log10(edges[i + 1]) - Math.log10(edges[i]))) : edges[i] + f * (edges[i + 1] - edges[i]); } acc += c; }
  return edges[edges.length - 1];
}
// N50: the length at which half of the bases are in longer reads; w = bases per bin
export function n50(edges, w, log = true) {
  let total = 0; for (const c of w) total += c; if (!(total > 0)) return null;
  let acc = 0; const want = total / 2;
  for (let i = w.length - 1; i >= 0; i--) { const c = w[i]; if (acc + c >= want && c > 0) { const f = (want - acc) / c; return log ? Math.pow(10, Math.log10(edges[i + 1]) - f * (Math.log10(edges[i + 1]) - Math.log10(edges[i]))) : edges[i + 1] - f * (edges[i + 1] - edges[i]); } acc += c; }
  return edges[0];
}
export const sum = a => { let t = 0; for (const v of a) t += v || 0; return t; };
export const cumsum = a => { const o = new Float64Array(a.length); let t = 0; for (let i = 0; i < a.length; i++) { t += a[i] || 0; o[i] = t; } return o; };

// ---- run health ----
export const flags = b => { const t = b.tiers || {}, rh = t.run_health || {}, rl = t.read_level || {}, dyn = b.modules && b.modules.dyn; return {
  stub: !b.modules || !b.modules.run || rh.state === 'absent',
  quick: rh.detail === 'quick',
  full: rh.detail === 'full' && !!dyn,
  // "so far": checks still wait for the read count. A run whose read table is cut off, or could not be read, is in
  // quick detail for good: it is not being counted, and its verdict is not "so far".
  soFar: !!(b.status && b.status.provisional),
  partlyRead: rh.detail === 'quick' && !!dyn && dyn.state === 'partial' && (dyn.cells || []).length > 0,
  cut: rh.detail === 'quick' && !!dyn && !(b.status && (b.status.pending || []).length),
  unread: !!dyn && dyn.state === 'error',
  t2: !!rl.state && rl.state !== 'absent',
  t2partial: rl.state === 'partial',
  running: ((b.identity || {}).protocol_run || {}).run_state === 'running',
}; };
export const runHours = b => summary(b, 'run').run_duration_h || 0;

// output: reads and estimated bases per 600 s bin
export function output(b) {
  return once(b, 'output', () => {
    const c = cell0(b, 'yield'); if (!c || !c.series || !c.series.yield_series) return null;
    const v = c.series.yield_series.values, n = v.reads.length, step = 600 / 3600;
    const tEnd = new Float64Array(n + 1), cr = new Float64Array(n + 1), cb = new Float64Array(n + 1);
    for (let i = 0; i < n; i++) { tEnd[i + 1] = (i + 1) * step; cr[i + 1] = cr[i] + v.reads[i]; cb[i + 1] = cb[i] + v.estimated_bases[i]; }
    // per hour: sums of six bins
    const nh = Math.ceil(n / 6), th = new Float64Array(nh), rh = new Float64Array(nh), bh = new Float64Array(nh);
    for (let i = 0; i < n; i++) { const k = (i / 6) | 0; rh[k] += v.reads[i]; bh[k] += v.estimated_bases[i]; }
    for (let k = 0; k < nh; k++) { th[k] = k + 0.5; const bins = Math.min(6, n - k * 6); if (bins < 6) { rh[k] *= 6 / bins; bh[k] *= 6 / bins; } }
    return { n, t: tEnd, reads: cr, bases: cb, th, readsPerH: rh, basesPerH: bh, raw: v };
  });
}
export function valueAt(t, cum, h) { if (!t || h > t[t.length - 1] + 1e-6) return null; const i = Math.min(t.length - 1, Math.round(h * 6)); return cum[i]; }

export const STATE_GROUPS = [
  ['sequencing', 'Sequencing', ['strand', 'adapter'], '--fc-stategroup-sequencing'],
  ['pore', 'Pore free', ['pore'], '--fc-stategroup-pore'],
  ['recovering', 'Recovering', ['unavailable', 'unblocking'], '--fc-stategroup-recovering'],
  ['inactive', 'No pore', ['no_pore', 'zero', 'disabled', 'unknown_positive', 'unknown_negative', 'multiple', 'saturated'], '--fc-stategroup-inactive'],
  ['unclassified', 'Unclassified', ['locked', 'membrane', 'unclassified', 'unclassified_following_reset', 'pending_manual_reset', 'pending_mux_change'], '--fc-stategroup-unclassified'],
];
// channel activity: channel-seconds per 600 s bin, in MinKNOW's five groups, plus occupancy and working channels
export function activity(b) {
  return once(b, 'activity', () => {
    const c = cell0(b, 'fc'); if (!c || !c.series || !c.series.fc_channel_state_series) return null;
    const v = c.series.fc_channel_state_series.values, n = v.strand.length, t = new Float64Array(n), total = new Float64Array(n), g = {};
    for (const [id, , members] of STATE_GROUPS) { const a = g[id] = new Float64Array(n); for (const m of members) if (v[m]) for (let i = 0; i < n; i++) a[i] += v[m][i]; }
    const occ = new Array(n), working = new Float64Array(n), share = {};
    for (let i = 0; i < n; i++) { t[i] = (i + 0.5) / 6; for (const k in g) total[i] += g[k][i]; const busy = v.strand[i] + v.adapter[i], w = busy + v.pore[i]; occ[i] = w > 0 ? busy / w : null; working[i] = w / 600; }
    for (const k in g) { const a = share[k] = new Float64Array(n); for (let i = 0; i < n; i++) a[i] = total[i] > 0 ? g[k][i] / total[i] : 0; }
    const seqSec = new Float64Array(n), workSec = new Float64Array(n); for (let i = 0; i < n; i++) { seqSec[i] = v.strand[i] + v.adapter[i]; workSec[i] = seqSec[i] + v.pore[i]; }
    return { n, t, share, occ, working, total, seqSec, workSec, raw: v };
  });
}
// Channel activity for drawing. A pore scan takes every channel away from sequencing for a few minutes, about every
// hour and a half, and a 10-minute step shows that as a notch in every band. A step counts as a scan step when it lies
// beside a recorded scan AND the states MinKNOW reports only while it is assessing channels (zero, saturated, multiple,
// unknown, locked, unclassified, pending) rise clearly above the steps around it; a fall that has another cause (pores
// blocking before a scan frees them) is left as it is. Scan steps take the shares of the steps beside them and the
// scans are ticked above the plot. Nothing is lost: "measured" keeps the stored shares, which the read-out and the
// table print, and one control draws them.
const SCAN_STATES = ['zero', 'saturated', 'multiple', 'unknown_positive', 'unknown_negative', 'locked', 'unclassified', 'unclassified_following_reset', 'pending_manual_reset', 'pending_mux_change'];
export function activityEven(b) {
  return once(b, 'activityEven', () => {
    const a = activity(b), sc = scans(b); if (!a) return null;
    const n = a.n, near = new Uint8Array(n), hit = new Uint8Array(n), busy = new Float64Array(n);
    for (const k of SCAN_STATES) { const v = a.raw[k]; if (v) for (let i = 0; i < n; i++) busy[i] += v[i]; }
    for (let i = 0; i < n; i++) busy[i] = a.total[i] > 0 ? busy[i] / a.total[i] : 0;
    if (sc) for (const t of sc.t) { const k = Math.floor(t * 6 + 1e-9); for (const q of [k - 1, k, k + 1]) if (q >= 0 && q < n) near[q] = 1; }
    for (let i = 0; i < n; i++) {
      if (!near[i]) continue;
      let l = i - 1, r = i + 1; while (l >= 0 && near[l]) l--; while (r < n && near[r]) r++;
      const ref = Math.min(l >= 0 ? busy[l] : Infinity, r < n ? busy[r] : Infinity);
      if (isFinite(ref) && busy[i] > ref + 0.02) hit[i] = 1;
    }
    const share = {};
    for (const key of Object.keys(a.share)) {
      const src = a.share[key], out = share[key] = Float64Array.from(src);
      for (let i = 0; i < n; i++) {
        if (!hit[i]) continue;
        let l = i - 1, r = i + 1; while (l >= 0 && hit[l]) l--; while (r < n && hit[r]) r++;
        out[i] = l >= 0 && r < n ? src[l] + (src[r] - src[l]) * (i - l) / (r - l) : l >= 0 ? src[l] : r < n ? src[r] : src[i];
      }
    }
    let scanSteps = 0; for (let i = 0; i < n; i++) scanSteps += hit[i];
    return { ...a, share, measured: a.share, scanStep: hit, scanSteps };
  });
}
export const SCAN_CLASSES = [
  ['single_pore', 'Available', '--fc-scan-single-pore'], ['reserved_pore', 'In reserve', '--fc-scan-reserved-pore'], ['unavailable', 'Unavailable', '--fc-scan-unavailable'],
  ['multiple', 'Multiple', '--fc-scan-multiple'], ['saturated', 'Saturated', '--fc-scan-saturated'], ['zero', 'Zero', '--fc-scan-zero'], ['other', 'Other', '--fc-scan-other'],
];
export function scans(b) {
  return once(b, 'scans', () => {
    const c = cell0(b, 'fc'); if (!c || !c.series || !c.series.fc_pore_scan_series) return null;
    const v = c.series.fc_pore_scan_series.values, n = v.single_pore.length, t = Array.from(v.t_s, s => s / 3600), pores = new Float64Array(n);
    for (let i = 0; i < n; i++) pores[i] = v.single_pore[i] + (v.reserved_pore ? v.reserved_pore[i] : 0);
    return { n, t, v, pores };
  });
}
export function temperature(b) {
  const c = cell0(b, 'fc'), s = c && c.series && c.series.fc_temperature_series; if (!s || !s.values || !s.values.flowcell_c) return null;
  const step = b.axes[s.axis].step / 3600, v = s.values, n = v.flowcell_c.length, p = (mod(b, 'fc').params || {});
  return { t: Array.from({ length: n }, (_, i) => (s.start || 0) * step + i * step), flowcell: v.flowcell_c, chamber: v.chamber_c || [], lo: v.target_min_c || [], hi: v.target_max_c || [], dropped: !!p.temperature_stop_sample_dropped };
}
export function voltage(b) {
  const c = cell0(b, 'fc'), s = c && c.series && c.series.fc_bias_voltage_series; if (!s) return null;
  return { t: s.values.t_s.map(x => x / 3600), mv: s.values.voltage_mv };
}
export function messages(b) {
  const c = cell0(b, 'run'), tb = c && c.tables && c.tables.run_messages; if (!tb) return [];
  const t0 = Date.parse(b.identity.protocol_run.time_origin || b.identity.protocol_run.started);
  return tb.rows.map(r => ({ time: r[0], h: (Date.parse(r[0]) - t0) / 3.6e6, severity: r[1], id: r[2], text: r[3] }));
}

export const END_REASONS = [
  ['signal_positive', 'Ended normally', '--fc-end-signal-positive'], ['unblock_mux_change', 'Ejected: strand stalled', '--fc-end-unblock-mux-change'],
  ['data_service_unblock_mux_change', 'Ejected by adaptive sampling', '--fc-end-data-service-unblock-mux-change'], ['signal_negative', 'Ended with a current drop', '--fc-end-signal-negative'],
  ['mux_change', 'Cut by a pore scan', '--fc-end-mux-change'], ['analysis_config_change', 'Settings changed', '--fc-end-analysis-config-change'],
  ['device_data_error', 'Device data error', '--fc-end-device-data-error'], ['api_request', 'Stopped by software', '--fc-end-api-request'], ['paused', 'Paused', '--fc-end-paused'], ['unknown', 'Not recorded', '--fc-end-unknown'],
];
export function endReasons(b) {
  return once(b, 'end', () => {
    const m = mod(b, 'dyn'); if (!m) return null;
    const by = {}; for (const c of m.cells) by[c.at[0]] = c;
    const rows = []; let reads = 0, basesT = 0;
    for (const [id, label, token] of END_REASONS) { const c = by[id]; if (!c) continue; rows.push({ id, label, token, reads: c.scalars.reads, bases: c.scalars.estimated_bases, ms: c.scalars.duration_ms, cell: c }); reads += c.scalars.reads; basesT += c.scalars.estimated_bases; }
    for (const id of Object.keys(by)) if (!END_REASONS.some(e => e[0] === id)) { const c = by[id]; rows.push({ id, label: id, token: '--fc-end-other', reads: c.scalars.reads, bases: c.scalars.estimated_bases, ms: c.scalars.duration_ms, cell: c }); reads += c.scalars.reads; basesT += c.scalars.estimated_bases; }
    // per hour of read start
    const ser = r => (r.cell.series && r.cell.series.dyn_end_reason_series && r.cell.series.dyn_end_reason_series.values) || { reads: [], estimated_bases: [] };
    const nh = Math.max(1, ...rows.map(r => ser(r).reads.length));
    const t = Array.from({ length: nh }, (_, i) => i + 0.5), total = new Float64Array(nh), totalB = new Float64Array(nh);
    for (const r of rows) { const v = ser(r); r.perH = new Float64Array(nh); r.perHB = new Float64Array(nh); for (let i = 0; i < v.reads.length; i++) { r.perH[i] = v.reads[i]; r.perHB[i] = v.estimated_bases[i]; total[i] += v.reads[i]; totalB[i] += v.estimated_bases[i]; } }
    return { rows, reads, bases: basesT, t, total, totalB, nh };
  });
}
// speed of reads that ended normally: median per hour of read start, and the whole-run histogram
export function speed(b) {
  return once(b, 'speed', () => {
    const m = mod(b, 'dyn'); if (!m) return null;
    const c = m.cells.find(x => x.at[0] === 'signal_positive'); if (!c || !c.hists2d || !c.hists2d.dyn_speed_est_series) return null;
    const h = c.hists2d.dyn_speed_est_series, ax = b.axes[h.y], ny = ax.n, nx = Math.max(...h.i) + 1, grid = Array.from({ length: nx }, () => new Float64Array(ny));
    for (let k = 0; k < h.i.length; k++) grid[h.i[k]][h.j[k]] += h.n[k];
    const e = axisEdges(ax), t = [], med = [], q1 = [], q3 = [];
    for (let i = 0; i < nx; i++) { const tot = sum(grid[i]); t.push(i + 0.5); if (tot < 200) { med.push(null); q1.push(null); q3.push(null); continue; } med.push(quantile(e, grid[i], 0.5)); q1.push(quantile(e, grid[i], 0.25)); q3.push(quantile(e, grid[i], 0.75)); }
    const hist = expand(c.hists && c.hists.dyn_speed_est_hist, ny);
    return { t, med, q1, q3, edges: e, hist, grid, nx, ny };
  });
}
const histAxis = (b, cells, name, fallback) => { for (const c of cells) { const h = c.hists && c.hists[name]; if (h && h.axis && b.axes[h.axis]) return h.axis; } return fallback; };
export function durationHist(b) {
  return once(b, 'dur', () => { const m = mod(b, 'dyn'); if (!m) return null; const name = histAxis(b, m.cells, 'dyn_duration_hist', 'duration_log40'), ax = b.axes[name]; if (!ax) return null;
    const n = new Float64Array(ax.n); for (const c of m.cells) if (c.hists && c.hists.dyn_duration_hist) addInto(n, expand(c.hists.dyn_duration_hist, ax.n)); return { edges: axisEdges(ax), n }; });
}
export function lenEst(b) {
  return once(b, 'lenest', () => {
    const m = mod(b, 'len.est'); if (!m) return null; const name = histAxis(b, m.cells, 'len_est_hist', 'read_length_log40'), ax = b.axes[name]; if (!ax) return null;
    const n = new Float64Array(ax.n), w = new Float64Array(ax.n); let outside = 0;
    for (const c of m.cells) { const h = c.hists && c.hists.len_est_hist; if (!h) continue; addInto(n, expand(h, ax.n, 'n')); addInto(w, expand(h, ax.n, 'w')); outside += (h.under || 0) + (h.over || 0); }
    return { edges: axisEdges(ax), n, w, outside };
  });
}
// the x domain that holds all but the outermost reads (default 0.05 % at each end), so a handful of extreme reads does not flatten the picture
export function spanQ(edges, counts, tail = 0.0005, pad = 3) { let tot = 0; for (const c of counts) tot += c; if (!(tot > 0)) return [edges[0], edges[edges.length - 1]]; let a = 0, acc = 0; while (a < counts.length - 1 && acc + counts[a] < tot * tail) acc += counts[a++]; let b = counts.length - 1; acc = 0; while (b > a && acc + counts[b] < tot * tail) acc += counts[b--]; a = Math.max(0, a - pad); b = Math.min(counts.length - 1, b + pad); return [edges[a], edges[b + 1]]; }
// trim leading and trailing empty bins: [lo, hi] of the x domain that holds data
export function span(edges, counts, pad = 2) { let a = 0, b = counts.length - 1; while (a < b && !(counts[a] > 0)) a++; while (b > a && !(counts[b] > 0)) b--; a = Math.max(0, a - pad); b = Math.min(counts.length - 1, b + pad); return [edges[a], edges[b + 1]]; }

// ---- channel grid (u6: 6 bits per value, base64url alphabet, square-root scale) ----
const CODE = new Int8Array(128).fill(-1); 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-_'.split('').forEach((c, i) => { CODE[c.charCodeAt(0)] = i; });
// MIME base64 uses + and /; accept both alphabets
CODE['+'.charCodeAt(0)] = 62; CODE['/'.charCodeAt(0)] = 63;
export function decodeU6(str) { const out = new Uint8Array(str.length); for (let i = 0; i < str.length; i++) { const v = CODE[str.charCodeAt(i)]; out[i] = v < 0 ? 0 : v; } return out; }
// sparkline points (0..1) of an index row; with the run's hours, the points after the run ended are left out
// (the index pads the 32 points of the fixed 0 to 96 h axis with zeros)
export function sparkPoints(str, hours) { const out = [], n = hours != null && hours === hours ? Math.max(1, Math.min(32, Math.ceil(hours / 3 - 1e-9))) : 32;
  for (let i = 0; i < (str || '').length && i < n; i++) { const v = CODE[str.charCodeAt(i)]; if (v < 0) break; out.push(v / 63); } return out; }

// ---- read data: sections with the dimensions sample and filter (pass, fail); numbers are exact unless the section says partial ----
export function samplesOf(b) { return (b.samples || []).filter(s => s.key !== 'unclassified'); }
// cells of a section as [{ sample, filter, cell }], whatever the order of its dimensions
export function readCells(b, id) {
  return once(b, 'cells:' + id, () => {
    const m = mod(b, id); if (!m) return [];
    const si = (m.dims || []).indexOf('sample'), fi = (m.dims || []).indexOf('filter');
    return (m.cells || []).map(c => ({ sample: si >= 0 ? c.at[si] : null, filter: fi >= 0 ? c.at[fi] : 'pass', cell: c }));
  });
}
// which read sets the files hold: { pass, fail }
export function filtersOf(b) { const out = { pass: false, fail: false }; for (const id of ['len', 'yield.basecalled', 'qual']) for (const c of readCells(b, id)) if (c.filter in out) out[c.filter] = true; return out; }
// the state of a read-data section: complete is the normal case and says nothing
export function sectionState(b, id) {
  const m = mod(b, id); if (!m) return null;
  const partial = m.state === 'partial', cov = m.coverage || null;
  return { partial, reason: partial ? (m.partial_reason || 'partly_read') : null, coverage: cov, population: m.population || null, error: m.state === 'error' };
}
// sample keys that have a cell in the section
export function sampleKeysWithData(b, id = 'len') { const s = new Set(); for (const c of readCells(b, id)) if (c.sample && c.sample !== 'unclassified') s.add(c.sample); return s; }
// length, quality, identity and mapping-quality histograms of the chosen samples and read sets, added up
export function readHists(b, keys, filters = ['pass']) {
  const want = new Set(keys), fset = new Set(filters), pick = id => readCells(b, id).filter(c => want.has(c.sample) && fset.has(c.filter)).map(c => c.cell);
  const L = pick('len'), Q = pick('qual'), A = pick('aln');
  const laxName = histAxis(b, readCells(b, 'len').map(c => c.cell), 'len_hist', 'read_length_log40'), lax = b.axes[laxName], qax = b.axes[histAxis(b, readCells(b, 'qual').map(c => c.cell), 'qual_hist', 'qscore_0p25')];
  const iax = b.axes[histAxis(b, readCells(b, 'aln').map(c => c.cell), 'aln_identity_hist', 'identity_q_0p1')], max = b.axes[histAxis(b, readCells(b, 'aln').map(c => c.cell), 'aln_mapq_hist', 'mapq_1')];
  const out = { len: lax ? { edges: axisEdges(lax), n: new Float64Array(lax.n), w: new Float64Array(lax.n) } : null, qual: qax ? { edges: axisEdges(qax), n: new Float64Array(qax.n), missing: 0 } : null, reads: 0, bases: 0,
    ident: iax ? { edges: axisEdges(iax), n: new Float64Array(iax.n), over: 0 } : null, mapq: max ? { edges: axisEdges(max), n: new Float64Array(max.n) } : null,
    mapped: 0, unmapped: 0, alnReads: 0, perfect: 0, secondary: 0, supplementary: 0, columns: 0, edits: 0 };
  for (const l of L) { const h = l.hists && l.hists.len_hist; if (h && out.len) { addInto(out.len.n, expand(h, lax.n, 'n')); addInto(out.len.w, expand(h, lax.n, 'w')); } const sc = l.scalars || {}; out.reads += sc.reads || 0; out.bases += sc.bases || 0; }
  for (const q of Q) { const h = q.hists && q.hists.qual_hist; if (h && out.qual) { addInto(out.qual.n, expand(h, qax.n)); out.qual.missing += h.missing || 0; } }
  for (const a of A) {
    const hs = a.hists || {}, sc = a.scalars || {};
    if (hs.aln_identity_hist && out.ident) { addInto(out.ident.n, expand(hs.aln_identity_hist, iax.n)); out.ident.over += hs.aln_identity_hist.over || 0; }
    if (hs.aln_mapq_hist && out.mapq) addInto(out.mapq.n, expand(hs.aln_mapq_hist, max.n));
    out.mapped += sc.mapped_reads || 0; out.unmapped += sc.unmapped_reads || 0; out.alnReads += sc.reads || 0; out.perfect += sc.perfect_reads || 0;
    out.secondary += sc.secondary_records || 0; out.supplementary += sc.supplementary_records || 0; out.columns += sc.alignment_columns || 0; out.edits += sc.edit_distance || 0;
  }
  return out;
}
// reads and bases per sample and read set, from yield.basecalled (else from len): { key: { pass: {reads, bases}, fail: {...} } }
export function sampleYield(b) {
  return once(b, 'sampleYield', () => {
    const out = {}, src = readCells(b, 'yield.basecalled').length ? readCells(b, 'yield.basecalled') : readCells(b, 'len');
    for (const c of src) { const o = out[c.sample] || (out[c.sample] = {}), sc = c.cell.scalars || {}; o[c.filter] = { reads: sc.reads || 0, bases: sc.bases || 0 }; }
    return out;
  });
}
// length x quality grid of the chosen samples and read sets, added up: { xe, ye, grid[ix][iy] reads, w[ix][iy] bases }
export function lenQualGrid(b, keys, filters = ['pass']) {
  const want = new Set(keys), fset = new Set(filters), cells = readCells(b, 'qual').filter(c => want.has(c.sample) && fset.has(c.filter)).map(c => c.cell.hists2d && c.cell.hists2d.qual_vs_length).filter(Boolean);
  if (!cells.length) return null;
  const xa = b.axes[cells[0].x], ya = b.axes[cells[0].y]; if (!xa || !ya) return null;
  const grid = Array.from({ length: xa.n }, () => new Float64Array(ya.n)), w = Array.from({ length: xa.n }, () => new Float64Array(ya.n)); let hasW = true, outside = 0;
  for (const h of cells) { if (!h.w) hasW = false; outside += h.outside || 0; for (let i = 0; i < h.i.length; i++) if (h.i[i] < xa.n && h.j[i] < ya.n) { grid[h.i[i]][h.j[i]] += h.n[i]; if (h.w) w[h.i[i]][h.j[i]] += h.w[i]; } }
  return { xe: axisEdges(xa), ye: axisEdges(ya), grid, w: hasW ? w : null, nx: xa.n, ny: ya.n, outside };
}
export const FUNNEL_LABELS = { sequenced: 'Sequenced on the flow cell', basecalled: 'Basecalled', pass: 'Passed the quality filter', classified: 'Given to a sample', full_length: 'Full-length', genome_mapped: 'Mapped to the genome', transcriptome_mapped: 'Mapped to transcripts', quantified: 'Counted towards a transcript' };
export function funnel(b) { const c = cell0(b, 'funnel'); const t = c && c.counts && c.counts.funnel_stage_table; if (!t) return null; return t.cats.map((id, i) => ({ id, label: FUNNEL_LABELS[id] || id, n: t.n[i], missing: (t.missing || []).includes(id), source: (t.sources || {})[id] })); }
export function funnelSamples(b) { const m = mod(b, 'funnel.samples'); if (!m) return {}; const by = {}; for (const c of m.cells) { const t = c.counts && c.counts.funnel_stage_table, o = {}; if (!t) continue; t.cats.forEach((id, i) => { o[id] = (t.missing || []).includes(id) ? null : t.n[i]; }); by[c.at[0]] = o; } return by; }
export function cdna(b) { const m = mod(b, 'cdna'); if (!m) return {}; const by = {}; for (const c of m.cells) by[c.at[0]] = c.scalars; return by; }

// ---- comparison: output = working channel-hours x sequencing share x estimated bases per sequencing second ----
export function factors(b, hoursEnd) {
  const a = activity(b), o = output(b); if (!a || !o) return null;
  const nb = Math.min(a.n, o.n, Math.round(hoursEnd * 6)); let work = 0, seq = 0;
  for (let i = 0; i < nb; i++) { work += a.workSec[i]; seq += a.seqSec[i]; }
  const basesOut = o.bases[nb];
  return { hours: nb / 6, bases: basesOut, reads: o.reads[nb], workHours: work / 3600, seqShare: work > 0 ? seq / work : null, basesPerSeqSec: seq > 0 ? basesOut / seq : null };
}
export { fin };
