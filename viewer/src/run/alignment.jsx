// Charts drawn from alignment records and from the tags the basecaller writes on every read: the summary of read
// alignment (a ring with its table), what the bases of the mapped reads are, depth along the genome, reads per
// chromosome, lengths, calibration of the quality scores, and under "More charts" strand, soft clipping and the
// 2-D grids; a section of tracks over run time on the time axis of the run timeline (one cursor for all of them);
// a section on the transcripts (poly(A) tail, coverage along the transcript body, saturation, reads per transcript);
// and the chart of reads and bases per sample that stands with the Samples table.
//
// Every chart has a title, one sentence written from its numbers, a read-out, "Show as table", the reads it describes,
// and follows the Samples control. Numbers are exact (every record of the files) unless the section says partial.
import { setParams } from '../core/state.js';
import { Tile } from '../shell/parts.jsx';
import { ChartPanel, Track, Empty } from '../charts/ChartPanel.jsx';
import { genomeTrack } from '../charts/kit.js';
import { tk } from '../charts/theme.js';
import * as F from '../core/format.js';
import * as D from '../core/derive.js';
import * as A from '../core/derive_align.js';
import { metric, panelTitle } from '../core/labels.js';
import { Section, Disclosure, AllNumbers, countNumbers } from './section.jsx';
import { scopeTag, partialInfo, readSetWords, sampleName, chosenWords, ofReads, variantNote } from './common.js';
import { Grid } from '../charts/Grid.jsx';

const ML = 58, MR = 18;                                  // the margins of the run timeline: one hairline through every track
const pid = panel => 'p-' + panel.replace('/', '-');
const how = id => metric(id).how;
const lower = s => s ? s.charAt(0).toLowerCase() + s.slice(1) : s;
const tagOf = (st, p) => [st.tag, p ? p.tag : null].filter(Boolean).join(' · ');
const q2p = q => 100 * (1 - Math.pow(10, -q / 10));      // identity in per cent from its place on the error-rate scale
const pctId = q => (q >= 30 ? q2p(q).toFixed(2) : q2p(q).toFixed(1)) + ' %';
// the same on an axis: no trailing zeros, so the ticks read 90 %, 95 %, 99 %, 99.5 %, 99.9 %, 99.99 %
const tickId = q => +(q >= 30 ? q2p(q).toFixed(2) : q2p(q).toFixed(1)) + ' %';
// a length on a logarithmic axis, as the x axis of the kit writes it: the unit at 1, 10, 100 ..., bare numbers at 2 and 5
const lenTicks = (d0, d1) => { const out = []; for (let e = Math.floor(Math.log10(d0)); e <= Math.ceil(Math.log10(d1)); e++) for (const m of [1, 2, 5]) { const v = m * Math.pow(10, e); if (v >= d0 * 0.999 && v <= d1 * 1.001) out.push(v); } return out; };
const lenTick = v => { const m = v / Math.pow(10, Math.floor(Math.log10(v) + 1e-9)), u = v >= 1e6 ? [1e6, ' Mb'] : v >= 1e3 ? [1e3, ' kb'] : [1, ' b']; return Math.abs(m - 1) < 1e-6 ? +(v / u[0]).toPrecision(3) + u[1] : String(+(v / u[0]).toPrecision(3)); };
// run time as hours and minutes: one way of writing a time on every track
const hm = h => { const m = Math.max(0, Math.round(h * 60)); return Math.floor(m / 60) + ':' + String(m % 60).padStart(2, '0'); };
const stepWords = (a, e) => 'Run time ' + hm(Math.max(0, a)) + ' to ' + hm(e);
const ID_TICKS = [10, 13, 20, 23, 30, 40];                // 90, 95, 99, 99.5, 99.9 and 99.99 % identity
const mids = (e, log) => Array.from({ length: e.length - 1 }, (_, i) => log ? Math.sqrt(e[i] * e[i + 1]) : (e[i] + e[i + 1]) / 2);
const slot = (s, i) => { const m = /(\d+)$/.exec((s && s.barcode) || ''); return ((m ? +m[1] - 1 : i) % 12) + 1; };
const colorOf = (c, key) => { const i = c.named.findIndex(q => q.key === key); return tk('--fc-cat-' + slot(c.named[i], i < 0 ? 11 : i)); };
const brush = r => { if (!r) return setParams({ t: '' }); let a = Math.max(0, Math.round(r[0] * 6) / 6), e = Math.round(r[1] * 6) / 6; if (e - a < 0.5) e = a + 0.5; setParams({ t: +a.toFixed(2) + '-' + +e.toFixed(2) }); };
const share = (v, t) => t > 0 ? v / t : 0;
const small = v => v >= 0.1 ? F.pct1(v) : v > 0 ? (v * 100).toPrecision(2) + ' %' : '0 %';
const depthX = v => (v >= 100 ? Math.round(v) : v >= 10 ? v.toFixed(1) : v >= 1 ? v.toFixed(2) : v > 0 ? v.toPrecision(2) : '0') + '×';
const Mb = v => v >= 1e6 ? +(v / 1e6).toFixed(1) + ' Mb' : v >= 1e3 ? +(v / 1e3).toFixed(1) + ' kb' : v + ' b';
// "only for the run" notes: a chart of one sample that is stored once for all samples
// a section that has numbers, but none for the samples chosen in the Samples control: said, not left out
const NoneFor = ({ c, what }) => <Empty title={'No ' + what + ' for ' + chosenWords(c)}>The files read so far hold no {what} of {c.sKeys.length > 1 ? 'these samples' : 'this sample'}: other sections of the page may. <button type="button" class="fc-linkbtn" onClick={() => setParams({ s: '' })}>Show all samples</button></Empty>;
const RunOnly = ({ title }) => <Empty title={title + ': all samples only'}>This run has more than twelve samples: this chart is kept once, for all samples added up. Choose “all samples” in the Samples control to see it.</Empty>;

// ---------- cells of a grid for the heat form: square root of the count, so that thin regions still show ----------
function heatCells(g) {
  const I = [], J = [], T = []; let mx = 0; for (let i = 0; i < g.nx; i++) for (let j = 0; j < g.ny; j++) if (g.grid[i][j] > mx) mx = g.grid[i][j];
  for (let i = 0; i < g.nx; i++) for (let j = 0; j < g.ny; j++) { const n = g.grid[i][j]; if (n > 0) { I.push(i); J.push(j); T.push(Math.sqrt(n / mx)); } }
  const xn = new Float64Array(g.nx), yn = new Float64Array(g.ny); for (let i = 0; i < g.nx; i++) for (let j = 0; j < g.ny; j++) { xn[i] += g.grid[i][j]; yn[j] += g.grid[i][j]; }
  return { I, J, T, xn, yn, total: D.sum(xn) };
}

// =====================================================================================================
// The Alignment section: what opens it
// =====================================================================================================
// `more`: the panels that follow the lengths in the same row of charts (identity, calibration), so that a wide
// window lays four of them side by side instead of two pairs with an empty third of the page beside them
export function AlignmentOpening({ b, c, p, more }) {
  const mod = D.mod(b, 'aln'); if (!mod) return null;
  const T = A.alnTotals(b, c.keys, c.filters); if (!T.cells) return null;
  const whole = A.wholeRun(c), pop = mod.population || {}, set = readSetWords(pop, ['pass']), st0 = scopeTag(c, ['samples', 'reads']), st = { scope: [set, st0.scope].filter(Boolean).join(' · '), tag: st0.tag }, tag = tagOf(st, p);
  const withSet = follows => [set, scopeTag(c, follows).scope].filter(Boolean).join(' · ');
  // ---- summary of read alignment: a ring of the records by kind, and its table ----
  const records = T.mapped + T.unmapped + T.secondary + T.supplementary || 1;
  const parts = [{ label: 'Primary, mapped', value: T.mapped, color: tk('--fc-cat-1'), unit: 'reads' }, { label: 'Not mapped', value: T.unmapped, color: tk('--fc-fate-unmapped'), unit: 'reads' },
    { label: 'Secondary records', value: T.secondary, color: tk('--fc-cat-3') }, { label: 'Supplementary records', value: T.supplementary, color: tk('--fc-cat-5') }];
  // a read class is labelled with its share of the reads (the number people ask for), a record class with its count only:
  // a share of all records beside "mapped" would be taken for the mapped share
  for (const q of parts) q.text = q.unit ? F.count(q.value) + ' · ' + small(share(q.value, T.reads)) + ' of reads' : F.count(q.value) + ' records, not reads';
  const ring = { form: 'ring', parts, center: { value: F.count(T.reads), label: 'reads' },
    readout: i => [[parts[i].label, ''], [F.int(parts[i].value), parts[i].unit || 'records']].concat(parts[i].unit ? [[small(share(parts[i].value, T.reads)), 'of the reads']] : ['counted, never read as a read']).concat([[small(parts[i].value / records), 'of all records']]) };
  const recRows = [['Reads (one primary record each)', F.int(T.reads), '100 %', small(T.reads / records)], ['Mapped', F.int(T.mapped), small(share(T.mapped, T.reads)), small(T.mapped / records)],
    ['Not mapped', F.int(T.unmapped), small(share(T.unmapped, T.reads)), small(T.unmapped / records)], ['Secondary records', F.int(T.secondary), F.DASH, small(T.secondary / records)],
    ['Supplementary records', F.int(T.supplementary), F.DASH, small(T.supplementary / records)], ['Reads with a supplementary alignment', F.int(T.withSupp), small(share(T.withSupp, T.reads)), F.DASH], ['All records', F.int(records), F.DASH, '100 %']];
  const recTable = () => ({ head: ['Kind', 'Count', 'Share of reads', 'Share of records'], rows: recRows });
  const ringReading = small(share(T.mapped, T.reads)) + ' of the ' + F.count(T.reads) + ' ' + ofReads(pop).replace(/^the /, '') + ' map to the reference; ' + small(share(T.withSupp, T.mapped)) + ' of the mapped reads continue in a supplementary alignment.' + (T.secondary ? ' The files also hold ' + F.count(T.secondary) + ' secondary records: other places a read could lie.' : '');
  // ---- what the bases of the mapped reads are ----
  let baseSpec = null, baseTable = null, baseReading = null, baseLegend = null;
  if (T.cigar) {
    const noNm = Math.max(0, T.mm - T.mmNm), match = Math.max(0, T.mmNm - T.mism);
    const bp = [{ name: 'Match', value: match, color: tk('--fc-cat-1') }, { name: 'Mismatch', value: T.mism, color: tk('--fc-cat-4') }, { name: 'Inserted', value: T.ins, color: tk('--fc-cat-2') },
      { name: 'Deleted', value: T.del, color: tk('--fc-cat-6') }, { name: 'Soft-clipped', value: T.soft, color: tk('--fc-cat-8') }, { name: 'Hard-clipped', value: T.hard, color: tk('--fc-cat-10') }]
      .concat(noNm ? [{ name: 'Match or mismatch, not told apart', value: noNm, color: tk('--fc-chart-none') }] : []);
    const all = bp.reduce((a, q) => a + q.value, 0) || 1, rest = all - match || 1;
    // two bars, each the whole of what it shows: every position, and below it the positions that are not a match
    // enlarged to the full width. A segment is labelled with its share of its own bar, so length and label agree.
    const rows = [{ label: 'All positions', parts: bp.map(q => ({ value: q.value, color: q.color, label: q.name + ' ' + small(q.value / all), labelShort: small(q.value / all) })) },
      { label: 'The ' + small(rest / all) + ' not a match', parts: bp.slice(1).map(q => ({ value: q.value, color: q.color, label: q.name + ' ' + small(q.value / rest), labelShort: small(q.value / rest) })) }];
    baseSpec = { form: 'hstack', norm: true, segLabels: true, rowH: 40, barH: 24, rows,
      readout: k => { const r = Math.floor(k / 1000), q = bp[(k % 1000) + r]; return [[q.name, ''], [F.int(q.value), 'positions'], [small(q.value / all), 'of all positions']].concat(r ? [[small(q.value / rest), 'of the positions that are not a match']] : []); } };
    baseLegend = bp.map(q => ({ name: q.name + (q.value > 0 ? '' : ' (none)'), color: q.color }));
    baseTable = () => ({ head: ['What the position is', 'Positions', 'Share', 'Per 100 alignment columns'], rows: bp.map(q => [q.name, F.int(q.value), small(q.value / all), /clipped/.test(q.name) ? F.DASH : (q.value / ((T.mm + T.ins + T.del) || 1) * 100).toFixed(3)]),
      note: 'Positions are the alignment columns (match, mismatch, insertion, deletion) and the clipped bases of the primary mapped records. Mismatches are the edit distance (NM) minus inserted and deleted bases.' });
    baseReading = small(match / all) + ' of all positions of the mapped reads are matches; ' + small(T.mism / all) + ' are mismatches, ' + small(T.ins / all) + ' inserted and ' + small(T.del / all) + ' deleted bases, and ' + small(T.soft / all) + ' are soft-clipped ends. The lower bar is the ' + small(rest / all) + ' that are not a match, drawn to the full width: its labels are shares of that bar.';
  }
  // ---- strand of the alignments ----
  const sp = [{ name: 'Forward strand', value: T.fwd, color: tk('--fc-cat-1') }, { name: 'Reverse strand', value: T.rev, color: tk('--fc-cat-3') }], mapped = T.fwd + T.rev || 1;
  const strand = { form: 'hstack', norm: true, segLabels: true, rowH: 40, barH: 24, rows: [{ label: '', parts: sp.map(q => ({ value: q.value, color: q.color, label: q.name + ' ' + F.pct1(q.value / mapped) })) }],
    readout: k => { const q = sp[k % 1000]; return [[q.name, ''], [F.int(q.value), 'mapped reads'], [F.pct1(q.value / mapped), 'of the mapped reads']]; } };
  // ---- depth along the genome ----
  // one track for what is chosen, or one per sample when several are chosen: a genome reads better long, so these
  // are tracks (the whole width of the page), not panels
  const ct = A.contigs(b, c.keys, c.filters);
  // common: { log, top, lo, hi } for tracks that must share one scale (one per sample): log when any of them needs it
  const depthStats = keys => { const dw = A.depthWindows(b, keys, c.filters, false); if (!dw || dw.unavailable) return null; const pos = Array.from(dw.depth).filter(v => v > 0).sort((x, y) => x - y); if (!pos.length) return null;
    const med = pos[pos.length >> 1], max = pos[pos.length - 1]; return { log: med > 0 && max / med > 40, top: Math.max(pos[Math.min(pos.length - 1, Math.floor(pos.length * 0.995))] * 1.15, med * 2), lo: pos[0], hi: max }; };
  const mkDepth = (keys, all, common) => { const dw = A.depthWindows(b, keys, c.filters, all); let depth = null; if (!dw) return null; if (dw.unavailable) return { unavailable: true };
    const ref = dw.ref, n = dw.depth.length, pos = Array.from(dw.depth).filter(v => v > 0).sort((x, y) => x - y), med = pos.length ? pos[pos.length >> 1] : 0, max = pos.length ? pos[pos.length - 1] : 0;
    const log = common ? common.log : c.logY || (med > 0 && max / med > 40), top = log ? null : common ? common.top : (pos.length ? Math.max(pos[Math.min(pos.length - 1, Math.floor(pos.length * 0.995))] * 1.15, med * 2) : 1);
    const strip = ref.contigs.every(q => !q.windows || q.windows < n * 0.012 || /^chr/i.test(q.name)), isBig = q => q.windows >= n * 0.012 || q.name === ref.mitochondrial;
    const blocks = []; for (const q of ref.contigs) { if (isBig(q)) blocks.push({ name: q.name, label: strip ? q.name.replace(/^chr/i, '') : q.name, first: q.first_window, windows: q.windows, big: true }); else { const last = blocks[blocks.length - 1]; if (last && !last.big) last.windows = q.first_window + q.windows - last.first; else blocks.push({ name: 'smaller contigs', label: 'other', first: q.first_window, windows: q.windows, big: false }); } }
    const meanOf = q => { let bs = 0, ln = 0; for (let i = q.first; i < q.first + q.windows; i++) { bs += dw.bases[i]; ln += dw.len[i]; } return ln ? bs / ln : 0; };
    const means = blocks.map(meanOf), total = D.sum(dw.bases), mean = ref.total_length ? total / ref.total_length : 0, zero = n - pos.length, over = top == null ? 0 : pos.filter(v => v > top).length;
    const where = i => { const q = ref.contigs.find(o => i >= o.first_window && i < o.first_window + o.windows), k = i - q.first_window, a = k * ref.window_bases; return { q, a, e: a + dw.len[i] }; };
    const spec = genomeTrack({ contigs: blocks, values: Array.from(dw.depth), means, top, log, title: 'Depth, ×' + (log ? ', logarithmic' : ''), fmt: v => v >= 1000 ? F.int(v) : String(+v.toPrecision(3)), name: 'Depth',
      readout: i => { const w = where(i), blk = blocks.findIndex(q => i >= q.first && i < q.first + q.windows); return [[w.q.name + ': ' + Mb(w.a) + ' to ' + Mb(w.e), ''], [depthX(dw.depth[i]), 'depth'], [F.bases(dw.bases[i]), 'aligned'], [depthX(means[blk]), 'mean of ' + blocks[blk].name]]; } });
    spec.ml = ML; spec.mr = MR;
    if (common && log && common.lo > 0) spec.y.domain = [Math.pow(10, Math.floor(Math.log10(common.lo))), Math.pow(10, Math.ceil(Math.log10(common.hi)))];
    depth = { spec, log, mean, med, zero, n, over, top, ref, blocks, means, run: dw.run,
      reading: 'The genome is read ' + depthX(mean) + ' on average; half of the windows that hold a read lie under ' + depthX(med) + (zero ? ', and ' + F.int(zero) + ' of ' + F.int(n) + ' windows hold none' : '') + '.' + (over ? ' ' + F.int(over) + ' windows rise above the top of the chart (' + depthX(top) + '); their depth is in the read-out and the table.' : '') + (log && !c.logY ? ' Drawn on a logarithmic scale, because a few windows are read far more often than the rest.' : ''),
      table: () => ({ head: ['Contig', 'From', 'To', 'Aligned bases', 'Depth, ×'], rows: Array.from(dw.depth, (v, i) => { const w = where(i); return [w.q.name, F.int(w.a), F.int(w.e), F.int(dw.bases[i]), v >= 10 ? v.toFixed(1) : v.toFixed(3)]; }),
        note: 'Windows of ' + F.int(ref.window_bases) + ' bases; the last window of a contig is shorter. Depth = aligned bases in the window / its length; positions under a match, mismatch or deletion count, skipped regions (introns) do not.' + (ref.other && ref.other.sequences ? ' ' + F.int(ref.other.sequences) + ' short contigs have no window of their own.' : '') }) };
    return depth; };
  let common = null;
  if (c.sKeys.length > 1) { const all = c.sKeys.map(k => depthStats([k])).filter(Boolean); if (all.length) common = { log: c.logY || all.some(q => q.log), top: Math.max(...all.map(q => q.top)), lo: Math.min(...all.map(q => q.lo)), hi: Math.max(...all.map(q => q.hi)) }; }
  const depths = c.sKeys.length > 1 ? c.sKeys.map(k => ({ key: k, d: mkDepth([k], false, common) })).filter(q => q.d) : [{ key: null, d: mkDepth(c.keys, whole) }].filter(q => q.d);
  const depth = depths.length && !depths[0].d.unavailable ? depths[0].d : null, depthNA = depths.length && depths[0].d.unavailable;
  // one sentence about every track when there are several: the mean depth of each sample, on the scale they share
  const depthsReading = depths.length > 1 ? 'One track per sample, all on one scale' + (common && common.log ? ' (logarithmic)' : '') + '. Mean depth: ' + depths.filter(q => !q.d.unavailable).map(q => c.short(q.key, 24) + ' ' + depthX(q.d.mean)).join(', ') + '.' : null;
  // ---- reads and bases per chromosome ----
  let contig = null;
  if (ct && ct.rows.length) {
    const byB = c.by === 'bases', tot = byB ? ct.bases : ct.reads, floor = ct.total ? ct.total * 0.005 : 1e6, big = ct.rows.filter(r => r.name !== 'other contigs' && (r.length == null || r.length >= floor || r.name === ct.mito || ct.rows.length <= 30));
    const restRows = ct.rows.filter(r => !big.includes(r)), rest = { name: 'other contigs', reads: restRows.reduce((a, r) => a + r.reads, 0), bases: restRows.reduce((a, r) => a + r.bases, 0), length: restRows.reduce((a, r) => a + (r.length || 0), 0) };
    const bars = big.concat(rest.reads || rest.bases ? [rest] : []), strip = big.every(r => /^chr/i.test(r.name)), val = r => share(byB ? r.bases : r.reads, tot) * 100;
    const x = bars.map((_, i) => i + 0.5), mi = bars.findIndex(r => r.name === ct.mito), mito = mi >= 0 ? bars[mi] : null, top = bars.slice().sort((p1, p2) => val(p2) - val(p1))[0];
    contig = { spec: { form: 'xy', x: { domain: [0, bars.length], tickMarks: false, ticks: bars.map((r, i) => ({ v: i + 0.5, label: r.name === 'other contigs' ? 'other' : strip ? r.name.replace(/^chr/i, '') : r.name, major: true })), title: strip ? 'Chromosome' : 'Contig' },
        y: { title: 'Share of the ' + (byB ? 'aligned bases' : 'mapped reads') + ', %' }, series: [{ name: 'Share', type: 'vbar', x, parts: [{ y: bars.map((r, i) => i === mi ? 0 : val(r)), color: tk('--fc-cat-1') }, { y: bars.map((r, i) => i === mi ? val(r) : 0), color: tk('--fc-cat-4') }] }], cursor: { x },
        readout: i => [[bars[i].name, ''], [F.count(bars[i].reads), 'reads'], [small(share(bars[i].reads, ct.reads)), 'of the mapped reads'], [F.bases(bars[i].bases), 'aligned']].concat(bars[i].length ? [[depthX(bars[i].bases / bars[i].length), 'mean depth']] : []) },
      legend: mito ? [{ name: 'Chromosomes', color: tk('--fc-cat-1') }, { name: 'Mitochondrial genome (' + mito.name + ')', color: tk('--fc-cat-4') }] : null,
      reading: (mito ? small(share(mito.reads, ct.reads)) + ' of the mapped reads lie on the mitochondrial genome (' + mito.name + '). ' : 'No mitochondrial genome was found among the reference names. ') + 'Most ' + (byB ? 'aligned bases' : 'reads') + ' are on ' + top.name + ' (' + small(val(top) / 100) + ').',
      table: () => ({ head: ['Contig', 'Length, bases', 'Reads', 'Share of mapped reads', 'Aligned bases', 'Mean depth, ×'], rows: ct.rows.map(r => [r.name, r.length ? F.int(r.length) : F.DASH, F.int(r.reads), small(share(r.reads, ct.reads)), F.int(r.bases), r.length ? (r.bases / r.length).toFixed(2) : F.DASH]) }) };
  }
  // ---- read length and aligned length, one over the other ----
  const H = D.readHists(b, c.keys, c.filters), ah = A.hist1(b, 'aln', 'aln_aligned_length_hist', c.keys, c.filters);
  let lens = null;
  if (ah && H.len && D.sum(ah.n) > 0) {
    const pct = (a) => { const t = D.sum(a) || 1; return Array.from(a, v => v / t * 100); }, dom = D.spanQ(H.len.edges, H.len.n), lpop = (D.mod(b, 'len') || {}).population || {};
    const amed = D.quantile(ah.edges, ah.n, 0.5, true), lmed = D.quantile(H.len.edges, H.len.n, 0.5, true), an50 = ah.w ? D.n50(ah.edges, ah.w) : null, inside = share(T.mm + T.ins, T.mm + T.ins + T.soft);
    const ser = [{ name: 'Read length', type: 'outline', edges: H.len.edges, y: pct(H.len.n), color: tk('--fc-chart-annotation') }, { name: 'Aligned length', type: 'outline', edges: ah.edges, y: pct(ah.n), color: tk('--fc-this-run'), width: 2 }];
    lens = { spec: { form: 'xy', x: { domain: dom, log: true, kind: 'length', title: 'Length, logarithmic' }, y: { title: 'Share of the reads per bin, %' }, series: ser, cursor: { x: mids(ah.edges, true) },
        vlines: [amed ? { x: amed, label: 'aligned, median ' + F.bases(amed), strong: true } : null, lmed ? { x: lmed, label: 'read, median ' + F.bases(lmed), side: 'left' } : null].filter(Boolean),
        readout: i => { const k = D.sum(H.len.n) && H.len.edges.length === ah.edges.length ? i : -1; return [[F.int(ah.edges[i]) + ' to ' + F.int(ah.edges[i + 1]) + ' bases', ''], [F.int(ah.n[i]), 'reads aligned over this length']].concat(k >= 0 ? [[F.int(H.len.n[k]), 'reads of this length']] : []); } },
      legend: [{ name: 'Read length (' + readSetWords(lpop, c.filters) + ')', kind: 'line', color: tk('--fc-chart-annotation') }, { name: 'Aligned length (mapped reads)', kind: 'line', color: tk('--fc-this-run') }],
      reading: 'The aligned part of a mapped read is typically ' + F.bases(amed) + ' long' + (an50 ? ' (N50 ' + F.bases(an50) + ')' : '') + (T.cigar ? '; ' + small(inside) + ' of the bases of the mapped reads lie inside their alignment.' : '.'),
      table: () => { const same = H.len.edges.length === ah.edges.length; return { head: ['Length, bases', 'Reads of this length', 'Reads aligned over this length', 'Aligned bases'],
        rows: Array.from(ah.n, (v, i) => [F.int(ah.edges[i]) + ' to ' + F.int(ah.edges[i + 1]), same ? F.int(H.len.n[i]) : F.DASH, F.int(v), ah.w ? F.int(ah.w[i]) : F.DASH]).filter((_, i) => ah.n[i] > 0 || (same && H.len.n[i] > 0)),
        note: 'Reads of this length: ' + readSetWords(lpop, c.filters) + '. Reads aligned over this length: the mapped reads, by the length of the part inside the alignment.' }; } };
  }
  const heat220 = 240;
  return <>
    <Grid>
      <ChartPanel id={pid('aln/records')} title={panelTitle('aln/records')} reading={ringReading} scope={st.scope} tag={tag} height={236} spec={ring} how={how('aln_record_counts')} table={recTable} />
      <figure class="fc-chart" id="p-aln-records-table">
        <figcaption><div class="fc-chart__head"><span class="fc-chart__title">The same, in numbers</span></div>
          <p class="fc-chart__reading">A read counts once, by its primary record. Secondary and supplementary records are counted, never read as reads.</p></figcaption>
        <div class="fc-tablewrap"><table class="fc-table fc-table--wrap"><caption class="fc-sr">Alignment records by kind</caption>
          <thead><tr><th scope="col">Kind</th><th scope="col" class="fc-num">Count</th><th scope="col" class="fc-num">Share of reads</th><th scope="col" class="fc-num">Share of records</th></tr></thead>
          <tbody>{recRows.map(r => <tr><th scope="row">{r[0]}</th><td class="fc-num">{r[1]}</td><td class="fc-num">{r[2]}</td><td class="fc-num">{r[3]}</td></tr>)}</tbody></table></div>
      </figure>
      {baseSpec ? <ChartPanel id={pid('aln/bases')} title={panelTitle('aln/bases')} reading={baseReading} scope={st.scope} tag={tag} height={104} spec={baseSpec} legend={baseLegend} how={how('aln_error_rates')} table={baseTable} /> : null}
      {baseSpec ? <ChartPanel id={pid('aln/strand')} title={panelTitle('aln/strand')} reading={F.pct1(T.fwd / mapped) + ' of the mapped reads lie on the forward strand of the reference and ' + F.pct1(T.rev / mapped) + ' on the reverse strand.'} scope={st.scope} tag={tag} height={64} spec={strand} how={how('aln_strand_frac')}
        table={() => ({ head: ['Strand', 'Mapped reads', 'Share'], rows: sp.map(q => [q.name, F.int(q.value), F.pct1(q.value / mapped)]) })} /> : null}
    </Grid>
    {depth ? <>
      <p class="fc-tracks__hint" id="p-aln-depth-reading"><b>{panelTitle('aln/depth')}.</b> {depthsReading || depth.reading}</p>
      <div class="fc-tracks">{depths.filter(q => !q.d.unavailable).map((q, k) => <Track id={pid('aln/depth') + (k ? '-' + k : '')} title={q.key ? sampleName(c, q.key) : panelTitle('aln/depth')}
        sub={'depth of each window of ' + Mb(q.d.ref.window_bases) + (q.d.log ? ', logarithmic' : '') + ' · line: mean of the chromosome' + (q.key ? ' · mean ' + depthX(q.d.mean) : '')} height={depths.length > 1 ? 130 : 190} spec={q.d.spec}
        tag={tagOf(scopeTag(c, ['samples', 'reads', 'y']), partialInfo(b, 'aln.windows') || p) || null} how={how('aln_depth_windows')} table={q.d.table} />)}</div></>
      : depthNA ? <RunOnly title={panelTitle('aln/depth')} /> : null}
    <Grid>
      {contig ? <ChartPanel id={pid('aln/contigs')} title={panelTitle('aln/contigs')} reading={contig.reading} scope={withSet(['samples', 'reads', 'by'])} tag={tagOf(scopeTag(c, ['samples', 'reads', 'by']), p)} height={heat220} spec={contig.spec} legend={contig.legend} how={how('aln_contig_table')} table={contig.table} /> : null}
      {lens ? <ChartPanel id={pid('aln/lengths')} title={panelTitle('aln/lengths')} reading={lens.reading} scope={st.scope} tag={tag} height={heat220} spec={lens.spec} legend={lens.legend} how={how('aln_aligned_length_hist')} table={lens.table} /> : null}
      {more}
    </Grid>
  </>;
}

// the basecaller's quality against the quality measured from the alignment: stands beside the identity chart
export function Calibration({ b, c, p }) {
  const whole = A.wholeRun(c), st0 = scopeTag(c, ['samples', 'reads']), st = { scope: [readSetWords((D.mod(b, 'aln') || {}).population || {}, ['pass']), st0.scope].filter(Boolean).join(' · '), tag: st0.tag };
  const cg = A.grid2(b, 'aln', 'aln_q_calibration', c.keys, c.filters, whole); if (!cg) return null;
  if (cg.unavailable) return <RunOnly title={panelTitle('aln/calibration')} />;
  const hc = heatCells(cg); if (!hc.total) return null;
  const qs = A.gridQuantiles(cg), mx = D.quantile(cg.xe, hc.xn, 0.5), my = D.quantile(cg.ye, hc.yn, 0.5), dx = D.spanQ(cg.xe, hc.xn, 0.0005, 2), dy = D.spanQ(cg.ye, hc.yn, 0.0005, 2), lo = Math.max(dx[0], dy[0]), hi = Math.min(dx[1], dy[1]), xm = mids(cg.xe);
  const med = new Array(cg.nx).fill(null); for (const q of qs) med[q.i] = q.med;
  const diff = my - mx;
  const spec = { form: 'xy', x: { domain: dx, title: 'Mean quality of the read, from the basecaller, Q', fmt: v => String(v) }, y: { domain: dy, title: 'Quality measured from the alignment, Q', fmt: v => String(v) }, heat: { xe: cg.xe, ye: cg.ye, i: hc.I, j: hc.J, t: hc.T },
    series: [{ name: 'Equal', x: [lo, hi], y: [lo, hi], color: tk('--fc-chart-annotation'), dash: [5, 3], width: 1.5 }, { name: 'Median measured quality', x: xm, y: med, color: tk('--fc-text'), width: 2 }], cursor: { x: xm },
    readout: i => [['Q ' + cg.xe[i] + ' to ' + cg.xe[i + 1] + ' from the basecaller', ''], [F.int(hc.xn[i]), 'reads']].concat(med[i] != null ? [['Q ' + med[i].toFixed(1), 'measured, median'], [pctId(med[i]), 'identity']] : ['too few reads for a median']) };
  const reading = 'A typical read is scored ' + F.q(mx) + ' by the basecaller and measures ' + F.q(my) + ' against the reference (' + pctId(my) + ' identity): the basecaller ' + (Math.abs(diff) < 0.75 ? 'tells the truth' : diff > 0 ? 'is cautious by ' + diff.toFixed(1) + ' Q' : 'promises ' + (-diff).toFixed(1) + ' Q more than the alignment finds') + '.' + (cg.outside ? ' ' + F.count(cg.outside) + ' reads without a single difference have no measured quality and are not on the grid.' : '');
  return <ChartPanel id={pid('aln/calibration')} title={panelTitle('aln/calibration')} reading={reading} scope={st.scope} tag={tagOf(scopeTag(c, ['samples', 'reads'], cg.run && !whole ? 'all samples' : null), p)} height={220} spec={spec} how={how('aln_q_calibration')}
    legend={[{ name: 'Median measured quality at each basecaller quality', kind: 'line', color: tk('--fc-text') }, { name: 'Where the two are equal', kind: 'line', dash: true, color: tk('--fc-chart-annotation') }]}
    table={() => ({ head: ['Quality from the basecaller', 'Reads', 'Measured quality, median', 'Middle half', 'Identity at the median, %'], rows: qs.map(q => ['Q ' + cg.xe[q.i] + ' to ' + cg.xe[q.i + 1], F.int(q.n), q.med.toFixed(1), q.q25.toFixed(1) + ' to ' + q.q75.toFixed(1), q2p(q.med).toFixed(2)]),
      note: 'Measured quality = -10 log10 of the share of alignment columns that differ from the reference. Rows with too few reads for a median (under 50; under 10 in a small run) are left out.' })} />;
}

// =====================================================================================================
// The Alignment section: "More charts"
// =====================================================================================================
// first: a panel of the section that goes in front of these (mapping quality)
export function AlignmentMore({ b, c, p, first }) {
  const mod = D.mod(b, 'aln'); if (!mod) return null;
  const T = A.alnTotals(b, c.keys, c.filters); if (!T.cells || !T.cigar) return <Grid wide={3}>{first}</Grid>;
  const whole = A.wholeRun(c), st0 = scopeTag(c, ['samples', 'reads']), st = { scope: [readSetWords(mod.population || {}, ['pass']), st0.scope].filter(Boolean).join(' · '), tag: st0.tag }, tag = tagOf(st, p);
  // soft clipping
  const sh = A.hist1(b, 'aln', 'aln_soft_clip_hist', c.keys, c.filters);
  let clip = null;
  if (sh && D.sum(sh.n) > 0) {
    const tot = D.sum(sh.n), none = sh.n[0] / tot; let over10 = 0; for (let i = 10; i < sh.n.length; i++) over10 += sh.n[i];
    clip = { spec: { form: 'xy', x: { domain: [0, 1], title: 'Soft-clipped share of the read', fmt: v => Math.round(v * 100) + ' %' }, y: { log: true, title: 'Reads per bin, logarithmic' }, series: [{ name: 'Reads', type: 'hist', edges: sh.edges, y: Array.from(sh.n) }], cursor: { x: mids(sh.edges) },
        readout: i => [[Math.round(sh.edges[i] * 100) + ' to ' + Math.round(sh.edges[i + 1] * 100) + ' % of the read soft-clipped', ''], [F.int(sh.n[i]), 'reads'], [small(sh.n[i] / tot), 'of the mapped reads']] },
      reading: F.pct1(none) + ' of the mapped reads have under 1 % of their bases soft-clipped; ' + small(over10 / tot) + ' have 10 % or more outside their alignment (adapter or primer left on the read, or two molecules joined).',
      table: () => ({ head: ['Soft-clipped share of the read, %', 'Reads'], rows: Array.from(sh.n, (v, i) => [Math.round(sh.edges[i] * 100) + ' to ' + Math.round(sh.edges[i + 1] * 100), F.int(v)]).filter((_, i) => sh.n[i] > 0) }) };
  }
  // the grids
  const heat = (name, o) => {
    const g = A.grid2(b, 'aln', name, c.keys, c.filters, whole); if (!g) return null; if (g.unavailable) return { unavailable: true };
    const hc = heatCells(g); if (!hc.total) return null;
    const xm = mids(g.xe, g.xlog), qs = o.median ? A.gridQuantiles(g) : [], med = new Array(g.nx).fill(null); for (const q of qs) med[q.i] = q.med;
    const dx = D.spanQ(g.xe, hc.xn, 0.0005, 2), dy = D.spanQ(g.ye, hc.yn, 0.0005, 2), yIsId = /identity/.test(g.ylog ? '' : o.y), ser = [];
    if (o.diagonal) { const lo = Math.max(dx[0], dy[0]), hi = Math.min(dx[1], dy[1]); ser.push({ name: 'Equal', x: [lo, hi], y: [lo, hi], color: tk('--fc-chart-annotation'), dash: [5, 3], width: 1.5 }); }
    if (o.median) ser.push({ name: 'Median', x: xm, y: med, color: tk('--fc-text'), width: 2 });
    const yTicks = yIsId ? ID_TICKS.filter(v => v >= dy[0] && v <= dy[1]) : null, xr = i => o.xr(g.xe[i], g.xe[i + 1]);
    // a logarithmic length axis keeps the range of the data (the kit would widen it to whole powers of ten) and is
    // written like the x axis: "100 b", "1 kb", "10 kb"
    const yLen = g.ylog ? lenTicks(dy[0], dy[1]) : null;
    return { g, run: g.run, spec: { form: 'xy', x: { domain: dx, log: g.xlog, kind: g.xlog ? 'length' : 'linear', title: o.xTitle, fmt: v => String(v) }, y: Object.assign({ domain: dy, log: g.ylog, title: o.yTitle }, yTicks && yTicks.length > 1 ? { ticks: yTicks, fmt: tickId } : yLen && yLen.length > 1 ? { ticks: yLen, fmt: lenTick } : { fmt: g.ylog ? undefined : (v => String(v)) }),
        heat: { xe: g.xe, ye: g.ye, i: hc.I, j: hc.J, t: hc.T }, series: ser, cursor: { x: xm },
        readout: i => [[xr(i), ''], [F.int(hc.xn[i]), 'reads']].concat(o.median ? (med[i] != null ? [[o.yv(med[i]), 'median']] : ['too few reads for a median']) : []) },
      legend: (o.median ? [{ name: 'Median at each ' + o.xWord, kind: 'line', color: tk('--fc-text') }] : []).concat(o.diagonal ? [{ name: 'Where the two are equal', kind: 'line', dash: true, color: tk('--fc-chart-annotation') }] : []),
      outside: g.outside, qs,
      table: () => { const rows = []; for (let k = 0; k < hc.I.length; k++) rows.push([o.xr(g.xe[hc.I[k]], g.xe[hc.I[k] + 1]), o.yr(g.ye[hc.J[k]], g.ye[hc.J[k] + 1]), F.int(g.grid[hc.I[k]][hc.J[k]])]); return { head: [o.xHead, o.yHead, 'Reads'], rows }; } };
  };
  const lenR = (a0, a1) => F.int(a0) + ' to ' + F.int(a1) + ' bases', idR = (a0, a1) => q2p(a0).toFixed(2) + ' to ' + q2p(a1).toFixed(2) + ' %', qR = (a0, a1) => 'Q ' + a0 + ' to ' + a1;
  const la = heat('aln_length_vs_aligned', { diagonal: true, xTitle: 'Read length, logarithmic', yTitle: 'Aligned length, logarithmic', xr: lenR, yr: lenR, xHead: 'Read length', yHead: 'Aligned length', y: 'length' });
  const ia = heat('aln_identity_vs_aligned', { median: true, xTitle: 'Aligned length, logarithmic', yTitle: 'Identity, spaced by error rate', xr: lenR, yr: idR, yv: pctId, xWord: 'aligned length', xHead: 'Aligned length', yHead: 'Identity', y: 'identity' });
  const qi = heat('aln_q_calibration', { median: true, xTitle: 'Mean quality of the read, Q', yTitle: 'Identity, spaced by error rate', xr: qR, yr: idR, yv: pctId, xWord: 'read quality', xHead: 'Read quality', yHead: 'Identity', y: 'identity' });
  // accuracy by read length: the median and the middle half of the identity at each aligned length
  let acc = null;
  if (ia && !ia.unavailable && ia.qs.length > 1) {
    const g = ia.g, xs = ia.qs.map(q => Math.sqrt(g.xe[q.i] * g.xe[q.i + 1])), med = ia.qs.map(q => q.med), lo = ia.qs.map(q => q.q25), hi = ia.qs.map(q => q.q75);
    const dom = [g.xe[ia.qs[0].i], g.xe[ia.qs[ia.qs.length - 1].i + 1]], y0 = Math.floor(Math.min(...lo) - 0.5), y1 = Math.ceil(Math.max(...hi) + 0.5), ticks = ID_TICKS.concat([17, 26]).sort((p1, p2) => p1 - p2).filter(v => v >= y0 && v <= y1);
    const best = ia.qs.reduce((p1, q) => q.n > p1.n ? q : p1, ia.qs[0]), first = ia.qs[0], last = ia.qs[ia.qs.length - 1];
    acc = { spec: { form: 'xy', x: { domain: dom, log: true, kind: 'length', title: 'Aligned length, logarithmic' }, y: { domain: [y0, y1], ticks: ticks.length > 1 ? ticks : undefined, fmt: tickId, title: 'Identity, spaced by error rate' },
        series: [{ name: 'Median identity', x: xs, y: med, band: { lo, hi }, color: tk('--fc-this-run') }], cursor: { x: xs },
        readout: i => [[lenR(g.xe[ia.qs[i].i], g.xe[ia.qs[i].i + 1]), ''], [pctId(med[i]), 'median identity'], [pctId(lo[i]) + ' to ' + pctId(hi[i]), 'middle half'], [F.int(ia.qs[i].n), 'reads']] },
      legend: [{ name: 'Median identity', kind: 'line', color: tk('--fc-this-run') }, { name: 'Middle half of the reads', kind: 'band', color: tk('--fc-this-run') }],
      reading: 'Reads aligned over about ' + F.bases(Math.sqrt(g.xe[best.i] * g.xe[best.i + 1])) + ' (the most common length) match the reference at ' + pctId(best.med) + '; the shortest alignments drawn reach ' + pctId(first.med) + ' and the longest ' + pctId(last.med) + '.',
      table: () => ({ head: ['Aligned length, bases', 'Reads', 'Median identity, %', 'Middle half, %'], rows: ia.qs.map(q => [F.int(g.xe[q.i]) + ' to ' + F.int(g.xe[q.i + 1]), F.int(q.n), q2p(q.med).toFixed(2), q2p(q.q25).toFixed(2) + ' to ' + q2p(q.q75).toFixed(2)]), note: 'Lengths with too few reads for a median (under 50; under 10 in a small run) are left out.' }) };
  }
  const gtag = h => tagOf(scopeTag(c, ['samples', 'reads'], h && h.run && !whole ? 'all samples' : null), p);
  const panel = (id, h, reading, hid) => !h ? null : h.unavailable ? <RunOnly title={panelTitle(id)} /> : <ChartPanel id={pid(id)} title={panelTitle(id)} reading={reading + (h.outside ? ' ' + F.count(h.outside) + ' reads fall outside the axes (a read without a single difference has no place on the identity axis).' : '')} scope={st.scope} tag={gtag(h)} height={240} spec={h.spec} legend={h.legend} how={how(hid)} table={h.table} />;
  return <>
    <Grid wide={3}>
      {first}
      {clip ? <ChartPanel id={pid('aln/soft_clip')} title={panelTitle('aln/soft_clip')} reading={clip.reading} scope={st.scope} tag={tag} height={240} spec={clip.spec} how={how('aln_soft_clip_hist')} table={clip.table} /> : null}
      {acc ? <ChartPanel id={pid('aln/identity_by_length')} title={panelTitle('aln/identity_by_length')} reading={acc.reading} scope={st.scope} tag={gtag(ia)} height={240} spec={acc.spec} legend={acc.legend} how={how('aln_identity_by_length')} table={acc.table} /> : null}
      {panel('aln/length_vs_aligned', la, 'Reads on the dashed line align from end to end; reads under it lost part of their length to clipping.', 'aln_length_vs_aligned')}
      {panel('aln/identity_vs_aligned', ia, 'Where the mapped reads sit by the length of their alignment and their identity with the reference. Under about 300 bases the cells form a chequered pattern: a short alignment can only have a few different numbers of differences, so only some identities are possible there.', 'aln_identity_vs_aligned')}
      {panel('aln/quality_vs_identity', qi, 'Where the mapped reads sit by the basecaller’s quality and their identity with the reference.', 'aln_quality_vs_identity')}
    </Grid>
  </>;
}

// =====================================================================================================
// Reads over run time: tracks on the time axis of the run timeline
// =====================================================================================================
// lines of one quantity (A.overTime) on one common list of hours
function onGrid(ot) {
  const ts = [...new Set(ot.lines.flatMap(l => l.t.map(v => +v.toFixed(4))))].sort((p1, p2) => p1 - p2), at = new Map(ts.map((t, i) => [t, i]));
  return { t: ts, lines: ot.lines.map(l => { const z = () => new Array(ts.length).fill(null), o = { key: l.key, n: z(), q25: z(), med: z(), q75: z() }; l.t.forEach((t, k) => { const i = at.get(+t.toFixed(4)); o.n[i] = l.n[k]; o.q25[i] = l.q25[k]; o.med[i] = l.med[k]; o.q75[i] = l.q75[k]; }); return o; }) };
}
export function ReadsOverTime({ b, c, evRule, evSec }) {
  const dynm = D.mod(b, 'dyn.basecalled'), alt = D.mod(b, 'aln.time'); if (!dynm && !alt) return null;
  const whole = A.wholeRun(c), clock = A.clockOf(b, 'dyn.basecalled') || A.clockOf(b, 'aln.time'), pop = (dynm || alt).population || {}, set = readSetWords(pop, c.filters), p = partialInfo(b, 'dyn.basecalled') || partialInfo(b, 'aln.time');
  const runEnd = D.runHours(b);
  const xAxis = labels => ({ domain: c.win, kind: 'time', labels });
  const names = k => k == null ? 'all samples' : sampleName(c, k);
  const tracks = [];
  // ---- output, from the start time of every read ----
  const byB = c.by === 'bases', unit = byB ? 1e9 : 1e6, uw = byB ? 'Gb' : 'M reads';
  if (dynm) {
    let t = null, cur = null, all = null, step = 'ten minutes';       // or 'hour'
    if (whole) { const o = A.outputFromReads(b); if (o) { t = c.out === 'cum' ? o.t : o.th; const main = o.hasPass ? o.pass : o.fail, k = byB ? 'bases' : 'reads'; cur = Array.from(main[k]); if (o.hasPass && o.hasFail && c.filters.length > 1) all = cur.map((v, i) => v + o.fail[k][i]); else if (!c.filters.includes('pass') && o.hasFail) cur = Array.from(o.fail[k]); } }
    else { const o = A.outputOfSamples(b, c.keys, c.filters); if (o) { t = o.t; cur = Array.from(byB ? o.bases : o.reads); step = 'hour'; } else tracks.push({ id: 'reads/output_over_time', na: 'Kept for all samples added up when a run has more than twelve samples: choose “all samples” to see it.' }); }
    if (t && cur) {
      // Steps are cut on the clock, so the first step begins before the run and the last may end after it: a rate per
      // hour is worked out over the part of the step the run covers, or the run would seem to start and end slowly.
      const stepH = step === 'hour' ? 1 : 1 / 6, off = clock ? (step === 'hour' ? clock.offset_3600s_s : clock.offset_600s_s) / 3600 : 0, known = !!(clock && clock.origin_known);
      const startOf = i => off + i * stepH, endOf = i => off + (i + 1) * stepH, covered = i => known ? Math.max(1 / 6, (Math.min(endOf(i), Math.max(runEnd, startOf(i))) - Math.max(0, startOf(i))) / stepH) : 1;
      const cum = a => { let s = 0; return a.map(v => (s += v) / unit); }, perH = a => a.map((v, i) => v / unit / stepH / Math.min(1, covered(i))), f = c.out === 'cum' ? cum : perH;
      const y = f(all || cur), y2 = all ? f(cur) : null, total = D.sum(all || cur), last = t[t.length - 1], raw = all || cur;
      let half = null; { let s = 0; const a = all || cur; for (let i = 0; i < a.length; i++) { s += a[i]; if (s >= total / 2) { half = t[i]; break; } } }
      const ser = [{ name: all ? 'All reads' : 'Reads', x: t, y, type: c.out === 'cum' ? 'line' : 'step' }].concat(y2 ? [{ name: 'Passed reads', x: t, y: y2, color: tk('--fc-chart-annotation'), dash: [5, 3], width: 1.5 }] : []);
      tracks.push({ id: 'reads/output_over_time', sub: (byB ? 'basecalled bases, Gb' : 'reads, millions') + (c.out === 'cum' ? ', added up' : ' per hour') + ' · steps of ' + (step === 'hour' ? 'one hour' : 'ten minutes') + (c.out === 'cum' ? '' : ' · the first and the last step are rated over the part the run covers'), height: 130, how: how('dyn_output_series'),
        legend: y2 ? [{ name: 'All reads', kind: 'line', color: tk('--fc-this-run') }, { name: 'Reads that passed the quality filter', kind: 'line', dash: true, color: tk('--fc-chart-annotation') }] : null,
        spec: { form: 'xy', sync: 'time', ml: ML, mr: MR, x: xAxis(false), y: {}, series: ser, cursor: { x: t }, onBrush: brush,
          readout: i => [[c.out === 'cum' ? 'Run time 0:00 to ' + hm(Math.min(endOf(i), Math.max(runEnd, endOf(i) - stepH))) : stepWords(startOf(i), endOf(i)), ''], [F.sig(y[i]) + ' ' + uw + (c.out === 'cum' ? '' : ' per hour'), ''], all ? 'all reads' : lower(set)].concat(y2 ? [[F.sig(y2[i]) + ' ' + uw, 'passed']] : []) },
        table: () => { let s1 = 0, s2 = 0; return { head: ['Run time, from', 'to (h:mm)', (byB ? 'Bases' : 'Reads') + ' in the step', 'Added up'].concat(y2 ? ['Passed, in the step', 'Passed, added up'] : []),
          rows: raw.map((v, i) => [hm(startOf(i)), hm(endOf(i)), F.int(v), F.int(s1 += v)].concat(y2 ? [F.int(cur[i]), F.int(s2 += cur[i])] : [])),
          note: 'Steps are cut on the clock: the first one begins before the run did' + (known ? ' (it is written from 0:00)' : '') + '.' }; } });
      tracks.half = half; tracks.total = total; tracks.last = last;
    }
  }
  // ---- a median with its middle half, hour by hour ----
  const quart = (id, sec, name, o) => {
    const ot = A.overTime(b, sec, name, c); if (!ot || !ot.lines.length) { if (ot && ot.unavailable) tracks.push({ id, na: 'Kept for all samples added up when a run has more than twelve samples: choose “all samples” to see it.' }); return null; }
    const G = onGrid(ot), many = G.lines.length > 1, fmt = o.fmt;
    if (!G.t.length) { tracks.push({ id, na: 'No hour of the run holds 30 reads of this selection: no median is drawn.' }); return null; }
    const ser = G.lines.map((l, k) => ({ name: names(l.key), x: G.t, y: l.med, color: many ? colorOf(c, l.key) : tk('--fc-this-run'), band: many ? null : { lo: l.q25, hi: l.q75 }, width: 2 })).concat(o.extra || []);
    const lo = Math.min(...G.lines.flatMap(l => l.q25.filter(v => v != null))), hi = Math.max(...G.lines.flatMap(l => l.q75.filter(v => v != null)));
    const first = G.lines[0], firstMed = first.med.find(v => v != null), lastMed = first.med.slice().reverse().find(v => v != null);
    tracks.push({ id, sub: o.sub + (many ? '' : ' · line: median, band: middle half') + (ot.run || many ? '' : ' · ' + names(first.key)), height: o.height || 112, how: how(o.how),
      legend: many ? ser.slice(0, G.lines.length).map(q => ({ name: q.name, label: q.name, kind: 'line', color: q.color })).concat((o.legend || []).filter(q => q.dash)) : o.legend || null,
      spec: { form: 'xy', sync: 'time', ml: ML, mr: MR, x: xAxis(false), y: Object.assign({ zero: false }, o.y ? o.y(lo, hi) : {}), series: ser, cursor: { x: G.t }, onBrush: brush,
        readout: i => [[stepWords(G.t[i] - 0.5, G.t[i] + 0.5), '']].concat(many ? G.lines.map(l => [l.med[i] == null ? 'too few reads' : fmt(l.med[i]), c.short(l.key, 16)])
          : first.med[i] == null ? ['too few reads'] : [[fmt(first.med[i]), 'median'], [fmt(first.q25[i]) + ' to ' + fmt(first.q75[i]), 'middle half'], [F.int(first.n[i]), 'reads']], o.readExtra ? o.readExtra(G.t[i]) : []) },
      table: () => ({ head: ['Run time, from', 'to (h:mm)'].concat(G.lines.flatMap(l => many ? [names(l.key) + ', median'] : ['Reads', 'Lower quartile', 'Median', 'Upper quartile'])),
        rows: G.t.map((t, i) => [hm(t - 0.5), hm(t + 0.5)].concat(G.lines.flatMap(l => many ? [l.med[i] == null ? F.DASH : fmt(l.med[i])] : [l.n[i] == null ? F.DASH : F.int(l.n[i]), l.q25[i] == null ? F.DASH : fmt(l.q25[i]), l.med[i] == null ? F.DASH : fmt(l.med[i]), l.q75[i] == null ? F.DASH : fmt(l.q75[i])]))),
        note: 'One row per hour of the clock; hours with fewer than 30 reads are left out.' }) });
    return { firstMed, lastMed };
  };
  const idAxis = (lo, hi) => { const a = Math.floor(lo - 0.3), e = Math.ceil(hi + 0.3), ticks = [3, 7, 10, 13, 17, 20, 23, 26, 30, 40].filter(v => v >= a && v <= e); return { domain: [a, e], ticks: ticks.length > 1 ? ticks : undefined, fmt: tickId }; };
  const len = dynm ? quart('reads/length_over_time', 'dyn.basecalled', 'len_over_time', { sub: 'read length, bases', fmt: v => F.bases(v), how: 'len_over_time' }) : null;
  const qual = dynm ? quart('reads/quality_over_time', 'dyn.basecalled', 'qual_over_time', { sub: 'mean read quality, Q', fmt: v => 'Q ' + v.toFixed(1), how: 'qual_over_time' }) : null;
  if (alt) {
    quart('reads/aligned_length_over_time', 'aln.time', 'aln_aligned_length_over_time', { sub: 'aligned part of the mapped reads, bases', fmt: v => F.bases(v), how: 'aln_aligned_length_over_time' });
    quart('reads/identity_over_time', 'aln.time', 'aln_identity_over_time', { sub: 'identity of the mapped reads with the reference', fmt: pctId, y: idAxis, how: 'aln_identity_over_time' });
  }
  if (dynm) {
    const est = D.speed(b), extra = est ? [{ name: 'Estimated from the signal', x: est.t, y: est.med, color: tk('--fc-chart-annotation'), dash: [5, 3], width: 1.5 }] : [];
    quart('reads/speed_over_time', 'dyn.basecalled', 'dyn_speed_over_time', { sub: 'bases per second, measured on the reads', fmt: v => Math.round(v) + ' bases/s', how: 'dyn_speed_over_time', extra,
      legend: est ? [{ name: 'Measured: bases of the read / seconds in the pore', kind: 'line', color: tk('--fc-this-run') }, { name: 'Estimated from the signal before basecalling (run timeline)', kind: 'line', dash: true, color: tk('--fc-chart-annotation') }] : null,
      readExtra: est ? (t => { const i = Math.max(0, Math.min(est.t.length - 1, Math.floor(t))); return est.med[i] == null ? [] : [[Math.round(est.med[i]) + ' bases/s', 'estimated']]; }) : null });
  }
  const shown = tracks.filter(Boolean);
  if (!shown.length) return c.sKeys.length ? <Section id="reads_over_time" title="Reads over run time" x={{ evRule, evSec }}><NoneFor c={c} what="reads with a start time" /></Section> : null;
  shown.forEach((t, i) => { if (t.spec) t.spec.x = xAxis(i === shown.length - 1); });
  const off = clock ? Math.round(-clock.offset_3600s_s / 60) : 0;
  const reading = (tracks.total ? 'Half of the ' + (byB ? 'bases' : 'reads') + ' of these files started in the first ' + F.hours(tracks.half) + '.' : '') + (len && len.firstMed && len.lastMed ? ' Median read length went from ' + F.bases(len.firstMed) + ' in the first hour drawn to ' + F.bases(len.lastMed) + ' in the last' : '') + (qual && qual.firstMed != null ? (len ? ', read quality from ' : ' Median read quality went from ') + F.q(qual.firstMed) + ' to ' + F.q(qual.lastMed) + '.' : len ? '.' : '');
  const speedOnTrimmed = dynm && dynm.params && dynm.params.speed_on_trimmed_reads;
  // the alignment was read, but its reads carry no start time: say why two tracks are not there
  const alnM = D.mod(b, 'aln'), alnRead = alnM && (alnM.cells || []).length > 0, alnTimed = alt && (alt.cells || []).some(q => (q.scalars || {}).reads_with_start_time > 0);
  const noAlnTime = alnRead && !alnTimed ? 'Aligned length and identity over run time are not drawn: the reads of the alignment files of this run carry no start time (no st tag on their records).' : null;
  const qWhole = dynm && dynm.params && dynm.params.qual_whole_numbers;
  return <Section id="reads_over_time" title="Reads over run time" reading={reading.trim() || null} x={{ evRule, evSec }}>
    {p ? <p class="fc-section__note fc-notice fc-partial" role="status">{p.text}</p> : null}
    {variantNote(b, ['dyn.basecalled', 'aln.time']) ? <p class="fc-section__note fc-notice fc-variant" role="status">{variantNote(b, ['dyn.basecalled', 'aln.time'])}</p> : null}
    <p class="fc-tracks__hint">The same time axis and the same cursor as the run timeline: hours since the start of sequencing, from the start time each read carries. These are {set}{whole ? '' : ' of ' + chosenWords(c)}{c.custom && (c.minq || c.minlen) ? ', not cut by the custom filter' : ''}.
      {clock && clock.origin_known ? ' Reads are grouped by the hour of the clock (ten minutes for output), so a step begins ' + (off ? off + ' min before' : 'on') + ' a whole hour of run time.' : ' The start of the run is not known for these files: hours count from the first read.'}
      {speedOnTrimmed ? ' Speed is measured on reads that were trimmed after basecalling, so it reads a little low.' : ''}
      {qWhole ? ' These files give read quality as whole numbers, so the quality of an hour moves in steps of 1 Q.' : ''}{noAlnTime ? ' ' + noAlnTime : ''}</p>
    <div class="fc-tracks">
      {shown.map(t => t.na ? <div class="fc-track fc-track--na" id={pid(t.id)}><div class="fc-track__label"><div class="fc-track__title">{panelTitle(t.id)}</div></div><div class="fc-track__plot fc-muted">{t.na}</div></div>
        : <Track id={pid(t.id)} title={panelTitle(t.id)} sub={t.sub} height={t.height} spec={t.spec} legend={t.legend} how={t.how} table={t.table} tag={p ? p.tag : null} />)}
    </div>
    <Disclosure id="n-time" c={c} title={'All numbers (' + countNumbers(b, ['dyn.basecalled', 'aln.time']) + ')'}><AllNumbers b={b} mods={['dyn.basecalled', 'aln.time']} /></Disclosure>
  </Section>;
}

// =====================================================================================================
// Transcripts: poly(A) tail, coverage along the transcript body, saturation, reads per transcript
// =====================================================================================================
// the classes of "reads per transcript" start at 1, 2, 4, 8 ...: every fourth one is labelled (1, 16, 256, 4,096 ...),
// so the labels are evenly spaced and none is rounded to a number that is not a class
const rptTicks = x => x.map((v, i) => ({ v: i + 0.5, label: i % 4 === 0 ? F.int(v) : '', major: i % 4 === 0 }));
const RPT_TITLE = 'Reads per transcript: a class starts at the number shown and is twice as wide as the one before';
export function Transcripts({ b, c, evRule, evSec }) {
  const pm = D.mod(b, 'rna.polya'), tx = D.mod(b, 'rna.tx'), cov = D.mod(b, 'rna.coverage'); if (!pm && !tx && !cov) return null;
  const whole = A.wholeRun(c), p = partialInfo(b, 'rna.tx') || partialInfo(b, 'rna.polya'), txSet = readSetWords(((tx || pm || cov) || {}).population || {}, c.filters);
  const st0 = scopeTag(c, ['samples', 'reads', 'overlay']), st = { scope: [txSet, st0.scope].filter(Boolean).join(' · '), tag: st0.tag }, tag = tagOf(st, p);
  const names = k => k == null ? 'All samples' : sampleName(c, k), overlay = c.sKeys.length > 1;
  // ---- poly(A) tail ----
  const P = pm ? A.polya(b, c.keys, c.filters) : null;
  let tail = null;
  if (P && P.hist && D.sum(P.hist.n) > 0) {
    const h = P.hist, tot = D.sum(h.n) + h.over, med = D.quantile(h.edges, h.n, Math.min(1, 0.5 * tot / (D.sum(h.n) || 1))), dom = D.spanQ(h.edges, h.n, 0.002, 4), xm = mids(h.edges), set = readSetWords(pm.population || {}, c.filters);
    let ser, legend = null;
    if (overlay) { ser = c.sKeys.map(k => { const q = A.polya(b, [k], c.filters); if (!q || !q.hist) return null; const t = D.sum(q.hist.n) || 1; return { name: names(k), type: 'outline', edges: q.hist.edges, y: Array.from(q.hist.n, v => v / t * 100), color: colorOf(c, k) }; }).filter(Boolean); legend = ser.map(q => ({ name: q.name, label: q.name, kind: 'line', color: q.color })); }
    else { const u = Math.max(...h.n) >= 1e6 ? 1e6 : 1; ser = [{ name: 'Reads', type: 'hist', edges: h.edges, y: Array.from(h.n, v => v / u), u }]; }
    tail = { med, spec: { form: 'xy', x: { domain: [0, Math.max(dom[1], 50)], title: 'Poly(A) tail length, nucleotides', fmt: v => String(v) }, y: overlay ? { title: 'Share of each sample’s reads with a tail, %' } : { title: 'Reads per bin' + (ser[0].u === 1e6 ? ', millions' : '') }, series: ser,
        vlines: overlay ? [] : [{ x: med, label: 'median ' + Math.round(med) + ' nt', strong: true }], cursor: { x: xm },
        readout: i => [[h.edges[i] + ' to ' + h.edges[i + 1] + ' nucleotides', '']].concat(overlay ? ser.map(q => [q.y[i].toFixed(2) + ' %', c.shorten(q.name, 18)]) : [[F.int(h.n[i]), 'reads']]) },
      legend, reading: 'A tail was found in ' + F.pct1(share(P.found, P.reads)) + ' of the ' + F.count(P.reads) + ' reads; its typical length is ' + Math.round(med) + ' nucleotides.' + (P.zero + P.none ? ' ' + F.count(P.zero) + ' reads carry a tail of zero and ' + F.count(P.none) + ' no estimate.' : '') + (P.reads > P.tagged ? ' ' + F.count(P.reads - P.tagged) + ' reads carry no estimate tag at all.' : '') + ' These are ' + set + '.',
      table: () => ({ head: ['Tail length, nucleotides', 'Reads'], rows: Array.from(h.n, (v, i) => [h.edges[i] + ' to ' + h.edges[i + 1], F.int(v)]).filter((_, i) => h.n[i] > 0).concat(h.over ? [['1,000 and longer', F.int(h.over)]] : []).concat([['Tail of zero', F.int(P.zero)], ['No tail found', F.int(P.none)], ['No estimate tag', F.int(P.reads - P.tagged)]]) }) };
  }
  // ---- coverage along the transcript body ----
  const BC = cov ? A.bodyCoverage(b, c.keys) : null;
  let body = null;
  if (BC && BC.some(q => q.reads > 0)) {
    // four length classes: four colours (no two of them green) and, for the two longer ones, a dashed line, so that
    // the classes are never told apart by colour alone
    const BODY_DASH = [null, null, [6, 3], [2, 2]];
    const xs = Array.from({ length: 100 }, (_, i) => i + 0.5), cols = ['--fc-cat-1', '--fc-cat-2', '--fc-cat-4', '--fc-cat-10'], cl = BC.filter(q => q.reads >= 20).map((q, k) => ({ ...q, color: tk(cols[BC.indexOf(q)]), dash: BODY_DASH[BC.indexOf(q)], pct: q.cover.map(v => v == null ? null : v * 100) }));
    const long = cl[cl.length - 1], short = cl[0], start = q => q.pct.slice(0, 5).reduce((a, v) => a + v, 0) / 5, end = q => q.pct.slice(95).reduce((a, v) => a + v, 0) / 5, ends = q => Math.round(start(q)) + ' % of the reads cover the 5′ end and ' + Math.round(end(q)) + ' % the 3′ end';
    body = { spec: { form: 'xy', x: { domain: [0, 100], ticks: [0, 25, 50, 75, 100].map(v => ({ v, label: v === 0 ? '5′ end' : v === 100 ? '3′ end' : v + ' %', major: true })), title: 'Position along the transcript' }, y: { domain: [0, 100], title: 'Share of the reads that cover the position, %' },
        series: cl.map(q => ({ name: q.label, x: xs, y: q.pct, color: q.color, dash: q.dash })), cursor: { x: xs }, readout: i => [[i + ' to ' + (i + 1) + ' % of the way from the 5′ end', '']].concat(cl.map(q => [q.pct[i].toFixed(1) + ' %', q.label])) },
      legend: cl.map(q => ({ name: 'Transcripts ' + q.label + ' (' + F.count(q.reads) + ' reads)', kind: 'line', color: q.color, dash: !!q.dash })),
      reading: cl.length ? 'On transcripts ' + short.label + ', ' + ends(short) + (long !== short ? '; on transcripts ' + long.label + ', ' + ends(long) : '') + ' (the first and the last twentieth of the transcript).' : null,
      table: () => ({ head: ['Position from the 5′ end, %'].concat(cl.map(q => q.label + ', % of reads')), rows: xs.map((_, i) => [i + ' to ' + (i + 1)].concat(cl.map(q => q.pct[i].toFixed(2)))), note: 'A read covers a position when its primary alignment covers the middle of that hundredth of the transcript. Classes with fewer than 20 reads are left out.' }) };
  }
  // ---- saturation, and reads per transcript ----
  const SC = tx ? A.txCurves(b, c, 'rna_saturation_curve') : [], RP = tx ? A.txCurves(b, c, 'rna_reads_per_transcript') : [];
  let sat = null, rpt = null;
  if (SC.length) {
    const one = SC.length === 1, cv = SC[0].cv, last = cv.x.length - 1, u = cv.x[last] >= 1e6 ? 1e6 : 1e3, uw = u === 1e6 ? 'millions' : 'thousands';
    const ser = one ? [{ name: 'Transcripts', x: cv.x.map(v => v / u), y: cv.y.transcripts, color: tk('--fc-this-run'), dots: true }].concat(cv.y.genes ? [{ name: 'Genes', x: cv.x.map(v => v / u), y: cv.y.genes, color: tk('--fc-cat-3'), dots: true }, { name: 'Genes with 10 reads or more', x: cv.x.map(v => v / u), y: cv.y.genes_10_reads, color: tk('--fc-cat-3'), dash: [5, 3], width: 1.5 }] : [])
      : SC.map(q => ({ name: names(q.key), x: q.cv.x.map(v => v / u), y: q.cv.y.transcripts, color: colorOf(c, q.key) }));
    const gain = k => { const y = cv.y[k]; return y && y[last] ? (y[last] - y[last - 1]) / y[last] : null; }, g = gain('transcripts');
    const xmax = Math.max(...SC.map(q => q.cv.x[q.cv.x.length - 1])) / u;
    sat = { spec: { form: 'xy', x: { domain: [0, xmax], title: 'Mapped reads sequenced, ' + uw, fmt: v => String(v) }, y: { title: one ? 'Found (at least one read)' : 'Transcripts found' }, series: ser, cursor: { x: ser[0].x },
        readout: i => [[F.count(SC[0].cv.x[i]) + ' reads', '']].concat(one ? [[F.int(cv.y.transcripts[i]), 'transcripts']].concat(cv.y.genes ? [[F.int(cv.y.genes[i]), 'genes'], [F.int(cv.y.genes_10_reads[i]), 'genes with 10 reads or more']] : []) : SC.map(q => [q.cv.y.transcripts[i] == null ? F.DASH : F.int(q.cv.y.transcripts[i]), c.short(q.key, 16)])) },
      legend: ser.map(q => ({ name: q.name, label: q.name, kind: 'line', dash: !!q.dash, color: q.color })),
      reading: one ? F.int(cv.y.transcripts[last]) + ' transcripts' + (cv.y.genes ? ' of ' + F.int(cv.y.genes[last]) + ' genes' : '') + ' have at least one of the ' + F.count(cv.x[last]) + ' mapped reads. The last tenth of the reads added ' + small(g) + ' more transcripts' + (cv.y.genes ? ' and ' + small(gain('genes')) + ' more genes' : '') + ': ' + (g > 0.02 ? 'more reads would still find more transcripts.' : 'the curve is close to flat.')
        : 'One line per sample: the transcripts each would show with fewer reads.',
      table: () => ({ head: ['Mapped reads sequenced'].concat(one ? ['Transcripts found'].concat(cv.y.genes ? ['Genes found', 'Genes with 10 reads or more'] : []) : SC.map(q => names(q.key))),
        rows: one ? cv.x.map((x, i) => [F.int(x), F.int(cv.y.transcripts[i])].concat(cv.y.genes ? [F.int(cv.y.genes[i]), F.int(cv.y.genes_10_reads[i])] : [])) : SC[0].cv.x.map((_, i) => [Math.round((i + 1) / SC[0].cv.x.length * 100) + ' % of the reads'].concat(SC.map(q => q.cv.y.transcripts[i] == null ? F.DASH : F.int(q.cv.y.transcripts[i])))),
        note: 'Whole numbers. Expected numbers, worked out exactly from the read count of every transcript: what a random part of the mapped reads of that size would have found. No read is drawn and the order of the reads in the run plays no part.' }) };
  }
  if (RP.length > 1) {
    const n = Math.max(...RP.map(q => q.cv.x.length)), xs = Array.from({ length: n }, (_, i) => i + 0.5), x0 = RP.reduce((p1, q) => q.cv.x.length > p1.length ? q.cv.x : p1, RP[0].cv.x), rng = i => i === 0 ? '1 read' : x0[i] + ' to ' + (x0[i] * 2 - 1) + ' reads';
    const ser = RP.map(q => { const t = D.sum(q.cv.y.transcripts) || 1; return { name: names(q.key), x: xs, y: xs.map((_, i) => (q.cv.y.transcripts[i] || 0) / t * 100), color: colorOf(c, q.key), dots: true }; });
    rpt = { spec: { form: 'xy', x: { domain: [0, n], tickMarks: false, ticks: rptTicks(x0), title: RPT_TITLE }, y: { title: 'Share of each sample’s transcripts, %' },
        series: ser, cursor: { x: xs }, readout: i => [[rng(i), '']].concat(ser.map(q => [q.y[i].toFixed(1) + ' %', c.shorten(q.name, 16)])) },
      legend: ser.map(q => ({ name: q.name, label: q.name, kind: 'line', color: q.color })), reading: 'One line per sample: the share of its transcripts that received so many reads.',
      table: () => ({ head: ['Reads per transcript'].concat(RP.map(q => names(q.key) + ', transcripts')), rows: xs.map((_, i) => [rng(i)].concat(RP.map(q => F.int(q.cv.y.transcripts[i] || 0)))) }) };
  }
  if (RP.length === 1) {
    const cv = RP[0].cv, n = cv.x.length, tr = cv.y.transcripts, rd = cv.y.reads, tt = D.sum(tr) || 1, tr2 = D.sum(rd) || 1, xs = cv.x.map((_, i) => i + 0.5), rng = i => i === 0 ? '1 read' : cv.x[i] + ' to ' + (cv.x[i] * 2 - 1) + ' reads';
    let top = 0, topReads = 0; for (let i = n - 1; i >= 0 && top < tt * 0.01; i--) { top += tr[i]; topReads += rd[i]; }
    rpt = { spec: { form: 'xy', x: { domain: [0, n], tickMarks: false, ticks: rptTicks(cv.x), title: RPT_TITLE }, y: { title: 'Transcripts' },
        series: [{ name: 'Transcripts', type: 'vbar', x: xs, parts: [{ y: tr, color: tk('--fc-this-run') }] }], cursor: { x: xs }, readout: i => [[rng(i), ''], [F.int(tr[i]), 'transcripts'], [small(tr[i] / tt), 'of the transcripts found'], [small(rd[i] / tr2), 'of the mapped reads']] },
      reading: F.pct1(tr[0] / tt) + ' of the transcripts found have a single read; the ' + F.int(top) + ' transcripts with the most reads (' + small(top / tt) + ') hold ' + F.pct1(topReads / tr2) + ' of the mapped reads.',
      table: () => ({ head: ['Reads per transcript', 'Transcripts', 'Reads they hold'], rows: cv.x.map((_, i) => [rng(i), F.int(tr[i]), F.int(rd[i])]) }) };
  }
  const S = D.summary(b, 'rna.tx'), TT = tx ? A.alnTotals(b, c.keys, c.filters, 'rna.tx') : null, cap = text => p ? text + ' · ' + p.word : text;
  // why the charts of the alignment to the transcripts are not there: its files wait to be read (nothing was counted
  // yet), or the counts of each transcript were not kept
  const txWaits = tx && !(tx.cells || []).length && (D.sectionState(b, 'rna.tx') || {}).reason === 'waiting';
  const satNote = txWaits ? 'The alignments of this run to the transcripts are linked but not read yet: coverage along the transcript, transcripts found and reads per transcript appear when a scan has read them.'
    : tx && !SC.length ? (c.sKeys.length ? 'Reads per transcript and the saturation curve are kept per sample for runs of up to twelve samples.' : 'Reads per transcript and the saturation curve are not shown: at least one file was read without the count of each transcript.') : null;
  const reading = (TT && TT.reads ? F.pct1(share(TT.mapped, TT.reads)) + ' of ' + ofReads((tx || {}).population) + (c.sKeys.length ? ' of ' + chosenWords(c) : '') + ' align to a transcript' + (whole && S.rna_tx_transcripts_detected != null ? '; ' + F.int(S.rna_tx_transcripts_detected) + ' transcripts' + (S.rna_tx_genes_detected != null ? ' of ' + F.int(S.rna_tx_genes_detected) + ' genes' : '') + ' have at least one read' : '') + '.' : '') + (tail ? ' The typical poly(A) tail is ' + Math.round(tail.med) + ' nucleotides.' : '');
  if (c.sKeys.length && !tail && !body && !sat && !rpt && !(TT && TT.reads) && !satNote) return <Section id="transcripts" title="Transcripts" x={{ evRule, evSec }}>{p ? <p class="fc-section__note fc-notice fc-partial" role="status">{p.text}</p> : null}<NoneFor c={c} what="transcript numbers" /></Section>;
  return <Section id="transcripts" title="Transcripts" reading={reading.trim() || null} x={{ evRule, evSec }}>
    {p ? <p class="fc-section__note fc-notice fc-partial" role="status">{p.text}</p> : null}
    {variantNote(b, ['rna.tx', 'rna.polya', 'rna.coverage']) ? <p class="fc-section__note fc-notice fc-variant" role="status">{variantNote(b, ['rna.tx', 'rna.polya', 'rna.coverage'])}</p> : null}
    <div class="fc-tiles fc-section__tiles">
      {TT && TT.reads ? <Tile label="Reads on a transcript" value={F.pctNum(share(TT.mapped, TT.reads))} unit="%" caption={cap(F.count(TT.mapped) + ' of ' + F.count(TT.reads) + ' reads')} title="rna_tx_aligned_frac" /> : null}
      {whole && S.rna_tx_transcripts_detected != null ? <Tile label="Transcripts found" value={F.countParts(S.rna_tx_transcripts_detected)[0]} unit={F.countParts(S.rna_tx_transcripts_detected)[1]} caption={cap('with at least one read')} title="rna_tx_transcripts_detected" /> : null}
      {whole && S.rna_tx_genes_detected != null ? <Tile label="Genes found" value={F.countParts(S.rna_tx_genes_detected)[0]} unit={F.countParts(S.rna_tx_genes_detected)[1]} caption={cap('with at least one read')} title="rna_tx_genes_detected" /> : null}
      {TT && TT.mapped ? <Tile label="On the opposite strand" value={F.pctNum(share(TT.rev, TT.mapped))} unit="%" caption={cap('of the reads on a transcript')} title="rna_antisense_frac" /> : null}
      {tail ? <Tile label="Median poly(A) tail, reads with a tail" value={String(Math.round(tail.med))} unit="nt" caption={cap('every read in which the basecaller found a tail')} title="rna_polya_median" /> : null}
      {P && P.reads ? <Tile label="Tail found" value={F.pctNum(share(P.found, P.reads))} unit="%" caption={cap(F.count(P.found) + ' of ' + F.count(P.reads) + ' reads')} title="rna_polya_detect_frac" /> : null}
    </div>
    <Grid>
      {tail ? <ChartPanel id={pid('rna/polya')} title={panelTitle('rna/polya')} reading={tail.reading} scope={st.scope} tag={tagOf(st, partialInfo(b, 'rna.polya'))} height={220} spec={tail.spec} legend={tail.legend} legendToggle={overlay} how={how('rna_polya_hist')} table={tail.table} /> : null}
      {body ? <ChartPanel id={pid('rna/body_coverage')} title={panelTitle('rna/body_coverage')} reading={body.reading} scope={[txSet, scopeTag(c, ['samples']).scope].filter(Boolean).join(' · ')} tag={tagOf(scopeTag(c, ['samples']), partialInfo(b, 'rna.coverage'))} height={220} spec={body.spec} legend={body.legend} how={how('rna_body_coverage')} table={body.table} /> : null}
      {sat ? <ChartPanel id={pid('rna/saturation')} title={panelTitle('rna/saturation')} reading={sat.reading} scope={st.scope} tag={tag} height={220} spec={sat.spec} legend={sat.legend} how={how('rna_saturation_curve')} table={sat.table} /> : null}
      {rpt ? <ChartPanel id={pid('rna/reads_per_transcript')} title={panelTitle('rna/reads_per_transcript')} reading={rpt.reading} scope={st.scope} tag={tag} height={220} spec={rpt.spec} legend={rpt.legend} how={how('rna_reads_per_transcript')} table={rpt.table} /> : null}
    </Grid>
    {satNote ? <p class="fc-section__note">{satNote}</p> : null}
    <Disclosure id="n-tx" c={c} title={'All numbers (' + countNumbers(b, ['rna.tx', 'rna.coverage', 'rna.polya']) + ')'}><AllNumbers b={b} mods={['rna.tx', 'rna.coverage', 'rna.polya']} /></Disclosure>
  </Section>;
}

// =====================================================================================================
// Reads and bases per sample: the chart that stands with the Samples table
// =====================================================================================================
export function SamplesChart({ b, c }) {
  const sy = D.sampleYield(b), rows = c.named.filter(s => sy[s.key] && sy[s.key].pass).map((s, i) => ({ s, i, y: sy[s.key].pass })); if (rows.length < 2 && !(sy.unclassified && sy.unclassified.pass)) return null;
  const un = sy.unclassified && sy.unclassified.pass, all = rows.map(r => ({ key: r.s.key, name: r.s.name, reads: r.y.reads, bases: r.y.bases, color: colorOf(c, r.s.key) })).concat(un ? [{ key: 'unclassified', name: 'No barcode found', reads: un.reads, bases: un.bases, color: tk('--fc-fate-no-barcode') }] : []);
  const pop = (D.mod(b, 'len') || {}).population, set = readSetWords(pop, ['pass']), p = partialInfo(b, 'len');
  const chart = (k, word, fmt) => { const tot = all.reduce((a, r) => a + r[k], 0) || 1, named = all.filter(r => r.key !== 'unclassified'), vals = named.map(r => r[k]).sort((p1, p2) => p1 - p2), lo = vals[0], hi = vals[vals.length - 1];
    return { spec: { form: 'hbars', rowH: 24, rows: all.map(r => ({ label: r.name, short: c.short(r.key, 22), value: r[k], color: r.color, text: fmt(r[k]), note: F.pct1(r[k] / tot) })), readout: i => [[c.short(all[i].key, 30), ''], [F.int(all[i].reads), 'reads'], [F.bases(all[i].bases), ''], [F.pct1(all[i][k] / tot), 'of the ' + word]] },
      reading: named.length > 1 ? 'The sample with the most ' + word + ' has ' + (lo > 0 ? (hi / lo).toFixed(1) + ' times' : 'far more than') + ' the ' + word + ' of the one with the fewest' + (un ? '; ' + F.pct1(un[k] / tot) + ' of the ' + word + ' carry no barcode.' : '.') : F.pct1((un ? un[k] : 0) / tot) + ' of the ' + word + ' carry no barcode.' }; };
  const R = chart('reads', 'reads', F.count), B = chart('bases', 'bases', F.bases), table = () => ({ head: ['Sample', 'Reads', 'Bases'], rows: all.map(r => [r.name, F.int(r.reads), F.int(r.bases)]) }), h = all.length * 24 + 10;
  return <Grid>
    <ChartPanel id={pid('yield/per_sample')} title="Reads per sample" reading={R.reading + ' Counts are ' + set + '.'} tag={p ? p.tag : null} height={h} spec={R.spec} how={how('yield_basecalled_per_sample')} table={table} />
    <ChartPanel id="p-yield-per_sample-bases" title="Bases per sample" reading={B.reading} tag={p ? p.tag : null} height={h} spec={B.spec} how={how('yield_basecalled_per_sample')} table={table} />
  </Grid>;
}
