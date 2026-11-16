// Pychopper: everything the statistics file of the cDNA pre-processor holds (section cdna.pychopper of a bundle).
// One table per sample; every chart shows the samples one by one, or the samples chosen in the Samples control added
// up. Collected numbers: exact, and each says what it is a count of (reads, or pieces of reads).
import { Tile } from '../shell/parts.jsx';
import { ChartPanel } from '../charts/ChartPanel.jsx';
import { tk } from '../charts/theme.js';
import * as F from '../core/format.js';
import * as D from '../core/derive.js';
import { Section, Disclosure, AllNumbers, countNumbers, ruleState } from './section.jsx';
import { scopeTag, sharedTables, sectionResult, ruleById } from './common.js';
import { sampleColor } from './readdata.jsx';
import { Grid } from '../charts/Grid.jsx';

export const PYCHOPPER = 'cdna.pychopper';
const ADD = ['input', 'pass_reads', 'qc_fail', 'len_fail', 'primers_found', 'rescued_reads', 'rescue_segments', 'unusable', 'oriented', 'written', 'strand_plus', 'strand_minus',
  'rescue_strand_plus', 'rescue_strand_minus', 'unusable_no_hit', 'unusable_multi_hit', 'hits_reads', 'primer_config_other', 'trimmed_reads', 'trimmed_bases', 'umi_detected', 'umi_detected_final'];
const HISTS = ['cdna_rescue_segment_hist', 'cdna_unusable_hits_hist', 'cdna_rescue_hits_hist', 'cdna_trimmed_bases_hist'];
const share = (a, b) => F.isNum(a) && b > 0 ? a / b : null;
// a share in words: one decimal, two for a share under one in a thousand (so that a few reads never read as "0.0 %")
const pc = f => !F.isNum(f) ? F.DASH : f > 0 && f < 0.00005 ? 'under 0.01 %' : f > 0 && f < 0.001 ? (f * 100).toFixed(2) + ' %' : F.pct1(f);
const pcNum = f => !F.isNum(f) ? F.DASH : f > 0 && f < 0.001 ? (f * 100).toFixed(2) : F.pctNum(f);
const centers = e => Array.from({ length: e.length - 1 }, (_, i) => (e[i] + e[i + 1]) / 2);

// a histogram of a cell as a full array on its axis, with what lies beyond the axis
function hist(b, h) {
  const edges = D.edgesOf(b, h.axis), n = D.expand(h, edges.length - 1), w = h.w ? D.expand(h, edges.length - 1, 'w') : null;
  return { edges, n, w, over: h.over || 0, wOver: h.w_over || 0, min: h.min, max: h.max };
}
// the chosen samples added up: every number that is a sum, the named primer arrangements, the four histograms.
// A number one of the samples does not have (UMI detection was off) is not known for the sum either.
function pool(b, cells) {
  const sc = {}, unknown = new Set();
  for (const k of ADD) sc[k] = 0;
  for (const c of cells) for (const k of ADD) { const v = c.scalars[k]; if (v == null) unknown.add(k); else sc[k] += v; }
  for (const k of unknown) sc[k] = null;
  const cfg = cells[0].counts.cdna_primer_config_table, cats = cfg.cats, n = cats.map(() => 0);
  for (const c of cells) c.counts.cdna_primer_config_table.n.forEach((v, i) => { n[i] += v; });
  const H = {};
  for (const id of HISTS) { let t = null; for (const c of cells) { const h = c.hists && c.hists[id]; if (!h) continue; const x = hist(b, h);
    if (!t) t = x; else { D.addInto(t.n, x.n); if (t.w && x.w) D.addInto(t.w, x.w); t.over += x.over; t.wOver += x.wOver; t.min = Math.min(t.min, x.min); t.max = Math.max(t.max, x.max); } } H[id] = t; }
  return { sc, cats, n, H, kinds: Math.max(...cells.map(c => c.scalars.primer_configs || 0)) };
}
// the cut-off of one sample: the point of its tuning curve that is marked as chosen
const cutoff = c => { const v = c.curves && c.curves.cdna_autotune_curve; if (!v) return null; const i = v.y.chosen ? v.y.chosen.indexOf(1) : -1; return i < 0 ? null : v.x[i]; };
const fmtQ = q => q == null ? 'not known' : q !== 0 && Math.abs(q) < 0.001 ? String(+q.toPrecision(2)) : String(+q.toFixed(3));

// one name per class, on every chart, tile and table of the page
const CLASSES = [['primers_found', 'Full-length, one cDNA', '--fc-cdna-full-length', 'one cDNA with both primers'], ['rescued_reads', 'Full-length, rescued', '--fc-cdna-rescued', 'two or more cDNAs joined together, cut apart'],
  ['unusable', 'Unusable', '--fc-cdna-unusable', 'no usable primer arrangement'], ['qc_fail', 'Dropped: low quality', '--fc-cdna-quality-fail', 'not classified']];
const HITS = [['No primer', '--fc-ord-1', h => h.n[0]], ['One primer', '--fc-ord-2', h => h.n[1]], ['Two primers', '--fc-ord-3', h => h.n[2]], ['Three primers', '--fc-ord-4', h => h.n[3]],
  ['Four or more', '--fc-ord-5', h => { let t = h.over; for (let i = 4; i < h.n.length; i++) t += h.n[i]; return t; }]];
const UMI = [['UMI found and kept', '--fc-cat-1'], ['UMI found, not kept', '--fc-cat-3'], ['No UMI found', '--fc-chart-none', 'none']];

export function PychopperSection({ b, c, evRule, evSec }) {
  const m = D.mod(b, PYCHOPPER); if (!m || !m.cells || !m.cells.length) return null;
  const by = {}; for (const cell of m.cells) by[cell.at[0]] = cell;
  const all = c.named.filter(s => by[s.key]), names = {}; all.forEach(s => { names[s.key] = s.name; });
  if (!all.length) return null;
  const chosen = all.filter(s => c.sKeys.includes(s.key)), sel = chosen.length ? chosen : all, cells = sel.map(s => by[s.key]);
  const P = pool(b, cells), sc = P.sc, one = sel.length === 1, whom = one ? sel[0].name : chosen.length ? F.plural(sel.length, 'sample') + ' chosen' : all.length === 1 ? all[0].name : 'all ' + all.length + ' samples';
  const st = scopeTag(c, ['samples']), missing = m.state === 'partial' ? (m.notes || []).filter(n => n.level === 'warning') : [];
  // a table that covers the sample on two flow cells: said at the top, in the sentence, under every count and on every row
  const shared = sharedTables(b, [PYCHOPPER], c, 'pychopper'), sharedSel = !!shared && sel.some(q => shared.keys.has(q.key)), capS = t => sharedSel ? t + ' · ' + shared.caption : t;
  const mark = k => shared ? shared.mark(k) : '', mid = (name, n = 26) => c.shorten(name, n), rule = ruleById(b, 'cdna_full_length_frac'), whole = sel.length === all.length;
  const tag = [st.tag, m.state === 'partial' ? 'PARTIAL · ' + all.length + ' of ' + c.named.length + ' samples' : null].filter(Boolean).join(' · ');
  const strong = k => chosen.length > 0 && chosen.length < all.length && c.sKeys.includes(k);
  // rows of the per-sample charts: the samples chosen added up first (when there is more than one), then every sample
  const rowsOf = partsOf => (sel.length > 1 ? [{ label: chosen.length ? sel.length + ' chosen' : 'All samples', strong: true, parts: partsOf(sc, P), sum: true }] : [])
    .concat(all.map(s => ({ label: s.name + mark(s.key), short: c.short(s.key, 22) + mark(s.key), strong: strong(s.key), parts: partsOf(by[s.key].scalars, null, by[s.key]), key: s.key })));
  const readStack = (rows, names2, unit, of) => k => { const r = rows[Math.floor(k / 1000)], i = k % 1000, t = r.parts.reduce((a, q) => a + q.value, 0);
    return [[r.short || r.label, ''], [pc(r.parts[i].value / (t || 1)), F.lower(names2[i])], [F.int(r.parts[i].value), unit], 'of ' + F.count(t) + ' ' + of]; };
  const height = rows => rows.length * 22 + 34;

  // ---- 1. what pychopper made of the reads ----
  const clsRows = rowsOf(s => CLASSES.map(([k, , token]) => ({ value: s[k] || 0, color: tk(token) })));
  const fl = share(sc.primers_found + sc.rescued_reads, sc.pass_reads), flBy = all.map(s => share(by[s.key].scalars.primers_found + by[s.key].scalars.rescued_reads, by[s.key].scalars.pass_reads));
  const range = a => { const v = a.filter(F.isNum); return v.length > 1 ? ' (' + F.pct1(Math.min(...v)) + ' to ' + F.pct1(Math.max(...v)) + ' across ' + (chosen.length ? 'all ' + all.length : 'the') + ' samples)' : ''; };
  // ---- 2. strand ----
  const strandRows = rowsOf(s => [{ value: s.strand_plus || 0, color: tk('--fc-cat-1') }, { value: s.strand_minus || 0, color: tk('--fc-cat-2') }]);
  const plus = share(sc.strand_plus, sc.strand_plus + sc.strand_minus), plusBy = all.map(s => share(by[s.key].scalars.strand_plus, by[s.key].scalars.strand_plus + by[s.key].scalars.strand_minus));
  // ---- 3. rescue: pieces per read that was cut apart ----
  const seg = P.H.cdna_rescue_segment_hist, segN = seg ? Array.from(seg.n) : [], segMax = seg ? Math.max(3, Math.min(seg.n.length - 1, seg.over ? seg.n.length - 1 : (F.isNum(seg.max) ? seg.max : 3))) : 3;
  const two = seg ? share(seg.n[2], sc.rescued_reads) : null;
  const countSpec = (h, lo, hi, xt, yt, unit) => ({ form: 'xy', x: { domain: [lo, hi + 1], kind: 'linear', title: xt, ticks: Array.from({ length: hi - lo + 1 }, (_, i) => lo + i).filter(v => hi - lo < 12 || v % 2 === 0).map(v => ({ v: v + 0.5, label: String(v), major: true })) },
    y: c.logY ? { log: true, title: yt + ', logarithmic' } : { title: yt }, series: [{ name: unit, type: 'hist', edges: h.edges, y: Array.from(h.n) }], cursor: { x: centers(h.edges).filter((_, i) => i >= lo && i <= hi) },
    readout: i => [[String(lo + i), xt.toLowerCase().replace(/ per .*$| on .*$/, '')], [F.int(h.n[lo + i]), unit]] });
  const countTable = (h, head, unit) => () => ({ head: [head, unit], rows: Array.from(h.n, (v, i) => [String(i), F.int(v)]).filter((_, i) => h.n[i] > 0).concat(h.over ? [[(h.n.length) + ' or more', F.int(h.over)]] : []) });
  // ---- 4. primer hits ----
  const hitRows = rowsOf((s, pooled, cell) => { const h = pooled ? pooled.H.cdna_unusable_hits_hist : hist(b, cell.hists.cdna_unusable_hits_hist); return HITS.map(([, token, f]) => ({ value: h ? f(h) : 0, color: tk(token) })); });
  const noHit = share(sc.unusable_no_hit, sc.unusable), multi = share(sc.unusable_multi_hit, sc.unusable), rh = P.H.cdna_rescue_hits_hist;
  const rhMax = rh ? Math.max(5, Math.min(rh.n.length - 1, F.isNum(rh.max) ? rh.max : 5)) : 5;
  // ---- 5. primer arrangements ----
  const cfgTotal = sc.hits_reads || 1, cfgRows = P.cats.map((name, i) => ({ label: name.replace(/,/g, ', '), value: P.n[i], text: F.count(P.n[i]), note: pc(P.n[i] / cfgTotal) })).filter(r => r.value > 0)
    .concat(sc.primer_config_other > 0 ? [{ label: 'All other arrangements', value: sc.primer_config_other, text: F.count(sc.primer_config_other), note: pc(sc.primer_config_other / cfgTotal), color: tk('--fc-chart-none') }] : []);
  const topTwo = P.n.slice().sort((a, x) => x - a).slice(0, 2).reduce((a, v) => a + v, 0) / cfgTotal;
  // ---- 6. bases cut off ----
  const tr = P.H.cdna_trimmed_bases_hist, trMed = tr ? D.quantile(tr.edges, tr.n, Math.min(1, 0.5 * (D.sum(tr.n) + tr.over) / (D.sum(tr.n) || 1))) : null, trMean = share(sc.trimmed_bases, sc.trimmed_reads),
    trExact = sel.length === m.cells.length && F.isNum((m.summary || {}).cdna_trimmed_bases_median) ? m.summary.cdna_trimmed_bases_median : null, trHalf = trExact != null ? trExact : trMed != null ? Math.round(trMed) : null;
  const overlay = chosen.length > 1 && chosen.length <= 8;
  let trSpec = null, trLegend = null, trTable = null;
  if (tr) {
    const dom = D.spanQ(tr.edges, tr.n, 0.002, 4);
    if (overlay) { const ser = chosen.map(s => { const h = hist(b, by[s.key].hists.cdna_trimmed_bases_hist), t = (D.sum(h.n) + h.over) || 1, i = c.named.findIndex(q => q.key === s.key); return { name: s.name, type: 'outline', edges: h.edges, y: Array.from(h.n, v => v / t * 100), color: sampleColor(s, i) }; });
      trSpec = { form: 'xy', x: { domain: dom, kind: 'linear', title: 'Bases cut off a read' }, y: { title: 'Share of each sample’s usable reads per bin, %' }, series: ser, cursor: { x: centers(tr.edges) }, readout: i => [[tr.edges[i] + ' to ' + (tr.edges[i + 1] - 1) + ' bases', '']].concat(ser.map(q => [q.y[i].toFixed(2) + ' %', mid(q.name, 18)])) };
      trLegend = ser.map(q => ({ name: q.name, label: q.name, kind: 'line', color: q.color }));
      trTable = () => ({ head: ['Bases cut off'].concat(ser.map(q => q.name + ', %')), rows: Array.from(tr.n, (_, i) => [tr.edges[i] + ' to ' + (tr.edges[i + 1] - 1)].concat(ser.map(q => q.y[i].toFixed(2)))).filter((_, i) => tr.n[i] > 0) });
    } else { const u = Math.max(...tr.n) >= 1e6 ? 1e6 : 1;
      trSpec = { form: 'xy', x: { domain: dom, kind: 'linear', title: 'Bases cut off a read' }, y: c.logY ? { log: true, title: 'Usable reads per bin' + (u === 1e6 ? ', millions' : '') + ', logarithmic' } : { title: 'Usable reads per bin' + (u === 1e6 ? ', millions' : '') },
        series: [{ name: 'Usable reads', type: 'hist', edges: tr.edges, y: Array.from(tr.n, v => v / u) }], vlines: trHalf != null ? [{ x: trHalf + 0.5, label: 'median ' + trHalf + ' b', strong: true }] : [], cursor: { x: centers(tr.edges) },
        readout: i => [[tr.edges[i] + ' to ' + (tr.edges[i + 1] - 1) + ' bases', ''], [F.int(tr.n[i]), 'usable reads'], [pc(tr.n[i] / (sc.trimmed_reads || 1)), 'of them']] };
      trTable = () => ({ head: ['Bases cut off', 'Usable reads'], rows: Array.from(tr.n, (v, i) => [tr.edges[i] + ' to ' + (tr.edges[i + 1] - 1), F.int(v)]).filter((_, i) => tr.n[i] > 0).concat(tr.over ? [[tr.edges[tr.edges.length - 1] + ' or more', F.int(tr.over)]] : []) }); }
  }
  // ---- 7. dropped, and UMI ----
  const dropRows = rowsOf(s => [{ value: s.qc_fail || 0, color: tk('--fc-cdna-quality-fail') }, { value: s.len_fail || 0, color: tk('--fc-cdna-length-fail') }]).filter(r => !r.sum);
  const umiKnown = all.filter(s => by[s.key].scalars.umi_detected != null), umiParts = s => { const d = s.umi_detected || 0, f = s.umi_detected_final == null ? d : s.umi_detected_final; return [f, Math.max(0, d - f), Math.max(0, (s.oriented || 0) - d)]; };
  const umiRows = umiKnown.length ? (sc.umi_detected != null && sel.length > 1 ? [{ label: chosen.length ? sel.length + ' chosen' : 'All samples', strong: true, parts: umiParts(sc).map((v, i) => ({ value: v, color: tk(UMI[i][1]) })) }] : [])
    .concat(umiKnown.map(s => ({ label: s.name + mark(s.key), short: c.short(s.key, 22) + mark(s.key), strong: strong(s.key), parts: umiParts(by[s.key].scalars).map((v, i) => ({ value: v, color: tk(UMI[i][1]) })) }))) : [];
  const umi = share(sc.umi_detected, sc.oriented), umiFinal = share(sc.umi_detected_final, sc.oriented);
  // ---- 8. the cut-off pychopper tuned itself to ----
  const curves = all.map((s, i) => ({ s, i: c.named.findIndex(q => q.key === s.key), v: by[s.key].curves && by[s.key].curves.cdna_autotune_curve, q: cutoff(by[s.key]) })).filter(o => o.v && o.v.x.length > 1);
  const cuts = [...new Set(all.map(s => cutoff(by[s.key])).filter(q => q != null))].sort((a, x) => a - x), cutSel = [...new Set(sel.map(s => cutoff(by[s.key])).filter(q => q != null))];
  let tuneSpec = null, tuneLegend = null, tuneTable = null;
  if (curves.length) {
    const few = curves.length <= 8, some = !few && curves.some(o => c.sKeys.includes(o.s.key)), lit = o => few || c.sKeys.includes(o.s.key), x = curves[0].v.x;
    const ser = curves.map(o => { const top = Math.max(...o.v.y.reads.filter(F.isNum), 1); return { name: o.s.name, type: 'line', x: o.v.x, y: o.v.y.reads.map(v => F.isNum(v) ? v / top * 100 : NaN), raw: o.v.y.reads, lit: lit(o),
      color: lit(o) ? sampleColor(o.s, o.i) : some ? tk('--fc-chart-none-line') : tk('--fc-this-run'), width: lit(o) ? 1.75 : 1, alpha: few || some ? 1 : 0.45 }; })
      .sort((a, q) => (a.lit > q.lit) - (a.lit < q.lit));
    const low = Math.min(...ser.flatMap(q => q.y.filter(F.isNum))), tally = new Map(); for (const o of curves) if (o.q != null) tally.set(o.q, (tally.get(o.q) || 0) + 1);
    const most = [...tally.entries()].sort((a, q) => q[1] - a[1] || a[0] - q[0]).map(e => e[0]);
    tuneSpec = { form: 'xy', x: { domain: [Math.min(...x), Math.max(...x)], kind: 'linear', title: 'Cut-off tried (q)' }, y: { domain: [Math.max(0, Math.floor(low / 5) * 5), 100], title: 'Reads classified, % of the sample’s best' }, series: ser,
      vlines: most.slice(0, 4).map((q, k) => ({ x: q, label: k ? '' : 'chosen: ' + fmtQ(q) + (most.length > 1 ? ' (' + tally.get(q) + ' of ' + curves.length + ' samples)' : ''), strong: k === 0 })), cursor: { x },
      readout: i => { const list = few ? ser : ser.filter(q => q.lit), shown = list.slice(0, 4);      // the unit of the axis; what does not fit is in the table, and said
        return [['q ' + fmtQ(x[i]), 'reads classified, % of the sample’s best']].concat(shown.map(q => [F.isNum(q.y[i]) ? q.y[i].toFixed(1) + ' %' : F.DASH, c.shorten(q.name, 14)])).concat(list.length > shown.length ? ['+ ' + (list.length - shown.length) + ' more in the table'] : []).concat(few ? [] : [ser.length + ' samples drawn, one line each']); } };
    tuneLegend = few ? ser.map(q => ({ name: q.name, label: q.name, kind: 'line', color: q.color })) : some ? ser.filter(q => q.lit).slice(0, 8).map(q => ({ name: q.name, label: q.name, kind: 'line', color: q.color })).concat([{ name: 'other', label: 'Other samples', kind: 'line', color: tk('--fc-chart-none-line') }]) : null;
    tuneTable = () => ({ head: ['Cut-off (q)'].concat(curves.map(o => o.s.name)), rows: x.map((q, i) => [fmtQ(q)].concat(curves.map(o => F.isNum(o.v.y.reads[i]) ? F.int(o.v.y.reads[i]) + (o.q != null && Math.abs(o.q - q) < 1e-12 ? ' (chosen)' : '') : F.DASH))) });
  }
  // one definition per number: the three classes are shares of the reads pychopper classified (the full-length share
  // with the rescued reads is the number the check judges); what was dropped before classifying is a share of the reads given
  const reading = 'Pychopper was given ' + F.count(sc.input) + ' reads' + (one ? ' of ' + sel[0].name : all.length > 1 ? ' (' + whom + ')' : '') + (sharedSel ? ', merged from ' + shared.both + ' of the sample,' : '')
    + ' and classified ' + F.count(sc.pass_reads) + ' of them. Of the reads classified, ' + pc(fl) + ' are full-length cDNA (' + pc(share(sc.primers_found, sc.pass_reads)) + ' as one cDNA with both primers, ' + pc(share(sc.rescued_reads, sc.pass_reads))
    + ' rescued by cutting apart a read that held several) and ' + pc(share(sc.unusable, sc.pass_reads)) + ' are unusable. ' + pc(share(sc.qc_fail, sc.input)) + ' of the reads given were dropped for low quality before classifying.'
    + (rule && rule.outcome === 'not_evaluated' ? ' ' + rule.message : '');
  const sampleTable = cols => () => ({ head: ['Sample'].concat(cols.map(x => x[0])), rows: (sel.length > 1 ? [[chosen.length ? 'Samples chosen' : 'All samples'].concat(cols.map(x => x[1](sc, P)))] : []).concat(all.map(s => [s.name + mark(s.key)].concat(cols.map(x => x[1](by[s.key].scalars, null, by[s.key]))))),
    note: shared ? shared.legend : undefined });
  const hOf = (pooled, cell, id) => pooled ? pooled.H[id] : hist(b, cell.hists[id]);
  // the check on the full-length share judges all samples together: its mark stands on the tile while all are shown
  const ruleSt = rule && whole && rule.outcome !== 'not_evaluated' ? ruleState(rule) : null;
  return <Section id="pychopper" title="Pychopper" result={sectionResult(b, 'cdna')} reading={reading} x={{ evRule, evSec }}>
    <span id="s-cdna" />
    {shared ? <p class="fc-section__note fc-notice fc-shared" role="status">{shared.text}</p> : null}
    {(m.notes || []).filter(n => n.code === 'pychopper_categories_absent').map(n => <p class="fc-section__note fc-notice fc-shared" role="status">{n.message}</p>)}
    {m.state === 'partial' ? <p class="fc-section__note fc-notice fc-partial" role="status">Partial: {all.length} of {c.named.length} samples have a pychopper statistics file linked to this run. {missing.length ? missing[0].message : 'The numbers here cover those samples only.'}</p> : null}
    <div class="fc-tiles fc-section__tiles">
      <Tile label="Full-length reads" glyph={ruleSt} state={ruleSt === 'warn' || ruleSt === 'fail' ? ruleSt : null} value={pc(fl).replace(/\s*%$/, '')} unit="%" caption={capS(F.count(sc.primers_found + sc.rescued_reads) + ' of ' + F.count(sc.pass_reads) + ' reads classified, with the rescued ones')} title="cdna_full_length_frac" />
      <Tile label="Full-length, rescued" value={pcNum(share(sc.rescued_reads, sc.pass_reads))} unit="%" caption={capS(F.int(sc.rescued_reads) + ' of ' + F.count(sc.pass_reads) + ' reads classified, cut into ' + F.int(sc.rescue_segments) + ' pieces')} title="cdna_rescued_frac" />
      <Tile label="Unusable" value={pcNum(share(sc.unusable, sc.pass_reads))} unit="%" caption={capS(F.count(sc.unusable) + ' of ' + F.count(sc.pass_reads) + ' reads classified; ' + F.pct1(noHit) + ' of them without any primer')} title="cdna_unusable_frac" />
      <Tile label="Dropped: low quality" value={pcNum(share(sc.qc_fail, sc.input))} unit="%" caption={capS(F.int(sc.qc_fail) + ' of ' + F.count(sc.input) + ' reads given, before classifying')} title="cdna_qc_fail_frac" />
      <Tile label="Dropped: too short" value={pcNum(share(sc.len_fail, sc.oriented))} unit="%" caption={capS(F.int(sc.len_fail) + ' of ' + F.count(sc.oriented) + ' pieces kept, after classifying')} title="cdna_len_fail_frac" />
      {sc.umi_detected != null ? <Tile label="With a UMI" value={F.pctNum(umi)} unit="%" caption={capS(F.count(sc.umi_detected) + ' of ' + F.count(sc.oriented) + ' pieces kept' + (umiFinal != null ? '; ' + F.pct1(umiFinal) + ' after pychopper’s own checks' : ''))} title="cdna_umi_detected_frac" />
        : <Tile label="With a UMI" state="unknown" value="not known" caption={umiKnown.length ? 'UMI detection was off for some of these samples' : 'UMI detection was off'} />}
      <Tile label="Reads written" value={F.countParts(sc.written)[0]} unit={F.countParts(sc.written)[1]} caption={capS('full-length reads and pieces that were long enough')} />
      <Tile label="Cut-off chosen" value={cutSel.length === 1 ? fmtQ(cutSel[0]) : cutSel.length ? cutSel.length + ' values' : F.DASH} caption={cutSel.length > 1 ? 'the samples tuned to different cut-offs' : 'pychopper tuned it on a sample of the reads'} title="cdna_pychopper_cutoff" />
    </div>
    {/* one bar per sample: side by side, as many as the window holds */}
    <Grid>
      <ChartPanel id="p-pych-class" title="What pychopper made of the reads" reading={F.pct1(fl) + ' of the reads classified are full-length (one cDNA, or rescued)' + range(flBy) + '. A bar is every read given to pychopper.'} scope={st.scope} tag={tag} height={height(clsRows)}
        spec={{ form: 'hstack', norm: true, rowH: 22, rows: clsRows, readout: readStack(clsRows, CLASSES.map(x => x[1]), 'reads', 'reads given to pychopper') }} legend={CLASSES.map(([, name, token]) => ({ name, color: tk(token) }))}
        how={'Each bar is one sample’s reads given to pychopper. Full-length, one cDNA: both primers, one molecule. Full-length, rescued: a read made of two or more cDNAs, cut apart. Unusable: no usable primer arrangement. Dropped: below pychopper’s quality limit, never classified.' + (shared ? ' ' + shared.legend + '.' : '')}
        table={sampleTable([['Given to pychopper', s => F.int(s.input)], ['Classified', s => F.int(s.pass_reads)], ['Full-length, one cDNA', s => F.int(s.primers_found)], ['Full-length, rescued', s => F.int(s.rescued_reads)], ['Unusable', s => F.int(s.unusable)], ['Dropped: low quality', s => F.int(s.qc_fail)], ['Full-length share of the reads classified', s => F.pct1(share(s.primers_found + s.rescued_reads, s.pass_reads))]])} />
      <ChartPanel id="p-pych-strand" title="Strand of the full-length reads" reading={F.pct1(plus) + ' of the full-length reads are on the plus strand' + range(plusBy) + '; about half is expected.'} scope={st.scope} tag={tag} height={height(strandRows)}
        spec={{ form: 'hstack', norm: true, rowH: 22, rows: strandRows, readout: readStack(strandRows, ['Plus strand', 'Minus strand'], 'reads', 'full-length reads') }} legend={[{ name: 'Plus strand', color: tk('--fc-cat-1') }, { name: 'Minus strand', color: tk('--fc-cat-2') }]}
        table={sampleTable([['Full-length, plus', s => F.int(s.strand_plus)], ['Full-length, minus', s => F.int(s.strand_minus)], ['Plus share', s => F.pct1(share(s.strand_plus, s.strand_plus + s.strand_minus))], ['Rescued pieces, plus', s => F.int(s.rescue_strand_plus)], ['Rescued pieces, minus', s => F.int(s.rescue_strand_minus)]])} />
      <ChartPanel id="p-pych-hits" title="Primers found on the unusable reads" reading={F.pct1(noHit) + ' of the unusable reads carry no primer at all, ' + F.pct1(multi) + ' three or more.'} scope={st.scope} tag={tag} height={height(hitRows)}
        spec={{ form: 'hstack', norm: true, rowH: 22, rows: hitRows, readout: readStack(hitRows, HITS.map(x => x[0]), 'reads', 'unusable reads') }} legend={HITS.map(([name, token]) => ({ name, color: tk(token) }))}
        how="No primer at all usually means the read is not cDNA, or is badly damaged. Three or more primers mean several cDNAs joined together in a way pychopper could not cut apart."
        table={sampleTable([['Unusable reads', s => F.int(s.unusable)]].concat(HITS.map(([name, , f]) => [name, (s, pooled, cell) => F.int(f(hOf(pooled, cell, 'cdna_unusable_hits_hist')))])))} />
      <ChartPanel id="p-pych-dropped" title="Dropped by pychopper" reading={F.int(sc.qc_fail) + ' reads were dropped for low quality before classifying, ' + F.int(sc.len_fail) + ' pieces for being too short after it.'} scope={st.scope} tag={tag} height={height(dropRows)}
        spec={{ form: 'hstack', rowH: 22, rows: dropRows, fmt: v => v >= 1e6 ? +(v / 1e6).toPrecision(3) + ' M' : F.int(v), readout: k => { const r = dropRows[Math.floor(k / 1000)], i = k % 1000; return [[r.short || r.label, ''], [F.int(r.parts[i].value), i ? 'pieces too short' : 'reads of low quality']]; } }}
        legend={[{ name: 'Dropped: low quality (reads)', color: tk('--fc-cdna-quality-fail') }, { name: 'Dropped: too short (pieces)', color: tk('--fc-cdna-length-fail') }]}
        table={sampleTable([['Dropped: low quality, reads', s => F.int(s.qc_fail)], ['Share of reads given', s => pc(share(s.qc_fail, s.input))], ['Dropped: too short, pieces', s => F.int(s.len_fail)], ['Share of pieces', s => pc(share(s.len_fail, s.oriented))], ['Reads written', s => F.int(s.written)]])} />
      {umiRows.length ? <ChartPanel id="p-pych-umi" title="UMI found, per sample" reading={sc.umi_detected != null ? 'A UMI was found in ' + F.pct1(umi) + ' of the pieces kept' + (umiFinal != null ? '; ' + F.pct1(umiFinal) + ' still count after pychopper’s own checks.' : '.') : 'UMI detection was on for ' + umiKnown.length + ' of ' + all.length + ' samples.'} scope={st.scope} tag={tag} height={height(umiRows)}
        spec={{ form: 'hstack', norm: true, rowH: 22, rows: umiRows, readout: readStack(umiRows, UMI.map(x => x[0]), 'pieces', 'pieces kept') }} legend={UMI.map(([name, token, kind]) => ({ name, color: tk(token), kind }))}
        table={() => ({ head: ['Sample', 'Pieces kept', 'UMI found', 'UMI found and kept', 'Share with a UMI'], rows: umiKnown.map(s => { const q = by[s.key].scalars; return [s.name + mark(s.key), F.int(q.oriented), F.int(q.umi_detected), F.int(q.umi_detected_final), F.pct1(share(q.umi_detected, q.oriented))]; }), note: shared ? shared.legend : undefined })} /> : null}
    </Grid>
    {/* the samples chosen, added up */}
    <Grid>
      {seg && sc.rescued_reads > 0 ? <ChartPanel id="p-pych-pieces" title="Pieces cut out of each rescued read" reading={F.int(sc.rescued_reads) + ' reads were cut into ' + F.int(sc.rescue_segments) + ' pieces; ' + F.pct1(two) + ' of them into two.'} scope={scopeTag(c, ['samples', 'y']).scope} tag={tag} height={200}
        spec={countSpec(seg, 2, segMax, 'Pieces per read', 'Reads', 'reads')} table={countTable(seg, 'Pieces', 'Reads')} /> : null}
      {rh && sc.rescue_segments > 0 ? <ChartPanel id="p-pych-rescuehits" title="Primers found on the rescued reads" reading={'Counted once per piece: a read cut into two pieces with four primers adds two at 4.'} scope={scopeTag(c, ['samples', 'y']).scope} tag={tag} height={200}
        spec={countSpec(rh, 2, rhMax, 'Primers on the read', 'Pieces', 'pieces')} table={countTable(rh, 'Primers on the read', 'Pieces')} /> : null}
      <ChartPanel id="p-pych-config" title="Primer arrangements found" reading={'The two most common arrangements hold ' + F.pct1(topTwo) + ' of the reads with a primer; up to ' + F.int(P.kinds) + ' different arrangements were seen in one sample.'} scope={st.scope} tag={tag} height={cfgRows.length * 24 + 10} lines={1}
        spec={{ form: 'hbars', rowH: 24, rows: cfgRows, readout: i => [[cfgRows[i].label, ''], [F.int(cfgRows[i].value), 'reads'], [cfgRows[i].note, 'of the reads with a primer']] }}
        how="The primers found on a read, in order; a minus sign means the primer was found on the other strand. A full-length cDNA has one primer at each end, in opposite directions."
        table={() => ({ head: ['Arrangement', 'Reads', 'Share'].concat(sel.length > 1 && sel.length <= 6 ? sel.map(s => s.name) : []),
          rows: P.cats.map((name, i) => [name.replace(/,/g, ', '), F.int(P.n[i]), pc(P.n[i] / cfgTotal)].concat(sel.length > 1 && sel.length <= 6 ? sel.map(s => F.int(by[s.key].counts.cdna_primer_config_table.n[i])) : []))
            .concat([['All other arrangements', F.int(sc.primer_config_other), pc(sc.primer_config_other / cfgTotal)].concat(sel.length > 1 && sel.length <= 6 ? sel.map(s => F.int(by[s.key].scalars.primer_config_other)) : [])]) })} />
      {trSpec ? <ChartPanel id="p-pych-trimmed" title="Bases cut off each usable read" reading={overlay ? 'Each sample is drawn as an outline scaled to its own total, so shapes compare.' : 'Half of the usable reads lost ' + (trHalf != null ? trHalf : F.DASH) + ' bases or fewer (mean ' + (trMean != null ? trMean.toFixed(1) : F.DASH) + '): the two primers and the adapter.'}
        scope={scopeTag(c, ['samples', 'y']).scope} tag={tag} height={220} spec={trSpec} legend={trLegend} legendToggle={overlay} table={trTable} /> : null}
      {tuneSpec ? <ChartPanel id="p-pych-tune" title="The cut-off pychopper tuned itself to" reading={cuts.length === 1 ? 'Every sample kept the cut-off ' + fmtQ(cuts[0]) + ', where most reads of its tuning sample were classified.' : 'The samples kept ' + cuts.length + ' different cut-offs (' + cuts.map(fmtQ).join(', ') + '): each kept the one where most reads of its tuning sample were classified.'}
        tag={tag} height={220} spec={tuneSpec} legend={tuneLegend} legendToggle={curves.length <= 8} table={tuneTable}
        how="Before it classifies, pychopper tries 30 cut-offs for a primer hit on a sample of the reads and keeps the one that classifies most. A curve is one sample; it is drawn relative to that sample’s best cut-off." /> : null}
    </Grid>
    <Disclosure id="n-pychopper" c={c} title={'All numbers (' + countNumbers(b, [PYCHOPPER]) + ')'}><AllNumbers b={b} mods={[PYCHOPPER]} /></Disclosure>
  </Section>;
}
