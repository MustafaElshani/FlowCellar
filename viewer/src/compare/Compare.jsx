// Compare 2 to 8 runs. Letters stay with their item and are the series
// markers; every single number is cut at the longest window every item reached; runs are aligned on hours
// since start, never on the clock. Opens with "Why the output differs": the exact factorisation
// output = working channel-hours x sequencing share x estimated bases per sequencing second.
import { index, useBundle, bundleState, rowById } from '../core/data.js';
import { useEffect } from 'preact/hooks';
import { setParams, param, go, selected, hrefRun } from '../core/state.js';
import { Badge, Mid, Seg, Crumbs, When, compareItems } from '../shell/parts.jsx';
import { ChartPanel } from '../charts/ChartPanel.jsx';
import { tk } from '../charts/theme.js';
import * as F from '../core/format.js';
import * as D from '../core/derive.js';
import { ASSAY_WORD, POPULATION_WORD } from '../core/labels.js';
import { SavePdf, PaperHead, PaperScope } from '../shell/paper.jsx';

const LETTERS = 'ABCDEFGH', DASH = [[], [7, 3], [2, 3], [9, 3, 2, 3], [12, 4], [4, 2, 1, 2], [1, 3], [14, 3, 4, 3]];
const colorOf = k => tk('--fc-cat-' + (k + 1));
const pctDiff = r => (r >= 1 ? '+' : '−') + Math.abs((r - 1) * 100).toFixed(Math.abs(r - 1) < 0.0995 ? 1 : 0) + ' %';

function words(f, ref, usual) {
  const lo = Math.log(f.bases / ref.bases), parts = [['the flow cell', Math.log(f.workHours / ref.workHours)], ['loading', Math.log(f.seqShare / ref.seqShare)], ['the library', Math.log(f.basesPerSeqSec / ref.basesPerSeqSec)]];
  if (Math.abs(lo) < usual) return ['No real difference.', 'The output differs by less than runs of this kind usually differ from each other in this lab (' + Math.round(usual * 100) + ' %).'];
  const same = parts.filter(p => Math.sign(p[1]) === Math.sign(lo)).sort((a, b) => Math.abs(b[1]) - Math.abs(a[1]));
  const dir = lo < 0 ? 'lower' : 'higher', exp = { 'the flow cell': lo < 0 ? 'fewer channels worked, or for less time' : 'more channels worked, or for longer', loading: lo < 0 ? 'working pores spent less of their time on a strand' : 'working pores spent more of their time on a strand', 'the library': lo < 0 ? 'fewer bases per second of sequencing: shorter or slower molecules' : 'more bases per second of sequencing' };
  if (!same.length) return ['Factors pull in different directions.', ''];
  if (Math.abs(same[0][1]) > Math.abs(lo) / 2 && (same.length < 2 || Math.abs(same[0][1]) > 1.6 * Math.abs(same[1][1]))) return ['Mostly ' + same[0][0] + ':', exp[same[0][0]] + '.'];
  if (same.length >= 2) return [(same[0][0][0].toUpperCase() + same[0][0].slice(1)) + ' and ' + same[1][0] + ':', exp[same[0][0]] + '; ' + exp[same[1][0]] + '.'];
  return ['Mostly ' + same[0][0] + ':', exp[same[0][0]] + '.'];
}
const Bar = ({ r }) => { const w = Math.min(60, Math.abs(Math.log(r)) * 60); return <span class={'fc-diffbar' + (r < 1 ? ' is-neg' : '')} aria-hidden="true"><i style={{ width: Math.max(1, w) + 'px' }} /></span>; };

export function Compare() {
  // the runs compared are the runs the address names: known ones, each once, eight at most; what was left out is said
  const ix = index.value, asked = compareItems(param('items', '')), ids = asked.ids;
  const rows = ids.map(rowById), bundles = ids.map(id => useBundle(id)), ready = bundles.every(Boolean);
  // the ticked runs follow this page, so that the Runs table and the application bar agree with what is compared
  useEffect(() => { if (ids.length && selected.value.join(',') !== ids.join(',')) selected.value = ids; }, [ids.join(',')]);
  const showAs = param('as', 'measured'), by = param('by', 'bases'), only = param('all', '') !== '1';
  const remove = id => { const n = ids.filter(x => x !== id); selected.value = n; if (n.length) setParams({ items: n.join(',') }); else go('/runs'); };
  const dropped = [asked.unknown.length ? asked.unknown.join(', ') + (asked.unknown.length === 1 ? ' is not a run of this library and was left out.' : ' are not runs of this library and were left out.') : null,
    asked.over.length ? 'Compare shows up to 8 runs: ' + asked.over.join(', ') + (asked.over.length === 1 ? ' was left out.' : ' were left out.') : null].filter(Boolean);
  const crumbs = <><Crumbs><a href="#/home">Library</a> {'›'} <a href="#/runs">Runs</a> {'›'} <b>Compare {F.plural(ids.length, 'run')}</b> <span class="fc-secondary">{'·'} this view is in the address bar: copy it to share</span></Crumbs>
    <h1 class="fc-sr">Compare {F.plural(ids.length, 'run')}</h1>
    {dropped.length ? <div class="fc-banner fc-compare__notices" role="status">{dropped.map(t => <span class="fc-notice">{t}</span>)}</div> : null}</>;
  const chipsOf = save => <div class="fc-items">{rows.map((r, k) => <div class="fc-item">
    <span class="fc-item__letter" style={{ color: colorOf(k) }}><b>{LETTERS[k]}</b><svg width="26" height="8" aria-hidden="true"><line x1="0" y1="4" x2="26" y2="4" stroke="currentColor" stroke-width="2" stroke-dasharray={DASH[k].join(' ') || undefined} /></svg></span>
    <span class="fc-item__body"><span><a class="fc-mono fc-id" href={hrefRun(r.id)}>{r.fc}</a> <span class="fc-secondary"><When iso={r.started} floating={r.stub} />{k === 0 ? ' · reference' : ''}</span> <Badge v={r.health} soFar={r.provisional} /></span>
      <span class="fc-item__name">{r.run_name ? <Mid text={r.run_name} tail={12} /> : 'no run name'}</span></span>
    <button type="button" class="fc-close" aria-label={'Remove ' + r.fc + ' from the comparison'} onClick={() => remove(r.id)}>Remove</button></div>)}
    {ids.length < 8 ? <a class="fc-btn fc-btn--ghost" href="#/runs">Add a run from the Runs table (up to 8)</a> : null}{save ? <span class="fc-items__actions"><SavePdf /></span> : null}</div>, chips = chipsOf(false);
  if (ids.length < 2) return <div class="fc-page">{crumbs}{chips}<div class="fc-empty fc-page__empty"><b>Compare needs two to eight runs.</b><span>Tick runs in the Runs table and choose Compare. <a href="#/runs">To the Runs table</a></span></div></div>;
  if (!ready) return <div class="fc-page">{crumbs}{chips}<p class="fc-loading">{ids.some(id => bundleState(id) === 'failed') ? 'One of these runs could not be read.' : 'Loading the runs …'}</p></div>;
  const withT1 = bundles.map(b => !D.flags(b).stub && !!D.output(b) && !!D.activity(b));
  if (!withT1.every(Boolean)) return <div class="fc-page">{crumbs}{chips}<div class="fc-empty fc-page__empty"><b>{rows.filter((_, k) => !withT1[k]).map(r => r.fc).join(', ')}: only raw data was found.</b><span>There is no run health to compare. Remove it above to compare the others.</span></div></div>;

  // ---- the common window: the longest every item reached ----
  const durs = bundles.map(D.runHours), T = Math.max(1 / 6, Math.floor(Math.min(...durs) * 6) / 6), Tw = F.hours(T).replace(' ', ' ');
  const fac = bundles.map(b => D.factors(b, T)), ref = fac[0];
  // "no real difference": below the lab's usual run-to-run spread of output for this kind of run, and never below 10 %
  const usual = (() => { const g = ix.groups.get(rows[0].group), m = g && g.metrics && g.metrics.yield_estimated_bases_at_48h; return Math.max(0.10, m && m.n >= 5 && m.spread ? m.spread : 0.10); })();
  const pores = bundles.map(b => D.summary(b, 'fc').fc_pores_first_scan);
  const norm = k => showAs === 'pores' ? 1000 / (pores[k] || 1) : showAs === 'chh' ? 1 / (fac[k].workHours || 1) : 1;
  const unitWord = showAs === 'pores' ? ' per 1,000 starting pores' : showAs === 'chh' ? ' per working channel-hour' : '';

  // ---- notices ----
  const notices = [];
  const assays = new Set(rows.map(r => r.assay_group)); if (assays.size > 1) notices.push('Different assays (' + [...assays].map(a => ASSAY_WORD[a] || a).join(', ') + '): length and library numbers are not comparable.');
  const mk = new Set(rows.map(r => (r.minknow || '').split('.').slice(0, 2).join('.'))); if (mk.size > 1) notices.push('MinKNOW versions are mixed (' + [...mk].join(', ') + '): MinKNOW’s speed estimate changed between versions, so speed and estimated bases are not on the same footing.');
  const reused = rows.filter(r => r.use > 1); if (reused.length && reused.length < rows.length) notices.push(reused.map(r => LETTERS[rows.indexOf(r)]).join(', ') + (reused.length > 1 ? ' are re-used flow cells' : ' is a re-used flow cell') + ': pores and total output are shown but say little about the library.');
  const grp = new Set(rows.map(r => r.group)); if (grp.size > 1) notices.push('These runs are not all the same device, flow cell type and kit, so the lab’s typical values differ between them.');
  const noT2 = rows.filter(r => r.read_data === 'absent'), someT2 = noT2.length < rows.length; if (someT2 && noT2.length) notices.push(noT2.map(r => LETTERS[rows.indexOf(r)]).join(', ') + (noT2.length > 1 ? ' have' : ' has') + ' no read data: read charts are drawn for the others only.');
  const pops = new Set(rows.filter(r => r.t2_population).map(r => r.t2_population)); if (pops.size > 1) notices.push('The read data of these runs describes different sets of reads (' + [...pops].map(k => POPULATION_WORD[k] || k.replace(/_/g, ' ')).join(', ') + '): read lengths are not like for like.');
  const part = rows.filter(r => r.read_data === 'partial'); if (part.length) notices.push(part.map(r => LETTERS[rows.indexOf(r)]).join(', ') + (part.length > 1 ? ' have' : ' has') + ' partial read data: its read numbers will still change.');
  const quick = rows.filter(r => r.detail_t1 === 'quick' && r.pending > 0); if (quick.length) notices.push(quick.map(r => LETTERS[rows.indexOf(r)]).join(', ') + (quick.length > 1 ? ' are' : ' is') + ' still counting reads: speed, why reads ended and estimated read length are missing there.');

  // ---- what differs ----
  const S = bundles.map(b => { const o = {}; for (const k of Object.keys(b.modules)) Object.assign(o, b.modules[k].summary || {}); return o; });
  const at = (b, k) => { const o = D.output(b); return o ? [o.bases[Math.round(T * 6)], o.reads[Math.round(T * 6)]] : [null, null]; };
  const defs = [
    ['Setup', [['MinKNOW', k => rows[k].minknow || F.DASH, null], ['Kit', k => rows[k].kit || F.DASH, null], ['Flow cell type', k => rows[k].product || F.DASH, null], ['Device position', k => rows[k].position || F.DASH, null],
      ['Run time', k => F.hours(durs[k]), k => durs[k]], ['Use of the flow cell', k => rows[k].use > 1 ? 're-used · use ' + rows[k].use : 'first use', null]]],
    ['Run health, first ' + Tw, [['Estimated bases at ' + Tw, k => F.bases(fac[k].bases), k => fac[k].bases], ['Reads at ' + Tw, k => F.count(fac[k].reads), k => fac[k].reads],
      ['Working channel-hours', k => F.int(fac[k].workHours), k => fac[k].workHours], ['Share of working time spent sequencing', k => F.pct(fac[k].seqShare), k => fac[k].seqShare], ['Estimated bases per sequencing second', k => F.sig(fac[k].basesPerSeqSec), k => fac[k].basesPerSeqSec]]],
    ['Run health, whole run', [['Pores at the first scan', k => F.int(pores[k]), k => pores[k]], ['Occupancy, first hour', k => durs[k] >= 1.5 ? F.pct(S[k].fc_occupancy_first_hour) : F.DASH, k => durs[k] >= 1.5 ? S[k].fc_occupancy_first_hour : null], ['Occupancy over the run', k => F.pct(S[k].fc_occupancy_median), k => S[k].fc_occupancy_median],
      ['Half of the pores gone after', k => S[k].fc_pore_half_life_h != null ? F.hours(S[k].fc_pore_half_life_h) : F.DASH, k => S[k].fc_pore_half_life_h], ['Blocked pores', k => F.pct(S[k].fc_recovering_frac), k => S[k].fc_recovering_frac],
      ['Speed through the pore, estimated', k => S[k].dyn_speed_est_median != null ? Math.round(S[k].dyn_speed_est_median) + ' bases/s' : 'being counted', k => S[k].dyn_speed_est_median], ['Reads that ended normally', k => S[k].dyn_signal_positive_frac != null ? F.pct1(S[k].dyn_signal_positive_frac) : 'being counted', k => S[k].dyn_signal_positive_frac],
      ['Estimated read N50', k => S[k].len_est_n50 != null ? F.bases(S[k].len_est_n50) : 'being counted', k => S[k].len_est_n50]]],
  ];
  // read data: FlowCellar's own numbers from the linked read files; a run without them shows a dash
  if (someT2) defs.push(['Read data' + (noT2.length ? ' · ' + noT2.map(r => LETTERS[rows.indexOf(r)]).join(', ') + ' ' + (noT2.length > 1 ? 'have' : 'has') + ' none' : ''), [
    ['Usable reads', k => rows[k].t2_usable != null ? F.pct(rows[k].t2_usable) : F.DASH, k => rows[k].t2_usable], ['Passed bases', k => rows[k].t2_passed_gb != null ? F.bases(rows[k].t2_passed_gb * 1e9) : F.DASH, k => rows[k].t2_passed_gb],
    ['Read N50', k => rows[k].t2_n50 != null ? F.bases(rows[k].t2_n50) : F.DASH, k => rows[k].t2_n50], ['Median read quality', k => rows[k].t2_q_med != null ? F.q(rows[k].t2_q_med) : F.DASH, k => rows[k].t2_q_med],
    ['Mapped reads', k => rows[k].t2_mapped_frac != null ? F.pct1(rows[k].t2_mapped_frac) : F.DASH, k => rows[k].t2_mapped_frac],
    ['Full-length reads (cDNA)', k => rows[k].t2_full_length != null ? F.pct1(rows[k].t2_full_length) : F.DASH, k => rows[k].t2_full_length]]]);
  let nRows = 0, nDiff = 0;
  const table = defs.map(([group, list]) => [group, list.map(([name, text, num]) => {
    const texts = ids.map((_, k) => text(k)); let spread = null, differs;
    if (num) { const vals = ids.map((_, k) => num(k)).filter(F.isNum); if (vals.length >= 2) { const s = vals.slice().sort((a, b) => a - b), med = s[s.length >> 1] || 1; spread = (s[s.length - 1] - s[0]) / Math.abs(med || 1); } differs = spread != null && spread >= 0.10; }
    else differs = new Set(texts).size > 1;
    nRows++; if (differs) nDiff++;
    return { name, texts, spread, differs, setup: !num };
  })]);

  // ---- overlays ----
  const xT = { domain: [0, T], kind: 'time' };
  const mk1 = (k, x, y) => ({ name: LETTERS[k], x, y, color: colorOf(k), dash: DASH[k], end: LETTERS[k], width: ids.length > 4 ? 1.5 : 2 });
  const legend = rows.map((r, k) => ({ name: LETTERS[k], label: LETTERS[k] + ' · ' + r.fc + ' ' + F.midCut(r.run_name || '', 28), kind: 'line', color: colorOf(k), dash: !!DASH[k].length }));
  const outs = bundles.map(D.output), acts = bundles.map(D.activity), scs = bundles.map(D.scans);
  const unit = by === 'reads' ? 1e6 : 1e9;
  const outSeries = bundles.map((b, k) => { const o = outs[k], n = Math.round(T * 6) + 1; return mk1(k, Array.from(o.t.slice(0, n)), Array.from((by === 'reads' ? o.reads : o.bases).slice(0, n), v => v / unit * norm(k) * (showAs === 'measured' ? 1 : 1e3))); });
  const yTitle = (by === 'reads' ? 'Reads' : 'Estimated bases') + (showAs === 'measured' ? (by === 'reads' ? ', millions' : ', Gb') : (by === 'reads' ? ', thousands' : ', Mb') + unitWord);
  const pick = (ser, i, fmt) => ser.map(q => { const j = Math.min(q.y.length - 1, i); return [fmt(q.y[j]), q.name]; });
  const longest = outs[outs.reduce((m, o, k) => (o.t.length > outs[m].t.length ? k : m), 0)].t;
  const outSpec = { form: 'xy', sync: 'cmp', mr: 26, x: xT, y: { title: yTitle }, series: outSeries, cursor: { x: Array.from(longest.slice(0, Math.round(T * 6) + 1)) }, readout: i => [['Hour ' + +(i / 6).toFixed(1), '']].concat(pick(outSeries, i, v => F.sig(v))) };
  // hourly median of the 10-minute steps: a pore scan takes channels away for a few minutes every hour and a half, and would draw as a comb
  const hourly = a => { const n = Math.min(a.n, Math.round(T * 6)), xs = [], ys = []; for (let i = 0; i < n; i += 6) { const w = Array.from(a.working.slice(i, Math.min(n, i + 6))).sort((p, q) => p - q); xs.push(i / 6 + 0.5); ys.push(w[w.length >> 1]); } return [xs, ys]; };
  const wSeries = bundles.map((b, k) => { const [xs, ys] = hourly(acts[k]), peak = Math.max(1, ...ys); return mk1(k, xs, ys.map(v => showAs === 'measured' ? v : v / peak * 100)); });
  const wSpec = { form: 'xy', sync: 'cmp', mr: 26, x: xT, y: { title: (showAs === 'measured' ? 'Working channels' : 'Working channels, % of each run’s peak') + ', hour by hour' }, series: wSeries, cursor: { x: wSeries.reduce((m, q) => q.x.length > m.length ? q.x : m, []) }, readout: i => [['Hour ' + i + ' to ' + (i + 1), '']].concat(pick(wSeries, i, v => F.int(v))) };
  const pSeries = bundles.map((b, k) => { const s = scs[k]; const xs = [], ys = []; if (s) for (let i = 0; i < s.n; i++) if (s.t[i] <= T + 1e-6) { xs.push(s.t[i]); ys.push(showAs === 'pores' ? s.pores[i] / (s.pores[0] || 1) * 100 : s.pores[i]); } return { ...mk1(k, xs, ys), type: 'step' }; });
  const scanX = pSeries.reduce((m, q) => q.x.length > m.length ? q.x : m, []);
  const pSpec = { form: 'xy', sync: 'cmp', mr: 26, x: xT, y: { title: showAs === 'pores' ? 'Working pores, % of the first scan' : 'Working pores at each scan' }, series: pSeries, cursor: { x: scanX }, readout: i => [['Scan ' + (i + 1) + ', about hour ' + (scanX[i] || 0).toFixed(1), '']].concat(pick(pSeries, i, v => F.int(v))) };
  const oSeries = bundles.map((b, k) => { const a = acts[k], n = Math.min(a.n, Math.round(T * 6)); const xs = [], ys = []; for (let i = 0; i < n; i += 3) { let s = 0, m = 0; for (let j = i; j < Math.min(n, i + 3); j++) if (a.occ[j] != null) { s += a.occ[j]; m++; } xs.push(a.t[Math.min(n - 1, i + 1)]); ys.push(m ? s / m : null); } return mk1(k, xs, ys); });
  const oSpec = { form: 'xy', sync: 'cmp', mr: 26, x: xT, y: { domain: [0, 1], ticks: [0, 0.25, 0.5, 0.75, 1], fmt: v => Math.round(v * 100) + (v === 1 ? ' %' : ''), title: 'Occupancy: busy share of working pores, in half-hour steps' }, series: oSeries, limits: [{ y: 0.7, kind: 'warn', label: 'ONT’s guide: 70 % or more', fit: false }], cursor: { x: oSeries[0].x }, readout: i => [['Hour ' + (oSeries[0].x[i] || 0).toFixed(1), '']].concat(pick(oSeries, i, v => v == null ? F.DASH : F.pct(v))) };
  // read length: basecalled where every item has it, else MinKNOW's estimate where every item was counted
  let lenPanel = null;
  // MinKNOW's estimate exists for every counted run and is like for like; basecalled read lengths are used when a run was not counted and every run has read data
  const allFull = bundles.every(b => D.flags(b).full && D.lenEst(b)), allT2 = !allFull && bundles.every(b => D.mod(b, 'len') && D.readCells(b, 'len').length);
  if (allT2 || allFull) {
    const hs = bundles.map(b => { if (allT2) { const keys = [...new Set(D.readCells(b, 'len').map(q => q.sample))], h = D.readHists(b, keys, ['pass']); return { edges: h.len.edges, a: by === 'reads' ? h.len.n : h.len.w }; } const l = D.lenEst(b); return { edges: l.edges, a: by === 'reads' ? l.n : l.w }; });
    const sameAxis = hs.every(h => h.edges.length === hs[0].edges.length);
    if (sameAxis) {
      const ser = hs.map((h, k) => { const t = D.sum(h.a) || 1; return { name: LETTERS[k], type: 'outline', edges: h.edges, y: Array.from(h.a, v => v / t * 100), color: colorOf(k), dash: DASH[k], width: 1.75 }; });
      // the axis ends where 99.9 % of each run's reads (or bases) lie: a handful of very long reads do not stretch it
      const lo = Math.min(...hs.map(h => D.spanQ(h.edges, h.a, 0.001, 2)[0])), hi = Math.max(...hs.map(h => D.spanQ(h.edges, h.a, 0.001, 2)[1])), e = hs[0].edges;
      lenPanel = <ChartPanel title={allT2 ? 'Read length' : 'Estimated read length'} reading={'Each run is drawn as an outline scaled to its own total, so a deep and a shallow run compare in shape.' + (allT2 ? ' Reads that passed the quality filter.' : ' MinKNOW’s estimate from the signal, before basecalling.')} tag={allT2 ? 'whole run' : 'ESTIMATED · whole run'} height={230}
        spec={{ form: 'xy', x: { domain: [lo, hi], log: true, kind: 'length', title: 'Read length, logarithmic' }, y: { title: 'Share of each run’s ' + (by === 'reads' ? 'reads' : 'bases') + ' per bin, %' }, series: ser, cursor: { x: Array.from({ length: e.length - 1 }, (_, i) => Math.sqrt(e[i] * e[i + 1])) }, readout: i => [[F.int(e[i]) + ' to ' + F.int(e[i + 1]) + ' bases', '']].concat(ser.map(q => [q.y[i].toFixed(2) + ' %', q.name])) }} legend={legend}
        table={() => ({ head: ['Read length, bases'].concat(ser.map(q => q.name + ', %')), rows: Array.from({ length: e.length - 1 }, (_, i) => [F.int(e[i]) + ' to ' + F.int(e[i + 1])].concat(ser.map(q => q.y[i].toFixed(2)))).filter((_, i) => ser.some(q => q.y[i] > 0)) })} />;
    }
  }
  return <div class="fc-page fc-compare">
    <PaperHead kind={'Comparison of ' + F.plural(ids.length, 'run')} bundles={bundles} file={'FlowCellar comparison ' + rows.map(r => r.fc).join(' ')} />
    {crumbs}{chipsOf(true)}
    {notices.length ? <div class="fc-banner fc-compare__notices">{notices.map(n => <span class="fc-notice">{n}</span>)}</div> : null}
    <div class="fc-controls">
      <span class="fc-ctrl fc-compare__time"><span class="fc-caps">Time</span><span class="fc-compare__timetext">first {Tw}, the length of the shortest run: every single number is cut there</span></span>
      <Seg name="cmp-by" label="Count by" value={by} onChange={v => setParams({ by: v === 'bases' ? '' : v })} opts={[{ v: 'reads', label: 'Reads' }, { v: 'bases', label: 'Bases (est.)' }]} />
      <Seg name="cmp-as" label="Show as" value={showAs} onChange={v => setParams({ as: v === 'measured' ? '' : v })} opts={[{ v: 'measured', label: 'Measured' }, { v: 'pores', label: 'Per 1,000 starting pores' }, { v: 'chh', label: 'Per working channel-hour' }]} />
    </div>
    <PaperScope items={[['Time', 'the first ' + Tw + ', the length of the shortest run: every single number is cut there'], ['Count by', by === 'reads' ? 'reads' : 'bases, estimated'],
      ['Show as', showAs === 'pores' ? 'per 1,000 starting pores' : showAs === 'chh' ? 'per working channel-hour' : 'measured'], ['Rows of “What differs”', only ? 'differences only' : 'all']]} />
    <div class="fc-grid fc-grid--why">
      <section class="fc-panel" aria-labelledby="h-why">
        <div class="fc-panel__head"><h2 class="fc-section__title" id="h-why">Why the output differs</h2><span class="fc-secondary">against A, first {Tw}</span></div>
        <div class="fc-panel__body">
          <div class="fc-tablewrap"><table class="fc-table fc-table--wrap fc-why-table">
            <thead><tr><th scope="col">Run</th><th scope="col" class="fc-num">Output by {Tw}<span class="fc-table__sub">estimated bases</span></th><th scope="col">Flow cell<span class="fc-table__sub">working channel-hours</span></th><th scope="col">Loading<span class="fc-table__sub">share of that time sequencing</span></th><th scope="col">Library<span class="fc-table__sub">est. bases per sequencing second</span></th><th scope="col">In words</th></tr></thead>
            <tbody>{fac.map((f, k) => { const w = k ? words(f, ref, usual) : null;
              return <tr><th scope="row"><b style={{ color: colorOf(k) }}>{LETTERS[k]}</b> <span class="fc-mono">{rows[k].fc}</span></th><td class="fc-num">{F.bases(f.bases)}{k ? <span class="fc-table__sub">{pctDiff(f.bases / ref.bases)}</span> : <span class="fc-table__sub">reference</span>}</td>
                {k ? [f.workHours / ref.workHours, f.seqShare / ref.seqShare, f.basesPerSeqSec / ref.basesPerSeqSec].map(r => <td class="fc-nowrap"><Bar r={r} /> {pctDiff(r)}</td>) : [F.int(f.workHours) + ' h', F.pct(f.seqShare), F.sig(f.basesPerSeqSec)].map(t => <td class="fc-secondary">{t}</td>)}
                <td>{k ? <><b>{w[0]}</b> {w[1]}</> : F.DASH}</td></tr>; })}</tbody></table></div>
          <p class="fc-section__note">Output = working channel-hours {'×'} share of that time sequencing {'×'} estimated bases per sequencing second. The three factors multiply to the output ratio exactly{fac.length > 1 ? ' (B: ' + [fac[1].workHours / ref.workHours, fac[1].seqShare / ref.seqShare, fac[1].basesPerSeqSec / ref.basesPerSeqSec].map(r => r.toFixed(3)).join(' × ') + ' = ' + (fac[1].bases / ref.bases).toFixed(3) + ')' : ''}; the words are a reading of them. Working time is time with a working pore; sequencing time is the part of it on a strand or its adapter.</p>
        </div>
      </section>
      <ChartPanel title="Output over time" reading={'Every run from hour 0, cut at ' + Tw + '.'} tag={(by === 'bases' ? 'ESTIMATED' : '') + (showAs !== 'measured' ? (by === 'bases' ? ' · ' : '') + unitWord.trim() : '')} height={250} spec={outSpec} legend={legend}
        table={() => ({ head: ['Hour'].concat(outSeries.map(q => q.name)), rows: outSpec.cursor.x.map((t, i) => [(+t).toFixed(2)].concat(outSeries.map(q => i < q.y.length ? F.sig(q.y[i]) : F.DASH))).filter((_, i) => i % 6 === 0), note: yTitle + '; one row per hour.' })} />
    </div>
    <section class="fc-panel" aria-labelledby="h-diff">
      <div class="fc-panel__head"><h2 class="fc-section__title" id="h-diff">What differs</h2>
        <Seg name="cmp-only" label="Rows" value={only ? '0' : '1'} onChange={v => setParams({ all: v === '1' ? '1' : '' })} opts={[{ v: '1', label: 'All rows' }, { v: '0', label: 'Differences only' }]} />
        <span class="fc-secondary">{nDiff} of {nRows} rows differ by 10 % or more, or in their setup</span></div>
      <div class="fc-tablewrap"><table class="fc-table fc-difftable"><thead><tr><th scope="col">Number</th>{rows.map((r, k) => <th scope="col" class="fc-num"><b style={{ color: colorOf(k) }}>{LETTERS[k]}</b> <span class="fc-mono">{r.fc}</span></th>)}<th scope="col" title="How far the runs lie apart: the largest value against the smallest. Small: under 10 %. Clear: 10 to 30 %. Large: 30 % or more.">Spread<span class="fc-table__sub">small under 10 % · clear · large from 30 %</span></th></tr></thead>
        <tbody>{table.map(([group, list]) => { const shown = list.filter(r => !only || r.differs); if (!shown.length) return null;
          return [<tr class="fc-table__group"><td colSpan={ids.length + 2}>{group}</td></tr>].concat(shown.map(r => <tr><th scope="row">{r.name}</th>{r.texts.map(t => <td class="fc-num">{t}</td>)}
            <td class="fc-nowrap">{r.setup ? (r.differs ? 'differs' : 'same') : r.spread == null ? F.DASH : <><span class="fc-spreadbar" aria-hidden="true"><i style={{ width: Math.min(40, Math.round(r.spread * 40)) + 'px' }} /></span> {r.spread >= 0.3 ? 'large' : r.spread >= 0.10 ? 'clear' : 'small'} <span class="fc-muted">{Math.round(r.spread * 100)} %</span></>}</td></tr>)); })}</tbody></table></div>
      {only && nRows > nDiff ? <p class="fc-section__note fc-panel__note">{nRows - nDiff} rows that differ by less than 10 % are hidden. The same 10 % is used for every number.</p> : null}
    </section>
    <div class="fc-compare__charts">
    <div class="fc-grid fc-grid--3">
      <ChartPanel title="Working channels" reading="Channels with a working pore, whether or not it held a strand." height={220} spec={wSpec} legend={legend}
        table={() => ({ head: ['Hour'].concat(wSeries.map(q => q.name)), rows: wSpec.cursor.x.map((t, i) => ['Hour ' + i + ' to ' + (i + 1)].concat(wSeries.map(q => i < q.y.length ? F.int(q.y[i]) : F.DASH))) })} />
      <ChartPanel title="Pores per scan" reading="Working wells found at each pore scan, about every hour and a half." height={220} spec={pSpec} legend={legend}
        table={() => ({ head: ['Scan'].concat(pSeries.map(q => q.name)), rows: scanX.map((t, i) => [String(i + 1)].concat(pSeries.map(q => i < q.y.length ? F.int(q.y[i]) : F.DASH))) })} />
      <ChartPanel title="Occupancy" reading="Low occupancy with plenty of pores points at the library or the loading, not at the flow cell." height={220} spec={oSpec} legend={legend}
        table={() => ({ head: ['Hour'].concat(oSeries.map(q => q.name + ', %')), rows: oSeries[0].x.map((t, i) => [(+t).toFixed(1)].concat(oSeries.map(q => q.y[i] == null ? F.DASH : (q.y[i] * 100).toFixed(1)))) })} />
    </div>
    {lenPanel ? <div class="fc-grid fc-grid--2">{lenPanel}</div> : null}
    </div>
    {lenPanel ? null : <p class="fc-section__note fc-notice">Read length is not drawn: it needs every run to have been counted (or every run to have read data).</p>}
  </div>;
}
