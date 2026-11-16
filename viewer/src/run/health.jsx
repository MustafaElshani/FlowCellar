// Run health: Run timeline, Flow cell, Reads through the pore. Everything here comes from the MinKNOW
// run folder alone and exists the moment the run ends.
import { setParams, param } from '../core/state.js';
import { Seg, Tile } from '../shell/parts.jsx';
import { ChartPanel, Track, Legend, ScaleLegend, Empty } from '../charts/ChartPanel.jsx';
import { tk, seqColor } from '../charts/theme.js';
import { layoutFor } from '../charts/layouts.js';
import * as F from '../core/format.js';
import * as D from '../core/derive.js';
import { TRACK_CLOSED, metric, panelTitle, label } from '../core/labels.js';
import { Section, Disclosure, SectionTiles, AllNumbers, countNumbers, ruleState, fmtValue } from './section.jsx';
import { sectionResult, scopeTag, ruleById } from './common.js';
import { Grid } from '../charts/Grid.jsx';

const ML = 58, MR = 18;
const pid = panel => 'p-' + panel.replace('/', '-');
const hourSpan = i => { const a = i * 10, f = m => Math.floor(m / 60) + ':' + String(m % 60).padStart(2, '0'); return 'Hour ' + f(a) + ' to ' + f(a + 10); };
const brush = r => { if (!r) return setParams({ t: '' }); let a = Math.max(0, Math.round(r[0] * 6) / 6), e = Math.round(r[1] * 6) / 6; if (e - a < 0.5) e = a + 0.5; setParams({ t: +a.toFixed(2) + '-' + +e.toFixed(2) }); };

export function RunTimeline({ b, c, evRule, evSec }) {
  const fl = D.flags(b), o = D.output(b), am = D.activity(b), sc = D.scans(b), y = D.summary(b, 'yield'), run = D.summary(b, 'run'), fc = D.summary(b, 'fc');
  if (!o || !am) return null;
  // channel activity is drawn with the pore-scan steps evened out unless the reader asks for the steps as measured
  const ae = D.activityEven(b), even = c.scans !== 'raw' && ae && ae.scanSteps > 0, a = even ? ae : am, meas = am.share;
  const evPanel = evRule && evRule.panel, isEv = p => evPanel === p;
  const x = labels => ({ domain: c.win, kind: 'time', labels });
  const msgs = D.messages(b), tmp = D.temperature(b), volt = D.voltage(b), sp = D.speed(b), er = D.endReasons(b), fcp = (D.mod(b, 'fc').params || {});
  const moreOpen = c.more.has('tracks') || (evPanel && TRACK_CLOSED.has(evPanel));
  // ---- output: the check's own unit while its evidence is shown ----
  const outRule = evRule && /^yield_/.test(evRule.id) ? evRule : null;
  const by = outRule ? (/reads/.test(outRule.id) ? 'reads' : 'bases') : c.by, mode = outRule ? 'cum' : c.out, unit = by === 'reads' ? 1e6 : 1e9, uw = by === 'reads' ? 'M reads' : 'Gb';
  const g = c.group, cur = c.typ && g && g.n_runs >= 5 && by === 'bases' ? g.curve_output_gb_per_3h : null;
  let outSpec, outLegend = [{ name: 'This run', kind: 'line', color: tk('--fc-this-run') }];
  if (mode === 'cum') {
    const cum = by === 'reads' ? o.reads : o.bases, ys = Array.from(cum, v => v / unit), pts = [];
    for (const h of [24, 48, 72]) { const v = D.valueAt(o.t, cum, h); if (v != null && h >= c.win[0] && h <= c.win[1]) pts.push({ x: h, y: v / unit, label: F.sig(v / unit) + (by === 'reads' ? ' M' : ' Gb') }); }
    let band = null;
    if (cur) { let acc = 0; const bx = [0], mid = [0]; cur.median.forEach((v, i) => { if (v == null) return; acc += v; bx.push((i + 1) * 3); mid.push(acc); }); band = { x: bx, lo: bx.map(() => null), hi: bx.map(() => null), mid }; outLegend.push({ name: 'typical for this lab today: median of ' + g.n_runs + ' runs, added up', kind: 'line', color: tk('--fc-chart-typical-line') }); }
    const limits = [];
    if (outRule && outRule.threshold) { const t = outRule.threshold, h = +/_at_(\d+)h/.exec(outRule.id)[1]; if (t.warn != null) limits.push({ y: t.warn / unit, kind: 'warn', label: (t.kind === 'baseline' ? 'this lab’s warning level at ' : 'fixed warning level at ') + h + ' h: ' + F.sig(t.warn / unit) + (by === 'reads' ? ' M' : ' Gb') }); if (t.fail != null) limits.push({ y: t.fail / unit, kind: 'fail', label: 'fail level at ' + h + ' h: ' + F.sig(t.fail / unit) + (by === 'reads' ? ' M' : ' Gb'), fit: false }); }
    outSpec = { form: 'xy', sync: 'time', ml: ML, mr: MR, x: x(false), y: {}, series: [{ name: 'This run', x: o.t, y: ys }], points: pts, band, limits, cursor: { x: o.t }, onBrush: brush,
      readout: i => [['Hour ' + +o.t[i].toFixed(1), ''], [F.sig(ys[i]) + ' ' + uw, by === 'bases' ? 'estimated' : '']] };
  } else {
    const per = by === 'reads' ? o.readsPerH : o.basesPerH, ys = Array.from(per, v => v / unit);
    const band = cur ? { x: cur.median.map((_, i) => i * 3 + 1.5), lo: cur.q25.map(v => v == null ? null : v / 3), hi: cur.q75.map(v => v == null ? null : v / 3), mid: cur.median.map(v => v == null ? null : v / 3) } : null;
    if (cur) outLegend.push({ name: 'typical for this lab today: middle half of ' + g.n_runs + ' runs', kind: 'band' });
    outSpec = { form: 'xy', sync: 'time', ml: ML, mr: MR, x: x(false), y: {}, series: [{ name: 'This run', x: o.th, y: ys }], band, cursor: { x: o.th }, onBrush: brush,
      readout: i => [['Hour ' + i + ' to ' + (i + 1), ''], [F.sig(ys[i]) + ' ' + uw + ' per hour', by === 'bases' ? 'estimated' : '']] };
  }
  // ---- channel activity ----
  const actSpec = { form: 'xy', sync: 'time', ml: ML, mr: MR, x: x(false), y: { domain: [0, 1], ticks: [0, 0.25, 0.5, 0.75, 1], fmt: v => Math.round(v * 100) + (v === 1 ? ' %' : '') }, stack: { x: a.t },
    series: D.STATE_GROUPS.map(([id, name, , token]) => ({ name, label: name, y: a.share[id], color: tk(token) })), cursor: { x: a.t }, onBrush: brush, events: even && sc ? sc.t.map(t => ({ x: t })) : null,
    readout: i => [[hourSpan(i) + (ae && ae.scanStep[i] ? ' · a pore scan ran in this step' : ''), '']].concat(D.STATE_GROUPS.map(([id, name]) => [F.pct(meas[id][i]), name.toLowerCase()])) };
  const actTable = () => ({ head: ['Time'].concat(D.STATE_GROUPS.map(g => g[1] + ', %'), ['Pore scan in this step']), rows: Array.from(am.t, (_, i) => [hourSpan(i)].concat(D.STATE_GROUPS.map(([id]) => (meas[id][i] * 100).toFixed(1)), [ae && ae.scanStep[i] ? 'yes' : ''])),
    note: 'As measured, in 10-minute steps. The chart evens out the steps in which a pore scan ran unless “as measured” is chosen.' });
  const occRule = ruleById(b, 'fc_occupancy_first_hour'), occLim = [{ y: 0.7, kind: 'warn', label: 'ONT’s guide: 70 % or more', fit: false }];
  if (isEv('fc/occupancy') && occRule && occRule.threshold.fail != null) occLim.push({ y: occRule.threshold.fail, kind: 'fail', label: 'fail below ' + F.pct(occRule.threshold.fail), fit: false });
  const occSpec = { form: 'xy', sync: 'time', ml: ML, mr: MR, x: x(false), y: { domain: [0, 1], ticks: [0, 0.5, 1], fmt: v => Math.round(v * 100) + (v === 1 ? ' %' : '') }, series: [{ name: 'Occupancy', x: a.t, y: a.occ }], limits: occLim, cursor: { x: a.t }, onBrush: brush,
    readout: i => [[hourSpan(i), ''], [a.occ[i] == null ? 'no working pore' : F.pct(a.occ[i]), a.occ[i] == null ? '' : 'of working pores busy']] };
  const poreRule = ruleById(b, 'fc_pores_first_scan');
  const scanSpec = sc ? { form: 'xy', sync: 'time', ml: ML, mr: MR, x: x(!moreOpen), y: {}, series: [{ name: 'Pores', type: 'vbar', x: sc.t, parts: [{ y: sc.v.single_pore, color: tk('--fc-scan-single-pore') }, { y: sc.v.reserved_pore, color: tk('--fc-scan-reserved-pore') }] }],
    limits: isEv('fc/pores_per_scan') && poreRule && poreRule.threshold.warn ? [{ y: poreRule.threshold.warn, kind: 'warn', label: 'limit at the first scan: ' + F.int(poreRule.threshold.warn) + ' or more' }] : [], cursor: { x: sc.t }, onBrush: brush,
    readout: i => [['Scan ' + (i + 1) + ', hour ' + sc.t[i].toFixed(1), ''], [F.int(sc.v.single_pore[i]), 'available'], [F.int(sc.v.reserved_pore[i]), 'in reserve']] } : null;
  // ---- events ----
  const evs = (sc ? sc.t.map((t, i) => ({ x: t, text: 'pore scan ' + (i + 1), kind: 'scan' })) : []).concat(msgs.map(m => ({ x: Math.max(0, m.h), text: m.text, kind: m.severity }))).sort((p, q) => p.x - q.x);
  const evSpec = { form: 'xy', sync: 'time', ml: ML, mr: MR, x: x(false), y: { domain: [0, 1], labels: false, ticks: [], grid: false }, series: [], strip: { scans: sc ? sc.t : [], marks: msgs.map(m => ({ x: Math.max(0, m.h), kind: m.severity })) }, cursor: { x: evs.map(e => e.x) }, onBrush: brush,
    readout: i => [['Hour ' + evs[i].x.toFixed(1), ''], evs[i].text] };
  // ---- closed tracks ----
  const tgt = fcp.target_speed_min_bps ? [Math.round(fcp.target_speed_min_bps), Math.round(fcp.target_speed_max_bps)] : null;
  const spSpec = sp ? { form: 'xy', sync: 'time', ml: ML, mr: MR, x: x(false), y: { zero: false }, series: [{ name: 'Median speed', x: sp.t, y: sp.med }], hband: tgt ? { y0: tgt[0], y1: tgt[1], label: 'target ' + tgt[0] + ' to ' + tgt[1] } : null, cursor: { x: sp.t }, onBrush: brush,
    readout: i => [['Hour ' + i + ' to ' + (i + 1), ''], [sp.med[i] == null ? 'too few reads' : Math.round(sp.med[i]) + ' bases per second', sp.med[i] == null ? '' : 'estimated, median of reads that ended normally']] } : null;
  const erSpec = er ? { form: 'xy', sync: 'time', ml: ML, mr: MR, x: x(false), y: { domain: [0, 1], ticks: [0, 0.5, 1], fmt: v => Math.round(v * 100) + (v === 1 ? ' %' : '') }, stack: { x: er.t },
    series: er.rows.map(r => ({ name: r.label, y: Array.from(by === 'bases' ? r.perHB : r.perH, (v, i) => { const t = by === 'bases' ? er.totalB[i] : er.total[i]; return t > 0 ? v / t : 0; }), color: tk(r.token) })), cursor: { x: er.t }, onBrush: brush,
    readout: i => [['Hour ' + i + ' to ' + (i + 1), '']].concat(er.rows.slice(0, 3).map(r => [F.pct1((by === 'bases' ? r.perHB[i] : r.perH[i]) / ((by === 'bases' ? er.totalB[i] : er.total[i]) || 1)), r.label.toLowerCase()])) } : null;
  const tmpSpec = tmp && tmp.t.length ? { form: 'xy', sync: 'time', ml: ML, mr: MR, x: x(false), y: { zero: false, fmt: v => v.toFixed(1) }, series: [{ name: 'Flow cell', x: tmp.t, y: tmp.flowcell, dots: tmp.t.length < 3 }], hband: tmp.lo[0] != null && tmp.hi[0] != null ? { y0: tmp.lo[0], y1: tmp.hi[0], label: 'target ' + tmp.lo[0] + ' to ' + tmp.hi[0] + ' °C' } : null, cursor: { x: tmp.t }, onBrush: brush,
    readout: i => [['Hour ' + +tmp.t[i].toFixed(1), ''], [tmp.flowcell[i] == null ? 'no value' : tmp.flowcell[i].toFixed(2) + ' °C', 'flow cell']] } : null;
  const voSpec = volt ? { form: 'xy', sync: 'time', ml: ML, mr: MR, x: x(true), y: { zero: false }, series: [{ name: 'Voltage', type: 'step', x: volt.t, y: volt.mv }], cursor: { x: volt.t }, onBrush: brush,
    readout: i => [['Hour ' + volt.t[i].toFixed(1), ''], [Math.round(volt.mv[i]) + ' mV', 'at this pore scan']] } : null;
  const half = y.yield_half_time_h, dur = run.run_duration_h;
  const reading = dur < 1 ? 'The run stopped after ' + F.hours(dur) + '; there is little to draw.' : half ? 'Half of the output came in the first ' + F.hours(half) + ' of ' + F.hours(dur) + '.' + (fc.fc_sequencing_duty != null ? ' Channels spent ' + F.pct(fc.fc_sequencing_duty) + ' of the run sequencing.' : '') : null;
  const h = (/_at_(\d+)h/.exec((b.status.rules.find(r => /^yield_.*_at_\d+h$/.test(r.id)) || {}).id || '') || [])[1];
  const tileIds = (h ? ['yield_reads_at_' + h + 'h', 'yield_estimated_bases_at_' + h + 'h'] : []).concat(['fc_occupancy_median', 'fc_recovering_frac', 'fc_recovering_last_third_frac', 'fc_temperature_max_dev_c']);
  const nMsg = msgs.length, nWarn = msgs.filter(m => m.severity === 'warning').length, nErr = msgs.filter(m => m.severity === 'error').length;
  return <Section id="run_timeline" title="Run timeline" result={sectionResult(b, 'run_timeline')} soFar={fl.soFar} reading={reading} x={{ evRule, evSec }}
    tools={<><Seg name="c-out" label="Output" value={mode} onChange={v => setParams({ out: v === 'cum' ? '' : v })} opts={[{ v: 'cum', label: 'Cumulative' }, { v: 'rate', label: 'Per hour' }]} />
      {c.timeOn ? <button type="button" class="fc-btn fc-btn--ghost" onClick={() => setParams({ t: '' })}>Whole run</button> : null}</>}>
    <SectionTiles b={b} ids={tileIds} />
    <p class="fc-tracks__hint fc-screen">One time axis and one cursor for every track: hours since the start of sequencing. Drag across a track, or use the Time control, to look at part of the run; double-click a track for the whole run.</p>
    <div class="fc-tracks">
      <Track id={pid('run/events')} evidence={isEv('run/events')} title="Events" how={metric('run_messages').how} table={() => ({ head: ['Hour', 'Event'], rows: evs.map(e => [e.x.toFixed(2), e.text]) })} sub={F.plural(sc ? sc.n : 0, 'pore scan') + ' · ' + (nErr ? F.plural(nErr, 'error') + ' · ' : '') + F.plural(nWarn, 'warning') + ' · ' + F.plural(nMsg - nWarn - nErr, 'other message')} height={46} spec={evSpec} />
      <Track id={pid('yield/output_over_time')} evidence={isEv('yield/output_over_time')} title="Output" sub={(by === 'reads' ? 'reads, millions' : 'estimated bases, Gb') + (mode === 'rate' ? ' per hour' : ', added up')} tag={by === 'bases' ? 'ESTIMATED' : (c.typ && g && g.n_runs >= 5 ? 'lab typical: estimated bases only' : null)} height={150} spec={outSpec} legend={outLegend.length > 1 ? outLegend : null} how={metric(outRule ? outRule.id : 'yield_reads_at_24h').how}
        table={() => ({ head: ['Hour', 'Reads, added up', 'Estimated bases, added up'], rows: Array.from(o.t, (t, i) => [t.toFixed(2), F.int(o.reads[i]), F.int(o.bases[i])]).filter((_, i) => i % 6 === 0 || i === o.t.length - 1), note: 'One row per hour; the chart has one point per 10 minutes.' })} />
      <Track id={pid('fc/channel_activity')} evidence={isEv('fc/channel_activity')} title="Channel activity" sub={'share of channels' + (even ? ' · pore-scan steps evened out, scans ticked above' : ae && ae.scanSteps ? ' · as measured: each pore scan shows as a notch' : '')} height={even ? 162 : 150} spec={actSpec}
        legend={D.STATE_GROUPS.map(([, name, , token]) => ({ name, color: tk(token) }))} how={metric('fc_channel_state_series').how} table={actTable}
        tools={ae && ae.scanSteps ? <button type="button" class="fc-linkbtn" aria-pressed={!even} onClick={() => setParams({ scans: even ? 'raw' : '' })} title="A pore scan takes every channel away from sequencing for a few minutes, about every hour and a half">{even ? 'Show as measured' : 'Even out pore scans'}</button> : null} />
      <Track id={pid('fc/occupancy')} evidence={isEv('fc/occupancy')} title="Occupancy" sub="busy share of working pores" height={96} spec={occSpec} how={metric('fc_occupancy_first_hour').how}
        table={() => ({ head: ['Time', 'Occupancy, %', 'Working channels'], rows: Array.from(am.t, (_, i) => [hourSpan(i), am.occ[i] == null ? 'no working pore' : (am.occ[i] * 100).toFixed(1), F.int(am.working[i])]) })} />
      {scanSpec ? <Track id={pid('fc/pores_per_scan')} evidence={isEv('fc/pores_per_scan')} title="Pores per scan" sub="working wells at each pore scan" height={moreOpen ? 110 : 128} spec={scanSpec} legend={[{ name: 'Available', color: tk('--fc-scan-single-pore') }, { name: 'In reserve', color: tk('--fc-scan-reserved-pore') }]} how={metric('fc_pores_first_scan').how}
        table={() => ({ head: ['Scan', 'Hour', 'Available', 'In reserve'], rows: sc.t.map((t, i) => [String(i + 1), t.toFixed(2), F.int(sc.v.single_pore[i]), F.int(sc.v.reserved_pore ? sc.v.reserved_pore[i] : 0)]) })} /> : null}
    </div>
    <Disclosure id="tracks" c={c} charts force={moreOpen && !c.more.has('tracks')} title={'4 more tracks on the same time axis: Speed, Why reads ended, Temperature, Voltage'}>
      <div class="fc-tracks">
        {spSpec ? <Track id={pid('dyn/speed_over_time')} evidence={isEv('dyn/speed_over_time')} title="Speed through the pore" sub="bases per second, estimated" tag="ESTIMATED" height={110} spec={spSpec} how={metric('dyn_speed_est_median').how}
          table={() => ({ head: ['Hour', 'Median speed, bases per second', 'Middle half'], rows: sp.t.map((t, i) => ['Hour ' + i + ' to ' + (i + 1), sp.med[i] == null ? 'too few reads' : String(Math.round(sp.med[i])), sp.q1[i] == null ? F.DASH : Math.round(sp.q1[i]) + ' to ' + Math.round(sp.q3[i])]) })} />
          : <div class="fc-track fc-track--na" id={pid('dyn/speed_over_time')}><div class="fc-track__label"><div class="fc-track__title">Speed through the pore</div></div><div class="fc-track__plot fc-muted">Being counted: this track appears when the read count finishes.</div></div>}
        {erSpec ? <Track id={pid('dyn/why_reads_ended')} evidence={isEv('dyn/why_reads_ended')} title="Why reads ended" sub={'share of ' + (by === 'bases' ? 'estimated bases' : 'reads') + ', by hour of read start'} height={110} spec={erSpec} legend={er.rows.map(r => ({ name: r.label, color: tk(r.token) }))} how={metric('dyn_signal_positive_frac').how}
          table={() => ({ head: ['Hour of read start'].concat(er.rows.map(r => r.label + ', reads')), rows: er.t.map((_, i) => ['Hour ' + i + ' to ' + (i + 1)].concat(er.rows.map(r => F.int(r.perH[i])))) })} />
          : <div class="fc-track fc-track--na" id={pid('dyn/why_reads_ended')}><div class="fc-track__label"><div class="fc-track__title">Why reads ended</div></div><div class="fc-track__plot fc-muted">Being counted: this track appears when the read count finishes.</div></div>}
        {tmpSpec ? <Track id={pid('fc/temperature')} evidence={isEv('fc/temperature')} title="Temperature" sub={'flow cell, °C' + (tmp.dropped ? ' · the sample taken as the run stopped is left out' : '')} height={96} spec={tmpSpec} how={metric('fc_temperature_max_dev_c').how}
          table={() => ({ head: ['Hour', 'Flow cell, °C', 'Target, °C'], rows: tmp.t.map((t, i) => [String(+t.toFixed(2)), tmp.flowcell[i] == null ? F.DASH : tmp.flowcell[i].toFixed(2), tmp.lo[i] != null ? tmp.lo[i] + ' to ' + tmp.hi[i] : F.DASH]) })} /> : null}
        {voSpec ? <Track id={pid('fc/voltage')} title="Voltage" sub="mV at each pore scan" height={110} spec={voSpec} table={() => ({ head: ['Scan', 'Hour', 'Voltage, mV'], rows: volt.t.map((t, i) => [String(i + 1), t.toFixed(2), String(Math.round(volt.mv[i]))]) })} />
          : <div class="fc-track fc-track--na"><div class="fc-track__label"><div class="fc-track__title">Voltage</div></div><div class="fc-track__plot fc-muted">Being counted: voltage per scan comes from the pore-scan table, which is read with the sequencing summary.</div></div>}
      </div>
    </Disclosure>
    <Disclosure id="n-timeline" c={c} title={'All numbers (' + countNumbers(b, ['yield', 'run']) + ')'}><AllNumbers b={b} mods={['yield', 'run']} /></Disclosure>
  </Section>;
}

// ---------------------------------------------------------------------------------------------
const MEASURES = [['bases', 'Estimated bases'], ['reads', 'Reads'], ['seq', 'Time sequencing'], ['wells', 'Working wells, first scan'], ['last', 'Last scan with a working well']];
export function FlowCellSection({ b, c, evRule, evSec }) {
  const fl = D.flags(b), fc = D.summary(b, 'fc'), cell = D.cell0(b, 'fc'), sc = D.scans(b), idn = b.identity.flow_cell;
  if (!cell || fl.stub) return null;
  const evPanel = evRule && evRule.panel, ch = cell.channels;
  const reading = F.int(fc.fc_pores_first_scan) + ' working pores at the first scan after loading' + (fc.fc_pore_half_life_h ? '; half were gone after ' + F.hours(fc.fc_pore_half_life_h) + (D.runHours(b) < fc.fc_pore_half_life_h ? ' at the rate seen in this run' : '') : '') + '.';
  // ---- the first pore scan ----
  const wells = sc ? D.SCAN_CLASSES.reduce((t, [id]) => t + (sc.v[id] ? sc.v[id][0] : 0), 0) : 0;
  const firstRows = sc ? D.SCAN_CLASSES.filter(([id]) => sc.v[id]).map(([id, label, token]) => ({ label, value: sc.v[id][0], text: F.int(sc.v[id][0]), note: F.pct(sc.v[id][0] / (wells || 1)) + ' of wells', color: tk(token) })) : [];
  const chRule = ruleById(b, 'fc_channels_with_pore_first_scan');
  // ---- map ----
  let map = null;
  if (ch && ch.fc_channel_map) {
    const m = param('map', c.by === 'reads' ? 'reads' : 'bases'), L = layoutFor(ch.fc_channel_map.layout, idn.channel_count || ch.fc_channel_map.values.reads.length), v = ch.fc_channel_map.values, wm = ch.fc_well_map && ch.fc_well_map.values;
    let values, max, colorOf, nodata, unitWord, fmt, scale;
    const p99 = arr => { const s = Array.from(arr).filter(q => q > 0).sort((p, q) => p - q); return s.length ? s[Math.min(s.length - 1, Math.floor(s.length * 0.99))] : 1; };
    if (m === 'reads') { values = v.reads; max = p99(values); fmt = q => F.int(q) + ' reads'; scale = [F.count(0), F.count(max) + ' reads']; }
    else if (m === 'seq' && v.sequencing_seconds) { values = v.sequencing_seconds; max = p99(values); fmt = q => F.hours(q / 3600) + ' sequencing'; scale = ['0', F.hours(max / 3600)]; }
    else if (m === 'wells' && wm && wm.good_wells_first_scan) { values = wm.good_wells_first_scan; max = 4; nodata = q => !(q > 0); colorOf = q => tk('--fc-ord-' + Math.min(5, q + 1)); fmt = q => q + ' working ' + (q === 1 ? 'well' : 'wells') + ' at the first scan'; }
    else if (m === 'last' && wm && wm.last_good_scan && sc) { values = wm.last_good_scan.map(i => i > 0 ? sc.t[Math.min(sc.n - 1, i - 1)] + 0.01 : 0); max = sc.t[sc.n - 1] || 1; fmt = q => q > 0 ? 'last working well at hour ' + q.toFixed(1) : 'no working well at any scan'; scale = ['hour 0', 'hour ' + Math.round(max)]; }
    else { values = v.estimated_bases; max = p99(values); fmt = q => F.bases(q) + ' estimated'; scale = ['0', F.bases(max)]; }
    const none = values.reduce((t, q) => t + ((nodata ? nodata(q) : !(q > 0)) ? 1 : 0), 0), st = scopeTag(c, [], L.physical === 'order' ? 'channel order, not positions' : L.physical === 'inferred' ? 'positions inferred' : null);
    map = <ChartPanel id={pid('fc/channel_map')} evidence={evPanel === 'fc/channel_map'} title="Flow cell map" cls="fc-chart--map"
      reading={F.int(fc.fc_active_channels) + ' of ' + F.int(idn.channel_count) + ' channels produced reads; ' + F.int(none) + (m === 'wells' ? ' had no working well at the first scan.' : m === 'last' ? ' never had a working well.' : ' are drawn hollow because they have no ' + (m === 'reads' ? 'reads' : m === 'seq' ? 'sequencing time' : 'output') + '.')}
      scope={st.scope} tag={st.tag} aspect={L.cols + ' / ' + L.rows} maxWidth={L.cols * 16}
      idle="Point at a cell to read its channel." how={metric('fc_active_channels').how}
      table={() => ({ head: ['Channel', 'Reads', 'Estimated bases', 'Time sequencing, hours'], rows: Array.from(v.reads, (q, k) => [F.int(k + 1), F.int(q), F.int(v.estimated_bases[k]), v.sequencing_seconds ? (v.sequencing_seconds[k] / 3600).toFixed(2) : F.DASH]) })}
      spec={{ form: 'map', cols: L.cols, rows: L.rows, pos: L.pos, at: (r, cc) => { const k = L.at(r, cc); return k != null && k < values.length ? k : null; }, values, max, colorOf, nodata,
        readout: k => [['Channel ' + F.int(k + 1), '']].concat(m === 'reads' || m === 'bases' ? [] : [[fmt(values[k]), '']], [[F.bases(v.estimated_bases[k]), 'estimated'], [F.int(v.reads[k]), 'reads']]) }}
      footer={<div class="fc-chart__foot">
        <label class="fc-ctrl"><span class="fc-caps">Show</span><select class="fc-selectbox" value={m} onChange={e => setParams({ map: e.currentTarget.value })}>{MEASURES.filter(([id]) => ((id !== 'wells' && id !== 'last') || wm) && (id !== 'seq' || v.sequencing_seconds) && (id !== 'last' || (wm && wm.last_good_scan))).map(([id, name]) => <option value={id}>{name}</option>)}</select></label>
        {m === 'wells' ? <Legend items={[1, 2, 3, 4].map(k => ({ name: k + (k === 1 ? ' well' : ' wells'), color: tk('--fc-ord-' + (k + 1)) })).concat([{ name: 'none', kind: 'none' }])} /> : <ScaleLegend lo={scale[0]} hi={scale[1]} none={m === 'last' ? 'never' : 'none'} />}
        <span class="fc-muted fc-chart__note">{L.note}</span>
      </div>} />;
  }
  // ---- more charts ----
  const grid = cell.grids && cell.grids.fc_channel_time_grid, mux = cell.counts && cell.counts.fc_mux_usage, opc = cell.hists2d && cell.hists2d.fc_open_pore_current;
  let gridPanel = null, muxPanel = null, opcPanel = null;
  if (c.more.has('fc')) {
    if (grid) {
      const gc = grid.channels.strand_seconds, codes = D.decodeU6(gc.data), rows = grid.rows, cols = grid.cols, pool = rows > 1000 ? 12 : 4, nr = Math.ceil(rows / pool), colS = b.axes[grid.cols_axis].step;
      const I = [], J = [], T = [];
      for (let r = 0; r < nr; r++) for (let k = 0; k < cols; k++) { let s = 0, n = 0; for (let q = r * pool; q < Math.min(rows, (r + 1) * pool); q++) { const v = codes[q * cols + k] * gc.scale; s += v * v; n++; } const t = s / n / colS; if (t > 0.004) { I.push(k); J.push(r); T.push(Math.min(1, t)); } }
      const xe = Array.from({ length: cols + 1 }, (_, i) => i * colS / 3600), ye = Array.from({ length: nr + 1 }, (_, i) => Math.min(rows, i * pool) + 1), st = scopeTag(c, ['time']);
      // what the read-out says of the step under the pointer: the channels that sequenced in it, and their mean share of the step
      const onCol = new Float64Array(cols), secCol = new Float64Array(cols); for (let q = 0; q < rows; q++) for (let k = 0; k < cols; k++) { const sv = codes[q * cols + k] * gc.scale, sec = sv * sv; if (sec > 0) onCol[k]++; secCol[k] += sec; }
      const hmm = h => { const m = Math.round(h * 60); return Math.floor(m / 60) + ':' + String(m % 60).padStart(2, '0'); };
      gridPanel = <ChartPanel id={pid('fc/channel_time')} title={panelTitle('fc/channel_time')} reading={'Share of each ' + F.hours(colS / 3600).replace(' ', ' ') + ' step that channels spent sequencing; channels in number order, ' + pool + ' per row.'} scope={st.scope} tag={st.tag} height={260}
        spec={{ form: 'xy', sync: 'time', x: { domain: c.win, kind: 'time' }, y: { domain: [1, rows], title: 'Channel number', ticks: F.niceTicks(0, rows, 5).map(v => Math.max(1, v)) }, heat: { xe, ye, i: I, j: J, t: T }, series: [],
          cursor: { x: Array.from({ length: cols }, (_, k) => (k + 0.5) * colS / 3600) }, readout: k => [['Run time ' + hmm(xe[k]) + ' to ' + hmm(xe[k + 1]), ''], [F.int(onCol[k]), 'of ' + F.int(rows) + ' channels sequencing'], [F.pct1(secCol[k] / rows / colS), 'of the step, mean of all channels']] }}
        table={() => ({ head: ['Time'].concat(['Channels sequencing', 'Mean share of the step spent sequencing, %']), rows: Array.from({ length: cols }, (_, k) => { let on = 0, tot = 0; for (let q = 0; q < rows; q++) { const sv = codes[q * cols + k] * gc.scale, sec = sv * sv; if (sec > 0) on++; tot += sec; } return ['Hour ' + +(k * colS / 3600).toFixed(1) + ' to ' + +((k + 1) * colS / 3600).toFixed(1), F.int(on), (tot / rows / colS * 100).toFixed(1)]; }),
          note: 'One row per time step, over all ' + F.int(rows) + ' channels; the chart shows every channel.' })}
        footer={<div class="fc-chart__foot"><ScaleLegend lo="0 %" hi="100 % of the time" none={null} /></div>} />;
    }
    if (mux) muxPanel = <ChartPanel id={pid('fc/well_map')} title={panelTitle('fc/well_map')} reading="Each channel has up to four wells and reads from one at a time; MinKNOW moves to the next when a pore is lost." height={mux.cats.length * 26 + 10}
      spec={{ form: 'hbars', rows: mux.cats.map((k, i) => ({ label: 'Well ' + k, value: c.by === 'bases' ? mux.w[i] : mux.n[i], text: c.by === 'bases' ? F.bases(mux.w[i]) : F.count(mux.n[i]), note: F.pct((c.by === 'bases' ? mux.w[i] : mux.n[i]) / D.sum(c.by === 'bases' ? mux.w : mux.n)) })),
        readout: i => [['Well ' + mux.cats[i], ''], [F.int(mux.n[i]), 'reads'], [F.pct1(mux.n[i] / (D.sum(mux.n) || 1)), 'of the reads']].concat(mux.w ? [[F.bases(mux.w[i]), 'estimated']] : []) }} scope={scopeTag(c, ['by']).scope} tag={scopeTag(c, ['by']).tag}
      table={() => ({ head: ['Well', 'Reads', 'Estimated bases'], rows: mux.cats.map((k, i) => ['Well ' + k, F.int(mux.n[i]), mux.w ? F.int(mux.w[i]) : F.DASH]) })} />;
    if (opc) {
      const ax = b.axes[opc.y], colMax = {}; for (let k = 0; k < opc.i.length; k++) colMax[opc.i[k]] = Math.max(colMax[opc.i[k]] || 0, opc.n[k]);
      const ns = Math.max(...opc.i) + 1, jmin = Math.min(...opc.j), jmax = Math.max(...opc.j), T = opc.n.map((n, k) => Math.sqrt(n / colMax[opc.i[k]]));
      // what the read-out says of the scan under the pointer: its wells, and the current most of them show
      const scanN = new Float64Array(ns), scanTop = new Array(ns).fill(null); for (let k = 0; k < opc.i.length; k++) { const si = opc.i[k]; scanN[si] += opc.n[k]; if (scanTop[si] == null || opc.n[k] > opc.n[scanTop[si]]) scanTop[si] = k; }
      opcPanel = <ChartPanel id={pid('fc/open_pore_current')} title={panelTitle('fc/open_pore_current')} reading="Current through open pores at each pore scan; each scan is scaled to its own most common value." tag={scopeTag(c, []).tag} height={220}
        spec={{ form: 'xy', x: { domain: [1, ns + 1], kind: 'linear', title: 'Pore scan number', fmt: v => String(Math.round(v)) }, y: { domain: [Math.max(0, (jmin - 2) * ax.step), (jmax + 3) * ax.step], title: 'Open-pore current, pA' }, heat: { xe: Array.from({ length: ns + 1 }, (_, i) => i + 1), ye: D.axisEdges(ax), i: opc.i, j: opc.j, t: T }, series: [],
          cursor: { x: Array.from({ length: ns }, (_, i) => i + 1.5) }, readout: i => [['Pore scan ' + (i + 1), ''], [F.int(scanN[i]), 'wells with a current']].concat(scanTop[i] != null ? [[(opc.j[scanTop[i]] * ax.step) + ' to ' + ((opc.j[scanTop[i]] + 1) * ax.step) + ' pA', 'the most common current']] : []) }}
        table={() => ({ head: ['Pore scan', 'Open-pore current, pA', 'Wells'], rows: opc.i.map((si, k) => [String(si + 1), (opc.j[k] * ax.step) + ' to ' + ((opc.j[k] + 1) * ax.step), F.int(opc.n[k])]), note: opc.outside ? F.int(opc.outside) + ' wells fell outside the axis and are not listed.' : null })}
        footer={<div class="fc-chart__foot"><ScaleLegend lo="few wells" hi="most wells of that scan" none={null} /></div>} />;
    }
  }
  const nMore = (grid ? 1 : 0) + (mux ? 1 : 0) + (opc ? 1 : 0);
  return <Section id="flow_cell" title="Flow cell" result={sectionResult(b, 'flow_cell')} soFar={fl.soFar} reading={reading} x={{ evRule, evSec }}>
    <SectionTiles b={b} ids={['fc_pores_first_scan', 'fc_pores_last', 'fc_half', 'fc_channels_with_pore_first_scan', 'fc_active_channels', 'fc_volt']} extra={{
      fc_pores_last: <Tile label="Pores at the last scan" value={F.int(fc.fc_pores_last_scan)} caption={sc ? 'scan ' + sc.n + ', hour ' + sc.t[sc.n - 1].toFixed(0) : null} />,
      fc_half: <Tile label="Half of the pores gone after" value={fc.fc_pore_half_life_h != null ? F.hours(fc.fc_pore_half_life_h).split(' ')[0] : F.DASH} unit={fc.fc_pore_half_life_h != null ? F.hours(fc.fc_pore_half_life_h).split(' ')[1] : ''} caption={fc.fc_pore_survival_24h_frac != null ? F.pct(fc.fc_pore_survival_24h_frac) + ' were left after 24 h' : null} />,
      fc_volt: fc.fc_bias_voltage_start_mv != null ? <Tile label="Voltage, start to end" value={String(Math.round(fc.fc_bias_voltage_start_mv)).replace('-', '−') + ' → ' + String(Math.round(fc.fc_bias_voltage_end_mv)).replace('-', '−')} unit="mV" /> : <Tile label="Voltage, start to end" state="unknown" value="being counted" caption="comes with the pore-scan table" />,
    }} />
    <div class="fc-grid fc-grid--map">
      {map || <Empty title="Flow cell map">{fl.cut ? 'Not available: the read table of this run is cut off or could not be read, and the map needs every read.' : fl.quick ? 'Being counted: the map of output per channel appears when the read count finishes; reload this page then.' : 'No per-channel numbers are in this run’s data.'}</Empty>}
      {firstRows.length ? <ChartPanel id={pid('fc/first_scan')} evidence={evPanel === 'fc/first_scan'} title={panelTitle('fc/first_scan')} cls="fc-chart--side"
        reading={'What the first scan found in the ' + F.int(wells) + ' wells: ' + F.int(fc.fc_channels_with_pore_first_scan) + ' of ' + F.int(idn.channel_count) + ' channels had at least one working pore' + (chRule && chRule.threshold.warn ? ' (the limit is ' + F.pct(chRule.threshold.warn) + ' or more).' : '.')}
        tag={scopeTag(c, []).tag} height={firstRows.length * 26 + 10} spec={{ form: 'hbars', rows: firstRows, readout: i => [[firstRows[i].text, firstRows[i].label.toLowerCase() + ' wells'], firstRows[i].note] }} how={metric('fc_channels_with_pore_first_scan').how}
        table={() => ({ head: ['Found in the well', 'Wells', 'Share of wells'], rows: firstRows.map(r => [r.label, r.text, r.note.replace(' of wells', '')]) })} /> : null}
    </div>
    {nMore ? <Disclosure id="fc" c={c} charts title={'More charts (' + nMore + '): channels over time, which wells were used, open-pore current'}><Grid>{gridPanel}{opcPanel}{muxPanel}</Grid></Disclosure>
      : fl.cut ? <p class="fc-section__note fc-notice">Channels over time and wells used need every read of the sequencing summary, and that file is cut off or could not be read.</p>
      : fl.quick ? <p class="fc-section__note fc-notice">Channels over time, wells used and open-pore current are being counted with the sequencing summary.</p> : null}
    <Disclosure id="n-fc" c={c} title={'All numbers (' + countNumbers(b, ['fc']) + ')'}><AllNumbers b={b} mods={['fc']} /></Disclosure>
  </Section>;
}

// ---------------------------------------------------------------------------------------------
export function ReadsThroughPore({ b, c, evRule, evSec }) {
  const fl = D.flags(b), dy = D.summary(b, 'dyn'), le = D.summary(b, 'len.est');
  if (fl.unread) return <Section id="reads_through_pore" title="Reads through the pore" x={{ evRule, evSec }} reading="Not known: the sequencing summary of this run could not be read.">
    <Empty title="The read table could not be read">Speed, why reads ended and estimated read length come from the sequencing summary, one row per read. The file is there, but it could not be read ({(((D.mod(b, 'dyn') || {}).notes || [])[0] || {}).message || 'see Details, Data sources'}). Another scan will not change this until the file does.</Empty></Section>;
  if (!fl.full && !fl.partlyRead) return <Section id="reads_through_pore" title="Reads through the pore" x={{ evRule, evSec }} reading="Being counted. Speed, why reads ended and estimated read length come from the sequencing summary, one row per read.">
    <Empty title="Being counted">This section appears when the read count finishes. {(b.status.pending || []).length ? <>{(b.status.pending || []).length} checks wait for it: {(b.status.pending || []).map(id => { const l = label(id); return l[0].toLowerCase() + l.slice(1); }).join(', ')}. </> : null}Reload this page then.</Empty></Section>;
  const evPanel = evRule && evRule.panel, L = D.lenEst(b), er = D.endReasons(b), du = D.durationHist(b), sp = D.speed(b);
  if (!L || !er) return null;
  const rl = (D.mod(b, 'dyn') || {}).run_level, poreTypes = rl && rl.counts && rl.counts.dyn_pore_type && rl.counts.dyn_pore_type.cats && rl.counts.dyn_pore_type.cats.length ? rl.counts.dyn_pore_type : null;
  const byB = c.by === 'bases', ly = byB ? L.w : L.n, dom = D.spanQ(L.edges, L.n), n50 = le.len_est_n50, med = le.len_est_median, n50r = ruleById(b, 'len_est_n50');
  const tot = D.sum(ly), cumFrom = (() => { const o = new Float64Array(ly.length); let acc = 0; for (let i = ly.length - 1; i >= 0; i--) { acc += ly[i]; o[i] = acc; } return o; })();
  const ycfg = (title) => c.logY ? { log: true, title: title + ', logarithmic' } : { title };
  const centers = e => Array.from({ length: e.length - 1 }, (_, i) => Math.sqrt(e[i] * e[i + 1]));
  const lenScale = byB ? 1e9 : 1e6, st1 = scopeTag(c, ['by', 'y'], 'ESTIMATED');
  const lenSpec = { form: 'xy', x: { domain: dom, log: true, kind: 'length', title: 'Estimated read length, logarithmic' }, y: ycfg(byB ? 'Estimated bases per bin, Gb' : 'Reads per bin, millions'), series: [{ name: 'All reads', type: 'hist', edges: L.edges, y: Array.from(ly, v => v / lenScale) }],
    vlines: [{ x: n50, label: 'N50 ' + F.bases(n50), strong: true }, { x: med, label: 'median ' + F.bases(med), side: 'left' }].concat(evPanel === 'len/estimated_length' && n50r && n50r.threshold.warn ? [{ x: n50r.threshold.warn, label: 'limit for N50: ' + F.bases(n50r.threshold.warn) }] : []), cursor: { x: centers(L.edges) },
    readout: i => [[F.int(L.edges[i]) + ' to ' + F.int(L.edges[i + 1]) + ' bases', ''], [byB ? F.bases(L.w[i]) : F.int(L.n[i]), byB ? 'estimated' : 'reads'], [F.pct1(ly[i] / tot), 'of the total'], [F.pct(cumFrom[i] / tot), 'in reads this long or longer']] };
  const erTot = byB ? er.bases : er.reads, erRule = ruleById(b, 'dyn_signal_positive_frac');
  const erRows = er.rows.map(r => ({ label: r.label, value: byB ? r.bases : r.reads, text: F.pct1((byB ? r.bases : r.reads) / erTot), note: byB ? F.bases(r.bases) : F.count(r.reads) + ' reads', color: tk(r.token) }));
  const st2 = scopeTag(c, ['by']);
  const dSpec = du ? { form: 'xy', x: { domain: D.spanQ(du.edges, du.n), log: true, kind: 'duration', title: 'Time in the pore, logarithmic' }, y: ycfg('Reads per bin, millions'), series: [{ name: 'Reads', type: 'hist', edges: du.edges, y: Array.from(du.n, v => v / 1e6) }],
    vlines: [{ x: dy.dyn_duration_median_s, label: 'median ' + F.sig(dy.dyn_duration_median_s) + ' s', strong: true }], cursor: { x: centers(du.edges) }, readout: i => [[F.sig(du.edges[i]) + ' to ' + F.sig(du.edges[i + 1]) + ' s', ''], [F.int(du.n[i]), 'reads']] } : null;
  const fcp = D.mod(b, 'fc').params || {}, sdom = sp ? D.spanQ(sp.edges, sp.hist, 0.0005, 4) : null;
  const sSpec = sp ? { form: 'xy', x: { domain: sdom, kind: 'linear', title: 'Bases per second, estimated' }, y: ycfg('Reads per bin, millions'), series: [{ name: 'Reads', type: 'hist', edges: sp.edges, y: Array.from(sp.hist, v => v / 1e6) }],
    vlines: [{ x: dy.dyn_speed_est_median, label: 'median ' + Math.round(dy.dyn_speed_est_median), strong: true }].concat(fcp.target_speed_min_bps ? [{ x: fcp.target_speed_min_bps, label: 'target from ' + Math.round(fcp.target_speed_min_bps), side: 'left' }, { x: fcp.target_speed_max_bps, label: 'to ' + Math.round(fcp.target_speed_max_bps) }] : []),
    cursor: { x: Array.from({ length: sp.edges.length - 1 }, (_, i) => (sp.edges[i] + sp.edges[i + 1]) / 2) }, readout: i => [[sp.edges[i] + ' to ' + sp.edges[i + 1] + ' bases per second', ''], [F.int(sp.hist[i]), 'reads']] } : null;
  const reading = 'Strands moved at about ' + Math.round(dy.dyn_speed_est_median) + ' bases per second and ' + F.pct(dy.dyn_signal_positive_frac) + ' of reads ended normally. Lengths here are MinKNOW’s estimate from the signal, before basecalling.';
  return <Section id="reads_through_pore" title="Reads through the pore" result={sectionResult(b, 'reads_through_pore')} reading={reading} x={{ evRule, evSec }}>
    <SectionTiles b={b} ids={['dyn_speed_est_median', 'dyn_signal_positive_frac', 'dyn_unblock_frac', 'len_est_n50', 'len_est_median', 'dyn_speed_drop_frac']} />
    {fl.partlyRead ? <p class="fc-section__note fc-notice">{((b.status.modifiers || []).find(m => m.code === 'reads_partly_counted') || {}).message || 'The sequencing summary of this run was read in part: the numbers here rest on the reads it holds.'} The file itself is cut off, so another scan will not change this.</p> : null}
    {(((D.mod(b, 'len.est') || {}).notes) || []).filter(n => n.code === 'quantiles_from_histogram').map(n => <p class="fc-section__note">{n.message}</p>)}
    <Grid>
      <ChartPanel id={pid('len/estimated_length')} evidence={evPanel === 'len/estimated_length'} title={panelTitle('len/estimated_length')} reading={'Half of the estimated bases are in reads of ' + F.bases(n50) + ' or longer.'} scope={st1.scope} tag={st1.tag} height={230} spec={lenSpec}
        how={metric('len_est_n50').how} table={() => ({ head: ['Estimated length, bases', 'Reads', 'Estimated bases'], rows: Array.from(L.n, (v, i) => [F.int(L.edges[i]) + ' to ' + F.int(L.edges[i + 1]), F.int(v), F.int(L.w[i])]).filter((_, i) => L.n[i] > 0) })} />
      <ChartPanel id={pid('dyn/why_reads_ended')+'-bars'} evidence={evPanel === 'dyn/why_reads_ended'} title="Why reads ended" reading={F.pct(dy.dyn_signal_positive_frac) + ' of reads ended normally; ' + F.pct(dy.dyn_unblock_frac) + ' were ejected because the strand stalled.' + (erRule && erRule.threshold.warn ? ' The limit for normal endings is ' + F.pct(erRule.threshold.warn) + ' or more.' : '')}
        scope={st2.scope} tag={st2.tag} height={erRows.length * 26 + 10} spec={{ form: 'hbars', rows: erRows, readout: i => [[erRows[i].text, erRows[i].label.toLowerCase()], [F.int(er.rows[i].reads), 'reads'], [F.bases(er.rows[i].bases), 'estimated']] }}
        how={metric('dyn_signal_positive_frac').how} table={() => ({ head: ['Why the read ended', 'Reads', 'Share of reads', 'Estimated bases'], rows: er.rows.map(r => [r.label, F.int(r.reads), F.pct1(r.reads / er.reads), F.int(r.bases)]) })} />
    </Grid>
    <Disclosure id="rtp" c={c} charts title="More charts (2): time each read spent in the pore, speed of reads that ended normally">
      <Grid>
        {dSpec ? <ChartPanel id={pid('dyn/duration')} title={panelTitle('dyn/duration')} reading={'A typical read spent ' + F.sig(dy.dyn_duration_median_s) + ' seconds in the pore.'} scope={scopeTag(c, ['y']).scope} tag={scopeTag(c, ['y']).tag} height={210} spec={dSpec}
          table={() => ({ head: ['Time in the pore, seconds', 'Reads'], rows: Array.from(du.n, (q, i) => [F.sig(du.edges[i]) + ' to ' + F.sig(du.edges[i + 1]), F.int(q)]).filter((_, i) => du.n[i] > 0) })} /> : null}
        {sSpec ? <ChartPanel id={pid('dyn/speed')} title={panelTitle('dyn/speed')} reading={'Reads of 1 second or longer; speed is estimated from MinKNOW’s signal events.'} scope={scopeTag(c, ['y']).scope} tag={scopeTag(c, ['y'], 'ESTIMATED').tag} height={210} spec={sSpec} how={metric('dyn_speed_est_median').how}
          table={() => ({ head: ['Bases per second, estimated', 'Reads'], rows: Array.from(sp.hist, (q, i) => [sp.edges[i] + ' to ' + sp.edges[i + 1], F.int(q)]).filter((_, i) => sp.hist[i] > 0) })} /> : null}
      </Grid>
    </Disclosure>
    <Disclosure id="n-rtp" c={c} title={'All numbers (' + countNumbers(b, ['dyn', 'len.est']) + ')'}><AllNumbers b={b} mods={['dyn', 'len.est']} />
      {poreTypes ? <p class="fc-section__note">Pore types MinKNOW recorded for the reads: {poreTypes.cats.map((k, i) => String(k).replace(/_/g, ' ') + ' ' + F.int(poreTypes.n[i])).join(' · ')}.</p> : null}
      {L.outside ? <p class="fc-section__note">{F.int(L.outside)} reads have an estimated length outside the axis of the length chart (under 10 bases or over 10 Mb); they are counted in the totals and not drawn.</p> : null}</Disclosure>
  </Section>;
}
