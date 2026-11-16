// Flow cell history: every run of one flow cell on one axis. No verdict
// for the flow cell: each run keeps its own. Runs are grouped into loadings; a restart is a short
// segment of the loading it belongs to.
import { index, useBundle, bundleState } from '../core/data.js';
import { setParams, param, hrefRun } from '../core/state.js';
import { Badge, Chip, Mid, Seg, Tile, Crumbs, When } from '../shell/parts.jsx';
import { ChartPanel, Track, Legend, ScaleLegend, Empty } from '../charts/ChartPanel.jsx';
import { tk } from '../charts/theme.js';
import { layoutFor } from '../charts/layouts.js';
import * as F from '../core/format.js';
import * as D from '../core/derive.js';
import { DEVICE_WORD } from '../core/labels.js';
import { SavePdf, PaperHead, PaperScope } from '../shell/paper.jsx';
import { printing } from '../core/print.js';

const ML = 58, MR = 18;
export function FlowCell({ id }) {
  const ix = index.value, rows = ix.byFc.get(id) || [];
  const crumbs = <Crumbs><a href="#/home">Library</a> {'›'} <a href="#/runs">Runs</a> {'›'} Flow cell <b class="fc-mono">{id}</b></Crumbs>;
  if (!rows.length) return <div class="fc-page">{crumbs}<div class="fc-empty fc-page__empty"><b>No run of flow cell {id} is in the library.</b><span><a href="#/runs">To the Runs table</a></span></div></div>;
  const bundles = rows.map(r => useBundle(r.id));
  if (!bundles.every(Boolean)) return <div class="fc-page">{crumbs}<div class="fc-runhead"><div class="fc-runhead__row"><h1 class="fc-runhead__id fc-mono">{id}</h1></div></div><p class="fc-loading">{rows.some(r => bundleState(r.id) === 'failed') ? 'One run of this flow cell could not be read.' : 'Loading the runs of this flow cell …'}</p></div>;
  if (bundles.every(b => D.flags(b).stub)) return <div class="fc-page">{crumbs}<div class="fc-runhead"><div class="fc-runhead__row"><h1 class="fc-runhead__id fc-mono">{id}</h1><span class="fc-runhead__chips"><Chip>{F.plural(rows.length, 'run')}</Chip><Chip>{F.date(rows[0].started)}</Chip></span></div></div>
    <div class="fc-empty fc-page__empty"><b>Only raw data was found for {rows.length === 1 ? 'the run' : 'the runs'} of this flow cell.</b><span>There is no MinKNOW report to draw its history from. {rows.map(r => <a class="fc-mono" href={hrefRun(r.id)}>{r.id}</a>).reduce((a, x, i) => a.concat(i ? [', '] : [], [x]), [])}</span></div></div>;
  const axis = param('axis', 'seq'), n = rows.length, last = n - 1;
  // loadings and restarts: the index says them for every run (use, restart); the newest run's own list is used when it has one
  const use = ((bundles[last].relations || {}).flow_cell_use) || {}, info = r => (use.runs || []).find(u => u.run_key === r.id) || { loading: r.use || 1, restart: !!r.restart };
  const t1 = bundles.map(b => !D.flags(b).stub), durs = bundles.map((b, k) => t1[k] ? D.runHours(b) : 0);
  const t0 = rows[0].startedMs, offs = []; let acc = 0;
  for (let k = 0; k < n; k++) { offs.push(axis === 'cal' ? (rows[k].startedMs - t0) / 3.6e6 : acc); acc += durs[k]; }
  const totalH = durs.reduce((a, v) => a + v, 0), span = axis === 'cal' ? Math.max(...offs.map((o, k) => o + durs[k])) : totalH;
  const gbs = rows.map(r => r.gb || 0), totalGb = gbs.reduce((a, v) => a + v, 0);
  const fcs = bundles.map(b => D.summary(b, 'fc')), firstT1 = t1.indexOf(true), poresNew = fcs[firstT1].fc_pores_first_scan, lastT1 = t1.lastIndexOf(true);
  const loadings = new Set(rows.map(r => info(r).loading)).size, restarts = rows.filter(r => info(r).restart).length;
  const idn = bundles[firstT1].identity, chem = idn.assay.chemistry;
  const many = n > 1, days = Math.round((rows[last].startedMs - t0) / 864e5);
  const projects = new Set(rows.flatMap(r => r.project_ids));
  const head = many
    ? 'Used ' + (n === 2 ? 'twice' : n + ' times') + (loadings < n ? ' in ' + loadings + ' loadings' : '') + (days > 0 ? ' over ' + (days + 1) + ' days' : ' on one day') + '. The first run gave ' + F.bases(gbs[0] * 1e9) + ' estimated; the ' + (n - 1 === 1 ? 'later run' : (n - 1) + ' later starts') + ' added ' + F.bases((totalGb - gbs[0]) * 1e9) + (n > 2 ? ' between them.' : '.')
    : 'Used once' + (t1[0] ? ': ' + F.bases(gbs[0] * 1e9) + ' estimated in ' + F.hours(durs[0]) + '.' : '. Only raw data was found for that run.');
  const detail = many && fcs[lastT1].fc_pores_first_scan && poresNew ? 'The last run started with ' + F.int(fcs[lastT1].fc_pores_first_scan) + ' working pores, ' + F.pct(fcs[lastT1].fc_pores_first_scan / poresNew) + ' of the ' + F.int(poresNew) + ' the flow cell had at its first scan when new.' : 'Every flow cell has this page; a second run of the same flow cell would appear here beside the first.';
  const x = labels => ({ domain: [0, Math.max(span, 0.5)], kind: 'time', labels });
  const vr = many ? offs.slice(1).map((o, i) => ({ x: o })) : [];
  // ---- joined series ----
  const px = [], py = [], firsts = [];
  const ox = [0], oy = [0]; let cum = 0;
  const ax = [], shares = D.STATE_GROUPS.map(() => []);
  bundles.forEach((b, k) => {
    if (!t1[k]) return; const sc = D.scans(b), o = D.output(b), a = D.activityEven(b) || D.activity(b), off = offs[k];
    if (sc) { for (let i = 0; i < sc.n; i++) { px.push(off + sc.t[i]); py.push(sc.pores[i]); } firsts.push({ x: off + sc.t[0], y: sc.pores[0], label: F.int(sc.pores[0]) }); px.push(off + durs[k]); py.push(null); }
    if (o) { if (axis === 'cal' && k) { ox.push(off); oy.push(cum / 1e9); } for (let i = 1; i < o.t.length; i++) { ox.push(off + o.t[i]); oy.push((cum + o.bases[i]) / 1e9); } cum += o.bases[o.t.length - 1]; }
    if (a) { for (let i = 0; i < a.n; i++) { ax.push(off + a.t[i]); D.STATE_GROUPS.forEach(([gid], gi) => shares[gi].push(a.share[gid][i])); }
      if (axis === 'cal' && k < n - 1) { ax.push(off + durs[k] + 0.01); D.STATE_GROUPS.forEach((_, gi) => shares[gi].push(0)); ax.push(offs[k + 1] - 0.01); D.STATE_GROUPS.forEach((_, gi) => shares[gi].push(0)); } }
  });
  const pSpec = { form: 'xy', sync: 'fch', ml: ML, mr: MR, x: x(false), y: {}, series: [{ name: 'Working pores', type: 'step', x: px, y: py, color: tk('--fc-scan-single-pore') }], points: firsts.map(p => ({ ...p, color: tk('--fc-scan-single-pore') })), vrules: vr, cursor: { x: px }, readout: i => { if (py[i] != null) return [['Hour ' + px[i].toFixed(1), ''], [F.int(py[i]), 'working pores']];
      let j = i; while (j > 0 && py[j] == null) j--;            // the end of a run: no scan there, the last scan before it is said
      return py[j] == null ? null : [['Hour ' + px[i].toFixed(1) + ', the end of a run', ''], [F.int(py[j]), 'working pores at its last scan, hour ' + px[j].toFixed(1)]]; } };
  const oSpec = { form: 'xy', sync: 'fch', ml: ML, mr: MR, x: x(false), y: {}, series: [{ name: 'Output', x: ox, y: oy }], vrules: vr, cursor: { x: ox }, readout: i => [['Hour ' + ox[i].toFixed(1), ''], [F.sig(oy[i]) + ' Gb', 'estimated, added up']] };
  const aSpec = { form: 'xy', sync: 'fch', ml: ML, mr: MR, x: x(true), y: { domain: [0, 1], ticks: [0, 0.5, 1], fmt: v => Math.round(v * 100) + (v === 1 ? ' %' : '') }, stack: { x: ax }, series: D.STATE_GROUPS.map(([, label, , token], gi) => ({ name: label, y: shares[gi], color: tk(token) })), vrules: vr, cursor: { x: ax }, readout: i => [['Hour ' + ax[i].toFixed(1), '']].concat(D.STATE_GROUPS.map(([, label], gi) => [F.pct(shares[gi][i]), label.toLowerCase()])) };
  // ---- maps on a shared scale ----
  const maps = bundles.map(b => { const c = D.cell0(b, 'fc'); return c && c.channels && c.channels.fc_channel_map ? c.channels.fc_channel_map : null; });
  let mapMax = 1; for (const m of maps) if (m) { const s = Array.from(m.values.estimated_bases).filter(v => v > 0).sort((p, q) => p - q); if (s.length) mapMax = Math.max(mapMax, s[Math.floor(s.length * 0.99)]); }
  const useWord = (r, k) => (k + 1) + ' · ' + (k === 0 ? 'new' : info(r).restart ? 'restart' : 're-use');
  const segW = k => Math.max(0.4, durs[k] / Math.max(span, 0.5) * 100), dom = Math.max(span, 0.5), paper = printing.value;      // on paper the table has no room for the reason of each run: they stand under it
  return <div class="fc-page fc-flowcell">
    <PaperHead kind="Flow cell history" bundles={bundles} file={'FlowCellar flow cell ' + id} />
    {crumbs}
    <header class="fc-runhead"><div class="fc-runhead__row"><h1 class="fc-runhead__id fc-mono">{id}</h1>
      <span class="fc-runhead__chips">{idn.flow_cell.product_code ? <Chip>{idn.flow_cell.product_code}{chem ? ' · ' + chem.toUpperCase() : ''}</Chip> : null}<Chip>{DEVICE_WORD[idn.device.type] || idn.device.type || 'device not known'} {'·'} {[...new Set(rows.map(r => r.position))].filter(Boolean).join(', ')}</Chip>
        <Chip>{F.plural(n, 'run')}</Chip><Chip>{F.date(rows[0].started, rows[0].stub)}{many ? ' to ' + F.date(rows[last].started, rows[last].stub) : ''}</Chip><Chip>{F.plural(projects.size, 'project')}</Chip></span>
      <span class="fc-runhead__actions"><SavePdf />{many && n <= 8 ? <a class="fc-btn" href={'#/compare?items=' + rows.filter((r, k) => t1[k]).map(r => r.id).join(',')}>Compare its {n} runs</a> : null}</span></div></header>
    <section class="fc-verdict fc-flowcell__summary"><p class="fc-verdict__headline">{head}</p><p class="fc-verdict__detail">{detail}</p></section>
    <div class="fc-keynumbers"><div class="fc-tiles">
      <Tile label="Runs" value={String(n)} caption={many ? F.plural(loadings, 'loading') + (restarts ? ' · ' + F.plural(restarts, 'restart') : '') : 'used once'} />
      <Tile label="Output over its life · estimated" value={F.basesParts(totalGb * 1e9)[0]} unit={F.basesParts(totalGb * 1e9)[1]} caption={many ? 'the first run: ' + F.bases(gbs[0] * 1e9) : null} />
      <Tile label="Hours sequenced" value={F.hours(totalH).split(' ')[0]} unit={F.hours(totalH).split(' ')[1]} caption={F.plural(n, 'run')} />
      <Tile label="Pores when new" value={F.int(poresNew)} caption="first scan of the first run" />
      <Tile label="Pores at the last scan" value={F.int(fcs[lastT1].fc_pores_last_scan)} caption={F.date(bundles[lastT1].identity.protocol_run.ended || rows[lastT1].started)} />
      <Tile label={'Estimated bases per starting pore, ' + (many ? 'last run' : 'this run')} value={F.sig(D.summary(bundles[lastT1], 'yield').yield_per_pore_mb)} unit="Mb" caption={many ? 'the first run: ' + F.sig(D.summary(bundles[0], 'yield').yield_per_pore_mb) + ' Mb' : 'estimated'} />
    </div></div>
    <PaperScope items={[['Time axis', axis === 'seq' ? 'hours sequenced, the runs joined end to end' : 'calendar: hours since the first start, the gaps between runs are real']]} />
    <div class="fc-controls"><Seg name="fch-axis" label="Time axis" value={axis} onChange={v => setParams({ axis: v === 'seq' ? '' : v })} opts={[{ v: 'seq', label: 'Hours sequenced' }, { v: 'cal', label: 'Calendar' }]} />
      <span class="fc-secondary">{axis === 'seq' ? 'Runs are joined end to end, so pore decay reads as one curve.' : 'Hours since the first start: the gaps between runs are real.'}</span></div>
    <section class="fc-section" aria-labelledby="h-life"><div class="fc-section__head"><h2 class="fc-section__title" id="h-life">Life of the flow cell</h2></div>
      <p class="fc-section__reading">One continuous axis; a vertical rule where a run ended and the next began.</p>
      <div class="fc-tracks">
        <div class="fc-track"><div class="fc-track__label"><div class="fc-track__title">Runs</div><div class="fc-track__sub">in the order they started</div></div>
          <div class="fc-track__plot"><div class="fc-runstrip" style={{ marginLeft: ML + 'px', marginRight: MR + 'px' }}>{rows.map((r, k) => <a class={'fc-runstrip__seg' + (info(r).restart ? ' is-restart' : '')} href={hrefRun(r.id)} style={{ left: offs[k] / dom * 100 + '%', width: segW(k) + '%' }} title={'Run ' + (k + 1) + ': ' + F.hours(durs[k]) + ', ' + (r.run_name || '')}><span>{(k + 1) + ' · ' + F.hours(durs[k])}</span></a>)}</div>
            {many ? <ol class="fc-runstrip__list">{rows.map((r, k) => <li><b>{k + 1}</b> <When iso={r.started} floating={r.stub} time /> {'·'} {F.hours(durs[k])} {'·'} loading {info(r).loading}{info(r).restart ? ' · restart' : ''}{k > 0 && axis === 'seq' ? ' · idle ' + F.hours(Math.max(0, (r.startedMs - rows[k - 1].startedMs) / 3.6e6 - durs[k - 1])) + ' before it' : ''}</li>)}</ol> : null}</div></div>
        <Track title="Working pores" sub="at each pore scan; the first scan of each run is marked" height={150} spec={pSpec}
          table={() => ({ head: ['Hour', 'Working pores'], rows: px.map((t, i) => py[i] == null ? null : [t.toFixed(2), F.int(py[i])]).filter(Boolean) })} />
        <Track title="Output" sub="estimated bases, Gb, added up over the runs" tag="ESTIMATED" height={130} spec={oSpec}
          table={() => ({ head: ['Hour', 'Estimated bases added up, Gb'], rows: ox.map((t, i) => [t.toFixed(2), F.sig(oy[i])]).filter((_, i) => i % 6 === 0 || i === ox.length - 1), note: 'One row per hour.' })} />
        <Track title="Channel activity" sub="share of channels · pore-scan steps evened out" height={130} spec={aSpec} legend={D.STATE_GROUPS.map(([, name, , token]) => ({ name, color: tk(token) }))}
          table={() => ({ head: ['Hour'].concat(D.STATE_GROUPS.map(g => g[1] + ', %')), rows: ax.map((t, i) => [t.toFixed(2)].concat(shares.map(sh => (sh[i] * 100).toFixed(1)))).filter((_, i) => i % 6 === 0), note: 'One row per hour, as drawn.' })} />
      </div>
    </section>
    <section class="fc-section" aria-labelledby="h-runs"><div class="fc-section__head"><h2 class="fc-section__title" id="h-runs">Its runs</h2></div>
      <div class="fc-tablewrap"><table class="fc-table fc-table--comfortable fc-table--wrap fc-fcruns"><caption class="fc-sr">Runs of this flow cell</caption>
        <thead><tr><th scope="col">Use</th><th scope="col">Verdict</th><th scope="col">Started<span class="fc-table__sub">{F.zone(rows[0].started)}</span></th><th scope="col">Run name {'·'} project</th><th scope="col" class="fc-num">Run time</th><th scope="col" class="fc-num">Pores at start</th><th scope="col" class="fc-num">Pores at end</th><th scope="col" class="fc-num">Output</th><th scope="col" class="fc-num">Per starting pore</th><th scope="col" class="fc-num">Occupancy</th>{paper ? null : <th scope="col">Reason</th>}</tr></thead>
        <tbody>{rows.map((r, k) => { const b = bundles[k], f = fcs[k], y = D.summary(b, 'yield');
          return <tr><td class="fc-nowrap">{useWord(r, k)}<span class="fc-table__sub">loading {info(r).loading}</span></td><td class="fc-nowrap"><Badge v={r.health} soFar={r.provisional} /></td><td class="fc-nowrap"><When iso={r.started} floating={r.stub} time /></td>
            <th scope="row"><a href={hrefRun(r.id)}><Mid text={r.run_name || r.id} tail={12} /></a><span class="fc-table__sub">{r.project || 'no project'}</span></th>
            <td class="fc-num fc-nowrap">{t1[k] ? F.hours(durs[k]) : F.DASH}</td><td class="fc-num">{F.int(f.fc_pores_first_scan)}</td><td class="fc-num">{durs[k] >= 1.4 ? F.int(f.fc_pores_last_scan) : F.DASH}</td><td class="fc-num fc-nowrap">{t1[k] ? F.bases(r.gb * 1e9) : F.DASH}</td>
            <td class="fc-num fc-nowrap">{y.yield_per_pore_mb != null && durs[k] >= 6 ? F.sig(y.yield_per_pore_mb) + ' Mb' : F.DASH}</td><td class="fc-num">{durs[k] >= 1.5 ? F.pct(f.fc_occupancy_first_hour) : F.DASH}</td>
            {paper ? null : <td class="fc-fcruns__reason">{(r.reasons && r.reasons[0]) || r.headline}</td>}</tr>; })}</tbody></table></div>
      {paper ? <ol class="fc-fcruns__reasons">{rows.map((r, k) => <li><b>{useWord(r, k)}.</b> {(r.reasons && r.reasons[0]) || r.headline}</li>)}</ol> : null}
      <p class="fc-empty fc-flowcell__rule"><b>How a re-used flow cell is judged.</b><span>Only on numbers that do not depend on how many pores were left: output per starting pore, occupancy, speed and read quality. Pore count, total output and pore life are shown beside the flow cell's earlier runs and carry no verdict. A re-use is never used to work out what is typical for the lab. A restart within minutes of a stopped run belongs to the same loading and is marked Incomplete when it is too short to judge.</span></p>
    </section>
    <section class="fc-section" aria-labelledby="h-where"><div class="fc-section__head"><h2 class="fc-section__title" id="h-where">Where on the flow cell</h2></div>
      <p class="fc-section__reading">Estimated bases per channel, one map per run, on one colour scale: a region that died in the first run shows as hollow in the later maps.</p>
      <div class="fc-grid fc-grid--maps">{rows.map((r, k) => { const m = maps[k]; if (!m) return <Empty title={'Run ' + (k + 1)}>{t1[k] ? 'The channel map is being counted.' : 'Only raw data was found.'}</Empty>;
        const L = layoutFor(m.layout, m.values.reads.length), v = m.values.estimated_bases, act = v.reduce((a, q) => a + (q > 0 ? 1 : 0), 0);
        return <ChartPanel cls="fc-chart--flat" title={'Run ' + (k + 1) + ' · ' + F.hours(durs[k])} reading={F.int(act) + ' channels produced reads; ' + F.bases(r.gb * 1e9) + ' estimated.'} aspect={L.cols + ' / ' + L.rows} maxWidth={L.cols * 12}
          spec={{ form: 'map', cols: L.cols, rows: L.rows, pos: L.pos, at: (rr, cc) => { const q = L.at(rr, cc); return q != null && q < v.length ? q : null; }, values: v, max: mapMax, readout: q => [['Channel ' + F.int(q + 1), ''], [F.bases(v[q]), 'estimated']] }}
          tag={L.physical === 'order' ? 'channel order, not positions' : L.physical === 'inferred' ? 'positions inferred' : null}
          table={() => ({ head: ['Channel', 'Reads', 'Estimated bases'], rows: Array.from(v, (q, ch) => [F.int(ch + 1), F.int(m.values.reads[ch]), F.int(q)]) })} />; })}</div>
      <div class="fc-chart__foot"><ScaleLegend lo="0" hi={F.bases(mapMax)} none="no reads" /><span class="fc-muted">{layoutFor(idn.flow_cell.channel_layout, idn.flow_cell.channel_count).note}</span></div>
    </section>
  </div>;
}
