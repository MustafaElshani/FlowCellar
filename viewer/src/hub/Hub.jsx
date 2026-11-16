// Runs: one row per run of the library, drawn from the index alone.
// A real <table>; two text lines and one fixed height per row; three tab stops per row (tick box, flow cell ID,
// Peek); header buttons sort. A column is either whole or hidden: at any width the columns that do not fit are
// hidden in a stated order and counted in the view bar, so no column is ever half visible.
import { useMemo, useRef, useState, useEffect } from 'preact/hooks';
import { index, useBundle, cachedBundle, bundleState } from '../core/data.js';
import { setParams, param, selected, toggleSelected, hrefRun, hrefFlowCell, refocus, textScale, useScrollMemory } from '../core/state.js';
import { makeEngine, hayOf, toggleTerm, hasTerm, removeTerm, verdictHint } from './query.js';
import { COLS, fitColumns, fitWords } from './columns.js';
import { Badge, Glyph, Mid, Spark, Differs, Tile, Updated, When, outputCell, runSoFar, dataWord, dataSoFar, basesLabel, loudNotices } from '../shell/parts.jsx';
import { ChartPanel } from '../charts/ChartPanel.jsx';
import { tk } from '../charts/theme.js';
import * as F from '../core/format.js';
import * as D from '../core/derive.js';
import { ASSAY_SHORT, ASSAY_WORD, DEVICE_WORD, MODIFIER_WORD, POPULATION_WORD, VERDICT_WORD, label } from '../core/labels.js';

const RANK = { fail: 0, warn: 1, incomplete: 2, running: 3, pass: 4, unknown: 5, not_available: 6 };
const SORTS = {
  health: r => RANK[r.health], data: r => RANK[r.data], fc: r => r.fc, started: r => r.startedMs,
  out: r => outputCell(r).sortValue, pores: r => r.pores, occ: r => r.occ_first, passed: r => r.t2_passed_gb, n50: r => r.t2_n50,
};
const DEFAULT_DIR = { health: 1, data: 1, fc: 1, started: -1, out: -1, pores: -1, occ: -1, passed: -1, n50: -1 };

// what stands under a Run health badge, as whole phrases
export function healthWords(r) {
  const out = [], mods = r.modifiers || [];
  if (r.run_state === 'running') out.push('sequencing now');
  if (r.use > 1) out.push('re-used · use ' + r.use);
  if (r.cause === 'raw_only') out.push('raw data only');
  else if (r.cause && r.hours != null) out.push('stopped after ' + F.hours(r.hours));
  for (const m of ['quick_detail', 'baseline_provisional', 'no_history']) if (mods.includes(m)) out.push(m === 'no_history' && out.length ? 'no history' : MODIFIER_WORD[m]);
  for (const m of mods) if (!['quick_detail', 'baseline_provisional', 'no_history', 'reused_flow_cell'].includes(m) && MODIFIER_WORD[m]) out.push(MODIFIER_WORD[m]);
  if (r.source_state === 'offline') out.push('folder offline'); else if (r.source_state === 'missing') out.push('folder gone');
  return out;
}
function dataSub(r) {
  if (r.t2_files_waiting) return F.plural(r.t2_files_waiting, 'file') + ' waiting';
  if (r.read_data === 'absent') return '';
  if (r.t2_samples && r.t2_samples_with_data != null && r.t2_samples_with_data < r.t2_samples) return r.t2_samples_with_data + ' of ' + r.t2_samples + ' samples';
  if (r.read_data === 'partial') return r.coverage != null && r.coverage < 1 ? F.pct(r.coverage) + ' read so far' : 'read in part';
  return POPULATION_WORD[r.t2_population] || 'reads linked';
}

function prepare(ix) {
  for (const r of ix.rows) {
    const names = []; for (const s of ix.samplesByRun.get(r.id) || []) for (const n of [s.name, s.alias]) if (n && !names.includes(n)) names.push(n);
    r.samples = names; r.hay = hayOf(r);
  }
  return makeEngine(ix.rows, Date.parse(ix.written));
}
// A facet reads one value of a row, or several: a run with samples of two projects counts under each of them
// (so the Project facet, "Group by project" and the Projects page all give the same numbers).
const projectWord = v => { const p = index.value && index.value.projectById.get(v); return p && p.name != null ? p.name : v; };
const FACETS = [
  ['health', 'Run health', r => r.health, v => VERDICT_WORD[v]], ['data', 'Read data', r => r.data, v => VERDICT_WORD[v]],
  ['assay', 'Assay', r => r.assay_group, v => ASSAY_WORD[v] || v], ['fctype', 'Flow cell type', r => r.product, v => v], ['kit', 'Kit', r => r.kit, v => v],
  ['device', 'Device', r => r.device_type, v => DEVICE_WORD[v] || v], ['project', 'Project', r => r.project_ids, projectWord],
];
const GROUPS = { none: 'None', project: 'Project', health: 'Run health', fctype: 'Flow cell type', kit: 'Kit', month: 'Month', fc: 'Flow cell' };
const GROUP_KEY = { project: r => r.project_ids.length ? r.project_ids.map(projectWord) : ['No project'], health: r => VERDICT_WORD[r.health], fctype: r => r.product || 'Flow cell type not known', kit: r => r.kit || 'Kit not known', month: r => F.monthKey(r.started, r.stub), fc: r => r.fc };
const fkey = (field, v) => field + ':' + String(v).replace(/[^A-Za-z0-9_.-]/g, '_');

function Facets({ eng, q, onQ, open, setOpen, all, skip }) {
  if (!open) return <aside class="fc-facets is-rail"><button type="button" class="fc-btn fc-facets__btn" data-fk="facets-open" aria-expanded="false" onClick={() => { setOpen(true); refocus('[data-fk="facets-hide"]'); }}>Filters</button></aside>;
  return <aside class="fc-facets" aria-label="Filters">
    <div class="fc-facets__head"><span class="fc-caps">Filters</span><button type="button" class="fc-btn fc-btn--ghost" data-fk="facets-hide" aria-expanded="true" onClick={() => { setOpen(false); refocus('[data-fk="facets-open"]', '.fc-viewbar__filters'); }}>Hide</button></div>
    {FACETS.filter(f => f[0] !== skip).map(([field, name, get, word]) => {
      const counts = eng.facet(q, field, get), seen = new Map(); for (const r of all) { const v = get(r); for (const x of Array.isArray(v) ? v : [v]) if (x != null) seen.set(x, 0); }
      const vals = [...seen.keys()].map(v => [v, counts.get(v) || 0]).sort((a, b) => b[1] - a[1] || String(word(a[0]) || a[0]).localeCompare(String(word(b[0]) || b[0])));
      if (!vals.length) return null;
      return <fieldset class="fc-facet" key={field}><legend class="fc-caps">{name}</legend>{vals.map(([v, c]) => {
        const on = hasTerm(q, field, v), k = fkey(field, v);
        return <label key={k} class={'fc-facet__row' + (c ? '' : ' is-zero')}><span class="fc-facet__name"><input type="checkbox" data-fk={k} checked={on} onChange={() => { onQ(toggleTerm(q, field, v)); refocus('[data-fk="' + k + '"]'); }} /> <Mid text={word(v) || String(v)} tail={8} /></span><span class="fc-facet__n">{c}</span></label>;
      })}</fieldset>;
    })}
  </aside>;
}

function HeadCell({ c, sort, onSort, width }) {
  const key = sort.replace(/^-/, ''), on = c.sort && key === c.sort, dir = on ? (sort[0] === '-' ? 'descending' : 'ascending') : undefined;
  const inner = <><span>{c.label}{on ? <span aria-hidden="true">{dir === 'descending' ? ' ▾' : ' ▴'}</span> : null}</span>{c.sub ? <span class="fc-table__sub">{c.sub}</span> : null}</>;
  return <th scope="col" aria-sort={dir} title={c.title} class={(c.num ? 'fc-num ' : '') + (on ? 'is-on' : '') + ' c-' + c.id} style={{ width: width + 'px' }}>
    {c.id === 'tick' || c.id === 'peek' ? <span class="fc-sr">{c.id === 'tick' ? 'Select' : 'Peek'}</span> : c.sort ? <button type="button" onClick={() => onSort(c.sort)}>{inner}</button> : inner}</th>;
}

function Row({ r, cols, pos, ticked, peeked, widths, scale = 1 }) {
  const out = outputCell(r), stub = r.stub, words = healthWords(r), partial = r.read_data === 'partial';
  const open = e => { if (e.target.closest('a,button,input,label')) return; setParams({ peek: peeked ? '' : r.id }); };
  const togglePeek = () => { setParams({ peek: peeked ? '' : r.id }); if (peeked) refocus('[data-peek="' + r.id + '"]'); else refocus('.fc-peek .fc-close'); };
  const basesTitle = r.t2_passed_gb == null ? 'Not available yet: appears when the reads of this run have been read' : basesLabel(r.t2_population);
  return <tr class={(ticked ? 'is-selected ' : '') + (peeked ? 'is-marked' : '')} aria-rowindex={pos + 2} onClick={open}>
    {cols.map(c => {
      switch (c.id) {
        case 'tick': return <td class="c-tick"><input type="checkbox" data-k="tick" checked={ticked} aria-label={'Select ' + r.fc + ' ' + (r.run_name || '')} onChange={() => toggleSelected(r.id)} /></td>;
        case 'health': return <td class="c-health"><Badge v={r.health} soFar={runSoFar(r)} /><span class="fc-cell__sub" title={words.length ? words.join(' · ') : undefined}>{fitWords(words, widths.health / scale) || ' '}</span></td>;
        case 'data': return <td><Badge v={r.data} word={dataWord(r)} soFar={dataSoFar(r)} /><span class="fc-cell__sub">{dataSub(r) || ' '}</span></td>;
        case 'fc': return <th scope="row" class="c-fc"><a class="fc-mono fc-id" data-k="id" href={hrefRun(r.id)}>{r.fc}</a>{r.use_of > 1 ? <a class="fc-reuse" href={hrefFlowCell(r.fc)} title="Re-used flow cell: open its history">history</a> : null}<span class="fc-cell__sub" title={r.run_name || undefined}>{r.run_name ? F.midCut(r.run_name, Math.max(14, Math.floor((widths.fc / scale - 14) / 6.5)), 10) : (stub ? 'no report in the folder' : ' ')}</span></th>;
        case 'started': return <td><When iso={r.started} floating={stub} /><span class="fc-cell__sub">{r.hours != null ? F.hours(r.hours >= 10 ? Math.round(r.hours) : r.hours) + (ASSAY_SHORT[r.assay_group] ? ' · ' + ASSAY_SHORT[r.assay_group] : '') : ' '}</span></td>;
        case 'spark': return <td class="c-spark">{r.hours != null ? <Spark s={r.spark_output} hours={r.hours} early={r.hours < (r.planned_h || 96) - 2} /> : <span class="fc-muted">{F.DASH}</span>}</td>;
        case 'out': return <td class="fc-num">{out.a}<span class="fc-cell__sub">{out.b || ' '}</span></td>;
        case 'pores': return <td class="fc-num">{F.int(r.pores)}</td>;
        case 'occ': return <td class="fc-num">{r.occ_first != null && r.hours >= 1.5 ? F.pct(r.occ_first) : F.DASH}</td>;
        case 'passed': return <td class="fc-num" title={basesTitle}>{r.t2_passed_gb != null ? F.bases(r.t2_passed_gb * 1e9) : F.DASH}{r.t2_passed_gb != null && partial ? <span class="fc-cell__sub">partial</span> : null}</td>;
        case 'n50': return <td class="fc-num" title={r.t2_n50 == null ? 'Not available yet: appears when the reads of this run have been read' : undefined}>{r.t2_n50 != null ? F.bases(r.t2_n50) : F.DASH}{r.t2_n50 != null && partial ? <span class="fc-cell__sub">partial</span> : null}</td>;
        case 'why': return <td class="c-why">{stub ? null : <Differs r={r} />}</td>;
        case 'peek': return <td class="c-peek"><button type="button" class="fc-btn fc-btn--icon fc-peekbtn" data-k="peek" data-peek={r.id} aria-label={'Peek at ' + r.fc} aria-pressed={peeked} onClick={togglePeek}><span aria-hidden="true">{'›'}</span></button></td>;
      }
      return null;
    })}
  </tr>;
}

// The peek drawer: everything except the pore-scan chart comes from the index and is instant. The run's own file
// is asked for only after the drawer has stayed on one row for a moment, so walking down the table loads nothing.
function Peek({ id }) {
  const ix = index.value, r = ix.rows[ix.rowOf.get(id)];
  const [armed, setArmed] = useState(null);
  useEffect(() => { const t = setTimeout(() => setArmed(id), 150); return () => clearTimeout(t); }, [id]);
  const b = armed === id ? useBundle(id) : cachedBundle(id), state = bundleState(id);
  if (!r) return null;
  const grp = ix.groups.get(r.group), stub = r.stub;
  const h = r.outcome_h, metricId = r.outcome_metric || '', isReads = metricId.includes('reads'), partial = r.read_data === 'partial';
  const pts = D.sparkPoints(r.spark_output, r.hours), x = pts.map((_, i) => i * 3 + 1.5), y = pts.map(v => v * (r.spark_output_max_gb_per_3h || 0) / 3);
  const cur = grp && grp.n_runs >= 5 ? grp.curve_output_gb_per_3h : null;
  const sc = b ? D.scans(b) : null, inSel = selected.value.includes(id), na = 'not available yet';
  // Escape closes the drawer and gives the keyboard back to the button that opened it
  const close = () => { setParams({ peek: '' }); refocus('[data-peek="' + id + '"]', '.fc-table--hub [data-k="id"]'); };
  const onKey = e => { if (e.key === 'Escape') { e.preventDefault(); e.stopPropagation(); close(); } };
  return <aside class="fc-peek" aria-label={'Peek at ' + r.fc} onKeyDown={onKey}>
    <div class="fc-peek__head">
      <div><a class="fc-mono fc-peek__id" href={hrefRun(r.id)}>{r.fc}</a> <span class="fc-secondary"><When iso={r.started} floating={stub} />{r.hours != null ? ' · ' + F.hours(r.hours) : ''} {'·'} <a href={hrefFlowCell(r.fc)}>History</a></span></div>
      <button type="button" class="fc-close" aria-label="Close the peek drawer (Escape)" onClick={close}>Close</button>
    </div>
    <div class="fc-peek__name">{r.run_name || (stub ? 'No run name: only raw data was found' : 'No run name')}</div>
    <div class="fc-peek__badges"><span><Badge v={r.health} soFar={runSoFar(r)} /> Run health</span><span><Badge v={r.data} word={dataWord(r)} soFar={dataSoFar(r)} /> Read data{r.t2_files_waiting ? <span class="fc-secondary"> {'·'} {F.plural(r.t2_files_waiting, 'file')} waiting</span> : null}</span></div>
    <p class="fc-peek__headline">{r.headline}</p>
    <ul class="fc-reasons">{r.reasons.slice(0, 3).map((t, k) => <li><Glyph v={r.reason_levels[k] === 'fail' ? 'fail' : 'warn'} /><span>{t}</span></li>)}
      {r.detail ? <li class="is-na"><Glyph v="unknown" label="Note" /><span>{r.detail}</span></li> : null}</ul>
    {healthWords(r).length ? <p class="fc-peek__note fc-secondary">{healthWords(r).join(' · ')}</p> : null}
    {!stub ? <div class="fc-tiles fc-peek__tiles">
      <Tile label={h ? label(metricId) : 'Reads'} value={h ? (isReads ? F.countParts(r.outcome_value)[0] : F.basesParts(r.outcome_value)[0]) : F.countParts(r.reads)[0]} unit={h ? (isReads ? F.countParts(r.outcome_value)[1] : F.basesParts(r.outcome_value)[1]) : F.countParts(r.reads)[1]} />
      <Tile label="Pores at the first scan" value={F.int(r.pores)} />
      <Tile label="Occupancy, 1st hour" value={r.hours >= 1.5 ? F.pctNum(r.occ_first) : F.DASH} unit={r.hours >= 1.5 ? '%' : ''} />
      <Tile label={basesLabel(r.t2_population)} value={r.t2_passed_gb != null ? F.basesParts(r.t2_passed_gb * 1e9)[0] : na} unit={r.t2_passed_gb != null ? F.basesParts(r.t2_passed_gb * 1e9)[1] : ''} state={r.t2_passed_gb == null ? 'unknown' : null} caption={r.t2_passed_gb != null && partial ? 'partial' : null} />
      <Tile label="Read N50" value={r.t2_n50 != null ? F.basesParts(r.t2_n50)[0] : na} unit={r.t2_n50 != null ? F.basesParts(r.t2_n50)[1] : ''} state={r.t2_n50 == null ? 'unknown' : null} caption={r.t2_n50 != null && partial ? 'partial' : null} />
      <Tile label="Read quality" value={r.t2_q_med != null ? F.q(r.t2_q_med) : na} state={r.t2_q_med == null ? 'unknown' : null} caption={r.t2_q_med != null ? (partial ? 'median · partial' : 'median') : null} />
      {r.t2_full_length != null ? <Tile label="Full-length reads" value={F.pctNum(r.t2_full_length)} unit="%" caption="cDNA, from the pychopper statistics" title="cdna_full_length_frac" /> : null}
    </div> : null}
    {!stub && pts.length > 1 ? <ChartPanel cls="fc-chart--flat" title="Output per hour" reading={cur ? 'This run against the middle half of ' + grp.n_runs + ' runs of the same device, flow cell type and kit in the library today.' : 'Estimated bases per hour, in 3-hour steps.'} height={120}
      spec={{ form: 'xy', x: { domain: [0, 96], kind: 'time' }, y: { title: 'Gb per hour, estimated', fmt: v => String(+v.toFixed(2)) }, series: [{ name: 'This run', x, y, color: tk('--fc-this-run') }],
        band: cur ? { x: cur.median.map((_, i) => i * 3 + 1.5), lo: cur.q25.map(v => v == null ? null : v / 3), hi: cur.q75.map(v => v == null ? null : v / 3), mid: cur.median.map(v => v == null ? null : v / 3) } : null,
        cursor: { x }, readout: i => [['Hour ' + (i * 3) + ' to ' + (i * 3 + 3), ''], [y[i].toFixed(2) + ' Gb per hour', '']] }}
      table={() => ({ head: ['Hours', 'Gb per hour, estimated'], rows: y.map((v, i) => [(i * 3) + ' to ' + (i * 3 + 3), v.toFixed(2)]) })}
      legend={cur ? [{ name: 'This run', kind: 'line', color: tk('--fc-this-run') }, { name: 'lab today, middle half', kind: 'band' }] : null} /> : null}
    {sc ? <ChartPanel cls="fc-chart--flat" title="Pores per scan" reading={F.int(sc.pores[0]) + ' working pores at the first scan; ' + F.int(sc.pores[sc.n - 1]) + ' at the last.'} height={110}
      spec={{ form: 'xy', x: { domain: [0, Math.max(1, Math.ceil(r.hours || 1))], kind: 'time' }, y: { title: 'Working pores' }, series: [{ name: 'Working pores', type: 'step', dots: sc.n < 4, x: sc.t, y: Array.from(sc.pores), color: tk('--fc-scan-single-pore') }], cursor: { x: sc.t }, readout: i => [['Scan ' + (i + 1) + ', hour ' + sc.t[i].toFixed(1), ''], [F.int(sc.pores[i]), 'working pores']] }}
      table={() => ({ head: ['Scan', 'Hour', 'Working pores'], rows: Array.from(sc.pores, (v, i) => [String(i + 1), sc.t[i].toFixed(1), F.int(v)]) })} />
      : !stub ? <p class="fc-peek__note fc-muted">{state === 'failed' ? 'The data file of this run could not be read: summary only.' : 'Loading the pore-scan chart …'}</p> : null}
    <div class="fc-peek__actions">
      <a class="fc-btn fc-btn--primary" href={hrefRun(r.id)}>Open run</a>
      <button type="button" class="fc-btn" aria-pressed={inSel} onClick={() => toggleSelected(id)}>{inSel ? 'In compare' : 'Add to compare'}</button>
      <a class="fc-btn" href={hrefFlowCell(r.fc)}>History of this flow cell</a>
    </div>
  </aside>;
}

const PAGE = 400;       // rows drawn at once; a longer list grows on request
// The Runs table with its bars, filters, peek drawer and selection tray. rowFilter restricts it to a part of the
// library (a project page); embedded makes it part of a page that scrolls instead of filling the window.
export function RunsView({ rowFilter, embedded, label: tableLabel }) {
  const ix = index.value, eng = useMemo(() => prepare(ix), [ix]);
  const q = param('q', ''), sort = param('sort', '-started'), group = param('group', 'none'), peek = param('peek', ''), allCols = param('cols', '') === 'all';
  const sel = selected.value, scale = textScale.value;
  const wrap = useRef(null), [width, setWidth] = useState(1300), [facetsOpen, setFacetsOpen] = useState(() => !embedded && window.innerWidth >= 1440), [shown, setShown] = useState(PAGE);
  useEffect(() => { const ro = new ResizeObserver(() => { if (wrap.current) setWidth(wrap.current.clientWidth); }); ro.observe(wrap.current); if (wrap.current) setWidth(wrap.current.clientWidth); return () => ro.disconnect(); }, []);
  useScrollMemory(wrap);                               // Back from a run returns to the row one left, not to the top
  const all = rowFilter ? ix.rows.filter(rowFilter) : ix.rows, inScope = rowFilter ? new Set(all.map(r => r.id)) : null;
  const res = eng.run(q), found = inScope ? res.rows.filter(r => inScope.has(r.id)) : res.rows, fit = fitColumns(width - 2, scale, allCols);
  const key = sort.replace(/^-/, ''), dir = sort[0] === '-' ? -1 : 1, get = SORTS[key] || SORTS.started;
  const rows = found.slice().sort((a, b) => { const x = get(a), y = get(b), xn = x == null || x !== x, yn = y == null || y !== y; if (xn || yn) return xn && yn ? b.startedMs - a.startedMs : xn ? 1 : -1; return x < y ? -dir : x > y ? dir : b.startedMs - a.startedMs; });
  const onSort = k => setParams({ sort: key === k ? (dir === 1 ? '-' + k : k) : (DEFAULT_DIR[k] === -1 ? '-' + k : k) });
  const onQ = v => setParams({ q: v });
  const page = rows.slice(0, shown);
  let body, multi = 0;
  if (group !== 'none' && GROUP_KEY[group]) {     // grouping, one level deep; a run of two projects stands in both groups
    const gs = new Map(); for (const r of page) { const ks = [].concat(GROUP_KEY[group](r)); if (ks.length > 1) multi++; for (const k of ks) { if (!gs.has(k)) gs.set(k, []); gs.get(k).push(r); } }
    body = []; let pos = 0;
    for (const [k, list] of gs) {
      const roll = ['pass', 'warn', 'fail', 'incomplete', 'running'].map(v => [v, list.filter(r => r.health === v).length]).filter(x => x[1]);
      body.push(<tr key={'g:' + k} class="fc-table__group"><td colSpan={fit.cols.length}>{k} <span class="fc-table__roll">{'·'} {F.plural(list.length, 'run')} {'·'} {roll.map(([v, n]) => n + ' ' + VERDICT_WORD[v]).join(' · ')}</span></td></tr>);
      for (const r of list) body.push(<Row key={k + '/' + r.id} r={r} cols={fit.cols} pos={pos++} ticked={sel.includes(r.id)} peeked={peek === r.id} widths={fit.widths} scale={scale} />);
    }
  } else body = page.map((r, pos) => <Row key={r.id} r={r} cols={fit.cols} pos={pos} ticked={sel.includes(r.id)} peeked={peek === r.id} widths={fit.widths} scale={scale} />);
  const nWarn = all.filter(r => r.verdict === 'warn').length, nFail = all.filter(r => r.verdict === 'fail').length;
  // Up and Down move focus to the same control of the next row
  const onKey = e => {
    if (e.key !== 'ArrowDown' && e.key !== 'ArrowUp') return; const k = e.target.getAttribute && e.target.getAttribute('data-k'); if (!k) return;
    let tr = e.target.closest('tr'); do { tr = e.key === 'ArrowDown' ? tr.nextElementSibling : tr.previousElementSibling; } while (tr && !tr.querySelector('[data-k]'));
    if (tr) { e.preventDefault(); const t = tr.querySelector('[data-k="' + k + '"]'); t.focus(); t.scrollIntoView({ block: 'nearest' }); }
  };
  // on this screen: what a running scan is doing, and notices about a folder that was skipped. Notices about single
  // runs are on Home and on the pages they concern.
  const lines = loudNotices(ix).map(n => n.message); if (ix.progress.state === 'scanning' && ix.progress.message) lines.unshift(ix.progress.message);
  const peekShown = peek && (!inScope || inScope.has(peek));
  const hint = !found.length && all.length ? verdictHint(q) : null, bad = res.terms.filter(t => !t.ok);
  const hintN = hint ? eng.run(hint.query).rows.filter(r => !inScope || inScope.has(r.id)).length : 0;
  // the next thing to focus when a filter chip is removed: the chip that takes its place, else the search box
  const dropTerm = (raw, k) => { onQ(removeTerm(q, raw)); refocus('.fc-filterbar .fc-filter:nth-of-type(' + (k + 1) + ') .fc-filter__x', res.terms.length > 1 ? '.fc-filterbar .fc-filter:last-of-type .fc-filter__x' : embedded ? '.fc-viewbar__filters' : '.fc-appbar__search input'); };
  return <div class={'fc-hub' + (embedded ? ' fc-hub--embedded' : '')}>
    {embedded ? null : <h1 class="fc-sr">Runs</h1>}
    <div class="fc-viewbar">
      <button type="button" class="fc-btn fc-viewbar__filters" aria-expanded={facetsOpen} onClick={() => { setFacetsOpen(!facetsOpen); refocus(facetsOpen ? '.fc-viewbar__filters' : '[data-fk="facets-hide"]', '.fc-viewbar__filters'); }}>Filters</button>
      <label class="fc-ctrl"><span class="fc-caps">Group by</span><select class="fc-selectbox" value={group} onChange={e => setParams({ group: e.currentTarget.value === 'none' ? '' : e.currentTarget.value })}>{Object.keys(GROUPS).map(k => <option value={k}>{GROUPS[k]}</option>)}</select></label>
      <span class="fc-ctrl"><span class="fc-caps">Columns</span>{allCols
        ? <button type="button" class="fc-linkbtn fc-viewbar__cols" data-fk="cols" title={fit.wouldHide ? 'Every column is shown and the table scrolls sideways. Press to hide the columns that do not fit: ' + fit.wouldHideNames.join(', ') : 'Every column is shown'} onClick={() => { setParams({ cols: '' }); refocus('[data-fk="cols"]'); }}>All{fit.wouldHide ? ' · fit to the window' : ''}</button>
        : fit.hidden ? <button type="button" class="fc-linkbtn fc-viewbar__cols" data-fk="cols" title={'Hidden at this width: ' + fit.hiddenNames.join(', ') + '. Press to show every column; the table then scrolls sideways.'} onClick={() => { setParams({ cols: 'all' }); refocus('[data-fk="cols"]'); }}>Default {'·'} {fit.hidden} hidden {'·'} show all</button>
        : <span class="fc-viewbar__cols" title="Every default column is shown">Default</span>}</span>
      <span class="fc-viewbar__count" role="status" aria-live="polite">{found.length === all.length ? F.plural(all.length, 'run') : found.length + ' of ' + F.plural(all.length, 'run')} {'·'} {nWarn} with a warning, {nFail ? nFail + ' failed' : 'none failed'}{multi ? ' · ' + F.plural(multi, 'run') + ' with samples of two projects ' + (multi === 1 ? 'stands' : 'stand') + ' under each' : ''} {'·'} library <Updated /></span>
    </div>
    <div class="fc-filterbar">
      <span class="fc-caps">Filters</span>
      {res.terms.length === 0 ? <span class="fc-muted">none. Typed terms and ticked filters appear here, for example <button type="button" class="fc-linkbtn" onClick={() => onQ('status:warn')}>status:warn</button>, <button type="button" class="fc-linkbtn" onClick={() => onQ('pores:<5000 is:first-use')}>pores:&lt;5000 is:first-use</button>, <button type="button" class="fc-linkbtn" onClick={() => onQ('no:reads')}>no:reads</button></span> : null}
      {res.terms.map((t, k) => { const i = t.raw.indexOf(':'), isField = t.field && i > 0;
        return <span key={t.raw} class={'fc-filter' + (t.ok ? '' : ' is-unknown') + (t.neg ? ' is-not' : '')} title={t.error || undefined}>
          <span class="fc-filter__field">{isField ? t.raw.slice(t.neg ? 1 : 0, i + 1) : 'text:'}</span><span class="fc-filter__value">{isField ? t.raw.slice(i + 1) : t.raw.replace(/^-/, '')}</span>
          <button type="button" class="fc-filter__x" aria-label={'Remove the filter ' + t.raw} onClick={() => dropTerm(t.raw, k)}>Remove</button></span>; })}
      {bad.length ? <span class="fc-filterbar__error" role="status">{bad.map(t => t.raw + ': ' + t.error).join(' ')} A term that is not understood is struck through and never applied.</span> : null}
      {res.terms.length ? <button type="button" class="fc-btn fc-btn--ghost" onClick={() => onQ('')}>Clear all</button> : null}
      {lines.length ? <span class="fc-filterbar__notice fc-notice">{lines.join(' ')}</span> : null}
    </div>
    <div class="fc-hub__main">
      <Facets eng={eng} q={q} onQ={onQ} open={facetsOpen} setOpen={setFacetsOpen} all={all} skip={rowFilter ? 'project' : null} />
      <div class="fc-tablewrap fc-hub__table" ref={wrap}>
        <a class="fc-skip" href="#after-table">Skip past the table</a>
        <table class="fc-table fc-table--hub" style={{ width: Math.max(fit.total, fit.min) + 'px' }} aria-rowcount={rows.length + 1} aria-label={tableLabel || 'Runs'} onKeyDown={onKey}>
          <thead><tr aria-rowindex="1">{fit.cols.map(c => <HeadCell c={c} sort={sort} onSort={onSort} width={fit.widths[c.id]} />)}</tr></thead>
          <tbody>{body}</tbody>
        </table>
        {rows.length > shown ? <p class="fc-hub__more"><button type="button" class="fc-btn" onClick={() => setShown(shown + PAGE)}>Show {Math.min(PAGE, rows.length - shown)} more of {rows.length} runs</button></p> : null}
        {rows.length === 0 ? <div class="fc-empty fc-hub__empty">{all.length ? <><b>No run matches these filters.</b>
          {hint ? <span>A plain word is looked for in names (flow cell, run, sample, project, kit), not in verdicts. For the verdict type <button type="button" class="fc-linkbtn" onClick={() => onQ(hint.query)}>{hint.term}</button>{' (' + F.plural(hintN, 'run') + ')'}.</span>
            : <span>Remove a filter above, or clear them all.{bad.length ? ' ' + F.plural(bad.length, 'term') + (bad.length === 1 ? ' was' : ' were') + ' not understood and ' + (bad.length === 1 ? 'is' : 'are') + ' struck through above: ' + bad.map(t => t.raw).join(', ') + '.' : ''}</span>}
          {(ix.analyses || []).length ? <span class="fc-hub__anhint">This search looks through flow cell runs. The {F.plural(ix.analyses.length, 'single-cell and spatial analysis', 'single-cell and spatial analyses')} of the library are listed under <a href="#/analyses">Analyses</a>.</span> : null}</>
          : <><b>No runs here yet.</b><span>Runs appear when a folder that holds them has been scanned: <span class="fc-mono">flowcellar scan /folder/with/runs -L {ix.library.path || '<this library>'}</span></span></>}</div> : null}
        <span id="after-table" tabIndex={-1} />
      </div>
      {peekShown ? <Peek id={peek} /> : null}
    </div>
    {sel.length ? <div class="fc-tray fc-surface-inverse" role="region" aria-label="Selected runs">
      <b>{sel.length} selected</b><span class="fc-mono fc-tray__ids">{sel.map(id => id.split('_')[0]).join(' · ')}</span><span class="fc-grow" />
      {sel.length >= 2 && sel.length <= 8 ? <a class="fc-btn fc-btn--primary" href={'#/compare?items=' + sel.join(',')}>Compare</a> : <span class="fc-tray__hint">{sel.length < 2 ? 'Tick one more run to compare' : 'Compare shows up to 8 runs'}</span>}
      <button type="button" class="fc-btn fc-btn--inverse" onClick={() => { selected.value = []; }}>Clear</button>
    </div> : null}
  </div>;
}
export const Hub = () => <RunsView />;
