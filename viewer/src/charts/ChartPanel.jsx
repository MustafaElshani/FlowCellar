// The shell around a chart: a <figure> whose title and reading are visible text, a meta row (scope line,
// corner tag), a read-out of reserved height, the plot, the legend, and "Show as table".
// The canvas is aria-hidden; the keyboard cursor is one focusable element with role="slider" whose
// aria-valuetext is the read-out sentence. "Show as table" gives the numbers of a chart as a real table.
//
// A read-out never changes the size or the place of anything. Its box has a fixed height (styles: .fc-chart__readout,
// .fc-track__value), text that does not fit ends in an ellipsis, and the whole sentence is always in aria-valuetext
// and in the table. The pointer changes the text inside that box and nothing else on the page.
import { useRef, useEffect, useLayoutEffect, useState } from 'preact/hooks';
import { Chart } from './kit.js';
import { tk } from './theme.js';
import { themeTick } from '../core/state.js';

// write a read-out without going through the framework: [[value, name], ...] or a string
function writeReadout(node, parts) {
  if (!node) return '';
  node.textContent = '';
  if (parts == null) return '';
  if (typeof parts === 'string') { node.textContent = parts; return parts; }
  let said = '';
  parts.forEach((p, k) => {
    if (k) node.append(' · ');
    if (typeof p === 'string') { node.append(p); said += (k ? ', ' : '') + p; return; }
    const b = document.createElement('b'); b.textContent = p[0]; node.append(b);
    if (p[1]) node.append(' ' + p[1]);
    said += (k ? ', ' : '') + p[0] + (p[1] ? ' ' + p[1] : '');
  });
  return said;
}

// ---- the room a read-out needs in the label column of a track ----
// The read-outs of bins spread over the whole chart are written out and the longest is measured on a canvas (no
// layout is asked for): the column keeps that many lines free, once, before the pointer comes.
const partText = p => typeof p === 'string' ? p : p[0] + (p[1] ? ' ' + p[1] : '');
// values with their names keep all the room they need; free text (a message of the sequencer) may be cut after about two lines
function longestReadout(spec, idle, textMax) {
  const flat = parts => parts == null ? '' : typeof parts === 'string' ? parts.slice(0, textMax) : parts.map(p => typeof p === 'string' ? p.slice(0, textMax) : partText(p)).join(' · ');
  let best = flat(idle); const n = spec && spec.readout ? (spec.form === 'hbars' ? spec.rows.length : spec.cursor ? spec.cursor.x.length : 0) : 0;
  const see = i => { let t = ''; try { t = flat(spec.readout(i)); } catch (e) { /* a bin without a read-out */ } if (t.length > best.length) best = t; };
  for (let i = 0, step = Math.max(1, Math.floor(n / 60)); i < n; i += step) see(i);
  if (n) see(n - 1);
  return best;
}
let ruler = null;
function textWidth(text, px) { if (!ruler) ruler = document.createElement('canvas').getContext('2d'); ruler.font = '600 ' + px + 'px ' + (tk('--fc-font-sans') || 'sans-serif'); return ruler.measureText(text).width; }
// lines: 2 to 8 in a narrow column (two is what a time and one value take), 1 or more where the label lies across the page
function readoutLines(node, spec, idle) {
  const w = node.clientWidth; if (!(w > 0)) return 0;
  const px = parseFloat(getComputedStyle(node).fontSize) || 11, per = w * 0.88;
  return Math.max(w < 320 ? 2 : 1, Math.min(8, Math.ceil(textWidth(longestReadout(spec, idle, Math.round(2 * per / (px * 0.56))), px) / per)));
}

// On a map the cell under the pointer may be far from the read-out above the plot, so the same words also stand in a small
// note beside the cell. The note lies over the plot (it takes no room and moves nothing), it never leaves the plot, and it
// holds nothing that the read-out, aria-valuetext and the table do not hold: the keyboard cursor needs no note.
function placeTip(host, chart, parts) {
  let tip = host.querySelector(':scope > .fc-tip'); const at = parts == null ? null : chart.anchor();
  if (!at) { if (tip) tip.hidden = true; return; }
  if (!tip) { tip = document.createElement('div'); tip.className = 'fc-tooltip fc-tip'; tip.setAttribute('aria-hidden', 'true'); host.append(tip); }
  writeReadout(tip, parts); tip.hidden = false;
  const w = tip.offsetWidth, h = tip.offsetHeight, W = host.clientWidth, H = host.clientHeight, gap = at[2] + 8;
  let x = at[0] + gap, y = at[1] - gap - h; if (x + w > W) x = at[0] - gap - w; if (y < 0) y = at[1] + gap;
  tip.style.transform = 'translate(' + Math.round(Math.max(0, Math.min(W - w, x))) + 'px,' + Math.round(Math.max(0, Math.min(H - h, y))) + 'px)';
}

export function useChart(spec, readoutRef, idle) {
  const host = useRef(null), chart = useRef(null), specRef = useRef(spec);
  specRef.current = spec;
  const tick = themeTick.value;
  const wrap = s => ({ ...s, onHover: (i, byPointer) => {
    const cur = specRef.current, parts = i < 0 ? (idle || null) : (cur.readout ? cur.readout(i) : null);
    const said = writeReadout(readoutRef && readoutRef.current, parts);
    if (host.current) host.current.setAttribute('aria-valuetext', said || 'no value selected');
    if (host.current && chart.current && cur.form === 'map') placeTip(host.current, chart.current, i >= 0 && byPointer ? parts : null);
    if (cur.onCursor) cur.onCursor(i);
  } });
  useEffect(() => {
    if (!chart.current) chart.current = new Chart(host.current, wrap(spec)); else chart.current.update(wrap(spec));
  });
  useEffect(() => { if (readoutRef && readoutRef.current && chart.current && chart.current.cur < 0) writeReadout(readoutRef.current, idle || null); });
  useEffect(() => () => { if (chart.current) { chart.current.destroy(); chart.current = null; } }, []);
  const onKey = e => {
    const c = chart.current; if (!c) return;
    const d = { ArrowRight: 1, ArrowDown: 1, ArrowLeft: -1, ArrowUp: -1, PageDown: 10, PageUp: -10, Home: -Infinity, End: Infinity }[e.key];
    if (d == null) { if (e.key === 'Escape') c.setCursor(-1, true); return; }
    e.preventDefault(); c.move(d);
  };
  return { host, onKey, tick };
}

export function Plot({ spec, height, label, readoutRef, idle, aspect, maxWidth }) {
  const { host, onKey } = useChart(spec, readoutRef, idle);
  const style = aspect ? { aspectRatio: aspect, maxWidth: maxWidth ? maxWidth + 'px' : 'none' } : { height: height + 'px' };
  // A chart without a read-out (a map of many small cells) is a picture: it takes no keyboard focus, since the
  // arrow keys would do nothing there, and it says where its numbers are.
  if (!spec.readout) return <div ref={host} class="fc-plot fc-plot--still" style={style} role="img" aria-label={label + '. A picture without a value read-out: the numbers are under “Show as table”.'} />;
  return <div ref={host} class="fc-plot" style={style} tabIndex={0} role="slider" aria-label={label} aria-valuemin={0} aria-valuemax={100} aria-valuenow={0} aria-valuetext="no value selected" onKeyDown={onKey} />;
}

export function Legend({ items, hidden, onToggle }) {
  if (!items || !items.length) return null;
  return <div class="fc-legend">{items.map(it => {
    const off = hidden && hidden.has(it.name);
    const sw = it.kind === 'line' ? <span class="fc-sw fc-sw--line" style={{ borderTopColor: it.color, borderTopStyle: it.dash ? 'dashed' : 'solid' }} />
      : it.kind === 'band' ? <span class="fc-sw fc-sw--band" style={it.color ? { background: 'color-mix(in srgb, ' + it.color + ' 22%, transparent)', boxShadow: 'none' } : null} /> : it.kind === 'none' ? <span class="fc-sw fc-sw--none" />
      : it.kind === 'limit' ? <span class="fc-sw fc-sw--line fc-sw--limit" style={{ borderTopColor: it.color }} /> : <span class="fc-sw" style={{ background: it.color }} />;
    return onToggle && it.toggle !== false
      ? <button type="button" class={'fc-legend__item fc-legend__btn' + (off ? ' is-hidden' : '')} aria-pressed={!off} onClick={() => onToggle(it.name)}>{sw}{it.label || it.name}</button>
      : <span class="fc-legend__item">{sw}{it.label || it.name}</span>;
  })}</div>;
}

export function ScaleLegend({ lo, hi, unit, none }) {
  const steps = [1, 2, 3, 4, 5, 6, 7, 8, 9];
  return <div class="fc-scale"><span>{lo}</span><span class="fc-scale__ramp">{steps.map(i => <i style={{ background: tk('--fc-seq-' + i) }} />)}</span><span>{hi}{unit ? ' ' + unit : ''}</span>
    {none ? <span class="fc-legend__item fc-scale__none"><span class="fc-sw fc-sw--none" />{none}</span> : null}</div>;
}

export function DataTable({ data, caption }) {
  if (!data) return null;
  const LIMIT = 3200, rows = data.rows.length > LIMIT ? data.rows.filter((_, i) => i % Math.ceil(data.rows.length / LIMIT) === 0) : data.rows;
  return <div class="fc-tablewrap fc-charttable" tabIndex={0} role="group" aria-label={caption}>
    <table class="fc-table"><caption class="fc-sr">{caption}</caption>
      <thead><tr>{data.head.map((h, k) => <th scope="col" class={k ? 'fc-num' : ''}>{h}</th>)}</tr></thead>
      <tbody>{rows.map(r => <tr>{r.map((c, k) => k ? <td class="fc-num">{c}</td> : <th scope="row">{c}</th>)}</tr>)}</tbody></table>
    {data.note ? <p class="fc-muted fc-charttable__note">{data.note}</p> : null}
    {rows.length < data.rows.length ? <p class="fc-muted fc-charttable__note">Every {Math.ceil(data.rows.length / LIMIT)}th row of {data.rows.length} is listed here.</p> : null}
  </div>;
}

// A stand-alone chart panel.
// lines: how many lines the read-out keeps free above the plot (2, or 1 for a chart whose read-out is always short).
export function ChartPanel({ id, title, reading, tag, scope, height = 220, spec, legend, table, idle, how, cls, footer, legendToggle, children, aspect, maxWidth, evidence, lines = 2 }) {
  const ro = useRef(null), [showTable, setShowTable] = useState(false), [showHow, setShowHow] = useState(false), [hidden, setHidden] = useState(() => new Set());
  const toggle = legendToggle ? name => { const n = new Set(hidden); n.has(name) ? n.delete(name) : n.add(name); setHidden(n); } : null;
  const shown = spec && hidden.size && spec.series ? { ...spec, series: spec.series.filter(q => !hidden.has(q.name)) } : spec;
  // --fc-plot-max: the most the plot is drawn wide (a map of cells of a fixed most size); wide layouts size the panel by it
  return <figure class={'fc-chart' + (cls ? ' ' + cls : '') + (evidence ? ' is-evidence' : '')} id={id} style={aspect && maxWidth ? { '--fc-plot-max': maxWidth + 'px' } : undefined}>
    <figcaption>
      <div class="fc-chart__head"><span class="fc-chart__title">{title}</span>
        <span class="fc-chart__menu">
          {how ? <button type="button" class="fc-btn" aria-expanded={showHow} onClick={() => setShowHow(!showHow)}>How to read</button> : null}
          {table ? <button type="button" class="fc-btn" aria-expanded={showTable} onClick={() => setShowTable(!showTable)}>{showTable ? 'Hide table' : 'Show as table'}</button> : null}
        </span></div>
      {reading ? <p class="fc-chart__reading">{reading}</p> : null}
      {scope || tag ? <div class="fc-chart__meta"><span class="fc-chart__scope">{scope || ''}</span>{tag ? <span class={'fc-tag fc-chart__tag' + (/^ESTIMATED/.test(tag) ? ' fc-tag--estimate' : '')}>{tag}</span> : null}</div> : null}
      {showHow ? <p class="fc-chart__how">{how}</p> : null}
    </figcaption>
    {spec ? <div class={'fc-chart__readout' + (lines === 1 ? '' : ' fc-chart__readout--2')} aria-hidden="true"><span class="fc-ro" ref={ro} /></div> : null}
    {spec ? <Plot spec={shown} height={height} label={title} readoutRef={ro} idle={idle} aspect={aspect} maxWidth={maxWidth} /> : null}
    {children}
    <Legend items={legend} hidden={hidden} onToggle={toggle} />
    {footer}
    {showTable && table ? <DataTable data={table()} caption={title + ', as a table'} /> : null}
  </figure>;
}

// One track of a timeline: a label column that prints the value under the cursor, and the plot.
export function Track({ id, title, sub, height = 110, spec, idle, legend, evidence, tag, table, how, tools }) {
  const ro = useRef(null), [showTable, setShowTable] = useState(false), [showHow, setShowHow] = useState(false);
  // the label column keeps the lines its longest read-out needs: set before the first paint, and again when the
  // window (and with it the column) changes width
  const latest = useRef(null); latest.current = { spec, idle };
  useLayoutEffect(() => {
    const fit = () => { const node = ro.current; if (!node) return; const n = readoutLines(node, latest.current.spec, latest.current.idle); if (n && node.dataset.lines !== String(n)) { node.dataset.lines = n; node.style.setProperty('--fc-ro-lines', n); } };
    fit(); addEventListener('resize', fit); return () => removeEventListener('resize', fit);
  });
  return <div class={'fc-track' + (evidence ? ' is-evidence' : '')} id={id} role="group" aria-label={title}>
    <div class="fc-track__label">
      <div class="fc-track__title">{title}</div>
      {sub ? <div class="fc-track__sub">{sub}</div> : null}
      {tag ? <div><span class="fc-tag">{tag}</span></div> : null}
      <div class="fc-track__value" ref={ro} aria-hidden="true" />
      {table || how || tools ? <div class="fc-track__menu">
        {tools}
        {how ? <button type="button" class="fc-linkbtn" aria-expanded={showHow} onClick={() => setShowHow(!showHow)}>How to read</button> : null}
        {table ? <button type="button" class="fc-linkbtn" aria-expanded={showTable} onClick={() => setShowTable(!showTable)}>{showTable ? 'Hide table' : 'Show as table'}</button> : null}
      </div> : null}
    </div>
    <div class="fc-track__plot">
      <Plot spec={spec} height={height} label={title} readoutRef={ro} idle={idle} />
      {legend ? <Legend items={legend} /> : null}
      {showHow && how ? <p class="fc-chart__how fc-track__how">{how}</p> : null}
      {showTable && table ? <DataTable data={table()} caption={title + ', as a table'} /> : null}
    </div>
  </div>;
}

export function Empty({ title, children }) { return <div class="fc-empty"><b>{title}</b><span>{children}</span></div>; }
