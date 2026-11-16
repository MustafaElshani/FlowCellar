// The spatial map of an analysis: its own canvas component (the chart kit is not involved).
//
// One picture with a switch ("Colour by"), not five pictures. Square bins (Visium HD) are painted as one pixel per
// bin on a hidden canvas and scaled without smoothing; spots (Visium) are drawn as discs on their staggered grid.
// The picture is painted when the data, the channel, the size or the theme changes, never for a pointer move: the
// pointer only moves a cursor on a second canvas and rewrites the text of a read-out box of fixed height, so
// nothing on the page changes size or place. The keyboard moves the same cursor (arrow keys); the read-out
// sentence is the element's aria-valuetext, and "Show as table" gives the numbers behind the picture.
import { useRef, useEffect, useMemo, useState } from 'preact/hooks';
import { tk, seqLUT } from '../charts/theme.js';
import { ScaleLegend, DataTable } from '../charts/ChartPanel.jsx';
import { themeTick } from '../core/state.js';
import { pixelRatio, onPaper, printing } from '../core/print.js';
import * as F from '../core/format.js';
import { decodeChannel, channelText, channelRange, CHANNEL_WORD, CHANNEL_ORDER, levelValue } from './sc.js';

const SQRT3_2 = Math.sqrt(3) / 2;
// every channel of a grid, decoded once per bundle
export function useChannels(grid) {
  return useMemo(() => { const out = {}; if (grid) for (const id of Object.keys(grid.channels)) out[id] = decodeChannel(grid.channels[id]); return out; }, [grid]);
}
export const channelIds = grid => CHANNEL_ORDER.filter(id => grid && grid.channels[id]);

// where a cell is drawn, in units of one cell: square bins tile the plane; spots sit on rows that are 0.866 apart,
// every second row half a spot to the right
function geometry(grid) {
  const hex = grid.layout && grid.layout.kind === 'hex', rows = grid.rows, cols = grid.cols;
  return hex ? { hex, rows, cols, w: cols + 0.5, h: (rows - 1) * SQRT3_2 + 1, x: (r, c) => c + (r % 2 ? 0.5 : 0) + 0.5, y: r => r * SQRT3_2 + 0.5 }
    : { hex, rows, cols, w: cols, h: rows, x: (r, c) => c + 0.5, y: r => r + 0.5 };
}
function colourOf(dec, ch, id, lut, i, none) {
  const q = dec.q[i];
  if (dec.nodata && q === 255) return none;
  if (id === 'image') return [q, q, q];
  return lut[Math.max(0, Math.min(255, Math.round(q / dec.hi * 255)))];
}

function paint(cv, grid, G, dec, ch, id, opts) {
  const dpr = pixelRatio(), W = cv.clientWidth, H = cv.clientHeight; if (!(W > 0 && H > 0)) return null;
  cv.width = Math.round(W * dpr); cv.height = Math.round(H * dpr);
  const g = cv.getContext('2d'); g.setTransform(dpr, 0, 0, dpr, 0, 0); g.clearRect(0, 0, W, H);
  const s = Math.min(W / G.w, H / G.h), ox = (W - G.w * s) / 2, oy = (H - G.h * s) / 2, lut = seqLUT();
  const none = tk('--fc-chart-none') || tk('--fc-border'), rgbNone = null, tissue = opts.tissue;
  if (G.hex) {
    const r = Math.max(1, s * 0.5 * (grid.layout.bin_um ? 55 / grid.layout.bin_um : 0.55) * 1.25);
    for (let row = 0; row < G.rows; row++) for (let c = 0; c < G.cols; c++) {
      const i = row * G.cols + c, q = dec.q[i], x = ox + G.x(row, c) * s, y = oy + G.y(row) * s;
      if (dec.nodata && q === 255) { g.fillStyle = none; g.globalAlpha = 0.35; } else { const col = colourOf(dec, ch, id, lut, i); g.fillStyle = 'rgb(' + col[0] + ',' + col[1] + ',' + col[2] + ')'; g.globalAlpha = 1; }
      g.beginPath(); g.arc(x, y, r, 0, Math.PI * 2); g.fill();
    }
    g.globalAlpha = 1;
    if (tissue && opts.outline) { g.strokeStyle = tk('--fc-text'); g.lineWidth = Math.max(1, r * 0.22);
      for (let row = 0; row < G.rows; row++) for (let c = 0; c < G.cols; c++) { const i = row * G.cols + c; if (!(levelValue(tissue.ch, tissue.dec.q[i]) >= 0.5)) continue; g.beginPath(); g.arc(ox + G.x(row, c) * s, oy + G.y(row) * s, r + g.lineWidth * 0.6, 0, Math.PI * 2); g.stroke(); } }
    return { s, ox, oy, r };
  }
  // square bins: one pixel per bin, scaled without smoothing
  const off = document.createElement('canvas'); off.width = G.cols; off.height = G.rows;
  const og = off.getContext('2d'), img = og.createImageData(G.cols, G.rows), d = img.data;
  for (let i = 0; i < G.rows * G.cols; i++) { const q = dec.q[i], k = i * 4;
    if (dec.nodata && q === 255) { d[k + 3] = 0; continue; }
    const col = colourOf(dec, ch, id, lut, i, rgbNone); d[k] = col[0]; d[k + 1] = col[1]; d[k + 2] = col[2]; d[k + 3] = 255; }
  og.putImageData(img, 0, 0);
  g.fillStyle = tk('--fc-chart-surface') || tk('--fc-surface'); g.fillRect(ox, oy, G.w * s, G.h * s);
  // On paper a browser may write this drawing into the PDF as it is, and a reader of the PDF would then smooth a picture of one pixel
  // per bin into a blur: there the bins are first made whole squares of device pixels (at most 4,096 a side), so they stay squares.
  let src = off; const k = printing.peek() ? Math.min(Math.ceil(s * dpr), Math.floor(4096 / Math.max(G.cols, G.rows))) : 1;
  if (k > 1) { src = document.createElement('canvas'); src.width = G.cols * k; src.height = G.rows * k; const sg = src.getContext('2d'); sg.imageSmoothingEnabled = false; sg.drawImage(off, 0, 0, src.width, src.height); }
  g.imageSmoothingEnabled = false; g.drawImage(src, ox, oy, G.w * s, G.h * s);
  if (tissue && opts.outline) {         // the edge of the tissue: between a bin that is mostly under tissue and one that is not
    const t = i => levelValue(tissue.ch, tissue.dec.q[i]) >= 0.5;
    g.strokeStyle = tk('--fc-text'); g.lineWidth = 1.25; g.beginPath();
    for (let r = 0; r < G.rows; r++) for (let c = 0; c < G.cols; c++) { const a = t(r * G.cols + c);
      if (c + 1 < G.cols && a !== t(r * G.cols + c + 1)) { const x = Math.round(ox + (c + 1) * s) + 0.5; g.moveTo(x, oy + r * s); g.lineTo(x, oy + (r + 1) * s); }
      if (r + 1 < G.rows && a !== t((r + 1) * G.cols + c)) { const y = Math.round(oy + (r + 1) * s) + 0.5; g.moveTo(ox + c * s, y); g.lineTo(ox + (c + 1) * s, y); } }
    g.stroke();
  }
  return { s, ox, oy };
}

// the cell under a point of the canvas, or -1
function hit(G, L, x, y) {
  const u = (x - L.ox) / L.s, v = (y - L.oy) / L.s;
  if (!G.hex) { const c = Math.floor(u), r = Math.floor(v); return r >= 0 && r < G.rows && c >= 0 && c < G.cols ? r * G.cols + c : -1; }
  const r = Math.round((v - 0.5) / SQRT3_2); if (r < 0 || r >= G.rows) return -1;
  const c = Math.round(u - 0.5 - (r % 2 ? 0.5 : 0)); return c >= 0 && c < G.cols ? r * G.cols + c : -1;
}
function drawCursor(cv, G, L, i) {
  const dpr = pixelRatio(), W = cv.clientWidth, H = cv.clientHeight;
  if (cv.width !== Math.round(W * dpr) || cv.height !== Math.round(H * dpr)) { cv.width = Math.round(W * dpr); cv.height = Math.round(H * dpr); }
  const g = cv.getContext('2d'); g.setTransform(dpr, 0, 0, dpr, 0, 0); g.clearRect(0, 0, W, H);
  if (i < 0 || !L) return;
  const r = Math.floor(i / G.cols), c = i % G.cols, x = L.ox + G.x(r, c) * L.s, y = L.oy + G.y(r) * L.s;
  g.strokeStyle = tk('--fc-text'); g.lineWidth = 2; g.fillStyle = 'transparent';
  if (G.hex) { g.beginPath(); g.arc(x, y, (L.r || L.s / 2) + 3, 0, Math.PI * 2); g.stroke(); }
  else { const h = Math.max(L.s, 6); g.strokeRect(Math.round(x - h / 2) - 1, Math.round(y - h / 2) - 1, Math.round(h) + 2, Math.round(h) + 2); }
}

// What one cell holds, as the parts of a read-out: every channel of the cell, the coloured one first.
export function cellParts(grid, chans, i, id, unit) {
  const r = Math.floor(i / grid.cols), c = i % grid.cols, um = grid.layout && grid.layout.bin_um;
  const where = grid.layout && grid.layout.kind === 'hex' ? 'spot row ' + (r + 1) + ', ' + (c + 1) : (um ? Math.round(um) + ' µm ' : '') + 'bin ' + (r + 1) + ', ' + (c + 1);
  const ids = [id].concat(channelIds(grid).filter(k => k !== id && k !== 'image'));
  return [where].concat(ids.map(k => channelText(k, chans[k].value(i), unit)));
}
function write(node, parts) { if (!node) return ''; node.textContent = ''; parts.forEach((p, k) => { if (k) node.append(' · '); if (k === 1) { const b = document.createElement('b'); b.textContent = p; node.append(b); } else node.append(p); }); return parts.join(', '); }

// The picture alone (also the small one on a card): no controls. `onCell` is told the cell under the cursor.
// `link`: a Set shared by the pictures of one grid. The cell pointed at in one of them is marked in the others too.
export function MapCanvas({ grid, chans, id, outline, unit, readoutRef, idle, label, still, maxWidth, link }) {
  const base = useRef(null), over = useRef(null), host = useRef(null), L = useRef(null), cur = useRef(-1), raf = useRef(0), ptr = useRef(null), follow = useRef(null);
  const G = useMemo(() => geometry(grid), [grid]);
  const tick = themeTick.value;
  const say = i => { const parts = i < 0 ? [idle || ''] : cellParts(grid, chans, i, id, unit); const said = write(readoutRef && readoutRef.current, parts); if (host.current && !still) host.current.setAttribute('aria-valuetext', i < 0 ? 'no cell selected' : said); };
  const setCursor = i => { if (i === cur.current) return; cur.current = i; if (over.current) drawCursor(over.current, G, L.current, i); say(i); if (link) for (const f of link) if (f !== follow.current) f(i); };
  // what another picture of the same grid calls: the cursor is drawn here, the read-out stays with the picture pointed at
  useEffect(() => { if (!link) return undefined; const f = i => { if (i === cur.current) return; cur.current = i; if (over.current) drawCursor(over.current, G, L.current, i); if (host.current) host.current.setAttribute('aria-valuetext', i < 0 ? 'no cell selected' : cellParts(grid, chans, i, id, unit).join(', ')); };
    follow.current = f; link.add(f); return () => { link.delete(f); }; }, [link, grid, chans, id, unit]);
  useEffect(() => {
    const dec = chans[id]; if (!dec || !base.current) return;
    // the cleanup of an effect runs after paint: an observer can still call once after the canvas has left the page
    const draw = () => { if (!base.current) return; L.current = paint(base.current, grid, G, dec, grid.channels[id], id, { outline, tissue: chans.tissue && id !== 'tissue' ? { dec: chans.tissue, ch: grid.channels.tissue } : null }); if (over.current) drawCursor(over.current, G, L.current, cur.current); say(cur.current); };
    draw();
    const ro = typeof ResizeObserver === 'function' ? new ResizeObserver(draw) : null; if (ro) ro.observe(base.current);
    const offPaper = onPaper(draw);                    // on its way to paper, and back, the picture is painted at once
    return () => { if (ro) ro.disconnect(); offPaper(); };
  }, [grid, id, outline, tick, chans]);
  useEffect(() => () => cancelAnimationFrame(raf.current), []);
  const style = { aspectRatio: String(G.w / G.h), maxWidth: maxWidth ? maxWidth + 'px' : undefined };
  if (still) return <div class="fc-spmap__plot fc-spmap__plot--still" style={style} role="img" aria-label={label}><canvas ref={base} class="fc-spmap__cv" aria-hidden="true" /></div>;
  // the pointer is read once per frame; only the cursor layer is drawn again
  const onMove = e => { ptr.current = [e.clientX, e.clientY]; if (raf.current) return; raf.current = requestAnimationFrame(() => { raf.current = 0; const p = ptr.current, el = over.current; if (!p || !el || !L.current) return; const b = el.getBoundingClientRect(); setCursor(hit(G, L.current, p[0] - b.left, p[1] - b.top)); }); };
  const onLeave = () => { ptr.current = null; setCursor(-1); };
  const onKey = e => {
    const step = e.shiftKey ? 5 : 1, d = { ArrowRight: [0, step], ArrowLeft: [0, -step], ArrowDown: [step, 0], ArrowUp: [-step, 0] }[e.key];
    if (e.key === 'Escape') { setCursor(-1); return; }
    if (e.key === 'Home' || e.key === 'End') { e.preventDefault(); setCursor(e.key === 'Home' ? 0 : G.rows * G.cols - 1); return; }
    if (!d) return; e.preventDefault();
    const i = cur.current < 0 ? Math.floor(G.rows / 2) * G.cols + Math.floor(G.cols / 2) : cur.current, r = Math.max(0, Math.min(G.rows - 1, Math.floor(i / G.cols) + d[0])), c = Math.max(0, Math.min(G.cols - 1, i % G.cols + d[1]));
    setCursor(r * G.cols + c);
  };
  return <div ref={host} class="fc-spmap__plot" style={style} tabIndex={0} role="slider" aria-label={label + '. Arrow keys move over the cells.'} aria-valuemin={0} aria-valuemax={G.rows * G.cols - 1} aria-valuenow={0}
    aria-valuetext="no cell selected" onKeyDown={onKey} onBlur={() => setCursor(-1)}>
    <canvas ref={base} class="fc-spmap__cv" aria-hidden="true" />
    <canvas ref={over} class="fc-spmap__cv fc-spmap__over" aria-hidden="true" onPointerMove={onMove} onPointerLeave={onLeave} />
  </div>;
}

// median, smallest and largest of a channel over the cells that have a value (and lie under tissue, where that is known)
function channelStats(grid, chans, id) {
  const dec = chans[id], t = chans.tissue && id !== 'tissue' ? chans.tissue : null, vals = [];
  for (let i = 0; i < dec.q.length; i++) { if (dec.nodata && dec.q[i] === 255) continue; if (t && !(levelValue(grid.channels.tissue, t.q[i]) >= 0.5)) continue; vals.push(dec.value(i)); }
  vals.sort((a, b) => a - b); const n = vals.length;
  return n ? { n, min: vals[0], med: vals[Math.floor((n - 1) / 2)], max: vals[n - 1] } : { n: 0 };
}

// The map with its switch, its legend, its read-out and its table. `exact`: rows [label, value] of the numbers that
// are counted from the full-resolution data (never read off the picture); `side`: what stands beside the picture.
export function SpatialMap({ grid, unit, title, reading, exact, side, note, initial }) {
  const chans = useChannels(grid), ids = channelIds(grid), ro = useRef(null), link = useMemo(() => new Set(), [grid]);
  const [id, setId] = useState(initial && ids.includes(initial) ? initial : ids[0]), [outline, setOutline] = useState(true), [showTable, setShowTable] = useState(false);
  const cur = ids.includes(id) ? id : ids[0], hasTissue = !!grid.channels.tissue, hex = grid.layout && grid.layout.kind === 'hex';
  // On a wide screen a second picture of the same slide stands beside the first, coloured by something else (the
  // tissue image where there is one): the two share the cursor, so one bin is read in both at once.
  const [id2, setId2] = useState(null), others = ids.filter(k => k !== cur), cur2 = others.includes(id2) ? id2 : others.includes('image') ? 'image' : others[0];
  const range2 = cur2 ? channelRange(cur2, grid.channels[cur2], chans[cur2]) : null;
  const range = channelRange(cur, grid.channels[cur], chans[cur]), um = grid.layout && grid.layout.bin_um;
  const what = hex ? F.int(grid.rows * grid.cols) + ' spots, each drawn at its place on the slide' : grid.rows + ' × ' + grid.cols + ' squares of ' + Math.round(um) + ' µm';
  const idle = 'Point at the picture, or use the arrow keys, to read one ' + (hex ? 'spot' : 'bin') + '.';
  const table = () => ({ head: ['What', 'Smallest', 'Median', 'Largest', (hex ? 'Spots' : 'Bins') + ' counted'],
    rows: ids.filter(k => k !== 'image').map(k => { const s = channelStats(grid, chans, k); return [CHANNEL_WORD[k] + (hasTissue && k !== 'tissue' ? ', under tissue' : ''), s.n ? channelText(k, s.min, unit) : F.DASH, s.n ? channelText(k, s.med, unit) : F.DASH, s.n ? channelText(k, s.max, unit) : F.DASH, F.int(s.n)]; })
      .concat((exact || []).map(r => [r[0], '', r[1], '', 'exact'])),
    note: 'The picture stores one of 255 levels per ' + (hex ? 'spot' : 'bin of ' + Math.round(um) + ' µm') + ': its rows are read off the picture and are approximate. Rows marked exact are counted from the full-resolution data.' });
  return <figure class={'fc-chart fc-chart--map fc-spmap' + (cur2 ? ' fc-spmap--two' : '')}>
    <figcaption>
      <div class="fc-chart__head"><span class="fc-chart__title">{title}</span>
        <span class="fc-chart__menu"><button type="button" class="fc-btn" aria-expanded={showTable} onClick={() => setShowTable(!showTable)}>{showTable ? 'Hide table' : 'Show as table'}</button></span></div>
      {reading ? <p class="fc-chart__reading">{reading}</p> : null}
      <div class="fc-chart__meta"><span class="fc-chart__scope">{what}{hasTissue && outline && cur !== 'tissue' ? ' · tissue outline drawn' : ''}{printing.value ? ' · coloured by ' + F.lower(CHANNEL_WORD[cur]) : ''}</span><span class="fc-tag fc-chart__tag">A PICTURE · exact totals are in the tiles</span></div>
    </figcaption>
    <div class="fc-spmap__body">
      <div class="fc-spmap__main fc-spmap__first">
        <div class="fc-chart__readout fc-chart__readout--2" aria-hidden="true"><span class="fc-ro" ref={ro} /></div>
        <MapCanvas grid={grid} chans={chans} id={cur} outline={outline} unit={unit} readoutRef={ro} idle={idle} link={link} label={title + ', coloured by ' + F.lower(CHANNEL_WORD[cur])} />
        <div class="fc-chart__foot">{cur === 'image' ? <span class="fc-muted">The tissue image the pipeline registered to the slide, in grey.</span>
          : <ScaleLegend lo={range.lo} hi={range.hi} unit={range.unit} none={chans[cur].nodata ? 'no data' : null} />}</div>
      </div>
      {cur2 ? <div class="fc-spmap__main fc-spmap__second">
        <div class="fc-chart__readout fc-chart__readout--2 fc-spmap__beside"><label>Beside it, the same slide coloured by <select class="fc-selectbox" value={cur2} onChange={e => setId2(e.currentTarget.value)}>{others.map(k => <option key={k} value={k}>{CHANNEL_WORD[k]}</option>)}</select></label></div>
        <MapCanvas grid={grid} chans={chans} id={cur2} outline={outline} unit={unit} readoutRef={ro} idle={idle} link={link} label={title + ', a second picture coloured by ' + F.lower(CHANNEL_WORD[cur2])} />
        <div class="fc-chart__foot">{cur2 === 'image' ? <span class="fc-muted">The tissue image the pipeline registered to the slide, in grey.</span>
          : <ScaleLegend lo={range2.lo} hi={range2.hi} unit={range2.unit} none={chans[cur2].nodata ? 'no data' : null} />}</div>
      </div> : null}
      <div class="fc-spmap__side">
        <div class="fc-spmap__pick" role="radiogroup" aria-label="Colour by">
          <span class="fc-caps">Colour by</span>
          {ids.map(k => <button key={k} type="button" role="radio" aria-checked={String(k === cur)} class={'fc-spmap__opt' + (k === cur ? ' is-on' : '')} onClick={() => setId(k)}>{CHANNEL_WORD[k]}</button>)}
        </div>
        {hasTissue ? <label class="fc-spmap__check"><input type="checkbox" checked={outline} onChange={e => setOutline(e.currentTarget.checked)} /> Tissue outline</label> : null}
        {side}
      </div>
    </div>
    {note ? <p class="fc-chart__note fc-muted">{note}</p> : null}
    {showTable ? <DataTable data={table()} caption={title + ', as a table'} /> : null}
  </figure>;
}
