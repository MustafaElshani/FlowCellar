// fc-charts: FlowCellar's own small canvas kit. No third-party chart library.
//
// One canvas per chart plus a hover layer; device pixel ratio capped at 2 on screen, and the scale of paper while the
// page is laid out for print (core/print.js: every chart paints itself again then). Colours come only from the
// tokens (theme.js). Read-outs, legends and tables are DOM elements built by the shell (ChartPanel.jsx);
// the kit reports the hovered bin through spec.onHover(i, byPointer) and never draws a floating box.
//
// Hover: the pointer is read once per frame (requestAnimationFrame), only the hover layer is drawn again, and of
// that layer only the part the last cursor covered is cleared. The chart itself is never painted for a pointer move.
// Anywhere on the canvas reads the nearest bin, row, segment or cell: the read-out goes back to its idle text
// only when the pointer leaves the chart.
//
// Forms: 'xy' (line, area, step, histogram, outline, stacked area, stacked bars at x positions, heat cells,
// with limit lines, bands and labelled marks), 'hbars' (ranked bars, also read fate), 'hstack' (one stacked bar
// per row), 'map' (flow cell map), 'strip' (events on a time axis), 'ring' (parts of a whole with their labels
// beside them).
//
// Added for the alignment charts, usable by any section:
//   'ring'            { form: 'ring', parts: [{ label, value, color, text }], center: { value, label }, readout: i => ... }
//                     at most six parts; each part is labelled beside the ring with its name and `text` (count and
//                     share), joined to it by a line, so no part is told by colour alone
//   a quartile band   on any line of 'xy': series: [{ x, y, band: { lo: [...], hi: [...] } }] fills between lo and hi
//                     in the colour of the line
//   a genome track    genomeTrack({ contigs: [{ name, first, windows }], values, means, ... }) returns an 'xy' spec:
//                     one bar per window, a hairline at every chromosome boundary, the names under the axis, the mean
//                     of each chromosome as a line (x.bounds and hsegs do this and can be used on any 'xy')
//   labels in a bar   'hstack' with segLabels: true writes part.label inside every segment that has room for it
import { tk, font, onData, alpha, seqColor, themeGen } from './theme.js';
import { pixelRatio, onPaper, printing } from '../core/print.js';
import { niceTicks, tickLabel } from '../core/format.js';

const TAU = Math.PI * 2;
export const bus = new Map();                                // sync key -> Set<Chart>
const el = (tag, cls) => { const e = document.createElement(tag); if (cls) e.className = cls; return e; };
const crisp = v => Math.round(v) + 0.5;
const fin = v => typeof v === 'number' && isFinite(v);
const clamp = (v, lo, hi) => Math.max(lo, Math.min(hi, v));
function lowerBound(a, v) { let lo = 0, hi = a.length; while (lo < hi) { const m = (lo + hi) >> 1; if (a[m] < v) lo = m + 1; else hi = m; } return lo; }
export function nearest(a, v) { if (!a || !a.length) return -1; let i = lowerBound(a, v); if (i > 0 && (i >= a.length || v - a[i - 1] < a[i] - v)) i--; return i; }

// ---- x ticks by kind ----
function timeTicks(d0, d1) {
  const span = d1 - d0, step = span > 60 ? 24 : span > 30 ? 12 : span > 14 ? 6 : span > 7 ? 2 : span > 2.5 ? 1 : span > 1.2 ? 0.5 : span > 0.4 ? 0.25 : 1 / 12;
  const out = [];
  for (let v = Math.ceil(d0 / step - 1e-9) * step; v <= d1 + 1e-9; v += step) {
    const r = Math.round(v * 12) / 12;
    out.push({ v: r, label: step >= 1 ? Math.round(r) + ' h' : (span <= 2.5 ? Math.round(r * 60) + ' min' : r + ' h'), major: true });
  }
  if (step === 24 && span <= 110) for (let v = Math.ceil(d0 / 12) * 12; v <= d1; v += 12) if (v % 24) out.push({ v, label: '', major: false });
  return out;
}
const SI = (v, unit) => v >= 1e6 ? v / 1e6 + ' M' + unit : v >= 1e3 ? v / 1e3 + ' k' + unit : v + ' ' + unit;
function logTicks(d0, d1, kind) {
  const out = [];
  for (let e = Math.floor(Math.log10(d0)); e <= Math.ceil(Math.log10(d1)); e++) for (const m of [1, 2, 5]) {
    const v = m * Math.pow(10, e); if (v < d0 * 0.999 || v > d1 * 1.001) continue;
    let label;
    if (kind === 'length') label = m === 1 ? SI(v, 'b') : String(v >= 1e6 ? v / 1e6 : v >= 1e3 ? v / 1e3 : v);
    else if (kind === 'duration') label = m === 1 ? (v >= 3600 ? v / 3600 + ' h' : v >= 60 && v % 60 === 0 ? v / 60 + ' min' : tickLabel(v) + ' s') : tickLabel(v);
    else label = m === 1 ? tickLabel(v) : '';
    out.push({ v, label, major: m === 1 });
  }
  return out;
}

export class Chart {
  constructor(host, spec) {
    this.host = host; host.classList.add('fck');
    this.cv = el('canvas', 'fck-base'); this.ov = el('canvas', 'fck-over');
    this.cv.setAttribute('aria-hidden', 'true'); this.ov.setAttribute('aria-hidden', 'true');
    host.append(this.cv, this.ov);
    this.g = this.cv.getContext('2d'); this.og = this.ov.getContext('2d');
    this.cur = -1; this.curX = null; this.spec = spec; this.W = 0; this.H = 0; this.drag = null;
    this.ptr = null; this.raf = 0; this.dirty = null;
    this.ov.addEventListener('pointermove', e => this.onPointer([e.clientX, e.clientY]));
    this.ov.addEventListener('pointerleave', () => this.onPointer(false));
    this.ov.addEventListener('pointerdown', e => this.onDown(e));
    this.ov.addEventListener('pointerup', e => this.onUp(e));
    this.ov.addEventListener('dblclick', () => { if (this.spec.onBrush) this.spec.onBrush(null); });
    this.ro = new ResizeObserver(() => { if (this.measure()) this.paint(); });
    this.ro.observe(host);
    this.offPaper = onPaper(() => this.refresh());
    this.join(); this.measure(); this.paint();
  }
  join() { const k = this.spec.sync; if (this.syncKey === k) return; this.leave(); this.syncKey = k; if (k) { if (!bus.has(k)) bus.set(k, new Set()); bus.get(k).add(this); } }
  leave() { if (this.syncKey && bus.has(this.syncKey)) bus.get(this.syncKey).delete(this); this.syncKey = null; }
  measure() {
    const r = this.host.getBoundingClientRect(), W = Math.max(40, Math.round(r.width)), H = Math.max(20, Math.round(r.height)), dpr = pixelRatio();
    if (W === this.W && H === this.H && dpr === this.dpr) return false;
    this.W = W; this.H = H; this.dpr = dpr;
    for (const c of [this.cv, this.ov]) { c.width = W * dpr; c.height = H * dpr; }
    return true;
  }
  // a new spec while a cursor stands (the keyboard cursor, or the pointer resting on the chart): the cursor stays
  // where it is when the new spec still has that bin, and the read-out is written again from the new numbers
  update(spec) {
    this.spec = spec; this.join(); this.measure(); this.paint();
    if (this.cur < 0) return;
    if (this.cur >= this.count()) { this.cur = -1; this.curX = null; }
    this.drawCursor(); if (spec.onHover) spec.onHover(this.cur);
  }
  count() { const s = this.spec; return s.form === 'hbars' ? s.rows.length : s.form === 'hstack' ? Infinity : s.form === 'ring' ? s.parts.length : s.form === 'map' ? s.values.length : s.cursor ? s.cursor.x.length : 0; }
  retheme() { this.paint(); }
  // paint again when the size, the scale or the theme is no longer the one of the last painting (the page went to paper, or came back)
  refresh() { if (this.measure() || this.gen !== themeGen()) { this.paint(); this.drawCursor(); } }
  destroy() { this.leave(); this.ro.disconnect(); this.offPaper(); if (this.raf) cancelAnimationFrame(this.raf); this.raf = 0; this.host.textContent = ''; }

  paint() {
    const g = this.g, s = this.spec; this.gen = themeGen();
    g.setTransform(this.dpr, 0, 0, this.dpr, 0, 0);
    g.clearRect(0, 0, this.W, this.H);
    g.fillStyle = tk('--fc-chart-surface'); g.fillRect(0, 0, this.W, this.H);
    g.font = font(); g.lineJoin = 'round'; g.lineCap = 'butt'; g.setLineDash([]); g.globalAlpha = 1;
    if (s.form === 'hbars') this.paintHBars(g);
    else if (s.form === 'hstack') this.paintHStack(g);
    else if (s.form === 'ring') this.paintRing(g);
    else if (s.form === 'map') this.paintMap(g);
    else this.paintXY(g);
    this.og.setTransform(this.dpr, 0, 0, this.dpr, 0, 0); this.og.clearRect(0, 0, this.W, this.H); this.dirty = null;
  }

  // =========================== cartesian ===========================
  prepXY() {
    const s = this.spec, g = this.g, W = this.W, H = this.H, x = s.x || {}, y = s.y || {};
    const d0 = x.domain[0], d1 = x.domain[1], xlog = !!x.log, fx = xlog ? Math.log10 : (v => v);
    // ---- y domain ----
    let y0, y1;
    const inX = v => v >= d0 - 1e-9 && v <= d1 + 1e-9;
    if (y.domain) { y0 = y.domain[0]; y1 = y.domain[1]; }
    else {
      let lo = Infinity, hi = -Infinity; const see = v => { if (fin(v)) { if (v < lo) lo = v; if (v > hi) hi = v; } };
      const ser = s.series || [];
      if (s.stack) { const n = s.stack.x.length; for (let i = 0; i < n; i++) { if (!inX(s.stack.x[i])) continue; let t = 0; for (const q of ser) t += q.y[i] || 0; see(t); } see(0); }
      else for (const q of ser) {
        if (q.parts) { for (let i = 0; i < q.x.length; i++) { if (!inX(q.x[i])) continue; let t = 0; for (const p of q.parts) t += p.y[i] || 0; see(t); } continue; }
        const n = q.y.length;
        for (let i = 0; i < n; i++) { const xv = q.edges ? q.edges[i + 1] : q.x[i]; const xa = q.edges ? q.edges[i] : xv; if (xv < d0 || xa > d1) continue; see(q.y[i]); if (q.band) { see(q.band.lo[i]); see(q.band.hi[i]); } }
      }
      for (const h of s.hsegs || []) see(h.y);
      if (s.band) for (let i = 0; i < s.band.x.length; i++) if (inX(s.band.x[i])) { see(s.band.hi[i]); see(s.band.lo[i]); }
      for (const l of s.limits || []) if (l.fit !== false) see(l.y);
      for (const p of s.points || []) if (inX(p.x)) see(p.y);
      if (s.hband) { see(s.hband.y0); see(s.hband.y1); }
      if (!fin(lo)) { lo = 0; hi = 1; }
      if (y.log) { y1 = hi > 0 ? hi : 1; let minPos = Infinity; for (const q of ser) for (const v of q.y || []) if (v > 0 && v < minPos) minPos = v; y0 = Math.max(Math.pow(10, Math.ceil(Math.log10(y1)) - 5), isFinite(minPos) ? Math.pow(10, Math.floor(Math.log10(minPos))) : 0); if (!(y1 > y0)) y1 = y0 * 10; }
      else if (y.zero === false) { const pad = (hi - lo || Math.abs(hi) || 1) * 0.12; y0 = lo - pad; y1 = hi + pad; }
      else { y0 = Math.min(0, lo); y1 = hi > 0 ? hi : 1; }
    }
    let yt;
    if (y.ticks) yt = y.ticks;
    else if (y.log) { yt = []; for (let e = Math.floor(Math.log10(y0)); e <= Math.ceil(Math.log10(y1)); e++) yt.push(Math.pow(10, e)); y0 = yt[0]; y1 = yt[yt.length - 1]; }
    else { yt = niceTicks(y0, y1, H < 130 ? 3 : H < 220 ? 4 : 5); if (!y.domain) { const st = yt.length > 1 ? yt[1] - yt[0] : 1; if (yt[yt.length - 1] < y1 - st * 1e-6) yt.push(+(yt[yt.length - 1] + st).toPrecision(12)); y1 = Math.max(y1, yt[yt.length - 1]); if (y.zero === false) y0 = Math.min(y0, yt[0]); } }
    const yfmt = y.fmt || (y.log ? (v => v >= 1e9 ? +(v / 1e9).toPrecision(3) + ' G' : v >= 1e6 ? +(v / 1e6).toPrecision(3) + ' M' : v >= 1e3 ? +(v / 1e3).toPrecision(3) + ' k' : String(+v.toPrecision(3))) : tickLabel);
    g.font = font();
    let lw = 0; if (y.labels !== false) for (const v of yt) lw = Math.max(lw, g.measureText(yfmt(v)).width);
    const ml = s.ml != null ? s.ml : Math.max(30, Math.ceil(lw) + 12), mr = s.mr != null ? s.mr : 14;
    // a title of an axis that is longer than the chart is wide takes two lines (cut at the space nearest its middle) instead of running out of it
    const twoLines = (t, room) => { if (g.measureText(t).width <= room) return [t]; const mid = t.length / 2; let cut = -1; for (let i = 0; i < t.length; i++) if (t[i] === ' ' && (cut < 0 || Math.abs(i - mid) < Math.abs(cut - mid))) cut = i; return cut > 0 ? [t.slice(0, cut), t.slice(cut + 1)] : [t]; };
    const yTitle = y.title ? twoLines(y.title, W - ml - 4) : null;
    const mt = (s.mt != null ? s.mt : (y.title ? 22 : 8)) + (yTitle ? 13 * (yTitle.length - 1) : 0) + (s.events ? 12 : 0) + (s.vlines && s.vlines.length ? 12 : 0);
    const xTitle = x.title && x.labels !== false ? twoLines(x.title, W - 8) : null;
    const mb = x.labels === false ? 5 : (xTitle ? 23 + 13 * xTitle.length : 22);
    const px = ml, py = mt, pw = Math.max(10, W - ml - mr), ph = Math.max(10, H - mt - mb);
    const fy = y.log ? Math.log10 : (v => v), fy0 = fy(y0), fy1 = fy(y1), fx0 = fx(d0), fx1 = fx(d1);
    const L = this.L = { px, py, pw, ph, d0, d1, y0, y1, yt, yfmt, xTitle, yTitle,
      X: v => px + (fx(v) - fx0) / (fx1 - fx0 || 1) * pw,
      Y: v => py + ph - (fy(y.log ? Math.max(v, y0) : v) - fy0) / (fy1 - fy0 || 1) * ph,
      invX: p => { const t = fx0 + (p - px) / pw * (fx1 - fx0); return xlog ? Math.pow(10, t) : t; } };
    L.xt = x.ticks || (x.kind === 'time' ? timeTicks(d0, d1) : xlog ? logTicks(d0, d1, x.kind) : niceTicks(d0, d1, pw < 360 ? 4 : 7).map(v => ({ v, label: (x.fmt || tickLabel)(v), major: true })));
    return L;
  }
  paintXY(g) {
    const s = this.spec, L = this.prepXY(), { px, py, pw, ph, X, Y } = L, x = s.x || {}, y = s.y || {};
    const surface = tk('--fc-chart-surface'), tick = tk('--fc-chart-tick'), by = py + ph;
    // gridlines: horizontal hairlines only, never on stacks and heat maps
    g.lineWidth = 1;
    if (!s.stack && !s.heat && y.grid !== false) { g.strokeStyle = tk('--fc-chart-grid'); g.beginPath(); for (const v of L.yt) { const yy = crisp(Y(v)); if (yy > by - 1) continue; g.moveTo(px, yy); g.lineTo(px + pw, yy); } g.stroke(); }
    if (x.bounds) { g.strokeStyle = tk('--fc-chart-grid'); g.lineWidth = 1; g.beginPath(); for (const v of x.bounds) { const xx = crisp(X(v)); if (xx <= px || xx >= px + pw) continue; g.moveTo(xx, py); g.lineTo(xx, by + 4); } g.stroke(); }
    g.save(); g.beginPath(); g.rect(px, py - 1, pw, ph + 2); g.clip();
    // lab-typical band: grey fill with a 1 px edge line, median hairline
    if (s.band) {
      const b = s.band, n = b.x.length; let started = false; g.beginPath();
      const idx = []; for (let i = 0; i < n; i++) if (fin(b.lo[i]) && fin(b.hi[i])) idx.push(i);
      if (idx.length > 1) {
        for (const i of idx) { const xx = X(b.x[i]), yy = Y(b.hi[i]); if (!started) { g.moveTo(xx, yy); started = true; } else g.lineTo(xx, yy); }
        for (let k = idx.length - 1; k >= 0; k--) g.lineTo(X(b.x[idx[k]]), Y(b.lo[idx[k]]));
        g.closePath(); g.fillStyle = tk('--fc-chart-typical'); g.fill();
        g.strokeStyle = tk('--fc-chart-typical-edge'); g.lineWidth = 1;
        for (const key of ['hi', 'lo']) { g.beginPath(); idx.forEach((i, k) => { const xx = X(b.x[i]), yy = Y(b[key][i]); k ? g.lineTo(xx, yy) : g.moveTo(xx, yy); }); g.stroke(); }
      }
      if (b.mid) { g.strokeStyle = tk('--fc-chart-typical-line'); g.lineWidth = 1; g.beginPath(); let st = false; for (let i = 0; i < n; i++) { if (!fin(b.mid[i])) continue; const xx = X(b.x[i]), yy = Y(b.mid[i]); st ? g.lineTo(xx, yy) : g.moveTo(xx, yy); st = true; } g.stroke(); }
    }
    // heat cells
    if (s.heat) { const h = s.heat; for (let k = 0; k < h.i.length; k++) { const x0 = X(h.xe[h.i[k]]), x1 = X(h.xe[h.i[k] + 1]), ya = Y(h.ye[h.j[k]]), yb = Y(h.ye[h.j[k] + 1]); g.fillStyle = seqColor(h.t[k]); g.fillRect(x0, yb, Math.max(1, x1 - x0 + 0.5), Math.max(1, ya - yb + 0.5)); } }
    // stacked area: full-strength fills, a 1 px line in the surface colour between bands
    if (s.stack) {
      const xs = s.stack.x, n = xs.length, ser = s.series, cum = new Float64Array(n), prev = new Float64Array(n), XP = new Float64Array(n);
      for (let i = 0; i < n; i++) XP[i] = X(xs[i]);
      this.bandsGeom = [];
      for (const q of ser) {
        prev.set(cum); for (let i = 0; i < n; i++) cum[i] += q.y[i] || 0;
        g.beginPath(); g.moveTo(XP[0], Y(cum[0])); for (let i = 1; i < n; i++) g.lineTo(XP[i], Y(cum[i]));
        for (let i = n - 1; i >= 0; i--) g.lineTo(XP[i], Y(prev[i])); g.closePath(); g.fillStyle = q.color; g.fill();
        if (n < pw / 3) { g.strokeStyle = surface; g.lineWidth = 1; g.beginPath(); g.moveTo(XP[0], Y(cum[0])); for (let i = 1; i < n; i++) g.lineTo(XP[i], Y(cum[i])); g.stroke(); }
        if (q.label) { // a label inside the band, where the band stays thick over a stretch (a single spike is not room)
          g.font = font(); const tw = g.measureText(q.label).width + 10, k = Math.max(1, Math.ceil(tw / Math.max(0.5, (XP[n - 1] - XP[0]) / Math.max(1, n - 1))));
          let best = -1, bt = 0, by0 = 0;
          for (let i = 0; i + k < n; i += Math.max(1, Math.floor(k / 3))) { if (XP[i] < px + 4 || XP[i + k] > px + pw - 4) continue; let top = -Infinity, bot = Infinity; for (let j = i; j <= i + k; j++) { top = Math.max(top, Y(cum[j])); bot = Math.min(bot, Y(prev[j])); } const t = bot - top; if (t > bt + 0.5) { bt = t; best = i; by0 = (top + bot) / 2; } }
          if (best >= 0 && bt >= 14) this.bandsGeom.push({ label: q.label, x: XP[best] + 5, y: by0, color: q.color });
        }
      }
      g.strokeStyle = alpha(surface, 0.35); g.lineWidth = 1; g.beginPath(); for (const v of L.yt) { if (v <= L.y0 || v >= L.y1) continue; const yy = crisp(Y(v)); g.moveTo(px, yy); g.lineTo(px + pw, yy); } g.stroke();
      g.font = font(); g.textBaseline = 'middle'; g.textAlign = 'left';
      for (const b of this.bandsGeom) { g.fillStyle = onData(b.color); g.fillText(b.label, b.x, b.y); }
    } else for (const q of s.series || []) this.paintSeries(g, q, L, !!s.heat);
    // a level over a stretch of the axis (the mean of one chromosome): a line on a ground of the surface colour
    for (const h of s.hsegs || []) { if (!fin(h.y)) continue; const xa = Math.max(px, X(h.x0)) + 1, xb = Math.min(px + pw, X(h.x1)) - 1, yy = crisp(Y(Math.min(h.y, L.y1))); if (xb - xa < 2) continue;
      g.strokeStyle = surface; g.lineWidth = 3.5; g.beginPath(); g.moveTo(xa, yy); g.lineTo(xb, yy); g.stroke(); g.strokeStyle = h.color || tk('--fc-chart-annotation'); g.lineWidth = 1.5; g.beginPath(); g.moveTo(xa, yy); g.lineTo(xb, yy); g.stroke(); }
    g.restore();
    // baseline
    g.strokeStyle = tk('--fc-chart-axis'); g.lineWidth = 1; g.beginPath(); g.moveTo(px, crisp(by)); g.lineTo(px + pw, crisp(by)); g.stroke();
    // events strip: pore scans as ticks on a hairline, messages as marks above it
    if (s.strip) {
      const cy = crisp(py + ph * 0.66);
      g.strokeStyle = tk('--fc-chart-axis'); g.lineWidth = 1; g.beginPath(); g.moveTo(px, cy); g.lineTo(px + pw, cy); g.stroke();
      g.strokeStyle = alpha(tk('--fc-chart-annotation'), 0.7); g.beginPath();
      for (const xv of s.strip.scans) { const xx = crisp(X(xv)); if (xx < px || xx > px + pw) continue; g.moveTo(xx, cy - 5); g.lineTo(xx, cy + 5); } g.stroke();
      g.textAlign = 'center'; g.textBaseline = 'middle';
      for (const m of s.strip.marks) {
        const xx = X(m.x); if (xx < px - 1 || xx > px + pw + 1) continue; const yy = cy - 14;
        if (m.kind === 'warning' || m.kind === 'error') { g.font = font('600'); g.fillStyle = m.kind === 'error' ? tk('--fc-status-fail-fg') : tk('--fc-status-warn-fg'); g.fillText(m.kind === 'error' ? '\u2715' : '\u25B2', xx, yy); g.font = font(); }
        else { g.strokeStyle = tk('--fc-chart-annotation'); g.lineWidth = 1.5; g.beginPath(); g.arc(xx, yy, 3, 0, TAU); g.stroke(); }
      }
    }
    // fixed target range: two hairlines, labelled once
    if (s.hband) { const hb = s.hband; g.strokeStyle = alpha(tk('--fc-chart-annotation'), 0.55); g.lineWidth = 1; g.beginPath(); for (const v of [hb.y0, hb.y1]) { const yy = crisp(Y(v)); g.moveTo(px, yy); g.lineTo(px + pw, yy); } g.stroke(); if (hb.label) this.halo(g, hb.label, px + pw - 4, Y(hb.y1) - 7, 'right', tk('--fc-chart-annotation')); }
    // limits: dashed lines in the status outline colour, labelled in words, haloed
    for (const l of s.limits || []) {
      const yy = crisp(Y(l.y)); if (yy < py - 1 || yy > by + 1) continue;
      const col = l.kind === 'fail' ? tk('--fc-status-fail-border') : l.kind === 'target' ? tk('--fc-chart-annotation') : tk('--fc-status-warn-border');
      g.strokeStyle = tk('--fc-chart-limit-halo'); g.lineWidth = 3; g.setLineDash([]); g.beginPath(); g.moveTo(px, yy); g.lineTo(px + pw, yy); g.stroke();
      g.strokeStyle = col; g.lineWidth = 1; g.setLineDash(l.kind === 'fail' ? [2, 2] : [5, 3]); g.beginPath(); g.moveTo(px, yy); g.lineTo(px + pw, yy); g.stroke(); g.setLineDash([]);
      if (l.label) {
        let above = 0, below = 0; for (const q of s.series || []) if (q.x && q.y && !q.parts) { const m = Math.min(q.y.length, 40); for (let i = 0; i < m; i++) { const v = q.y[i]; if (!fin(v) || q.x[i] < L.d0 || q.x[i] > L.d1) continue; if (v >= l.y) above++; else below++; } }
        const under = (above > below && by - yy > 16) || yy - py < 14;
        // the words of a limit stand on an opaque ground: a line that crosses them does not go through the letters
        g.font = font(); const lw2 = g.measureText(l.label).width, ty = yy + (under ? 8 : -7); g.fillStyle = tk('--fc-chart-surface'); g.fillRect(px + 3, ty - 7, lw2 + 6, 14);
        this.halo(g, l.label, px + 6, yy + (under ? 8 : -7), 'left', l.kind === 'fail' ? tk('--fc-status-fail-fg') : l.kind === 'target' ? tk('--fc-chart-annotation') : tk('--fc-status-warn-fg'));
      }
    }
    // vertical rules (run boundaries on the flow cell page)
    for (const v of s.vrules || []) { const xx = crisp(X(v.x)); if (xx < px || xx > px + pw) continue; g.strokeStyle = tk('--fc-chart-annotation'); g.lineWidth = 1; g.beginPath(); g.moveTo(xx, py - 4); g.lineTo(xx, by); g.stroke(); if (v.label) this.halo(g, v.label, xx + 4, py + 7, 'left', tk('--fc-chart-annotation')); }
    // N50, median: a 1 px line from the baseline to above the plot, label at the top
    let vl = 0;
    for (const v of s.vlines || []) {
      const xx = crisp(X(v.x)); if (xx < px || xx > px + pw) continue;
      g.strokeStyle = v.strong ? tk('--fc-chart-annotation') : alpha(tk('--fc-chart-annotation'), 0.55); g.lineWidth = 1; g.beginPath(); g.moveTo(xx, v.strong ? py - 12 : py - 2); g.lineTo(xx, by); g.stroke();
      g.font = font(v.strong ? '600' : ''); const right = v.side !== 'left'; const tw = g.measureText(v.label).width;
      const flip = right ? xx + 5 + tw > px + pw : xx - 5 - tw < px;
      this.halo(g, v.label, (right !== flip) ? xx + 5 : xx - 5, v.strong ? py - 6 : py + 7 + 12 * (vl++ % 2), (right !== flip) ? 'left' : 'right', tk('--fc-chart-annotation'), v.strong ? '600' : '');
    }
    // marked points with haloed value labels (24, 48, 72 h); a label that would sit on another one is left out
    const placed = [], free = (x0, y0, w, h) => { if (x0 < 0 || x0 + w > this.W || y0 < 0 || y0 + h > by + 2) return false; for (const r of placed) if (x0 < r[0] + r[2] + 7 && x0 + w + 7 > r[0] && y0 < r[1] + r[3] + 3 && y0 + h + 3 > r[1]) return false; return true; };
    for (const p of s.points || []) {
      const xx = X(p.x), yy = Y(p.y); if (xx < px - 1 || xx > px + pw + 1) continue;
      g.fillStyle = surface; g.beginPath(); g.arc(xx, yy, 6, 0, TAU); g.fill(); g.fillStyle = p.color || tk('--fc-this-run'); g.beginPath(); g.arc(xx, yy, 4, 0, TAU); g.fill();
      if (!p.label) continue;
      g.font = font('600'); const tw = g.measureText(p.label).width, th = 12;
      // left of and above the point first (the line of a curve that rises comes from below left); the label stands
      // on an opaque ground, so where no free place is clear of the line the line does not go through the letters
      for (const [lx, ly] of [[xx - tw - 6, yy - 18], [xx + 6, yy - 18], [xx + 8, yy - 6], [xx - tw - 6, yy + 8], [xx + 6, yy + 8], [xx - tw - 8, yy - 6]]) { if (!free(lx, ly, tw, th)) continue; placed.push([lx, ly, tw, th]);
        g.fillStyle = alpha(surface, 0.9); g.fillRect(lx - 3, ly - 1, tw + 6, th + 2); this.halo(g, p.label, lx, ly + th / 2, 'left', tk('--fc-text'), '600'); break; }
    }
    // line-end letters (Compare): kept at least a line apart
    const ends = [];
    for (const q of s.series || []) if (q.end && q.x) { let i = q.y.length - 1; while (i >= 0 && (!fin(q.y[i]) || q.x[i] > L.d1)) i--; if (i < 0) continue; const xx = X(q.x[i]), yy = Y(q.y[i]); g.fillStyle = surface; g.beginPath(); g.arc(xx, yy, 6, 0, TAU); g.fill(); g.fillStyle = q.color; g.beginPath(); g.arc(xx, yy, 4, 0, TAU); g.fill(); ends.push({ xx, yy, ly: yy, text: q.end }); }
    ends.sort((p, q) => p.yy - q.yy);
    for (let k = 1; k < ends.length; k++) if (ends[k].ly - ends[k - 1].ly < 12) ends[k].ly = ends[k - 1].ly + 12;
    const over = ends.length ? ends[ends.length - 1].ly - (by - 4) : 0; if (over > 0) for (const e of ends) e.ly -= over;
    for (const e of ends) this.halo(g, e.text, Math.min(e.xx + 9, this.W - 9), e.ly, 'left', tk('--fc-text'), '600');
    // event ticks above the plot
    if (s.events) { g.strokeStyle = alpha(tk('--fc-chart-annotation'), 0.7); g.lineWidth = 1; g.beginPath(); for (const e of s.events) { const xx = crisp(X(e.x)); if (xx < px || xx > px + pw) continue; g.moveTo(xx, py - 11); g.lineTo(xx, py - 7); } g.stroke(); }
    // ---- axis labels ----
    g.font = font(); g.fillStyle = tick; g.textBaseline = 'middle'; g.textAlign = 'right';
    if (y.labels !== false) { let last = Infinity; for (const v of L.yt) { const yy = Y(v); if (yy < py - 2 || yy > by + 2 || last - yy < 11) continue; last = yy; g.fillText(L.yfmt(v), px - 6, yy); } }
    if (L.yTitle) { g.textAlign = 'left'; g.textBaseline = 'top'; L.yTitle.forEach((t, k) => g.fillText(t, Math.max(0, px), 2 + 13 * k)); }
    if (x.labels !== false) {
      g.textBaseline = 'top'; let lastR = -Infinity;
      const draw = major => { for (const t of L.xt) { if (!!t.major !== major || !t.label) continue; const xx = X(t.v); if (xx < px - 1 || xx > px + pw + 1) continue; const w = g.measureText(t.label).width; let lx = xx - w / 2; if (lx < 0) lx = 0; if (lx + w > this.W) lx = this.W - w; if (major) { if (lx < lastR + 8) continue; lastR = lx + w; } else { let clash = false; for (const o of taken) if (lx < o[1] + 6 && lx + w > o[0] - 6) clash = true; if (clash) continue; } taken.push([lx, lx + w]); g.textAlign = 'left'; g.fillText(t.label, lx, by + 6); } };
      const taken = []; draw(true); draw(false);
      if (x.tickMarks !== false) { g.strokeStyle = tk('--fc-chart-axis'); g.beginPath(); for (const t of L.xt) { const xx = crisp(X(t.v)); if (xx < px - 1 || xx > px + pw + 1) continue; g.moveTo(xx, by); g.lineTo(xx, by + (t.major ? 4 : 2)); } g.stroke(); }
      if (L.xTitle) { g.textAlign = 'center'; g.fillStyle = tick; L.xTitle.forEach((t, k) => { const half = g.measureText(t).width / 2 + 2; g.fillText(t, Math.max(half, Math.min(this.W - half, px + pw / 2)), by + 21 + 13 * k); }); }
    }
  }
  halo(g, text, x, y, align, color, weight) {
    g.font = font(weight || ''); g.textAlign = align; g.textBaseline = 'middle';
    g.lineWidth = 3; g.strokeStyle = tk('--fc-chart-surface'); g.lineJoin = 'round'; g.strokeText(text, x, y);
    g.fillStyle = color; g.fillText(text, x, y); g.font = font();
  }
  // cased: the line lies over heat cells, whose scale runs to the colour of the text in either theme. It is drawn on a
  // wider line in the colour of the chart surface (dark under a light line, light under a dark one), so it never
  // vanishes where the cells are densest.
  paintSeries(g, q, L, cased) {
    const { X, Y, py, ph } = L, by = py + ph, type = q.type || 'line', color = q.color || tk('--fc-this-run');
    g.setLineDash(q.dash || []); g.lineWidth = q.width || 2; g.strokeStyle = color; g.fillStyle = color; g.globalAlpha = q.alpha || 1;
    if (type === 'hist' || type === 'outline') {
      const e = q.edges, n = q.y.length, w = n > 1 ? Math.abs(X(e[Math.min(n, Math.floor(n / 2) + 1)]) - X(e[Math.floor(n / 2)])) : 10;
      if (type === 'hist' && w <= 3.2) {
        // bars narrower than a few pixels are one filled outline: separate rectangles that end between two device pixels
        // overlap or leave a seam there, which shows as darker and lighter stripes that are not in the counts
        g.beginPath(); let open = false;
        for (let i = 0; i < n; i++) { const v = q.y[i], x0 = X(e[i]), x1 = X(e[i + 1]); if (!(v > 0)) { if (open) { g.lineTo(x0, by); g.closePath(); open = false; } continue; } const yy = Math.min(Y(v), by - 0.6); if (!open) { g.moveTo(x0, by); open = true; } g.lineTo(x0, yy); g.lineTo(x1, yy); if (i === n - 1) { g.lineTo(x1, by); g.closePath(); open = false; } }
        g.fill();
      } else if (type === 'hist') { for (let i = 0; i < n; i++) { const v = q.y[i]; if (!(v > 0)) continue; const x0 = X(e[i]), x1 = X(e[i + 1]), yy = Y(v); g.fillRect(x0, yy, Math.max(0.6, x1 - x0 - 1), by - yy); } }
      else { g.lineWidth = q.width || 1.5; g.beginPath(); let pen = false; for (let i = 0; i < n; i++) { const v = q.y[i] || 0, x0 = X(e[i]), x1 = X(e[i + 1]), yy = Y(v); if (!pen) { g.moveTo(x0, yy); pen = true; } else g.lineTo(x0, yy); g.lineTo(x1, yy); } g.stroke(); }
    } else if (type === 'vbar') {
      // stacked bars at x positions: at most 24 px thick, 1 px gap between segments
      const n = q.x.length, step = n > 1 ? Math.abs(X(q.x[1]) - X(q.x[0])) : 24, bw = Math.max(1.5, Math.min(24, step - (step > 6 ? 2 : 0.6)));
      for (let i = 0; i < n; i++) { let base = 0; const xx = X(q.x[i]) - bw / 2; for (const p of q.parts) { const v = p.y[i] || 0; if (v <= 0) continue; const ya = Y(base), yb = Y(base + v); g.fillStyle = p.color; const hgt = ya - yb - (base > 0 && ya - yb > 2 ? 1 : 0); if (hgt > 0.5) g.fillRect(xx, yb, bw, hgt); base += v; } }
    } else {
      const xs = q.x, ys = q.y, n = ys.length;
      // the band between two quartiles: a fill in the colour of the line, in stretches where both are known
      if (q.band) { const lo = q.band.lo, hi = q.band.hi; g.globalAlpha = q.band.alpha || 0.2; let i = 0;
        while (i < n) { while (i < n && !(fin(lo[i]) && fin(hi[i]))) i++; let j = i; while (j < n && fin(lo[j]) && fin(hi[j])) j++; if (j - i === 1) { const xx = X(xs[i]); g.fillRect(xx - 1.5, Y(hi[i]), 3, Math.max(1, Y(lo[i]) - Y(hi[i]))); } else if (j > i) { g.beginPath(); g.moveTo(X(xs[i]), Y(hi[i])); for (let k = i + 1; k < j; k++) g.lineTo(X(xs[k]), Y(hi[k])); for (let k = j - 1; k >= i; k--) g.lineTo(X(xs[k]), Y(lo[k])); g.closePath(); g.fill(); } i = j; }
        g.globalAlpha = q.alpha || 1; }
      const path = () => { g.beginPath(); let pen = false, px0 = 0, py0 = 0; for (let i = 0; i < n; i++) { const v = ys[i]; if (!fin(v)) { pen = false; continue; } const xx = X(xs[i]), yy = Y(v); if (!pen) { g.moveTo(xx, yy); pen = true; } else { if (type === 'step') g.lineTo(xx, py0); g.lineTo(xx, yy); } px0 = xx; py0 = yy; } };
      if (type === 'area') { g.globalAlpha = 0.14; g.beginPath(); let first = -1, last = -1; for (let i = 0; i < n; i++) if (fin(ys[i])) { if (first < 0) first = i; last = i; } if (first >= 0) { g.moveTo(X(xs[first]), by); for (let i = first; i <= last; i++) if (fin(ys[i])) g.lineTo(X(xs[i]), Y(ys[i])); g.lineTo(X(xs[last]), by); g.closePath(); g.fill(); } g.globalAlpha = q.alpha || 1; }
      if (cased) { g.setLineDash([]); g.strokeStyle = tk('--fc-chart-surface'); g.lineWidth = (q.width || 2) + 3; g.lineCap = 'round'; path(); g.stroke(); g.setLineDash(q.dash || []); g.strokeStyle = color; g.lineWidth = q.width || 2; }
      g.lineCap = 'round'; path(); g.stroke(); g.lineCap = 'butt';
      if (q.dots) for (let i = 0; i < n; i++) if (fin(ys[i])) { g.beginPath(); g.arc(X(xs[i]), Y(ys[i]), 2.5, 0, TAU); g.fill(); }
    }
    g.setLineDash([]); g.globalAlpha = 1;
  }

  // =========================== ranked bars (also read fate) ===========================
  paintHBars(g) {
    const s = this.spec, rows = s.rows, W = this.W, rowH = s.rowH || 26, th = Math.min(s.barH || 14, 24);
    g.font = font(); let lw = 0; for (const r of rows) lw = Math.max(lw, g.measureText(r.label).width);
    const most = s.labelW || Math.round(W * (printing.peek() ? 0.5 : 0.42)), useShort = lw + 10 > most && rows.some(r => r.short), text = r => useShort && r.short ? r.short : r.label;
    if (useShort) { lw = 0; for (const r of rows) lw = Math.max(lw, g.measureText(text(r)).width); }
    lw = Math.min(Math.ceil(lw) + 10, most);
    g.font = font('600'); let vw = 0; for (const r of rows) vw = Math.max(vw, g.measureText(r.text || '').width + (r.note ? g.measureText('  ' + r.note).width : 0));
    const x0 = lw, pw = Math.max(20, W - lw - Math.min(vw + 12, W * 0.45)), max = s.max || Math.max(1e-9, ...rows.map(r => r.value || 0));
    const f = s.log ? (v => v > 0 ? Math.log10(v) / Math.log10(max) : 0) : (v => v / max);
    this.L = { rows: rows.map((r, i) => ({ y: 4 + i * rowH, h: rowH })), x0, pw };
    g.textBaseline = 'middle';
    rows.forEach((r, i) => {
      const cy = 4 + i * rowH + rowH / 2;
      g.font = font(); g.textAlign = 'right'; g.fillStyle = tk('--fc-chart-label'); g.fillText(this.fit(g, text(r), lw - 8), x0 - 8, cy);
      if (r.unknown) { // a stage with no source: a dashed outline labelled "not known", never zero
        const uh = Math.min(rowH - 6, 20); g.strokeStyle = tk('--fc-fate-unknown-line'); g.lineWidth = 1; g.setLineDash([3, 3]); g.strokeRect(x0 + 0.5, Math.round(cy - uh / 2) + 0.5, pw - 1, uh - 1); g.setLineDash([]);
        g.textAlign = 'left'; g.fillStyle = tk('--fc-chart-tick'); g.fillText('not known', x0 + 8, cy);
        return;
      }
      const w = Math.max(r.value > 0 ? 1.5 : 0, f(r.value || 0) * pw);
      g.fillStyle = r.color || tk('--fc-this-run'); this.rrect(g, x0, cy - th / 2, w, th, 2); g.fill();
      g.textAlign = 'left'; g.font = font('600'); g.fillStyle = tk('--fc-text'); g.fillText(r.text || '', x0 + w + 6, cy);
      if (r.note) { const tw = g.measureText(r.text || '').width; g.font = font(); g.fillStyle = tk('--fc-chart-tick'); g.fillText(r.note, x0 + w + 6 + tw + 6, cy); }
    });
    g.strokeStyle = tk('--fc-chart-axis'); g.lineWidth = 1; g.beginPath(); g.moveTo(crisp(x0) - 1, 2); g.lineTo(crisp(x0) - 1, 4 + rows.length * rowH); g.stroke();
  }
  rrect(g, x, y, w, h, r) { r = Math.min(r, w / 2, h / 2); g.beginPath(); g.moveTo(x, y); g.lineTo(x + w - r, y); g.arcTo(x + w, y, x + w, y + r, r); g.lineTo(x + w, y + h - r); g.arcTo(x + w, y + h, x + w - r, y + h, r); g.lineTo(x, y + h); g.closePath(); }
  fit(g, text, w) { if (g.measureText(text).width <= w) return text; let a = Math.ceil(text.length * 0.55), b = Math.floor(text.length * 0.35); while (a + b > 4 && g.measureText(text.slice(0, a) + '…' + text.slice(-b)).width > w) { if (a > b) a--; else b--; } return text.slice(0, a) + '…' + text.slice(-b); }

  // =========================== one stacked bar per row ===========================
  paintHStack(g) {
    const s = this.spec, rows = s.rows, W = this.W, rowH = s.rowH || 24, th = Math.min(s.barH || 16, rowH - 6);
    // each label is measured in the weight it is drawn in: a chosen sample is bold and needs a little more room
    let lw = 0; for (const r of rows) { g.font = font(r.strong ? '600' : ''); lw = Math.max(lw, g.measureText(r.label).width); }
    // when the full names need more than their share of the width, every row is written with its short name (the words
    // the rows share left out), so that two rows never end up with the same cut text
    const most = s.labelW || Math.round(W * (printing.peek() ? 0.5 : 0.44)), useShort = lw + 10 > most && rows.some(r => r.short), text = r => useShort && r.short ? r.short : r.label;
    if (useShort) { lw = 0; for (const r of rows) { g.font = font(r.strong ? '600' : ''); lw = Math.max(lw, g.measureText(text(r)).width); } }
    lw = lw ? Math.min(Math.ceil(lw) + 10, most) : 2;
    const x0 = lw, pw = Math.max(20, W - lw - 12), max = s.norm ? 1 : Math.max(1e-9, ...rows.map(r => r.parts.reduce((a, p) => a + (p.value || 0), 0)));
    const top = 4, by = top + rows.length * rowH;
    this.L = { x0, pw, top, rowH, hits: [] };
    g.textBaseline = 'middle';
    rows.forEach((r, i) => {
      const cy = top + i * rowH + rowH / 2, total = r.parts.reduce((a, p) => a + (p.value || 0), 0) || 1;
      if (r.label) { g.font = font(r.strong ? '600' : ''); g.textAlign = 'right'; g.fillStyle = r.strong ? tk('--fc-text') : tk('--fc-chart-label'); g.fillText(this.fit(g, text(r), lw - 8), x0 - 8, cy); }
      let x = x0;
      r.parts.forEach((p, k) => {
        const w = (p.value || 0) / (s.norm ? total : max) * pw; if (w < 0.3) return;
        if (w >= 1.2) { g.fillStyle = p.color; g.fillRect(x, cy - th / 2, Math.max(0.5, w - (w > 3 ? 1 : 0)), th); }
        if (s.segLabels && p.label && w > 22) { g.font = font(); let lab = p.label; if (g.measureText(lab).width + 8 > w && p.labelShort) lab = p.labelShort; if (g.measureText(lab).width + 8 <= w) { g.textAlign = 'center'; g.fillStyle = onData(p.color); g.fillText(lab, x + w / 2, cy + 0.5); } }
        this.L.hits.push({ id: i * 1000 + k, x0: x, x1: x + w, y0: cy - rowH / 2, y1: cy + rowH / 2 });
        x += w;
      });
    });
    // share axis
    g.font = font(); g.fillStyle = tk('--fc-chart-tick'); g.textBaseline = 'top'; g.strokeStyle = tk('--fc-chart-axis'); g.lineWidth = 1;
    g.beginPath(); g.moveTo(x0, crisp(by + 2)); g.lineTo(x0 + pw, crisp(by + 2)); g.stroke();
    const ticks = s.norm ? [0, 0.25, 0.5, 0.75, 1] : niceTicks(0, max, 4);
    ticks.forEach((v, k) => { const xx = x0 + v / max * pw; if (xx > x0 + pw + 1) return; g.textAlign = k === 0 ? 'left' : (xx > x0 + pw - 16 ? 'right' : 'center'); g.fillText(s.norm ? Math.round(v * 100) + (v === 1 ? ' %' : '') : (s.fmt || tickLabel)(v), xx, by + 7); });
  }

  // =========================== parts of a whole: a ring with its labels beside it ===========================
  // Every part is named beside the ring, with its count and share (part.text), and joined to its arc by a line: nothing is
  // told by colour alone. A part too thin to see is drawn a little under half a percent wide; its label holds its true share.
  paintRing(g) {
    const s = this.spec, parts = s.parts, W = this.W, H = this.H, surface = tk('--fc-chart-surface');
    const total = parts.reduce((a, p) => a + (p.value > 0 ? p.value : 0), 0) || 1, rowH = 30;
    let lw = 0; for (const p of parts) { g.font = font('600'); lw = Math.max(lw, g.measureText(p.label).width); g.font = font(); lw = Math.max(lw, g.measureText(p.text || '').width); }
    const R = Math.max(26, Math.min(s.maxR || 104, (H - 12) / 2, (W - 2 * (Math.min(lw, W * 0.3) + 30)) / 2)), r0 = R * (s.inner || 0.6), cx = W / 2, cy = H / 2, room = cx - R - 26;
    const minA = TAU * 0.0045, raw = parts.map(p => (p.value > 0 ? p.value : 0) / total * TAU);
    let add = 0, big = 0; for (const a of raw) { if (a > 0 && a < minA) add += minA - a; else big += a; }
    const k = big > 0 ? Math.max(0, (big - add) / big) : 1, geo = []; let a0 = -Math.PI / 2;
    raw.forEach(a => { const w = a > 0 && a < minA ? minA : a * k; geo.push({ a0, a1: a0 + w, mid: a0 + w / 2 }); a0 += w; });
    geo.forEach((q, i) => { if (!(parts[i].value > 0)) return; g.beginPath(); g.arc(cx, cy, R, q.a0, q.a1); g.arc(cx, cy, r0, q.a1, q.a0, true); g.closePath(); g.fillStyle = parts[i].color; g.fill(); g.strokeStyle = surface; g.lineWidth = 1; g.stroke(); });
    if (s.center) { g.textAlign = 'center'; g.textBaseline = 'middle'; g.fillStyle = tk('--fc-text'); g.font = '600 ' + Math.round(Math.max(13, Math.min(22, r0 * 0.36))) + 'px ' + (tk('--fc-font-sans') || 'sans-serif'); g.fillText(s.center.value, cx, cy - (s.center.label ? 7 : 0)); if (s.center.label) { g.font = font(); g.fillStyle = tk('--fc-chart-tick'); g.fillText(s.center.label, cx, cy + 11); } }
    // labels: on the side of the ring their part is on, in the order of the parts, never closer than a row
    const labels = parts.map((p, i) => { const q = geo[i], c = Math.cos(q.mid), sn = Math.sin(q.mid); return { i, side: c >= 0 ? 1 : -1, ax: cx + c * (R + 2), ay: cy + sn * (R + 2), y: cy + sn * (R + 10) }; });
    for (const side of [1, -1]) { const col = labels.filter(l => l.side === side).sort((a, b) => a.y - b.y); if (!col.length) continue;
      const lo = rowH / 2 + 1, hi = H - rowH / 2 - 1; col.forEach((l, j) => { l.y = Math.max(l.y, j ? col[j - 1].y + rowH : lo); });
      for (let j = col.length - 1; j >= 0; j--) col[j].y = Math.min(col[j].y, j < col.length - 1 ? col[j + 1].y - rowH : hi);
      if (col[0].y < lo) col.forEach((l, j) => { l.y = lo + j * Math.min(rowH, (hi - lo) / Math.max(1, col.length - 1)); }); }
    g.textBaseline = 'middle'; g.lineWidth = 1;
    for (const l of labels) { const p = parts[l.i], ex = cx + l.side * (R + 12), tx = cx + l.side * (R + 22);
      g.strokeStyle = tk('--fc-chart-axis'); g.beginPath(); g.moveTo(l.ax, l.ay); g.lineTo(ex, crisp(l.y) - 0.5); g.lineTo(tx - l.side * 4, crisp(l.y) - 0.5); g.stroke();
      g.textAlign = l.side > 0 ? 'left' : 'right'; g.font = font('600'); g.fillStyle = tk('--fc-text'); g.fillText(this.fit(g, p.label, Math.max(30, room)), tx, l.y - 7);
      g.font = font(); g.fillStyle = tk('--fc-chart-label'); g.fillText(this.fit(g, p.text || '', Math.max(30, room)), tx, l.y + 7); l.x0 = l.side > 0 ? tx - 12 : 0; l.x1 = l.side > 0 ? W : tx + 12; }
    this.L = { ring: { cx, cy, R, r0, geo }, labels, rowH };
  }

  // =========================== flow cell map ===========================
  paintMap(g) {
    const s = this.spec, { cols, rows, pos, values, max } = s, W = this.W, H = this.H;
    const cell = Math.max(2, Math.min(s.maxCell || 16, Math.floor((W - 2) / cols * 100) / 100, Math.floor((H - 2) / rows * 100) / 100));
    const gap = cell >= 6 ? 1 : 0.6, ox = Math.floor((W - cols * cell) / 2), oy = 1;
    this.L = { cell, ox, oy, cols, rows };
    const none = tk('--fc-chart-none'), noneLine = tk('--fc-chart-none-line'), hatch = tk('--fc-chart-none-hatch');
    const n = values.length;
    for (let ch = 0; ch < n; ch++) {
      const [r, c] = pos(ch), x = ox + c * cell, y = oy + r * cell, v = values[ch], sz = cell - gap;
      if (s.nodata ? s.nodata(v) : !(v > 0)) {
        g.fillStyle = none; g.fillRect(x, y, sz, sz); g.strokeStyle = noneLine; g.lineWidth = cell >= 6 ? 0.8 : 0.6; g.strokeRect(x + 0.4, y + 0.4, sz - 0.8, sz - 0.8);
        if (cell >= 6) { g.strokeStyle = hatch; g.lineWidth = 0.8; g.beginPath(); g.moveTo(x + 1, y + sz - 1); g.lineTo(x + sz - 1, y + 1); g.stroke(); }
        continue;
      }
      g.fillStyle = s.colorOf ? s.colorOf(v) : seqColor(Math.min(1, v / (max || 1)));
      g.fillRect(x, y, sz, sz);
    }
  }

  // =========================== interaction ===========================
  // The last pointer position of a frame is the one that counts: [clientX, clientY], or false when the pointer left.
  onPointer(p) {
    this.ptr = p; if (this.raf) return;
    this.raf = requestAnimationFrame(() => { this.raf = 0; const q = this.ptr; this.ptr = null; if (q === false) { if (!this.drag) this.setCursor(-1, true); } else if (q) this.onMove(q[0], q[1]); });
  }
  local(e) { const r = this.ov.getBoundingClientRect(); return [e.clientX - r.left, e.clientY - r.top]; }
  onMove(cx, cy) {
    const s = this.spec, L = this.L; if (!L) return;
    const [x, y] = this.local({ clientX: cx, clientY: cy }); this.byPointer = true;
    if (this.drag) { this.drag.x1 = clamp(x, L.px, L.px + L.pw); this.drawCursor(); return; }
    if (s.form === 'hbars') { const n = s.rows.length; return this.setCursor(n ? clamp(Math.floor((y - 4) / (s.rowH || 26)), 0, n - 1) : -1); }
    if (s.form === 'hstack') {
      // the row under the pointer, and in it the segment under the pointer or the nearest one
      const row = clamp(Math.floor((y - L.top) / L.rowH), 0, s.rows.length - 1); let hit = -1, far = Infinity;
      for (const h of L.hits) { if (Math.floor(h.id / 1000) !== row) continue; const d = x < h.x0 ? h.x0 - x : x >= h.x1 ? x - h.x1 + 0.01 : 0; if (d < far) { far = d; hit = h.id; } }
      return this.setCursor(hit);
    }
    if (s.form === 'map') { const ch = s.at(clamp(Math.floor((y - L.oy) / L.cell), 0, L.rows - 1), clamp(Math.floor((x - L.ox) / L.cell), 0, L.cols - 1)); if (ch != null) this.setCursor(ch); return; }
    if (s.form === 'ring') {
      const r = L.ring, dx = x - r.cx, dy = y - r.cy; let hit = -1;
      if (Math.hypot(dx, dy) <= r.R + 8) { let a = Math.atan2(dy, dx); if (a < -Math.PI / 2) a += TAU; hit = r.geo.findIndex((q, i) => s.parts[i].value > 0 && a >= q.a0 && a < q.a1); if (hit < 0) hit = r.geo.reduce((b, q, i) => s.parts[i].value > 0 && Math.abs(q.mid - a) < Math.abs(r.geo[b].mid - a) ? i : b, 0); }
      else { let far = Infinity; for (const l of L.labels) if ((dx >= 0) === (l.side > 0) && Math.abs(y - l.y) < far) { far = Math.abs(y - l.y); hit = l.i; } }
      return this.setCursor(hit);
    }
    if (!s.cursor || !s.cursor.x.length) return;
    const xv = L.invX(clamp(x, L.px, L.px + L.pw)); this.setCursor(nearest(s.cursor.x, xv), true, xv);
  }
  onDown(e) { const s = this.spec, L = this.L; if (!s.onBrush || !L || L.px == null) return; const [x] = this.local(e); if (x < L.px || x > L.px + L.pw) return; this.drag = { x0: x, x1: x }; try { this.ov.setPointerCapture(e.pointerId); } catch (err) { /* synthetic events */ } }
  onUp(e) { const d = this.drag; if (!d) return; this.drag = null; this.ptr = null; try { this.ov.releasePointerCapture(e.pointerId); } catch (err) { /* nothing held */ } const L = this.L; d.x1 = clamp(this.local(e)[0], L.px, L.px + L.pw);      // where the button was let go, not where the last frame saw it
    if (Math.abs(d.x1 - d.x0) > 5) { const a = L.invX(Math.min(d.x0, d.x1)), b = L.invX(Math.max(d.x0, d.x1)); this.spec.onBrush([a, b]); } this.drawCursor(); }
  // i: the bin that is read (-1: none). xv: where on the x axis the pointer is; the hairline stands there, so that
  // charts that share a cursor draw one straight line. Without xv (keyboard) the hairline stands on the bin.
  setCursor(i, broadcast, xv) {
    const x = i < 0 || xv == null ? null : xv, moved = i !== this.cur;
    if (!moved && x === this.curX && !broadcast) return;
    this.cur = i; this.curX = x; this.drawCursor();
    if (moved && this.spec.onHover) this.spec.onHover(i, i >= 0 && !!this.byPointer);
    if (broadcast && this.syncKey && bus.has(this.syncKey)) { const at = i < 0 ? null : (xv != null ? xv : this.spec.cursor.x[i]); for (const c of bus.get(this.syncKey)) if (c !== this) c.follow(at); }
  }
  follow(xv) { const s = this.spec; if (!s.cursor || !this.L) return; const i = xv == null || !s.cursor.x.length ? -1 : nearest(s.cursor.x, xv); if (i >= 0 && (xv < this.L.d0 - 1e-9 || xv > this.L.d1 + 1e-9)) return this.setCursor(-1); this.setCursor(i, false, xv); }
  // where the cursor mark stands, in the plot's own pixels: [x, y, half its size]. Only a map says it (the shell puts a note beside the cell).
  anchor() { const s = this.spec, L = this.L, i = this.cur; if (i < 0 || !L || s.form !== 'map') return null; const [r, c] = s.pos(i); return [L.ox + (c + 0.5) * L.cell, L.oy + (r + 0.5) * L.cell, L.cell / 2]; }
  move(delta) { this.byPointer = false; if (this.spec.form === 'hstack') { const hs = this.L.hits; if (!hs.length) return; let k = hs.findIndex(q => q.id === this.cur); k = k < 0 ? (delta > 0 ? 0 : hs.length - 1) : Math.max(0, Math.min(hs.length - 1, delta === -Infinity ? 0 : delta === Infinity ? hs.length - 1 : k + delta)); return this.setCursor(hs[k].id, true); }
    const s = this.spec, n = s.form === 'hbars' ? s.rows.length : s.form === 'ring' ? s.parts.length : s.form === 'map' ? s.values.length : s.cursor ? s.cursor.x.length : 0; if (!n) return;
    let lo = 0, hi = n - 1;                                  // on an axis, the keyboard stays inside what is drawn
    if (s.cursor && this.L && this.L.d0 != null) { while (lo < hi && s.cursor.x[lo] < this.L.d0) lo++; while (hi > lo && s.cursor.x[hi] > this.L.d1) hi--; }
    let i = this.cur < 0 ? (delta > 0 ? lo : hi) : this.cur + delta; if (delta === -Infinity) i = lo; if (delta === Infinity) i = hi; i = Math.max(lo, Math.min(hi, i)); this.setCursor(i, true); }
  // Where a line stands at x: between two of its points on the line that joins them (a step holds its last value),
  // in pixels. null where the line has a gap or does not reach.
  lineY(q, x) {
    const xs = q.x, ys = q.y, L = this.L, n = Math.min(xs.length, ys.length); if (!n || x < xs[0] || x > xs[n - 1]) return null;
    const b = Math.min(n - 1, lowerBound(xs, x)), a = xs[b] === x ? b : b - 1; if (a < 0) return null;
    if (a === b) return fin(ys[a]) ? L.Y(ys[a]) : null;
    if (!fin(ys[a]) || !fin(ys[b])) return null;
    if (q.type === 'step') return L.Y(ys[a]);
    const xa = L.X(xs[a]), xb = L.X(xs[b]), ya = L.Y(ys[a]), yb = L.Y(ys[b]);
    return xb > xa ? ya + (yb - ya) * (L.X(x) - xa) / (xb - xa) : ya;
  }
  drawCursor() {
    const g = this.og, s = this.spec, L = this.L, i = this.cur;
    g.setTransform(this.dpr, 0, 0, this.dpr, 0, 0);
    // only what the last cursor covered is cleared, and what this one covers is remembered
    if (this.dirty) { g.clearRect(this.dirty[0], this.dirty[1], this.dirty[2], this.dirty[3]); this.dirty = null; }
    const mark = (x, y, w, h) => { const d = this.dirty, x0 = Math.floor(x) - 3, y0 = Math.floor(y) - 3, x1 = Math.ceil(x + w) + 3, y1 = Math.ceil(y + h) + 3; this.dirty = d ? [Math.min(d[0], x0), Math.min(d[1], y0), Math.max(d[0] + d[2], x1) - Math.min(d[0], x0), Math.max(d[1] + d[3], y1) - Math.min(d[1], y0)] : [x0, y0, x1 - x0, y1 - y0]; };
    if (this.drag && L) { const a = Math.min(this.drag.x0, this.drag.x1), b = Math.max(this.drag.x0, this.drag.x1); g.fillStyle = alpha(tk('--fc-text'), 0.10); g.fillRect(a, L.py, b - a, L.ph); g.strokeStyle = tk('--fc-evidence-line'); g.lineWidth = 1; g.beginPath(); g.moveTo(crisp(a), L.py); g.lineTo(crisp(a), L.py + L.ph); g.moveTo(crisp(b), L.py); g.lineTo(crisp(b), L.py + L.ph); g.stroke(); mark(a, L.py, b - a, L.ph); return; }
    if (i < 0 || !L) return;
    if (s.form === 'hbars') { const r = L.rows[i]; if (!r) return; g.strokeStyle = tk('--fc-text'); g.lineWidth = 1; g.strokeRect(0.5, r.y + 0.5, this.W - 1, r.h - 1); mark(0, r.y, this.W, r.h); return; }
    if (s.form === 'hstack') { const h = L.hits.find(q => q.id === i); if (!h) return; g.strokeStyle = tk('--fc-text'); g.lineWidth = 1; g.strokeRect(crisp(h.x0), h.y0 + 3.5, Math.max(2, h.x1 - h.x0 - 1), h.y1 - h.y0 - 7); mark(h.x0, h.y0, Math.max(3, h.x1 - h.x0), h.y1 - h.y0); return; }
    if (s.form === 'map') { const [r, c] = s.pos(i); g.strokeStyle = tk('--fc-text'); g.lineWidth = 1.5; g.strokeRect(L.ox + c * L.cell - 1, L.oy + r * L.cell - 1, L.cell + 1, L.cell + 1); mark(L.ox + c * L.cell - 2, L.oy + r * L.cell - 2, L.cell + 3, L.cell + 3); return; }
    if (s.form === 'ring') {       // the part that is read: its arc outlined, and a rule under its label
      const r = L.ring, q = r.geo[i], l = L.labels.find(o => o.i === i); if (!q || !l) return;
      g.strokeStyle = tk('--fc-text'); g.lineWidth = 1.5; g.lineJoin = 'round';
      if (s.parts[i].value > 0) { g.beginPath(); g.arc(r.cx, r.cy, r.R + 1, q.a0, q.a1); g.arc(r.cx, r.cy, r.r0 - 1, q.a1, q.a0, true); g.closePath(); g.stroke(); }
      const tx = r.cx + l.side * (r.R + 22), far = l.side > 0 ? Math.min(this.W - 2, tx + 120) : Math.max(2, tx - 120); g.lineWidth = 1; g.beginPath(); g.moveTo(tx, crisp(l.y + 15)); g.lineTo(far, crisp(l.y + 15)); g.stroke();
      mark(r.cx - r.R - 3, r.cy - r.R - 3, 2 * r.R + 6, 2 * r.R + 6); mark(Math.min(tx, far) - 1, l.y + 13, Math.abs(far - tx) + 2, 4); return;
    }
    if (!s.cursor) return;
    const xv = this.curX != null ? this.curX : s.cursor.x[i], xx = crisp(L.X(xv)); if (xx < L.px - 1 || xx > L.px + L.pw + 1) return;
    const by = L.py + L.ph;
    g.strokeStyle = tk('--fc-chart-cursor'); g.lineWidth = 1; g.beginPath(); g.moveTo(xx, L.py); g.lineTo(xx, by); g.stroke(); mark(xx - 1, L.py, 2, L.ph);
    if (s.stack) return;
    for (const q of s.series || []) {
      if (q.type === 'vbar') {       // bars stand apart: the one that is read is outlined
        if (!q.x || i >= q.x.length) continue; let tot = 0; for (const p of q.parts) tot += p.y[i] > 0 ? p.y[i] : 0; if (!(tot > 0)) continue;
        const n = q.x.length, step = n > 1 ? Math.abs(L.X(q.x[1]) - L.X(q.x[0])) : 24, bw = Math.max(1.5, Math.min(24, step - (step > 6 ? 2 : 0.6))), x0 = L.X(q.x[i]) - bw / 2, y0 = Math.max(L.py, L.Y(tot));
        if (x0 + bw < L.px || x0 > L.px + L.pw) continue;
        g.strokeStyle = tk('--fc-text'); g.lineWidth = 1; g.strokeRect(Math.round(x0) - 0.5, Math.round(y0) - 0.5, Math.round(bw) + 1, Math.round(by - y0) + 1); mark(x0 - 1, y0 - 1, bw + 2, by - y0 + 2); continue;
      }
      if (q.type === 'hist' || q.type === 'outline' || !q.x) continue;
      // a dot on each line where the hairline crosses it; on the bin itself where the line has a gap beside it
      let yy = this.lineY(q, xv), dx = xx;
      if (yy == null && q.x.length === s.cursor.x.length && fin(q.y[i])) { yy = L.Y(q.y[i]); dx = L.X(q.x[i]); }
      if (yy == null || dx < L.px - 1 || dx > L.px + L.pw + 1) continue;
      g.fillStyle = tk('--fc-chart-surface'); g.beginPath(); g.arc(dx, yy, 5, 0, TAU); g.fill(); g.fillStyle = q.color || tk('--fc-this-run'); g.beginPath(); g.arc(dx, yy, 3.5, 0, TAU); g.fill(); mark(dx - 5, yy - 5, 10, 10);
    }
  }
}

// A genome track as an 'xy' spec: one bar per window along the reference, a hairline at every chromosome boundary, the
// chromosome names under the axis (the ones that have room), and the mean of each chromosome as a line.
//   contigs: [{ name, first (index of its first window), windows }]   values: one number per window (depth)
//   means: one number per contig (optional)   top: the top of the y axis (optional: bars above it are cut off)
export function genomeTrack({ contigs, values, means, top, title, fmt, log, color, readout, name }) {
  const n = values.length, edges = new Float64Array(n + 1), mid = new Float64Array(n);
  for (let i = 0; i <= n; i++) edges[i] = i;
  for (let i = 0; i < n; i++) mid[i] = i + 0.5;
  const y = { title, fmt, grid: true }; if (log) y.log = true; else if (top != null) y.domain = [0, top];
  return { form: 'xy', mr: 6,
    x: { domain: [0, n], tickMarks: false, bounds: contigs.slice(1).map(c => c.first), ticks: contigs.map(c => ({ v: c.first + c.windows / 2, label: c.label != null ? c.label : c.name, major: true })) },
    y, series: [{ name: name || 'Depth', type: 'hist', edges, y: values, color }],
    hsegs: means ? contigs.map((c, k) => ({ x0: c.first, x1: c.first + c.windows, y: means[k] })) : [],
    cursor: { x: mid }, readout };
}
