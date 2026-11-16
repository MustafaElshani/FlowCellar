// Paper. A page of the viewer is saved as a PDF through the browser's own print dialog, with a layout made for
// paper: "Save as PDF" on a page and the browser's Print command (Ctrl+P) go the same way and give the same pages.
//
// What happens between the browser's `beforeprint` and `afterprint` events, all of it before the handler returns
// (a browser takes its picture of the page the moment `beforeprint` is done, so nothing here may wait for a frame):
//   1. <html> gets the class "fc-print" and the theme "print". The style sheet lays the page out at the width of
//      the paper (--fc-paper-w, styles/app.css) with no navigation and no controls, on screen as on paper, so a
//      chart is measured at the width it is printed at.
//   2. Every screen is drawn again, with the colours of paper: the specs of the charts hold colours, and those were
//      read from the theme of the screen (main.jsx reads `printing` for this: the colours are read again, and the
//      application draws itself again, whenever it changes). Closed groups of further charts are opened when that
//      was asked for.
//   3. Every chart paints itself again at its new size, PRINT_SCALE device pixels to a CSS pixel.
//   4. A heading is tied to the first block under it, so that it never stands alone at the foot of a sheet.
//   5. A table that is wider than the sheet is scaled down until it fits; links lose their addresses, so that no
//      address of this computer travels with the PDF.
// `afterprint` undoes all five: the same layout, theme, open groups and scroll position as before.
import { options } from 'preact';
import { signal, batch } from '@preact/signals';

export const printing = signal(false);        // true while the page is laid out for paper
export const printAll = signal(false);        // closed groups of further charts are opened for this print
// Device pixels per CSS pixel of a chart on paper. Firefox writes a chart into the PDF as a picture: at 3 it is still
// sharp when the sheet is read at 300 %. Chromium writes a chart as a drawing when it was painted during `beforeprint`,
// which is why the painting happens there and not before.
export const PRINT_SCALE = 3;
export const pixelRatio = () => printing.peek() ? PRINT_SCALE : Math.min(2, window.devicePixelRatio || 1);

// ---- charts: each one says how it paints itself again at the size and scale of now ----
const painters = new Set();
export function onPaper(paint) { painters.add(paint); return () => { painters.delete(paint); }; }
function repaint() { for (const paint of [...painters]) { try { paint(); } catch (e) { /* a chart that has left the page */ } } }

// ---- the name a browser offers for the PDF file is the title of the page ----
let fileTitle = '';
export function setPrintTitle(text) { fileTitle = String(text || '').replace(/[^\w .()+-]+/g, ' ').replace(/\s+/g, ' ').trim().slice(0, 120); }

// Draw the page again now: the renders a change asks for, and the effects of those renders (a chart is made and
// painted in an effect), before this function returns. The framework would do both later, after the next frame;
// its two scheduling hooks are held for the length of the call, as its own test helper does.
function renderNow(change) {
  const debounce = options.debounceRendering, frame = options.requestAnimationFrame;
  let drain = null, effects = [];
  options.debounceRendering = cb => { drain = cb; };
  options.requestAnimationFrame = cb => { effects.push(cb); };
  const render = () => { const d = drain; drain = null; if (d) d(); };
  try {
    change(); render();
    for (let turn = 0; effects.length && turn < 25; turn++) { const run = effects; effects = []; for (const f of run) f(); render(); }
  } finally {
    options.debounceRendering = debounce; options.requestAnimationFrame = frame;
    // whatever is still waiting goes the usual way
    if (drain) setTimeout(drain, 0);
    for (const f of effects) setTimeout(f, 0);
  }
}

// ---- a heading stays with what it heads ----
// Browsers differ in what they honour of "do not break after this" (Firefox honours none of it), so the heading is
// made one block that may not be cut, as tall as itself and the first block under it: the style sheet ends the
// heading in an empty box of the height --fc-keep and takes that height back with a margin. The height is measured
// here, in the layout of the paper.
const HEADS = '.fc-section__head, .fc-panel__head, .fc-disclosure.is-open, .fc-allchecks__head';
const BLOCKS = '.fc-chart, .fc-tiles, .fc-track, .fc-empty, .fc-ancard, .fc-readcard, tbody > tr, .fc-runstrip';
let kept = [];
function keepHeads() {
  const found = [];
  for (const head of document.querySelectorAll(HEADS)) {
    const top = head.getBoundingClientRect();
    if (!top.height) continue;
    // the first block under the heading: among what follows it inside its section, panel or group
    const scope = head.closest('.fc-section, .fc-panel, .fc-allchecks, .fc-verdict') || head.parentElement;
    let first = null;
    for (const el of scope.querySelectorAll(BLOCKS)) { if (!(head.compareDocumentPosition(el) & Node.DOCUMENT_POSITION_FOLLOWING)) continue; const r = el.getBoundingClientRect(); if (r.height > 0) { first = r; break; } }
    // down to the foot of that block and a little further (under a table: its first row and most of the second), at most two thirds of a sheet
    const keep = first ? Math.max(0, Math.min(640, Math.round(first.bottom - top.bottom) + 30)) : 60;
    found.push([head, keep]);
  }
  for (const [head, keep] of found) head.style.setProperty('--fc-keep', keep + 'px');
  kept = found.map(f => f[0]);
}
function freeHeads() { for (const head of kept) head.style.removeProperty('--fc-keep'); kept = []; }

// what is not in the PDF because it was closed on the page, said in the line under the selection (shell/paper.jsx)
function sayLeftOut() {
  const line = document.querySelector('.fc-paper__left'); if (!line) return;
  const closed = [...document.querySelectorAll('.fc-disclosure:not(.is-open)')], charts = closed.filter(d => d.dataset.group === 'charts').length, tables = closed.length - charts;
  const n = (k, one, many) => k + ' ' + (k === 1 ? one : many), every = line.dataset.every === '1';
  const out = [charts ? n(charts, 'group', 'groups') + ' of further charts' : null, tables ? n(tables, 'group', 'groups') + ' of numbers and details' : null].filter(Boolean).join(' and ');
  line.textContent = (every ? 'With every chart: the groups of further charts that were closed on the page are opened here. ' : 'The page as it was shown. ')
    + (out ? 'Left out, because closed on the page: ' + out + '.' : every ? 'Nothing that the page holds is left out.' : 'No group was closed on the page: nothing is left out.');
}

// ---- a table that is wider than the sheet is made smaller until it fits ----
// A comparison of eight runs has ten columns, and a sheet is narrower than any screen such a table was made for. A
// table that does not fit would make the browser shrink every sheet (Chromium, Firefox) or cut the table off. The
// table alone is scaled down instead (the `zoom` property: its rows still break between sheets), to half its size at
// most; below that, the browser's own shrinking is the lesser evil.
let fitted = [];
function fitTables() {
  fitted = [];
  for (const t of document.querySelectorAll('#app table')) {
    const box = t.parentElement, room = box ? box.clientWidth : 0; if (!(room > 0)) continue;
    for (let turn = 0, z = 1; turn < 3; turn++) {
      const w = t.getBoundingClientRect().width; if (w <= room + 0.5) break;
      z = Math.max(0.5, Math.floor(z * room / w * 1000) / 1000 - 0.004); t.style.zoom = String(z); if (!fitted.includes(t)) fitted.push(t);
      if (z <= 0.5) break;
    }
  }
}
function freeTables() { for (const t of fitted) t.style.zoom = ''; fitted = []; }

// ---- no link leaves the page with the PDF ----
// A browser writes every link of the page into the PDF with its whole address (Firefox does), and the address of a page of
// the viewer is where it lies: a path on this computer, or the address of a server that only this computer reaches. On
// paper the links are plain text: their addresses are taken off for the print and put back afterwards.
let links = [];
function unlink() { links = []; for (const a of document.querySelectorAll('a[href]')) { links.push([a, a.getAttribute('href')]); a.removeAttribute('href'); } }
function relink() { for (const [a, href] of links) if (!a.hasAttribute('href')) a.setAttribute('href', href); links = []; }

let was = null;
export function enterPrint() {
  if (was) return;
  const root = document.documentElement;
  was = { theme: root.getAttribute('data-theme'), x: window.scrollX, y: window.scrollY, title: document.title };
  root.classList.add('fc-print'); root.setAttribute('data-theme', 'print');
  renderNow(() => { printing.value = true; });
  repaint();
  sayLeftOut();
  fitTables();
  keepHeads();
  unlink();
  if (fileTitle) document.title = fileTitle;
}
export function leavePrint() {
  if (!was) return;
  const root = document.documentElement, s = was; was = null;
  relink(); freeHeads(); freeTables();
  root.classList.remove('fc-print');
  if (s.theme == null) root.removeAttribute('data-theme'); else root.setAttribute('data-theme', s.theme);
  document.title = s.title;
  renderNow(() => batch(() => { printing.value = false; printAll.value = false; }));
  repaint();
  window.scrollTo(s.x, s.y);
}
if (typeof addEventListener === 'function') {
  addEventListener('beforeprint', enterPrint);
  addEventListener('afterprint', leavePrint);
  // A browser that never says `afterprint` would leave the page in the layout of the paper: the first key or pointer
  // press brings the screen back (while a print dialog is open the page receives neither).
  for (const type of ['pointerdown', 'keydown']) addEventListener(type, () => { if (was) leavePrint(); }, true);
}

// "Save as PDF": the layout for paper, then the print dialog of the browser. everyChart: closed groups of further
// charts are opened for the PDF, and closed again afterwards.
// The layout itself is made when the browser says `beforeprint`, inside window.print(): a chart painted then, and
// not before, is what a browser can write into the PDF as a drawing instead of a picture.
export function savePdf(everyChart) {
  printAll.value = !!everyChart;
  try { window.print(); } catch (e) { /* a page that may not print */ }
  // a browser that has no dialog (or printed without saying so) must not leave the groups open for the next print
  setTimeout(() => { if (!was) printAll.value = false; }, 0);
}
