// "Save as PDF": a run page, an analysis page, a comparison of runs and the history of a flow cell go to paper, in
// Chromium and in Firefox, on A4 and on Letter, and come back to the screen as they were.
//
//   node verify/print.mjs --fixture FILE.json --out DIR [--papers a4,letter] [--engines chromium,firefox] [--quick]
//   node verify/print.mjs --out DIR --page NAME=KIND=URL [--page ...]      any page of any hub or report (KIND: run,
//                                                                          analysis, compare, flowcell, project)
//
//   --fixture  the line of JSON that tests/viewer/viewer_print_fixture.py prints, in a file: one-run reports with
//              3 and 13 samples, the hub of the synthetic library, a hub and a report with spatial analyses
//   --quick    A4 only and the first page of each kind
//   --only     NAME[,NAME]: these pages only (the names are printed with every check)
//
// For every page, in every browser:
//   in the page   the control is there and opens the print dialog; in the layout for paper nothing is wider than
//                 the sheet, no navigation or control is left, the head and the selection line are there, the theme
//                 is the light one of paper whatever the screen shows, and every chart is painted again at the width
//                 it has on the sheet with at least twice as many pixels; afterwards every element of the page is
//                 where it was, with the same theme, title, open groups and scroll position
//   in the PDF    (Chromium through Playwright, Firefox through its own print pipeline, verify/firefox_pdf.mjs) the
//                 sheets have the size asked for, their number is within bounds, none is blank, the text is text
//                 (the verdict sentence, the head, the selection line) and holds no control and no link, no heading stands
//                 alone at the foot of a sheet, pictures have at least twice the pixels of their size, and the file
//                 stays under 5 MB
// Exit code 0 when every check passed and nothing was logged, 1 otherwise, 2 when no Chromium is installed.
// Writes OUT/print_report.json, the PDFs and one picture of each page in the layout for paper.
import fs from 'node:fs';
import path from 'node:path';
import zlib from 'node:zlib';
import { execFileSync } from 'node:child_process';
import { launch, fileUrl, settle, setTheme, hasBrowser, viewerDir } from './browser.mjs';
import { firefoxSession } from './firefox_pdf.mjs';

const args = process.argv.slice(2), opt = n => { const i = args.indexOf('--' + n); return i >= 0 ? args[i + 1] : null; }, flag = n => args.includes('--' + n);
const outDir = path.resolve(opt('out') || 'verify-out'), quick = flag('quick');
const papers = (opt('papers') || (quick ? 'a4' : 'a4,letter')).split(','), engines = (opt('engines') || 'chromium,firefox').split(',');
const PAPER = { a4: { format: 'A4', pts: [595, 842], mm: [210, 297] }, letter: { format: 'Letter', pts: [612, 792], mm: [215.9, 279.4] } };
const MARGIN_MM = 10, PAPER_W = 714, MAX_BYTES = 5 * 1024 * 1024;
if (engines.includes('chromium') && !hasBrowser('chromium')) { console.error('No Playwright Chromium found. Install one with:\n  PLAYWRIGHT_BROWSERS_PATH=' + path.join(viewerDir, 'pw-browsers') + ' npx playwright install chromium'); process.exit(2); }
fs.mkdirSync(outDir, { recursive: true });

// ---- the pages ----
const KIND = {
  run: { ready: '.fc-run .fc-track canvas', word: /Run report/, pages: [5, 60], shows: /Samples: |Time: /, button: true },
  analysis: { ready: '.fc-analysis .fc-runmain .fc-section', word: /Analysis report/, pages: [2, 14], shows: /Pipeline: |Scale/, button: false },
  compare: { ready: '.fc-compare .fc-plot canvas', word: /Comparison of \d+ runs/, pages: [2, 10], shows: /Count by: /, button: false },
  flowcell: { ready: '.fc-flowcell .fc-fcruns', word: /Flow cell history/, pages: [2, 10], shows: /Time axis: /, button: false },
  // the page of a project: on paper a report of its runs, drawn from the index alone, so it holds tables and no chart
  project: { ready: '.fc-projectpage .fc-hub', word: /Project report/, pages: [1, 6], shows: /Runs: all \d+ runs? of the project/, button: false, charts: false },
};
const targets = [];
if (opt('fixture')) {
  const fx = JSON.parse(fs.readFileSync(opt('fixture'), 'utf8').trim().split('\n').pop());
  const gun = b64 => JSON.parse(zlib.gunzipSync(Buffer.from(b64, 'base64')).toString('utf8'));
  const ix = gun(/^__FC_INDEX__\("([A-Za-z0-9+/=]*)"\);\n$/.exec(fs.readFileSync(path.join(fx.hub, 'index.js'), 'utf8'))[1]), C = ix.runs.columns;
  const col = (k, i) => { const c = C[k]; return c == null ? undefined : (!Array.isArray(c) && c.dict ? c.dict[c.codes[i]] : c[i]); };
  const rows = Array.from({ length: ix.runs.n }, (_, i) => ({ id: col('id', i), fc: col('fc', i), detail: col('detail_t1', i), use_of: col('use_of', i), read_data: col('read_data', i), started: col('started', i) }));
  const full = rows.filter(r => r.detail === 'full').sort((a, b) => Date.parse(b.started) - Date.parse(a.started)), hub = fileUrl(path.join(fx.hub, 'index.html'));
  const reused = (rows.find(r => r.use_of > 1) || {}).fc;
  for (const n of Object.keys(fx.runs)) targets.push({ name: 'run_' + n + '_samples', kind: 'run', url: fileUrl(fx.runs[n]), every: n === Object.keys(fx.runs)[0] });
  if (fx.spatial) { const sh = fileUrl(path.join(fx.spatial.hub, 'index.html'));
    targets.push({ name: 'analysis_spatial', kind: 'analysis', url: sh + '#/analysis/' + fx.spatial.analysis, map: true });
    targets.push({ name: 'analysis_report', kind: 'analysis', url: fileUrl(fx.spatial.report), map: true }); }
  if (full.length >= 3) targets.push({ name: 'compare_3_runs', kind: 'compare', url: hub + '#/compare?items=' + full.slice(0, 3).map(r => r.id).join(',') });
  // eight runs: the table of differences has ten columns and is scaled down to the width of the sheet
  if (full.length >= 8) targets.push({ name: 'compare_8_runs', kind: 'compare', url: hub + '#/compare?items=' + full.slice(0, 8).map(r => r.id).join(','), fits: true });
  if (reused) targets.push({ name: 'flow_cell_history', kind: 'flowcell', url: hub + '#/flowcell/' + reused });
  const project = (ix.projects || []).filter(q => q.id !== '_none').sort((a, b) => (b.runs || 0) - (a.runs || 0))[0];
  if (project) targets.push({ name: 'project_report', kind: 'project', url: hub + '#/project/' + encodeURIComponent(project.id) });
  // the screen in the dark theme, scrolled down, with a group open: paper is light all the same, and the screen comes back as it was
  const withData = full.find(r => r.read_data && r.read_data !== 'absent') || full[0];
  if (withData) targets.push({ name: 'run_from_a_dark_screen', kind: 'run', url: hub + '#/run/' + withData.id + '?more=fc', dark: true, scroll: 900 });
}
for (let i = 0; i < args.length; i++) if (args[i] === '--page') { const m = /^([^=]+)=([a-z]+)=(.+)$/.exec(args[i + 1] || ''); if (!m || !KIND[m[2]]) { console.error('--page NAME=KIND=URL, KIND one of ' + Object.keys(KIND).join(', ')); process.exit(2); } targets.push({ name: m[1], kind: m[2], url: m[3] }); }
if (!targets.length) { console.error('usage: node verify/print.mjs --fixture FILE.json --out DIR   (or --page NAME=KIND=URL)'); process.exit(2); }
const only = (opt('only') || '').split(',').filter(Boolean);
const chosen = (quick ? targets.filter((t, k) => targets.findIndex(q => q.kind === t.kind) === k) : targets).filter(t => !only.length || only.includes(t.name));

const report = { started: new Date().toISOString(), papers, engines, checks: [], pdfs: [], shots: [], errors: [], requests: [] };
const check = (name, ok, detail) => { report.checks.push({ name, ok: !!ok, detail: detail == null ? '' : String(detail).slice(0, 400) }); console.log((ok ? 'ok   ' : 'FAIL ') + name + (detail != null && detail !== '' ? '  [' + String(detail).slice(0, 220) + ']' : '')); };

// ---- what a PDF holds, read with poppler's tools ----
const NO_CONTROL = ['Save as PDF', 'PDF with every chart', 'Show as table', 'Hide table', 'How to read', 'Add to compare', 'In compare', 'Skip to the content', 'Back to summary', 'Open analysis', 'Tissue outline', 'Colour by', 'Whole run ·', 'Remove', 'Add a run from'];
function pdfFacts(file) {
  const run = (cmd, a) => execFileSync(cmd, a, { maxBuffer: 256 * 1024 * 1024 }).toString('utf8');
  const info = run('pdfinfo', [file]), pages = +(/Pages:\s+(\d+)/.exec(info) || [0, 0])[1], size = /Page size:\s+([\d.]+) x ([\d.]+) pts/.exec(info);
  const texts = run('pdftotext', ['-layout', file, '-']).split('\f').slice(0, pages);
  const links = run('pdfinfo', ['-url', file]).split('\n').slice(1).filter(l => /\S/.test(l)).map(l => l.trim().split(/\s+/).slice(2).join(' '));
  const images = run('pdfimages', ['-list', file]).split('\n').slice(2).filter(Boolean).map(l => l.trim().split(/\s+/)).map(c => ({ page: +c[0], w: +c[3], h: +c[4], ppi: Math.min(+c[12], +c[13]) }));
  // ink: every sheet as a small grey picture; the share of pixels that are not white
  const tmp = fs.mkdtempSync(path.join(outDir, 'ink-')); run('pdftoppm', ['-r', '18', '-gray', file, path.join(tmp, 'p')]);
  const ink = fs.readdirSync(tmp).sort().map(n => { const b = fs.readFileSync(path.join(tmp, n)), head = /^P5\s+\d+\s+\d+\s+\d+\s/.exec(b.subarray(0, 48).toString('latin1')), px = b.subarray(head ? head[0].length : 0);
    let dark = 0; for (const v of px) if (v < 215) dark++; return dark / (px.length || 1); });
  fs.rmSync(tmp, { recursive: true, force: true });
  // Firefox writes "fi" and "fl" as one sign each; a reader of the PDF finds the two letters all the same
  const plain = t => t.normalize('NFKC');
  return { pages, pts: size ? [Math.round(+size[1]), Math.round(+size[2])] : null, bytes: fs.statSync(file).size, texts: texts.map(plain), text: plain(texts.join('\n')).replace(/\s+/g, ' '), images, ink, links };
}
function checkPdf(tag, file, t, paper, page) {
  const f = pdfFacts(file), K = KIND[t.kind], want = PAPER[paper].pts;
  report.pdfs.push({ file: path.basename(file), pages: f.pages, bytes: f.bytes, images: f.images.length });
  check(tag + ': a PDF of ' + f.pages + ' sheets, ' + Math.round(f.bytes / 1024) + ' kB', f.pages > 0 && f.bytes > 20000, path.basename(file));
  check(tag + ': the sheets are ' + PAPER[paper].format + ', upright', f.pts && Math.abs(f.pts[0] - want[0]) <= 2 && Math.abs(f.pts[1] - want[1]) <= 2, f.pts && f.pts.join(' x '));
  const bounds = t.pages || K.pages; check(tag + ': the number of sheets is within bounds (' + bounds.join(' to ') + ')', f.pages >= bounds[0] && f.pages <= bounds[1], f.pages);
  const blank = f.ink.map((v, k) => v < 0.002 ? k + 1 : 0).filter(Boolean); check(tag + ': no sheet is blank', blank.length === 0 && f.ink.length === f.pages, blank.length ? 'sheets ' + blank.join(', ') : '');
  check(tag + ': the file stays under 5 MB', f.bytes < MAX_BYTES, Math.round(f.bytes / 1024) + ' kB');
  const first = (f.texts[0] || '').replace(/\s+/g, ' ');
  check(tag + ': the first sheet opens with what this is, when its numbers were computed and by which version', K.word.test(first) && /Numbers computed .{4,60} with FlowCellar \d/.test(first) && /PDF made \d/.test(first), first.slice(0, 200));
  // the selection line stands under the summary: on the first sheet, or at the head of the second when the verdict has much to say
  const lead = f.texts.slice(0, 2).join(' ').replace(/\s+/g, ' ');
  check(tag + ': the selection that is shown is said right after the summary, on the first or the second sheet', /This PDF shows/.test(lead) && K.shows.test(lead), (/This PDF shows.{0,160}/.exec(lead) || [''])[0]);
  if (page.headline) { const head = page.headline.normalize('NFKC').replace(/\s+/g, ' ').trim(); check(tag + ': the verdict sentence is text in the PDF, to select and to search', f.text.includes(head.slice(0, 60)), head.slice(0, 70)); }
  const left = NO_CONTROL.filter(w => f.text.includes(w)); check(tag + ': no control and no navigation is in the text', left.length === 0, left.join(' | '));
  check(tag + ': no link is in the PDF (a link would carry the address of the page: a path or a server of this computer)', f.links.length === 0, f.links.length + ' links, first ' + String(f.links[0]).replace(/[^/#?]{12,}/g, '…').slice(0, 80));
  // a heading at the foot of a sheet: the last line of a sheet is the title of a section
  const alone = []; f.texts.forEach((tx, k) => { const lines = tx.split('\n').map(l => l.trim()).filter(Boolean), last = lines[lines.length - 1] || ''; for (const title of page.titles) if (last.startsWith(title) && last.length <= title.length + 16) alone.push('sheet ' + (k + 1) + ': ' + title); });
  check(tag + ': no heading stands alone at the foot of a sheet', alone.length === 0, alone.join(' · '));
  const coarse = f.images.filter(im => im.w > 60 && im.ppi < 190); check(tag + ': pictures have at least twice the pixels of their size on the sheet' + (f.images.length ? '' : ' (charts are drawings here: no picture at all)'), coarse.length === 0, coarse.slice(0, 3).map(im => 'sheet ' + im.page + ' ' + im.w + 'x' + im.h + ' at ' + im.ppi + ' ppi').join(', '));
  return f;
}

// ---- in the page ----
const STATE = () => { const box = new Map(); for (const e of document.querySelectorAll('#app, #app *')) { const r = e.getBoundingClientRect(); box.set(e, [r.left, r.top + scrollY, r.width, r.height]); } window.__box = box;
  return { n: box.size, theme: document.documentElement.getAttribute('data-theme'), cls: document.documentElement.className, title: document.title, y: scrollY, tall: document.documentElement.scrollHeight,
    open: document.querySelectorAll('.fc-disclosure.is-open').length, canvases: [...document.querySelectorAll('canvas')].map(c => c.width + 'x' + c.height).join(','), bg: getComputedStyle(document.body).backgroundColor }; };
const MOVED = () => { let n = 0, first = ''; const now = [...document.querySelectorAll('#app, #app *')];
  for (const e of now) { const b = window.__box.get(e), r = e.getBoundingClientRect(), v = [r.left, r.top + scrollY, r.width, r.height]; if (b && !b.some((x, i) => Math.abs(x - v[i]) > 0.5)) continue; n++; if (!first) first = e.tagName.toLowerCase() + '.' + String(e.className).split(' ')[0] + (b ? ' ' + b.map(Math.round).join(',') + ' -> ' + v.map(Math.round).join(',') : ' is new'); }
  if (now.length !== window.__box.size) { n++; first = first || now.length + ' elements, ' + window.__box.size + ' before'; } return { n, first }; };
const PAPER_FACTS = W => { const app = document.getElementById('app').getBoundingClientRect(), shown = e => e.getClientRects().length > 0, wide = [];
  for (const e of document.querySelectorAll('#app *')) { if (!shown(e)) continue; const r = e.getBoundingClientRect(); if (r.width && (r.right > app.right + 0.6 || r.left < app.left - 0.6)) wide.push(e.tagName.toLowerCase() + '.' + String(e.className).split(' ')[0] + ' ' + Math.round(r.left - app.left) + '..' + Math.round(r.right - app.left)); }
  const chrome = [...document.querySelectorAll('.fc-appbar, .fc-rail, .fc-controls, .fc-crumbs, .fc-chart__menu, .fc-track__menu, .fc-chart__readout, .fc-track__value, .fc-tip, .fc-tray, .fc-facets, .fc-peek, select, input, a.fc-btn, button:not(.fc-legend__btn):not(.fc-samples .fc-linkbtn):not(.fc-disclosure), canvas.fck-over, canvas.fc-spmap__over')].filter(shown)
    .map(e => e.tagName.toLowerCase() + '.' + String(e.className).split(' ')[0] + ' "' + (e.textContent || '').trim().slice(0, 24) + '"');
  const charts = [...document.querySelectorAll('canvas.fck-base, canvas.fc-spmap__cv:not(.fc-spmap__over)')].filter(shown).map(c => { const r = c.getBoundingClientRect(); return { w: c.width, h: c.height, cw: r.width, ch: r.height }; });
  // the colour of the ground of the first chart, read through a canvas of one pixel that is made for reading
  const pix = (() => { const c = [...document.querySelectorAll('canvas.fck-base')].filter(shown)[0]; if (!c) return null; const p = document.createElement('canvas'); p.width = p.height = 1; const g = p.getContext('2d', { willReadFrequently: true }); g.drawImage(c, 1, 1, 1, 1, 0, 0, 1, 1); const d = g.getImageData(0, 0, 1, 1).data; return [d[0], d[1], d[2]]; })();
  // the colours a screen read from the theme while it was drawn: the swatches of the legends that name a state of the channels
  // ("No pore" is light blue on paper and dark slate on a dark screen) against the colour the theme of now gives that state
  const want = (() => { const p = document.createElement('i'); p.style.color = getComputedStyle(document.documentElement).getPropertyValue('--fc-stategroup-inactive'); document.body.append(p); const c = getComputedStyle(p).color; p.remove(); return c; })();
  const noPore = [...document.querySelectorAll('.fc-legend__item')].filter(e => shown(e) && e.textContent.trim() === 'No pore').map(e => getComputedStyle(e.querySelector('.fc-sw')).backgroundColor);
  const head = document.querySelector('.fc-paper__head'), scope = document.querySelector('.fc-paper__scope'), left = document.querySelector('.fc-paper__left');
  return { W, want, noPore, zoomed: [...document.querySelectorAll('#app table')].map(e => +e.style.zoom).filter(z => z > 0), app: [Math.round(app.left), Math.round(app.width)], sideways: document.documentElement.scrollWidth - document.documentElement.clientWidth, wide, chrome, charts, pix,
    cls: document.documentElement.className, theme: document.documentElement.getAttribute('data-theme'), bg: getComputedStyle(document.body).backgroundColor, color: getComputedStyle(document.body).color,
    head: head ? head.innerText.replace(/\s+/g, ' ') : null, scope: scope ? scope.innerText.replace(/\s+/g, ' ') : null, left: left ? left.textContent : null,
    closedCharts: document.querySelectorAll('.fc-disclosure[data-group="charts"]:not(.is-open)').length, title: document.title,
    // every open group of further charts with the number of charts drawn in it
    groups: [...document.querySelectorAll('.fc-disclosure[data-group="charts"].is-open')].map(d => { const body = d.nextElementSibling; return [d.textContent.slice(0, 24), body && body.classList.contains('fc-disclosure__body') ? body.querySelectorAll('canvas.fck-base').length : 0]; }),
    headline: (document.querySelector('.fc-verdict__headline') || {}).textContent || '', titles: [...document.querySelectorAll('.fc-section__title')].map(e => e.textContent.trim()).filter(t => t.length > 3) }; };
const contentPx = paper => Math.floor((PAPER[paper].mm[0] - 2 * MARGIN_MM) / 25.4 * 96);

async function inBrowser(name) {
  const { browser, log, newPage, version } = await launch(name);
  const facts = new Map();
  for (const t of chosen) {
    const K = KIND[t.kind], tag = name + ' ' + t.name, origin = t.url.replace(/[#?].*$/, '').replace(/[^/]*$/, '');
    const page = await newPage({ width: 1366, height: 900, scheme: t.dark ? 'dark' : 'light', allowed: u => u.startsWith(origin) });
    // the print dialog cannot be opened by a script that is watched: window.print is replaced by what a browser does first, the `beforeprint` event
    await page.addInitScript(() => { window.__prints = 0; window.print = () => { window.__prints++; window.dispatchEvent(new Event('beforeprint')); }; });
    await page.goto(t.url); await page.waitForSelector(K.ready, { timeout: 40000 }); await settle(page, 1200);
    if (t.dark) { await setTheme(page, 'Dark'); await settle(page, 400);
      // on the screen, too, the colours a screen read from the theme follow the switch (they did not before paper needed them to)
      const d = await page.evaluate(PAPER_FACTS, 0); if (d.noPore.length) check(tag + ': on screen, the dark theme reaches the colours of the charts', d.noPore.every(c => c === d.want), d.noPore.join(' ') + ' against ' + d.want); }
    // ---- the control ----
    const btn = page.locator('[data-fk="save-pdf"]');
    check(tag + ': the page has a control named "Save as PDF"', await btn.count() === 1 && (await btn.innerText()).trim() === 'Save as PDF' && await btn.isVisible(), await btn.count());
    check(tag + ': "PDF with every chart" is offered on a run page, and only there', (await page.locator('[data-fk="save-pdf-all"]').count()) === (K.button ? 1 : 0));
    await btn.click(); await settle(page, 200);
    const pressed = await page.evaluate(() => ({ prints: window.__prints, cls: document.documentElement.className }));
    check(tag + ': the control opens the print dialog with the layout for paper', pressed.prints === 1 && /fc-print/.test(pressed.cls), JSON.stringify(pressed));
    await page.evaluate(() => window.dispatchEvent(new Event('afterprint'))); await settle(page, 500);
    // ---- the screen as a person left it (scrolled down, perhaps dark), then the browser's own Print command: `beforeprint` without the control ----
    if (t.scroll) { await page.evaluate(y => window.scrollTo(0, y), t.scroll); await settle(page, 300); }
    const before = await page.evaluate(STATE);
    await page.evaluate(() => window.dispatchEvent(new Event('beforeprint')));
    // ---- the layout for paper, at the width of each sheet, as the browser lays it out for print ----
    for (const paper of papers) {
      const W = contentPx(paper);
      await page.emulateMedia({ media: 'print' }); await page.setViewportSize({ width: W, height: 1000 }); await settle(page, 250);
      const f = await page.evaluate(PAPER_FACTS, W), ptag = tag + ' ' + paper;
      if (paper === papers[0]) facts.set(t.name, f);
      check(ptag + ': the page is laid out for paper, ' + PAPER_W + ' px wide inside a sheet of ' + W + ' px', /fc-print/.test(f.cls) && f.app[1] === PAPER_W && f.app[1] <= W, JSON.stringify(f.app));
      check(ptag + ': nothing is wider than the sheet', f.wide.length === 0 && f.sideways <= 0, f.wide.slice(0, 3).join(' · ') || f.sideways);
      check(ptag + ': no navigation, control, button or read-out is left', f.chrome.length === 0, f.chrome.slice(0, 4).join(' · '));
      check(ptag + ': the head says what this is and the line under the summary says which selection is shown', !!f.head && K.word.test(f.head) && /Numbers computed .* with FlowCellar \d/.test(f.head) && !!f.scope && /This PDF shows/.test(f.scope) && K.shows.test(f.scope), (f.head || 'no head') + ' | ' + (f.scope || 'no selection line').slice(0, 120));
      if (K.button) check(ptag + ': the line says what was left out because it was closed', /^The page as it was shown\./.test(f.left || ''), f.left);
      const m = /(\d+), (\d+), (\d+)/.exec(f.bg), c = /(\d+), (\d+), (\d+)/.exec(f.color);
      check(ptag + ': the theme of paper, light, whatever the screen shows', f.theme === 'print' && m && +m[1] + +m[2] + +m[3] >= 750 && c && +c[1] + +c[2] + +c[3] < 200 && (!f.pix || f.pix[0] + f.pix[1] + f.pix[2] >= 720), f.bg + ' / ' + f.color + ' / chart ' + f.pix);
      if (f.noPore.length) check(ptag + ': the colours of the charts are those of paper too (a legend names “No pore” in the colour the paper theme gives it)', f.noPore.every(c => c === f.want), f.noPore.join(' ') + ' against ' + f.want);
      const off = f.charts.filter(q => q.w < 2 * q.cw - 2 || Math.abs(q.w / 3 - q.cw) > 1.5 || Math.abs(q.h / 3 - q.ch) > 1.5);
      check(ptag + ': every chart (' + f.charts.length + ') is painted again at the width it has on the sheet, three pixels to one', (K.charts === false || f.charts.length > 0) && off.length === 0, off.slice(0, 3).map(q => q.w + 'x' + q.h + ' for ' + Math.round(q.cw) + 'x' + Math.round(q.ch)).join(', '));
      if (t.fits) check(ptag + ': a table that is wider than the sheet is scaled down, to no less than half its size', f.zoomed.length > 0 && f.zoomed.every(z => z >= 0.5 && z < 1), f.zoomed.join(' '));
      if (t.map) check(ptag + ': the spatial map is on the sheet and says what it is coloured by', await page.locator('.fc-spmap__first canvas').count() >= 1 && /coloured by/.test(await page.locator('.fc-spmap .fc-chart__scope').innerText()));
      await page.screenshot({ path: path.join(outDir, name + '_' + t.name + '_' + paper + '_paper.png'), fullPage: true }); report.shots.push(name + '_' + t.name + '_' + paper + '_paper.png');
    }
    await page.emulateMedia({ media: null }); await page.setViewportSize({ width: 1366, height: 900 });
    await page.evaluate(() => window.dispatchEvent(new Event('afterprint'))); await settle(page, 600);
    const same = async what => { const now = await page.evaluate(STATE), moved = await page.evaluate(MOVED), diff = ['theme', 'cls', 'title', 'y', 'tall', 'open', 'canvases', 'bg'].filter(k => now[k] !== before[k]);
      check(tag + ': ' + what + ', every element of the page is where it was, with the same theme, title, open groups, charts and scroll position', moved.n === 0 && diff.length === 0, moved.n ? moved.n + ' moved, first ' + moved.first : diff.map(k => k + ' ' + String(before[k]).slice(0, 40) + ' -> ' + String(now[k]).slice(0, 40)).join(' · ')); };
    await same('after the dialog');
    // ---- the PDFs (Chromium): the browser says `beforeprint` and `afterprint` itself ----
    const f0 = facts.get(t.name);
    if (name === 'chromium') for (const paper of papers) {
      const file = path.join(outDir, name + '_' + t.name + '_' + paper + '.pdf');
      await page.pdf({ path: file, format: PAPER[paper].format, margin: { top: MARGIN_MM + 'mm', bottom: MARGIN_MM + 'mm', left: MARGIN_MM + 'mm', right: MARGIN_MM + 'mm' }, printBackground: false });
      checkPdf(name + ' ' + t.name + ' ' + paper, file, t, paper, f0);
      await settle(page, 400); await same('after the PDF on ' + PAPER[paper].format);
    }
    // ---- with every chart: the closed groups of further charts are opened for the print, and closed again ----
    if (t.every && K.button) {
      const closed = await page.locator('.fc-disclosure[data-group="charts"]:not(.is-open)').count();
      await page.locator('[data-fk="save-pdf-all"]').click(); await settle(page, 300);
      const e = await page.evaluate(PAPER_FACTS, 0);
      check(tag + ': "PDF with every chart" opens the ' + closed + ' closed groups of further charts for the print', closed > 0 && e.closedCharts === 0 && /^With every chart/.test(e.left || '') && e.charts.length > f0.charts.length, e.closedCharts + ' still closed · ' + e.charts.length + ' charts, ' + f0.charts.length + ' as shown · ' + e.left);
      check(tag + ': every group that was opened holds its charts (also a group whose charts are made only while it is open)', e.groups.length === closed + f0.groups.length && e.groups.every(g => g[1] > 0), e.groups.map(g => g[0] + ': ' + g[1]).join(' · '));
      const off = e.charts.filter(q => q.w < 2 * q.cw - 2 || Math.abs(q.w / 3 - q.cw) > 1.5);
      check(tag + ': the charts of the opened groups are painted at the width of the sheet too', off.length === 0 && e.wide.length === 0, off.length + ' of ' + e.charts.length + ' · ' + e.wide.slice(0, 2).join(' · '));
      if (name === 'chromium') { const file = path.join(outDir, name + '_' + t.name + '_every_chart_a4.pdf');
        await page.pdf({ path: file, format: 'A4', margin: { top: MARGIN_MM + 'mm', bottom: MARGIN_MM + 'mm', left: MARGIN_MM + 'mm', right: MARGIN_MM + 'mm' } });
        const g = checkPdf(name + ' ' + t.name + ' with every chart', file, t, 'a4', e), asShown = report.pdfs.find(p => p.file === name + '_' + t.name + '_a4.pdf');
        check(tag + ': the PDF with every chart is longer than the page as shown', asShown && g.pages > asShown.pages, g.pages + ' against ' + (asShown && asShown.pages)); }
      else await page.evaluate(() => window.dispatchEvent(new Event('afterprint')));
      await settle(page, 500); await same('after the PDF with every chart');
    }
    // a key press brings the screen back when a browser never says `afterprint`
    await page.evaluate(() => window.dispatchEvent(new Event('beforeprint'))); await page.keyboard.press('Shift'); await settle(page, 300);
    check(tag + ': a key press returns the screen when the browser never says that printing is over', !/fc-print/.test(await page.evaluate(() => document.documentElement.className)));
    await page.finish();
  }
  await browser.close();
  for (const e of log.errors) report.errors.push({ engine: name, ...e }); for (const r of log.requests) report.requests.push(r);
  return { version, facts };
}

let chromiumFacts = null;
for (const name of engines) { const r = await inBrowser(name); report[name] = r.version; if (name === 'chromium' || !chromiumFacts) chromiumFacts = r.facts; }

// ---- Firefox's own PDFs ----
if (engines.includes('firefox')) {
  let ff = null; try { ff = await firefoxSession(); } catch (e) { check('firefox: its print pipeline can be driven', false, e.message); }
  if (ff) { try {
    report.firefox_pdf = await ff.version();
    for (const t of chosen) {
      const K = KIND[t.kind];
      await ff.open(t.url, 'document.querySelector(' + JSON.stringify(K.ready) + ')', 1500);
      await ff.js('window.__ev = []; if (!window.__evOn) { window.__evOn = 1; addEventListener("beforeprint", () => window.__ev.push("before"), true); addEventListener("afterprint", () => window.__ev.push(document.documentElement.className), true); } return 1');
      for (const paper of papers) {
        const file = path.join(outDir, 'firefox_' + t.name + '_' + paper + '.pdf'); await ff.pdf(file, paper);
        const f = checkPdf('firefox ' + t.name + ' ' + paper, file, t, paper, chromiumFacts.get(t.name));
        check('firefox ' + t.name + ' ' + paper + ': the charts are pictures of three pixels to one', (K.charts === false || f.images.length > 0) && f.images.filter(im => im.w > 60).every(im => im.ppi >= 270), f.images.length + ' pictures, least ' + Math.min(...f.images.map(im => im.ppi)) + ' ppi');
      }
      const after = JSON.parse(await ff.js('return JSON.stringify({ ev: window.__ev, cls: document.documentElement.className, theme: document.documentElement.getAttribute("data-theme") })'));
      check('firefox ' + t.name + ': Firefox says `beforeprint` and `afterprint`, the layout for paper is on between them and off afterwards', after.ev.length === 2 * papers.length && after.ev.filter(v => /fc-print/.test(v)).length === papers.length && !/fc-print/.test(after.cls) && after.theme !== 'print', JSON.stringify(after));
    }
  } finally { await ff.close(); } } else if (ff === null) check('firefox: its print pipeline can be driven', false, 'no Playwright Firefox found');
}

report.finished = new Date().toISOString();
const failed = report.checks.filter(c => !c.ok);
report.summary = { checks: report.checks.length, failed: failed.length, errors: report.errors.length, requests_to_another_origin: report.requests.length, pdfs: report.pdfs.length };
fs.writeFileSync(path.join(outDir, 'print_report.json'), JSON.stringify(report, null, 1));
console.log('\n' + report.checks.length + ' checks, ' + failed.length + ' failed · ' + report.pdfs.length + ' PDFs in ' + outDir + ' · page and console errors and policy violations: ' + report.errors.length + ' · requests to another origin: ' + report.requests.length);
for (const c of failed.slice(0, 30)) console.log('  FAILED ' + c.name + '  [' + c.detail.slice(0, 200) + ']');
for (const e of report.errors.slice(0, 12)) console.log('  ' + e.engine + ' ' + e.kind + ': ' + e.text.slice(0, 300));
process.exit(failed.length || report.errors.length || report.requests.length ? 1 : 0);
