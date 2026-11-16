// The Pychopper and RNA composition sections of one-run reports, used like a person would, in headless Firefox:
// both sections are there, every chart has a title, a reading, a read-out under the pointer and a table with rows,
// the page does not move while the pointer crosses a chart and never scrolls sideways, at two widths and in both
// themes, with all samples and with some samples chosen. Nothing may be logged.
//
//   node verify/libprep.mjs REPORT.html [REPORT.html ...] [--shots DIR] [--widths 1366,2560] [--pick 3]
//
//   --shots  also write one picture of each section per width and theme
//   --pick   how many samples to choose for the second pass (0 = skip it)
// Exit code 0 when every check passed, 1 otherwise.
import fs from 'node:fs';
import path from 'node:path';
import { launch, fileUrl, settle, open, setTheme, sideways } from './browser.mjs';

const args = process.argv.slice(2), opt = (n, d) => { const i = args.indexOf('--' + n); return i >= 0 ? args.splice(i, 2)[1] : d; };
const widths = opt('widths', '1366,2560').split(',').map(Number), shots = opt('shots', ''), pickN = +opt('pick', 3), reports = args;
if (!reports.length) { console.error('usage: node verify/libprep.mjs REPORT.html [...] [--shots DIR] [--widths 1366,2560] [--pick 3]'); process.exit(2); }
if (shots) fs.mkdirSync(shots, { recursive: true });
const SECTIONS = ['pychopper', 'rna_composition'], CHARTS = SECTIONS.map(s => '#s-' + s + ' .fc-chart').join(', '), PLOTS = SECTIONS.map(s => '#s-' + s + ' .fc-plot').join(', ');
const { browser, log, newPage } = await launch();
let failed = 0;
const check = (name, ok, detail) => { if (!ok) failed++; console.log((ok ? 'ok   ' : 'FAIL ') + name + (detail != null && detail !== '' ? '  [' + String(detail).slice(0, 200) + ']' : '')); };
// where every chart and both sections are: nothing of this may change while the pointer moves
const boxes = page => page.evaluate(sel => [...document.querySelectorAll(sel)].map(e => { const r = e.getBoundingClientRect(); return [Math.round(r.left), Math.round(r.top + scrollY), Math.round(r.width), Math.round(r.height)].join(','); }).join(';'), CHARTS + ', #s-pychopper, #s-rna_composition');

async function look(page, tag, width, pass) {
  const name = tag + ' ' + width + (pass ? ' ' + pass : '');
  for (const sec of SECTIONS) check(name + ': section ' + sec, await page.locator('#s-' + sec).count() === 1);
  const charts = await page.locator(CHARTS).all(), plots = await page.locator(PLOTS).all();
  check(name + ': charts drawn', charts.length >= 6, charts.length);
  let moved = 0, said = 0;
  for (const p of plots) {
    await p.scrollIntoViewIfNeeded(); const before = await boxes(page), bb = await p.boundingBox(); if (!bb) continue;
    for (const f of [0.3, 0.55, 0.8]) { await page.mouse.move(bb.x + bb.width * f, bb.y + bb.height * 0.5); await settle(page, 50); if ((await boxes(page)) !== before) moved++; }
    const text = await p.getAttribute('aria-valuetext'); if (text && text !== 'no value selected') said++;
    await page.mouse.move(2, 2);
  }
  check(name + ': the page does not move under the pointer', moved === 0, moved + ' moves over ' + plots.length + ' charts');
  check(name + ': every chart reads out the value under the pointer', said === plots.length, said + ' of ' + plots.length);
  let noTitle = 0, noReading = 0, noRows = 0;
  for (const ch of charts) {
    if (!(await ch.locator('.fc-chart__title').innerText()).trim()) noTitle++;
    if (!(await ch.locator('.fc-chart__reading').count()) || !(await ch.locator('.fc-chart__reading').innerText()).trim()) noReading++;
    const btn = ch.locator('button:has-text("Show as table")'); if (!(await btn.count())) { noRows++; continue; }
    await btn.click(); await settle(page, 40); if (!(await ch.locator('.fc-charttable tbody tr').count())) noRows++;
    await ch.locator('button:has-text("Hide table")').click();
  }
  check(name + ': every chart has a title and a reading written from its numbers', noTitle === 0 && noReading === 0, noTitle + ' without a title, ' + noReading + ' without a reading');
  check(name + ': every chart opens as a table with rows', noRows === 0, noRows);
  const words = await page.evaluate(sel => [...document.querySelectorAll(sel)].map(e => e.innerText).join(' '), SECTIONS.map(s => '#s-' + s).join(', '));
  check(name + ': no number reads as NaN, undefined or null', !/\bNaN\b|undefined|\bnull\b|Infinity/.test(words), (/.{0,30}(NaN|undefined|null|Infinity).{0,30}/.exec(words) || [''])[0]);
}

for (const rep of reports) {
  const tag = path.basename(rep).replace(/\.flowcellar\.html$/, '');
  for (const width of widths) {
    const page = await newPage({ width, height: 1000 });
    await open(page, fileUrl(rep)); await page.waitForSelector('.fc-runmain', { timeout: 20000 }); await settle(page, 300);
    for (const theme of ['Light', 'Dark']) {
      await setTheme(page, theme);
      check(tag + ' ' + width + ' ' + theme + ': no sideways scroll', (await sideways(page)) <= 0, await sideways(page));
      if (shots) {  // for the pictures only: bars that stick to the top of the window would lie over a section taller than the window
        const stuck = await page.evaluate(() => { const o = []; for (const e of document.querySelectorAll('*')) { const p = getComputedStyle(e).position; if (p === 'sticky' || p === 'fixed') { o.push(e); e.dataset.fcWas = e.style.position; e.style.position = 'static'; } } return o.length; });
        for (const sec of SECTIONS) { const el = page.locator('#s-' + sec); if (await el.count()) await el.screenshot({ path: path.join(shots, [tag, sec, width, theme.toLowerCase()].join('_') + '.png') }); }
        if (stuck) await page.evaluate(() => { for (const e of document.querySelectorAll('[data-fc-was]')) { e.style.position = e.dataset.fcWas; delete e.dataset.fcWas; } });
      }
    }
    await look(page, tag, width, '');
    if (pickN > 0) {  // some samples chosen, through the address like the Samples control writes it
      const keys = await page.evaluate(n => { const seen = new Set(); for (const m of document.querySelectorAll('#s-samples tbody tr button[aria-pressed]')) seen.add(m); return seen.size >= n ? n : 0; }, pickN);
      if (keys) {
        await page.evaluate(n => { location.hash = location.hash.split('?')[0] + '?s=' + Array.from({ length: n }, (_, i) => 's' + (i + 1)).join(','); }, Math.min(pickN, 3)); await settle(page, 400);
        await look(page, tag, width, 'with ' + Math.min(pickN, 3) + ' samples chosen');
        check(tag + ' ' + width + ': no sideways scroll with samples chosen', (await sideways(page)) <= 0);
      }
    }
    await page.finish();
  }
}
await browser.close();
check('nothing was logged (console, page errors, policy)', log.errors.length === 0, JSON.stringify(log.errors).slice(0, 400));
check('no request to another origin', log.requests.length === 0, log.requests.slice(0, 3).join(' '));
process.exit(failed ? 1 : 0);
