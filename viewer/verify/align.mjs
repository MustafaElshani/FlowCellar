// The charts drawn from alignment records and from the basecaller's tags, used like a person would, in headless
// Firefox: the Alignment section opens with the ring and its table, every chart has a title, a reading, a read-out
// under the pointer and under the arrow keys, and a table with rows; the page does not move while the pointer
// crosses a chart and never scrolls sideways; at two widths, in both themes, with all samples, with one sample and
// with two samples chosen. Nothing may be logged.
//
//   node verify/align.mjs REPORT.html [REPORT.html ...] [--shots DIR] [--widths 1366,2560]
//
// Exit code 0 when every check passed, 1 otherwise.
import fs from 'node:fs';
import path from 'node:path';
import { launch, fileUrl, settle, open, setTheme, sideways } from './browser.mjs';

const args = process.argv.slice(2), opt = (n, d) => { const i = args.indexOf('--' + n); return i >= 0 ? args.splice(i, 2)[1] : d; };
const widths = opt('widths', '1366,2560').split(',').map(Number), shots = opt('shots', ''), reports = args;
if (!reports.length) { console.error('usage: node verify/align.mjs REPORT.html [...] [--shots DIR] [--widths 1366,2560]'); process.exit(2); }
if (shots) fs.mkdirSync(shots, { recursive: true });
const SECTIONS = ['alignment', 'reads_over_time', 'transcripts'];
const OPENING = ['p-aln-records', 'p-aln-records-table', 'p-aln-bases', 'p-aln-strand', 'p-aln-depth', 'p-aln-contigs', 'p-aln-lengths', 'p-ident', 'p-aln-calibration'];
const MORE = ['p-mapq', 'p-aln-soft_clip', 'p-aln-identity_by_length', 'p-aln-length_vs_aligned', 'p-aln-identity_vs_aligned', 'p-aln-quality_vs_identity'];
const TRACKS = ['p-reads-output_over_time', 'p-reads-length_over_time', 'p-reads-quality_over_time', 'p-reads-aligned_length_over_time', 'p-reads-identity_over_time', 'p-reads-speed_over_time'];
const RNA = ['p-rna-polya', 'p-rna-body_coverage', 'p-rna-saturation', 'p-rna-reads_per_transcript'];
const MINE = SECTIONS.map(s => '#s-' + s).concat(['#s-samples']);
const CHARTS = MINE.map(s => s + ' .fc-chart, ' + s + ' .fc-track').join(', '), PLOTS = MINE.map(s => s + ' .fc-plot').join(', ');
const { browser, log, newPage } = await launch();
let failed = 0;
const check = (name, ok, detail) => { if (!ok) failed++; console.log((ok ? 'ok   ' : 'FAIL ') + name + (detail != null && detail !== '' ? '  [' + String(detail).slice(0, 300) + ']' : '')); };
// where every chart and track is: nothing of this may change while the pointer moves
const boxes = page => page.evaluate(sel => [...document.querySelectorAll(sel)].map(e => { const r = e.getBoundingClientRect(); return [Math.round(r.left), Math.round(r.top + scrollY), Math.round(r.width), Math.round(r.height)].join(','); }).join(';'), CHARTS + ', .fc-section');
const has = (page, id) => page.evaluate(i => !!document.getElementById(i), id);

async function look(page, name, { samples, coarse }) {
  for (const sec of SECTIONS) check(name + ': section ' + sec, await page.locator('#s-' + sec).count() === 1);
  const missing = []; for (const id of OPENING.concat(MORE, TRACKS, RNA)) if (!(await has(page, id))) missing.push(id);
  // with more than twelve samples, a chart of one sample that is kept for the run only says so instead of drawing
  const allowed = coarse && samples ? new Set(['p-aln-depth', 'p-aln-calibration', 'p-aln-length_vs_aligned', 'p-aln-identity_vs_aligned', 'p-aln-quality_vs_identity', 'p-aln-identity_by_length', 'p-rna-saturation', 'p-rna-reads_per_transcript']) : new Set();
  check(name + ': every chart is there', missing.filter(id => !allowed.has(id)).length === 0, missing.filter(id => !allowed.has(id)).join(' '));
  // the section opens with the ring, and its table stands beside it
  const order = await page.evaluate(() => [...document.querySelectorAll('#s-alignment .fc-chart')].slice(0, 2).map(e => e.id));
  check(name + ': the Alignment section opens with the ring and its table', order.join(' ') === 'p-aln-records p-aln-records-table', order.join(' '));
  const ringRows = await page.locator('#p-aln-records-table tbody tr').count();
  check(name + ': the table of the ring has its seven rows', ringRows === 7, ringRows);
  if (!samples && !coarse) check(name + ': reads and bases per sample stand with the Samples table', await has(page, 'p-yield-per_sample'));
  const plots = await page.locator(PLOTS).all();
  check(name + ': charts drawn', plots.length >= (coarse && samples ? 10 : 20), plots.length);
  let moved = 0, said = 0, keyed = 0, sliders = 0;
  for (const p of plots) {
    await p.scrollIntoViewIfNeeded(); const before = await boxes(page), bb = await p.boundingBox(); if (!bb) continue;
    const slider = (await p.getAttribute('role')) === 'slider'; if (slider) sliders++;
    let got = false;
    for (const f of [0.3, 0.55, 0.8]) { await page.mouse.move(bb.x + bb.width * f, bb.y + bb.height * 0.5); await settle(page, 40); if ((await boxes(page)) !== before) moved++; const t = await p.getAttribute('aria-valuetext'); if (t && t !== 'no value selected') got = true; }
    if (got) said++;
    await page.mouse.move(2, 2); await settle(page, 20);
    if (slider) { await p.focus(); await page.keyboard.press('ArrowRight'); await settle(page, 30); const a = await p.getAttribute('aria-valuetext'); await page.keyboard.press('End'); await settle(page, 30); const z = await p.getAttribute('aria-valuetext'); if (a && a !== 'no value selected' && z && z !== 'no value selected') keyed++; if ((await boxes(page)) !== before) moved++; await page.keyboard.press('Escape'); }
  }
  check(name + ': the page does not move under the pointer or the arrow keys', moved === 0, moved + ' moves over ' + plots.length + ' charts');
  check(name + ': every chart reads out the value under the pointer', said === sliders, said + ' of ' + sliders);
  check(name + ': every chart reads out under the arrow keys', keyed === sliders, keyed + ' of ' + sliders);
  // the ring: at most six parts, each read with its count and its share
  const ring = page.locator('#p-aln-records .fc-plot'); await ring.focus(); const parts = new Set();
  for (let i = 0; i < 9; i++) { await page.keyboard.press(i ? 'ArrowRight' : 'Home'); await settle(page, 20); parts.add(await ring.getAttribute('aria-valuetext')); }
  check(name + ': the ring has at most six parts, each with count and share', parts.size >= 2 && parts.size <= 6 && [...parts].every(t => /records|reads/.test(t) && /%/.test(t)), [...parts].join(' | '));
  await page.keyboard.press('Escape');
  // title, reading, table
  let noTitle = 0, noReading = 0, noRows = 0, tables = 0;
  for (const c of await page.locator(MINE.map(s => s + ' figure.fc-chart').join(', ')).all()) {
    if (!((await c.locator('.fc-chart__title').innerText()).trim())) noTitle++;
    if (!(await c.locator('.fc-chart__reading').count()) || (await c.locator('.fc-chart__reading').innerText()).trim().length < 20) noReading++;
    const btn = c.locator('button:has-text("Show as table")'); if (!(await btn.count())) continue;
    await btn.click(); await settle(page, 30); tables++; if (!(await c.locator('.fc-charttable tbody tr').count())) noRows++; await c.locator('button:has-text("Hide table")').click();
  }
  for (const t of await page.locator('#s-reads_over_time .fc-track:not(.fc-track--na), #s-alignment .fc-track').all()) { const btn = t.locator('button:has-text("Show as table")'); if (!(await btn.count())) { noRows++; continue; } await btn.click(); await settle(page, 30); tables++; if (!(await t.locator('.fc-charttable tbody tr').count())) noRows++; await t.locator('button:has-text("Hide table")').click(); }
  check(name + ': every chart has a title and a reading written from its numbers', noTitle === 0 && noReading === 0, noTitle + ' without a title, ' + noReading + ' without a reading');
  check(name + ': “Show as table” gives rows for every chart', noRows === 0 && tables >= (coarse && samples ? 8 : 18), noRows + ' empty of ' + tables);
  check(name + ': no sideways scrolling', (await sideways(page)) <= 0, await sideways(page));
}

for (const report of reports) {
  const tag = path.basename(report).replace(/\..*$/, '');
  for (const width of widths) for (const scheme of ['light', 'dark']) {
    const page = await newPage({ width, height: 900, scheme });
    await open(page, fileUrl(report));
    const key = await page.evaluate(() => JSON.parse(document.getElementById('fc-boot').textContent).run);
    const base = fileUrl(report) + '#/run/' + key + '?more=aln-more';
    await page.goto(base); await settle(page, 600);
    const info = await page.evaluate(() => ({ samples: [...document.querySelectorAll('#s-samples tbody tr button.fc-samples__name')].length }));
    const coarse = info.samples > 12;
    await look(page, tag + ' ' + width + ' ' + scheme, { samples: 0, coarse });
    if (shots) for (const sec of SECTIONS.concat(['samples'])) { const el = await page.$('#s-' + sec); if (el) await el.screenshot({ path: path.join(shots, `${tag}_${sec}_${width}_${scheme}.png`) }); }
    if (scheme === 'light' && info.samples > 1) {
      // one sample, then two: the charts follow the Samples control
      await page.goto(base + '&s=s1'); await settle(page, 500);
      await look(page, tag + ' ' + width + ' one sample', { samples: 1, coarse });
      if (shots) for (const sec of SECTIONS) { const el = await page.$('#s-' + sec); if (el) await el.screenshot({ path: path.join(shots, `${tag}_${sec}_${width}_one.png`) }); }
      await page.goto(base + '&s=s1,s2'); await settle(page, 500);
      await look(page, tag + ' ' + width + ' two samples', { samples: 2, coarse });
      if (!coarse) check(tag + ' ' + width + ' two samples: one depth track per sample', await page.locator('#s-alignment .fc-tracks .fc-track').count() === 2, await page.locator('#s-alignment .fc-tracks .fc-track').count());
      if (shots) for (const sec of SECTIONS) { const el = await page.$('#s-' + sec); if (el) await el.screenshot({ path: path.join(shots, `${tag}_${sec}_${width}_two.png`) }); }
    }
    await page.finish();
  }
}
check('nothing was logged (console errors, page errors, policy violations)', log.errors.length === 0, JSON.stringify(log.errors).slice(0, 600));
check('no request left the page', log.requests.length === 0, log.requests.slice(0, 3).join(' '));
await browser.close();
console.log(failed ? failed + ' check(s) FAILED' : 'all checks passed');
process.exit(failed ? 1 : 0);
