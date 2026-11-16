// The screens of the single-cell and spatial add-on, used like a person would, in headless Firefox: the Analyses
// tab, the page of an analysis (both pipelines of a sample), the panel on a run page and the block on Home, at
// several widths and in both themes, with screenshots. It works on any hub that holds analyses (it reads the
// hub's own index.js to know what to expect), so the same script checks the synthetic fixture and a real library.
//
//   node verify/analysis.mjs --hub DIR --out DIR [--widths 1366,2560] [--quick] [--report FILE.html]
//
// Exit code 0 when every check passed and nothing was logged, 1 otherwise. Writes OUT/analysis_report.json.
import fs from 'node:fs';
import path from 'node:path';
import zlib from 'node:zlib';
import { launch, fileUrl, settle, open, sideways } from './browser.mjs';

const args = process.argv.slice(2), opt = n => { const i = args.indexOf('--' + n); return i >= 0 ? args[i + 1] : null; }, flag = n => args.includes('--' + n);
const hub = opt('hub'), outDir = path.resolve(opt('out') || 'verify-out'), quick = flag('quick'), reportFile = opt('report');
if (!hub) { console.error('usage: node verify/analysis.mjs --hub DIR --out DIR [--widths 1366,2560] [--quick] [--report FILE]'); process.exit(2); }
const hubDir = path.resolve(hub), base = fileUrl(path.join(hubDir, 'index.html')), origin = fileUrl(hubDir);
const widths = (opt('widths') || '1366,2560').split(',').map(Number), themes = quick ? ['light'] : ['light', 'dark'];
fs.mkdirSync(outDir, { recursive: true });

const text0 = fs.readFileSync(path.join(hubDir, 'index.js'), 'utf8'), m0 = /^__FC_INDEX__\("([A-Za-z0-9+/=]*)"\);\n$/.exec(text0);
const ix = JSON.parse(zlib.gunzipSync(Buffer.from(m0[1], 'base64')).toString('utf8')), analyses = ix.analyses || [];
if (!analyses.length) { console.error('this hub holds no analysis'); process.exit(2); }
const runIds = new Set((() => { const c = ix.runs.columns.id; return Array.isArray(c) ? c : c.codes.map(k => c.dict[k]); })());
const of = (pipe, assay) => analyses.find(a => a.pipeline === pipe && a.assay === assay);
// one page per kind that the hub holds, a pair first
const kinds = [['spaceranger', 'visium_hd'], ['wf-single-cell', 'visium_hd'], ['spaceranger', 'visium'], ['wf-single-cell', 'visium'], ['wf-single-cell', 'singlecell']].map(k => of(...k)).filter(Boolean);
const linkedRun = (analyses.flatMap(a => a.links || []).find(l => l.run && runIds.has(l.run)) || {}).run || null;
// the flow cell runs of a sample: those of the analysis and of the other pipeline of the same sample, each once
const runsOf = a => new Set([a].concat((a.siblings || []).map(id => analyses.find(x => x.id === id)).filter(Boolean)).flatMap(r => (r.links || []).map(l => l.run || l.fc))).size;

const { browser, log, newPage, version } = await launch();
const report = { hub, browser: version, started: new Date().toISOString(), analyses: analyses.length, checks: [], shots: [], errors: log.errors, requests: log.requests };
const check = (name, ok, detail) => { report.checks.push({ name, ok: !!ok, detail: detail == null ? '' : String(detail).slice(0, 300) }); console.log((ok ? 'ok   ' : 'FAIL ') + name + (detail != null && detail !== '' ? '  [' + String(detail).slice(0, 170) + ']' : '')); };
const shot = async (page, name, full) => { if (full) { await page.evaluate(() => window.scrollTo(0, 0)); await settle(page, 150); } await page.screenshot({ path: path.join(outDir, name + '.png'), fullPage: !!full }); report.shots.push(name + '.png'); };
// the left and right edge of what a page shows, against the application bar: the same two lines on every screen
const edges = page => page.evaluate(() => { const bar = document.querySelector('.fc-appbar__in'), pg = document.querySelector('.fc-page'); if (!bar || !pg) return null; const cs = getComputedStyle(pg), b = bar.getBoundingClientRect(), p = pg.getBoundingClientRect(), bs = getComputedStyle(bar);
  return { barL: Math.round(b.left + parseFloat(bs.paddingLeft)), barR: Math.round(b.right - parseFloat(bs.paddingRight)), pageL: Math.round(p.left + parseFloat(cs.paddingLeft)), pageR: Math.round(p.right - parseFloat(cs.paddingRight)), w: innerWidth }; });
// every box of the page that could move when a read-out changes: position and size of the figures, sections and tiles
const boxes = page => page.evaluate(() => [...document.querySelectorAll('.fc-chart, .fc-section, .fc-tile, .fc-spmap__plot, .fc-chart__readout, .fc-spmap__side')].map(e => { const r = e.getBoundingClientRect(); return [Math.round(r.left), Math.round(r.top + scrollY), Math.round(r.width), Math.round(r.height)].join(','); }).join(';'));
// how much of a canvas is painted. It is read through a small copy made for reading: reading the page's own canvas again and again makes
// Chromium write a warning to the console (it would keep that canvas off the graphics card from then on), and a check must not change the page
const inked = (page, sel) => page.evaluate(s => { const c = document.querySelector(s); if (!c || !c.width) return 0; const p = document.createElement('canvas'); p.width = Math.min(240, c.width); p.height = Math.min(240, c.height);
  const g = p.getContext('2d', { willReadFrequently: true }); g.drawImage(c, 0, 0, p.width, p.height); const d = g.getImageData(0, 0, p.width, p.height).data; let n = 0; for (let i = 3; i < d.length; i += 4 * 97) if (d[i]) n++; return n; }, sel);

for (const width of widths) for (const scheme of themes) {
  const tag = width + '_' + scheme, page = await newPage({ width, height: width >= 2000 ? 1300 : 768, scheme, allowed: u => u.startsWith(origin) });
  // ---- the Analyses tab ----
  await open(page, base + '#/analyses'); await page.waitForSelector('.fc-analyses tbody tr', { timeout: 20000 }); await settle(page);
  const nRows = await page.locator('.fc-analyses tbody tr').count(), tab = await page.locator('.fc-tabs .fc-tab.is-on').innerText();
  check(tag + ' Analyses: the tab is lit and lists every sample once', /Analyses/.test(tab) && nRows >= 1 && nRows <= analyses.length, nRows + ' rows for ' + analyses.length + ' analyses');
  check(tag + ' Analyses: nothing is wider than the window', (await sideways(page)) <= 0, await sideways(page));
  const e1 = await edges(page); check(tag + ' Analyses: the page starts and ends on the lines of the application bar', e1 && e1.barL === e1.pageL && e1.barR === e1.pageR, JSON.stringify(e1));
  const head = (await page.locator('.fc-analyses thead').innerText()).replace(/\s+/g, ' ');
  check(tag + ' Analyses: the columns are sample, kind, chemistry, units, genes, UMIs, saturation, verdict, runs', ['Sample', 'Kind', 'Units', 'Genes', 'UMIs', 'Saturation', 'Verdict', 'Flow cell runs'].every(w => head.includes(w)), head);
  if (analyses.some(a => (a.links || []).some(l => !l.run || !runIds.has(l.run)))) check(tag + ' Analyses: a run that is not in the library is named and says so', await page.locator('.fc-runchip.is-missing').count() > 0);
  await shot(page, 'analyses_' + tag, true);
  // ---- Home ----
  await open(page, base + '#/home'); await settle(page);
  if (ix.runs.n) { check(tag + ' Home: has a short Analyses block', await page.locator('.fc-home__block--analyses').count() === 1 && await page.locator('.fc-home__block--analyses tbody tr').count() <= 4);
    check(tag + ' Home: nothing is wider than the window', (await sideways(page)) <= 0); await page.locator('.fc-home__block--analyses').scrollIntoViewIfNeeded(); await shot(page, 'home_analyses_' + tag); }
  // ---- one page per kind ----
  for (const a of kinds) {
    const name = a.pipeline.replace(/[^a-z]/g, '') + '_' + a.assay;
    await open(page, base + '#/analysis/' + a.id); await page.waitForSelector('.fc-analysis .fc-verdict__headline', { timeout: 20000 });
    await page.waitForSelector('.fc-analysis .fc-runmain .fc-section', { timeout: 30000 }); await settle(page, 500);
    check(tag + ' ' + name + ': header, verdict and six key numbers', await page.locator('.fc-analysis .fc-keynumbers .fc-tile').count() === 6 && (await page.locator('.fc-analysis .fc-verdict__detail').first().innerText()).includes('does not change the verdict'), await page.locator('.fc-analysis .fc-keynumbers .fc-tile').count());
    check(tag + ' ' + name + ': nothing is wider than the window', (await sideways(page)) <= 0, await sideways(page));
    const e2 = await edges(page); check(tag + ' ' + name + ': equal page edges', e2 && e2.barL === e2.pageL && e2.barR === e2.pageR, JSON.stringify(e2));
    await shot(page, 'analysis_' + name + '_' + tag, true);
    // an analysis made from several flow cells is one analysis of their merged reads: the header and the sentence under the verdict say so
    { const n = runsOf(a), line = (await page.locator('.fc-analysis .fc-anruns').innerText()).replace(/\s+/g, ' '), detail = (await page.locator('.fc-analysis .fc-verdict__detail').first().innerText()).replace(/\s+/g, ' ');
      if (n > 1) check(tag + ' ' + name + ': made from ' + n + ' flow cells: the header and the verdict say that it is one analysis of their merged reads', line.startsWith('Reads from ' + n + ' flow cell runs, merged before the analysis:') && detail.includes('one analysis of the reads of ' + (n === 2 ? 'both' : 'all ' + n) + ' flow cells together'), line.slice(0, 70) + ' | ' + detail);
      else check(tag + ' ' + name + ': made from one flow cell: nothing speaks of merged reads', !/merged/.test(line) && !/together/.test(detail), line.slice(0, 70) + ' | ' + detail); }
    const hasMap = await page.locator('.fc-spmap').count();
    if (a.assay !== 'singlecell') check(tag + ' ' + name + ': a spatial map, or a line that says why there is none', hasMap === 1 || await page.locator('#s-map .fc-empty').count() === 1);
    if (hasMap) {
      await page.locator('.fc-spmap').scrollIntoViewIfNeeded(); await settle(page, 300);
      check(tag + ' ' + name + ': the map is drawn', (await inked(page, '.fc-spmap__first .fc-spmap__plot canvas')) > 20, await inked(page, '.fc-spmap__first .fc-spmap__plot canvas'));
      // pointer over the map: the read-out changes its words and nothing changes its place or size
      const before = await boxes(page), idle = await page.locator('.fc-spmap .fc-ro').innerText(), bb = await page.locator('.fc-spmap__first .fc-spmap__plot').boundingBox();
      const t0 = Date.now(); let moved = 0;
      for (let k = 0; k < 40; k++) { await page.mouse.move(bb.x + bb.width * (0.2 + 0.6 * k / 40), bb.y + bb.height * (0.3 + 0.4 * ((k * 7) % 10) / 10)); moved++; }
      await settle(page, 80); const over = await page.locator('.fc-spmap .fc-ro').innerText(), during = await boxes(page), ms = Date.now() - t0;
      check(tag + ' ' + name + ': the read-out names the cell under the pointer', over !== idle && /(bin|spot row) \d+, \d+/.test(over), over);
      check(tag + ' ' + name + ': no layout shift on hover', before === during, before === during ? '' : 'boxes moved');
      check(tag + ' ' + name + ': the map stays responsive (40 pointer moves)', ms < 4000, ms + ' ms');
      await shot(page, 'analysis_' + name + '_map_' + tag);
      await page.mouse.move(2, 2); await settle(page, 80);
      check(tag + ' ' + name + ': the read-out goes back when the pointer leaves', (await page.locator('.fc-spmap .fc-ro').innerText()) === idle);
      // keyboard: the map takes focus, the arrow keys move a cursor, the sentence is the value of the slider
      await page.locator('.fc-spmap__first .fc-spmap__plot[role="slider"]').focus(); await page.keyboard.press('ArrowRight'); await page.keyboard.press('ArrowDown'); await settle(page, 60);
      const said = await page.locator('.fc-spmap__first .fc-spmap__plot[role="slider"]').getAttribute('aria-valuetext');
      check(tag + ' ' + name + ': arrow keys read a cell', /(bin|spot row) \d+, \d+/.test(said || ''), said);
      check(tag + ' ' + name + ': the keyboard cursor moves nothing either', (await boxes(page)) === before);
      // the switch: every channel can be chosen, and the picture is painted again
      const opts = page.locator('.fc-spmap__opt'), nOpt = await opts.count(); let painted = 0;
      for (let k = 0; k < nOpt; k++) { await opts.nth(k).click(); await settle(page, 120); if ((await inked(page, '.fc-spmap__first .fc-spmap__plot canvas')) > 20) painted++; if (k === 1) await shot(page, 'analysis_' + name + '_map2_' + tag); }
      check(tag + ' ' + name + ': colour by ' + nOpt + ' things, each drawn', nOpt >= 2 && painted === nOpt, painted + ' of ' + nOpt);
      await opts.nth(0).click();
      // a wide window: a second picture of the same slide beside the first, coloured by something else, with the same cursor
      const second = page.locator('.fc-spmap__second .fc-spmap__plot'), secondShown = (await second.count()) === 1 && await second.isVisible();
      if (width >= 1800 && nOpt >= 2) {
        check(tag + ' ' + name + ': a second picture stands beside the first and is drawn', secondShown && (await inked(page, '.fc-spmap__second .fc-spmap__plot canvas')) > 20);
        if (secondShown) {
          const b1 = await page.locator('.fc-spmap__first .fc-spmap__plot').boundingBox(), b2 = await second.boundingBox(), was = await boxes(page);
          check(tag + ' ' + name + ': the two pictures have one size and stand on one line', Math.abs(b1.width - b2.width) <= 1 && Math.abs(b1.y - b2.y) <= 1 && b2.x > b1.x + b1.width - 1, [b1.width, b2.width, b1.y, b2.y].map(Math.round).join(' '));
          await page.mouse.move(b1.x + b1.width * 0.5, b1.y + b1.height * 0.5); await settle(page, 80);
          const there = await second.getAttribute('aria-valuetext'), here = await page.locator('.fc-spmap__first .fc-spmap__plot').getAttribute('aria-valuetext');
          check(tag + ' ' + name + ': the cell pointed at in one picture is marked in the other', /(bin|spot row) \d+, \d+/.test(there || '') && (there || '').split(',').slice(0, 2).join(',') === (here || '').split(',').slice(0, 2).join(','), there);
          check(tag + ' ' + name + ': the second picture moves nothing either', (await boxes(page)) === was);
          await page.mouse.move(2, 2); await settle(page, 80);
        }
      } else check(tag + ' ' + name + ': one picture with its switch at this width', !secondShown);
      await page.locator('.fc-spmap .fc-chart__menu button').click(); await settle(page, 100);
      check(tag + ' ' + name + ': Show as table gives the numbers behind the picture', await page.locator('.fc-spmap table tbody tr').count() >= 2 && /exact/.test(await page.locator('.fc-spmap table').innerText()));
      await page.locator('.fc-spmap .fc-chart__menu button').click();
    }
    // a chart of the kit on this page: hover must not shift it either
    const plot = page.locator('.fc-analysis .fc-chart:not(.fc-spmap) .fc-plot').first();
    if (await plot.count()) { await plot.scrollIntoViewIfNeeded(); await settle(page, 200); const b0 = await boxes(page), pb = await plot.boundingBox(); await page.mouse.move(pb.x + pb.width * 0.5, pb.y + pb.height * 0.5); await settle(page, 120);
      check(tag + ' ' + name + ': no layout shift on hover of a chart', (await boxes(page)) === b0); await page.mouse.move(2, 2); }
    // the other pipeline of the same sample is a tab of the same page
    if ((a.siblings || []).length) {
      const seg = page.locator('.fc-ancontrols [role="radio"]'); check(tag + ' ' + name + ': the two pipelines are two tabs of one page', await seg.count() === 2);
      await page.locator('.fc-ancontrols [role="radio"][aria-checked="false"]').click(); await page.waitForSelector('.fc-analysis .fc-runmain .fc-section', { timeout: 30000 }); await settle(page, 500);
      check(tag + ' ' + name + ': the other tab draws', await page.locator('.fc-analysis .fc-keynumbers .fc-tile').count() === 6 && (await sideways(page)) <= 0);
      check(tag + ' ' + name + ': two pipelines compared', await page.locator('#s-pair').count() === 1);
    }
  }
  // ---- the run page: one panel that lists the analyses made from the run ----
  if (linkedRun) {
    await open(page, base + '#/run/' + linkedRun + '?at=analyses'); await page.waitForSelector('#s-analyses .fc-ancard', { timeout: 30000 }); await settle(page, 600);
    const t = (await page.locator('#s-analyses').innerText()).replace(/\s+/g, ' ');
    check(tag + ' run page: a panel lists the analyses made from the run and says whose verdict it is', /never change this flow cell’s verdict/.test(t) && /Open analysis/.test(t), t.slice(0, 120));
    check(tag + ' run page: the rail names the panel', /Single-cell and spatial/.test(await page.locator('.fc-runrail').innerText()));
    if (analyses.some(a => (a.links || []).some(l => l.run === linkedRun) && runsOf(a) > 1)) check(tag + ' run page: an analysis made from several flow cells says that it is of all of them together', /The analysis is of all of them together, not of this flow cell alone\./.test(t), t.slice(0, 260));
    check(tag + ' run page: nothing is wider than the window', (await sideways(page)) <= 0);
    await page.locator('#s-analyses').scrollIntoViewIfNeeded(); await settle(page, 300); await shot(page, 'run_analyses_' + tag);
  }
  await page.finish();
}
// ---- a one-file report of an analysis ----
if (reportFile) {
  const url = fileUrl(reportFile), page = await newPage({ allowed: u => u.startsWith(url) });
  await page.goto(url); await page.waitForSelector('.fc-analysis .fc-verdict__headline', { timeout: 20000 }); await page.waitForSelector('.fc-analysis .fc-runmain .fc-section', { timeout: 30000 }); await settle(page, 500);
  check('report: opens on its analysis with no tab of a library', await page.locator('.fc-tabs').count() === 0 && /Analysis report/.test(await page.locator('.fc-appbar__what').innerText()));
  check('report: draws the key numbers', await page.locator('.fc-analysis .fc-keynumbers .fc-tile').count() === 6);
  await shot(page, 'analysis_report_1366_light', true); await page.finish();
}
await browser.close();
report.finished = new Date().toISOString();
fs.writeFileSync(path.join(outDir, 'analysis_report.json'), JSON.stringify(report, null, 1));
const bad = report.checks.filter(c => !c.ok).length;
console.log('\n' + report.checks.length + ' checks, ' + bad + ' failed · errors: ' + log.errors.length + ' · requests to another origin: ' + log.requests.length);
for (const e of log.errors.slice(0, 12)) console.log('  ' + e.kind + ': ' + e.text.slice(0, 300));
process.exit(bad || log.errors.length || log.requests.length ? 1 : 0);
