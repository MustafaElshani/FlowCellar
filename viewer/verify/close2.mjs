// Regressions held as checks on any hub (the synthetic fixture and a real library):
//
//   A  a table that covers a sample on two flow cells is said so wherever its counts are shown: no section prints a
//      read total larger than the reads the flow cell sequenced without that sentence
//   B  a section whose files are linked but not read shows the notice and the scan command, and nothing that reads
//      as a measurement (no dash, no zero, no empty chart)
//   C  the Samples control offers every sample that a section of the page has numbers for
//   D  names of samples: a table heading holds the whole name, and the short names in read-outs stay different
//   E  the skip link moves the focus and leaves the address alone
//   F  a mistyped run address asks for no data file and logs nothing
//   G  the application bar keeps everything inside the window between 761 and 1,000 px
//   H  a row of the Analyses table opens the page of the pipeline whose verdict the row shows
//   I  tables of pages are as wide as the page (Projects, the samples of a run), like Runs and Analyses
//   J  a row of charts has columns of equal width: a chart left over in the last row is as wide as the others
//
//   node verify/close2.mjs --hub DIR --out DIR [--need shared,waiting,control,names,failrow]
//
// --need: the cases the hub must hold (the test of the fixture names them all); a case the hub does not hold is
// reported as "not on this hub" and is no failure otherwise. Exit code 0 when every check passed and nothing was logged.
import fs from 'node:fs';
import path from 'node:path';
import zlib from 'node:zlib';
import { launch, fileUrl, settle, open, sideways } from './browser.mjs';

const args = process.argv.slice(2), opt = n => { const i = args.indexOf('--' + n); return i >= 0 ? args[i + 1] : null; };
const hub = opt('hub'), outDir = path.resolve(opt('out') || 'verify-out'), need = (opt('need') || '').split(',').filter(Boolean), maxRuns = +(opt('runs') || 12);
if (!hub) { console.error('usage: node verify/close2.mjs --hub DIR --out DIR [--need shared,waiting,control,names,failrow]'); process.exit(2); }
const hubDir = path.resolve(hub), base = fileUrl(path.join(hubDir, 'index.html')), origin = fileUrl(hubDir);
fs.mkdirSync(outDir, { recursive: true });
const gun = b64 => JSON.parse(zlib.gunzipSync(Buffer.from(b64, 'base64')).toString('utf8'));
const ix = gun(/^__FC_INDEX__\("([A-Za-z0-9+/=]*)"\);\n$/.exec(fs.readFileSync(path.join(hubDir, 'index.js'), 'utf8'))[1]);
const bundleOf = id => { const f = path.join(hubDir, 'data', id + '.js'); if (!fs.existsSync(f)) return null; const m = /^__FCB_SIDECAR__\("[^"]*","([A-Za-z0-9+/=]*)","[0-9a-f]*"\);\n$/.exec(fs.readFileSync(f, 'utf8')); return m ? gun(m[1]) : null; };
const C = ix.runs.columns, N = ix.runs.n, col = (k, i) => { const c = C[k]; return c == null ? undefined : (!Array.isArray(c) && c.dict ? c.dict[c.codes[i]] : c[i]); };
const runs = Array.from({ length: N }, (_, i) => ({ id: col('id', i), fc: col('fc', i), read_data: col('read_data', i) }));
const analyses = ix.analyses || [];

const { browser, log, newPage, version } = await launch();
const report = { hub, browser: version, started: new Date().toISOString(), checks: [], cases: {}, shots: [] };
const check = (name, ok, detail) => { report.checks.push({ name, ok: !!ok, detail: detail == null ? '' : String(detail).slice(0, 400) }); console.log((ok ? 'ok   ' : 'FAIL ') + name + (detail != null && detail !== '' ? '  [' + String(detail).slice(0, 220) + ']' : '')); };
const seen = k => { report.cases[k] = (report.cases[k] || 0) + 1; };
const page0 = o => newPage({ width: 1366, height: 800, ...o, allowed: u => u.startsWith(origin) });
const txt = (page, sel) => page.evaluate(s => { const e = document.querySelector(s); return e ? e.innerText.replace(/\s+/g, ' ').trim() : null; }, sel);
// "130 M", "54,008,496", "9,000" at the start of a count in a sentence -> a number
const amount = t => { const m = /([\d][\d,.]*)\s*(M\b)?/.exec(t || ''); if (!m) return null; const v = parseFloat(m[1].replace(/,/g, '')); return m[2] ? v * 1e6 : v; };

// ---------------------------------------------------------------------------------------------------------
// A, B, C, D: the run pages that hold a case
// ---------------------------------------------------------------------------------------------------------
const SHARED = [['rna.biotype', 'rna_composition', /Of the ([\d][\d,.]*\s*M?) reads given to a transcript/], ['cdna.pychopper', 'pychopper', /Pychopper was given ([\d][\d,.]*\s*M?) reads/], ['cdna', 'cdna', /of the ([\d][\d,.]*\s*M?) reads pychopper classified/]];
const withData = runs.filter(r => r.read_data && r.read_data !== 'absent');
// the runs that hold a case, those with a shared table first, then those with pychopper tables; at most --runs of them
const hasNote = (b, id) => (((b.modules || {})[id] || {}).notes || []).some(n => n.code === 'covers_several_runs');
const weight = b => (['rna.biotype', 'cdna.pychopper', 'cdna'].some(id => hasNote(b, id)) ? 4 : 0) + ((((b.modules || {})['cdna.pychopper'] || {}).cells || []).length >= 2 ? 2 : 0) + ((b.modules || {}).aln && (b.modules.aln.cells || []).length ? 1 : 0);
const ordered = withData.map(r => ({ r, b: bundleOf(r.id) })).filter(q => q.b).map(q => ({ ...q, w: weight(q.b) })).sort((a, z) => z.w - a.w);
let looked = 0;
const page = await page0();
for (const { r, b } of ordered) {
  const mods = b.modules || {}, has = (id, code) => ((mods[id] || {}).notes || []).some(n => n.code === code);
  const shared = SHARED.filter(([id]) => has(id, 'covers_several_runs') && (mods[id].cells || []).length && !(id === 'cdna' && mods['cdna.pychopper']));
  const waiting = [['len', 'length_quality'], ['aln', 'alignment']].filter(([id]) => mods[id] && mods[id].state === 'partial' && mods[id].partial_reason === 'waiting' && !(mods[id].cells || []).length);
  const pych = (mods['cdna.pychopper'] || {}).cells || [];
  const interesting = shared.length || waiting.length || pych.length >= 2;
  if (!interesting || looked >= maxRuns) continue;
  looked++;
  const reads = ((mods.yield || {}).summary || {}).yield_reads;
  await open(page, base + '#/run/' + r.id); await page.waitForSelector('.fc-verdict__headline', { timeout: 30000 }); await settle(page, 700);
  // ---- A ----
  for (const [id, sec, re] of shared) {
    seen('shared');
    const said = await txt(page, '#s-' + sec + ' .fc-shared'), reading = await txt(page, '#s-' + sec + ' .fc-section__reading');
    check(r.fc + ' ' + sec + ': the table that covers two flow cells is said at the top of the section, as one analysis of merged reads', !!said && /^One (pychopper analysis|read table) for (two|\d+) flow cells: /.test(said) && /were merged before/.test(said) && /not this flow cell alone/.test(said), said || 'no notice');
    const m = re.exec(reading || ''), total = m ? amount(m[1]) : null;
    check(r.fc + ' ' + sec + ': the sentence gives its read total and says that it is of both flow cells', total != null && /merged from (both|all \d+) flow cells of the sample/.test(reading), reading);
    if (total != null && reads) check(r.fc + ' ' + sec + ': a read total above the reads of the flow cell (' + reads + ') never stands without that sentence', !(total > reads) || (/flow cells/.test(reading) && !!said), total + ' against ' + reads);
    const caps = await page.evaluate(s => [...document.querySelectorAll('#s-' + s + ' .fc-tile')].map(e => e.innerText.replace(/\s+/g, ' ')).filter(t => /(\d|M) of [\d.,]+( M)? (reads|pieces)|rows of the annotated/.test(t)), sec);
    check(r.fc + ' ' + sec + ': every tile that prints a count says that the flow cells were merged before the tool ran', caps.length > 0 && caps.every(t => /flow cells, merged before/.test(t)), caps.find(t => !/flow cells, merged before/.test(t)) || caps.length + ' tiles');
    await page.locator('#s-' + sec).scrollIntoViewIfNeeded(); await settle(page, 200);
    await page.locator('#s-' + sec).screenshot({ path: path.join(outDir, 'shared_' + r.fc + '_' + sec + '.png') }).catch(() => {}); report.shots.push('shared_' + r.fc + '_' + sec + '.png');
  }
  // ---- B ----
  for (const [id, sec] of waiting) {
    seen('waiting');
    const s = await page.evaluate(q => { const e = document.querySelector('#s-' + q); return e ? { text: e.innerText.replace(/\s+/g, ' '), tiles: e.querySelectorAll('.fc-tile').length, canvases: e.querySelectorAll('canvas').length, charts: e.querySelectorAll('.fc-chart').length } : null; }, sec);
    check(r.fc + ' ' + sec + ': files that wait are said, with the command that reads them', !!s && /wait(s)? for the next scan/.test(s.text) && /flowcellar scan -L \S+ --only /.test(s.text), s && s.text.slice(0, 200));
    check(r.fc + ' ' + sec + ': nothing that reads as a measurement while no file was read (no tile, no chart)', !!s && s.tiles === 0 && s.canvases === 0 && s.charts === 0, s && JSON.stringify({ tiles: s.tiles, canvases: s.canvases }));
    check(r.fc + ' ' + sec + ': no unfinished sentence, no zero printed as a count', !!s && !/of – or longer|These are reads\.| 0 · 0 b|No alignment numbers for this selection/.test(s.text), s && s.text.slice(0, 160));
  }
  // ---- C ----
  {
    const keys = new Set(); for (const [id, m] of Object.entries(mods)) { const si = (m.dims || []).indexOf('sample'); if (si >= 0) for (const c of m.cells || []) if (c.at[si] !== 'unclassified') keys.add(c.at[si]); }
    for (const c of (mods['funnel.samples'] || {}).cells || []) keys.add(c.at[0]);
    const named = (b.samples || []).filter(s => s.key !== 'unclassified' && keys.has(s.key));
    const opts = await page.evaluate(() => { const sel = [...document.querySelectorAll('.fc-runcontrols select')].find(s => [...s.options].some(o => /added up|All reads/.test(o.textContent))); return sel ? [...sel.options].map(o => [o.value, o.textContent]) : null; });
    if (opts && named.length) {
      seen('control');
      const offered = opts.filter(o => /^s\d+$/.test(o[0]));
      check(r.fc + ' Samples control: offers every sample a section has numbers for', offered.length === named.length && new RegExp('All ' + named.length + ' added up').test(opts[0][1]), offered.length + ' offered, ' + named.length + ' with numbers; first option: ' + opts[0][1]);
      check(r.fc + ' Samples control: the names offered are all different', new Set(offered.map(o => o[1])).size === offered.length, offered.map(o => o[1]).join(' | '));
      // the rows of the Samples table: a sample with numbers can be chosen by its row
      const rowsOn = await page.evaluate(() => [...document.querySelectorAll('#s-samples tbody tr')].filter(tr => tr.querySelector('button.fc-samples__name')).length);
      check(r.fc + ' Samples table: every sample with numbers can be chosen by its row', rowsOn === named.length, rowsOn + ' of ' + named.length);
      const first = named[0];
      await page.evaluate(k => { const sel = [...document.querySelectorAll('.fc-runcontrols select')].find(s => [...s.options].some(o => /added up/.test(o.textContent))); sel.value = k; sel.dispatchEvent(new Event('change', { bubbles: true })); }, first.key); await settle(page, 500);
      check(r.fc + ' Samples control: a choice is in the address', new RegExp('[?&]s=' + first.key + '\\b').test(page.url()), page.url().slice(-60));
      const none = await page.evaluate(() => [...document.querySelectorAll('.fc-runmain .fc-section[id]')].filter(s => s.id && !s.querySelector('.fc-chart, .fc-track, .fc-tablewrap, .fc-empty, .fc-notice, .fc-readcard')).map(s => s.id));
      check(r.fc + ' one sample chosen: no section is left as a bare title', none.length === 0, none.join(', '));
      await page.evaluate(() => { const sel = [...document.querySelectorAll('.fc-runcontrols select')].find(s => [...s.options].some(o => /added up/.test(o.textContent))); sel.value = ''; sel.dispatchEvent(new Event('change', { bubbles: true })); }); await settle(page, 400);
    }
  }
  // ---- D ----
  if (pych.length >= 2 && pych.length <= 6) {
    seen('names');
    const names = (b.samples || []).filter(s => pych.some(c => c.at[0] === s.key)).map(s => s.name);
    await page.locator('#p-pych-config').scrollIntoViewIfNeeded();
    await page.locator('#p-pych-config button:has-text("Show as table")').click(); await settle(page, 300);
    const heads = await page.evaluate(() => [...document.querySelectorAll('#p-pych-config table thead th')].map(e => e.textContent.trim()).slice(3));
    check(r.fc + ' a table heading holds the whole name of each sample', heads.length === names.length && heads.every((h, i) => h === names[i]), heads.join(' | '));
    // the read-out of every sample row of a per-sample chart: the names shown are all different
    await page.setViewportSize({ width: 1024, height: 800 }); await settle(page, 500);
    await page.locator('#p-pych-class').scrollIntoViewIfNeeded(); await settle(page, 200);
    const box = await page.evaluate(() => { const e = document.querySelector('#p-pych-class .fc-plot'), r = e.getBoundingClientRect(); return [r.left, r.top, r.width, r.height]; });
    const before = await page.evaluate(() => [...document.querySelectorAll('.fc-chart, .fc-section')].map(e => { const q = e.getBoundingClientRect(); return [Math.round(q.left), Math.round(q.top + scrollY), Math.round(q.width), Math.round(q.height)]; }));
    const rowH = 22, first = names.length > 1 ? 1 : 0, got = [];
    for (let k = 0; k < names.length; k++) { await page.mouse.move(box[0] + box[2] * 0.7, box[1] + 4 + (first + k) * rowH + rowH / 2); await settle(page, 80); got.push(await page.evaluate(() => { const b = document.querySelector('#p-pych-class .fc-ro b'); return b ? b.textContent : ''; })); }
    const after = await page.evaluate(() => [...document.querySelectorAll('.fc-chart, .fc-section')].map(e => { const q = e.getBoundingClientRect(); return [Math.round(q.left), Math.round(q.top + scrollY), Math.round(q.width), Math.round(q.height)]; }));
    check(r.fc + ' read-outs of the per-sample bars name ' + names.length + ' different samples', new Set(got).size === names.length && got.every(Boolean), got.join(' | '));
    check(r.fc + ' the pointer on the per-sample bars moves nothing', JSON.stringify(before) === JSON.stringify(after), '');
    await page.mouse.move(2, 2); await page.setViewportSize({ width: 1366, height: 800 }); await settle(page, 300);
  }
}
for (const k of ['shared', 'waiting', 'control', 'names']) if (!report.cases[k]) check('case "' + k + '": ' + (need.includes(k) ? 'this hub must hold one' : 'not on this hub'), !need.includes(k), '');

// ---------------------------------------------------------------------------------------------------------
// E: the skip link
// ---------------------------------------------------------------------------------------------------------
{
  const kb = await page0();                       // a page that was never clicked: the first Tab goes to its first stop
  await open(kb, base + '#/home'); await settle(kb, 300);
  await kb.keyboard.press('Tab');
  const on = await kb.evaluate(() => { const a = document.activeElement; return a ? a.className + '|' + a.textContent : ''; });
  await kb.keyboard.press('Enter'); await settle(kb, 300);
  const st = await kb.evaluate(() => ({ hash: location.hash, active: document.activeElement && document.activeElement.id, text: document.body.innerText.slice(0, 4000) }));
  await kb.finish();
  check('skip link: the first Tab stop is "Skip to the content"', /fc-skip/.test(on) && /Skip to the content/.test(on), on);
  check('skip link: Enter moves the focus to the content and leaves the address alone', st.hash === '#/home' && st.active === 'main' && !/not a page of the library/.test(st.text), JSON.stringify({ hash: st.hash, active: st.active }));
}
// ---------------------------------------------------------------------------------------------------------
// F: a mistyped run address
// ---------------------------------------------------------------------------------------------------------
{
  const asked = []; const onReq = r => { if (/XXX99999_deadbeef/.test(r.url())) asked.push(r.url()); }; page.on('request', onReq);
  const errs = log.errors.length;
  await open(page, base + '#/run/XXX99999_deadbeef'); await settle(page, 800);
  const t = await txt(page, '#main');
  check('a mistyped run address says that the run is not in the library', /not in the library/.test(t || ''), (t || '').slice(0, 120));
  check('a mistyped run address asks for no data file and logs nothing', asked.length === 0 && log.errors.length === errs, asked.concat(log.errors.slice(errs).map(e => e.text)).join(' ; '));
  page.off('request', onReq);
}
// ---------------------------------------------------------------------------------------------------------
// G: the application bar in a narrow window
// ---------------------------------------------------------------------------------------------------------
for (const w of [1000, 940, 900, 820, 768]) {
  await page.setViewportSize({ width: w, height: 600 }); await open(page, base + '#/runs'); await settle(page, 300);
  const kids = await page.evaluate(() => [...document.querySelectorAll('.fc-appbar__in > *')].map(e => { const r = e.getBoundingClientRect(); return { cls: String(e.className), l: Math.round(r.left), r: Math.round(r.right), shown: r.width > 0 }; }));
  const out = kids.filter(k => k.shown && (k.r > w || k.l < 0)), theme = kids.filter(k => /fc-theme/.test(k.cls) && k.shown), compare = kids.find(k => /fc-appbar__compare/.test(k.cls));
  check('application bar at ' + w + ' px: nothing lies outside the window', out.length === 0 && (await sideways(page)) <= 1, JSON.stringify(out));
  check('application bar at ' + w + ' px: a theme control and Compare are there', theme.length === 1 && !!compare && compare.shown, JSON.stringify({ theme: theme.length, compare: compare && compare.shown }));
}
await page.setViewportSize({ width: 1366, height: 800 });
// ---------------------------------------------------------------------------------------------------------
// H: a row of the Analyses table opens the pipeline that decides its verdict
// ---------------------------------------------------------------------------------------------------------
if (analyses.length) {
  const RANK = { fail: 3, warn: 2, pass: 1 }, byId = new Map(analyses.map(a => [a.id, a]));
  const pairs = analyses.filter(a => a.pipeline === 'spaceranger' && (a.siblings || []).some(s => byId.has(s))).map(a => [a, byId.get(a.siblings.find(s => byId.has(s)))]);
  const differ = pairs.filter(([a, s]) => (RANK[s.verdict] || 0) > (RANK[a.verdict] || 0));
  const WORD = { fail: 'Fail', warn: 'Warning', pass: 'Pass' };
  for (const [a, s] of differ.slice(0, 3)) {
    seen('failrow');
    await open(page, base + '#/analyses'); await page.waitForSelector('.fc-analyses tbody tr'); await settle(page, 300);
    const row = page.locator('.fc-analyses tbody tr', { has: page.locator('a[href*="' + a.id + '"]') }).first();
    const badge = (await row.locator('.fc-badge').first().innerText()).trim();
    await row.locator('a.fc-projects__name').click(); await page.waitForSelector('.fc-anverdict .fc-badge', { timeout: 30000 }); await settle(page, 500);
    const top = (await page.locator('.fc-anverdict .fc-badge').first().innerText()).trim(), both = await txt(page, '.fc-ancontrols__verdicts');
    check('Analyses: the row of ' + (a.label || a.id) + ' says ' + badge + ' and opens a page that says ' + badge, badge === WORD[s.verdict] && top === badge && /[?&]p=/.test(page.url()), 'row ' + badge + ', page ' + top + ', ' + page.url().slice(-50));
    check('Analyses: the page of a pair shows the verdict of each pipeline beside the switch', !!both && both.includes(WORD[a.verdict]) && both.includes(WORD[s.verdict]), both);
  }
  if (!differ.length) check('case "failrow": ' + (need.includes('failrow') ? 'this hub must hold a pair whose two verdicts differ' : 'not on this hub'), !need.includes('failrow'), pairs.length + ' pairs');
  // the sort of the Analyses table is in the address
  await open(page, base + '#/analyses'); await page.waitForSelector('.fc-analyses thead button'); await page.locator('.fc-analyses thead button', { hasText: 'Verdict' }).click(); await settle(page, 300);
  check('Analyses: the sort is kept in the address', /[?&]sort=-?verdict/.test(page.url()), page.url().slice(-40));
  // an empty search of Runs says where the analyses are
  await open(page, base + '#/runs?q=zzzz_no_such_run'); await settle(page, 400);
  const empty = await txt(page, '.fc-hub__empty');
  check('Runs: an empty search says that analyses are listed under Analyses', /listed under Analyses/.test(empty || ''), (empty || '').slice(0, 160));
  await open(page, base + '#/home'); await settle(page, 300);
  const look = await txt(page, '.fc-home__anlook');
  check('Home: the Analyses block says how many need a look', !!look && /(failed|with a warning|No analysis has a warning)/.test(look), look);
}
// ---------------------------------------------------------------------------------------------------------
// I, J: wide windows
// ---------------------------------------------------------------------------------------------------------
for (const w of [1920, 2560, 3440]) {
  const wide = await page0({ width: w, height: 1300 });
  const edge = async sel => wide.evaluate(s => { const pg = document.querySelector('.fc-page'), e = document.querySelector(s); if (!pg || !e) return null; const cs = getComputedStyle(pg), p = pg.getBoundingClientRect(), r = e.getBoundingClientRect(); return { right: Math.round(r.right), page: Math.round(p.right - parseFloat(cs.paddingRight)) }; }, sel);
  await open(wide, base + '#/projects'); await wide.waitForSelector('.fc-projectspage table'); await settle(wide, 300);
  const pr = await edge('.fc-projectspage .fc-home__block');
  check('Projects at ' + w + ' px: the table reaches the right edge of the page, like Runs', pr && Math.abs(pr.right - pr.page) <= 2, JSON.stringify(pr));
  const pick = withData.find(r => { const b = bundleOf(r.id); return b && (b.samples || []).length > 1; });
  if (pick) {
    await open(wide, base + '#/run/' + pick.id); await wide.waitForSelector('.fc-verdict__headline', { timeout: 30000 }); await settle(wide, 800);
    const st = await wide.evaluate(() => { const main = document.querySelector('.fc-runmain'), t = document.querySelector('#s-samples > .fc-tablewrap'); if (!main || !t) return null; const sec = t.parentElement, cs = getComputedStyle(sec), s = sec.getBoundingClientRect(); return { right: Math.round(t.getBoundingClientRect().right), inner: Math.round(s.right - parseFloat(cs.paddingRight)) }; });
    if (st) check('run at ' + w + ' px: the Samples table is as wide as its section', Math.abs(st.right - st.inner) <= 2, JSON.stringify(st));
    // every row of charts: columns of one width, and a row that is not full starts at the left
    const grids = await wide.evaluate(() => [...document.querySelectorAll('.fc-grid--n')].map(g => { const kids = [...g.children].map(e => e.getBoundingClientRect()).filter(r => r.width > 0); const ws = kids.map(r => Math.round(r.width)), tops = [...new Set(kids.map(r => Math.round(r.top)))]; return { n: kids.length, min: Math.min(...ws), max: Math.max(...ws), rows: tops.length, perRow: tops.map(t => kids.filter(r => Math.round(r.top) === t).length), left: kids.length ? Math.round(kids[0].left - g.getBoundingClientRect().left) : 0 }; }).filter(g => g.n > 0));
    const uneven = grids.filter(g => g.max - g.min > 2 || g.left !== 0), lonely = grids.filter(g => g.n >= 4 && g.perRow.length > 1 && g.perRow[g.perRow.length - 1] === 1 && g.perRow[0] > 2);
    check('run at ' + w + ' px: the charts of a row have one width (' + grids.length + ' rows of charts)', grids.length > 0 && uneven.length === 0, JSON.stringify(uneven.slice(0, 3)));
    check('run at ' + w + ' px: no row of three or more leaves a single chart over', lonely.length === 0, JSON.stringify(lonely.slice(0, 3)));
    check('run at ' + w + ' px: nothing is wider than the window', (await sideways(wide)) <= 1, await sideways(wide));
    if (w === 1920) { await wide.screenshot({ path: path.join(outDir, 'run_1920.png'), fullPage: true }); report.shots.push('run_1920.png'); }
  }
  await wide.finish();
}
await page.finish();
await browser.close();
report.finished = new Date().toISOString(); report.errors = log.errors; report.requests = log.requests;
const failed = report.checks.filter(c => !c.ok);
report.summary = { checks: report.checks.length, failed: failed.length, errors: log.errors.length, requests_to_another_origin: log.requests.length, cases: report.cases };
fs.writeFileSync(path.join(outDir, 'close2_report.json'), JSON.stringify(report, null, 1));
console.log('\n' + report.checks.length + ' checks, ' + failed.length + ' failed · page and console errors and policy violations: ' + log.errors.length + ' · requests to another origin: ' + log.requests.length + ' · cases: ' + JSON.stringify(report.cases));
for (const e of log.errors.slice(0, 15)) console.log('  ' + e.kind + ' @ ' + e.url.slice(-60) + ': ' + e.text.slice(0, 300));
for (const r of log.requests.slice(0, 10)) console.log('  request: ' + r);
process.exit(failed.length || log.errors.length || log.requests.length ? 1 : 0);
