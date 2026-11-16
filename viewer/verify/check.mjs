// Use a hub of the viewer like a person would, in headless Firefox: every screen in both themes at three widths,
// with screenshots and a report of what was checked. It works on any hub (it reads the hub's own index.js to know
// what to expect), so the same script checks the synthetic fixture and a real library.
//
//   node verify/check.mjs --hub DIR_OR_URL --out DIR [--fixture] [--report FILE.html] [--quick] [--label NAME]
//
//   --hub      a hub folder (opened from file://) or the address of a server that serves one
//   --fixture  also run the checks that only hold for the synthetic library of tests/viewer/viewer_synth.py
//              (names that try to be markup, the Reload control, a missing data file)
//   --report   a one-run report file to check as well
//   --quick    one theme and one width (for a fast look while developing)
// Exit code 0 when every check passed and nothing was logged, 1 otherwise. Writes OUT/report.json.
import fs from 'node:fs';
import path from 'node:path';
import zlib from 'node:zlib';
import { launch, fileUrl, settle, open, setTheme, sideways } from './browser.mjs';

const args = process.argv.slice(2), opt = n => { const i = args.indexOf('--' + n); return i >= 0 ? args[i + 1] : null; }, flag = n => args.includes('--' + n);
const hub = opt('hub'), outDir = path.resolve(opt('out') || 'verify-out'), fixture = flag('fixture'), quick = flag('quick'), label = opt('label') || '';
if (!hub) { console.error('usage: node verify/check.mjs --hub DIR_OR_URL --out DIR [--fixture] [--report FILE] [--quick]'); process.exit(2); }
const http = /^https?:/.test(hub), hubDir = http ? null : path.resolve(hub);
const base = http ? hub.replace(/\/?$/, '/') : fileUrl(path.join(hubDir, 'index.html')), origin = http ? new URL(base).origin : fileUrl(hubDir);
fs.mkdirSync(outDir, { recursive: true });

// ---- what the hub holds, read the way the viewer reads it ----
async function readIndex() {
  const text = http ? await (await fetch(base + 'index.js')).text() : fs.readFileSync(path.join(hubDir, 'index.js'), 'utf8');
  const m = /^__FC_INDEX__\("([A-Za-z0-9+/=]*)"\);\n$/.exec(text); if (!m) throw new Error('index.js is not one call of __FC_INDEX__');
  return JSON.parse(zlib.gunzipSync(Buffer.from(m[1], 'base64')).toString('utf8'));
}
// the read-data sections of one run that its own bundle marks partial
async function partialSections(id) {
  const text = http ? await (await fetch(base + 'data/' + id + '.js')).text() : fs.readFileSync(path.join(hubDir, 'data', id + '.js'), 'utf8');
  const m = /^__FCB_SIDECAR__\("[^"]*","([A-Za-z0-9+/=]*)","[0-9a-f]*"\);\n$/.exec(text); if (!m) throw new Error('data/' + id + '.js is not one call of __FCB_SIDECAR__');
  const b = JSON.parse(zlib.gunzipSync(Buffer.from(m[1], 'base64')).toString('utf8'));
  return Object.entries(b.modules || {}).filter(([, s]) => s.tier === 'read_level' && s.state === 'partial').map(([k]) => k);
}
const ix = await readIndex(), C = ix.runs.columns, N = ix.runs.n;
const col = (k, i) => { const c = C[k]; return c == null ? undefined : (!Array.isArray(c) && c.dict ? c.dict[c.codes[i]] : c[i]); };
const rows = Array.from({ length: N }, (_, i) => { const r = {}; for (const k of Object.keys(C)) r[k] = col(k, i); return r; });
const newest = rows.slice().sort((a, b) => Date.parse(b.started) - Date.parse(a.started));
const pick = f => (newest.find(f) || {}).id || null;
const kinds = {
  complete: pick(r => r.detail_t1 === 'full' && r.read_data === 'complete'), partial: pick(r => r.read_data === 'partial'), quickRun: pick(r => r.detail_t1 === 'quick' && r.run_state !== 'running' && r.provisional),
  // quick detail for good: the read table is cut off (or could not be read), so nothing waits
  cut: pick(r => r.detail_t1 === 'quick' && r.run_state !== 'running' && !r.provisional && (r.modifiers || []).some(m => m === 'reads_partly_counted' || m === 'reads_not_counted')),
  full: pick(r => r.detail_t1 === 'full' && r.read_data === 'absent'), stub: pick(r => r.detail_t1 === 'none' || r.detail_t1 == null), running: pick(r => r.run_state === 'running'),
  short: pick(r => r.cause === 'stopped_early' || r.cause === 'output_not_judged'), many: pick(r => (r.t2_samples || 0) >= 13),
};
const reused = (newest.find(r => r.use_of > 1) || {}).fc || null;
const fullRuns = newest.filter(r => r.detail_t1 === 'full').map(r => r.id);
const attention = newest.filter(r => r.verdict === 'warn' || r.verdict === 'fail');
const projects = ix.projects || [], firstProject = projects[0] || null;

const { browser, log, newPage, version } = await launch();
const report = { hub, label, browser: version, started: new Date().toISOString(), runs: N, kinds, checks: [], shots: [], errors: log.errors, requests: log.requests, notices: 0 };
const check = (name, ok, detail) => { report.checks.push({ name, ok: !!ok, detail: detail == null ? '' : String(detail).slice(0, 300) }); console.log((ok ? 'ok   ' : 'FAIL ') + name + (detail != null && detail !== '' ? '  [' + String(detail).slice(0, 160) + ']' : '')); };
const page0 = opts => newPage({ ...opts, allowed: u => u.startsWith(origin) });
const shot = async (page, name, full) => { if (full) { await page.evaluate(() => window.scrollTo(0, 0)); await settle(page, 200); } await page.screenshot({ path: path.join(outDir, name + '.png'), fullPage: !!full }); report.shots.push(name + '.png'); };
const text = (page, sel) => page.locator(sel).first().innerText().then(t => t.replace(/\s+/g, ' ').trim());
const count = (page, sel) => page.locator(sel).count();
const hubRows = page => count(page, '.fc-table--hub tbody tr:not(.fc-table__group)');
const pwned = page => page.evaluate(() => window.__pwned);
// no header cell of the Runs table reaches outside its frame, and the frame does not scroll sideways
const halfColumn = page => page.evaluate(() => { const w = document.querySelector('.fc-hub__table'); if (!w) return 'no table'; const r = w.getBoundingClientRect(); let worst = 0; for (const th of w.querySelectorAll('thead th')) worst = Math.max(worst, th.getBoundingClientRect().right - r.right); return Math.max(worst, w.scrollWidth - w.clientWidth); });
// every control a pointer can use can be reached with the Tab key
const unreachable = page => page.evaluate(() => { const bad = []; for (const el of document.querySelectorAll('a[href], button, input, select, textarea, [role="button"], [role="radio"], [role="slider"], [role="checkbox"], [onclick]')) {
  const cs = getComputedStyle(el), r = el.getBoundingClientRect(); if (cs.display === 'none' || cs.visibility === 'hidden' || (!r.width && !r.height)) continue;
  if (el.disabled || el.getAttribute('aria-hidden') === 'true') continue; if (el.tabIndex < 0) bad.push(el.tagName.toLowerCase() + '.' + el.className + ' "' + (el.textContent || '').trim().slice(0, 30) + '"'); } return bad; });

const themes = quick ? ['Light'] : ['Light', 'Dark'], widths = quick ? [1366] : [1366, 1280, 820];
const H = w => (w === 1366 ? 768 : w === 1280 ? 720 : 1180);

// =====================================================================================================
// 1. Home: the default address, drawn from the index, numbers equal to the index
// =====================================================================================================
{
  const page = await page0({ width: 1366, height: 768 });
  await open(page, base);
  check('the default address is Home', await count(page, '.fc-home') === 1 && /#?$|#\/home/.test(new URL(page.url()).hash), page.url().slice(-30));
  const lib = ix.library || {}, c = lib.counts || {}, line = await text(page, '.fc-home__line');
  check('Home header: runs, flow cells and projects equal the index', line.includes(N + ' run') && line.includes(c.flow_cells + ' flow cell') && (!c.projects || line.includes(c.projects + ' project')), line);
  check('Home header names the library', (await text(page, '.fc-home__name')) === (lib.name || 'FlowCellar library'), await text(page, '.fc-home__name'));
  const tb = lib.estimated_bases / 1e12, shown = /([\d.,]+)\s(Tb|Gb|Mb) estimated/.exec(line);
  check('Home header: estimated output equals the index', shown && Math.abs(parseFloat(shown[1].replace(/,/g, '')) * (shown[2] === 'Tb' ? 1 : shown[2] === 'Gb' ? 1e-3 : 1e-6) - tb) <= Math.max(0.06 * tb, 1e-6), shown && shown[0]);
  // needs attention
  const att = await page.locator('#h-attention').locator('xpath=ancestor::section').locator('.fc-homerow').evaluateAll(els => els.map(e => ({ id: e.querySelector('a.fc-id').getAttribute('href'), text: e.querySelector('.fc-homerow__text').textContent, badge: e.querySelector('.fc-badge').textContent })));
  check('Needs attention: the newest ' + Math.min(6, attention.length) + ' runs with a warning or a fail, newest first', att.length === Math.min(6, attention.length) && att.every((a, k) => a.id.endsWith(encodeURIComponent(attention[k].id))), att.map(a => a.id.slice(-18)).join(' '));
  check('Needs attention: each row shows the stored headline sentence', att.every((a, k) => a.text === attention[k].headline), att[0] && att[0].text);
  if (attention.length) { const all = await text(page, '#h-attention ~ .fc-home__all'); check('Needs attention: the link counts all of them', all.includes('All ' + attention.length + ' '), all); }
  else check('Needs attention says that nothing needs attention', /Nothing needs attention/.test(await text(page, '.fc-home__calm')));
  // what still waits is the index's own progress (one count): the runs whose reads are not counted yet and the files not read yet
  const pgr = ix.progress || {}, counting = (pgr.tier1 || {}).runs_waiting || 0, filesWaiting = (pgr.read_data || {}).files_waiting || 0;
  if (counting || filesWaiting) {
    const aside = (await page.locator('.fc-home__progress').allInnerTexts()).join(' ').replace(/\s+/g, ' ');
    if (counting) check('runs still being counted are a line of their own, not warnings', aside.includes(counting + ' run'), aside.slice(0, 160));
    check('Home says whether a scan is running, and how to go on when none is', pgr.state === 'scanning' ? /A scan is running/.test(aside) : /Waiting for the next scan/.test(aside) && /flowcellar scan -L /.test(aside), aside.slice(0, 200));
    if (pgr.message) check('Home shows the progress sentence of the index', aside.includes(pgr.message), pgr.message);
  } else check('nothing waits: Home shows no progress block', await count(page, '.fc-home__progress') === 0);
  const cutRuns = rows.filter(r => (r.modifiers || []).some(m => m === 'reads_partly_counted' || m === 'reads_not_counted')).length;
  if (cutRuns) check('a run whose read table is cut off is not called "counting"', (await page.locator('.fc-home__aside').allInnerTexts()).some(t => /cut off or cannot be read/.test(t) && t.includes(cutRuns + ' run')), cutRuns);
  // latest runs
  const latest = await page.locator('#h-latest').locator('xpath=ancestor::section').locator('.fc-homerow a.fc-id').evaluateAll(els => els.map(e => e.getAttribute('href')));
  check('Latest runs: the six newest, newest first', latest.length === Math.min(6, N) && latest.every((h, k) => h.endsWith(encodeURIComponent(newest[k].id))), latest.length);
  check('Latest runs: each has a badge, a sparkline or none, and an output', await count(page, '.fc-homerow--latest .fc-badge') === latest.length && await count(page, '.fc-homerow--latest .fc-homerow__out') === latest.length);
  // projects
  const prow = await page.locator('.fc-home .fc-projects:not(.fc-analyses) tbody tr').evaluateAll(trs => trs.map(tr => [...tr.children].map(td => td.innerText.replace(/\s+/g, ' ').trim())));
  check('Projects: one row per project (at most eight on Home)', prow.length === Math.min(8, projects.length), prow.length + ' of ' + projects.length);
  if (firstProject) check('Projects: the first row equals the index (runs, flow cells, where it is)', prow[0][1] === String(firstProject.runs) && prow[0][2] === String(firstProject.flow_cells) && firstProject.roots.every(r => prow[0][7].includes(r)), prow[0].join(' | '));
  const two = projects.find(p => (p.roots || []).length > 1);
  if (two) { const i = projects.indexOf(two); if (i < prow.length) check('a project found under two folders shows both', two.roots.every(r => prow[i][7].includes(r)), prow[i][7]); }
  // roots and commands
  const rrow = await page.locator('.fc-roots tbody tr').evaluateAll(trs => trs.map(tr => tr.innerText.replace(/\s+/g, ' ').trim()));
  check('Where this library comes from: one row per scanned folder, with its path and run count', rrow.length === (ix.roots || []).length && (ix.roots || []).every((r, k) => rrow[k].includes(r.path) && rrow[k].includes(String(r.runs))), rrow.length);
  const off = (ix.roots || []).find(r => r.state === 'offline'); if (off) check('an offline folder is said to be offline', /Offline/.test(rrow[ix.roots.indexOf(off)]), rrow[ix.roots.indexOf(off)]);
  const cmd = await text(page, '.fc-howto .fc-cmd');
  check('the commands to add more carry this library\'s own path', cmd.includes('flowcellar scan /another/folder -L') && cmd.includes('flowcellar run /one/run_folder -L') && (!lib.path || cmd.split(lib.path).length === 3), cmd);
  check('one sentence says that nothing is copied', /Nothing is copied/.test(await text(page, '.fc-howto')));
  // activity
  check('Activity: two charts are drawn', await count(page, '.fc-home .fc-chart canvas.fck-base') === 2);
  await page.locator('.fc-home .fc-chart button:has-text("Show as table")').first().click(); await settle(page, 200);
  const months = await page.locator('.fc-home .fc-charttable tbody tr').evaluateAll(trs => trs.map(tr => +tr.children[1].textContent.replace(/,/g, '')));
  check('Activity: "Show as table" lists every month and the runs add up to the library', months.length === ix.activity.months.length && months.reduce((a, v) => a + v, 0) === N, months.length + ' months, ' + months.reduce((a, v) => a + v, 0) + ' runs');
  const bottom = await page.evaluate(() => document.querySelector('.fc-home__top').getBoundingClientRect().bottom + window.scrollY);
  check('the first three blocks need no scrolling at 1,366 x 768', bottom <= 768, Math.round(bottom) + ' px');
  // tabs, logo, search, Back
  const tabs = await page.locator('.fc-appbar .fc-tab').allInnerTexts();
  // a hub that holds single-cell or spatial analyses has a fourth tab for them; one without has three
  const nAn = (ix.analyses || []).length;
  check('tabs Home, Runs, Projects' + (nAn ? ', Analyses' : ''), tabs.length === (nAn ? 4 : 3) && /^Home/.test(tabs[0]) && /^Runs/.test(tabs[1]) && /^Projects/.test(tabs[2]) && (!nAn || /^Analyses/.test(tabs[3])), tabs.join(' | ').replace(/\s+/g, ' '));
  await page.click('.fc-appbar .fc-tab:has-text("Runs")'); await page.waitForSelector('.fc-table--hub'); await settle(page, 200);
  check('the Runs tab lists every run', await hubRows(page) === N, await hubRows(page));
  await page.click('.fc-appbar__logo'); await page.waitForSelector('.fc-home'); check('the logo goes Home', new URL(page.url()).hash === '#/home');
  await page.goBack(); await page.waitForSelector('.fc-table--hub'); check('Back returns to Runs', new URL(page.url()).hash.startsWith('#/runs'));
  await page.goForward(); await page.waitForSelector('.fc-home');
  await page.fill('.fc-home__search input', 'status:warn'); await page.keyboard.press('Enter'); await page.waitForSelector('.fc-table--hub'); await settle(page, 200);
  check('the search box of Home goes to Runs with the query', /q=status:warn/.test(page.url()) && await hubRows(page) === rows.filter(r => r.verdict === 'warn').length, (await hubRows(page)) + ' rows · ' + new URL(page.url()).hash);
  await page.goBack(); await page.waitForSelector('.fc-home'); check('Back returns Home from a search', await count(page, '.fc-home') === 1);
  const bad = await unreachable(page); check('Home: every control can be reached with the keyboard', bad.length === 0, bad.join('; '));
  if (fixture) check('names that try to be markup did not run on Home', await pwned(page) === undefined, await pwned(page));
  await page.finish();
}

// =====================================================================================================
// 2. Runs: table, search, sort, filters, peek, selection, compare
// =====================================================================================================
{
  const page = await page0({ width: 1366, height: 768 });
  await open(page, base + '#/runs');
  const ids = await page.locator('.fc-table--hub tbody tr a.fc-id').evaluateAll(els => els.map(e => e.getAttribute('href')));
  check('Runs: newest first', ids.every((h, k) => h.endsWith(encodeURIComponent(newest[k].id))) || ids.length === N, ids.length);
  check('Runs at 1,366 px: no column is half visible', await halfColumn(page) <= 1, await halfColumn(page));
  const nWarn = rows.filter(r => r.verdict === 'warn').length;
  for (const [q, expect] of [['status:warn', nWarn], ['no:reads', rows.filter(r => r.read_data === 'absent').length], ['wibble:1', N], ['is:reused', rows.filter(r => r.use > 1).length], ['fc:' + newest[0].fc, rows.filter(r => r.fc.toLowerCase().startsWith(newest[0].fc.toLowerCase())).length]]) {
    await page.fill('.fc-appbar input[type=search]', q); await settle(page, 250);
    check('search "' + q.replace(/fc:.*/, 'fc:<newest flow cell>') + '"', await hubRows(page) === expect, (await hubRows(page)) + ' rows, expected ' + expect);
  }
  check('a term that is not understood is shown struck through and never applied', await count(page, '.fc-filter.is-unknown') === 0);
  await page.fill('.fc-appbar input[type=search]', 'wibble:1'); await settle(page, 200); check('an unknown field is marked', await count(page, '.fc-filter.is-unknown') === 1);
  await page.fill('.fc-appbar input[type=search]', ''); await settle(page, 200);
  await page.click('.fc-table--hub th button:has-text("Pores")'); await settle(page, 200);
  const top = rows.filter(r => r.pores != null).sort((a, b) => b.pores - a.pores)[0];
  if (top) check('sorting by pores puts the largest first', (await page.locator('.fc-table--hub tbody tr a.fc-id').first().getAttribute('href')).endsWith(encodeURIComponent(top.id)) || rows.filter(r => r.pores === top.pores).length > 1, top.pores);
  check('the sorted header says so (aria-sort)', await count(page, '.fc-table--hub th[aria-sort="descending"]') === 1);
  await page.click('.fc-table--hub th button:has-text("Started")'); await settle(page, 200);
  await page.click('.fc-facets__btn'); await settle(page, 200);
  if (nWarn) { await page.locator('.fc-facet:has-text("Run health") label:has-text("Warning") input').check(); await settle(page, 250); check('filter Run health = Warning', await hubRows(page) === rows.filter(r => r.health === 'warn').length, await hubRows(page)); await page.click('button:has-text("Clear all")'); }
  await page.click('.fc-facets__head button'); await settle(page, 200);
  // peek from a click outside the ID; the pore-scan chart arrives from the run's own file
  const target = fullRuns[0] || newest[0].id;
  await page.locator('.fc-table--hub tbody tr', { has: page.locator('a[href$="' + encodeURIComponent(target) + '"]') }).locator('td').nth(3).click(); await settle(page, 900);
  check('a click on a row opens the peek drawer', await count(page, '.fc-peek') === 1 && /peek=/.test(page.url()));
  check('the peek drawer loads the run\'s own data on demand (pore-scan chart)', await count(page, '.fc-peek .fc-chart canvas.fck-base') >= 1, await count(page, '.fc-peek .fc-chart canvas.fck-base'));
  await shot(page, 'runs_peek_1366_light');
  await page.click('.fc-peek .fc-close'); await settle(page, 200);
  // keyboard: Down moves to the next ID, Enter opens that run
  await page.locator('.fc-table--hub tbody tr a.fc-id').first().focus(); await page.keyboard.press('ArrowDown'); const f = await page.evaluate(() => document.activeElement.textContent);
  await page.keyboard.press('Enter'); await page.waitForSelector('.fc-runhead__id'); check('keyboard: Down moves to the next flow cell ID, Enter opens that run', (await text(page, '.fc-runhead__id')) === f, f);
  await page.goBack(); await page.waitForSelector('.fc-table--hub'); await settle(page, 200);
  // selection and compare
  if (fullRuns.length >= 2) {
    for (const id of fullRuns.slice(0, 2)) await page.locator('.fc-table--hub tbody tr', { has: page.locator('a[href$="' + encodeURIComponent(id) + '"]') }).locator('input[type=checkbox]').check();
    await settle(page, 200); await page.click('.fc-tray a:has-text("Compare")'); await page.waitForSelector('.fc-compare', { timeout: 20000 }); await settle(page, 900);
    const note = await text(page, '.fc-compare .fc-section__note'), m = /\(B: ([\d.]+) × ([\d.]+) × ([\d.]+) = ([\d.]+)\)/.exec(note);
    check('compare from the tray: the three factors multiply to the output ratio', m && Math.abs(m[1] * m[2] * m[3] - m[4]) < 0.003, m ? m.slice(1).join(' ') : note.slice(0, 100));
    check('compare: four overlays and the differences table are drawn', await count(page, '.fc-compare .fc-chart canvas.fck-base') >= 4 && await count(page, '.fc-difftable tbody tr') >= 2, await count(page, '.fc-compare .fc-chart canvas.fck-base') + ' charts, ' + await count(page, '.fc-difftable tbody tr') + ' rows');
    const bad = await unreachable(page); check('Compare: every control can be reached with the keyboard', bad.length === 0, bad.join('; '));
  }
  if (fixture) check('names that try to be markup did not run on Runs and Compare', await pwned(page) === undefined, await pwned(page));
  await page.finish();
}

// =====================================================================================================
// 2b. Regressions: each defect once found is held here, on whatever hub is checked
// =====================================================================================================
{
  const page = await page0({ width: 1366, height: 768 });
  const active = () => page.evaluate(() => { const a = document.activeElement; return a ? a.tagName + '.' + String(a.className).split(' ')[0] + '[' + (a.getAttribute('data-fk') || a.getAttribute('data-k') || '') + ']' : 'none'; });
  await open(page, base + '#/runs');
  // --- the header row sticks; the row header of a run does not pile up under it ---
  const sticky = await page.evaluate(() => { const wrap = document.querySelector('.fc-hub__table'); wrap.scrollTop = 600; const head = document.querySelector('.fc-table--hub thead th.c-fc'), hr = head.getBoundingClientRect(), wr = wrap.getBoundingClientRect();
    return { scrolled: wrap.scrollTop, headTop: Math.round(hr.top - wr.top), bodySticky: [...document.querySelectorAll('.fc-table--hub tbody th')].filter(e => getComputedStyle(e).position === 'sticky').length, bg: getComputedStyle(head).backgroundColor }; });
  check('Runs: no row header is sticky (they piled up in the header cell)', sticky.bodySticky === 0, JSON.stringify(sticky));
  if (sticky.scrolled > 0) check('Runs: the header row stays at the top of the table and is opaque', sticky.headTop === 0 && !/rgba\(0, 0, 0, 0\)|transparent/.test(sticky.bg), JSON.stringify(sticky));
  // --- Back returns to the place in the table ---
  if (sticky.scrolled > 100) {
    const before = sticky.scrolled; await page.locator('.fc-table--hub tbody tr a.fc-id').nth(Math.min(N - 1, 16)).click(); await page.waitForSelector('.fc-runhead', { timeout: 20000 }); await settle(page, 300);
    await page.goBack(); await page.waitForSelector('.fc-table--hub'); await settle(page, 500);
    const after = await page.evaluate(() => document.querySelector('.fc-hub__table').scrollTop);
    check('Runs: Back from a run returns to the scroll position', Math.abs(after - before) < 80, before + ' -> ' + after);
  }
  // --- search: a verdict value that is none is an error, a verdict word as plain text gets a hint ---
  await open(page, base + '#/runs?q=status:wran');
  const bar = (await text(page, '.fc-filterbar'));
  check('search "status:wran": struck through, explained, never applied', await count(page, '.fc-filter.is-unknown') === 1 && /is not a value of status/.test(bar) && /Did you mean warn/.test(bar) && await hubRows(page) === Math.min(N, 400), bar.slice(0, 160));
  await open(page, base + '#/runs?q=Warning');
  if (await count(page, '.fc-hub__empty')) { const e = await text(page, '.fc-hub__empty'); const nW = rows.filter(r => r.verdict === 'warn').length;
    check('search "Warning" (the word on a badge): points to status:warn and says how many', /status:warn/.test(e) && e.includes('(' + nW + ' run'), e);
    await page.click('.fc-hub__empty .fc-linkbtn'); await settle(page, 300); check('following the hint shows the runs with a warning', await hubRows(page) === nW, await hubRows(page)); }
  // --- keyboard: the peek drawer takes focus, Escape gives it back; a ticked filter keeps focus ---
  await open(page, base + '#/runs');
  await page.focus('.fc-table--hub tbody tr [data-k="peek"]'); await page.keyboard.press('Enter'); await settle(page, 500);
  check('keyboard: opening the peek drawer moves focus into it', await page.evaluate(() => !!document.activeElement.closest('.fc-peek')), await active());
  await page.keyboard.press('Escape'); await settle(page, 500);
  check('keyboard: Escape closes the drawer and focus returns to its button', await count(page, '.fc-peek') === 0 && /\[peek\]/.test(await active()), await active());
  if (await count(page, '[data-fk="facets-open"]')) { await page.focus('[data-fk="facets-open"]'); await page.keyboard.press('Enter'); await settle(page, 400); check('keyboard: opening Filters keeps focus', /facets-hide/.test(await active()), await active()); }
  const fk = await page.evaluate(() => { const i = document.querySelector('.fc-facet input[type=checkbox]'); return i ? i.getAttribute('data-fk') : null; });
  if (fk) { await page.focus('[data-fk="' + fk + '"]'); await page.keyboard.press('Space'); await settle(page, 500); check('keyboard: ticking a filter keeps focus on it', (await active()).includes(fk), await active());
    await page.focus('.fc-filter__x'); const ring = await page.evaluate(() => { const st = getComputedStyle(document.activeElement); return st.outlineStyle + ' ' + st.outlineColor; });
    await page.keyboard.press('Enter'); await settle(page, 500); check('keyboard: removing a filter chip leaves focus on a control', !/^BODY/.test(await active()), await active()); }
  for (const sel of ['.fc-table--hub thead th button', '.fc-theme button']) { await page.focus(sel); await page.keyboard.press('Shift+Tab'); await page.keyboard.press('Tab');
    const ring = await page.evaluate(() => { const st = getComputedStyle(document.activeElement); return st.outlineStyle + ' ' + st.outlineWidth + ' ' + st.outlineColor; });
    check('a focus ring on ' + (sel.includes('thead') ? 'a sort header' : 'a segmented button'), /^solid/.test(ring) && !/rgba\(0, 0, 0, 0\)|transparent/.test(ring), ring); }
  // --- the Project filter, Group by project and the Projects page give the same numbers ---
  await open(page, base + '#/runs?group=project');
  const groups = await page.locator('.fc-table__group').evaluateAll(els => els.map(e => e.textContent));
  check('Group by project: one group per project of the index', groups.length === projects.length, groups.length + ' groups, ' + projects.length + ' projects');
  const bad = projects.filter(p => { const want = p.name == null ? 'No project' : p.name; const g = groups.find(t => t.startsWith(want + ' ')); return !g || !g.includes(' ' + p.runs + ' run'); });
  check('Group by project: every group holds the runs the Projects table counts', bad.length === 0, bad.slice(0, 3).map(p => p.id).join(' '));
  for (const p of projects.filter(p => p.name != null).slice(0, 4)) { await open(page, base + '#/runs?q=' + encodeURIComponent('project:' + p.id)); check('project:' + p.id.slice(0, 24) + ' matches that project only', await hubRows(page) === p.runs, (await hubRows(page)) + ' rows, ' + p.runs + ' runs'); }
  // --- addresses that are no page ---
  await open(page, base + '#/nonsense');
  check('an unknown address: the title says Not found and no tab is lit', /^Not found/.test(await page.title()) && await count(page, '.fc-tab.is-on') === 0, await page.title());
  const lone = [...new Set(rows.map(r => r.fc))].find(fc => rows.filter(r => r.fc === fc).length === 1);
  if (lone) { await open(page, base + '#/run/' + lone); await settle(page, 600); check('a flow cell id where a run key belongs opens its run', page.url().includes('#/run/' + rows.find(r => r.fc === lone).id), page.url().slice(-40)); }
  // --- compare: the bar follows the page; an unknown id is named; one run is "1 run" ---
  if (fullRuns.length >= 2) {
    await open(page, base + '#/compare?items=' + fullRuns.slice(0, 2).join(',') + ',XXX00000_00000000'); await page.waitForSelector('canvas', { timeout: 20000 }); await settle(page, 500);
    check('compare: an id that is no run is named, not dropped silently', /XXX00000_00000000 is not a run of this library/.test(await text(page, '.fc-page')));
    check('compare: the application bar counts the runs of the page', /Compare \(2\)/.test(await text(page, '.fc-appbar__compare')), await text(page, '.fc-appbar__compare'));
    await page.locator('.fc-item .fc-close').first().click(); await settle(page, 500);
    check('compare: removing a run updates the bar; one run reads "1 run"', /Compare \(1\)/.test(await text(page, '.fc-appbar__compare')) && /Compare 1 run( |$)/.test(await text(page, '.fc-crumbs')), await text(page, '.fc-crumbs'));
  }
  // --- a run: All checks across the card; numbers in hours; times name their zone; the section is in the address ---
  const rid = kinds.complete || kinds.full || fullRuns[0];
  if (rid) {
    const wide = await page0({ width: 1920, height: 1080 });
    await openRun(wide, rid); await wide.click('[data-fk="allchecks-run_health"]'); await settle(wide, 400);
    const m = await wide.evaluate(() => { const box = document.querySelector('.fc-allchecks'), card = document.querySelector('.fc-verdict'), tb = box.querySelector('table'), st = getComputedStyle(box);
      return { box: Math.round(box.getBoundingClientRect().width), card: Math.round(card.getBoundingClientRect().width), perRow: Math.round(tb.getBoundingClientRect().height / tb.querySelectorAll('tbody tr').length), maxH: st.maxHeight, of: st.overflowY }; });
    check('All checks: across the whole summary card, as tall as its rows, no row taller than three lines', m.box > m.card - 40 && m.maxH === 'none' && m.of === 'visible' && m.perRow < 70, JSON.stringify(m));
    await shot(wide, 'run_all_checks_1920_light');
    await openRun(wide, rid, '?more=n-timeline,n-fc,n-rtp,src');
    const body = (await wide.evaluate(() => document.body.innerText)).replace(/\u00a0/g, ' ');
    const hoursRows = body.split('\n').filter(l => /\t[a-z0-9_.]+_h$/.test(l) && !/_per_h$/.test(l));
    check('All numbers: every number whose id ends in _h is written in hours', hoursRows.length > 0 && hoursRows.every(l => /\t[\d.,]+ (h|min)\t/.test(l)), hoursRows.find(l => !/\t[\d.,]+ (h|min)\t/.test(l)) || hoursRows.length + ' rows');
    const src = await wide.evaluate(() => { const el = document.querySelector('#s-data_sources'); return el ? el.innerText : ''; });
    check('Data sources: sections and file kinds in the words of the page, no internal name', src.length > 50 && !/\b(funnel|len\.est|yield\.basecalled|quant\.[a-z]+|minknow\.[a-z_]+|dorado\.bam|align\.[a-z_]+)\b/.test(src), (src.match(/\b(funnel\S*|quant\.\S*|minknow\.\S*|align\.\S*)/g) || []).slice(0, 4).join(' '));
    check('Data sources: a date in a file name is not taken for a run id', !/run id 20\d{6},/.test(src));
    check('the start of a run names its time zone', /Started \d+ \w+ \d{4}, \d\d:\d\d UTC/.test(await text(wide, '.fc-runhead__line')), (await text(wide, '.fc-runhead__line')).slice(0, 60));
    await wide.locator('.fc-runrail a', { hasText: 'Flow cell' }).first().click(); await settle(wide, 300);
    check('side navigation: the section is in the address', /[?&]at=flow_cell/.test(wide.url()), wide.url().slice(-60));
    const still = await wide.evaluate(() => [...document.querySelectorAll('.fc-plot')].filter(p => p.tabIndex === 0 && !p.getAttribute('role')).length + [...document.querySelectorAll('.fc-plot--still')].filter(p => p.tabIndex === 0).length);
    check('a chart without a read-out takes no keyboard focus', still === 0, still);
    await wide.finish();
  }
  // --- text shown twice as large (every font size doubled): the table keeps whole columns, the bar grows, nothing pushes the page wider ---
  for (const [name, hash] of [['runs', '#/runs'], ['home', '#/home']].concat(rid ? [['run', '#/run/' + rid]] : [])) {
    const z = await page0({ width: 1366, height: 768 }); await open(z, base + hash); if (name === 'run') { await z.waitForSelector('.fc-verdict__headline', { timeout: 20000 }); await settle(z, 800); }
    const double = () => z.evaluate(() => { const els = [...document.querySelectorAll('body *')].filter(e => !e.dataset.z), sizes = els.map(e => parseFloat(getComputedStyle(e).fontSize)); els.forEach((e, i) => { e.style.fontSize = (sizes[i] * 2) + 'px'; e.dataset.z = 1; }); });
    await double(); await settle(z, 700); await double(); await settle(z, 500);
    const m = await z.evaluate(() => { const tabs = [...document.querySelectorAll('.fc-tab')].map(t => t.getBoundingClientRect()); let overlap = 0; for (let i = 1; i < tabs.length; i++) if (tabs[i].left < tabs[i - 1].right - 1 && Math.abs(tabs[i].top - tabs[i - 1].top) < 5) overlap++;
      return { sideways: document.documentElement.scrollWidth - document.documentElement.clientWidth, cut: [...document.querySelectorAll('.fc-table--hub tbody td.fc-num')].filter(td => td.scrollWidth > td.clientWidth + 1).length, overlap, cols: document.querySelectorAll('.fc-table--hub thead th').length }; });
    check('text twice as large, ' + name + ': no sideways scroll, no number cut, tabs apart', m.sideways <= 1 && m.cut === 0 && m.overlap === 0, JSON.stringify(m));
    await shot(z, 'zoom_text_' + name + '_1366_light'); await z.finish();
  }
  // --- times: the clock of the computer the page is read on, named ---
  if (rid) for (const tz of ['Asia/Tokyo', 'America/Los_Angeles']) {
    const ctxPage = await page0({ width: 1366, height: 768, timezoneId: tz }); await openRun(ctxPage, rid);
    const line = await text(ctxPage, '.fc-runhead__line'), started = Date.parse(rows.find(r => r.id === rid).started), want = new Intl.DateTimeFormat('en-GB', { timeZone: tz, hour: '2-digit', minute: '2-digit', hour12: false }).format(started);
    check('times in ' + tz + ': shown on that clock, with its zone', line.includes(', ' + want + ' UTC' + (tz === 'Asia/Tokyo' ? '+9' : '−')), line.slice(0, 50) + ' (want ' + want + ')');
    await ctxPage.finish();
  }
  await page.finish();
}

// =====================================================================================================
// 3. Projects and a project page
// =====================================================================================================
if (firstProject) {
  const page = await page0({ width: 1366, height: 768 });
  await open(page, base + '#/projects');
  check('Projects tab: one row per project', await count(page, '.fc-projects tbody tr') === projects.length, await count(page, '.fc-projects tbody tr'));
  await page.locator('.fc-projects tbody tr a').first().click(); await page.waitForSelector('.fc-projectpage'); await settle(page, 300);
  const head = await text(page, '.fc-projectpage .fc-runhead'), tiles = await text(page, '.fc-projectpage .fc-keynumbers');
  check('project page: header and numbers equal the index', head.includes(firstProject.runs + ' run') && head.includes(firstProject.flow_cells + ' flow cell') && tiles.includes(String(firstProject.runs)), head.slice(0, 120));
  check('project page: where it is, one row per folder', await count(page, '.fc-projectpage .fc-panel table tbody tr') >= Math.max(1, firstProject.roots.length), await count(page, '.fc-projectpage .fc-panel table tbody tr'));
  check('project page: the Runs table is filtered to the project', await hubRows(page) === firstProject.runs, (await hubRows(page)) + ' rows, project has ' + firstProject.runs);
  await page.goBack(); await page.waitForSelector('.fc-projectspage'); check('Back returns to Projects', new URL(page.url()).hash === '#/projects');
  if (fixture) check('names that try to be markup did not run on the project pages', await pwned(page) === undefined, await pwned(page));
  await page.finish();
}

// =====================================================================================================
// 4. Run views by state
// =====================================================================================================
async function openRun(page, id, extra = '') { await open(page, base + '#/run/' + id + extra); await page.waitForSelector('.fc-verdict__headline', { timeout: 20000 }); await settle(page, 900); }
for (const [kind, id] of Object.entries(kinds)) {
  if (!id) continue;
  const page = await page0({ width: 1366, height: 768 }), row = rows.find(r => r.id === id);
  await openRun(page, id);
  check(kind + ' run: the stored headline is shown word for word', (await text(page, '.fc-verdict__headline')) === row.headline.replace(/\s+/g, ' ').trim(), row.headline.slice(0, 80));
  const body = await text(page, '.fc-run');
  if (kind === 'stub') check('raw data only: says so, and draws no chart', /Only raw data was found/.test(body) && await count(page, '.fc-run canvas') === 0);
  else {
    await page.waitForSelector('.fc-track canvas', { timeout: 20000 });
    check(kind + ' run: the timeline is drawn', await count(page, '.fc-tracks canvas.fck-base') >= 4, await count(page, '.fc-tracks canvas.fck-base'));
    const parts = await page.locator('.fc-verdict__parthead .fc-badge').allInnerTexts();
    check(kind + ' run: two parts, Run health and Read data', parts.length === 2, parts.join(' | '));
    if (row.read_data === 'absent') check(kind + ' run: Read data reads "Not available yet" and one card says what will appear', parts[1] === 'Not available yet' && await count(page, '#s-read_data') === 1, parts[1]);
    if (kind === 'complete') {
      // every linked file was read, so nothing is marked partial, except a section the bundle itself calls partial
      // (read fate, when a stage such as "basecalled" has no file at all): that one says so, and only that one
      const own = await partialSections(id), marks = await count(page, '.fc-partial');
      check('complete read data carries no "partial" mark, except on a section that lacks a file', own.length ? marks === own.length && !/partial:/.test(body) : marks === 0 && !/PARTIAL|partial:/.test(body), own.join(', ') || 'none');
      check('complete read data: key numbers are FlowCellar\'s own (passed bases, read N50, quality)', /(Passed bases|Basecalled bases|Bases given to samples|Bases in full-length reads|Bases in aligned reads) [\d.,]+ (Gb|Mb)/.test(body) && /Read N50 [\d.,]+ k?b/.test(body) && /Median read quality Q [\d.]+/.test(body), (/Bases[^A-Z]*/.exec(body) || [''])[0]);
      const stage = row.t2_population; check('complete read data: "Passed" names only reads that passed a quality filter', stage === 'basecalled_pass' ? /Passed bases/.test(body) : !/Passed bases/.test(body), stage); }
    if (kind === 'partial') check('partial read data is marked: a line, tags on the charts and a word in the tile captions', await count(page, '.fc-partial') >= 1 && /PARTIAL ·/.test(body) && /partial:/.test(body), await count(page, '.fc-partial'));
    if (kind === 'quickRun') check('quick detail: "so far" and what waits for the read count', /so far/.test(parts[0]) && /being counted/i.test(body) && /wait for the read count/.test(body), parts[0]);
    if (kind === 'cut') check('a read table that is cut off: no "so far" on Run health, nothing "being counted", and the page says why', !/so far/.test(parts[0]) && !/being counted/i.test(body) && /(cut off|could not be read|read count is not complete)/.test(body), parts[0]);
    if (kind === 'running') check('a run in progress wears the Running badge', /Running/.test(parts[0]), parts[0]);
    if (kind === 'many') check('fourteen samples: the Samples table lists them all', await count(page, '.fc-samples tbody tr') >= row.t2_samples, await count(page, '.fc-samples tbody tr'));
    // open every drawer, then look for anything that reads like an identifier or a missing label
    for (let n = 0; n < 40 && await count(page, '.fc-runmain .fc-disclosure:not(.is-open)'); n++) { await page.locator('.fc-runmain .fc-disclosure:not(.is-open)').first().click({ timeout: 3000 }); await settle(page, 80); }
    await settle(page, 500);
    const all = await text(page, '.fc-run');
    check(kind + ' run: no "label missing", nothing imported or sampled', !/label missing|imported|sampled/i.test(all), (/label missing|imported|sampled/i.exec(all) || [''])[0]);
    const raw = await page.locator('.fc-runmain tr[title] th[scope=row]').evaluateAll(ths => ths.map(t => t.textContent.trim()).filter((t, k) => /^[a-z]+(_[a-z0-9]+)+$/.test(t) || t === ths[k].parentElement.getAttribute('title')));
    check(kind + ' run: no number is named by a bare identifier', raw.length === 0, raw.slice(0, 5).join(', '));
    const bad = await unreachable(page); check(kind + ' run: every control can be reached with the keyboard', bad.length === 0, bad.slice(0, 4).join('; '));
    check(kind + ' run: every chart offers "Show as table"', await page.evaluate(() => [...document.querySelectorAll('.fc-run .fc-chart, .fc-run .fc-track:not(.fc-track--na)')].filter(c => c.querySelector('canvas') && ![...c.querySelectorAll('button')].some(b => /table/i.test(b.textContent))).map(c => (c.querySelector('.fc-chart__title, .fc-track__title') || {}).textContent)).then(l => l.length ? l.join(', ') : '') === '', '');
    await shot(page, 'run_' + kind + '_all_1366_light', true);
  }
  if (fixture) check(kind + ' run: names that try to be markup did not run', await pwned(page) === undefined, await pwned(page));
  await page.finish();
}
// a reason leads to its evidence; hover and keyboard on charts; a drag sets the Time control
{
  const withReason = newest.find(r => (r.reason_checks || []).length && r.detail_t1 !== 'none');
  if (withReason) {
    const page = await page0({ width: 1366, height: 768 }); await openRun(page, withReason.id);
    await page.locator('.fc-reasons__go').first().click(); await settle(page, 900);
    check('a reason leads to its evidence and is a history entry', (await count(page, '.fc-section.is-evidence') === 1 || await count(page, '.is-evidence') >= 1) && /ev=/.test(page.url()), new URL(page.url()).hash);
    if (await count(page, 'button:has-text("Back to summary")')) { await page.click('button:has-text("Back to summary")'); await settle(page, 300); check('"Back to summary" returns to the top', await page.evaluate(() => window.scrollY) < 5); }
    await page.finish();
  }
  const id = kinds.complete || kinds.partial || kinds.full;
  if (id) {
    const page = await page0({ width: 1366, height: 768 }); await openRun(page, id); await page.waitForSelector('.fc-track canvas');
    const out = page.locator('#p-yield-output_over_time .fc-plot'); await out.scrollIntoViewIfNeeded(); await settle(page, 200);
    const box = await out.boundingBox(); await page.mouse.move(box.x + box.width * 0.5, box.y + box.height * 0.5); await settle(page, 250);
    const v1 = await text(page, '#p-yield-output_over_time .fc-track__value'), v2 = await text(page, '#p-fc-occupancy .fc-track__value');
    check('hover: the Output track prints its value and the Occupancy track follows the same cursor', /Hour/.test(v1) && /Hour/.test(v2), v1 + ' | ' + v2);
    await page.mouse.move(box.x + box.width * 0.30, box.y + box.height * 0.5); await page.mouse.down(); await page.mouse.move(box.x + box.width * 0.62, box.y + box.height * 0.5, { steps: 6 }); await page.mouse.up(); await settle(page, 500);
    check('a drag across a track sets the Time control', /[?&]t=\d/.test(page.url()), new URL(page.url()).hash);
    await out.dblclick(); await settle(page, 400); check('a double-click returns to the whole run', !/[?&]t=\d/.test(page.url()));
    const len = page.locator('#p-len-estimated_length .fc-plot'); await len.scrollIntoViewIfNeeded(); await len.focus(); await page.keyboard.press('ArrowRight'); await page.keyboard.press('ArrowRight'); await settle(page, 150);
    const said = await len.getAttribute('aria-valuetext'); check('keyboard: arrow keys move a chart cursor and the read-out is its spoken text', /bases/.test(said || ''), said);
    const evened = await count(page, '#p-fc-channel_activity button:has-text("Show as measured")');      // offered only when a pore scan left its mark on a 10-minute step
    if (evened) { await page.locator('#p-fc-channel_activity button:has-text("Show as measured")').click(); await settle(page, 300);
      check('channel activity: one control shows the pore-scan steps as measured', /scans=raw/.test(page.url()) && await count(page, '#p-fc-channel_activity button:has-text("Even out pore scans")') === 1); }
    await page.locator('#p-fc-channel_activity button:has-text("Show as table")').click(); await settle(page, 200);
    check('channel activity: the table lists the measured steps' + (evened ? ' and marks the scan steps' : ''), await count(page, '#p-fc-channel_activity .fc-charttable tbody tr') >= 1 && (!evened || (await page.locator('#p-fc-channel_activity .fc-charttable tbody').innerText()).includes('yes')), await count(page, '#p-fc-channel_activity .fc-charttable tbody tr') + ' rows');
    const map = page.locator('#p-fc-channel_map .fc-plot'); await map.scrollIntoViewIfNeeded(); const mb = await map.boundingBox(); await page.mouse.move(mb.x + mb.width * 0.4, mb.y + mb.height * 0.5); await settle(page, 200);
    check('flow cell map: pointing at a cell names its channel, and the map says how positions are known', /Channel [\d,]+/.test(await text(page, '#p-fc-channel_map .fc-chart__readout')) && /positions inferred|channel order/.test(await text(page, '#p-fc-channel_map')), await text(page, '#p-fc-channel_map .fc-chart__readout'));
    await page.finish();
    // the controls drive the read-data charts
    const row = rows.find(r => r.id === id);
    if (row.read_data !== 'absent') {
      const p2 = await page0({ width: 1366, height: 768 }); await openRun(p2, id, '?by=bases&y=log');
      const samples = p2.locator('.fc-runcontrols select').first(), opt1 = await samples.locator('option').nth(1).getAttribute('value'); await samples.selectOption(opt1); await settle(p2, 400);
      check('the Samples control narrows the read-data charts to one sample', /[?&]s=/.test(p2.url()) && (await text(p2, '#p-len .fc-chart__scope')).length > 0, await text(p2, '#p-len .fc-chart__scope'));
      await p2.click('.fc-runcontrols button:has-text("Custom")'); await settle(p2, 300);
      if (await count(p2, '.fc-custom input')) { await p2.fill('.fc-custom input >> nth=1', '1000'); await settle(p2, 500);
        check('the Custom filter reports the reads it keeps', await count(p2, '.fc-evidence:has-text("Custom filter")') === 1, (await p2.locator('.fc-evidence').first().innerText().catch(() => '')).slice(0, 110)); }
      if (kind_has_fail(row)) { await p2.click('.fc-runcontrols button:has-text("All")'); await settle(p2, 400); check('Reads: All adds the reads that failed the filter', /reads=all/.test(p2.url()) && /all reads, passed and failed/.test(await text(p2, '#p-len .fc-chart__scope')), await text(p2, '#p-len .fc-chart__scope')); }
      await p2.evaluate(() => document.getElementById('s-length_quality').scrollIntoView()); await settle(p2, 300); await shot(p2, 'run_controls_1366_light');
      await p2.finish();
    }
  }
}
function kind_has_fail(row) { return fixture && row.t2_population === 'basecalled_all'; }

// =====================================================================================================
// 5. Flow cell history
// =====================================================================================================
if (reused) {
  const page = await page0({ width: 1366, height: 768 }), n = rows.filter(r => r.fc === reused).length;
  await open(page, base + '#/flowcell/' + reused); await page.waitForSelector('.fc-fcruns', { timeout: 20000 }); await settle(page, 900);
  const t = await text(page, '.fc-flowcell');
  check('flow cell history: every run of the flow cell, its loadings, one map per run that has one', t.includes(n + ' runs') && await count(page, '.fc-fcruns tbody tr') === n && /loading/.test(t), (/Used[^.]*\./.exec(t) || [''])[0]);
  const bad = await unreachable(page); check('Flow cell history: every control can be reached with the keyboard', bad.length === 0, bad.join('; '));
  await page.finish();
}

// =====================================================================================================
// 6. Every screen, both themes, three widths: screenshots, no sideways scroll, no half column
// =====================================================================================================
const screens = [['home', '#/home', '.fc-home'], ['runs', '#/runs', '.fc-table--hub'], ['projects', '#/projects', '.fc-projectspage']];
if (firstProject) screens.push(['project', '#/project/' + encodeURIComponent(firstProject.id), '.fc-projectpage']);
for (const k of ['complete', 'partial', 'quickRun', 'full', 'stub']) if (kinds[k]) screens.push(['run_' + k, '#/run/' + kinds[k], '.fc-verdict__headline']);
if (fullRuns.length >= 2) screens.push(['compare', '#/compare?items=' + fullRuns.slice(0, 2).join(','), '.fc-compare .fc-plot canvas']);
if (reused) screens.push(['flowcell', '#/flowcell/' + reused, '.fc-fcruns']);
for (const w of widths) for (const theme of themes) {
  const page = await page0({ width: w, height: H(w), scheme: 'light' });
  await open(page, base + '#/home'); await setTheme(page, theme);
  check('theme ' + theme + ' at ' + w + ' px is set on the page', await page.evaluate(() => document.documentElement.getAttribute('data-theme')) === theme.toLowerCase());
  for (const [name, hash, ready] of screens) {
    await page.goto(base + hash); await page.waitForSelector(ready, { timeout: 20000 }); await settle(page, name.startsWith('run_') || name === 'compare' || name === 'flowcell' ? 900 : 350);
    const wide = await sideways(page);
    check(name + ' at ' + w + ' px, ' + theme.toLowerCase() + ': the page does not scroll sideways', wide <= 1, wide);
    if (name === 'runs' || name === 'project') { const hc = await halfColumn(page); check(name + ' at ' + w + ' px, ' + theme.toLowerCase() + ': every column of the Runs table is whole', hc <= 1, hc); }
    if (name === 'run_complete' || name === 'run_partial' || name === 'run_full') { const cols = await page.evaluate(() => { const b = document.querySelector('.fc-runbody'); return b ? getComputedStyle(b).gridTemplateColumns.split(' ').length : 0; }); check(name + ' at ' + w + ' px: ' + (w < 900 ? 'one column' : 'rail and sections'), w < 900 ? cols === 1 : cols === 2, cols); }
    await shot(page, name + '_' + w + '_' + theme.toLowerCase(), name === 'home' || name === 'project' || name === 'projects');
    if (name.startsWith('run_') && name !== 'run_stub' && theme === 'Light' && w !== 1280) { await page.evaluate(() => { const s = document.getElementById('s-flow_cell'); if (s) s.scrollIntoView(); }); await settle(page, 400); await shot(page, name + '_flowcell_' + w + '_light'); }
  }
  await page.finish();
}
// =====================================================================================================
// 6b. One width for every screen, and a pointer on a chart that moves nothing
// =====================================================================================================
{
  // --- a wide window: the content of the application bar and of every screen start and end on the same two lines ---
  const page = await page0({ width: 2560, height: 1440 }), found = [];
  for (const [name, hash, ready] of screens) {
    await page.goto(base + hash); await page.waitForSelector(ready, { timeout: 20000 }); await settle(page, name.startsWith('run_') || name === 'compare' || name === 'flowcell' ? 900 : 350);
    found.push({ name, ...await page.evaluate(() => { const W = document.documentElement.clientWidth, shown = el => el.getBoundingClientRect().width > 0;
      const box = el => { const r = el.getBoundingClientRect(), cs = getComputedStyle(el); return Math.round(r.left + parseFloat(cs.paddingLeft) + parseFloat(cs.borderLeftWidth)) + '/' + Math.round(W - r.right + parseFloat(cs.paddingRight) + parseFloat(cs.borderRightWidth)); };
      const wide = sel => Math.max(0, ...[...document.querySelectorAll(sel)].map(e => Math.round(e.getBoundingClientRect().width)));
      return { bar: box(document.querySelector('.fc-appbar__in')), pages: [...document.querySelectorAll('#main .fc-page, #main .fc-hub > .fc-viewbar, #main .fc-hub > .fc-filterbar')].filter(shown).map(box),
        sideways: document.documentElement.scrollWidth - W, chart: wide('#main .fc-section .fc-chart:not(.fc-chart--map):not(.fc-chart--flat) > .fc-plot, #main .fc-compare .fc-chart > .fc-plot'), tile: wide('#main .fc-tile') }; }) });
  }
  const lines = new Set(found.flatMap(f => [f.bar, ...f.pages])), gutter = [...lines][0].split('/').map(Number);
  check('at 2,560 px the application bar and the content of every screen start and end on the same two lines', lines.size === 1 && found.every(f => f.pages.length >= 1), lines.size === 1 ? found.length + ' screens, ' + [...lines][0] + ' px' : found.map(f => f.name + ' ' + f.bar + ' | ' + f.pages.join(' ')).join(' · '));
  check('at 2,560 px no screen is a centred column: the content reaches both gutters', lines.size === 1 && gutter[0] === gutter[1] && gutter[0] <= 24, [...lines].join(' vs '));
  check('at 2,560 px no screen scrolls sideways', found.every(f => f.sideways <= 1), found.filter(f => f.sideways > 1).map(f => f.name + ' ' + f.sideways).join(' '));
  // A row of charts fills its section in columns of one width (charts/Grid.jsx): a chart is at most 1,200 px wide
  // (--fc-chart-wide), so that a row of two reaches the right edge at this width like the rows of three and more.
  check('at 2,560 px the width is used for columns: no chart of a run or of Compare is wider than 1,200 px, no key number wider than 420 px', found.every(f => f.chart <= 1201 && f.tile <= 421), found.filter(f => f.chart > 1201 || f.tile > 421).map(f => f.name + ' chart ' + f.chart + ' tile ' + f.tile).join(' · '));
  await page.finish();
}
// --- a pointer that crosses a chart writes the read-out and draws the cursor, and nothing else changes: the box of
//     every element of the page is compared with its box before the pointer came (apart from the text inside a read-out
//     and the note beside the cell of a map, which lies over the plot) ---
if (!quick) {
  const SNAP = () => { const box = () => { const m = new Map(); for (const e of document.querySelectorAll('.fc-appbar, .fc-appbar *, #main, #main *')) { const ro = e.closest('.fc-chart__readout, .fc-track__value'); if ((ro && ro !== e) || e.closest('.fc-tip')) continue; const r = e.getBoundingClientRect(); m.set(e, [r.left, r.top, r.width, r.height]); } return m; };
    window.__box = box; window.__was = box(); window.__tall = document.documentElement.scrollHeight; };
  const MOVED = async k => { await new Promise(r => requestAnimationFrame(() => requestAnimationFrame(r))); const now = window.__box(); let n = 0, first = '';
    for (const [e, r] of now) { const b = window.__was.get(e); if (b && !b.some((v, i) => Math.abs(v - r[i]) > 0.01)) continue; n++; if (!first) first = e.tagName.toLowerCase() + '.' + String(e.className).split(' ')[0] + (b ? ' ' + b.map(Math.round).join(',') + ' -> ' + r.map(Math.round).join(',') : ' is new'); }
    if (now.size !== window.__was.size) n++; if (document.documentElement.scrollHeight !== window.__tall) { n++; first = first || 'the page is ' + (document.documentElement.scrollHeight - window.__tall) + ' px taller'; }
    const ro = window.__charts[k].querySelector('.fc-chart__readout, .fc-track__value'); return { n, first, text: ro ? ro.textContent : '' }; };
  const rid = kinds.complete || kinds.partial || kinds.full || fullRuns[0], hovered = [['Home', '#/home', '.fc-home .fc-chart canvas', 1]];
  if (rid) hovered.push(['a run', '#/run/' + rid + '?more=tracks,fc,rtp', '.fc-track canvas', 6]);
  if (fullRuns.length >= 2) hovered.push(['Compare', '#/compare?items=' + fullRuns.slice(0, 2).join(','), '.fc-compare .fc-plot canvas', 3]);
  const page = await page0({ width: 1366, height: 900 });
  for (const [name, hash, ready, least] of hovered) {
    await page.goto(base + hash); await page.waitForSelector(ready, { timeout: 20000 }); await settle(page, 1000);
    const n = await page.evaluate(() => { window.__charts = [...document.querySelectorAll('#main .fc-chart, #main .fc-track')].filter(c => c.querySelector('.fc-plot canvas.fck-over')); return window.__charts.length; });
    const bad = []; let steps = 0, wrote = 0;
    for (let k = 0; k < n; k++) {
      const c = await page.evaluate(async k => { const el = window.__charts[k]; el.scrollIntoView({ block: 'center' }); await new Promise(r => requestAnimationFrame(() => requestAnimationFrame(r))); const r = el.querySelector('.fc-plot').getBoundingClientRect(), t = el.querySelector('.fc-chart__title, .fc-track__title'); return { x: r.left, y: r.top, w: r.width, h: r.height, title: t ? t.textContent : '' }; }, k);
      await page.mouse.move(1, 1); await settle(page, 80); await page.evaluate(SNAP);
      const idle = (await page.evaluate(MOVED, k)).text; let worst = null, said = false;
      // from outside the chart, across it on a diagonal, and out again: 24 steps
      for (let i = 0; i <= 24; i++) { await page.mouse.move(c.x - 20 + (c.w + 40) * i / 24, Math.max(2, Math.min(898, c.y - 10 + (c.h + 20) * i / 24))); const m = await page.evaluate(MOVED, k); steps++; if (m.text !== idle) said = true; if (m.n && !worst) worst = c.title + ': ' + m.n + ' moved, ' + m.first; }
      await page.mouse.move(1, 1); const back = await page.evaluate(MOVED, k);
      if (back.n && !worst) worst = c.title + ': ' + back.n + ' still moved after the pointer left, ' + back.first; if (back.text !== idle && !worst) worst = c.title + ': the read-out did not go back to its idle text';
      if (worst) bad.push(worst); if (said) wrote++;
    }
    check('hover on ' + name + ': the pointer crosses ' + n + ' charts in ' + steps + ' steps and no element of the page changes its size or its place', n > 0 && bad.length === 0, bad.length ? bad.length + ' charts: ' + bad[0] : '');
    check('hover on ' + name + ': the read-outs are written while the pointer is on the chart, and go back when it leaves', wrote >= least, wrote + ' of ' + n + ' charts');
  }
  await page.finish();
}
// the system's dark setting with the switch on Auto
{ const page = await newPage({ width: 1366, height: 768, scheme: 'dark', allowed: u => u.startsWith(origin) }); await open(page, base + '#/home'); const bg = await page.evaluate(() => getComputedStyle(document.body).backgroundColor); const m = /(\d+), (\d+), (\d+)/.exec(bg); check('Auto follows the system\'s dark setting', m && +m[1] + +m[2] + +m[3] < 200, bg); await page.finish(); }

// =====================================================================================================
// 7. Fixture only: names that try to be markup, Reload, a missing data file, a missing index, an empty library
// =====================================================================================================
if (fixture) {
  const hostile = rows.find(r => /<script>/.test(r.run_name || ''));
  if (hostile) {
    const page = await page0({ width: 1366, height: 768 }); await open(page, base + '#/run/' + hostile.id); await page.waitForSelector('.fc-track canvas', { timeout: 20000 }); await settle(page, 600);
    check('a run name that is markup is shown as text', (await page.locator('.fc-runhead__name').innerText()).includes('<script>') && await count(page, '.fc-runhead__name script, .fc-runhead__name img') === 0);
    for (let n = 0; n < 40 && await count(page, '.fc-runmain .fc-disclosure:not(.is-open)'); n++) { await page.locator('.fc-runmain .fc-disclosure:not(.is-open)').first().click({ timeout: 3000 }); await settle(page, 60); }
    check('sample, project and message texts that are markup are shown as text and nothing ran', await pwned(page) === undefined && (await text(page, '.fc-samples')).includes('<svg onload') && await count(page, '.fc-run svg[onload], .fc-run img[onerror], #main script') === 0, await pwned(page));
    await page.fill('.fc-appbar input[type=search]', '<img src=x onerror=window.__pwned=9>'); await page.keyboard.press('Enter'); await page.waitForSelector('.fc-table--hub'); await settle(page, 300);
    check('a search text that is markup is a filter chip, not markup', await pwned(page) === undefined && await count(page, '.fc-filter img') === 0);
    await page.finish();
    const p2 = await page0({ width: 1366, height: 768 }); await open(p2, base + '#/run/' + encodeURIComponent('"><img src=x onerror=window.__pwned=7>')); await settle(p2, 600);
    check('an address that is markup shows "not in the library" and runs nothing', await pwned(p2) === undefined && /not in the library/.test(await text(p2, '#main')));
    await p2.goto(base + '#/project/' + encodeURIComponent('<svg onload=window.__pwned=8>')); await settle(p2, 400);
    check('a project address that is markup runs nothing', await pwned(p2) === undefined && /not in the library/.test(await text(p2, '#main')));
    await p2.finish();
  }
  if (!http) {
    // Reload: a newer index.js beside the page replaces the library on request
    const next = path.join(hubDir, 'index.next.js'), cur = path.join(hubDir, 'index.js');
    if (fs.existsSync(next)) {
      const page = await page0({ width: 1366, height: 768 }); await open(page, base + '#/home');
      const keep = fs.readFileSync(cur);
      try { fs.copyFileSync(next, cur); await page.click('.fc-updated button:has-text("Reload")'); await page.waitForSelector('.fc-libnotices', { timeout: 8000 }).catch(() => {}); await settle(page, 400);
        check('Reload reads a newer index and shows what the scan is doing', /Counting reads/.test(await page.locator('.fc-home').innerText()), (await page.locator('.fc-libnotices').innerText().catch(() => '')).slice(0, 120)); }
      finally { fs.writeFileSync(cur, keep); }
      await page.finish();
    }
    // a data file that is missing leaves the run on its summary
    const victim = kinds.full || kinds.complete, file = path.join(hubDir, 'data', victim + '.js'), keep = fs.readFileSync(file);
    try { fs.unlinkSync(file); const page = await newPage({ width: 1366, height: 768, allowed: u => u.startsWith(origin) }); await page.goto(base + '#/run/' + victim); await page.waitForSelector('.fc-runfail', { timeout: 20000 });
      check('a missing data file leaves the run on "summary only"', /could not be read/.test(await text(page, '.fc-runfail')) && await count(page, '.fc-verdict__headline') === 1, await text(page, '.fc-runfail'));
      log.errors.splice(0, log.errors.length, ...log.errors.filter(e => !(e.url.includes(victim))));       // the browser reports the missing file; that is the point of this check
      await page.finish(); } finally { fs.writeFileSync(file, keep); }
    // a damaged data file is refused
    try { const t = keep.toString('utf8'), i = t.indexOf('","') - 40; fs.writeFileSync(file, t.slice(0, i) + (t[i] === 'A' ? 'B' : 'A') + t.slice(i + 1));
      const page = await newPage({ width: 1366, height: 768, allowed: u => u.startsWith(origin) }); await page.goto(base + '#/run/' + victim); await page.waitForSelector('.fc-runfail', { timeout: 20000 });
      check('a damaged data file is refused, with a banner', /damaged|could not be read/.test(await text(page, '.fc-runfail')), await text(page, '.fc-runfail'));
      log.errors.splice(0, log.errors.length, ...log.errors.filter(e => !(e.url.includes(victim)))); await page.finish(); } finally { fs.writeFileSync(file, keep); }
  }
}

// =====================================================================================================
// 8. A one-run report
// =====================================================================================================
if (opt('report')) {
  const file = path.resolve(opt('report')), url = fileUrl(file);
  const page = await newPage({ width: 1366, height: 768, allowed: u => u.startsWith(url) });
  await page.goto(url); await page.waitForSelector('.fc-appbar', { timeout: 20000 }); await page.waitForSelector('.fc-verdict__headline'); await settle(page, 900);
  check('report: opens on its run, with no library tabs', await count(page, '.fc-run') === 1 && await count(page, '.fc-appbar .fc-tab') === 0 && /Run report/.test(await text(page, '.fc-appbar')), await text(page, '.fc-appbar'));
  check('report: the charts are drawn from the data inside the file', await count(page, '.fc-tracks canvas.fck-base') >= 4);
  check('report: names of things outside the file are not links', await count(page, '.fc-crumbs a') === 0 && await count(page, '.fc-runhead__chips a') === 0);
  if (await count(page, '.fc-reasons__go')) { await page.locator('.fc-reasons__go').first().click(); await settle(page, 600); check('report: a reason leads to its evidence', /ev=/.test(page.url()) && await count(page, '.fc-run') === 1); }
  for (const w of widths) { await page.setViewportSize({ width: w, height: H(w) }); await page.goto(url + '#/'); await page.waitForSelector('.fc-track canvas'); await settle(page, 700); check('report at ' + w + ' px: no sideways scroll', await sideways(page) <= 1, await sideways(page)); await shot(page, 'report_' + w + '_light'); }
  await setTheme(page, 'Dark'); await shot(page, 'report_820_dark');
  if (fixture) check('report: names that try to be markup did not run', await pwned(page) === undefined, await pwned(page));
  await page.finish();
  // with scripts off, the verdict and the key numbers are readable as plain text
  const ctx = await browser.newContext({ viewport: { width: 1366, height: 768 }, javaScriptEnabled: false }), p2 = await ctx.newPage();
  await p2.goto(url); const plain = (await p2.locator('body').innerText()).replace(/\s+/g, ' ');
  check('report with scripts off: verdict, sentence and key numbers are readable', /FlowCellar run report/i.test(plain) && /(Pass|Warning|Fail|Incomplete|Running)/.test(plain) && /Reads|Pores|raw data/i.test(plain), plain.slice(0, 200));
  await p2.screenshot({ path: path.join(outDir, 'report_noscript.png') }); report.shots.push('report_noscript.png'); await ctx.close();
}

await browser.close();
report.notices = log.notices; report.finished = new Date().toISOString();
const failed = report.checks.filter(c => !c.ok);
report.summary = { checks: report.checks.length, failed: failed.length, errors: log.errors.length, requests_to_another_origin: log.requests.length, screenshots: report.shots.length, browser_notices: log.notices };
fs.writeFileSync(path.join(outDir, 'report.json'), JSON.stringify(report, null, 1));
console.log('\n' + report.checks.length + ' checks, ' + failed.length + ' failed · page and console errors and policy violations: ' + log.errors.length + ' · requests to another origin: ' + log.requests.length + ' · ' + report.shots.length + ' screenshots in ' + outDir + ' · browser notices: ' + log.notices);
for (const e of log.errors.slice(0, 15)) console.log('  ' + e.kind + ' @ ' + e.url.slice(-60) + ': ' + e.text.slice(0, 300));
for (const r of log.requests.slice(0, 10)) console.log('  request: ' + r);
process.exit(failed.length || log.errors.length || log.requests.length ? 1 : 0);
