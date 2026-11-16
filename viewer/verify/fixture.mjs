// The browser checks on synthetic data, end to end: make the fixture with the package's Python (no Node at run
// time is involved in writing pages), then use it from file://, over HTTP, as a report, empty and without an index.
//
//   node verify/fixture.mjs [--out DIR] [--quick] [--keep]
//
// Python: FLOWCELLAR_PYTHON, else <repo>/.venv/bin/python, else python3. Screenshots and report.json land in
// DIR/shots-file, DIR/shots-http (default DIR: viewer/.cache/verify).
import { spawn, spawnSync } from 'node:child_process';
import fs from 'node:fs';
import path from 'node:path';
import { launch, fileUrl, settle, repo, viewerDir, here } from './browser.mjs';

const args = process.argv.slice(2), opt = n => { const i = args.indexOf('--' + n); return i >= 0 ? args[i + 1] : null; }, flag = n => args.includes('--' + n);
const out = path.resolve(opt('out') || path.join(viewerDir, '.cache', 'verify'));
const venv = path.join(repo, '.venv', 'bin', 'python'), python = process.env.FLOWCELLAR_PYTHON || (fs.existsSync(venv) ? venv : 'python3');
const py = (...a) => { const r = spawnSync(python, a, { cwd: repo, encoding: 'utf8' }); if (r.status !== 0) { console.error(r.stderr || r.stdout); process.exit(2); } return r.stdout; };
const node = (...a) => spawnSync(process.execPath, a, { cwd: viewerDir, stdio: 'inherit' }).status;
const made = JSON.parse(py(path.join(repo, 'tests', 'viewer', 'viewer_make_fixture.py'), path.join(out, 'fixture')).trim().split('\n').pop());
const common = flag('quick') ? ['--quick'] : [];
let failed = 0;

console.log('\n=== the hub of the synthetic library, opened from file:// ===');
failed += node(path.join(here, 'check.mjs'), '--hub', made.hub, '--out', path.join(out, 'shots-file'), '--fixture', '--label', 'file', '--report', made.reports.TST00001_0a1b2c3d, ...common) ? 1 : 0;

console.log('\n=== the same library over HTTP (serve_file behind a test server) ===');
const srv = spawn(python, [path.join(repo, 'tests', 'viewer', 'viewer_serve_stub.py'), made.library], { cwd: repo, stdio: ['ignore', 'pipe', 'inherit'] });
const address = await new Promise((resolve, reject) => { let buf = ''; srv.stdout.on('data', d => { buf += d; const m = /http:\/\/127\.0\.0\.1:\d+\//.exec(buf); if (m) resolve(m[0]); }); srv.on('exit', () => reject(new Error('the test server stopped'))); setTimeout(() => reject(new Error('the test server did not start')), 15000); });
try { failed += node(path.join(here, 'check.mjs'), '--hub', address, '--out', path.join(out, 'shots-http'), '--fixture', '--label', 'http', '--quick') ? 1 : 0; }
finally { srv.kill(); }

console.log('\n=== an empty library, a folder without an index, reports of other states ===');
{
  const { browser, log, newPage } = await launch(); let bad = 0;
  const check = (name, ok, detail) => { if (!ok) bad++; console.log((ok ? 'ok   ' : 'FAIL ') + name + (detail ? '  [' + String(detail).slice(0, 160) + ']' : '')); };
  const shots = path.join(out, 'shots-file'); fs.mkdirSync(shots, { recursive: true });
  for (const scheme of ['light', 'dark']) {
    const base = fileUrl(path.join(made.empty, 'index.html')), page = await newPage({ scheme, allowed: u => u.startsWith(fileUrl(made.empty)) });
    await page.goto(base); await page.waitForSelector('.fc-home__empty', { timeout: 20000 }); await settle(page);
    const t = (await page.locator('.fc-home').innerText()).replace(/\s+/g, ' ');
    check('empty library (' + scheme + '): says what a library is and shows the scan command, with no empty table', /A library is a folder/.test(t) && /flowcellar scan \/folder\/with\/runs/.test(t) && await page.locator('.fc-home table').count() === 0, t.slice(0, 110));
    await page.screenshot({ path: path.join(shots, 'home_empty_1366_' + scheme + '.png') });
    await page.goto(base + '#/runs'); await page.waitForSelector('.fc-hub__empty'); check('empty library (' + scheme + '): Runs says there are none yet', /No runs here yet/.test(await page.locator('.fc-hub__empty').innerText()));
    await page.goto(base + '#/projects'); await page.waitForSelector('.fc-projectspage'); check('empty library (' + scheme + '): Projects says there is none yet', /No project yet/.test(await page.locator('.fc-projectspage').innerText()));
    await page.finish();
  }
  {
    const page = await newPage({ allowed: u => u.startsWith(fileUrl(made.noindex)) });
    await page.goto(fileUrl(path.join(made.noindex, 'index.html'))); await page.waitForSelector('.fc-boot--failed', { timeout: 20000 });
    check('a folder without an index says so', (await page.locator('.fc-boot--failed').innerText()) === 'This folder holds no library index. Run flowcellar hub.', await page.locator('.fc-boot--failed').innerText());
    await page.screenshot({ path: path.join(shots, 'noindex_1366_light.png') });
    log.errors.splice(0, log.errors.length, ...log.errors.filter(e => !/index\.js/.test(e.text) && !e.url.includes('noindex')));
    await page.finish();
  }
  for (const [key, file] of Object.entries(made.reports)) {
    const url = fileUrl(file), page = await newPage({ allowed: u => u.startsWith(url) });
    await page.goto(url); await page.waitForSelector('.fc-verdict__headline', { timeout: 20000 }); await settle(page, 700);
    const raw = fs.readFileSync(file, 'utf8');
    check('report of ' + key + ': opens, names its run, runs nothing it should not', (await page.locator('.fc-runhead__id').innerText()) === key.split('_')[0] && await page.evaluate(() => window.__pwned) === undefined && !/<script>window\.__pwned/.test(raw));
    await page.screenshot({ path: path.join(shots, 'report_' + key.split('_')[0] + '_1366_light.png') });
    await page.finish();
  }
  await browser.close();
  console.log('errors: ' + log.errors.length + ' · requests to another origin: ' + log.requests.length);
  for (const e of log.errors.slice(0, 10)) console.log('  ' + e.kind + ': ' + e.text.slice(0, 200));
  if (bad || log.errors.length || log.requests.length) failed++;
}
if (!flag('keep')) fs.rmSync(path.join(out, 'fixture'), { recursive: true, force: true });
console.log(failed ? '\nFAILED: ' + failed + ' part(s) of the browser checks' : '\nAll browser checks passed. Screenshots: ' + out);
process.exit(failed ? 1 : 0);
