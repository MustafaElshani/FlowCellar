// Screenshots of a hub for the README: whole regions of a page and single charts.
//
//   node tools/readme_shots.mjs --hub LIBRARY/hub --out DIR --spec SPEC.json
//
// SPEC.json is a list of pages; tools/make_readme_pictures.py writes it and calls this script:
//
//   { "route": "#/home", "width": 1120, "height": 900, "scale": 1.6, "scheme": "light",
//     "click": ["text=Space Ranger"],                 optional: things to press first
//     "shots": [ { "name": "home", "from": null, "to": ".fc-home__block--projects", "pad": 12 },    a region of the page
//                { "name": "len", "selector": "#p-len" } ] }                                        one element
//
// A region with "tight": true is only as wide as its "from" element. { "name": "x", "pdf": true } writes x.pdf, the
// PDF the browser makes of the page (what "Save as PDF" gives; Chromium only).
//
// A route may hold {run:TST00106} (the run of that flow cell; {run:TST00103:2} its second one) and
// {analysis:visium_hd_b:spaceranger}: they are looked up in the hub's index, so the script needs no run key.
// The browser is Playwright's, from viewer/node_modules, with the browsers of viewer/pw-browsers (Chromium when it is
// there, else Firefox; FLOWCELLAR_BROWSER_ENGINE chooses). Nothing is fetched: a hub page reads only its own files.
import fs from 'node:fs';
import path from 'node:path';
import zlib from 'node:zlib';
import { createRequire } from 'node:module';
import { fileURLToPath, pathToFileURL } from 'node:url';

const repo = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const args = process.argv.slice(2), opt = n => { const i = args.indexOf('--' + n); return i >= 0 ? args[i + 1] : null; };
const hub = opt('hub'), outDir = opt('out'), specPath = opt('spec');
if (!hub || !outDir || !specPath) { console.error('usage: node tools/readme_shots.mjs --hub LIBRARY/hub --out DIR --spec SPEC.json'); process.exit(2); }
const browsers = process.env.PLAYWRIGHT_BROWSERS_PATH || path.join(repo, 'viewer', 'pw-browsers');
process.env.PLAYWRIGHT_BROWSERS_PATH = browsers;
const has = name => fs.existsSync(browsers) && fs.readdirSync(browsers).some(n => n.startsWith(name + '-'));
const engine = (process.env.FLOWCELLAR_BROWSER_ENGINE || (has('chromium') ? 'chromium' : 'firefox')).toLowerCase();
if (!has(engine)) { console.error(`No Playwright ${engine} under ${browsers}. Install one: cd viewer && PLAYWRIGHT_BROWSERS_PATH=pw-browsers npx playwright install ${engine}`); process.exit(2); }
const pw = createRequire(path.join(repo, 'viewer', 'package.json'))('playwright');

// ---- the hub's index: which run a flow cell id means, which analysis a sample name
const hubDir = path.resolve(hub);
const m = /^__FC_INDEX__\("([A-Za-z0-9+/=]*)"\);\n$/.exec(fs.readFileSync(path.join(hubDir, 'index.js'), 'utf8'));
if (!m) { console.error(`${hubDir}/index.js is not the index of a hub`); process.exit(2); }
const ix = JSON.parse(zlib.gunzipSync(Buffer.from(m[1], 'base64')).toString('utf8'));
const C = ix.runs.columns, col = (k, i) => { const c = C[k]; return !Array.isArray(c) && c.dict ? c.dict[c.codes[i]] : c[i]; };
const runs = Array.from({ length: ix.runs.n }, (_, i) => ({ id: col('id', i), fc: col('fc', i), started: col('started', i) })).sort((a, b) => String(a.started).localeCompare(String(b.started)));
const resolve = route => route.replace(/\{run:([A-Z0-9]+)(?::(\d+))?\}/g, (_, fc, nth) => {
  const mine = runs.filter(r => r.fc === fc); const r = mine[(+nth || 1) - 1];
  if (!r) throw new Error(`the hub holds no run ${nth || 1} of flow cell ${fc}`);
  return r.id;
}).replace(/\{analysis:([^:}]+):([^}]+)\}/g, (_, label, pipeline) => {
  const a = (ix.analyses || []).find(x => x.label === label && String(x.pipeline).replace(/[^a-z]/g, '') === pipeline.replace(/[^a-z]/g, ''));
  if (!a) throw new Error(`the hub holds no ${pipeline} analysis of ${label}`);
  return a.id;
});

const spec = JSON.parse(fs.readFileSync(specPath, 'utf8'));
fs.mkdirSync(outDir, { recursive: true });
const browser = await pw[engine].launch();
const report = { browser: engine + ' ' + browser.version(), shots: [], errors: [] };
for (const job of spec) {
  const ctx = await browser.newContext({ viewport: { width: job.width || 1120, height: job.height || 900 }, deviceScaleFactor: job.scale || 1,
    colorScheme: job.scheme || 'light', timezoneId: 'UTC', locale: 'en-GB' });
  const page = await ctx.newPage();
  page.on('pageerror', e => report.errors.push(String(e.message || e).slice(0, 300)));
  const route = resolve(job.route);
  await page.goto(pathToFileURL(path.join(hubDir, 'index.html')).href + route);
  await page.waitForSelector('.fc-appbar', { timeout: 30000 });
  await page.waitForTimeout(job.wait || 1500);
  if (job.scheme === 'dark') {                                   // the page follows the system by default; press Dark so that it also says so
    const dark = page.locator('.fc-appbar [role="radiogroup"] button:has-text("Dark")');
    if (await dark.count()) { await dark.first().click(); await page.waitForTimeout(300); }
  }
  for (const sel of job.click || []) { await page.locator(sel).first().click(); await page.waitForTimeout(500); }
  for (const sel of job.hide || []) await page.addStyleTag({ content: `${sel} { display: none !important; }` });
  await page.mouse.move(0, 0);
  for (const shot of job.shots) {
    const file = path.join(outDir, shot.name + '.png');
    try {
      if (shot.pdf) {
        // what "Save as PDF" gives: the browser's own PDF of the page, on A4 with the margins the page asks for (Chromium only)
        if (engine !== 'chromium') throw new Error('a PDF is made with Chromium only');
        await page.pdf({ path: path.join(outDir, shot.name + '.pdf'), format: 'A4', margin: { top: '10mm', bottom: '10mm', left: '10mm', right: '10mm' } });
        await page.waitForTimeout(300);
        report.shots.push({ name: shot.name, route, file: path.join(outDir, shot.name + '.pdf') });
        continue;
      }
      if (shot.selector) {
        const el = page.locator(shot.selector).nth(shot.nth || 0);
        await el.scrollIntoViewIfNeeded();
        await page.waitForTimeout(250);
        await el.screenshot({ path: file });
      } else {
        // a region: the full width of the page, from the top (or the top of "from") to the bottom of "to"
        // ("tight": only as wide as "from", for a part of a page without the list of sections beside it)
        const box = await page.evaluate(([from, to]) => {
          const y = s => { const e = s && document.querySelector(s); if (!e) return null; const b = e.getBoundingClientRect(); return [b.top + window.scrollY, b.bottom + window.scrollY, b.left, b.right]; };
          const a = y(from), z = y(to);
          return { top: a ? a[0] : 0, bottom: z ? z[1] : document.documentElement.scrollHeight, width: document.documentElement.clientWidth,
            left: a ? a[2] : 0, right: a ? a[3] : document.documentElement.clientWidth, missing: (from && !a) || (to && !z) };
        }, [shot.from || null, shot.to || null]);
        if (box.missing) throw new Error(`nothing matches ${shot.from || ''} ${shot.to || ''}`);
        const pad = shot.pad == null ? 12 : shot.pad;
        const top = Math.max(0, box.top - (shot.from ? pad : 0));
        const x = shot.tight ? Math.max(0, box.left - pad) : 0, width = shot.tight ? Math.min(box.width - x, box.right - box.left + 2 * pad) : box.width;
        await page.screenshot({ path: file, fullPage: true, clip: { x, y: top, width, height: Math.min(box.bottom + pad - top, shot.max || 4000) } });
      }
      report.shots.push({ name: shot.name, route, file });
    } catch (e) { report.errors.push(`${shot.name}: ${String(e.message || e).slice(0, 300)}`); }
  }
  await ctx.close();
}
await browser.close();
console.log(JSON.stringify(report));
process.exit(report.errors.length ? 1 : 0);
