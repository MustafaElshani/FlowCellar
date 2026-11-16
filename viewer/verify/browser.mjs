// Shared by the browser checks: a headless browser (Playwright's own build), pages that record every console
// error, page error, policy violation and request, and small helpers.
//
// Which browser: Firefox, or Chromium when FLOWCELLAR_BROWSER_ENGINE=chromium (a script may also ask for one by name:
// launch('chromium')). The browser is used by path. It is looked for in PLAYWRIGHT_BROWSERS_PATH, then in
// viewer/pw-browsers, then in the folder the development machine keeps it in; nothing is downloaded by these scripts.
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath, pathToFileURL } from 'node:url';

export const here = path.dirname(fileURLToPath(import.meta.url));
export const viewerDir = path.resolve(here, '..');
export const repo = path.resolve(viewerDir, '..');
export const ENGINES = ['firefox', 'chromium'];
export const engine = (() => { const e = (process.env.FLOWCELLAR_BROWSER_ENGINE || 'firefox').toLowerCase(); if (!ENGINES.includes(e)) { console.error('FLOWCELLAR_BROWSER_ENGINE must be firefox or chromium, not ' + e); process.exit(2); } return e; })();
const candidates = [process.env.PLAYWRIGHT_BROWSERS_PATH, path.join(viewerDir, 'pw-browsers'), process.env.FLOWCELLAR_PW_BROWSERS].filter(Boolean);
// the folder that holds a Playwright build of the browser, or null
export const browsersDir = name => candidates.find(d => fs.existsSync(d) && fs.readdirSync(d).some(n => n.startsWith(name + '-'))) || null;
export const hasBrowser = name => !!browsersDir(name);
const Name = n => n[0].toUpperCase() + n.slice(1);
const found = browsersDir(engine);
if (!found) { console.error('No Playwright ' + Name(engine) + ' found. Set PLAYWRIGHT_BROWSERS_PATH (or FLOWCELLAR_PW_BROWSERS) to a folder that holds one, or run:\n  PLAYWRIGHT_BROWSERS_PATH=' + path.join(viewerDir, 'pw-browsers') + ' npx playwright install ' + engine); process.exit(2); }
process.env.PLAYWRIGHT_BROWSERS_PATH = found;
const pw = (await import('playwright')).default;

// Firefox prints this once per page load whenever DecompressionStream is used, on any page (also on an empty test
// page that only inflates gzip). It is counted apart and is not an error of the viewer.
const BROWSER_NOTICE = /The execution of microtasks on this page was reordered/;

export async function launch(name = engine) {
  if (!ENGINES.includes(name)) throw new Error('no such browser: ' + name);
  if (name !== engine && !hasBrowser(name)) throw new Error('No Playwright ' + Name(name) + ' found');
  if (name !== engine) process.env.PLAYWRIGHT_BROWSERS_PATH = browsersDir(name);
  const browser = await pw[name].launch();
  const log = { errors: [], notices: 0, requests: [], pages: 0 };
  // allowed(url): the page's own origin. Everything else that is asked for is recorded as a request to another origin.
  async function newPage({ width = 1366, height = 768, scheme = 'light', allowed = () => true, timezoneId } = {}) {
    const ctx = await browser.newContext({ viewport: { width, height }, colorScheme: scheme, deviceScaleFactor: 1, ...(timezoneId ? { timezoneId } : {}) });
    const page = await ctx.newPage(); log.pages++;
    page.on('console', m => { if (m.type() !== 'error' && m.type() !== 'warning') return; const text = m.text(); if (BROWSER_NOTICE.test(text)) { log.notices++; return; } log.errors.push({ kind: 'console.' + m.type(), text: text.slice(0, 500), url: page.url() }); });
    page.on('pageerror', e => log.errors.push({ kind: 'pageerror', text: String(e.stack || e.message).slice(0, 800), url: page.url() }));
    page.on('request', r => { const u = r.url(); if (u.startsWith('data:') || u.startsWith('blob:') || u.startsWith('about:')) return; if (!allowed(u)) log.requests.push(u); });
    await page.addInitScript(() => { window.__csp = []; document.addEventListener('securitypolicyviolation', e => window.__csp.push(e.violatedDirective + ' ' + (e.blockedURI || '') + ' ' + (e.sourceFile || '') + ':' + e.lineNumber)); });
    page.finish = async () => { const v = await page.evaluate(() => window.__csp || []).catch(() => []); for (const x of v) log.errors.push({ kind: 'csp', text: x, url: page.url() }); await ctx.close(); };
    return page;
  }
  return { browser, log, newPage, name, version: name + ' ' + browser.version() };
}
export const fileUrl = p => pathToFileURL(path.resolve(p)).href;
export async function settle(page, ms = 300) { await page.waitForTimeout(ms); await page.evaluate(() => new Promise(r => requestAnimationFrame(() => requestAnimationFrame(r)))); }
// open a page of the viewer and wait until it has drawn its application bar
export async function open(page, url) { await page.goto(url); await page.waitForSelector('.fc-appbar', { timeout: 20000 }); await settle(page); }
// choose a theme: by its button of the three, or (a window of 1,000 px or less shows one button that steps through
// the three choices) by pressing that button until it names the theme
export async function setTheme(page, name) {
  const seg = page.locator(`.fc-appbar [role="radiogroup"] button:has-text("${name}")`);
  if (await seg.isVisible()) await seg.click();
  else { const one = page.locator('.fc-appbar .fc-theme1'); for (let k = 0; k < 3 && (await one.innerText()).trim() !== name; k++) { await one.click(); await page.waitForTimeout(80); } }
  await settle(page, 200);
}
// true when the page scrolls sideways (something is wider than the window)
export const sideways = page => page.evaluate(() => document.documentElement.scrollWidth - document.documentElement.clientWidth);
