// Open one page of the viewer in headless Firefox, print what went wrong, and write one screenshot.
//
//   node verify/smoke.mjs FILE_OR_URL [--hash '#/runs'] [--out shot.png] [--width 1366] [--height 768] [--dark] [--full]
//                         [--click SELECTOR] [--scroll SELECTOR] [--wait MS] [--text N]
import { launch, fileUrl, settle } from './browser.mjs';

const args = process.argv.slice(2), opt = n => { const i = args.indexOf('--' + n); return i >= 0 ? args[i + 1] : null; }, flag = n => args.includes('--' + n);
const target = args.find(a => !a.startsWith('--') && args[args.indexOf(a) - 1]?.startsWith('--') !== true) || args[0];
if (!target) { console.error('usage: node verify/smoke.mjs FILE_OR_URL [--hash #/runs] [--out shot.png] [--width 1366] [--dark] [--full]'); process.exit(2); }
const base = /^https?:/.test(target) ? target : fileUrl(target), origin = /^https?:/.test(base) ? new URL(base).origin : base.slice(0, base.lastIndexOf('/') + 1);
const { log, newPage, browser, version } = await launch();
const page = await newPage({ width: +(opt('width') || 1366), height: +(opt('height') || 768), scheme: flag('dark') ? 'dark' : 'light', allowed: u => u.startsWith(origin) });
await page.goto(base + (opt('hash') || ''));
await page.waitForSelector('.fc-appbar, .fc-boot--failed', { timeout: 20000 }).catch(() => {});
await settle(page, +(opt('wait') || 900));
if (opt('click')) { for (const sel of opt('click').split(';;')) { await page.locator(sel).first().click(); await settle(page, 400); } }
if (opt('scroll')) { await page.locator(opt('scroll')).first().evaluate(el => el.scrollIntoView({ block: 'start' })); await page.evaluate(() => window.scrollBy(0, -90)); await settle(page, 500); }
if (opt('out')) await page.screenshot({ path: opt('out'), fullPage: flag('full') });
const text = (await page.evaluate(() => document.body.innerText)).replace(/\s+/g, ' ').slice(0, +(opt('text') || 300));
await page.finish(); await browser.close();
console.log(version + ' · ' + (await Promise.resolve(text)));
console.log('errors: ' + log.errors.length + ' · requests to another origin: ' + log.requests.length + ' · browser notices: ' + log.notices);
for (const e of log.errors.slice(0, 12)) console.log('  ' + e.kind + ': ' + e.text);
for (const r of log.requests.slice(0, 12)) console.log('  request: ' + r);
process.exit(log.errors.length || log.requests.length ? 1 : 0);
