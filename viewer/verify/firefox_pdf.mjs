// A real PDF from Firefox, made by Firefox's own print pipeline (the one behind its print dialog).
//
// Playwright can ask only Chromium for a PDF. Firefox is therefore driven here over Marionette, its WebDriver back
// end: a headless Firefox with a throw-away profile, a few commands over one local socket, and "WebDriver:Print",
// which fires `beforeprint` and `afterprint` on the page like the dialog does. Nothing is installed or downloaded:
// the Firefox is the Playwright build the other checks use.
//
//   const ff = await firefoxSession();            // null when no Firefox is found
//   await ff.open(url); await ff.js('return document.title'); await ff.pdf('out.pdf', 'a4'); await ff.close();
//
//   node verify/firefox_pdf.mjs URL OUT.pdf [a4|letter]
import fs from 'node:fs';
import net from 'node:net';
import os from 'node:os';
import path from 'node:path';
import { spawn } from 'node:child_process';
import { fileURLToPath, pathToFileURL } from 'node:url';

const here = path.dirname(fileURLToPath(import.meta.url)), viewerDir = path.resolve(here, '..');
export const PAPERS = { a4: { width: 21.0, height: 29.7 }, letter: { width: 21.59, height: 27.94 } };      // centimetres
const sleep = ms => new Promise(r => setTimeout(r, ms));

export function firefoxBinary() {
  if (process.env.FLOWCELLAR_FIREFOX && fs.existsSync(process.env.FLOWCELLAR_FIREFOX)) return process.env.FLOWCELLAR_FIREFOX;
  for (const dir of [process.env.PLAYWRIGHT_BROWSERS_PATH, path.join(viewerDir, 'pw-browsers'), process.env.FLOWCELLAR_PW_BROWSERS].filter(Boolean)) {
    if (!fs.existsSync(dir)) continue;
    for (const n of fs.readdirSync(dir).filter(n => n.startsWith('firefox-')).sort().reverse()) { const bin = path.join(dir, n, 'firefox', 'firefox'); if (fs.existsSync(bin)) return bin; }
  }
  return null;
}

export async function firefoxSession({ width = 1366, height = 900 } = {}) {
  const bin = firefoxBinary(); if (!bin) return null;
  const port = 21000 + Math.floor(Math.random() * 20000), profile = fs.mkdtempSync(path.join(os.tmpdir(), 'fcff-'));
  const prefs = { 'marionette.port': port, 'browser.shell.checkDefaultBrowser': false, 'datareporting.policy.dataSubmissionEnabled': false, 'toolkit.telemetry.enabled': false, 'app.update.enabled': false,
    'browser.aboutwelcome.enabled': false, 'extensions.update.enabled': false, 'app.normandy.enabled': false, 'network.captive-portal-service.enabled': false, 'network.connectivity-service.enabled': false,
    'browser.startup.homepage_override.mstone': 'ignore', 'startup.homepage_welcome_url': '', 'browser.safebrowsing.malware.enabled': false, 'browser.safebrowsing.phishing.enabled': false, 'network.prefetch-next': false };
  fs.writeFileSync(path.join(profile, 'user.js'), Object.entries(prefs).map(([k, v]) => 'user_pref(' + JSON.stringify(k) + ', ' + JSON.stringify(v) + ');').join('\n') + '\n');
  const proc = spawn(bin, ['--headless', '--marionette', '--no-remote', '--profile', profile, 'about:blank'], { stdio: ['ignore', 'ignore', 'pipe'] });
  let said = ''; proc.stderr.on('data', d => { said = (said + d).slice(-2000); });
  const close = async () => { try { proc.kill('SIGKILL'); } catch (e) { /* gone */ } await sleep(250); fs.rmSync(profile, { recursive: true, force: true }); };
  let sock = null;
  for (let i = 0; i < 160 && !sock; i++) { sock = await new Promise(res => { const s = net.connect(port, '127.0.0.1'); s.once('connect', () => res(s)); s.once('error', () => res(null)); }); if (!sock) await sleep(250); }
  if (!sock) { await close(); throw new Error('Firefox did not open its Marionette port: ' + said.slice(-300)); }
  // messages are "<length>:<json>"; a command is [0, id, name, params] and its answer [1, id, error, result]
  let buf = Buffer.alloc(0), id = 0; const waiting = new Map();
  sock.on('data', d => { buf = Buffer.concat([buf, d]);
    for (;;) { const i = buf.indexOf(58); if (i < 0) return; const n = +buf.subarray(0, i).toString(); if (buf.length < i + 1 + n) return;
      const msg = JSON.parse(buf.subarray(i + 1, i + 1 + n).toString('utf8')); buf = buf.subarray(i + 1 + n);
      if (!Array.isArray(msg)) continue; const w = waiting.get(msg[1]); if (!w) continue; waiting.delete(msg[1]); if (msg[2]) w.rej(new Error(String(msg[2].message || JSON.stringify(msg[2])).slice(0, 500))); else w.res(msg[3]); } });
  const send = (name, params = {}) => new Promise((res, rej) => { const m = JSON.stringify([0, ++id, name, params]); waiting.set(id, { res, rej }); sock.write(Buffer.byteLength(m) + ':' + m); });
  try { await send('WebDriver:NewSession', { capabilities: { alwaysMatch: {} } }); await send('WebDriver:SetWindowRect', { width, height }).catch(() => {}); } catch (e) { await close(); throw e; }
  const js = (script, args = []) => send('WebDriver:ExecuteScript', { script, args, newSandbox: false }).then(r => r.value);
  return {
    js,
    version: () => js('return navigator.userAgent').then(ua => 'firefox ' + ((/Firefox\/([\d.]+)/.exec(ua) || [])[1] || '?')),
    // open a page of the viewer and wait until `ready` (a JavaScript expression) holds
    async open(url, ready = '!!document.querySelector(".fc-appbar")', settle = 1200) {
      await send('WebDriver:Navigate', { url: /^[a-z]+:/.test(url) ? url : pathToFileURL(path.resolve(url)).href });
      let ok = false; for (let i = 0; i < 120 && !ok; i++) { ok = await js('return !!(' + ready + ')').catch(() => false); if (!ok) await sleep(250); }
      await sleep(settle); return ok;
    },
    async pdf(out, paper = 'a4') {
      const r = await send('WebDriver:Print', { page: PAPERS[paper], margin: { top: 1, bottom: 1, left: 1, right: 1 }, scale: 1, shrinkToFit: true, background: false, orientation: 'portrait' });
      fs.writeFileSync(out, Buffer.from(r.value, 'base64')); return out;
    },
    async close() { await send('Marionette:Quit', { flags: ['eForceQuit'] }).catch(() => {}); await sleep(200); sock.destroy(); await close(); },
  };
}

if (process.argv[1] && path.resolve(process.argv[1]) === fileURLToPath(import.meta.url)) {
  const [url, out, paper = 'a4'] = process.argv.slice(2);
  if (!url || !out || !PAPERS[paper]) { console.error('usage: node verify/firefox_pdf.mjs URL OUT.pdf [a4|letter]'); process.exit(2); }
  const ff = await firefoxSession(); if (!ff) { console.error('No Playwright Firefox found'); process.exit(2); }
  let code = 0;
  try { await ff.open(url); await ff.pdf(out, paper); console.log((await ff.version()) + ' wrote ' + out + ' (' + fs.statSync(out).size + ' bytes)'); } catch (e) { console.error(e.message); code = 1; }
  await ff.close(); process.exit(code);
}
