// Data. The boot block of the page says how it gets its data (two ways, one start-up path):
//   report  the index and one run are inside the page, in script elements the browser does not run;
//   hub     the page holds no data: "index.js" and "data/<run key>.js" beside it are loaded as scripts, the
//           index at start-up and each run when a screen first needs it. The same files work from file://
//           and over "flowcellar serve".
// A payload is base64 of gzip bytes. Nothing is ever fetched: no fetch, no XHR, no import(), no worker.
import { signal } from '@preact/signals';

export const RUN_KEY = /^[A-Z]{2,3}[0-9]{4,6}_[0-9a-f]{8,12}$/;
// a single-cell or spatial analysis: a bundle of its own kind, loaded the same way (data/<analysis id>.js)
export const ANALYSIS_KEY = /^an_[a-z]{2,16}_[0-9a-f]{12}$/;
function readBoot() {
  try { const b = JSON.parse(document.getElementById('fc-boot').textContent); if (b && (b.mode === 'report' || b.mode === 'hub')) return b; } catch (e) { /* no boot block: decide from the page */ }
  return { mode: document.getElementById('fci') ? 'report' : 'hub', generator: null, schema: null, run: null };
}
export const boot = readBoot();
export const isReport = boot.mode === 'report';
export const canInflate = typeof DecompressionStream === 'function' && typeof Response === 'function';

function bytesOf(b64) { const bin = atob(b64), out = new Uint8Array(bin.length); for (let i = 0; i < bin.length; i++) out[i] = bin.charCodeAt(i); return out; }
function inflate(bytes) { return new Response(new Blob([bytes]).stream().pipeThrough(new DecompressionStream('gzip'))).json(); }
async function sha256(bytes) {
  if (!(globalThis.crypto && crypto.subtle)) return null;
  const d = new Uint8Array(await crypto.subtle.digest('SHA-256', bytes)); let s = '';
  for (const v of d) s += v.toString(16).padStart(2, '0'); return s;
}
// decode one payload; a checksum that does not match is an error named "damaged"
async function decode(b64, sha) {
  const bytes = bytesOf(String(b64).trim()), data = await inflate(bytes);
  if (sha) { const got = await sha256(bytes); if (got && got !== sha) { const e = new Error('the checksum of the data does not match'); e.damaged = true; throw e; } }
  return data;
}

// ---- sibling files, loaded as scripts (hub) ----
const MAX_IN_FLIGHT = 4;
let inFlight = 0, generation = 0; const queue = [], waiting = new Map();
function inject(src, done) {
  const s = document.createElement('script'); s.async = true; s.src = src;
  const end = ok => { s.onload = s.onerror = null; s.remove(); done(ok); };
  s.onload = () => end(true); s.onerror = () => end(false);
  document.head.appendChild(s);
}
function pump() {
  while (inFlight < MAX_IN_FLIGHT && queue.length) {
    const job = queue.shift(); inFlight++;
    inject(job.src, ok => { inFlight--; job.done(ok); pump(); });
  }
}
function sidecar(id) {
  return new Promise((resolve, reject) => {
    if (!RUN_KEY.test(id) && !ANALYSIS_KEY.test(id)) return reject(new Error('not a run key'));
    const w = { payload: null }; waiting.set(id, w);
    queue.push({ src: 'data/' + id + '.js' + (generation ? '?t=' + generation : ''), done: ok => {
      waiting.delete(id);
      if (!ok || !w.payload) return reject(new Error(ok ? 'the data file does not hold this run' : 'the data file of this run was not found'));
      resolve(w.payload);
    } });
    pump();
  });
}
let indexWait = null;
if (!isReport) {
  // defined before any sibling file loads
  window.__FCB_SIDECAR__ = (id, b64, sha) => { const w = waiting.get(String(id)); if (w) w.payload = { b64: String(b64), sha: sha ? String(sha) : null }; };
  window.__FC_INDEX__ = b64 => { if (indexWait) indexWait.payload = String(b64); };
}
function indexSidecar() {
  return new Promise((resolve, reject) => {
    const w = indexWait = { payload: null };
    inject('index.js' + (generation ? '?t=' + generation : ''), ok => {
      if (indexWait === w) indexWait = null;
      if (!ok || !w.payload) { const e = new Error('no index'); e.missing = true; return reject(e); }
      resolve(w.payload);
    });
  });
}

// ---- runs: decoded on demand, the most recently used kept ----
export const bundleTick = signal(0);        // bumps when a run finishes decoding
const cache = new Map(), KEEP = 24; let clock = 0;
function evict() {
  if (cache.size <= KEEP) return;
  const old = [...cache.entries()].filter(([, e]) => e.state !== 'loading' && clock - e.used > 40).sort((a, b) => a[1].used - b[1].used);
  for (const [id] of old) { if (cache.size <= KEEP) break; cache.delete(id); }
}
export function loadBundle(id) {
  let e = cache.get(id);
  if (e) { e.used = ++clock; return e.promise; }
  e = { state: 'loading', data: null, error: null, damaged: false, used: ++clock };
  let source;
  if (isReport) { const el = document.getElementById('fcb:' + id); source = el ? Promise.resolve({ b64: el.textContent, sha: el.getAttribute('data-sha256') }) : Promise.reject(new Error('this run is not in the file')); }
  else { // a key the index does not hold (a mistyped address) has no data file: none is asked for, so nothing fails in the console
    const ix = index.value, known = !ix || ix.rowOf.get(id) != null || (ix.analysisById && ix.analysisById.has(id));
    source = known ? sidecar(id) : Promise.reject(new Error('this run is not in the library')); }
  e.promise = source.then(p => decode(p.b64, p.sha)).then(
    d => { e.state = 'ready'; e.data = d; bundleTick.value++; evict(); return d; },
    err => { e.state = 'failed'; e.error = String(err && err.message || err); e.damaged = !!(err && err.damaged); bundleTick.value++; return null; });
  cache.set(id, e);
  return e.promise;
}
// Synchronous read for components: the run when decoded, null while loading (and starts the load).
export function useBundle(id) {
  bundleTick.value;                         // subscribe
  if (!id) return null;
  const e = cache.get(id);
  if (!e) { loadBundle(id); return null; }
  e.used = ++clock;
  return e.state === 'ready' ? e.data : null;
}
// the run when it is already decoded; never starts a load
export function cachedBundle(id) { bundleTick.value; const e = cache.get(id); return e && e.state === 'ready' ? e.data : null; }
export function bundleState(id) { bundleTick.value; const e = cache.get(id); return e ? e.state : 'absent'; }
export function bundleError(id) { bundleTick.value; const e = cache.get(id); return e && e.state === 'failed' ? { message: e.error, damaged: e.damaged } : null; }
export function retryBundle(id) { const e = cache.get(id); if (e && e.state === 'failed') { cache.delete(id); generation++; bundleTick.value++; } }

// ---- the index: columnar -> rows and look-up tables ----
export const index = signal(null);
export const indexState = signal({ state: 'loading', error: null });        // loading | ready | missing | failed | old_browser
const column = c => (c && !Array.isArray(c) && c.dict ? i => c.dict[c.codes[i]] : i => c[i]);
const arr = v => Array.isArray(v) ? v : [];
function hydrate(ix) {
  const C = (ix.runs && ix.runs.columns) || {}, n = (ix.runs && ix.runs.n) || 0, names = Object.keys(C), get = names.map(k => column(C[k]));
  const rows = [], rowOf = new Map();
  for (let i = 0; i < n; i++) {
    const r = { i }; for (let k = 0; k < names.length; k++) r[names[k]] = get[k](i);
    r.startedMs = Date.parse(r.started) || 0;
    r.stub = r.detail_t1 == null || r.detail_t1 === 'none';
    r.health = r.health || 'not_available'; r.data = r.data || 'not_available';
    r.modifiers = arr(r.modifiers); r.reasons = arr(r.reasons); r.reason_levels = arr(r.reason_levels); r.reason_checks = arr(r.reason_checks); r.reason_tiers = arr(r.reason_tiers);
    r.roots = arr(r.roots); r.project_ids = arr(r.project_ids); r.read_data = r.read_data || 'absent';
    r.sig = [r.verdict, r.detail_t1, r.read_data, r.checks_done, r.bundle_gzip_bytes, r.headline, r.t2_reads].join('|');
    rows.push(r); rowOf.set(r.id, i);
  }
  const samplesByRun = new Map();
  for (const s of arr(ix.samples)) { if (!samplesByRun.has(s.run)) samplesByRun.set(s.run, []); samplesByRun.get(s.run).push(s); }
  const groups = new Map(); for (const g of arr(ix.baselines && ix.baselines.groups)) groups.set(g.group, g);
  const byFc = new Map(); for (const r of rows) { if (!byFc.has(r.fc)) byFc.set(r.fc, []); byFc.get(r.fc).push(r); }
  for (const l of byFc.values()) l.sort((a, b) => a.startedMs - b.startedMs);
  const projects = arr(ix.projects), projectById = new Map(); for (const p of projects) projectById.set(p.id, p);
  const roots = arr(ix.roots), rootById = new Map(); for (const r of roots) rootById.set(r.id, r);
  const lib = ix.library || {};
  // single-cell and spatial analyses: each link is joined to the runs of this index here, so a link resolves by itself once its run is in the library
  const analyses = arr(ix.analyses), analysisById = new Map(), analysesByRun = new Map();
  for (const a of analyses) { analysisById.set(a.id, a); a.links = arr(a.links); a.siblings = arr(a.siblings); a.reasons = arr(a.reasons); a.reason_levels = arr(a.reason_levels);
    for (const l of a.links) { l.resolved = l.run != null && rowOf.has(l.run); if (l.resolved) { if (!analysesByRun.has(l.run)) analysesByRun.set(l.run, []); analysesByRun.get(l.run).push(a); } } }
  return { analyses, analysisById, analysesByRun, raw: ix, n, rows, rowOf, samplesByRun, groups, byFc, written: ix.written, notices: arr(ix.notices), library: lib, projects, projectById, roots, rootById,
    activity: ix.activity || { months: [], runs: [], estimated_bases: [], reads: [] }, progress: ix.progress || { state: 'idle' } };
}
function install(ix) {
  const next = hydrate(ix), prev = index.value;
  if (prev) for (const r of next.rows) { const o = prev.rows[prev.rowOf.get(r.id)]; if (o && o.sig !== r.sig) cache.delete(r.id); }     // a run that changed is read again
  if (prev) for (const a of next.analyses) { const o = prev.analysisById.get(a.id); if (o && (o.bundle_gzip_bytes !== a.bundle_gzip_bytes || o.verdict !== a.verdict || o.headline !== a.headline)) cache.delete(a.id); }
  index.value = next; indexState.value = { state: 'ready', error: null };
  return next;
}
export async function loadIndex() {
  if (!canInflate) { indexState.value = { state: 'old_browser', error: null }; throw new Error('old browser'); }
  try {
    const el = isReport ? document.getElementById('fci') : null;
    const b64 = isReport ? (el ? el.textContent : null) : await indexSidecar();
    if (b64 == null) { const e = new Error('no index'); e.missing = true; throw e; }
    return install(await decode(b64));
  } catch (e) { indexState.value = { state: e.missing ? 'missing' : 'failed', error: String(e && e.message || e) }; throw e; }
}
// Read index.js again (hub only). Returns true when a newer index was installed.
export const reloading = signal(false);
export async function reloadIndex() {
  if (isReport || reloading.value) return false;
  reloading.value = true; generation++;
  try { const ix = await decode(await indexSidecar()); if (index.value && ix.written === index.value.written) return false; install(ix); bundleTick.value++; return true; }
  catch (e) { return false; }
  finally { reloading.value = false; }
}
// while a scan is filling the library, look for a newer index once a minute
if (!isReport) setInterval(() => { const ix = index.value; if (ix && ix.progress && ix.progress.state === 'scanning' && !document.hidden) reloadIndex(); }, 60000);
export const rowById = id => { const ix = index.value; const i = ix && ix.rowOf.get(id); return i == null ? null : ix.rows[i]; };
