// State: signals and nothing else; all view state lives in the address fragment.
import { signal, computed } from '@preact/signals';
import { useEffect } from 'preact/hooks';

// ---- routing: fragment only ("#/run/TST00001_0a1b2c3d?t=0-24&by=bases") ----
export let HOME = '/home';                    // a one-run report sets this to its run (setHome)
export function setHome(path) { HOME = path; route.value = parseHash(location.hash); }
export function parseHash(hash) {
  const h = (hash || '').replace(/^#/, ''), i = h.indexOf('?'), path = (i < 0 ? h : h.slice(0, i)) || HOME, q = {};
  if (i >= 0) for (const kv of h.slice(i + 1).split('&')) {
    if (!kv) continue; const j = kv.indexOf('='), k = j < 0 ? kv : kv.slice(0, j);
    let v = j < 0 ? '' : kv.slice(j + 1).replace(/\+/g, '%20'); try { v = decodeURIComponent(v); } catch (e) { /* keep raw */ }
    q[k] = v;
  }
  const seg = path.split('/').filter(Boolean); let id = null;
  if (seg[1]) { try { id = decodeURIComponent(seg[1]); } catch (e) { id = seg[1]; } }
  return { screen: seg[0] || 'home', id, q, path };
}
const SAFE = /%(3A|2C|2F|3C|3E|22)/g;       // keep : , / < > " readable; a space is +, a literal + stays %2B
export function makeHash(path, q) {
  const parts = [];
  for (const k of Object.keys(q || {})) {
    const v = q[k]; if (v == null || v === '' || v === false) continue;
    parts.push(k + '=' + encodeURIComponent(String(v)).replace(SAFE, m => decodeURIComponent(m)).replace(/%20/g, '+'));
  }
  return '#' + path + (parts.length ? '?' + parts.join('&') : '');
}
export const route = signal(parseHash(location.hash));
addEventListener('hashchange', () => { route.value = parseHash(location.hash); });
export function go(path, q) { const h = makeHash(path, q); if (h === location.hash) return; location.hash = h; }       // a new history entry
export function setParams(patch, push) {                                                                              // filters and controls replace the entry
  const r = route.value, q = { ...r.q, ...patch }, h = makeHash(r.path, q);
  if (push) { location.hash = h; return; }
  try { history.replaceState(history.state, '', h); } catch (e) { /* file:// quirks: fall through to the signal */ }
  route.value = parseHash(h);
}
export const param = (k, d) => { const v = route.value.q[k]; return v == null || v === '' ? d : v; };
// Keep the keyboard where it was: after a control re-draws the page (a filter, a drawer), focus the element that
// matches `selector` once the page has been drawn again. Without this, focus falls back to the page body.
export function refocus(selector, fallback) {
  const run = () => { const el = document.querySelector(selector) || (fallback ? document.querySelector(fallback) : null); if (el && el !== document.activeElement) { try { el.focus({ preventScroll: false }); } catch (e) { /* not focusable any more */ } } };
  if (typeof requestAnimationFrame === 'function') requestAnimationFrame(() => requestAnimationFrame(run)); else setTimeout(run, 30);
}
export const hrefRun = id => '#/run/' + encodeURIComponent(id);
export const hrefProject = id => '#/project/' + encodeURIComponent(id);
export const hrefFlowCell = id => '#/flowcell/' + encodeURIComponent(id);

// ---- theme: Auto, Light, Dark on <html data-theme>; remembered per browser and per library as a convenience ----
let prefix = 'fc:_:';
function stored(k) { try { return localStorage.getItem(prefix + k); } catch (e) { return null; } }
function store(k, v) { try { localStorage.setItem(prefix + k, v); } catch (e) { /* private window: fine */ } }
export const themeChoice = signal('auto');
export const themeTick = signal(0);          // charts re-read the tokens when this changes
export function applyTheme(choice, remember = true) {
  themeChoice.value = choice; if (remember) store('theme', choice);
  const el = document.documentElement;
  if (choice === 'auto') el.removeAttribute('data-theme'); else el.setAttribute('data-theme', choice);
  themeTick.value++;
}
// storage keys are namespaced "fc:<library id>:" once the library is known
export function useLibraryStorage(libraryId) {
  prefix = 'fc:' + (libraryId || '_') + ':';
  const t = stored('theme'); if (t === 'light' || t === 'dark' || t === 'auto') { if (t !== themeChoice.value) applyTheme(t, false); }
}
try { matchMedia('(prefers-color-scheme: dark)').addEventListener('change', () => { themeTick.value++; }); } catch (e) { /* old engines */ }
export const isDark = computed(() => { themeTick.value; const c = themeChoice.value; return c === 'dark' || (c === 'auto' && typeof matchMedia === 'function' && matchMedia('(prefers-color-scheme: dark)').matches); });

// ---- text size: 1 at the designed size, 2 when the browser shows text twice as large (text-only zoom, a
// larger default font). Measured by a probe the application draws once (shell/parts.jsx TextProbe); layouts whose
// widths are stated in pixels for text of the designed size multiply them by it.
export const textScale = signal(1);

// ---- scroll memory: Back returns to where one was. The position of a scrolling container (and of the window)
// is kept in the history entry of the page, so a new visit starts at the top and Back does not.
function keep(patch) { try { history.replaceState({ ...(history.state || {}), ...patch }, ''); } catch (e) { /* file:// quirks */ } }
export function useScrollMemory(ref) {
  useEffect(() => {
    const el = ref.current, st = history.state || {};
    const apply = () => { if (el && st.fcWrap) el.scrollTop = st.fcWrap; if (st.fcWin) window.scrollTo(0, st.fcWin); };
    apply(); const again = setTimeout(apply, 90);         // once more when the rows have their height
    let timer = 0;
    const now = () => keep({ fcWrap: el ? el.scrollTop : 0, fcWin: window.scrollY });
    const later = () => { clearTimeout(timer); timer = setTimeout(now, 120); };
    if (el) { el.addEventListener('scroll', later, { passive: true }); el.addEventListener('click', now, true); }
    addEventListener('scroll', later, { passive: true });
    return () => { clearTimeout(again); clearTimeout(timer); if (el) { el.removeEventListener('scroll', later); el.removeEventListener('click', now, true); } removeEventListener('scroll', later); };
  }, []);
}

// ---- compare selection (ticked rows of the Runs table; kept across screens) ----
export const selected = signal([]);
export function toggleSelected(id) { const s = selected.value; selected.value = s.includes(id) ? s.filter(x => x !== id) : s.concat(id); }
