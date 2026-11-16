// FlowCellar viewer: start-up, routing, shell. One script for the one-run report and for the hub of a library;
// the boot block of the page says which it is (core/data.js).
import { render } from 'preact';
import { effect } from '@preact/signals';
import { loadIndex, index, indexState, rowById, boot, isReport, canInflate } from './core/data.js';
import { route, setParams, setHome, themeTick, themeChoice, applyTheme, useLibraryStorage, param, textScale, makeHash, parseHash } from './core/state.js';
import { useEffect, useRef } from 'preact/hooks';
import { resetTheme } from './charts/theme.js';
import { printing } from './core/print.js';
import { AppBar } from './shell/parts.jsx';
import { Home } from './home/Home.jsx';
import { Projects, Project } from './projects/Projects.jsx';
import { Hub } from './hub/Hub.jsx';
import { Run } from './run/Run.jsx';
import { Compare } from './compare/Compare.jsx';
import { FlowCell } from './flowcell/FlowCell.jsx';
import { Analyses } from './analysis/Analyses.jsx';
import { Analysis } from './analysis/Analysis.jsx';
import { VERDICT_WORD } from './core/labels.js';

// charts read the tokens again whenever the theme changes, and when the page goes to paper (the theme "print") or comes back
effect(() => { themeTick.value; printing.value; resetTheme(); });

const SCREENS = ['home', 'runs', 'projects', 'project', 'run', 'compare', 'flowcell', 'analyses', 'analysis'];
function NotFound({ r }) {
  return <div class="fc-page"><div class="fc-empty fc-page__empty"><h1 class="fc-empty__title">This address is not a page of the library.</h1><span><span class="fc-mono">#{r.path}</span> names no screen. <a href="#/home">Home</a> {'·'} <a href="#/runs">Runs</a> {'·'} <a href="#/projects">Projects</a></span></div></div>;
}
// The size text is really drawn at, against the size it was designed at. A hidden box 10em wide at a font size of
// 10px is 100px wide; with text-only zoom, or a larger minimum font size in the browser, it is wider.
function TextProbe() {
  const ref = useRef(null);
  useEffect(() => {
    const el = ref.current; if (!el || typeof ResizeObserver !== 'function') return;
    const read = () => { const w = el.getBoundingClientRect().width, s = w > 0 ? Math.max(1, Math.min(3, Math.round(w / 100 * 20) / 20)) : 1; if (s !== textScale.value) textScale.value = s; };
    const ro = new ResizeObserver(read); ro.observe(el); read();
    return () => ro.disconnect();
  }, []);
  return <span ref={ref} class="fc-textprobe" aria-hidden="true" />;
}
// "#/run/TST00001": a flow cell id where a run key is expected. One run of that flow cell: that run; several: its history.
function resolveRun(id, ix) {
  if (!id || ix.rowOf.get(id) != null) return null;
  const list = ix.byFc.get(id); if (!list || !list.length) return null;
  return list.length === 1 ? '/run/' + encodeURIComponent(list[0].id) : '/flowcell/' + encodeURIComponent(id);
}

function App() {
  const r = route.value, ix = index.value;
  // Every screen draws itself again when the theme changes, and on the way to paper and back: a screen reads colours from
  // the theme while it is drawn (the series of its charts, its legends). A screen that reads signals is otherwise drawn
  // again only when one of its own signals or props changes, so the change is handed to it as a prop.
  const paint = themeTick.value + (printing.value ? ':paper' : '');
  if (!ix) return null;
  let screen, runs = false;
  if (isReport) screen = boot.analysis ? <Analysis key={boot.analysis} id={boot.analysis} paint={paint} /> : <Run key={boot.run} id={boot.run} paint={paint} />;        // a report is its one run (or its one analysis), whatever the address says
  else switch (r.screen) {
    case 'home': screen = <Home paint={paint} />; break;
    case 'runs': screen = <Hub paint={paint} />; runs = true; break;
    case 'projects': screen = <Projects paint={paint} />; break;
    case 'project': screen = r.id ? <Project key={r.id} paint={paint} /> : <Projects paint={paint} />; break;
    case 'run': screen = r.id ? <Run key={r.id} id={r.id} paint={paint} /> : <Hub paint={paint} />; break;
    case 'compare': screen = <Compare paint={paint} />; break;
    case 'flowcell': screen = r.id ? <FlowCell key={r.id} id={r.id} paint={paint} /> : <Hub paint={paint} />; break;
    case 'analyses': screen = <Analyses paint={paint} />; break;
    case 'analysis': screen = r.id ? <Analysis key={r.id} id={r.id} paint={paint} /> : <Analyses paint={paint} />; break;
    default: screen = <NotFound r={r} />;
  }
  return <>
    {/* the address of a page is its hash: the skip link must not change it, so it moves the focus itself */}
    <a class="fc-skip" href="#main" onClick={e => { e.preventDefault(); const el = document.getElementById('main'); if (el) { el.focus({ preventScroll: true }); el.scrollIntoView({ block: 'start' }); } }}>Skip to the content</a>
    <AppBar q={runs ? param('q', '') : undefined} onQuery={runs ? v => setParams({ q: v }) : undefined} />
    <TextProbe />
    <div id="main" tabIndex={-1}>{screen}</div>
  </>;
}
// a flow cell id where a run key belongs is led to its run, or to the history of the flow cell
effect(() => {
  const r = route.value, ix = index.value; if (!ix || isReport || r.screen !== 'run') return;
  const to = resolveRun(r.id, ix); if (!to) return;
  const h = makeHash(to, r.q); try { history.replaceState(history.state, '', h); } catch (e) { /* file:// quirks */ }
  route.value = parseHash(h);
});
// the title of the tab says where one is; a report keeps naming its run and verdict
effect(() => {
  const r = route.value, ix = index.value; if (!ix) return;
  const id = isReport ? boot.run : r.screen === 'run' ? r.id : null, row = id ? rowById(id) : null, lib = ix.library && ix.library.name;
  if (isReport && boot.analysis) { const a = ix.analysisById.get(boot.analysis); document.title = (a ? (a.label || a.id) + ' · ' + (VERDICT_WORD[a.verdict] || '') + ' · ' : '') + 'FlowCellar analysis report'; return; }
  if (!isReport && (r.screen === 'analyses' || r.screen === 'analysis')) { const a = r.id ? ix.analysisById.get(r.id) : null; document.title = (a ? (a.label || a.id) + ' · ' + (VERDICT_WORD[a.verdict] || '') : r.id ? 'Not found' : 'Analyses') + ' · ' + (lib ? lib + ' · ' : '') + 'FlowCellar'; return; }
  if (isReport) { document.title = (row ? row.fc + ' · ' + (VERDICT_WORD[row.verdict] || '') + (row.provisional && (row.verdict === 'pass' || row.verdict === 'warn') ? ' so far' : '') + ' · ' : '') + 'FlowCellar run report'; return; }
  const p = r.screen === 'project' && r.id ? ix.projectById.get(r.id) : null;
  const where = row ? row.fc + ' · ' + (VERDICT_WORD[row.verdict] || '') + (row.provisional ? ' so far' : '') : r.screen === 'compare' ? 'Compare' : r.screen === 'flowcell' && r.id ? r.id + ' · flow cell history'
    : r.screen === 'runs' ? 'Runs' : r.screen === 'projects' ? 'Projects' : p ? (p.name == null ? 'No project' : p.name) : r.screen === 'home' ? 'Home'
    : !SCREENS.includes(r.screen) || (r.screen === 'run' && r.id) || (r.screen === 'project' && r.id) ? 'Not found' : 'Home';
  // Home is the page that gets bookmarked, and its title is what the bookmark is called
  document.title = where === 'Home' ? 'FlowCellar' + (lib ? ' - ' + lib : '') : where + ' · ' + (lib ? lib + ' · ' : '') + 'FlowCellar';
});
// a new screen starts at its top (a reason that leads to its evidence scrolls there by itself)
let lastPlace = '';
effect(() => { const r = route.value, place = r.screen + '/' + (r.id || ''); if (place !== lastPlace) { lastPlace = place; if (!r.q.ev) window.scrollTo(0, 0); } });

const root = document.getElementById('app');
function say(text, cls) { const p = document.createElement('p'); p.className = cls || 'fc-boot'; p.textContent = text; return p; }
function fail() {
  const st = indexState.value;
  const text = st.state === 'old_browser' ? 'This browser is too old to draw the charts of FlowCellar (it cannot unpack the data). Firefox 113, Chrome or Edge 80 and Safari 16.4, or newer, can.'
    : st.state === 'missing' ? 'This folder holds no library index. Run flowcellar hub.'
    : isReport ? 'FlowCellar could not read the data inside this file: ' + (st.error || 'unknown reason') + '.'
    : 'FlowCellar could not read the library index: ' + (st.error || 'unknown reason') + '. Run flowcellar hub to write it again.';
  // the text summary a report carries stays on the page: it is all that can be shown
  if (!isReport || st.state === 'missing') root.textContent = '';
  const p = say(text, 'fc-boot fc-boot--failed'); p.setAttribute('role', 'alert'); root.appendChild(p);
}
if (isReport && boot.run) setHome('/run/' + boot.run);
if (isReport && boot.analysis) setHome('/analysis/' + boot.analysis);
applyTheme(themeChoice.value, false);
if (!canInflate) { indexState.value = { state: 'old_browser', error: null }; fail(); }
else {
  if (!isReport) { root.textContent = ''; root.appendChild(say('Reading the library …')); }
  loadIndex().then(ix => { useLibraryStorage(isReport ? 'report' : ix.library && ix.library.id); root.textContent = ''; render(<App />, root); }).catch(fail);
}
