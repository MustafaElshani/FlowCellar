// Small shared components. Class names are those of styles/components.css.
import brand from '../gen/brand.js';
import { VERDICT_WORD } from '../core/labels.js';
import { route, go, themeChoice, applyTheme, selected, refocus, param } from '../core/state.js';
import { index, isReport, reloadIndex, reloading, boot } from '../core/data.js';
import { sparkPoints } from '../core/derive.js';
import * as F from '../core/format.js';
import { useRef, useEffect } from 'preact/hooks';

// A verdict is a shape and a word; fill density and colour repeat them. "not_available" (a part without data)
// and "unknown" (data, but no check could be assessed) share the grey dash and differ in their word.
const BADGE_CLASS = { pass: 'pass', warn: 'warn', fail: 'fail', incomplete: 'incomplete', running: 'running' };
export function Badge({ v, lg, soFar, word }) {
  const k = v || 'unknown';
  return <span class={'fc-badge fc-badge--' + (BADGE_CLASS[k] || 'unknown') + (lg ? ' fc-badge--lg' : '')}>{word || (VERDICT_WORD[k] || k) + (soFar && (k === 'pass' || k === 'warn') ? ' so far' : '')}</span>;
}
export const Glyph = ({ v, label }) => { const k = v === 'ok' ? 'pass' : v; return <span class={'fc-glyph fc-glyph--' + (BADGE_CLASS[k] && k !== 'incomplete' && k !== 'running' ? k : 'unknown')} role="img" aria-label={label || VERDICT_WORD[k] || k} />; };
// A date or a time with its zone in the tooltip (and beside it when `zone` is asked for). `floating`: a time read
// from a folder name, shown as written (see core/format.js).
export function When({ iso, floating, time, zone }) {
  if (!iso) return <span class="fc-muted">{F.DASH}</span>;
  return <time class="fc-when" dateTime={iso} title={F.timeTitle(iso, floating)}>{time ? F.dateTime(iso, floating) : F.date(iso, floating)}{zone && !floating ? ' ' + F.zone(iso) : ''}</time>;
}
export const Chip = ({ children, href, title, note }) => href
  ? <a class={'fc-chip' + (note ? ' fc-chip--note' : '')} href={href} title={title}>{children}</a>
  : <span class={'fc-chip' + (note ? ' fc-chip--note' : '')} title={title}>{children}</span>;
// verdict tally as small counts: "✓ 4  ▲ 1  ✕ 2", each with its word for a screen reader
const TALLY = ['fail', 'warn', 'pass', 'incomplete', 'running', 'unknown'];
export function Tally({ counts }) {
  const c = counts || {}, items = TALLY.filter(k => c[k] > 0);
  if (!items.length) return <span class="fc-muted">{F.DASH}</span>;
  return <span class="fc-tally">{items.map(k => <span class="fc-tally__item" title={c[k] + ' ' + VERDICT_WORD[k]}>{k === 'incomplete' || k === 'running'
    ? <span class={'fc-tally__g fc-tally__g--' + k} role="img" aria-label={VERDICT_WORD[k]} /> : <Glyph v={k} />}<span>{c[k]}</span></span>)}</span>;
}

// Middle truncation: the head shrinks with an ellipsis, the last 10 characters are kept.
export function Mid({ text, tail = 10, cls }) {
  if (!text) return null;
  if (text.length <= tail + 6) return <span class={cls} title={text}>{text}</span>;
  return <span class={'fc-trunc-mid' + (cls ? ' ' + cls : '')} title={text}><span>{text.slice(0, -tail)}</span><span>{text.slice(-tail)}</span></span>;
}

// Left and Right (Up and Down) move through the choices of a radio group and choose, as in any radio group.
export function radioKeys(e) {
  const d = e.key === 'ArrowRight' || e.key === 'ArrowDown' ? 1 : e.key === 'ArrowLeft' || e.key === 'ArrowUp' ? -1 : 0; if (!d) return;
  const all = [...e.currentTarget.querySelectorAll('button[role="radio"]')], i = all.indexOf(document.activeElement); if (i < 0 || all.length < 2) return;
  e.preventDefault(); const next = all[(i + d + all.length) % all.length]; next.focus(); next.click();
}
// A radio group drawn as a segmented control. opts: [{ v, label, off: reason }]
export function Seg({ label, value, opts, onChange, name }) {
  const pick = v => { onChange(v); refocus('[data-seg="' + name + ':' + v + '"]'); };
  return <span class="fc-ctrl"><span class="fc-caps" id={name + '-l'}>{label}</span>
    <span class="fc-seg" role="radiogroup" aria-labelledby={name + '-l'} onKeyDown={radioKeys}>{opts.map(o => o.off
      ? <span key={o.v} role="radio" aria-checked="false" aria-disabled="true" class="is-off" tabIndex={0} title={o.off}>{o.label}</span>
      : <button key={o.v} type="button" role="radio" data-seg={name + ':' + o.v} aria-checked={String(value === o.v)} class={value === o.v ? 'is-on' : ''} onClick={() => pick(o.v)}>{o.label}</button>)}</span></span>;
}
export function Select({ label, value, opts, onChange, name }) {
  return <label class="fc-ctrl"><span class="fc-caps">{label}</span>
    <select class="fc-selectbox" name={name} value={value} onChange={e => onChange(e.currentTarget.value)}>{opts.map(o => <option value={o.v} disabled={!!o.off}>{o.label}</option>)}</select></label>;
}

// Key-number tile: label with the check glyph, value and unit, the strip, a caption.
// strip: { lo, hi, band: [a, b], mark, limits: [v], better: 'more' | 'less' }
export function Tile({ label, glyph, value, unit, caption, strip, state, href, onClick, title, sub }) {
  const cls = 'fc-tile' + (state === 'warn' ? ' is-triggered' : state === 'fail' ? ' is-triggered is-fail' : state === 'unknown' ? ' is-unknown' : '');
  const pos = v => Math.max(0, Math.min(100, (v - strip.lo) / (strip.hi - strip.lo || 1) * 100)) + '%';
  const body = <>
    <span class="fc-tile__label"><span>{label}</span>{glyph ? <Glyph v={glyph} /> : null}</span>
    <span class="fc-tile__value">{value}{unit ? <small> {unit}</small> : null}</span>
    {strip ? <span class="fc-strip" aria-hidden="true">
      {strip.band ? <span class="fc-strip__band" style={{ left: pos(strip.band[0]), width: `calc(${pos(strip.band[1])} - ${pos(strip.band[0])})` }} /> : null}
      {(strip.limits || []).map(l => <span class="fc-strip__limit" style={{ left: pos(l) }} />)}
      <span class="fc-strip__mark" style={{ left: pos(strip.mark) }} />
    </span> : null}
    {caption ? <span class="fc-tile__caption">{caption}</span> : null}
    {sub ? <span class="fc-tile__caption">{sub}</span> : null}
  </>;
  if (href) return <a class={cls} href={href} title={title} onClick={onClick}>{body}</a>;
  if (onClick) return <button type="button" class={cls} title={title} onClick={onClick}>{body}</button>;
  return <div class={cls} title={title}>{body}</div>;
}

// Sparkline: one SVG path, fixed 0 to 96 h axis, a filled dot where a run stopped early.
export function Spark({ s, w = 84, h = 18, early, typ, hours }) {
  const pts = sparkPoints(s, hours), n = 32, sx = w / (n - 1);
  if (!pts.length) return <svg class="fc-spark" width={w} height={h} aria-hidden="true"><line x1="0" y1={h - 0.5} x2={w} y2={h - 0.5} /></svg>;
  let d = ''; pts.forEach((v, i) => { d += (i ? 'L' : 'M') + (i * sx).toFixed(1) + ' ' + (h - 1.5 - v * (h - 4)).toFixed(1); });
  const lx = (pts.length - 1) * sx, ly = h - 1.5 - pts[pts.length - 1] * (h - 4);
  return <svg class="fc-spark" width={w} height={h} viewBox={`0 0 ${w} ${h}`} aria-hidden="true">
    {typ ? <path class="typ" d={typ} /> : null}
    <line x1="0" y1={h - 0.5} x2={w} y2={h - 0.5} /><path d={d} />
    {early && pts.length < n ? <circle cx={lx.toFixed(1)} cy={ly.toFixed(1)} r="2.2" /> : null}
  </svg>;
}

export function Logo() {
  return <svg viewBox={brand.viewBox} role="img" aria-label="FlowCellar" focusable="false">
    {brand.parts.map(p => <path fill={p.fill} d={p.d} />)}
    <g transform={brand.wordTransform}>{brand.word.map(p => <path fill={p.fill} d={p.d} />)}</g>
  </svg>;
}

// Three choices side by side; in a narrow window one button that steps through them (it never leaves the bar).
const THEMES = ['auto', 'light', 'dark'], themeWord = t => t[0].toUpperCase() + t.slice(1);
const ThemeSwitch = () => <>
  <span class="fc-seg fc-seg--inverse fc-theme" role="radiogroup" aria-label="Theme" onKeyDown={radioKeys}>{THEMES.map(t =>
    <button key={t} type="button" role="radio" data-seg={'theme:' + t} aria-checked={String(themeChoice.value === t)} class={themeChoice.value === t ? 'is-on' : ''} onClick={() => { applyTheme(t); refocus('[data-seg="theme:' + t + '"]'); }}>{themeWord(t)}</button>)}</span>
  <button type="button" class="fc-btn fc-btn--inverse fc-theme1" title="Theme: Auto, Light or Dark. Press to change." aria-label={'Theme: ' + themeWord(themeChoice.value) + '. Press to change.'}
    onClick={() => applyTheme(THEMES[(THEMES.indexOf(themeChoice.value) + 1) % 3])}>{themeWord(themeChoice.value)}</button>
</>;

// One application bar on every screen: logo (goes Home), tabs, search, compare, theme. A one-run report has no
// library behind it, so its bar carries only the logo, what the file is, and the theme.
// The bar may grow (larger text, a narrow window): what sticks under it follows its real height.
function useBarHeight() {
  const ref = useRef(null);
  useEffect(() => {
    const el = ref.current; if (!el || typeof ResizeObserver !== 'function') return;
    const set = () => { const h = Math.round(el.getBoundingClientRect().height); if (h > 0) document.documentElement.style.setProperty('--fc-appbar-h', h + 'px'); };
    const ro = new ResizeObserver(set); ro.observe(el); set();
    return () => ro.disconnect();
  }, []);
  return ref;
}
export function AppBar({ q, onQuery }) {
  const ix = index.value, r = route.value, sel = selected.value, bar = useBarHeight();
  if (isReport) return <header ref={bar} class="fc-appbar fc-surface-inverse"><div class="fc-appbar__in">
    <a class="fc-appbar__logo" href="#/" aria-label="FlowCellar: the top of this report"><Logo /></a>
    <span class="fc-appbar__what">{boot.analysis ? 'Analysis report' : 'Run report'} {'·'} a snapshot{ix && ix.written ? ' made ' + F.date(ix.written) : ''}. Self-contained: works offline, sends nothing.</span>
    <ThemeSwitch />
  </div></header>;
  // the tab that is lit says where one is; an address that is no page of the library lights none
  const tab = r.screen === 'home' ? 'home' : r.screen === 'projects' || r.screen === 'project' ? 'projects' : r.screen === 'analyses' || r.screen === 'analysis' ? 'analyses' : ['runs', 'run', 'compare', 'flowcell'].includes(r.screen) ? 'runs' : null;
  const submit = e => { e.preventDefault(); const v = new FormData(e.currentTarget).get('q') || ''; if (r.screen !== 'runs') { go('/runs', { q: v }); refocus('.fc-appbar__search input'); } };
  const Tab = ({ id, href, children }) => <a class={'fc-tab' + (tab === id ? ' is-on' : '')} href={href} aria-current={tab === id && r.screen === id ? 'page' : undefined}>{children}</a>;
  // on the compare page the count and the link follow the page itself (its address), elsewhere the ticked runs
  const onCompare = r.screen === 'compare', items = onCompare ? compareItems(param('items', '')).ids : sel.slice(0, 8), nSel = onCompare ? items.length : sel.length;
  // the bar is the same on every screen: its content starts and ends on the gutter lines of the page (styles: --fc-gutter)
  return <header ref={bar} class="fc-appbar fc-surface-inverse"><div class="fc-appbar__in">
    <a class="fc-appbar__logo" href="#/home" aria-label="FlowCellar: Home"><Logo /></a>
    <nav class="fc-tabs" aria-label="Library">
      <Tab id="home" href="#/home">Home</Tab>
      <Tab id="runs" href="#/runs">Runs <span class="fc-tab__n">{ix ? ix.n : ''}</span></Tab>
      <Tab id="projects" href="#/projects">Projects <span class="fc-tab__n">{ix ? ix.projects.filter(p => p.id !== '_none').length : ''}</span></Tab>
      {ix && ix.analyses.length ? <Tab id="analyses" href="#/analyses">Analyses <span class="fc-tab__n">{ix.analyses.length}</span></Tab> : null}
    </nav>
    {r.screen === 'home' ? <span class="fc-grow" /> : <form class="fc-search fc-search--inverse fc-appbar__search" role="search" onSubmit={submit}>
      <input class="fc-search__input" type="search" name="q" aria-label="Search runs" autocomplete="off" spellcheck={false}
        placeholder="Search flow cell, run, sample, project … or type kit: status: after:" value={q != null ? q : ''} onInput={onQuery ? e => onQuery(e.currentTarget.value) : undefined} />
    </form>}
    <a class={'fc-btn fc-btn--inverse fc-appbar__compare' + (nSel < 2 ? ' is-quiet' : '')} href={nSel >= 2 ? '#/compare?items=' + items.join(',') : '#/runs'} aria-current={onCompare ? 'page' : undefined}
      title={nSel < 2 ? 'Tick two to eight runs in the Runs table to compare them' : onCompare ? 'The runs compared on this page' : 'Compare the ticked runs'}>Compare ({nSel})</a>
    <ThemeSwitch />
  </div></header>;
}
// The runs a compare address names: known runs only, each once, eight at most; what was left out is named.
export function compareItems(text) {
  const ix = index.value, asked = String(text || '').split(',').map(s => s.trim()).filter(Boolean), ids = [], unknown = [], over = [];
  for (const id of asked) { if (ids.includes(id)) continue; if (!ix || ix.rowOf.get(id) == null) unknown.push(id); else if (ids.length >= 8) over.push(id); else ids.push(id); }
  return { ids, unknown, over };
}
export const Crumbs = ({ children }) => <nav class="fc-crumbs" aria-label="Breadcrumb">{children}</nav>;

// "updated 13 Nov 2026, 19:47" with a control that reads the library again (a scan may have added to it)
export function Updated({ word = 'updated' }) {
  const ix = index.value; if (!ix) return null;
  return <span class="fc-updated">{word} <When iso={ix.written} time zone />{isReport ? null : <> <button type="button" class="fc-linkbtn" disabled={reloading.value} onClick={() => reloadIndex()} title="Read the library again: a scan may have added runs since this page was opened">{reloading.value ? 'Reading …' : 'Reload'}</button></>}</span>;
}
// what the last scan said about the library, and what a running scan is doing: lines, never verdicts
// Notices that concern single runs (two candidate sources, a file that failed) are said once, on Home under
// "Where this library comes from", and on the pages of the runs and projects they concern: not on every screen.
const QUIET_NOTICES = new Set(['two_candidate_sources', 'inputs_failed', 'judged_before_older_runs', 'runs_not_walked', 'config_not_used', 'baseline_spread_guard']);
export const loudNotices = ix => ix.notices.filter(n => n.message && !QUIET_NOTICES.has(n.code));
export const quietNotices = ix => ix.notices.filter(n => n.message && QUIET_NOTICES.has(n.code));
export function LibraryNotices({ cls, noProgress }) {
  const ix = index.value; if (!ix) return null;
  const lines = loudNotices(ix).map(n => n.message);
  if (!noProgress && ix.progress && ix.progress.state === 'scanning' && ix.progress.message) lines.unshift(ix.progress.message);
  if (!lines.length) return null;
  return <div class={'fc-libnotices' + (cls ? ' ' + cls : '')} role="status">{lines.map(t => <p class="fc-notice">{t}</p>)}</div>;
}
// runs whose read data holds files that are not counted because two candidates describe the same reads
export const runsInConflict = ix => ix.rows.filter(r => (r.t2_conflicts || 0) > 0);

// "Differs from typical": words only for what is off typical, lows first; no verdict shape.
export function differs(r) {
  const parts = [['flow cell', r.use > 1 ? null : r.mark_flow_cell], ['loading', r.mark_loading], ['library', r.mark_library]];
  if (parts.every(p => p[1] == null)) return { none: true, text: 'no lab history yet' };
  const lo = parts.filter(p => p[1] === -1).map(p => p[0]), hi = parts.filter(p => p[1] === 1).map(p => p[0]);
  const join = a => a.length === 1 ? a[0] : a.length === 2 ? a[0] + ' and ' + a[1] : a.slice(0, -1).join(', ') + ' and ' + a[a.length - 1];
  const out = [];
  if (lo.length) out.push({ dir: 'lo', text: join(lo) + ' low' });
  if (hi.length) out.push({ dir: 'hi', text: join(hi) + ' high' });
  return { parts: out, text: out.map(o => o.text).join(' · ') };
}
export function Differs({ r }) {
  const d = differs(r);
  if (d.none) return <span class="fc-why"><span class="none">{d.text}</span></span>;
  if (!d.parts.length) return null;
  return <span class="fc-why" title={d.text}>{d.parts.map((p, k) => <>{k ? ' · ' : ''}<span class={p.dir}>{p.text}</span></>)}</span>;
}
// "Output at 48 h" of a run in the judged unit, with the other unit or the horizon the run reached under it
export function outputCell(r) {
  if (r.hours == null) return { a: F.DASH, b: '', sortValue: null };
  const rna = r.assay_group !== 'gdna';
  let h = 48; if (r.reads_48h == null) h = r.reads_24h != null ? 24 : r.reads_12h != null ? 12 : r.reads_6h != null ? 6 : 0;
  if (!h) return { a: F.DASH, b: 'ran ' + F.hours(r.hours), sortValue: null };
  const reads = r['reads_' + h + 'h'], gbv = r['gb_' + h + 'h'], ran = F.hours(r.hours >= 10 ? Math.round(r.hours) : r.hours);
  const a = rna ? F.count(reads) + ' reads' : F.bases(gbv * 1e9) + ' est.', b = h === 48 ? (rna ? F.bases(gbv * 1e9) + ' est.' : F.count(reads) + ' reads') : 'at ' + h + ' h · ran ' + ran;
  // sorted by the figure that is shown: a run that ended before 48 h is ranked by what it reached, and says "at 24 h"
  return { a, b, sortValue: rna ? reads : gbv, h };
}
export const runSoFar = r => !!r.provisional;
// Read data is a verdict "so far" while files of the run are linked but not read: one check assessed from a
// small table must not read as a plain Pass.
export const dataSoFar = r => (r.t2_files_waiting || 0) > 0 && (r.data === 'pass' || r.data === 'warn');
// the word on a Read data badge
export const dataWord = r => r.data === 'not_available' ? 'Not available yet' : r.data === 'unknown' ? 'Not assessed' : undefined;
// The name of a count of bases after the reads it counts: "Passed" only where the reads passed a quality filter.
export const BASES_LABEL = { basecalled_pass: 'Passed bases', basecalled_all: 'Basecalled bases', classified: 'Bases given to samples', full_length_trimmed: 'Bases in full-length reads', aligned_primary: 'Bases in aligned reads' };
export const READS_LABEL = { basecalled_pass: 'Reads passed', basecalled_all: 'Reads basecalled', classified: 'Reads given to samples', full_length_trimmed: 'Full-length reads', aligned_primary: 'Aligned reads' };
export const basesLabel = stage => BASES_LABEL[stage] || 'Bases read';
