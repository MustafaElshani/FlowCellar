// Fixed section anatomy: title, result of the section, one sentence written from the data, tiles, charts, then
// closed drawers "More charts" and "All numbers".
import { useState } from 'preact/hooks';
import { Badge, Glyph, Tile } from '../shell/parts.jsx';
import { setParams, param, makeHash, route } from '../core/state.js';
import { checkName, metric, labelOf, VERDICT_WORD } from '../core/labels.js';
import * as F from '../core/format.js';
import * as D from '../core/derive.js';
import { fmtValue } from '../core/values.js';

// what one stored check shows: a triggered check warns or fails by its severity (a key check may fail, a supporting one
// warns at most); a triggered note (severity info or none) changes nothing and reads as passed
export const ruleState = r => r.outcome === 'triggered' ? (r.severity === 'fail' ? (r.role === 'key' || r.role == null ? 'fail' : 'warn') : r.severity === 'warn' ? 'warn' : 'pass') : r.outcome === 'not_evaluated' ? 'unknown' : 'pass';

export function Section({ id, title, result, reading, x, tools, children, soFar }) {
  const ev = x && x.evSec === id ? x.evRule : null;
  return <section class={'fc-section' + (ev ? ' is-evidence' : '')} id={'s-' + id} data-sec={id} aria-labelledby={'h-' + id}>
    <div class="fc-section__head"><h2 class="fc-section__title" id={'h-' + id}>{title}</h2>
      {result ? <Badge v={result} soFar={soFar} /> : null}{tools ? <span class="fc-section__tools">{tools}</span> : null}</div>
    {reading ? <p class="fc-section__reading">{reading}</p> : null}
    {ev ? <div class="fc-evidence" role="status"><Glyph v={ruleState(ev)} /><span>Showing the evidence for <b>{checkName(ev)}</b>: {ev.message}</span>
      <button type="button" class="fc-btn" onClick={() => { setParams({ ev: '' }); window.scrollTo(0, 0); }}>Back to summary</button></div> : null}
    {children}
  </section>;
}

// open and closed drawers are in the address (more=tracks,fc,...). charts: a group of further charts (not a table of
// numbers): "PDF with every chart" opens these for the length of a print, whatever the address says (run/common.js
// puts them into c.more then; a group marked `charts` must be listed there in CHART_GROUPS).
export function Disclosure({ id, title, c, children, force, charts }) {
  const open = force || c.more.has(id);
  const toggle = () => { const m = new Set(c.more); m.has(id) ? m.delete(id) : m.add(id); setParams({ more: [...m].join(',') }); };
  return <><button type="button" class={'fc-disclosure' + (open ? ' is-open' : '')} data-group={charts ? 'charts' : 'numbers'} aria-expanded={open} onClick={toggle}>{title}</button>{open ? <div class="fc-disclosure__body">{children}</div> : null}</>;
}

export function allSummary(b) { const o = {}; for (const k of Object.keys(b.modules || {})) Object.assign(o, b.modules[k].summary || {}); return o; }
export { fmtValue };
const split = s => { const m = /^(.*?)[\s ]([^\s \d][^\d]*)$/.exec(s); return m ? [m[1], m[2]] : [s, '']; };

// supporting numbers of a section, in the catalogue's tile order
export function SectionTiles({ b, ids, extra = {} }) {
  const S = allSummary(b), pending = new Set(b.status.pending || []), rules = b.status.rules || [];
  const href = r => makeHash(route.value.path, { ...route.value.q, ev: r.id });
  return <div class="fc-tiles fc-section__tiles">{ids.map(id => {
    if (extra[id]) return extra[id];
    const r = rules.find(q => q.id === id), m = metric(id), key = r && r.metric && S[r.metric] != null ? r.metric : id, v = S[key], name = labelOf(id);
    if (pending.has(id)) return <Tile label={name.label} state="unknown" value="being counted" caption="waits for the read count" title={id} />;
    if (v == null) return null;
    const st = r ? ruleState(r) : null, [val, unit] = split(fmtValue(key, v));
    const i = r && r.message ? r.message.indexOf('. ') : -1, cap = r && i > 0 && st !== 'unknown' ? r.message.slice(i + 2).replace(/\.$/, '') : null;
    return <Tile label={name.label} glyph={st && st !== 'unknown' ? st : null} state={st === 'warn' || st === 'fail' ? st : null} value={val} unit={unit} caption={cap} href={r ? href(r) : undefined} title={m.how || id} />;
  })}</div>;
}

// Every stored number of the given sections. A number without display copy is named from its id in words; the id
// is always in its own column and is the row's title, so nothing is ever shown as a bare identifier.
export function AllNumbers({ b, mods }) {
  const rows = []; for (const id of mods) { const m = b.modules[id]; if (!m) continue; for (const [k, v] of Object.entries(m.summary || {})) if (!k.startsWith('_')) rows.push([k, v, id]); }
  return <div class="fc-tablewrap fc-allnumbers"><table class="fc-table"><caption class="fc-sr">All numbers of this section</caption><thead><tr><th scope="col">Number</th><th scope="col" class="fc-num">Value</th><th scope="col">Id</th></tr></thead>
    <tbody>{rows.map(([k, v]) => <tr title={k}><th scope="row">{labelOf(k).label}</th><td class="fc-num">{fmtValue(k, v)}</td><td class="fc-mono fc-muted">{k}</td></tr>)}</tbody></table></div>;
}
export const countNumbers = (b, mods) => mods.reduce((a, id) => a + Object.keys((b.modules[id] || {}).summary || {}).filter(k => !k.startsWith('_')).length, 0);
export { VERDICT_WORD };
