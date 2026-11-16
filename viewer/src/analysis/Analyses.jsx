// The Analyses tab (the single-cell and spatial samples of a library), the short block of them on Home, and the
// panel on a run page that lists the analyses made from that run. All three are drawn from the index alone; only
// the small map on a run's card reads the analysis itself.
import { useState } from 'preact/hooks';
import { index, useBundle, isReport } from '../core/data.js';
import { param, setParams } from '../core/state.js';
import { Badge, Chip, Mid, Tile, Updated, When, LibraryNotices } from '../shell/parts.jsx';
import { Section } from '../run/section.jsx';
import * as F from '../core/format.js';
import { VERDICT_WORD } from '../core/labels.js';
import { RunChips } from './Analysis.jsx';
import { MapCanvas, useChannels } from './SpatialMap.jsx';
import * as S from './sc.js';

const pipelines = rows => rows.map(r => (S.PIPELINE_WORD[r.pipeline] || r.pipeline) + (r.version ? ' ' + r.version : '')).join(' · ');
const chemistry = rows => { const r = rows.find(x => x.chemistry) || rows.find(x => x.kit); return r ? r.chemistry || r.kit : null; };
const pick = (rows, f) => { for (const r of rows) { const v = f(r); if (v != null) return v; } return null; };
const flags = rows => { const out = []; const alarms = rows.reduce((t, r) => t + (r.alarms || 0), 0);
  if (alarms) out.push([F.plural(alarms, 'alarm'), 'Alarms the pipeline itself raised: quoted on the page of the analysis']);
  if (rows.some(r => r.from_report)) out.push(['from report', 'Some numbers exist only in a report’s HTML and were read from it']);
  if (rows.some(r => r.partial)) out.push(['partial', 'Part of the analysis was not read: its page says what']);
  if (rows.some(r => r.source_state === 'missing')) out.push(['folder gone', 'The folder of this analysis was not found by the last scan; what was read before is kept']);
  if (rows.some(r => (r.copies || 1) > 1)) out.push(['copies', 'The same output lies in more than one folder: it is one analysis']);
  return out; };
// the sentence under the badge of a row: the first reason of the analysis that decides the verdict (of the two
// pipelines of a pair, the one with the worse verdict), so that a Fail is never explained by a warning
const RANK = { fail: 3, warn: 2, incomplete: 1 };
const decides = rows => { const by = rows.slice().sort((a, b) => (RANK[b.verdict] || 0) - (RANK[a.verdict] || 0)); return by.find(x => (x.reasons || []).length) || by[0]; };
const first = rows => { const r = decides(rows); return (r.reasons && r.reasons[0]) || r.headline || ''; };
// The address a row opens: the page of the sample, on the pipeline whose verdict the row shows, so that a row that
// says Fail never opens a page that says Warning (the lead of a pair is the Space Ranger output; the other is "?p=").
export const hrefOfRow = rows => { const lead = rows[0], d = decides(rows); return S.hrefAnalysis(lead.id) + (d !== lead && (RANK[d.verdict] || 0) > (RANK[lead.verdict] || 0) ? '?p=' + encodeURIComponent(d.pipeline) : ''); };

// what one entry (a sample: one or two pipelines) shows in the table
function cells(e) {
  const rows = e.rows, lead = rows[0], u = S.unitWords(lead.assay);
  const units = pick(rows, r => r.spatial_units_under_tissue != null ? r.spatial_units_under_tissue : null), unitsAlt = pick(rows, r => S.unitsOf(r));
  const genes = pick(rows, r => S.perUnit(r, 'genes')), umis = pick(rows, r => S.perUnit(r, 'umis')), sat = pick(rows, r => S.saturationOf(r));
  return { rows, lead, u, units: units != null ? units : unitsAlt, unitsWord: units != null ? u.plural + ' under tissue' : lead.assay === 'singlecell' ? 'cells' : u.plural + ' with data',
    genes, umis, sat, under: pick(rows, r => r.spatial_fraction_reads_under_tissue), verdict: S.worstVerdict(rows), links: S.linksOf(rows), date: pick(rows, r => r.date) };
}
const SORT = { sample: c => (c.lead.label || '').toLowerCase(), kind: c => c.lead.assay || '', units: c => c.units, genes: c => c.genes && c.genes.v, umis: c => c.umis && c.umis.v, sat: c => c.sat && c.sat.v,
  under: c => c.under, verdict: c => ({ fail: 3, warn: 2, pass: 1 }[c.verdict] || 0), date: c => c.date || '' };
const DESC = new Set(['units', 'genes', 'umis', 'sat', 'under', 'verdict', 'date']);

export function AnalysesTable({ entries, limit, sortable }) {
  // the sort of the Analyses page is in its address, like the sort of Runs (?sort=verdict, ?sort=-date); the short
  // table on Home is not sorted by the reader
  const [local, setLocal] = useState(''), sort = sortable ? param('sort', '') : local, setSort = v => sortable ? setParams({ sort: v }) : setLocal(v);
  const key = sort.replace(/^-/, ''), dir = sort[0] === '-' ? -1 : 1, get = SORT[key];
  let list = entries.map(cells); const twice = S.repeatedLabels(entries);
  if (get) list = list.slice().sort((a, b) => { const x = get(a), y = get(b); if (x == null || y == null) return x == null && y == null ? 0 : x == null ? 1 : -1; return x < y ? -dir : x > y ? dir : 0; });
  if (limit) list = list.slice(0, limit);
  const Th = ({ k, cls, children, title }) => { const on = sortable && key === k, d = on ? (dir === -1 ? 'descending' : 'ascending') : undefined;
    return <th scope="col" class={(cls || '') + (on ? ' is-on' : '')} aria-sort={d} title={title}>{sortable ? <button type="button" onClick={() => setSort(key === k ? (dir === 1 ? '-' + k : k) : (DESC.has(k) ? '-' + k : k))}><span>{children}{on ? <span aria-hidden="true">{dir === -1 ? ' ▾' : ' ▴'}</span> : null}</span></button> : children}</th>; };
  return <div class="fc-tablewrap"><table class="fc-table fc-table--comfortable fc-projects fc-analyses">
    <caption class="fc-sr">Single-cell and spatial analyses of this library{get ? '' : ', the newest first'}</caption>
    <thead><tr><Th k="sample">Sample</Th><Th k="kind">Kind</Th><th scope="col" class="c-opt">Chemistry, slide</th><Th k="units" cls="fc-num">Units</Th>
      <Th k="genes" cls="fc-num" title="Genes per unit as the pipeline states them: a median for spots and cells, a mean for Visium HD bins">Genes</Th>
      <Th k="umis" cls="fc-num" title="UMIs per unit as the pipeline states them: a median for spots and cells, a mean for Visium HD bins">UMIs</Th>
      <Th k="sat" cls="fc-num">Saturation</Th><Th k="under" cls="fc-num c-opt">Reads under tissue</Th><Th k="verdict">Verdict</Th><th scope="col">Flow cell runs</th></tr></thead>
    <tbody>{list.map(c => { const l = c.lead, href = hrefOfRow(c.rows), fl = flags(c.rows);
      return <tr key={l.id} class="fc-projects__row" onClick={e => { if (!e.target.closest('a')) location.hash = href; }}>
        <th scope="row"><a class="fc-projects__name" href={href}><Mid text={l.label || l.id} tail={12} /></a>
          <span class="fc-table__sub">{pipelines(c.rows)}{c.date ? <> {'·'} <When iso={c.date} /></> : null}</span>
          {twice.has(l.label || '') && l.folder ? <span class="fc-table__sub fc-analyses__folder" title={'Analysed more than once: this one lies in ' + l.folder}>in <Mid text={l.folder.split('/').pop()} tail={14} /></span> : null}</th>
        <td class="fc-nowrap"><span class="fc-chip fc-analyses__kind" title={S.ASSAY_WORD[l.assay] || ''}>{S.ASSAY_BADGE[l.assay] || '?'}</span> {S.ASSAY_WORD[l.assay] || F.DASH}</td>
        <td class="c-opt">{chemistry(c.rows) || F.DASH}{pick(c.rows, r => r.slide) ? <span class="fc-table__sub fc-mono">{pick(c.rows, r => r.slide)}</span> : null}</td>
        <td class="fc-num">{F.int(c.units)}<span class="fc-table__sub">{c.unitsWord}</span></td>
        <td class="fc-num">{c.genes ? S.num1(c.genes.v) : F.DASH}<span class="fc-table__sub">{c.genes ? c.genes.stat : ' '}</span></td>
        <td class="fc-num">{c.umis ? S.num1(c.umis.v) : F.DASH}<span class="fc-table__sub">{c.umis ? c.umis.stat : ' '}</span></td>
        <td class="fc-num">{c.sat ? F.pct(c.sat.v) : F.DASH}<span class="fc-table__sub">{c.sat ? c.sat.from : ' '}</span></td>
        <td class="fc-num c-opt">{F.isNum(c.under) ? F.pct(c.under) : F.DASH}</td>
        <td><Badge v={c.verdict} /><span class="fc-table__sub fc-analyses__why" title={first(c.rows)}>{first(c.rows)}</span>{fl.length ? <span class="fc-analyses__flags">{fl.map(f => <span class="fc-tag" title={f[1]}>{f[0]}</span>)}</span> : null}</td>
        <td><RunChips links={c.links} small /></td>
      </tr>; })}</tbody></table></div>;
}

export function Analyses() {
  const ix = index.value, rows = ix.analyses || [], list = S.entries(rows);
  const n = { visium: 0, visium_hd: 0, singlecell: 0 }; for (const k of Object.keys(n)) n[k] = S.sampleCount(rows.filter(r => r.assay === k));
  const kinds = Object.keys(n).filter(k => n[k]).map(k => n[k] + ' ' + S.ASSAY_WORD[k]).join(', '), look = list.filter(e => ['warn', 'fail'].includes(S.worstVerdict(e.rows))).length;
  const missing = new Set(); for (const r of rows) for (const ln of r.links || []) if (!ln.resolved && ln.fc) missing.add(ln.fc);
  return <div class="fc-page fc-projectspage fc-analysespage">
    <header class="fc-home__head"><div class="fc-home__title"><h1 class="fc-home__name">Analyses</h1>
      <p class="fc-home__line">{F.plural(rows.length, 'analysis', 'analyses')} of {F.plural(S.sampleCount(rows), 'single-cell or spatial sample')} in {ix.library.name || 'this library'}{kinds ? ' (' + kinds + ')' : ''} {'·'} <Updated /></p></div></header>
    <LibraryNotices cls="fc-home__notices" />
    {list.length ? <section class="fc-panel fc-home__block" aria-label="Analyses"><AnalysesTable entries={list} sortable /></section>
      : <div class="fc-empty fc-page__empty"><b>No analysis yet.</b><span>An analysis is one output folder of wf-single-cell or Space Ranger. A scan finds them under the scanned folders; one folder is added with <span class="fc-mono">flowcellar analysis DIR -L {ix.library.path || '<this library>'}</span>.</span></div>}
    {list.length ? <p class="fc-section__note">{look ? F.plural(look, 'row') + (look === 1 ? ' needs' : ' need') + ' a look. ' : 'No analysis has a warning or a failed check. '}A verdict here is about the tissue, the library and the workflow: it never changes the verdict of a flow cell run.
      {missing.size ? ' ' + F.plural(missing.size, 'flow cell') + ' named by the reads ' + (missing.size === 1 ? 'is' : 'are') + ' not in the library (' + [...missing].slice(0, 6).join(', ') + (missing.size > 6 ? ', …' : '') + '): the link appears when the run folder is scanned.' : ''}</p> : null}
    {list.length ? <p class="fc-section__note">The ONT workflow analysis and the Space Ranger output of the same reads are one row; its page has both. A sample that was analysed again is another row, named by its folder. Genes and UMIs are per unit as the pipeline states them: a median for spots and cells, a mean for Visium HD bins.</p> : null}
  </div>;
}

// Home: the newest samples, one line each. Shown only when the library holds an analysis.
export function HomeAnalyses() {
  const ix = index.value, rows = (ix && ix.analyses) || []; if (!rows.length) return null;
  const list = S.entries(rows), nFail = list.filter(e => S.worstVerdict(e.rows) === 'fail').length, nWarn = list.filter(e => S.worstVerdict(e.rows) === 'warn').length;
  // "Needs attention" above lists flow cell runs only: what the analyses need is said here, with the way to them
  const look = [nFail ? nFail + ' failed' : null, nWarn ? nWarn + ' with a warning' : null].filter(Boolean).join(', ');
  return <section class="fc-panel fc-home__block fc-home__block--analyses" aria-labelledby="h-analyses">
    <div class="fc-panel__head"><h2 class="fc-section__title" id="h-analyses">Analyses</h2><span class="fc-secondary">single-cell and spatial samples, with the flow cell runs they came from</span>
      <a class="fc-home__all" href="#/analyses">All {F.plural(rows.length, 'analysis', 'analyses')}</a></div>
    <p class="fc-home__anlook">{look ? <><a href="#/analyses?sort=-verdict">{look}</a>, among {F.plural(list.length, 'sample')}. </> : 'No analysis has a warning or a failed check. '}<span class="fc-secondary">The four newest are listed; “Needs attention” above is about flow cell runs only, and the search box finds runs, not analyses.</span></p>
    <AnalysesTable entries={list} limit={4} />
  </section>;
}

// ---------- the run view: the analyses made from this run ----------
function Thumb({ id, label }) {
  const b = useBundle(id), grid = b ? (S.cellOf(b).grids || {}).spatial_map : null, chans = useChannels(grid);
  if (!grid || !chans.umis) return <div class="fc-anthumb fc-anthumb--none" aria-hidden="true">{b ? 'no map' : ''}</div>;
  return <div class="fc-anthumb"><MapCanvas grid={grid} chans={chans} id="umis" outline still label={'Map of ' + label + ', coloured by UMIs'} /></div>;
}
function Card({ fam, run, folder }) {
  const lead = fam.find(r => r.pipeline === 'spaceranger') || fam[0], c = cells({ rows: fam }), u = c.u, mine = c.links.find(l => l.run === run), others = c.links.filter(l => l.run !== run);
  const verdict = S.worstVerdict(fam), why = first(fam.slice().sort((a, b) => (b.reasons || []).length - (a.reasons || []).length));
  return <div class="fc-ancard">
    <Thumb id={lead.id} label={lead.label || lead.id} />
    <div class="fc-ancard__main">
      <div class="fc-ancard__head"><a class="fc-ancard__name" href={S.hrefAnalysis(lead.id)}><Mid text={lead.label || lead.id} tail={12} /></a><Chip>{S.ASSAY_WORD[lead.assay] || 'Assay not known'}</Chip>{folder && lead.folder ? <Chip title={'Analysed more than once: this one lies in ' + lead.folder}>in {F.midCut(lead.folder.split('/').pop(), 34, 14)}</Chip> : null}{fam.map(r => <Chip>{(S.PIPELINE_WORD[r.pipeline] || r.pipeline) + (r.version ? ' ' + r.version : '')}</Chip>)}
        <Badge v={verdict} /><span class="fc-ancard__why">{verdict === 'pass' ? '' : why}</span></div>
      <p class="fc-ancard__share">{mine && (mine.method === 'sample_name_chain' || mine.confidence === 'low') ? <>Linked to this flow cell by sample name only: nothing in the files says how many of its reads came from it{others.length ? '. Also linked: ' : '.'}</>
        : mine && F.isNum(mine.share) ? <><b>{S.cap(S.aboutShare(mine.share))}</b> of this analysis’s reads came from this flow cell (of the records that were looked at){others.length ? '; the rest from ' : '.'}</> : <>Made from this flow cell{others.length ? ' and ' : '.'}</>}
        {others.map((l, k) => <>{k ? ', ' : ''}<span class="fc-mono">{l.fc}</span>{l.resolved ? '' : ' (not in library)'}</>)}{others.length ? '. The analysis is of all of them together, not of this flow cell alone.' : ''} <span class="fc-muted">{mine ? S.cap(S.LINK_WORD[mine.method] || '') + (mine.confidence && mine.confidence !== 'high' ? ' (' + mine.confidence + ' confidence)' : '') + '.' : ''}</span></p>
      <div class="fc-tiles fc-ancard__tiles">
        <Tile label={S.cap(c.unitsWord)} value={F.int(c.units)} />
        <Tile label={'UMIs per ' + u.name.replace('8 µm ', '') + (c.umis ? ', ' + c.umis.stat : '')} value={c.umis ? S.num1(c.umis.v) : F.DASH} />
        {F.isNum(c.under) ? <Tile label="Reads under tissue" value={F.pctNum(c.under)} unit="%" /> : <Tile label={'Genes per ' + u.name.replace('8 µm ', '') + (c.genes ? ', ' + c.genes.stat : '')} value={c.genes ? S.num1(c.genes.v) : F.DASH} />}
        <Tile label="Saturation" value={c.sat ? F.pctNum(c.sat.v) : F.DASH} unit={c.sat ? '%' : ''} caption={c.sat ? c.sat.from : null} />
      </div>
    </div>
    <a class="fc-btn fc-ancard__open" href={S.hrefAnalysis(lead.id)}>Open analysis</a>
  </div>;
}
// rows of the analyses whose reads came from one run, grouped by sample
export function analysesOfRun(id) { const ix = index.value, rows = (ix && ix.analysesByRun && ix.analysesByRun.get(id)) || []; return rows.length ? S.entries(ix.analyses).filter(e => e.rows.some(r => rows.includes(r))) : []; }
export function RunAnalyses({ b }) {
  if (isReport) return null;
  const id = b.identity.run_key, list = analysesOfRun(id); if (!list.length) return null;
  const worst = S.worstVerdict(list.flatMap(e => e.rows)), twice = S.repeatedLabels(list.map(e => ({ rows: [e.rows.find(r => r.pipeline === 'spaceranger') || e.rows[0]] })));
  return <Section id="analyses" title="Single-cell and spatial" reading="Analysis verdicts are about the tissue and the workflow. They never change this flow cell’s verdict."
    tools={<span class="fc-secondary">{F.plural(list.length, 'analysis', 'analyses')} from this flow cell {'·'} {VERDICT_WORD[worst] || worst}</span>}>
    {list.map(e => <Card fam={e.rows} run={id} folder={twice.has((e.rows.find(r => r.pipeline === 'spaceranger') || e.rows[0]).label || '')} />)}
  </Section>;
}
