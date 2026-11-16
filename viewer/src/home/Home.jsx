// Home: the landing page of a library, drawn from the index alone (no run file is read for it).
// It answers, for someone who has never read the documentation: what is in this library, where it came from,
// what needs a look, and how to add more.
import { index } from '../core/data.js';
import { go, hrefRun, hrefProject, refocus } from '../core/state.js';
import { useState } from 'preact/hooks';
import { Badge, Glyph, Mid, Spark, Tally, Updated, When, LibraryNotices, quietNotices, runsInConflict, outputCell, runSoFar } from '../shell/parts.jsx';
import { HomeAnalyses } from '../analysis/Analyses.jsx';
import { ChartPanel } from '../charts/ChartPanel.jsx';
import { tk } from '../charts/theme.js';
import * as F from '../core/format.js';
import { ASSAY_SHORT, ROOT_STATE_WORD } from '../core/labels.js';

const ATTENTION = 6, LATEST = 6, PROJECTS_ON_HOME = 8;
const needsLook = r => r.verdict === 'fail' || r.verdict === 'warn';
const quote = p => /^[A-Za-z0-9_@%+=:,./~-]+$/.test(p) ? p : "'" + String(p).replace(/'/g, "'\\''") + "'";

function SearchBox() {
  const submit = e => { e.preventDefault(); go('/runs', { q: new FormData(e.currentTarget).get('q') || '' }); refocus('.fc-appbar__search input'); };
  return <form class="fc-search fc-home__search" role="search" onSubmit={submit}>
    <input class="fc-search__input" type="search" name="q" aria-label="Search the runs of this library" autocomplete="off" spellcheck={false}
      placeholder="Search flow cell, run, sample, project … or type kit: status: after:" />
    <button type="submit" class="fc-btn fc-btn--primary">Search runs</button>
  </form>;
}

// the commands that add to this library, with its own path filled in
export function AddCommands({ path }) {
  const lib = path ? quote(path) : '<this library>';
  return <div class="fc-howto">
    <p><b>To add to this library</b>, scan another folder, or add one run:</p>
    <pre class="fc-cmd fc-mono" tabIndex={0} aria-label="Commands that add to this library">{'flowcellar scan /another/folder -L ' + lib + '\nflowcellar run /one/run_folder -L ' + lib}</pre>
    <p class="fc-secondary">Nothing is copied or moved: FlowCellar only reads the folders where the runs are, and keeps its own small files in the library. Scanning again reads only what is new.</p>
  </div>;
}

function EmptyLibrary({ lib, roots }) {
  const path = lib.path ? quote(lib.path) : '<this library>';
  return <div class="fc-page fc-home">
    <header class="fc-home__head"><div class="fc-home__title"><h1 class="fc-home__name">{lib.name || 'FlowCellar library'}</h1><p class="fc-home__line">No runs yet {'·'} <Updated /></p></div></header>
    <LibraryNotices cls="fc-home__notices" />
    <section class="fc-empty fc-home__empty" aria-label="An empty library">
      <b>This library has no runs yet.</b>
      <span>A library is a folder that FlowCellar fills with one small file per sequencing run: its numbers, its verdict and these pages.
        It is built by scanning the folders where your runs are; nothing is copied or moved, the folders are only read.
        Scan again whenever there are new runs: only what is new is read.</span>
      {roots.length ? <div class="fc-home__scanned"><span>{roots.length === 1 ? 'The folder that was scanned held no run folder' : 'The folders that were scanned held no run folder'} (a run folder is one MinKNOW wrote: it holds a final summary or a report, or raw data):</span>
        <ul>{roots.map(r => <li><span class="fc-mono fc-break">{r.path}</span> {'·'} {r.state && r.state !== 'reachable' ? <span class="fc-notice">{(ROOT_STATE_WORD[r.state] || r.state).toLowerCase()}</span> : '0 runs found'}{r.last_scanned ? <> {'·'} scanned <When iso={r.last_scanned} time /></> : ' · not scanned yet'}</li>)}</ul>
        <span>If that is not where the runs are, scan the right folder:</span></div> : null}
      <pre class="fc-cmd fc-mono" tabIndex={0} aria-label="The command that fills this library">{'flowcellar scan /folder/with/runs /another/folder -L ' + path}</pre>
      <span class="fc-secondary">Then reload this page. A single run can be added with <span class="fc-mono">flowcellar run /one/run_folder -L {path}</span>.</span>
    </section>
    <HomeAnalyses />
  </div>;
}

function AttentionRow({ r }) {
  return <li class="fc-homerow">
    <span class="fc-homerow__badge"><Badge v={r.verdict} soFar={runSoFar(r)} /></span>
    <span class="fc-homerow__main">
      <span class="fc-homerow__line"><a class="fc-mono fc-id" href={hrefRun(r.id)}>{r.fc}</a>{r.run_name ? <Mid text={r.run_name} tail={10} cls="fc-homerow__name" /> : null}<span class="fc-homerow__date"><When iso={r.started} floating={r.stub} /></span></span>
      <span class="fc-homerow__text" title={r.headline}>{r.headline}</span>
    </span>
  </li>;
}
function LatestRow({ r }) {
  const out = outputCell(r);
  return <li class="fc-homerow fc-homerow--latest">
    <span class="fc-homerow__badge"><Badge v={r.verdict} soFar={runSoFar(r)} /></span>
    <span class="fc-homerow__main">
      <span class="fc-homerow__line"><a class="fc-mono fc-id" href={hrefRun(r.id)}>{r.fc}</a>{r.run_name ? <Mid text={r.run_name} tail={10} cls="fc-homerow__name" /> : <span class="fc-homerow__name fc-muted">{r.stub ? 'raw data only' : 'no run name'}</span>}</span>
      <span class="fc-homerow__text"><When iso={r.started} floating={r.stub} />{r.hours != null ? ' · ' + F.hours(r.hours >= 10 ? Math.round(r.hours) : r.hours) : ''}{ASSAY_SHORT[r.assay_group] ? ' · ' + ASSAY_SHORT[r.assay_group] : ''}</span>
    </span>
    <span class="fc-homerow__spark">{r.hours != null ? <Spark s={r.spark_output} hours={r.hours} early={r.hours < (r.planned_h || 96) - 2} /> : null}</span>
    <span class="fc-homerow__out fc-num">{out.a}<span class="fc-cell__sub">{out.b || ' '}</span></span>
  </li>;
}

// The projects table: also the Projects tab, and the header numbers of a project page read the same rows.
const P_SORT = { name: p => p.name == null ? null : p.name.toLowerCase(), runs: p => p.runs, fcs: p => p.flow_cells, samples: p => p.samples || null, period: p => p.last || null,
  verdicts: p => { const v = p.verdicts || {}; return (v.fail || 0) * 1e6 + (v.warn || 0) * 1e3 + (v.incomplete || 0); }, output: p => p.estimated_bases || null };
const P_DESC = new Set(['runs', 'fcs', 'samples', 'period', 'verdicts', 'output']);
export function ProjectsTable({ projects, roots, limit, sortable }) {
  const [sort, setSort] = useState('');                 // '' = as the index lists them: the newest project first
  const key = sort.replace(/^-/, ''), dir = sort[0] === '-' ? -1 : 1, get = P_SORT[key];
  const sorted = get ? projects.slice().sort((a, b) => { const x = get(a), y = get(b); if (x == null || y == null) return x == null && y == null ? 0 : x == null ? 1 : -1; return x < y ? -dir : x > y ? dir : 0; }) : projects;
  const shown = limit ? sorted.slice(0, limit) : sorted, rootTitle = id => { const r = roots.get(id); return r ? r.path + (r.state && r.state !== 'reachable' ? ' (' + (ROOT_STATE_WORD[r.state] || r.state).toLowerCase() + ')' : '') : undefined; };
  const Th = ({ k, cls, children }) => { const on = sortable && key === k, d = on ? (dir === -1 ? 'descending' : 'ascending') : undefined;
    return <th scope="col" class={(cls || '') + (on ? ' is-on' : '')} aria-sort={d}>{sortable ? <button type="button" onClick={() => setSort(key === k ? (dir === 1 ? '-' + k : k) : (P_DESC.has(k) ? '-' + k : k))}><span>{children}{on ? <span aria-hidden="true">{dir === -1 ? ' ▾' : ' ▴'}</span> : null}</span></button> : children}</th>; };
  return <div class="fc-tablewrap"><table class="fc-table fc-table--comfortable fc-projects">
    <caption class="fc-sr">Projects of this library{get ? '' : ', the newest first'}</caption>
    <thead><tr><Th k="name">Project</Th><Th k="runs" cls="fc-num">Runs</Th><Th k="fcs" cls="fc-num c-opt">Flow cells</Th><Th k="samples" cls="fc-num c-opt">Samples</Th><Th k="period">Period</Th>
      <Th k="verdicts">Verdicts</Th><Th k="output" cls="fc-num">Output, estimated</Th><th scope="col">Where it is</th></tr></thead>
    <tbody>{shown.map(p => <tr key={p.id} class="fc-projects__row" onClick={e => { if (!e.target.closest('a')) go('/project/' + encodeURIComponent(p.id)); }}>
      <th scope="row"><a class="fc-projects__name" href={hrefProject(p.id)}>{p.name == null ? 'No project' : <Mid text={p.name} tail={10} />}</a></th>
      <td class="fc-num">{F.int(p.runs)}</td><td class="fc-num c-opt">{F.int(p.flow_cells)}</td><td class="fc-num c-opt">{p.samples ? F.int(p.samples) : F.DASH}</td>
      <td class="fc-nowrap">{F.period(p.first, p.last)}</td><td><Tally counts={p.verdicts} /></td><td class="fc-num fc-nowrap">{p.estimated_bases ? F.basesBig(p.estimated_bases) : F.DASH}</td>
      <td class="fc-projects__where">{(p.roots || []).length ? p.roots.map(id => <span class="fc-chip" title={rootTitle(id)}>{id}</span>) : F.DASH}</td>
    </tr>)}</tbody></table></div>;
}

export function RootsTable({ roots, counts }) {
  return <div class="fc-tablewrap"><table class="fc-table fc-table--comfortable fc-table--wrap fc-roots">
    <caption class="fc-sr">Folders this library was built from</caption>
    <thead><tr><th scope="col">Folder</th><th scope="col" class="fc-num">Runs found</th><th scope="col">Last scanned</th><th scope="col">State</th>
      <th scope="col" class="fc-num">Full detail</th><th scope="col" class="fc-num">Quick detail</th><th scope="col" class="fc-num">Raw data only</th></tr></thead>
    <tbody>{roots.map(r => { const d = r.detail || {}, off = r.state && r.state !== 'reachable';
      return <tr>
        <th scope="row"><span class="fc-mono fc-break">{r.path}</span><span class="fc-table__sub">{r.id}{r.live ? ' · the sequencer writes here' : ''}</span></th>
        <td class="fc-num">{counts ? F.int(counts.get(r.id) || 0) + ' of ' + F.int(r.runs) : F.int(r.runs)}</td><td class="fc-nowrap">{r.last_scanned ? <When iso={r.last_scanned} time /> : 'not yet'}</td>
        <td>{off ? <span class="fc-notice">{r.state === 'offline' ? <>Offline{r.offline_since ? <> since <When iso={r.offline_since} time /></> : null}</> : (ROOT_STATE_WORD[r.state] || r.state)}</span> : 'Reachable'}</td>
        <td class="fc-num">{F.int(d.full || 0)}</td><td class="fc-num">{F.int(d.quick || 0)}</td><td class="fc-num">{F.int(d.raw_only || 0)}</td></tr>; })}</tbody></table></div>;
}

function Activity({ a }) {
  const n = a.months.length; if (!n) return null;
  const xs = a.months.map((_, i) => i + 0.5), color = tk('--fc-this-run');
  // month labels: every month of a short span; quarters, half-years or Januaries of a longer one. The year is said under each January.
  const step = n <= 12 ? 1 : n <= 24 ? 3 : n <= 60 ? 6 : 12;
  const ticks = a.months.map((m, i) => ({ v: i + 0.5, label: (+m.slice(5) - 1) % step === 0 ? (m.slice(5) === '01' || (n <= 12 && i === 0) ? F.monthLabel(m) : F.monthLabel(m, true)) : '', major: true }));
  const spec = (ys, unit, fmt) => ({ form: 'xy', x: { domain: [0, n], kind: 'linear', ticks }, y: { title: unit }, series: [{ name: unit, type: 'vbar', x: xs, parts: [{ y: ys, color }] }], cursor: { x: xs },
    readout: i => [[F.monthLabel(a.months[i]), ''], [fmt(ys[i], i), '']] });
  const gb = a.estimated_bases.map(v => (v || 0) / 1e9);
  const total = a.runs.reduce((t, v) => t + v, 0), busiest = a.runs.indexOf(Math.max(...a.runs));
  return <section class="fc-home__block fc-home__block--activity" aria-labelledby="h-activity">
    <div class="fc-section__head"><h2 class="fc-section__title" id="h-activity">Activity</h2><span class="fc-secondary">by the month a run started</span></div>
    <div class="fc-grid fc-grid--2">
      <ChartPanel title="Runs per month" reading={F.plural(total, 'run') + ' over ' + F.plural(n, 'month') + '; the most in ' + F.monthLabel(a.months[busiest]) + ' (' + a.runs[busiest] + ').'} height={150}
        spec={spec(a.runs, 'Runs', v => F.plural(v, 'run'))} table={() => ({ head: ['Month', 'Runs'], rows: a.months.map((m, i) => [F.monthLabel(m), F.int(a.runs[i])]) })} />
      <ChartPanel title="Estimated output per month" reading="Estimated bases of the runs that started in each month: MinKNOW's estimate from the signal, before basecalling." height={150}
        spec={spec(gb, 'Gb, estimated', (v, i) => F.basesBig(a.estimated_bases[i] || 0) + ' estimated')}
        table={() => ({ head: ['Month', 'Estimated bases', 'Reads'], rows: a.months.map((m, i) => [F.monthLabel(m), F.basesBig(a.estimated_bases[i] || 0), F.count((a.reads || [])[i] || 0)]) })} />
    </div>
  </section>;
}

export function Home() {
  const ix = index.value, lib = ix.library || {};
  if (!ix.n) return <EmptyLibrary lib={lib} roots={ix.roots} />;
  const c = lib.counts || {}, rows = ix.rows.slice().sort((a, b) => b.startedMs - a.startedMs);
  const attention = rows.filter(needsLook), nFail = attention.filter(r => r.verdict === 'fail').length;
  // What still waits comes from the index's own progress (one count, the one the scan prints): runs whose reads
  // are not counted yet and read-data files that are not read yet. A run whose read table is cut off does not wait.
  const pg = ix.progress || {}, t1 = pg.tier1 || {}, rdp = pg.read_data || {}, scanning = pg.state === 'scanning';
  const isCut = r => r.modifiers.includes('reads_partly_counted') || r.modifiers.includes('reads_not_counted');
  const counting = t1.runs_waiting != null ? t1.runs_waiting : rows.filter(r => r.detail_t1 === 'quick' && r.pending > 0 && r.run_state !== 'running' && !isCut(r)).length;
  const cut = rows.filter(isCut).length, running = rows.filter(r => r.run_state === 'running').length;
  const waiting = rdp.files_waiting != null ? rdp.files_waiting : 0;
  const scanCmd = 'flowcellar scan -L ' + (lib.path ? quote(lib.path) : '<this library>');
  const conflict = runsInConflict(ix), quiet = quietNotices(ix);
  const nFc = c.flow_cells != null ? c.flow_cells : ix.byFc.size, projects = ix.projects, named = projects.filter(p => p.id !== '_none');
  const nProj = c.projects != null ? c.projects : named.length, bases = lib.estimated_bases != null ? lib.estimated_bases : rows.reduce((t, r) => t + (r.gb || 0) * 1e9, 0);
  const offline = ix.roots.filter(r => r.state && r.state !== 'reachable').length;
  return <div class="fc-page fc-home">
    <header class="fc-home__head">
      <div class="fc-home__title">
        <h1 class="fc-home__name">{lib.name || 'FlowCellar library'}</h1>
        <p class="fc-home__line">{F.plural(ix.n, 'run')} on {F.plural(nFc, 'flow cell')}{nProj ? ' in ' + F.plural(nProj, 'project') : ''} {'·'} {F.basesBig(bases)} estimated {'·'} <Updated /></p>
      </div>
      <SearchBox />
    </header>
    <LibraryNotices cls="fc-home__notices" noProgress />
    <div class="fc-grid fc-grid--2 fc-home__top">
      <section class="fc-panel" aria-labelledby="h-attention">
        <div class="fc-panel__head"><h2 class="fc-section__title" id="h-attention">Needs attention</h2>
          {attention.length ? <a class="fc-home__all" href="#/runs?q=status:warn,fail">All {attention.length} with a warning{nFail ? ' or a fail' : ''}</a> : null}</div>
        <div class="fc-panel__body">
          {attention.length ? <ul class="fc-homerows">{attention.slice(0, ATTENTION).map(r => <AttentionRow r={r} />)}</ul>
            : <p class="fc-home__calm"><Glyph v="pass" label="Pass" /> Nothing needs attention: no run has a warning or a failed check.</p>}
          {counting || waiting ? <div class="fc-home__aside fc-home__progress" role="status">
            <p class="fc-notice">{scanning ? <b>A scan is running. </b> : <b>Waiting for the next scan. </b>}{pg.message || [counting ? 'The reads of ' + F.plural(counting, 'run') + ' are not counted yet.' : '', waiting ? 'Read data: ' + F.plural(waiting, 'file') + (waiting === 1 ? ' waits' : ' wait') + ' to be read.' : ''].filter(Boolean).join(' ')}</p>
            <p>{counting ? <>{F.plural(counting, 'run')} {counting === 1 ? 'reads' : 'read'} “so far” until {counting === 1 ? 'its reads are' : 'their reads are'} counted (<a href="#/runs?q=is:counting+-is:running">show {counting === 1 ? 'it' : 'them'}</a>). </> : null}{waiting ? <>Read data of {F.plural(rows.filter(r => r.t2_files_waiting > 0).length, 'run')} is a verdict so far (<a href="#/runs?q=is:waiting">show</a>). </> : null}</p>
            {scanning ? <p class="fc-secondary">This page looks for news once a minute.</p> : <p>No scan is running. To go on: <span class="fc-mono fc-break">{scanCmd}</span></p>}
          </div> : null}
          {cut ? <p class="fc-home__aside fc-notice">{F.plural(cut, 'run')} {cut === 1 ? 'has' : 'have'} a read table that is cut off or cannot be read: {cut === 1 ? 'its' : 'their'} reads were counted as far as the file goes, and another scan will not change that. <a href="#/runs?q=is:cut">Show {cut === 1 ? 'it' : 'them'}</a></p> : null}
          {running ? <p class="fc-home__aside fc-notice">{F.plural(running, 'run')} {running === 1 ? 'is' : 'are'} sequencing now. <a href="#/runs?q=is:running">Show {running === 1 ? 'it' : 'them'}</a></p> : null}
        </div>
      </section>
      <section class="fc-panel" aria-labelledby="h-latest">
        <div class="fc-panel__head"><h2 class="fc-section__title" id="h-latest">Latest runs</h2><span class="fc-secondary">output in the first 48 h</span><a class="fc-home__all" href="#/runs">All {ix.n} runs</a></div>
        <div class="fc-panel__body"><ul class="fc-homerows">{rows.slice(0, LATEST).map(r => <LatestRow r={r} />)}</ul></div>
      </section>
    </div>
    <section class="fc-panel fc-home__block fc-home__block--projects" aria-labelledby="h-projects">
      <div class="fc-panel__head"><h2 class="fc-section__title" id="h-projects">Projects</h2><span class="fc-secondary">named from the folders the runs were found in; a project found in two places is one project</span>
        {projects.length > PROJECTS_ON_HOME ? <a class="fc-home__all" href="#/projects">All {named.length} projects</a> : null}</div>
      {projects.length ? <ProjectsTable projects={projects} roots={ix.rootById} limit={PROJECTS_ON_HOME} sortable /> : <p class="fc-panel__body fc-secondary">No project folder was recognised: every run sits directly under a scanned folder.</p>}
    </section>
    <section class="fc-panel fc-home__block fc-home__block--roots" aria-labelledby="h-roots">
      <div class="fc-panel__head"><h2 class="fc-section__title" id="h-roots">Where this library comes from</h2>
        <span class="fc-secondary">{F.plural(ix.roots.length, 'folder')} scanned{offline ? ' · ' + offline + ' not reachable at the last scan' : ''}</span></div>
      {ix.roots.length ? <RootsTable roots={ix.roots} /> : <p class="fc-panel__body fc-secondary">The index names no scanned folder: these runs were added one by one.</p>}
      {quiet.length ? <div class="fc-panel__body fc-home__quiet" role="status">{quiet.map(n => <p class="fc-notice">{n.message}{n.code === 'two_candidate_sources' ? <> {conflict.length ? <><a href="#/runs?q=is:conflict">The {F.plural(conflict.length, 'run')} this concerns</a>. </> : null}
        Each run’s page names its two files under Details, Data sources.</> : null}</p>)}</div> : null}
      <div class="fc-panel__body"><AddCommands path={lib.path} /></div>
    </section>
    <HomeAnalyses />
    <Activity a={ix.activity} />
  </div>;
}
