// The Projects tab (every project of the library) and the page of one project: its numbers, where it is,
// and the Runs table filtered to it. Both are drawn from the index alone.
// On paper ("Save as PDF") the page of a project is a report of the project: its numbers, every run of it in the
// order they were sequenced with what each one gave, why a run was not a pass, and the single-cell and spatial
// analyses made from its runs. The folders it lies in and the Runs table of the screen are left out.
import { index } from '../core/data.js';
import { route, hrefProject } from '../core/state.js';
import { Badge, Chip, Crumbs, Tile, Tally, Updated, When, LibraryNotices, quietNotices, runsInConflict, runSoFar } from '../shell/parts.jsx';
import { SavePdf, PaperHead, PaperScope } from '../shell/paper.jsx';
import { printing } from '../core/print.js';
import * as S from '../analysis/sc.js';
import { ProjectsTable } from '../home/Home.jsx';
import { RunsView } from '../hub/Hub.jsx';
import * as F from '../core/format.js';
import { ROOT_STATE_WORD, ASSAY_WORD, VERDICT_WORD } from '../core/labels.js';

export function Projects() {
  const ix = index.value, named = ix.projects.filter(p => p.id !== '_none'), none = ix.projects.find(p => p.id === '_none');
  const multi = ix.rows.filter(r => r.project_ids.length > 1).length, sum = ix.projects.reduce((t, p) => t + (p.runs || 0), 0);
  return <div class="fc-page fc-projectspage">
    <header class="fc-home__head"><div class="fc-home__title"><h1 class="fc-home__name">Projects</h1>
      <p class="fc-home__line">{F.plural(named.length, 'project')} in {ix.library.name || 'this library'}{none ? ' · ' + F.plural(none.runs, 'run') + ' without a project' : ''} {'·'} <Updated /></p></div></header>
    <LibraryNotices cls="fc-home__notices" />
    {ix.projects.length ? <section class="fc-panel fc-home__block" aria-label="Projects"><ProjectsTable projects={ix.projects} roots={ix.rootById} sortable /></section>
      : <div class="fc-empty fc-page__empty"><b>No project yet.</b><span>Projects appear when runs have been scanned.</span></div>}
    {multi ? <p class="fc-section__note fc-notice">{F.plural(multi, 'run')} {multi === 1 ? 'carries' : 'carry'} samples of two projects and {multi === 1 ? 'is' : 'are'} counted under each of them: the rows above add up to {F.int(sum)} runs, the library holds {F.int(ix.n)}. <a href="#/runs?group=project">Runs grouped by project</a></p> : null}
    <p class="fc-section__note">A project is named from the folder its runs were found in (the first folder under each scanned folder, unless the library's configuration says otherwise). The same name under two scanned folders is one project, and a flow cell that carries samples of two projects is listed under both.</p>
  </div>;
}

// The runs of a project on paper: oldest first, one row each, the widest output as the full bar.
function PaperRuns({ rows, analyses }) {
  const runs = [...rows].sort((a, b) => String(a.started || '').localeCompare(String(b.started || ''))), max = Math.max(1e-9, ...runs.map(r => r.gb || 0));
  const why = runs.filter(r => r.verdict !== 'pass');
  const flow = a => a.links.map(l => (l.fc || l.run || 'flow cell not known') + (l.share != null && a.links.length > 1 ? ' ' + F.pct(l.share) : '') + (l.resolved ? '' : ' (run not in library)')).join(', ');
  return <>
    <section class="fc-section fc-paperproject" aria-labelledby="h-paper-runs"><div class="fc-section__head"><h2 class="fc-section__title" id="h-paper-runs">Runs of this project</h2><span class="fc-secondary">oldest first; output is estimated from the signal, over the whole run</span></div>
      <div class="fc-tablewrap"><table class="fc-table fc-table--wrap fc-paperruns"><caption class="fc-sr">Runs of this project</caption>
        <thead><tr><th scope="col">Started</th><th scope="col">Flow cell {'·'} run name</th><th scope="col">Kind</th><th scope="col">Verdict</th><th scope="col" class="fc-num">Ran</th><th scope="col">Output</th><th scope="col" class="fc-num">Reads</th><th scope="col" class="fc-num">Pores<span class="fc-table__sub">first scan</span></th><th scope="col" class="fc-num">Occupancy<span class="fc-table__sub">first hour</span></th><th scope="col" class="fc-num">Read N50</th></tr></thead>
        <tbody>{runs.map(r => <tr>
          <td class="fc-nowrap">{F.date(r.started, r.stub)}</td>
          <th scope="row"><span class="fc-nowrap"><span class="fc-mono">{r.fc}</span>{r.use > 1 ? <span class="fc-muted"> use {r.use}</span> : null}</span><span class="fc-table__sub">{r.run_name || r.id}</span></th>
          <td class="fc-nowrap">{ASSAY_WORD[r.assay] || ASSAY_WORD[r.assay_group] || F.DASH}<span class="fc-table__sub">{r.kit || ''}</span></td>
          <td class="fc-nowrap"><Badge v={r.verdict} soFar={runSoFar(r)} /></td>
          <td class="fc-num fc-nowrap">{F.hours(r.hours)}</td>
          <td class="fc-nowrap">{r.gb != null ? <><span class="fc-cellbar" style={{ width: Math.max(1, Math.round(r.gb / max * 44)) + 'px' }} />{F.bases(r.gb * 1e9)}</> : F.DASH}</td>
          <td class="fc-num fc-nowrap">{F.count(r.reads)}</td><td class="fc-num">{F.int(r.pores)}</td><td class="fc-num">{F.pct(r.occ_first)}</td>
          <td class="fc-num fc-nowrap">{r.t2_n50 != null ? F.bases(r.t2_n50) : r.n50_est != null ? <>{F.bases(r.n50_est)}<span class="fc-table__sub">estimated</span></> : F.DASH}</td></tr>)}</tbody></table></div>
      {why.length ? <ol class="fc-fcruns__reasons">{why.map(r => <li><b><span class="fc-mono">{r.fc}</span> {r.run_name || ''} {'·'} {VERDICT_WORD[r.verdict] || r.verdict}.</b> {[r.headline, ...(r.reasons || [])].filter(Boolean).filter((x, k, all) => all.indexOf(x) === k).join(' ')}</li>)}</ol>
        : <p class="fc-section__note">Every run of this project passed.</p>}
    </section>
    {analyses.length ? <section class="fc-section fc-paperproject" aria-labelledby="h-paper-an"><div class="fc-section__head"><h2 class="fc-section__title" id="h-paper-an">Single-cell and spatial analyses</h2><span class="fc-secondary">made from runs of this project; an analysis of several flow cells is one analysis of their merged reads</span></div>
      <div class="fc-tablewrap"><table class="fc-table fc-table--wrap fc-paperruns"><caption class="fc-sr">Analyses made from runs of this project</caption>
        <thead><tr><th scope="col">Analysis</th><th scope="col">Kind</th><th scope="col">Pipeline</th><th scope="col">Verdict</th><th scope="col">Reads from</th><th scope="col">In a sentence</th></tr></thead>
        <tbody>{analyses.map(a => <tr><th scope="row">{a.label || a.id}</th><td>{S.ASSAY_WORD[a.assay] || F.DASH}</td><td class="fc-nowrap">{(S.PIPELINE_WORD[a.pipeline] || a.pipeline) + (a.version ? ' ' + a.version : '')}</td>
          <td class="fc-nowrap"><Badge v={a.verdict} /></td><td>{flow(a)}</td><td>{a.headline || F.DASH}</td></tr>)}</tbody></table></div></section> : null}
  </>;
}

export function Project() {
  const ix = index.value, id = route.value.id, p = ix.projectById.get(id);
  const crumbs = <Crumbs><a href="#/home">Library</a> {'›'} <a href="#/projects">Projects</a> {'›'} <b>{p ? (p.name == null ? 'No project' : p.name) : id}</b></Crumbs>;
  if (!p) return <div class="fc-page">{crumbs}<div class="fc-empty fc-page__empty"><b>This project is not in the library.</b><span>No project has the id {id}. <a href="#/projects">All projects</a></span></div></div>;
  const inProject = id === '_none' ? (r => !r.project_ids.length) : (r => r.project_ids.includes(id));
  const rows = ix.rows.filter(inProject), perRoot = new Map();
  for (const r of rows) for (const k of r.roots) perRoot.set(k, (perRoot.get(k) || 0) + 1);
  const rootIds = [...new Set([...(p.roots || []), ...perRoot.keys()])], v = p.verdicts || {};
  const look = (v.warn || 0) + (v.fail || 0), reads = p.reads != null ? p.reads : rows.reduce((t, r) => t + (r.reads || 0), 0);
  // one sentence on the verdicts of the project: what needs a look, and what could not be judged
  const nInc = v.incomplete || 0, nRun = v.running || 0, shared = rows.filter(r => r.project_ids.length > 1).length;
  const verdictLine = (look ? F.plural(look, 'run') + (look === 1 ? ' needs' : ' need') + ' a look' : 'No run of this project has a warning or a failed check')
    + (nInc ? '; ' + F.plural(nInc, 'run') + (nInc === 1 ? ' is' : ' are') + ' incomplete (too short to judge, or raw data only)' : '') + (nRun ? '; ' + nRun + ' sequencing now' : '') + '.';
  const paper = printing.value, name = p.name == null ? 'No project' : p.name;
  const made = [{ provenance: { updated: ix.written, generator: (ix.raw.meta || {}).generator } }];       // the page is drawn from the index: when it was written, and by which version
  const mine = new Set(rows.map(r => r.id)), conflicts = runsInConflict(ix).filter(r => mine.has(r.id)), note = conflicts.length ? quietNotices(ix).find(n => n.code === 'two_candidate_sources') : null;
  return <div class="fc-projectpage">
    <div class="fc-page fc-projectpage__top">
      <PaperHead kind="Project report" bundles={made} file={'FlowCellar project ' + name} />
      {crumbs}
      <header class="fc-runhead"><div class="fc-runhead__row">
        <h1 class="fc-runhead__id">{p.name == null ? 'No project' : p.name}</h1>
        <span class="fc-runhead__chips"><Chip>{F.plural(p.runs, 'run')}</Chip><Chip>{F.plural(p.flow_cells, 'flow cell')}</Chip>{p.samples ? <Chip>{F.plural(p.samples, 'sample')}</Chip> : null}<Chip>{F.period(p.first, p.last)}</Chip></span>
        <span class="fc-runhead__actions"><SavePdf /></span></div>
        <p class="fc-runhead__line">{p.name == null ? 'Runs that sit directly under a scanned folder, so no project folder names them.' : 'Named from the folder its runs were found in.'} {verdictLine}{shared ? ' ' + F.plural(shared, 'run') + (shared === 1 ? ' carries' : ' carry') + ' samples of another project too and ' + (shared === 1 ? 'is' : 'are') + ' counted in both.' : ''}</p>
      </header>
      {note ? <div class="fc-libnotices fc-home__notices" role="status"><p class="fc-notice">Read data of {F.plural(conflicts.length, 'run')} of this project holds files that are not counted, because two candidate sources describe the same reads ({conflicts.slice(0, 4).map(r => r.fc).join(', ')}{conflicts.length > 4 ? ', …' : ''}). To choose, put <span class="fc-mono">primary = true</span> on the <span class="fc-mono">[[inputs]]</span> entry of the files to count, in the library’s <span class="fc-mono">flowcellar.toml</span>, and scan again.</p></div> : null}
      <div class="fc-keynumbers"><div class="fc-tiles">
        <Tile label="Runs" value={F.int(p.runs)} caption={F.plural(p.flow_cells, 'flow cell')} />
        <Tile label="Samples" value={p.samples ? F.int(p.samples) : F.DASH} state={p.samples ? null : 'unknown'} caption={p.samples ? 'named in the read data' : 'known when read data is linked'} />
        <Tile label="Estimated output" value={F.basesParts(p.estimated_bases)[0]} unit={F.basesParts(p.estimated_bases)[1]} caption="from signal, all runs added up" />
        <Tile label="Reads" value={F.countParts(reads)[0]} unit={F.countParts(reads)[1]} caption="sequenced, all runs added up" />
        <Tile label="First run" value={F.date(p.first)} title={F.timeTitle(p.first)} caption={p.last && F.date(p.last) !== F.date(p.first) ? 'last run ' + F.date(p.last) : 'one day'} />
        <div class="fc-tile"><span class="fc-tile__label"><span>Verdicts</span></span><span class="fc-tile__value fc-tile__value--tally"><Tally counts={v} /></span><span class="fc-tile__caption">one per run; the Runs table below says why</span></div>
      </div></div>
      <PaperScope items={[['Runs', 'all ' + F.plural(rows.length, 'run') + ' of the project, oldest first'], ['Left out', 'the folders the project lies in, and the search and filters of the Runs table on screen']]} />
      {paper ? <PaperRuns rows={rows} analyses={ix.analyses.filter(a => a.links.some(l => l.resolved && mine.has(l.run)))} /> : null}
      <section class="fc-panel fc-home__block fc-projectpage__where" aria-labelledby="h-where">
        <div class="fc-panel__head"><h2 class="fc-section__title" id="h-where">Where it is</h2><span class="fc-secondary">{rootIds.length > 1 ? 'found under ' + rootIds.length + ' scanned folders: one project' : 'the scanned folder its runs were found under'}</span></div>
        {rootIds.length ? <div class="fc-tablewrap"><table class="fc-table fc-table--comfortable fc-table--wrap"><caption class="fc-sr">Scanned folders that hold this project</caption>
          <thead><tr><th scope="col">Folder</th><th scope="col" class="fc-num">Runs of this project</th><th scope="col">State</th><th scope="col">Last scanned</th></tr></thead>
          <tbody>{rootIds.map(k => { const r = ix.rootById.get(k) || { id: k };
            return <tr><th scope="row"><span class="fc-mono fc-break">{r.path || k}</span><span class="fc-table__sub">{k}</span></th><td class="fc-num">{perRoot.has(k) ? F.int(perRoot.get(k)) : <span class="fc-muted" title="Read-data files of this project were found here">files only</span>}</td>
              <td>{r.state && r.state !== 'reachable' ? <span class="fc-notice">{ROOT_STATE_WORD[r.state] || r.state}{r.offline_since ? <> since <When iso={r.offline_since} time /></> : null}</span> : 'Reachable'}</td><td class="fc-nowrap">{r.last_scanned ? <When iso={r.last_scanned} time /> : F.DASH}</td></tr>; })}</tbody></table></div>
          : <p class="fc-panel__body fc-secondary">The index names no folder for this project.</p>}
      </section>
      <h2 class="fc-section__title fc-projectpage__runs" id="h-project-runs">Runs of this project</h2>
    </div>
    <RunsView key={id} rowFilter={inProject} embedded label={'Runs of ' + (p.name == null ? 'no project' : p.name)} />
  </div>;
}
export { hrefProject };
