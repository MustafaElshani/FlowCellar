// The run view: header, two-part verdict with reasons that jump to their evidence, six key numbers,
// the controls bar that drives every chart, a rail that is a status column, and the sections in the order people
// ask their questions. The same view is the one-run report.
import { useEffect, useRef, useState } from 'preact/hooks';
import { useBundle, bundleState, bundleError, retryBundle, rowById, isReport, index } from '../core/data.js';
import { route, setParams, param, selected, toggleSelected, makeHash, hrefFlowCell, hrefProject, refocus } from '../core/state.js';
import { ruleState, fmtValue } from './section.jsx';
import { controls, sectionResult, ruleById, partialInfo, readSetWords } from './common.js';
import { Badge, Glyph, Chip, Mid, Tile, Seg, Crumbs, When, basesLabel } from '../shell/parts.jsx';
import * as F from '../core/format.js';
import * as D from '../core/derive.js';
import { ASSAY_WORD, DEVICE_WORD, END_WORD, SECTIONS, checkName, sectionOf, label, ROLE_WORD } from '../core/labels.js';
import { RunTimeline, FlowCellSection, ReadsThroughPore } from './health.jsx';
import { ReadDataCard, ReadFate, SamplesSection, LengthQuality, Alignment, CdnaLibrary } from './readdata.jsx';
import { PychopperSection, PYCHOPPER } from './pychopper.jsx';
import { RnaComposition, BIOTYPE } from './rnacomp.jsx';
import { ReadsOverTime, Transcripts } from './alignment.jsx';
import { Details } from './details.jsx';
import { RunAnalyses, analysesOfRun } from '../analysis/Analyses.jsx';
import { SavePdf, PaperHead, PaperScope } from '../shell/paper.jsx';

const hub = href => isReport ? undefined : href;        // in a report, names of things outside the file are plain text
const lower = s => s ? s.charAt(0).toLowerCase() + s.slice(1) : s;

function limitText(r) {
  const t = r.threshold || {}, f = v => fmtValue(r.metric || r.id, v), bl = t.baseline || {};
  if ((t.applied === 'baseline' || t.kind === 'baseline') && (bl.warn_level != null || bl.fail_level != null))
    return 'this lab’s warning level: ' + f(bl.warn_level) + (bl.fail_level != null ? '; fail level ' + f(bl.fail_level) : '') + (t.warn != null ? ' (fixed level ' + f(t.warn) + ')' : '');
  if (t.kind === 'min' || (t.kind === 'baseline' && t.warn != null)) return 'warn below ' + f(t.warn) + (t.fail != null ? ', fail below ' + f(t.fail) : '');
  if (t.kind === 'max') return 'warn above ' + f(t.warn) + (t.fail != null ? ', fail above ' + f(t.fail) : '');
  if (t.kind === 'range' && t.pass_band) return 'pass from ' + t.pass_band.map(f).join(' to ');
  if (t.kind === 'run_config' && t.pass_band) return 'target ' + t.pass_band.map(f).join(' to ');
  return 'a rule:';
}
const evHref = (b, r) => makeHash(route.value.path, { ...route.value.q, ev: r.id });
// The check on the full-length share is shown in the Pychopper section when the bundle holds the whole pychopper table
// (one section for the cDNA step, not two that count the same reads): its evidence lives there, under that name.
const secOf = (b, r) => { const sec = sectionOf(r); return sec === 'cdna' && D.mod(b, PYCHOPPER) ? 'pychopper' : sec; };
const secName = (b, r) => { const sec = secOf(b, r); return sec === 'pychopper' ? 'Pychopper' : SECTIONS[sec]; };

function Reason({ b, r }) {
  const st = ruleState(r), name = checkName(r), ci = (r.message || '').indexOf(':'), starts = ci > 0 && ci < 64 && !/^Not assessed/.test(r.message);
  return <li class={st === 'unknown' ? 'is-na' : ''}><Glyph v={st} label={st === 'unknown' ? 'Not assessed' : undefined} /><span>{starts ? null : <b>{name}. </b>}{r.message}</span>
    <a class="fc-reasons__go" href={evHref(b, r)}>{secName(b, r) || 'Evidence'}</a></li>;
}

// The table behind a verdict. It opens across the whole summary card and is as tall as its rows (no box that
// scrolls inside the page); where a limit comes from stands under the limit, so that no column is a sliver.
function AllChecks({ b, rules, title, tier }) {
  return <div class="fc-allchecks" id={'allchecks-' + tier} role="group" aria-label={'All checks of ' + title}>
    <div class="fc-allchecks__head"><b>All checks: {title}</b><span class="fc-secondary">{F.plural(rules.length, 'check')}{tier === 'run_health' && (b.status.pending || []).length ? ' · ' + b.status.pending.length + ' more wait for the read count' : ''}</span></div>
    <div class="fc-tablewrap"><table class="fc-table fc-table--wrap fc-allchecks__table"><caption class="fc-sr">All checks of {title}</caption>
    <thead><tr><th scope="col">Result</th><th scope="col">Check</th><th scope="col" class="fc-num">Value</th><th scope="col">Limit, and where it comes from</th><th scope="col">Can</th><th scope="col">Evidence</th></tr></thead>
    <tbody>{rules.map(r => { const st = ruleState(r);
      return <tr><td class="fc-nowrap"><Glyph v={st} /> {st === 'pass' ? 'Passed' : st === 'warn' ? 'Warning' : st === 'fail' ? 'Fail' : 'Not assessed'}</td><th scope="row" title={r.id}>{checkName(r)}</th>
        <td class="fc-num fc-nowrap">{st === 'unknown' ? F.DASH : fmtValue(r.metric || r.id, r.value)}</td><td class="fc-allchecks__limit">{limitText(r)}{r.threshold && r.threshold.basis ? <span class="fc-allchecks__from">{r.threshold.basis}</span> : null}</td><td class="fc-nowrap">{ROLE_WORD[r.role] || r.role}</td>
        <td class="fc-nowrap"><a href={evHref(b, r)}>{secName(b, r)}</a></td></tr>; })}
      {tier === 'run_health' ? (b.status.pending || []).map(id => <tr class="is-dim"><td class="fc-nowrap"><Glyph v="unknown" label="Waiting" /> Waiting</td><th scope="row" title={id}>{label(id)}</th><td class="fc-num">{F.DASH}</td><td colSpan={3}>waits for the read count</td></tr>) : null}
    </tbody></table></div></div>;
}
const rulesOf = (b, tier) => ((b.status || {}).rules || []).filter(r => r.tier === tier);
const filesWaiting = b => ((b.provenance || {}).inputs || []).filter(i => (i.notes || []).some(n => n.code === 'waiting')).length;

function Part({ b, tier, title, open, onToggle }) {
  const st = b.status, v = (st.by_tier || {})[tier] || 'not_available', rules = rulesOf(b, tier), pend = tier === 'run_health' ? (st.pending || []) : [];
  const waitingFiles = tier === 'read_level' ? filesWaiting(b) : 0, judged = rules.filter(r => r.role !== 'info');
  const trig = rules.filter(r => r.outcome === 'triggered' && (r.severity === 'warn' || r.severity === 'fail')).sort((a, b2) => (a.severity === 'fail' ? 0 : 1) - (b2.severity === 'fail' ? 0 : 1) || (a.role === 'key' ? 0 : 1) - (b2.role === 'key' ? 0 : 1));
  const ok = rules.filter(r => r.outcome !== 'not_evaluated' && !trig.includes(r)), na = rules.filter(r => r.outcome === 'not_evaluated');
  // nothing triggered: up to three passed checks that can fail a run, so a Pass shows the numbers it rests on
  const pick = ['fc_pores_first_scan', 'fc_occupancy_first_hour'].map(id => ok.find(r => r.id === id)).concat(ok.filter(r => /^yield_/.test(r.id) && r.role === 'key'), ok.filter(r => r.role === 'key' && r.tier === 'read_level')).filter(Boolean);
  const shown = trig.length ? trig.slice(0, 3) : [...new Set(pick)].slice(0, 3), rest = trig.length > 3 ? trig.length - 3 : 0;
  const nWarn = trig.filter(r => ruleState(r) === 'warn').length, nFail = trig.length - nWarn;
  const count = !rules.length ? (v === 'running' ? 'checks are made as the numbers arrive' : 'no check was applied: the run is too short to judge') : pend.length ? `${rules.length} checks done, ${pend.length} wait for the read count`
    : waitingFiles ? `${judged.length - na.filter(r => r.role !== 'info').length} of ${judged.length} checks assessed · ${F.plural(waitingFiles, 'file')} of this run not read yet`
    : !trig.length && !na.length ? `${ok.length} of ${rules.length} checks passed`
    : [ok.length + ' passed', nFail ? nFail + ' failed' : null, nWarn ? F.plural(nWarn, 'warning') : null, na.length ? na.length + ' not assessed' : null].filter(Boolean).join(' · ');
  const tiers = b.tiers || {}, absent = tier === 'read_level' && (!tiers.read_level || tiers.read_level.state === 'absent'), stub = tier === 'run_health' && D.flags(b).stub;
  // a gap in read fate (a stage without a source) is said here only when a check could not be assessed because of it
  const fateGap = rules.some(r => (r.module === 'funnel' || r.module === 'funnel.samples') && r.outcome === 'not_evaluated') ? partialInfo(b, 'funnel') : null;
  const pr = b.identity.protocol_run, part = tier === 'read_level' && !absent ? partialInfo(b, 'len') || partialInfo(b, 'yield.basecalled') || fateGap : null;
  return <div class="fc-verdict__part">
    <div class="fc-verdict__parthead"><Badge v={v} soFar={tier === 'run_health' ? st.provisional : waitingFiles > 0} /><b>{title}</b>
      {!absent && !stub ? <span class="fc-secondary">{count}</span> : null}
      {rules.length ? <button type="button" class="fc-btn fc-btn--ghost fc-verdict__all" data-fk={'allchecks-' + tier} aria-expanded={open} aria-controls={'allchecks-' + tier} onClick={() => { onToggle(); refocus('[data-fk="allchecks-' + tier + '"]'); }}>All checks</button> : null}</div>
    {absent ? <div class="fc-verdict__absent"><p>No basecalled reads are linked to this run. FlowCellar looks for basecalled or aligned read files that name flow cell {b.identity.flow_cell.id} and run {pr.short_id} inside them.</p>
      <p class="fc-secondary">{b.identity.assay.barcoded ? 'The kit carries barcodes; MinKNOW recorded no sample names. They arrive with the read data. ' : ''}<a href="#s-read_data" onClick={e => { e.preventDefault(); document.getElementById('s-read_data')?.scrollIntoView(); }}>What will appear</a></p></div>
      : stub ? <p class="fc-secondary">Only raw data was found in this folder: there is no MinKNOW report to judge the flow cell from.</p>
      : <ul class="fc-reasons">{shown.map(r => <Reason b={b} r={r} />)}{(trig.length ? [] : na).slice(0, Math.max(0, 3 - shown.length)).concat(trig.length ? na.slice(0, 1) : []).map(r => <Reason b={b} r={r} />)}
        {rest ? <li><span /><span class="fc-secondary">+ {rest} more in All checks</span></li> : null}</ul>}
    {part ? <p class="fc-verdict__note fc-notice">{part.text}</p> : null}
    {pend.length ? <p class="fc-verdict__note">Waiting for the read count: {pend.map(id => lower(label(id))).join(', ')}.</p> : null}
  </div>;
}

function Verdict({ b }) {
  const st = b.status, [open, setOpen] = useState('');         // which part shows all its checks: '', run_health or read_level
  const toggle = tier => setOpen(open === tier ? '' : tier);
  return <section class="fc-verdict" id="s-summary" data-sec="summary" aria-label="Summary">
    <div class="fc-verdict__head"><Badge lg v={st.verdict} soFar={st.provisional} />
      <div><p class="fc-verdict__headline">{st.headline}</p>{st.detail ? <p class="fc-verdict__detail">{st.detail}</p> : null}</div></div>
    {(st.modifiers || []).filter(m => m.code !== 'read_data_so_far').map(m => <p class="fc-verdict__note fc-notice">{m.message}{m.code === 'quick_detail' && !isReport ? ' Reload this page then.' : ''}{m.code === 'reused_flow_cell' && !isReport ? <> <a href={hrefFlowCell(b.identity.flow_cell.id)}>History of this flow cell</a></> : null}</p>)}
    <div class="fc-verdict__parts"><Part b={b} tier="run_health" title="Run health" open={open === 'run_health'} onToggle={() => toggle('run_health')} /><Part b={b} tier="read_level" title="Read data" open={open === 'read_level'} onToggle={() => toggle('read_level')} /></div>
    {open ? <AllChecks b={b} rules={rulesOf(b, open)} tier={open} title={open === 'run_health' ? 'Run health' : 'Read data'} /> : null}
  </section>;
}

// ---------- key numbers ----------
function stripFor(value, band, limits, lo0) {
  const vals = [value].concat(band || [], limits || []).filter(F.isNum); if (!vals.length) return null;
  let lo = Math.min(...vals), hi = Math.max(...vals); const pad = (hi - lo || Math.abs(hi) || 1) * 0.18; lo = lo0 != null ? lo0 : Math.max(0, lo - pad); hi += pad;
  if (!band && !(limits || []).some(F.isNum)) return null;
  return { lo, hi, band, mark: value, limits: (limits || []).filter(F.isNum) };
}
function KeyNumbers({ b, c }) {
  const st = b.status, y = D.summary(b, 'yield'), fc = D.summary(b, 'fc'), le = D.summary(b, 'len.est'), fl = D.flags(b), tiles = [];
  const jump = id => { const r = ruleById(b, id); return r ? evHref(b, r) : undefined; };
  const g = id => { const r = ruleById(b, id); return r ? (ruleState(r) === 'unknown' ? null : ruleState(r)) : null; };
  const trig = id => { const s = g(id); return s === 'warn' || s === 'fail' ? s : null; };
  const out = (st.rules || []).find(r => /^yield_(reads|estimated_bases)_at_\d+h$/.test(r.id) && r.role === 'key') || (st.rules || []).find(r => /^yield_.*_at_\d+h$/.test(r.id));
  const h = out ? +/_at_(\d+)h$/.exec(out.id)[1] : null;
  const outTile = id => {
    const r = ruleById(b, id), isReads = /reads/.test(id), v = y[id], parts = isReads ? F.countParts(v) : F.basesParts(v), t = r && r.threshold, bl = t && t.baseline;
    const fmt = x => isReads ? F.count(x) : F.bases(x), hasBand = bl && bl.n >= 5 && bl.band, band = hasBand && c.typ ? bl.band : null, level = bl && bl.warn_level != null ? bl.warn_level : t && t.warn;
    const cap = hasBand ? 'typical ' + (isReads ? F.countParts(bl.band[0])[0] : F.basesParts(bl.band[0])[0]) + '–' + fmt(bl.band[1]) + ' · more is better' : t && t.warn != null ? 'fixed level ' + fmt(t.warn) + ' · no lab history yet' : '';
    return <Tile label={label(id)} glyph={r ? g(id) : null} state={trig(id)} value={parts[0]} unit={parts[1]} href={jump(id)} title={id}
      strip={r && (band || level != null) ? stripFor(v, band, [level]) : null} caption={cap} />;
  };
  const poresTile = () => { const r = ruleById(b, 'fc_pores_first_scan'), ex = (st.explanations || []).find(e => e.id === 'fc_pores_first_scan'), band = ex && ex.baseline && ex.baseline.band && c.typ ? ex.baseline.band : null, lim = r && r.threshold && r.threshold.warn;
    return <Tile label="Pores at the first scan" glyph={g('fc_pores_first_scan')} state={trig('fc_pores_first_scan')} value={F.int(fc.fc_pores_first_scan)} href={jump('fc_pores_first_scan')} strip={stripFor(fc.fc_pores_first_scan, band, [lim])}
      caption={(band ? (ex.direction === 'high' ? 'above typical · ' : ex.direction === 'low' ? 'below typical · ' : '') + 'typical ' + F.int(band[0]) + '–' + F.int(band[1]) + ' · ' : '') + (lim ? 'limit ' + F.int(lim) : '')} />; };
  const occTile = () => { const r = ruleById(b, 'fc_occupancy_first_hour'), short = D.runHours(b) < 1.5 || fc.fc_occupancy_first_hour == null;
    return <Tile label="Occupancy, first hour" glyph={short ? null : g('fc_occupancy_first_hour')} state={short ? 'unknown' : trig('fc_occupancy_first_hour')} value={short ? 'too short' : F.pctNum(fc.fc_occupancy_first_hour)} unit={short ? '' : '%'} href={jump('fc_occupancy_first_hour')}
      strip={short ? null : { lo: 0, hi: 1, mark: fc.fc_occupancy_first_hour, limits: r && r.threshold && F.isNum(r.threshold.warn) ? [r.threshold.warn] : [] }} caption={short ? 'measured from minute 30 to minute 90' : 'guide 70 %'} />; };
  if (fl.stub) return null;
  if (!fl.t2) {
    const reads = <Tile label="Reads" value={F.countParts(y.yield_reads)[0]} unit={F.countParts(y.yield_reads)[1]} caption={'in ' + F.hours(D.runHours(b))} />;
    const basesT = <Tile label="Estimated bases · from signal" value={F.basesParts(y.yield_estimated_bases)[0]} unit={F.basesParts(y.yield_estimated_bases)[1]} caption={h && y['yield_estimated_bases_at_' + h + 'h'] ? F.bases(y['yield_estimated_bases_at_' + h + 'h']) + ' at ' + h + ' h' : 'in ' + F.hours(D.runHours(b))} />;
    const atH = out ? outTile(out.id) : <Tile label="Output at the judged horizon" state="unknown" value="too short to judge" caption={'the run ' + (fl.running ? 'has run for ' : 'ended after ') + F.hours(D.runHours(b))} />;
    if (c.rna) tiles.push(reads, atH, basesT); else tiles.push(basesT, atH, reads);
    tiles.push(poresTile(), occTile());
    const n50r = ruleById(b, 'len_est_n50');
    const approx = ((D.mod(b, 'len.est') || {}).params || {}).quantile_source === 'histogram';
    tiles.push(!fl.full && le.len_est_n50 == null ? <Tile label="Estimated read N50" state="unknown" value={fl.cut ? 'not known' : 'being counted'} caption={fl.cut ? 'the read table of this run could not be read' : 'appears when the read count finishes'} />
      : <Tile label="Estimated read N50" glyph={g('len_est_n50')} state={trig('len_est_n50')} value={(approx ? '≈ ' : '') + F.basesParts(le.len_est_n50)[0]} unit={F.basesParts(le.len_est_n50)[1]} href={jump('len_est_n50')}
        title={approx ? 'Read from a histogram with 40 steps per tenfold length: good to a few per cent. Counting the reads of this run again (flowcellar scan --force --tier 1 --only ' + b.identity.run_key + ') makes it exact.' : 'Exact: counted read by read in the sequencing summary'}
        caption={(n50r && n50r.threshold && n50r.threshold.warn ? 'estimated from signal · limit ' + F.bases(n50r.threshold.warn) : 'estimated from signal') + (approx ? ' · approximate' : '')} />);
  } else {
    // Read data: FlowCellar's own numbers from every read of the linked files, or "not available yet"
    const fn = D.summary(b, 'funnel'), fu = D.funnel(b) || [], last = fu.filter(s => !s.missing).pop(), rl = b.tiers.read_level, pop = (D.mod(b, 'len') || {}).population || rl.population;
    const all = c.samples.map(s => s.key).concat(c.hasUnclassified ? ['unclassified'] : []), H = D.readHists(b, all, ['pass']), sy = D.sampleYield(b);
    const passBases = Object.values(sy).reduce((t, o) => t + ((o.pass || {}).bases || 0), 0), set = readSetWords(pop, ['pass']);
    const pl = partialInfo(b, 'len'), py = partialInfo(b, 'yield.basecalled') || pl, withP = (text, p) => p ? text + ' · ' + p.word : text;
    const n50 = H.len ? D.n50(H.len.edges, H.len.w) : null, qmed = H.qual ? D.quantile(H.qual.edges, H.qual.n, 0.5) : null;
    const usable = fn.funnel_usable_frac != null && last && last.id !== 'sequenced', na = 'not available yet';
    tiles.push(usable ? <Tile label="Usable reads" glyph={g('funnel_usable_frac')} state={trig('funnel_usable_frac')} value={F.pctNum(fn.funnel_usable_frac)} unit="%" href={jump('funnel_usable_frac')} caption={F.count(last.n) + ' of ' + F.count(fn.funnel_reads_sequenced) + ' sequenced'} sub={lower(last.label)} />
      : <Tile label="Usable reads" state="unknown" value={partialInfo(b, 'funnel') ? na : 'not known'} caption={partialInfo(b, 'funnel') ? 'needs the read data of this run to be complete' : 'no count of the reads that reached a sample is linked'} />);
    // the name of the number says which reads it counts: "Passed" only where the reads passed a quality filter
    const bl = basesLabel(pop && pop.stage), soon = 'appears when the read files of this run have been read';
    tiles.push(passBases > 0 ? <Tile label={bl} value={F.basesParts(passBases)[0]} unit={F.basesParts(passBases)[1]} caption={withP(set, py)} />
      : <Tile label={bl} state="unknown" value={na} caption={soon} />);
    tiles.push(n50 ? <Tile label="Read N50" value={F.basesParts(n50)[0]} unit={F.basesParts(n50)[1]} caption={withP(set, pl)} /> : <Tile label="Read N50" state="unknown" value={na} caption={soon} />);
    tiles.push(qmed ? <Tile label="Median read quality" value={F.q(qmed)} caption={withP(set, partialInfo(b, 'qual'))} /> : <Tile label="Median read quality" state="unknown" value={na} caption={soon} />);
    tiles.push(<Tile label="Samples with read data" value={rl.samples_with_data != null ? String(rl.samples_with_data) : F.DASH} unit={rl.samples_expected != null ? 'of ' + rl.samples_expected : ''}
      state={rl.samples_expected && rl.samples_with_data < rl.samples_expected ? 'unknown' : null} caption={rl.samples_expected && rl.samples_with_data < rl.samples_expected ? 'the others have no files yet' : 'every sample of the sample sheet'} />);
    tiles.push(<Tile label="Flow cell" glyph={g('fc_pores_first_scan')} state={trig('fc_pores_first_scan')} value={F.int(fc.fc_pores_first_scan)} unit="pores" href={jump('fc_pores_first_scan')} caption={F.pct(fc.fc_occupancy_first_hour) + ' occupancy, first hour'} />);
  }
  return <div class="fc-keynumbers">
    <div class="fc-tiles">{tiles}</div>
    <p class="fc-keynumbers__note">{c.peers >= 5 ? <>Lab typical <span class="fc-sw fc-sw--band" aria-hidden="true" />: compared with {c.peers} runs of this device, flow cell type and kit that started before this one.{h ? ' The judged horizon is ' + h + ' h on every part of this page.' : ''}</> : <>No lab history for this device, flow cell type and kit yet ({c.peers ? c.peers + ' comparable runs; comparisons start at 5' : 'no comparable run before this one'}): built-in limits only.{h ? ' The judged horizon is ' + h + ' h.' : ''}</>}</p>
  </div>;
}

// ---------- controls bar ----------
function Controls({ b, c, stuck }) {
  const fl = D.flags(b), noRead = !fl.t2, pop = (D.mod(b, 'len') || {}).population || b.tiers.read_level.population, q = pop && pop.min_qscore;
  const whyNone = 'No read data yet', whyOnly = 'These files hold only reads that passed' + (q ? ' Q ' + q : ' the quality filter');
  const presets = [['', 'Whole run · 0–' + F.hours(c.full[1])], ['0-1', 'First hour'], ['0-24', 'First 24 h'], ['24-48', '24 to 48 h'], [Math.max(0, c.full[1] - 24) + '-' + c.full[1], 'Last 24 h']].filter(p => { if (!p[0]) return true; const [a, e] = p[0].split('-').map(Number); return e <= c.full[1] + 1e-6 && a < e; });
  const tv = param('t', ''), known = presets.some(p => p[0] === tv);
  const off = c.timeOn || c.by !== c.byDefault || c.logY || (c.typOk && !c.typ) || c.sKeys.length || c.reads !== 'pass';
  const setT = (a, e) => { if (!(e > a)) return; setParams({ t: +a.toFixed(2) + '-' + +e.toFixed(2) }); };
  const offAll = noRead ? whyNone : c.bothFilters ? null : whyOnly;
  return <div class={'fc-controls fc-runcontrols' + (stuck ? ' is-stuck' : '')} role="group" aria-label="Controls for every chart">
    <span class="fc-runcontrols__id"><Glyph v={b.status.verdict} /> <span class="fc-mono">{b.identity.flow_cell.id}</span></span>
    <Seg name="c-reads" label="Reads" value={noRead ? '' : c.reads} onChange={v => setParams({ reads: v === 'pass' ? '' : v, minq: v === 'custom' ? (c.minq || 10) : '', minlen: v === 'custom' ? (c.minlen || '') : '' })}
      opts={[{ v: 'all', label: 'All', off: offAll }, { v: 'pass', label: 'Pass', off: noRead ? whyNone : null }, { v: 'fail', label: 'Fail', off: offAll }, { v: 'custom', label: 'Custom', off: noRead ? whyNone : null }]} />
    {c.custom ? <span class="fc-ctrl fc-custom" role="group" aria-label="Custom read filter" title="Cuts the charts of Read length and quality; the other sections say that they are not cut by it"><label>min Q <input type="number" min="0" max="40" step="1" value={c.minq} aria-label="Custom filter: lowest mean read quality, Q" onInput={e => setParams({ minq: e.currentTarget.value })} /></label><label>min length <input type="number" min="0" step="100" value={c.minlen} aria-label="Custom filter: shortest read length, bases" onInput={e => setParams({ minlen: e.currentTarget.value })} /> b</label></span> : null}
    <label class="fc-ctrl"><span class="fc-caps">Samples</span>
      {noRead ? <span class="fc-select is-off" tabIndex={0} title="No read data yet: MinKNOW recorded no sample names">not known yet</span>
        : <select class="fc-selectbox" value={c.sKeys.length === 1 ? c.sKeys[0] : c.sKeys.length ? '__multi' : ''} onChange={e => setParams({ s: e.currentTarget.value === '__multi' ? c.sKeys.join(',') : e.currentTarget.value })}>
          <option value="">{c.samples.length ? 'All ' + c.samples.length + ' added up' : 'All reads'}</option>{c.sKeys.length > 1 ? <option value="__multi">{c.sKeys.length} samples chosen</option> : null}{c.samples.map(s => <option value={s.key} title={s.name}>{c.short(s.key, 34)}</option>)}
          {c.hasUnclassified ? <option value="unclassified">No barcode found</option> : null}</select>}</label>
    <Seg name="c-by" label="Count by" value={c.by} onChange={v => setParams({ by: v === c.byDefault ? '' : v })} opts={[{ v: 'reads', label: 'Reads' }, { v: 'bases', label: noRead ? 'Bases (est.)' : 'Bases' }]} />
    <span class="fc-ctrl"><label class="fc-ctrl"><span class="fc-caps">Time</span>
      <select class="fc-selectbox" value={known ? tv : '__range'} onChange={e => { const v = e.currentTarget.value; if (v === '__set') setParams({ tf: '1' }); else setParams({ t: v === '__range' ? tv : v, tf: v ? '1' : '' }); }}>{presets.map(p => <option value={p[0]}>{p[1]}</option>)}{!known ? <option value="__range">Hours {+c.win[0].toFixed(1)} to {+c.win[1].toFixed(1)}</option> : null}<option value="__set">From and to, in hours {'…'}</option></select></label>
      {c.timeOn || param('tf', '') ? <><label class="fc-timefield">from <input type="number" min="0" max={c.full[1]} step="1" value={+c.win[0].toFixed(2)} aria-label="Time window from, hours" onChange={e => setT(+e.currentTarget.value, c.win[1])} /></label>
      <label class="fc-timefield">to <input type="number" min="0" max={c.full[1]} step="1" value={+c.win[1].toFixed(2)} aria-label="Time window to, hours" onChange={e => setT(c.win[0], Math.min(c.full[1], +e.currentTarget.value))} /> h</label></> : null}</span>
    <Seg name="c-y" label="Scale" value={c.logY ? 'log' : 'lin'} onChange={v => setParams({ y: v === 'lin' ? '' : v })} opts={[{ v: 'lin', label: 'Linear' }, { v: 'log', label: 'Log' }]} />
    <Seg name="c-typ" label="Lab typical" value={c.typ ? '1' : '0'} onChange={v => setParams({ typ: v === '1' ? '' : '0' })} opts={[{ v: '1', label: 'Shown', off: c.typOk ? null : 'Fewer than 5 comparable runs before this one: comparisons start at 5' }, { v: '0', label: 'Hidden' }]} />
    {off ? <button type="button" class="fc-btn fc-btn--ghost" onClick={() => setParams({ t: '', tf: '', by: '', y: '', typ: '', s: '', reads: '', minq: '', minlen: '' })}>Reset</button> : null}
  </div>;
}

// ---------- rail: a status column ----------
function Rail({ items, current }) {
  return <nav class="fc-rail fc-runrail" aria-label="Sections of this run">
    {items.map(it => it.group ? <div class="fc-rail__group fc-caps">{it.group}</div>
      : it.off ? <div class="fc-rail__item is-off"><span>{it.label}</span><span class="fc-rail__mark">{'–'}</span></div>
      : <a class={'fc-rail__item' + (current === it.id ? ' is-on' : '')} href={makeHash(route.value.path, { ...route.value.q, at: it.id === 'summary' ? '' : it.id, ev: '' })} aria-current={current === it.id ? 'location' : undefined}
          onClick={e => { e.preventDefault(); setParams({ at: it.id === 'summary' ? '' : it.id, ev: '' }); (it.id === 'summary' ? window.scrollTo(0, 0) : document.getElementById('s-' + it.id)?.scrollIntoView({ behavior: 'instant', block: 'start' })); }}>
        <span>{it.label}</span><span class="fc-rail__mark">{it.result ? <Glyph v={it.result} /> : it.mark || ''}</span></a>)}
  </nav>;
}

function Header({ b }) {
  const id = b.identity, pr = id.protocol_run, fl = D.flags(b), run = D.summary(b, 'run'), use = (b.relations || {}).flow_cell_use || { index: 1, of: 1 }, inSel = selected.value.includes(id.run_key);
  const n = D.samplesOf(b).length, bc = (id.basecalling || []).find(x => x.model);
  const ended = fl.running ? 'sequencing now' : END_WORD[run.run_stop_cause] || (run.run_stop_cause ? String(run.run_stop_cause).replace(/_/g, ' ').toLowerCase() : null);
  const planned = pr.planned_duration_s ? pr.planned_duration_s / 3600 : null, dev = id.device || {};
  // a folder with raw data only says nothing about the use of its flow cell: no chip claims a first use
  const useText = fl.stub ? null : use.of > 1 ? (use.index > 1 ? 're-used · use ' + use.index + ' of ' + use.of : 'first use · re-used later') : use.source === 'unknown' && !index.value ? 'use of the flow cell not known' : 'first use of this flow cell';
  return <header class="fc-runhead">
    <div class="fc-runhead__row">
      <h1 class="fc-runhead__id fc-mono">{id.flow_cell.id}</h1>
      <span class="fc-runhead__name" title="The MinKNOW sample ID: the name of the whole flow cell load">{pr.sample_id ? <Mid text={pr.sample_id} tail={12} /> : <span class="fc-muted">no run name</span>}</span>
      <span class="fc-runhead__chips">
        {id.assay.type !== 'unknown' ? <Chip href={hub('#/runs?q=assay:' + id.assay.group)} title={(id.assay.evidence || []).join('; ')}>{ASSAY_WORD[id.assay.type] || id.assay.type}{n ? ' · ' + F.plural(n, 'sample') : id.assay.barcoded ? ' · barcoding kit' : ''}</Chip> : null}
        {dev.type ? <Chip href={hub('#/runs?q=pos:' + dev.position)}>{DEVICE_WORD[dev.type] || dev.type}{/-[AB]$/.test(dev.position || '') ? ' · position ' + dev.position.slice(-1) : ''}</Chip> : dev.position ? <Chip>{dev.position}</Chip> : null}
        {id.flow_cell.product_code ? <Chip href={hub('#/runs?q=fctype:' + id.flow_cell.product_code)}>{id.flow_cell.product_code}</Chip> : null}
        {id.kit.sequencing_kit ? <Chip href={hub('#/runs?q=kit:' + id.kit.sequencing_kit)}>{id.kit.sequencing_kit}</Chip> : null}
        {bc ? <Chip title={bc.model}>{(bc.model_tier || 'model') + ' ' + ((/@(v[\d.]+)/.exec(bc.model) || [])[1] || '')}</Chip> : null}
        {useText ? <Chip href={hub(hrefFlowCell(id.flow_cell.id))} title={isReport ? undefined : 'Open the history of this flow cell'}>{useText}</Chip> : null}
      </span>
      <span class="fc-runhead__actions"><SavePdf charts={!fl.stub} />{isReport ? null : <button type="button" class="fc-btn" aria-pressed={inSel} onClick={() => toggleSelected(id.run_key)}>{inSel ? 'In compare ✓' : 'Add to compare'}</button>}</span>
    </div>
    <p class="fc-runhead__line">Started <When iso={pr.started} floating={fl.stub} time zone />{fl.stub ? ' (from the folder name, the sequencer’s clock)' : <> {'·'} {fl.running ? 'has run ' : 'ran '}{F.hours(run.run_duration_h)}{planned ? ' of ' + F.hours(planned) + ' planned' : ''}{ended ? ' · ' + ended : ''}</>}
      {' · '}{fl.stub ? 'no MinKNOW report in the folder' : <>run health computed <When iso={b.provenance.updated || b.provenance.created} time /> from {fl.full ? 'MinKNOW’s files and the sequencing summary' : 'MinKNOW’s report files'}</>}</p>
  </header>;
}

function RunCrumbs({ b, row, id }) {
  const fc = b ? b.identity.flow_cell.id : row ? row.fc : id, projects = b ? ((b.relations || {}).projects || []) : [];
  if (isReport) return <Crumbs>Run report {'›'} {projects.length ? <>{projects.map(p => p.name).join(', ')} {'›'} </> : null}<b class="fc-mono">{fc}</b></Crumbs>;
  return <Crumbs><a href="#/home">Library</a> {'›'} <a href="#/runs">Runs</a> {'›'} {projects.length ? <>{projects.map((p, k) => <>{k ? ', ' : ''}<a href={hrefProject(p.id)}>{p.name}</a></>)} {'›'} </> : null}<b class="fc-mono">{fc}</b></Crumbs>;
}

// The selection a run page shows, as the controls bar states it: on paper there is no bar to read it from.
function scopeItems(b, c) {
  const fl = D.flags(b), pop = (D.mod(b, 'len') || {}).population || b.tiers.read_level.population, items = [];
  if (fl.t2) {
    const names = c.sKeys.map(k => k === 'unclassified' ? 'No barcode found' : ((c.named.find(s => s.key === k) || {}).name || k));
    items.push(['Samples', c.sKeys.length ? names.slice(0, 6).join(', ') + (names.length > 6 ? ' and ' + (names.length - 6) + ' more' : '') + (names.length > 1 ? ' (' + names.length + ' of ' + c.samples.length + ')' : '') : c.samples.length ? 'all ' + c.samples.length + ' added up' : 'all reads']);
    items.push(['Reads', c.custom ? 'custom filter' + (c.minq ? ', Q ' + c.minq + ' or better' : '') + (c.minlen ? ', ' + F.int(c.minlen) + ' bases or longer' : '') + (c.minq || c.minlen ? '' : ', nothing cut')
      : c.reads === 'all' ? 'all, passed and failed' : c.reads === 'fail' ? 'failed the quality filter' : c.bothFilters ? 'passed the quality filter' : readSetWords(pop, ['pass'])]);
  }
  items.push(['Time', c.timeOn ? 'hours ' + +c.win[0].toFixed(2) + ' to ' + +c.win[1].toFixed(2) + ' of ' + F.hours(c.full[1]) : 'the whole run, 0–' + F.hours(c.full[1])]);
  items.push(['Count by', c.by === 'bases' ? (fl.t2 ? 'bases' : 'bases, estimated') : 'reads']);
  items.push(['Scale', c.logY ? 'logarithmic' : 'linear']);
  items.push(['Lab typical', c.typOk ? (c.typ ? 'shown' : 'hidden') : 'none yet: fewer than 5 comparable runs']);
  return items;
}

export function Run({ id }) {
  const row = rowById(id), b = useBundle(id), state = bundleState(id), err = bundleError(id);
  const [current, setCurrent] = useState('summary'), [stuck, setStuck] = useState(false), sentinel = useRef(null);
  const ev = param('ev', '');
  // follow a reason to its evidence: scroll to the panel, which is opened and outlined by the section
  useEffect(() => {
    if (!b || !ev) return; const r = ruleById(b, ev); if (!r) return;
    const t = setTimeout(() => { const el = (r.panel && document.getElementById('p-' + r.panel.replace('/', '-'))) || (secOf(b, r) === 'pychopper' && document.getElementById('p-pych-class')) || document.getElementById('s-' + secOf(b, r)); if (el) el.scrollIntoView({ behavior: 'instant', block: 'center' }); }, 60);
    return () => clearTimeout(t);
  }, [b, ev]);
  // a section named in the address (?at=flow_cell): a reload or a shared link opens there, not at the top
  const at = param('at', ''), atDone = useRef('');
  useEffect(() => { if (!b || !at || ev || atDone.current === id + '/' + at) return; atDone.current = id + '/' + at;
    const t = setTimeout(() => { const el = document.getElementById('s-' + at); if (el) el.scrollIntoView({ behavior: 'instant', block: 'start' }); }, 80); return () => clearTimeout(t); }, [b, at, id]);
  useEffect(() => { if (!ev && !at) window.scrollTo(0, 0); }, [id]);
  useEffect(() => {
    if (!b) return;
    const io = new IntersectionObserver(es => { for (const e of es) if (e.target === sentinel.current) setStuck(!e.isIntersecting); }, { rootMargin: '-41px 0px 0px 0px' });
    if (sentinel.current) io.observe(sentinel.current);
    const secs = [...document.querySelectorAll('[data-sec]')];
    const spy = () => { let cur = 'summary'; for (const s of secs) if (s.getBoundingClientRect().top < 180) cur = s.getAttribute('data-sec'); setCurrent(cur); };
    addEventListener('scroll', spy, { passive: true });
    return () => { io.disconnect(); removeEventListener('scroll', spy); };
  }, [b, id]);
  if (!row && !b && (state === 'failed' || state === 'absent')) return <div class="fc-page"><RunCrumbs id={id} /><div class="fc-empty fc-page__empty"><b>This run is not in the library.</b><span>No run has the key {id}. {isReport ? null : <a href="#/runs">All runs</a>}</span></div></div>;
  if (!b) return <div class="fc-page"><RunCrumbs row={row} id={id} />
    <div class="fc-runhead"><div class="fc-runhead__row"><h1 class="fc-runhead__id fc-mono">{row ? row.fc : id}</h1><span class="fc-runhead__name">{row && row.run_name}</span></div></div>
    {row ? <section class="fc-verdict"><div class="fc-verdict__head"><Badge lg v={row.verdict} soFar={row.provisional} /><div><p class="fc-verdict__headline">{row.headline}</p>{row.detail ? <p class="fc-verdict__detail">{row.detail}</p> : null}</div></div>
      {row.reasons.length ? <ul class="fc-reasons">{row.reasons.map((t, k) => <li><Glyph v={row.reason_levels[k] === 'fail' ? 'fail' : 'warn'} /><span>{t}</span></li>)}</ul> : null}</section> : null}
    {state === 'failed' ? <div class={'fc-banner fc-runfail' + (err && err.damaged ? ' fc-banner--strong' : '')} role="alert">
      <span>{err && err.damaged ? 'The data of this run is damaged: its checksum does not match, so nothing of it is drawn.' : 'The data file of this run could not be read (' + (err ? err.message : 'unknown reason') + '). The summary above comes from the library index.'}</span>
      {isReport ? null : <button type="button" class="fc-btn" onClick={() => retryBundle(id)}>Try again</button>}</div>
      : <p class="fc-loading">Loading the charts of this run …</p>}</div>;
  const c = controls(b), fl = D.flags(b), evRule = ev ? ruleById(b, ev) : null;
  const x = { b, c, row, evRule, evSec: evRule ? secOf(b, evRule) : null };
  const res = s => sectionResult(b, s), nIn = ((b.provenance || {}).inputs || []).length;
  const items = [{ id: 'summary', label: 'Summary', result: ['incomplete', 'running'].includes(b.status.verdict) ? null : b.status.verdict }];
  if (!fl.stub) {
    items.push({ group: 'Run health' }, { id: 'run_timeline', label: 'Run timeline', result: res('run_timeline') }, { id: 'flow_cell', label: 'Flow cell', result: res('flow_cell') });
    items.push(fl.full ? { id: 'reads_through_pore', label: 'Reads through the pore', result: res('reads_through_pore') } : { id: 'reads_through_pore', label: 'Reads through the pore', mark: fl.cut ? 'cut off' : fl.partlyRead ? 'in part' : 'counting' });
    if (fl.t2) {
      items.push({ group: 'Read data' + (fl.t2partial ? ' · partial' : '') });
      if (D.funnel(b)) items.push({ id: 'read_fate', label: 'Read fate', result: res('read_fate') });
      if (c.named.length) items.push({ id: 'samples', label: 'Samples', mark: String(c.named.length) });
      if (D.mod(b, 'len')) items.push({ id: 'length_quality', label: 'Read length and quality', result: res('length_quality') });
      if (D.mod(b, 'aln')) items.push({ id: 'alignment', label: 'Alignment', result: res('alignment') });
      if (D.mod(b, 'dyn.basecalled') || D.mod(b, 'aln.time')) items.push({ id: 'reads_over_time', label: 'Reads over run time' });
      if (D.mod(b, 'rna.tx') || D.mod(b, 'rna.polya') || D.mod(b, 'rna.coverage')) items.push({ id: 'transcripts', label: 'Transcripts' });
      { // library preparation tables: every category of the pychopper statistics (with the check on the full-length
        // share), and the annotated read tables. A bundle without the whole pychopper table shows its classes as "cDNA library".
        const py = D.mod(b, PYCHOPPER), rna = D.mod(b, BIOTYPE), cd = D.mod(b, 'cdna') && !py;
        if (py || rna || cd) items.push({ group: 'Library' });
        if (cd) items.push({ id: 'cdna', label: 'cDNA library', result: res('cdna') });
        if (py) items.push({ id: 'pychopper', label: 'Pychopper', result: res('cdna'), mark: py.state === 'partial' ? 'partial' : undefined });
        if (rna) items.push({ id: 'rna_composition', label: 'RNA composition', mark: rna.state === 'partial' ? 'partial' : undefined }); }
    } else items.push({ group: 'Read data · not available yet' }, { id: 'read_data', label: 'What will appear' });
    { const an = isReport ? [] : analysesOfRun(id); if (an.length) items.push({ group: 'Analyses' }, { id: 'analyses', label: 'Single-cell and spatial', mark: String(an.length) }); }
    items.push({ group: 'Details' }, { id: 'run_log', label: 'Run log', result: res('run_log') }, { id: 'setup', label: 'Setup and software' }, { id: 'data_sources', label: 'Data sources', mark: String(nIn) });
  }
  const projects = ((b.relations || {}).projects || []).map(p => p.name).filter(Boolean);
  return <div class="fc-page fc-run">
    <PaperHead kind="Run report" bundles={[b]} about={projects.length ? (projects.length > 1 ? 'projects ' : 'project ') + projects.join(', ') : ''} file={'FlowCellar run ' + b.identity.flow_cell.id + ' ' + String(b.identity.protocol_run.started || '').slice(0, 10)} />
    <RunCrumbs b={b} row={row} id={id} />
    <Header b={b} />
    <div class="fc-runtop"><Verdict b={b} /><KeyNumbers b={b} c={c} /></div>
    {fl.stub ? null : <PaperScope items={scopeItems(b, c)} groups />}
    {fl.stub ? <div class="fc-empty fc-page__empty"><b>Only raw data was found for this run.</b><span>{((b.provenance.warnings || [])[0] || {}).message} Run health appears when the folder holds MinKNOW's report files; read data when basecalled reads are linked. Folder: <span class="fc-mono fc-break">{(b.provenance.run_dir || {}).path}</span></span></div> : <>
      <div ref={sentinel} class="fc-sentinel" />
      <Controls b={b} c={c} stuck={stuck} />
      <div class="fc-runbody">
        <Rail items={items} current={current} />
        <main class="fc-runmain">
          <RunTimeline {...x} />
          <FlowCellSection {...x} />
          <ReadsThroughPore {...x} />
          {fl.t2 ? <><ReadFate {...x} /><SamplesSection {...x} /><LengthQuality {...x} /><Alignment {...x} /><ReadsOverTime {...x} /><Transcripts {...x} /><CdnaLibrary {...x} /><PychopperSection {...x} /><RnaComposition {...x} /></> : <ReadDataCard {...x} />}
          <RunAnalyses {...x} />
          <Details {...x} />
        </main>
      </div></>}
  </div>;
}
