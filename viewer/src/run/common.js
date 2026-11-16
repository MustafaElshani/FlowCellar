// Shared by the sections of the run view: the controls (all of them in the address), the scope line and corner
// tag of a chart, the status of a section, and the words for read data that is not complete.
import { index } from '../core/data.js';
import { param } from '../core/state.js';
import { printAll } from '../core/print.js';
import * as F from '../core/format.js';
import * as D from '../core/derive.js';
import { sectionOf, POPULATION_LONG } from '../core/labels.js';
import { ruleState } from './section.jsx';

// The groups of further charts of a run page (the ids of their drawers). "PDF with every chart" opens these for the
// length of a print: they are then open for every part of the page that asks, also for a section that builds the
// charts of its group only while the group is open.
export const CHART_GROUPS = ['tracks', 'fc', 'rtp', 'aln-more'];

// ---------- controls: all in the address ----------
export function controls(b) {
  const dur = Math.max(D.runHours(b), 1 / 6), full = [0, dur > 2 ? Math.ceil(dur) : dur];
  const t = param('t', ''); let win = full, timeOn = false;
  const m = /^(\d+(?:\.\d+)?)-(\d+(?:\.\d+)?)$/.exec(t); if (m && +m[2] > +m[1]) { win = [Math.max(0, +m[1]), Math.min(full[1], +m[2])]; timeOn = win[1] > win[0]; if (!timeOn) win = full; }
  const idn = b.identity, rna = idn.assay.group !== 'gdna', by = param('by', rna ? 'reads' : 'bases') === 'bases' ? 'bases' : 'reads';
  const group = idn.device.type + ' | ' + idn.flow_cell.product_code + ' | ' + idn.kit.sequencing_kit, g = index.value ? index.value.groups.get(group) : null;
  const judged = (b.status.rules || []).find(r => r.threshold && r.threshold.baseline);
  const peers = judged ? judged.threshold.baseline.n : 0;
  const typOk = peers >= 5, typ = typOk && param('typ', '1') !== '0';
  // samples: every sample that any section of the page holds a number for (reads, an alignment, a pychopper table, an
  // annotated read table ...), so the control never offers fewer samples than a section below it shows;
  // "No barcode found" is offered when the files hold such reads
  const withData = new Set(Object.keys(D.funnelSamples(b)));
  for (const id of Object.keys(b.modules || {})) if ((b.modules[id].dims || []).includes('sample')) for (const k of D.sampleKeysWithData(b, id)) withData.add(k);
  const named = D.samplesOf(b), samples = named.filter(s => withData.has(s.key)), shorten = F.shortNamer(named.map(s => s.name));
  const hasUnclassified = D.readCells(b, 'len').some(c => c.sample === 'unclassified') || D.readCells(b, 'yield.basecalled').some(c => c.sample === 'unclassified');
  const pickable = samples.map(s => s.key).concat(hasUnclassified ? ['unclassified'] : []);
  const sKeys = param('s', '').split(',').filter(k => pickable.includes(k));
  const keys = sKeys.length ? sKeys : pickable;                     // what a chart adds up
  // reads: Pass by default; All and Fail exist only when the files hold reads that failed the filter
  const has = D.filtersOf(b), both = has.pass && has.fail; let reads = param('reads', 'pass');
  if (!['all', 'pass', 'fail', 'custom'].includes(reads)) reads = 'pass';
  if (!both && (reads === 'all' || reads === 'fail')) reads = 'pass';
  const only = has.pass || !has.fail ? 'pass' : 'fail';
  const existing = both ? ['pass', 'fail'] : [only];
  const filters = reads === 'all' || reads === 'custom' ? existing : both ? [reads] : [only];
  const minq = +param('minq', 0) || 0, minlen = +param('minlen', 0) || 0, custom = reads === 'custom';
  return { full, win, timeOn, by, byDefault: rna ? 'reads' : 'bases', logY: param('y', 'lin') === 'log', typ, typOk, peers, group: g, groupName: group,
    samples, named, sKeys, keys, hasUnclassified, reads, filters, bothFilters: both, custom, minq, minlen, rna,
    // a sample's name where room is short (a read-out, a legend, a row label): the words all samples share are left out
    short: (k, max) => k === 'unclassified' ? 'No barcode found' : shorten((named.find(s => s.key === k) || {}).name || k, max), shorten,
    out: param('out', 'cum'), ev: param('ev', ''), more: new Set(param('more', '').split(',').filter(Boolean).concat(printAll.value ? CHART_GROUPS : [])), scans: param('scans', '') === 'raw' ? 'raw' : 'even' };
}
export const sampleName = (c, k) => k === 'unclassified' ? 'No barcode found' : ((c.named.find(s => s.key === k) || {}).name || k);
// "sample A", "A + B", "3 samples": the chosen samples in a sentence
export const chosenWords = c => c.sKeys.length > 2 ? c.sKeys.length + ' samples' : c.sKeys.map(k => sampleName(c, k)).join(' + ');
// scope line and corner tag: the tag is computed as active controls minus what the chart follows
export function scopeTag(c, follows, extra) {
  const f = new Set(follows), scope = [], tag = [];
  if (c.custom && (c.minq || c.minlen)) (f.has('custom') ? scope : tag).push(f.has('custom') ? 'custom filter' + (c.minq ? ', Q ≥ ' + c.minq : '') + (c.minlen ? ', ≥ ' + F.int(c.minlen) + ' b' : '') : 'not cut by the Custom filter');
  else if (c.reads !== 'pass' && c.bothFilters) (f.has('reads') ? scope : tag).push(f.has('reads') ? (c.reads === 'all' || c.reads === 'custom' ? 'all reads, passed and failed' : 'reads that failed the quality filter') : 'not affected by Reads');
  // several samples chosen: a chart either adds them up or draws them one over the other, and says which
  if (c.sKeys.length) (f.has('samples') ? scope : tag).push(f.has('samples') ? c.sKeys.map(k => c.short ? c.short(k, 30) : sampleName(c, k)).join(' + ') + (c.sKeys.length > 1 ? (f.has('overlay') ? ', one over the other' : ', added up') : '') : 'all samples');
  if (c.by !== c.byDefault && f.has('by')) scope.push('by ' + c.by);
  if (c.timeOn) (f.has('time') ? scope : tag).push(f.has('time') ? 'hours ' + +c.win[0].toFixed(1) + ' to ' + +c.win[1].toFixed(1) : 'whole run only');
  if (c.logY && f.has('y')) scope.push('log scale');
  if (extra) tag.unshift(extra);
  return { scope: scope.join(' · '), tag: tag.join(' · ') };
}

export function sectionResult(b, sec) {
  let worst = null; const order = { pass: 0, unknown: 0, warn: 1, fail: 2 };
  for (const r of b.status.rules || []) { if (sectionOf(r) !== sec) continue; const s = ruleState(r); if (s === 'unknown') { if (worst == null) worst = 'unknown'; continue; } if (worst == null || worst === 'unknown' || order[s] > order[worst]) worst = s; }
  return worst;
}
export const ruleById = (b, id) => (b.status.rules || []).find(r => r.id === id);

// ---------- read data that is not complete ----------
// Complete (every read of every file) is the normal case and carries no mark. Anything else is said three times:
// a line under the section title, a tag on every chart, and a word in the captions of its tiles.
export function partialInfo(b, id) {
  const st = D.sectionState(b, id); if (!st || !st.partial) return null;
  const rl = (b.tiers || {}).read_level || {}, cov = st.coverage || rl.coverage || {}, frac = F.isNum(cov.fraction) ? cov.fraction : null;
  switch (st.reason) {
    case 'samples_missing': return { tag: 'PARTIAL · ' + rl.samples_with_data + ' of ' + rl.samples_expected + ' samples', word: 'partial: ' + rl.samples_with_data + ' of ' + rl.samples_expected + ' samples',
      text: 'Partial: read data exists for ' + rl.samples_with_data + ' of ' + rl.samples_expected + ' samples. The numbers and charts here cover those samples only; the others have no files yet.' };
    case 'partly_read': return { tag: 'PARTIAL · ' + (frac != null ? F.pct(frac) + ' read' : 'read in part'), word: 'partial: ' + (frac != null ? F.pct(frac) + ' read' : 'read in part'),
      text: 'Partial: ' + (F.isNum(cov.records) && F.isNum(cov.records_total) ? F.int(cov.records) + ' of ' + F.int(cov.records_total) + ' records' : 'part of the file') + (frac != null ? ' (' + F.pct(frac) + ')' : '') + ' had been read when this page was made. Counts will grow and shapes may still shift; the rest is read by the next scan.' };
    case 'waiting': return { tag: 'PARTIAL · files waiting', word: 'partial: files waiting', text: 'Partial: some files of this run are linked but have not been read yet. The numbers here will change when the next scan has read them.' };
    default: {
      // read fate: say which stage is not known and why, instead of "a file was not found"
      if (id === 'funnel' || id === 'funnel.samples') {
        const miss = (D.funnel(b) || []).filter(s => s.missing), why = s => (s.source || 'not known').replace(/^not (known|assessed)( yet)?( for [^:]+)?: /, '');
        if (miss.length) return { tag: 'PARTIAL · ' + (miss.length === 1 ? miss[0].label.toLowerCase() + ' not known' : miss.length + ' stages not known'), word: 'partial: ' + (miss.length === 1 ? miss[0].label.toLowerCase() + ' not known' : miss.length + ' stages not known'),
          text: (miss.length === 1 ? 'One stage of read fate is not known: ' : miss.length + ' stages of read fate are not known: ') + miss.map(s => s.label + ' (' + why(s) + ')').join('; ') + '. Such a stage is drawn as an outline, never as zero; the other stages are whole.' };
      }
      return { tag: 'PARTIAL · files missing', word: 'partial: files missing', text: 'Partial: a file this section needs was not found, so some of its numbers are missing.' };
    }
  }
}
// "99.8 % of the full-length, trimmed reads": the read population inside a sentence
const OF_READS = { basecalled_all: 'the basecalled reads', basecalled_pass: 'the reads that passed the quality filter', classified: 'the reads given to a sample', full_length_trimmed: 'the full-length, trimmed reads', aligned_primary: 'the aligned reads' };
export const ofReads = pop => (pop && OF_READS[pop.stage]) || 'the reads';
// The reads a number describes, in words: the read population of the files, narrowed by the Reads control.
export function readSetWords(pop, filters) {
  const q = pop && pop.min_qscore ? ', Q ≥ ' + pop.min_qscore : '', f = filters || ['pass'];
  if (!pop || !pop.stage) return f.length > 1 ? 'all reads' : f[0] === 'fail' ? 'reads that failed the quality filter' : 'reads';
  if (pop.stage === 'basecalled_all') return f.length > 1 ? 'every read the basecaller wrote' : f[0] === 'fail' ? 'reads that failed the quality filter' + (pop.min_qscore ? ', below Q ' + pop.min_qscore : '') : 'reads that passed the quality filter' + q;
  return (POPULATION_LONG[pop.stage] || pop.stage.replace(/_/g, ' ')) + q;
}

// ---------- notes of a section that the page must say, not only store ----------
// A table without identifiers on its reads (a pychopper statistics table, an annotated read table) that covers a
// sample on two flow cells is ONE analysis of reads that were merged before the tool ran: the page of each flow cell
// shows that one analysis, and says so. { keys: Set of sample keys, names, n: flow cells, others: [flow cell ids],
// text, caption, legend, both, mark(key) } or null when no table of the sections is shared. A read total of such a
// section is never shown without `text`. kind: 'pychopper' (a pychopper statistics table) or 'table' (an annotated
// read table).
const SHARED_WORDS = {
  pychopper: { one: 'pychopper analysis', it: 'analysis', when: 'before pychopper ran', short: 'merged before pychopper', legend: 'were merged before pychopper: one analysis for the sample' },
  table: { one: 'read table', it: 'table', when: 'before the table was made', short: 'merged before the read table was made', legend: 'were merged before the read table was made: one table for the sample' },
};
const listed = a => a.length < 2 ? a.join('') : a.slice(0, -1).join(', ') + ' and ' + a[a.length - 1];
export function sharedTables(b, ids, c, kind) {
  const inputs = new Map((((b.provenance || {}).inputs) || []).map(i => [i.id, i])), keys = new Set(), others = new Set(); let n = 0;
  for (const id of ids) for (const note of ((D.mod(b, id) || {}).notes || [])) {
    if (note.code !== 'covers_several_runs') continue;
    const inp = inputs.get(note.input); n = Math.max(n, +note.value || ((inp && (inp.shared_with || []).length) || 1) + 1);
    for (const k of (inp && inp.samples) || []) keys.add(k);
    for (const r of (inp && inp.shared_with) || []) others.add(String(r).split('_')[0]);
    if (!inp || !(inp.samples || []).length) { const name = String(note.message || '').split(':')[0], hit = c.named.find(q => q.name === name); if (hit) keys.add(hit.key); }
  }
  if (!keys.size) return null;
  const mine = c.named.filter(q => keys.has(q.key)), all = (D.mod(b, ids[0]) || {}).cells ? new Set(D.sampleKeysWithData(b, ids[0])) : keys, every = [...all].every(k => keys.has(k));
  const w = SHARED_WORDS[kind] || SHARED_WORDS.table, own = b.identity.flow_cell.id, other = [...others].filter(id => id !== own), both = n > 2 ? 'all ' + n + ' flow cells' : 'both flow cells';
  const cells = n > 2 ? n + ' flow cells' : 'two flow cells', single = every ? all.size === 1 : mine.length === 1;
  const who = every ? (all.size === 1 ? 'this sample' : 'every sample here') : mine.length + ' of the ' + all.size + ' samples here (' + mine.slice(0, 4).map(q => q.name).join(', ') + (mine.length > 4 ? ', …' : '') + ')';
  // the other flow cells by their ids where the files name them all; otherwise in words that stay true ("another flow cell")
  const named = other.length === n - 1, cellsMerged = named ? listed(['this flow cell'].concat(other)) : 'this flow cell and ' + (n > 2 ? (n - 1) + ' other flow cells' : 'another flow cell');
  const pages = named ? 'the page' + (other.length > 1 ? 's' : '') + ' of ' + listed(other) : n > 2 ? 'the pages of the other flow cells' : 'the page of the other flow cell';
  return { keys, names: mine.map(q => q.name), n, others: other, every, both,
    text: 'One ' + w.one + ' for ' + cells + ': the reads of ' + cellsMerged + ' were merged ' + w.when + ', so ' + who + ' ' + (every || single ? 'has' : 'have') + ' a single ' + w.one + ' of ' + both
      + '. Every count of ' + (every ? 'this section' : 'those samples') + ' is that merged ' + w.it + ', not this flow cell alone; the same ' + w.it + ' is shown on ' + pages + '.',
    caption: both + ', ' + w.short, mark: k => keys.has(k) ? ' *' : '', legend: '* reads of ' + both + ' ' + w.legend };
}
// "Read from the files of version X": a section whose source is not the primary version of its files
export const variantNote = (b, ids) => { for (const id of ids) { const n = ((D.mod(b, id) || {}).notes || []).find(q => q.code === 'source_variant'); if (n) return n.message; } return null; };
// the files behind the given sections that are linked to the run but not read yet: how many, their size, and the
// command that reads them (with the path of this library when the page knows it)
export function waitingFiles(b, ids) {
  const mine = new Set(); for (const id of ids) for (const k of ((D.mod(b, id) || {}).inputs || [])) mine.add(k);
  const of = (((b.provenance || {}).inputs) || []).filter(i => mine.has(i.id) && (i.notes || []).some(n => n.code === 'waiting'));
  const lib = index.value && index.value.library && index.value.library.path;
  return { n: of.length, bytes: of.reduce((a, i) => a + (i.size || 0), 0), command: 'flowcellar scan -L ' + (lib || '<library>') + ' --only ' + b.identity.run_key };
}
// a section whose files all wait and that holds no number yet: only this notice is shown (no sentence, tile or chart
// that would print dashes and zeros as if they were measured)
export function waitingOnly(b, ids, what) {
  const st = D.sectionState(b, ids[0]), m = D.mod(b, ids[0]); if (!m || !st || !st.partial || st.reason !== 'waiting' || (m.cells || []).length) return null;
  const w = waitingFiles(b, ids);
  return { n: w.n, command: w.command, text: (w.n ? F.plural(w.n, 'file') + ' of this run ' + (w.n === 1 ? 'is' : 'are') + ' linked and wait' + (w.n === 1 ? 's' : '') + ' for the next scan' + (w.bytes ? ' (' + F.bytes(w.bytes) + ')' : '') : 'The files of this run are linked and wait for the next scan')
    + ': ' + what + ' appear when they have been read.' };
}
