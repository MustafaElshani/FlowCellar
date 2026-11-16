// Interface words. Never shown: tier, rule, threshold, baseline, z-score, bundle, card, partial, sidecar.
import cat from '../gen/catalogue.js';

// The overall verdict has six values; each of its two parts may also be "not_available" (that part has no data).
export const VERDICT_WORD = { pass: 'Pass', warn: 'Warning', fail: 'Fail', incomplete: 'Incomplete', running: 'Running', unknown: 'Not assessed', not_available: 'Not available yet' };
export const ASSAY_WORD = { cdna_pcr: 'PCR-cDNA', cdna_direct: 'Direct cDNA', cdna: 'cDNA', gdna: 'Genomic DNA', direct_rna: 'Direct RNA', drna: 'Direct RNA', amplicon: 'Amplicon', unknown: 'Assay not known' };
export const ASSAY_SHORT = { cdna: 'cDNA', gdna: 'gDNA', drna: 'direct RNA', unknown: '' };
export const DEVICE_WORD = { P2_SOLO: 'PromethION 2 Solo', P2I: 'PromethION 2 Integrated', MINION_MK1D: 'MinION Mk1D', MINION_MK1B: 'MinION Mk1B', MINION_MK1C: 'MinION Mk1C', P2: 'PromethION 2', P24: 'PromethION 24', P48: 'PromethION 48', GRIDION: 'GridION' };
export const POPULATION_WORD = { basecalled_all: 'all basecalled reads', basecalled_pass: 'passed reads', classified: 'classified reads', full_length_trimmed: 'full-length reads', aligned_primary: 'aligned reads' };
export const POPULATION_LONG = { basecalled_all: 'every read the basecaller wrote', basecalled_pass: 'reads that passed the quality filter', classified: 'reads given to a sample', full_length_trimmed: 'full-length reads, trimmed', aligned_primary: 'aligned reads' };

export const metric = id => cat.metrics[id] || {};

// A readable name for a number the catalogue has no display copy for: the id in words. The id itself stays
// available as the title of whatever shows the label, so it can still be searched and quoted.
const DROP = new Set(['run', 'fc', 'yield', 'dyn', 'len', 'qual', 'aln', 'cdna', 'funnel', 'bc', 'rna', 'drna', 'mod', 'dup', 'as']);
const WORD = { est: 'estimated', frac: 'share', n50: 'N50', n90: 'N90', mapq: 'mapping quality', mux: 'pore scan', umi: 'UMI', gc: 'GC', polya: 'poly(A)', minknow: 'MinKNOW', dev: 'off target', max: 'largest', min: 'smallest', pct: 'per cent', num: 'number of' };
const UNIT = { h: 'hours', s: 'seconds', ms: 'milliseconds', c: '°C', mv: 'mV', mb: 'Mb', gb: 'Gb', kb: 'kb', bps: 'bases per second', pa: 'pA', frac: 'share' };
export function readableId(id) {
  const t = String(id).split('_').filter(Boolean); if (t.length > 1 && DROP.has(t[0])) t.shift();
  let unit = ''; const last = t[t.length - 1];
  if (t.length > 1 && UNIT[last] && !/^\d/.test(last)) { unit = UNIT[last]; t.pop(); }
  const words = t.map(w => /^\d+h$/.test(w) ? w.slice(0, -1) + ' h' : WORD[w] || w);
  let s = words.join(' ').replace(/\bat (\d+ h)\b/, 'at $1').replace(/\bper (\w+)/g, 'per $1');
  s = s.charAt(0).toUpperCase() + s.slice(1);
  return unit ? s + ', ' + unit : s;
}
// { label, known }: known is false when the label was derived from the id
export function labelOf(id) {
  const m = cat.metrics[id];
  if (m && m.label) return { label: m.label, known: true };
  return { label: readableId(id), known: false };
}
export const label = id => labelOf(id).label;
// the name of a check: its label, else the words before the colon of its sentence
export function checkName(rule) {
  const l = labelOf(rule.id); if (l.known) return l.label;
  const i = (rule.message || '').indexOf(':'); return i > 0 && i < 48 && !/^Not assessed/.test(rule.message) ? rule.message.slice(0, i) : l.label;
}
// where the evidence of a check lives: panel id -> section
export const SECTIONS = {
  summary: 'Summary', run_timeline: 'Run timeline', flow_cell: 'Flow cell', reads_through_pore: 'Reads through the pore',
  read_fate: 'Read fate', samples: 'Samples', length_quality: 'Read length and quality', alignment: 'Alignment', cdna: 'cDNA library',
  reads_over_time: 'Reads over run time', transcripts: 'Transcripts',
  run_log: 'Run log', setup: 'Setup and software', data_sources: 'Data sources',
};
const PANEL_SECTION = {}; for (const p of cat.panels || []) PANEL_SECTION[p.id] = p.section;
const MODULE_SECTION = { funnel: 'read_fate', 'funnel.samples': 'read_fate', cdna: 'cdna', aln: 'alignment', len: 'length_quality', qual: 'length_quality', 'yield.basecalled': 'length_quality',
  'aln.windows': 'alignment', 'aln.time': 'reads_over_time', 'dyn.basecalled': 'reads_over_time', 'rna.polya': 'transcripts', 'rna.tx': 'transcripts', 'rna.coverage': 'transcripts' };
export function sectionOf(rule) { const s = (rule.panel && PANEL_SECTION[rule.panel]) || MODULE_SECTION[rule.module]; return s && SECTIONS[s] ? s : 'summary'; }
export const panelTitle = id => { const p = (cat.panels || []).find(x => x.id === id); return p ? p.title : id; };
export const TRACK_CLOSED = new Set(['dyn/speed_over_time', 'dyn/why_reads_ended', 'fc/temperature', 'fc/voltage']);

export const MODIFIER_WORD = { reused_flow_cell: 're-used', flow_cell_use_unknown: 'use not known', baseline_provisional: 'early comparison', baseline_spread_guard: 'wide lab spread', quick_detail: 'counting reads', no_history: 'no lab history yet',
  judged_before_older_runs: 'judged before older runs', custom_profile: 'custom limits', reads_partly_counted: 'reads partly counted' };
export const CAUSE_WORD = { raw_only: 'raw data only', stopped_early: 'stopped early', output_not_judged: 'too short to judge', no_final_summary: 'no final summary', inputs_missing: 'files missing' };
export const END_WORD = { PROTOCOL_COMPLETED: 'ended as planned', PROTOCOL_STOPPED_BY_USER: 'stopped by a person', PROTOCOL_FINISHED_WITH_ERROR: 'ended with an error' };
export const ROLE_WORD = { key: 'fail the run', supporting: 'warn', info: 'shown only' };
export const ROOT_STATE_WORD = { reachable: 'Reachable', offline: 'Offline', skipped: 'Skipped in the last scan' };
// why read data is not complete, in words (section.partial_reason)
export const PARTIAL_WORD = { inputs_missing: 'files missing', partly_read: 'read in part', samples_missing: 'not every sample', waiting: 'files waiting' };
