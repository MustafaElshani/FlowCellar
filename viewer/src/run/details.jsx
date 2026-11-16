// Details: can I trust this page? Run log, setup and software, data sources. Closed by default, with a
// title line that summarises the content.
import * as F from '../core/format.js';
import * as D from '../core/derive.js';
import { Disclosure } from './section.jsx';
import { printing } from '../core/print.js';
import { DEVICE_WORD, ASSAY_WORD } from '../core/labels.js';

const STATE_WORD = { ok: 'read in full', partially_read: 'read in part', missing: 'missing', unreadable: 'could not be read', skipped: 'not read' };
const STATE_MORE = { growing: 'still being written', truncated: 'ends early: read up to the break', error: 'could not be read', represented: 'the same data as another file: counted once' };
const ROLE_FEEDS = { 'minknow.final_summary': 'Identity', 'minknow.sample_sheet': 'Identity, samples', 'minknow.report_json': 'Run timeline, Flow cell, Run log', 'minknow.throughput': 'Output', 'minknow.pore_activity': 'Channel activity, Occupancy',
  'minknow.sequencing_summary': 'Reads through the pore, Flow cell map', 'minknow.pore_scan_data': 'Flow cell map, open-pore current, voltage', 'minknow.barcode_alignment': 'Samples', 'lab.sample_sheet': 'Samples', 'dorado.sample_sheet': 'Samples',
  'pychopper.stats': 'cDNA library, Read fate', 'align.genome_bam': 'Alignment, Read length and quality', 'align.transcriptome_bam': 'Transcripts', 'polya.read_table': 'RNA composition', 'align.genome_bai': 'Read fate', 'quant.salmon_meta': 'Read fate', 'salmon.meta_info': 'Read fate',
  'quant.oarfish_meta': 'Read fate', 'quant.kallisto_run_info': 'Read fate', 'dorado.bam': 'Read length and quality, Samples' };
// A file feeds sections of the page: always in the words of the page, never as an internal name.
const SECTION_WORD = { run: 'Run timeline, Run log', yield: 'Output', fc: 'Flow cell', dyn: 'Reads through the pore', 'len.est': 'Reads through the pore', 'yield.basecalled': 'Read length and quality', len: 'Read length and quality',
  qual: 'Read length and quality', aln: 'Alignment', cdna: 'cDNA library', funnel: 'Read fate', 'funnel.samples': 'Read fate, Samples',
  'aln.windows': 'Alignment', 'aln.time': 'Reads over run time', 'dyn.basecalled': 'Reads over run time', 'rna.polya': 'Transcripts', 'rna.tx': 'Transcripts', 'rna.coverage': 'Transcripts',
  'cdna.pychopper': 'Pychopper', 'rna.biotype': 'RNA composition', 'rna.genes': 'RNA composition' };
const ROLE_WORD = { 'minknow.final_summary': 'MinKNOW final summary', 'minknow.sample_sheet': 'MinKNOW sample sheet', 'minknow.report_json': 'MinKNOW report', 'minknow.throughput': 'MinKNOW throughput table', 'minknow.pore_activity': 'MinKNOW pore activity table',
  'minknow.sequencing_summary': 'sequencing summary', 'minknow.pore_scan': 'pore-scan table', 'minknow.pore_scan_data': 'pore-scan table', 'minknow.barcode_alignment': 'MinKNOW barcode table', 'minknow.output_hash': 'MinKNOW list of output files', 'minknow.report_html': 'MinKNOW report page', 'minknow.report_md': 'MinKNOW report text',
  'basecaller.sample_sheet': 'basecaller sample sheet', 'dorado.bam': 'basecalled reads (BAM)', 'dorado.summary': 'basecaller summary', 'align.genome_bam': 'genome alignment (BAM)', 'align.transcriptome_bam': 'transcriptome alignment (BAM)', 'pychopper.stats': 'pychopper statistics',
  'quant.salmon': 'salmon result', 'quant.oarfish': 'oarfish result', 'quant.kallisto': 'kallisto result', 'polya.read_table': 'annotated read table', 'modkit.summary': 'modkit summary' };
export const roleWord = role => ROLE_WORD[role] || String(role || '').replace(/^[a-z]+\./, '').replace(/_/g, ' ');
function feeds(i) {
  const words = []; for (const m of i.used_by || []) for (const w of (SECTION_WORD[m] || '').split(', ')) if (w && !words.includes(w)) words.push(w);
  if (words.length) return words.join(', ');
  if (ROLE_FEEDS[i.role]) return ROLE_FEEDS[i.role];
  return <span class="fc-muted" title={i.not_counted_reason || (i.link && i.link.not_counted) || undefined}>{i.not_counted_reason ? 'not counted: two candidate sources' : 'not used on this page'}</span>;
}
// a note that names a file kind by its internal name says it in words
const roleInWords = text => String(text || '').replace(/\b(minknow|dorado|align|quant|pychopper|basecaller|polya|modkit)\.[a-z_]+\b/g, m => roleWord(m));
const size = F.bytes;
const LINK_WORD = { run_dir: 'inside the run folder', bam_header_rg: 'named in the file’s header', read_tag: 'named on its reads', summary_column: 'named in its rows', sample_sheet: 'named in the sample sheet', sample_name_chain: 'by sample name', path_rule: 'by a configured folder pattern' };
const coverageWords = i => { const cv = i.coverage; if (!cv || cv.kind === 'full') return ''; return ' · ' + (F.isNum(cv.fraction) ? F.pct(cv.fraction) + ' read' : 'read in part') + (cv.kind === 'truncated' ? ' (the file ends early)' : ''); };

export function Details({ b, c }) {
  const msgs = D.messages(b), run = D.cell0(b, 'run'), sc = D.scans(b), id = b.identity, pr = id.protocol_run, sw = id.software || {}, inputs = b.provenance.inputs || [];
  const nWarn = msgs.filter(m => m.severity === 'warning').length, nErr = msgs.filter(m => m.severity === 'error').length;
  const note = ((D.mod(b, 'run') || {}).notes || []).find(n => n.code === 'mux_scan_messages_moved'), rel = b.relations || {};
  const s = run ? run.scalars : {}, bc = (id.basecalling || [])[0] || {};
  const setup = [
    ['Flow cell', id.flow_cell.id + ' · ' + (id.flow_cell.product_code || 'type not known') + ' · ' + F.int(id.flow_cell.channel_count) + ' channels'],
    ['Device', (DEVICE_WORD[id.device.type] || id.device.type || F.DASH) + ' · position ' + (id.device.position || F.DASH)],
    ['Kit', (id.kit.sequencing_kit || F.DASH) + (id.kit.speed_bps ? ' · ' + id.kit.speed_bps + ' bases per second' : '')],
    ['Assay', (ASSAY_WORD[id.assay.type] || id.assay.type) + ' · found from ' + (id.assay.evidence || []).join(', ')],
    ['Protocol', pr.protocol || F.DASH],
    ['Run', (pr.id || pr.short_id) + ' · experiment ' + (pr.experiment_id || F.DASH) + ' · key ' + id.run_key],
    ['Acquisitions', (id.acquisitions || []).map(a => a.purpose + ' ' + a.short_id + (a.read_count ? ' (' + F.int(a.read_count) + ' reads)' : '')).join(' · ') || F.DASH],
    ['Planned run time', pr.planned_duration_s ? F.hours(pr.planned_duration_s / 3600) : F.DASH],
    ['MinKNOW', [sw.minknow && 'MinKNOW ' + sw.minknow, sw.bream && 'Bream ' + sw.bream, sw.distribution && 'distribution ' + sw.distribution].filter(Boolean).join(' · ') || F.DASH],
    ['Basecalling during the run', bc.live ? 'on' : 'off'],
    ['Basecalling afterwards', bc.model ? 'Dorado ' + (bc.version || '') + ' · ' + bc.model + (bc.min_qscore ? ' · kept at Q ' + bc.min_qscore + ' or better' : '') : 'not linked'],
    ['Raw data written', s.raw_bytes_completed ? size(s.raw_bytes_completed) + ' in ' + F.int(s.pod5_files || s.fast5_files) + ' files' + (s.files_in_fallback ? ' · ' + s.files_in_fallback + ' in a fall-back folder' : '') : F.DASH],
    ['Run folder', (b.provenance.run_dir || {}).path || F.DASH],
    ['Found under', (rel.paths || []).length ? rel.paths.map(p => p.root + ': ' + p.relpath).join(' · ') : F.DASH],
    ['Project', (rel.projects || []).length ? rel.projects.map(p => p.name).join(' · ') : 'none: the run folder sits directly under a scanned folder'],
    ['Made by', 'FlowCellar ' + ((b.provenance.generator || {}).version || '') + ' · data format ' + b.schema_version + ' · ' + F.dateTime(b.provenance.updated || b.provenance.created)],
  ];
  const forced = c.ev === 'run_messages' || c.ev === 'run_stop_cause' || c.ev === 'run_output_files';
  if (printing.value && !forced && !['log', 'setup', 'src'].some(k => c.more.has(k))) return null;      // on paper, a heading with three closed groups under it is a heading alone
  return <section class="fc-section fc-details" aria-labelledby="h-details">
    <div class="fc-section__head"><h2 class="fc-section__title" id="h-details">Details</h2></div>
    <div id="s-run_log" data-sec="run_log"><Disclosure id="log" c={c} force={forced} title={'Run log · ' + F.plural(nErr, 'error') + ' · ' + F.plural(nWarn, 'warning') + ' · ' + F.plural(sc ? sc.n : 0, 'pore scan') + ' (shown as one row) · ' + F.plural(msgs.length, 'other message')}>
      <div class="fc-tablewrap"><table class="fc-table fc-table--wrap"><caption class="fc-sr">Messages from MinKNOW</caption><thead><tr><th scope="col">Time<span class="fc-table__sub">{F.zone(pr.started)}</span></th><th scope="col" class="fc-num">Hour</th><th scope="col">Kind</th><th scope="col">Message</th></tr></thead>
        <tbody>{msgs.slice().sort((p, q) => (p.severity === 'error' ? 0 : p.severity === 'warning' ? 1 : 2) - (q.severity === 'error' ? 0 : q.severity === 'warning' ? 1 : 2) || p.h - q.h).map(m => <tr><td class="fc-nowrap" title={F.timeTitle(m.time)}>{F.dateTime(m.time)}</td><td class="fc-num">{Math.max(0, m.h).toFixed(1)}</td><td class="fc-nowrap">{m.severity === 'warning' ? 'Warning' : m.severity === 'error' ? 'Error' : 'Information'}</td><td>{m.text}</td></tr>)}
          {sc ? <tr><td class="fc-nowrap">every {((D.mod(b, 'fc') || {}).params || {}).mux_scan_period_h || 1.5} h</td><td class="fc-num">{F.DASH}</td><td>Information</td><td>{sc.n} pore scans, from hour {sc.t[0].toFixed(1)} to hour {sc.t[sc.n - 1].toFixed(1)}. {note ? note.message : ''}</td></tr> : null}</tbody></table></div>
    </Disclosure></div>
    <div id="s-setup" data-sec="setup"><Disclosure id="setup" c={c} title={'Setup and software · ' + (sw.minknow ? 'MinKNOW ' + sw.minknow : 'MinKNOW version not known') + (pr.planned_duration_s ? ' · run limit ' + F.hours(pr.planned_duration_s / 3600) : '') + ' · basecalling during the run ' + (bc.live ? 'on' : 'off')}>
      <div class="fc-tablewrap"><table class="fc-table fc-table--wrap fc-kv"><caption class="fc-sr">Setup and software</caption><tbody>{setup.map(([k, v]) => <tr><th scope="row">{k}</th><td class={/folder|Protocol|Run$|Found under/.test(k) ? 'fc-mono fc-break' : ''}>{v}</td></tr>)}</tbody></table></div>
    </Disclosure></div>
    <div id="s-data_sources" data-sec="data_sources"><Disclosure id="src" c={c} title={'Data sources · ' + F.plural(inputs.length, 'file') + ' · ' + (inputs.every(i => !i.link || i.link.confidence === 'certain' || i.link.confirmed !== false) ? 'all matched to this run' : 'some matches are not confirmed')}>
      <div class="fc-tablewrap"><table class="fc-table fc-table--wrap"><caption class="fc-sr">Files this page was made from</caption><thead><tr><th scope="col">File</th><th scope="col">Feeds</th><th scope="col" class="fc-num">Size</th><th scope="col">Match to this run</th><th scope="col">State</th></tr></thead>
        <tbody>{inputs.map(i => <tr><td class="fc-mono fc-break" title={i.path}>{i.relpath || i.path}</td><td>{feeds(i)}<span class="fc-table__sub">{roleWord(i.role)}</span></td><td class="fc-num fc-nowrap">{size(i.size)}</td>
          <td title={i.link && i.link.evidence ? i.link.evidence : undefined}>{i.link ? (LINK_WORD[i.link.method] || (i.link.method || '').replace(/_/g, ' ')) + ' · ' + i.link.confidence + (i.scope === 'shared' ? ' · the file also holds other runs' : '') : F.DASH}</td>
          <td>{STATE_WORD[i.state] || STATE_MORE[i.state] || i.state}{i.records != null ? ' · ' + F.int(i.records) + ' records' : ''}{coverageWords(i)}</td></tr>)}</tbody></table></div>
      {(b.provenance.warnings || []).map(w => <p class="fc-section__note fc-notice">{roleInWords(w.message)}</p>)}
    </Disclosure></div>
  </section>;
}
