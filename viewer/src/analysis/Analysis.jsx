// The analysis view (wireframe 08b): one single-cell or spatial sample. It has the skeleton of the run view: identity,
// the runs it came from, a verdict with its reasons, six key numbers, a rail, sections. The ONT workflow analysis and
// the Space Ranger output of the same reads are two tabs of one page. The verdict is about the analysis: it never
// changes the verdict of a flow cell run, and the page says so.
import { useEffect, useState } from 'preact/hooks';
import { useBundle, bundleState, bundleError, retryBundle, rowById, isReport, index } from '../core/data.js';
import { route, setParams, param, makeHash, hrefRun } from '../core/state.js';
import { Badge, Glyph, Chip, Tile, Seg, Crumbs, When, Mid } from '../shell/parts.jsx';
import { ChartPanel, DataTable, Empty } from '../charts/ChartPanel.jsx';
import { tk } from '../charts/theme.js';
import { Section, ruleState } from '../run/section.jsx';
import * as F from '../core/format.js';
import { VERDICT_WORD } from '../core/labels.js';
import { SpatialMap } from './SpatialMap.jsx';
import { Grid } from '../charts/Grid.jsx';
import * as S from './sc.js';
import { SavePdf, PaperHead, PaperScope } from '../shell/paper.jsx';

const hub = href => isReport ? undefined : href;
const FromReport = () => <span class="fc-tag fc-fromreport" title="This number exists only in the workflow's own HTML report, not in a data file: it was read from the report.">from report</span>;
const withTag = (b, key, text) => S.fromReport(b, key) ? <>{text} <FromReport /></> : text;
const centers = (e, log) => Array.from({ length: e.length - 1 }, (_, i) => log ? Math.sqrt(e[i] * e[i + 1]) : (e[i] + e[i + 1]) / 2);
const span = h => { let a = 0, b = h.n.length - 1; while (a < b && !(h.n[a] > 0)) a++; while (b > a && !(h.n[b] > 0)) b--; return [Math.max(0, a - 1), Math.min(h.n.length - 1, b + 1)]; };
const total = h => { let t = 0; for (const v of h.n) t += v; return t; };

// ---------- the runs an analysis came from ----------
// A link whose run is in the library is a link to that run with its verdict; one whose run is not is a greyed chip
// that says so. The index joins on the run key, so the chip becomes a link by itself once that run is scanned.
export function RunChips({ links, small }) {
  if (!links || !links.length) return <span class="fc-runchips"><span class="fc-runchip is-missing" title="Nothing in the files of this analysis names a flow cell run. State it with: flowcellar analysis DIR --run FLOWCELL">no flow cell run named</span></span>;
  return <span class="fc-runchips">{links.map(ln => {
    const row = ln.run ? rowById(ln.run) : null, how = S.LINK_WORD[ln.method] || '', conf = ln.confidence && ln.confidence !== 'high' ? ' (' + ln.confidence + ' confidence)' : '';
    // A share comes from the records that were probed, not from all of them: it is "about", to the whole per cent. A link
    // made by the sample name alone knows nothing about shares: it says how it was made, where it can be seen.
    const byName = ln.method === 'sample_name_chain' || ln.confidence === 'low', share = byName ? ' · by sample name only' : F.isNum(ln.share) ? ' ' + S.aboutShare(ln.share) : '';
    if (row && !isReport) return <a class="fc-runchip" href={hrefRun(row.id)} title={cap1(how) + conf + (F.isNum(ln.share) && !byName ? '. The share is of the records that were looked at, not of all records' : '') + '. Open the run.'}><span class="fc-mono fc-id">{row.fc}</span>{share}{small ? null : <Badge v={row.verdict} soFar={row.provisional} />}</a>;
    return <span class="fc-runchip is-missing" title={cap1(how) + conf + (isReport ? '.' : '. This run is not in the library; the chip becomes a link when it is scanned.')}><span class="fc-mono fc-id">{ln.fc || 'run'}</span>{share}{isReport ? '' : ' · run not in library'}</span>;
  })}</span>;
}
const cap1 = s => s ? s.charAt(0).toUpperCase() + s.slice(1) : '';

// ---------- verdict ----------
const SECTION_OF = { sp_pipeline_alarms: 'summary', sp_registration: 'map', sp_units_under_tissue: 'map', sp_fraction_reads_under_tissue: 'map', sp_segmentation_summary: 'bins',
  sr_usable_reads_fraction: 'fate', sr_valid_barcodes: 'fate', sr_valid_umis: 'fate', sr_mapping: 'mapping', sr_antisense_fraction: 'mapping', sc_full_length_fraction: 'fate', sc_strand_plus_fraction: 'fate',
  sc_valid_barcode_fraction: 'fate', sc_valid_umi_fraction: 'fate', sc_gene_assigned_fraction: 'fate', sc_transcript_assigned_fraction: 'fate', sc_reads_chrm_fraction: 'units', sc_mito_fraction_per_barcode: 'units',
  sc_cells: 'rank', sp_wf_sr_umi_ratio: 'pair', sc_wf_execution: 'workflow', sc_link_runs: 'summary' };
const jump = id => e => { e.preventDefault(); setParams({ at: id }); const el = document.getElementById('s-' + id); if (el) el.scrollIntoView({ behavior: 'instant', block: 'start' }); };
// An analysis made from the reads of several flow cells is one analysis of all of them, merged before the pipeline ran:
// the sentence under the verdict says so (the stored sentence is the same for every analysis and does not).
const ABOUT_ANALYSIS = 'This verdict is about the analysis. It does not change the verdict of any flow cell.';
const aboutAnalysis = (stored, runs) => runs > 1 && (!stored || stored === ABOUT_ANALYSIS) ? 'This verdict is about the analysis: one analysis of the reads of ' + (runs === 2 ? 'both' : 'all ' + runs) + ' flow cells together. It does not change the verdict of any flow cell.' : stored || ABOUT_ANALYSIS;
function Verdict({ b, runs }) {
  const st = b.status || {}, rules = st.rules || [], [all, setAll] = useState(false);
  const trig = rules.filter(r => r.outcome === 'triggered' && (r.severity === 'warn' || r.severity === 'fail') && r.role !== 'info').sort((a, x) => (a.id === 'sp_pipeline_alarms' ? 0 : 1) - (x.id === 'sp_pipeline_alarms' ? 0 : 1));
  const notes = rules.filter(r => r.outcome === 'triggered' && r.role === 'info'), ok = rules.filter(r => r.outcome === 'ok' && r.role !== 'info').length, skipped = rules.filter(r => r.outcome === 'not_evaluated' && r.role !== 'info').length;
  return <section class="fc-verdict fc-anverdict" aria-label="Verdict of this analysis">
    <div class="fc-verdict__head"><Badge lg v={st.verdict || 'unknown'} />
      <div><p class="fc-verdict__headline">{st.headline}</p><p class="fc-verdict__detail">{aboutAnalysis(st.detail, runs)}</p></div></div>
    {trig.length ? <ul class="fc-reasons">{trig.map(r => <li><Glyph v={ruleState(r)} /><span>{r.message}</span>
      {SECTION_OF[r.id] && SECTION_OF[r.id] !== 'summary' ? <a class="fc-reasons__go" href={makeHash(route.value.path, { ...route.value.q, at: SECTION_OF[r.id] })} onClick={jump(SECTION_OF[r.id])}>Evidence</a> : null}</li>)}</ul> : null}
    {(st.modifiers || []).map(m => <p class="fc-verdict__note">{m.message}</p>)}
    {notes.map(r => <p class="fc-verdict__note">{r.message}</p>)}
    <p class="fc-verdict__note"><button type="button" class="fc-linkbtn" aria-expanded={all} onClick={() => setAll(!all)}>{all ? 'Hide the checks' : 'All ' + rules.filter(r => r.role !== 'info').length + ' checks'}</button>
      {' '}{ok} passed{trig.length ? ' · ' + trig.length + ' triggered' : ''}{skipped ? ' · ' + skipped + ' not assessed (no source for the number)' : ''}</p>
    {all ? <div class="fc-tablewrap fc-allchecks__wrap"><table class="fc-table fc-table--wrap"><caption class="fc-sr">Every check of this analysis</caption>
      <thead><tr><th scope="col">Check</th><th scope="col">Result</th><th scope="col">What was found</th><th scope="col">Limit</th></tr></thead>
      <tbody>{rules.map(r => { const t = r.threshold || {}, lim = t.pass_band ? 'outside ' + t.pass_band.join(' to ') : [t.warn != null ? 'warn ' + (t.kind === 'max' ? 'above ' : 'below ') + fmtLimit(r, t.warn) : null, t.fail != null ? 'fail ' + (t.kind === 'max' ? 'above ' : 'below ') + fmtLimit(r, t.fail) : null].filter(Boolean).join(', ');
        return <tr><th scope="row">{S.metricLabel(r.id)}</th><td class="fc-nowrap"><Glyph v={ruleState(r)} /> {r.outcome === 'not_evaluated' ? 'Not assessed' : r.outcome === 'ok' ? 'Passed' : r.role === 'info' ? 'Noted' : VERDICT_WORD[ruleState(r)]}</td><td>{r.message}</td>
          <td>{lim || (t.kind === 'rule' ? 'a rule, not a number' : F.DASH)}{t.basis ? <span class="fc-table__sub">{t.basis}</span> : null}</td></tr>; })}</tbody></table></div> : null}
  </section>;
}
const fmtLimit = (r, v) => (r.metric || '').includes('units') || (r.metric || '').includes('cells') ? F.int(v) : r.id === 'sp_wf_sr_umi_ratio' ? String(v) : F.pct(v);

// ---------- six key numbers ----------
function KeyNumbers({ b, other }) {
  const s = S.summaryOf(b), a = S.analysisOf(b), u = S.unitWords(a.assay), rules = (b.status || {}).rules || [];
  const st = id => { const r = rules.find(q => q.id === id); const k = r ? ruleState(r) : null; return k === 'warn' || k === 'fail' ? k : null; }, gl = id => { const r = rules.find(q => q.id === id); return r && r.outcome !== 'not_evaluated' ? ruleState(r) : null; };
  const o = other ? S.summaryOf(other) : {}, tiles = [], stat = s.unit_stat || 'median';
  const none = (label, why) => <Tile label={label} state="unknown" value="not available" caption={why} />;
  if (a.pipeline === 'spaceranger') {
    tiles.push(<Tile label={S.cap(u.plural) + ' under tissue'} glyph={gl('sp_units_under_tissue')} state={st('sp_units_under_tissue')} value={F.int(s.spatial_units_under_tissue)} caption={F.isNum(s.spatial_fraction_under_tissue) ? F.pct(s.spatial_fraction_under_tissue) + ' of the capture area' : null} />);
    tiles.push(<Tile label={'UMIs per ' + u.name.replace('8 µm ', '') + ', ' + stat} value={S.num1(s.spatial_umis_per_unit)} caption="under tissue, as Space Ranger states it" />);
    tiles.push(<Tile label={'Genes per ' + u.name.replace('8 µm ', '') + ', ' + stat} value={S.num1(s.spatial_genes_per_unit)} caption={F.isNum(s.spatial_genes_detected) ? F.int(s.spatial_genes_detected) + ' genes detected' : null} />);
    tiles.push(<Tile label="Reads under tissue" glyph={gl('sp_fraction_reads_under_tissue')} state={st('sp_fraction_reads_under_tissue')} value={F.pctNum(s.spatial_fraction_reads_under_tissue)} unit="%" caption="of the usable reads; warn below 50 %" />);
    tiles.push(<Tile label="Saturation · Space Ranger" value={F.pctNum(s.sr_saturation)} unit="%" caption={F.isNum(o.saturation_value) ? 'ONT workflow: ' + F.pct(o.saturation_value) : 'reads that repeat a molecule'} />);
    tiles.push(F.isNum(s.sr_usable_reads_fraction) ? <Tile label="Usable reads" glyph={gl('sr_usable_reads_fraction')} state={st('sr_usable_reads_fraction')} value={F.pctNum(s.sr_usable_reads_fraction)} unit="%" caption={'of ' + F.count(s.sr_reads) + ' reads; pass every filter, under tissue'} />
      : <Tile label="Reads" value={F.countParts(s.sr_reads)[0]} unit={F.countParts(s.sr_reads)[1]} caption="given to Space Ranger" />);
  } else {
    const units = s.barcodes_detected != null && a.assay !== 'singlecell' ? s.barcodes_detected : s.barcodes_cells;
    tiles.push(F.isNum(s.reads_input) ? <Tile label={withTag(b, 'reads_input', 'Input reads')} value={F.countParts(s.reads_input)[0]} unit={F.countParts(s.reads_input)[1]} caption="handed to the workflow" /> : none('Input reads', 'only in the workflow’s report, which was not read'));
    // more full-length pieces than reads (a read that holds two molecules is split): a share of the reads cannot be
    // stated and would read "100.0 %"; the tile shows the pieces instead
    const split = F.isNum(s.reads_full_length) && F.isNum(s.reads_input) && s.reads_full_length > s.reads_input;
    tiles.push(split ? <Tile label={withTag(b, 'reads_full_length', 'Full-length pieces')} glyph={gl('sc_full_length_fraction')} state={st('sc_full_length_fraction')} value={F.countParts(s.reads_full_length)[0]} unit={F.countParts(s.reads_full_length)[1]}
        caption={'from ' + F.count(s.reads_input) + ' reads: a read that holds two molecules end to end is split, so there are more pieces than reads'} />
      : F.isNum(s.reads_full_length_fraction) ? <Tile label={withTag(b, 'reads_full_length_fraction', 'Full-length reads')} glyph={gl('sc_full_length_fraction')} state={st('sc_full_length_fraction')} value={F.pctNum(s.reads_full_length_fraction)} unit="%"
      caption={(F.isNum(s.reads_full_length) && F.isNum(s.reads_input) ? F.count(s.reads_full_length) + ' of ' + F.count(s.reads_input) + ' reads · ' : '') + 'both adapters, in order'} /> : none('Full-length reads', 'only in the workflow’s report'));
    tiles.push(F.isNum(s.reads_gene_assigned_fraction) ? <Tile label={withTag(b, 'reads_gene_assigned_fraction', 'Assigned to a gene')} glyph={gl('sc_gene_assigned_fraction')} state={st('sc_gene_assigned_fraction')} value={F.pctNum(s.reads_gene_assigned_fraction)} unit="%" caption="of the full-length reads" /> : none('Assigned to a gene', 'needs the per-read table'));
    tiles.push(<Tile label={a.assay === 'singlecell' ? 'Cells' : S.cap(u.plural) + ' with data'} value={F.int(units)} caption={F.isNum(s.barcodes_reads_per_cell_mean) ? F.int(s.barcodes_reads_per_cell_mean) + ' input reads per ' + u.name : null} />);
    tiles.push(<Tile label={'UMIs per ' + u.name + ', median'} value={S.num1(s.barcodes_umis_median)} caption={F.isNum(s.barcodes_genes_median) ? S.num1(s.barcodes_genes_median) + ' genes, median' : null} />);
    tiles.push(F.isNum(s.saturation_value) ? <Tile label={withTag(b, 'saturation_value', 'Saturation · ONT workflow')} value={F.pctNum(s.saturation_value)} unit="%" caption={F.isNum(o.sr_saturation) ? 'Space Ranger: ' + F.pct(o.sr_saturation) : 'gene-assigned reads that repeat a molecule'} /> : none('Saturation · ONT workflow', 'needs the per-read table'));
  }
  return <div class="fc-keynumbers"><div class="fc-tiles">{tiles}</div></div>;
}

// ---------- sections ----------
function MapSection({ b, label }) {
  const grid = (S.cellOf(b).grids || {}).spatial_map, s = S.summaryOf(b), a = S.analysisOf(b), u = S.unitWords(a.assay);
  const noMap = S.notesOf(b, ['map_needs_positions', 'map_not_made', 'h5py_missing', 'input_unreadable'])[0];
  const reg = S.tableMap(b, 'spatial_registration'), rule = ((b.status || {}).rules || []).find(r => r.id === 'sp_registration' && r.outcome === 'triggered');
  const exact = [];
  if (F.isNum(s.spatial_units_under_tissue)) exact.push([S.cap(u.plural) + ' under tissue', F.int(s.spatial_units_under_tissue)]);
  if (F.isNum(s.spatial_fraction_reads_under_tissue)) exact.push(['Reads under tissue', F.pct1(s.spatial_fraction_reads_under_tissue)]);
  if (F.isNum(s.spatial_fraction_umis_under_tissue)) exact.push(['UMIs under tissue', F.pct1(s.spatial_fraction_umis_under_tissue)]);
  if (F.isNum(s.mito_umi_fraction)) exact.push(['Mitochondrial share of all UMIs', F.pct1(s.mito_umi_fraction)]);
  const side = <>
    {Object.keys(reg).length ? <div class="fc-spmap__box"><b>Image alignment</b><span>{reg.alignment === 'manual alignment file' ? 'Manual alignment file' + (reg.alignment_file ? ' (' + F.midCut(reg.alignment_file, 34, 14) + ')' : '') : reg.alignment ? 'Automatic alignment' : 'As the pipeline recorded it'}.
      {rule ? ' ' + rule.message : reg.is_gaussian_fit_outside_range === 'false' ? ' The pipeline’s check of molecules against the image found no offset.' : ''}</span>{S.fromReport(b, 'spatial_registration') ? <FromReport /> : null}</div> : null}
    <div class="fc-tiles fc-spmap__tiles">
      {F.isNum(s.spatial_fraction_reads_under_tissue) ? <Tile label="Reads under tissue" value={F.pctNum(s.spatial_fraction_reads_under_tissue)} unit="%" /> : null}
      {F.isNum(s.spatial_fraction_umis_under_tissue) ? <Tile label="UMIs under tissue" value={F.pctNum(s.spatial_fraction_umis_under_tissue)} unit="%" /> : null}
    </div></>;
  const reading = a.pipeline === 'spaceranger' ? 'The picture shows where the molecules are and where the tissue is. Exact totals are in the tiles, not read off the picture.'
    : grid && grid.channels.tissue ? 'Molecules and genes per ' + u.name + ' as the ONT workflow counted them. The workflow does not know where the ' + u.plural + ' lie: their places and the tissue outline are taken from the Space Ranger output of the same reads.'
    : 'Molecules and genes per ' + u.name.replace(/^\d+ µm /, '') + ' as the ONT workflow counted them. The workflow does not know where the tissue lies: the tissue outline and image are on the Space Ranger tab.';
  return <Section id="map" title="Spatial map" result={sectionResult(b, 'map')} reading={grid ? reading : null}>
    {grid ? <SpatialMap grid={grid} unit={(a.unit || {}).name} title={'Spatial map of ' + (label || a.label || 'this sample')} exact={exact} side={side} />
      : <Empty title="No spatial map for this analysis.">{noMap ? noMap.message : 'The files that carry the positions of the ' + u.plural + ' were not found.'}</Empty>}
  </Section>;
}

function SaturationSection({ b, other }) {
  const h = S.histOf(b, 'saturation_reads_per_molecule'), s = S.summaryOf(b), a = S.analysisOf(b), sr = a.pipeline === 'spaceranger';
  const stated = sr ? s.sr_saturation : s.saturation_value, word = sr ? 'Space Ranger' : 'the ONT workflow';
  if (!h || !h.w) return <Section id="saturation" title="Saturation" reading={F.isNum(stated) ? F.pct1(stated) + ' of the reads repeat a molecule already seen, as ' + word + ' counts it.' : null}>
    <Empty title="No saturation curve.">{sr ? 'The curve is drawn from molecule_info.h5, which was not found or could not be read.' : 'The curve is drawn from the workflow’s per-read table (one pass over it, a few minutes). It was not read: scan the library again, or run flowcellar analysis without --no-read-pass.'}</Empty></Section>;
  const c = S.saturationCurve(h), reads = c.reads, color = tk('--fc-this-run'), M = v => v >= 1e6 ? F.sig(v / 1e6) + ' M' : F.int(v);
  const end = c.at1, differs = F.isNum(stated) && Math.abs(stated - end) > 0.004;
  const reading = F.pct1(end) + ' of ' + (sr ? 'the reads counted for a gene' : 'the gene-assigned reads with a usable UMI') + ' repeat a molecule already seen.'
    + (c.gain2 != null ? ' Doubling the reads would add about ' + F.pctNum(c.gain2) + ' % more molecules.' : '');
  const ro = k => i => [[M(reads[i]), 'reads'], k === 'sat' ? [F.pct1(c.saturation[i]), 'saturation'] : [M(c.molecules[i]), 'molecules'], [F.pct(reads[i] / c.total), 'of today’s depth']];
  const x = { domain: [0, c.total], kind: 'linear', title: 'Reads, drawn at random from the ' + M(c.total) + ' that were counted', fmt: v => v >= 1e6 ? +(v / 1e6).toPrecision(3) + ' M' : F.int(v) };
  const table = () => ({ head: ['Share of the reads', 'Reads', 'Saturation', 'Unique molecules'], rows: reads.map((r, i) => [F.pct(r / c.total), F.int(r), F.pct1(c.saturation[i]), F.int(c.molecules[i])]).filter((_, i) => (i + 1) % 4 === 0),
    note: 'Computed from how many molecules were seen once, twice, three times and so on: the exact expectation when a random share of the reads is kept.' });
  return <Section id="saturation" title="Saturation" reading={reading}>
    {differs ? <p class="fc-section__note fc-notice">{sr ? 'Space Ranger states ' + F.pct1(stated) + '; the curve ends at ' + F.pct1(end) + ' because it is drawn from the molecule table.'
      : 'The workflow’s count matrix holds more molecules than the reads with a usable UMI resolve (it also counts reads whose UMI could not be corrected), so the saturation on the tile, ' + F.pct1(stated) + ', is lower than the end of this curve, ' + F.pct1(end) + '.'}
      {F.isNum(s.saturation_reported) ? ' The workflow’s own report shows ' + F.pct1(s.saturation_reported) + ', counted over all reads with a valid barcode.' : ''}</p> : null}
    <Grid>
      <ChartPanel title="Saturation against reads" reading={'The curve flattens when more reads mostly repeat molecules. Today’s depth is the right end: ' + F.pct1(end) + '.'} height={200} lines={1}
        spec={{ form: 'xy', x, y: { domain: [0, 1], ticks: [0, 0.25, 0.5, 0.75, 1], fmt: v => Math.round(v * 100) + (v === 1 ? ' %' : ''), title: 'Reads that repeat a molecule' }, series: [{ name: 'Saturation', x: reads, y: c.saturation, color }], cursor: { x: reads }, readout: ro('sat') }} table={table} />
      <ChartPanel title="Unique molecules against reads" reading={'With ' + M(c.total) + ' reads, ' + M(c.molecules1) + ' molecules were found' + (c.gain2 != null ? '; twice the reads would find about ' + M(c.molecules1 * (1 + c.gain2)) + '.' : '.')} height={200} lines={1}
        spec={{ form: 'xy', x, y: { title: 'Unique molecules', fmt: v => v >= 1e6 ? +(v / 1e6).toPrecision(3) + ' M' : F.int(v) }, series: [{ name: 'Molecules', x: reads, y: c.molecules, color }], cursor: { x: reads }, readout: ro('mol') }} table={table} />
    </Grid>
  </Section>;
}

const FATE_WORD = { input: 'Handed to the workflow', full_length: 'Full length (both adapters)', valid_barcode: 'Matched to a barcode', gene_assigned: 'Assigned to a gene', transcript_assigned: 'Assigned to one isoform', unique_molecules: 'Unique molecules',
  raw: 'Given to Space Ranger', feature: 'Counted for a gene', usable: 'Usable: every filter passed, under tissue' };
const ADAPTER_WORD = { full_length: 'Full length', single_adapter1: 'Adapter 1 only', single_adapter2: 'Adapter 2 only', double_adapter1: 'Adapter 1 twice', double_adapter2: 'Adapter 2 twice', no_adapters: 'No adapter', other: 'Other' };
function FateSection({ b }) {
  const c = S.cellOf(b), a = S.analysisOf(b), s = S.summaryOf(b), sr = a.pipeline === 'spaceranger', f = (c.counts || {})[sr ? 'sr_read_funnel' : 'reads_funnel'];
  const ad = (c.counts || {}).reads_adapter_config;
  if (!f) return <Section id="fate" title="Read fate" result={sectionResult(b, 'fate')}><Empty title="Read fate is not known.">{sr ? 'The read funnel comes from molecule_info.h5, which was not found or could not be read.' : 'None of the files that carry the read counts was found.'}</Empty></Section>;
  const missing = new Set(f.missing || []), first = f.n.find((v, i) => !missing.has(f.cats[i])) || 1, src = f.sources || {};
  const rows = f.cats.map((k, i) => ({ label: FATE_WORD[k] || k, value: missing.has(k) ? 0 : f.n[i], text: missing.has(k) ? 'not known' : F.count(f.n[i]), note: missing.has(k) ? '' : k === 'unique_molecules' ? 'molecules' : F.pct(f.n[i] / first) + ' of the first',
    color: k === 'unique_molecules' ? tk('--fc-cat-3') : k === 'usable' || k === 'gene_assigned' ? tk('--fc-fate-usable') : tk('--fc-this-run'), from: src[k] }));
  const last = sr ? 'usable' : 'gene_assigned', iLast = f.cats.indexOf(last);
  const reading = iLast >= 0 && !missing.has(last) && first > 1 ? 'Of ' + F.count(first) + ' reads, ' + F.pct(f.n[iLast] / first) + (sr ? ' passed every Space Ranger filter and lie under tissue.' : ' were assigned to a gene.')
    + (!sr && F.isNum(s.reads_full_length) && F.isNum(s.reads_input) && s.reads_full_length > s.reads_input ? ' Full-length reads can exceed the input: a read that holds two molecules end to end is split.' : '') : null;
  const anyReport = Object.values(src).some(v => v === 'workflow report');
  const adRows = ad ? [{ label: 'Reads', parts: ad.cats.map((k, i) => ({ value: ad.n[i], color: k === 'full_length' ? tk('--fc-fate-usable') : k === 'other' || k === 'no_adapters' ? tk('--fc-fate-no-barcode') : tk('--fc-cat-' + (i + 1)) })) }] : null, adTotal = ad ? ad.n.reduce((t, v) => t + v, 0) : 0;
  return <Section id="fate" title="Read fate" result={sectionResult(b, 'fate')} reading={reading}>
    <Grid>
      <ChartPanel title={sr ? 'Reads through Space Ranger' : 'Reads through the workflow'} tag={anyReport ? 'SOME COUNTS FROM THE REPORT' : null} height={rows.length * 26 + 10} lines={1}
        reading={sr ? 'Each bar is a subset of the one above it.' : 'Each bar of reads is a subset of the one above it; the last bar counts molecules, not reads.'}
        spec={{ form: 'hbars', rows, readout: i => [[rows[i].text, F.lower(rows[i].label)], rows[i].note, rows[i].from ? 'from the ' + rows[i].from : ''].filter(Boolean) }}
        table={() => ({ head: ['Stage', 'Count', 'Share of the first stage', 'Read from'], rows: rows.map((r, i) => [r.label, missing.has(f.cats[i]) ? 'not known' : F.int(f.n[i]), missing.has(f.cats[i]) || f.cats[i] === 'unique_molecules' ? F.DASH : F.pct1(f.n[i] / first), r.from || 'data file']) })} />
      {adRows ? <ChartPanel title="Adapter configuration" tag={S.fromReport(b, 'reads_adapter_config') ? 'FROM THE REPORT' : null} height={22 + 34} lines={1}
        reading={F.pct(ad.n[ad.cats.indexOf('full_length')] / (adTotal || 1)) + ' of the classified reads carry both adapters in the right order.' + (F.isNum(s.reads_plus_fraction) ? ' ' + F.pct(s.reads_plus_fraction) + ' of the stranded reads are in plus orientation.' : '')}
        spec={{ form: 'hstack', norm: true, rowH: 22, rows: adRows, readout: k => { const i = k % 1000; return [[F.pct1(ad.n[i] / (adTotal || 1)), (ADAPTER_WORD[ad.cats[i]] || ad.cats[i]).toLowerCase()], [F.int(ad.n[i]), 'reads']]; } }}
        legend={ad.cats.map((k, i) => ({ name: ADAPTER_WORD[k] || k, color: adRows[0].parts[i].color }))}
        table={() => ({ head: ['Adapter layout', 'Reads', 'Share'], rows: ad.cats.map((k, i) => [ADAPTER_WORD[k] || k, F.int(ad.n[i]), F.pct1(ad.n[i] / (adTotal || 1))]),
          note: S.fromReport(b, 'reads_adapter_config') ? 'Rebuilt from the percentages shown in the workflow’s report: the adapter statistics file itself was not kept.' : undefined })} /> : null}
    </Grid>
    {sr ? <div class="fc-tiles fc-section__tiles">
      <CheckTile b={b} id="sr_valid_barcodes" label="Valid barcodes" v={s.sr_valid_barcodes} caption="reads with a barcode of the slide" />
      <CheckTile b={b} id="sr_valid_umis" label="Valid UMIs" v={s.sr_valid_umis} caption="reads with a usable UMI" />
      {F.isNum(s.sr_corrected_bc_frac) ? <Tile label={withTag(b, 'sr_corrected_bc_frac', 'Barcodes corrected')} value={F.pctNum(s.sr_corrected_bc_frac)} unit="%" caption="needed a correction to match" /> : null}
      {F.isNum(s.sr_low_support_umi_reads_frac) ? <Tile label={withTag(b, 'sr_low_support_umi_reads_frac', 'Low-support UMI reads')} value={F.pctNum(s.sr_low_support_umi_reads_frac)} unit="%" caption="reads in UMIs seen too rarely" /> : null}
    </div> : <div class="fc-tiles fc-section__tiles">
      <CheckTile b={b} id="sc_valid_barcode_fraction" label="Valid barcode rate" v={s.reads_valid_barcode_fraction} k="reads_valid_barcode_fraction" caption="of the full-length reads" />
      <CheckTile b={b} id="sc_valid_umi_fraction" label="Valid UMI rate" v={s.reads_valid_umi_fraction} k="reads_valid_umi_fraction" caption="of the barcoded reads" />
      <CheckTile b={b} id="sc_transcript_assigned_fraction" label="Resolved to one isoform" v={s.reads_transcript_assigned_fraction} k="reads_transcript_assigned_fraction" caption="of the full-length reads" />
      {F.isNum(s.reads_report_unmapped_fraction) ? <Tile label={withTag(b, 'reads_report_unmapped_fraction', 'Unmapped')} value={F.pctNum(s.reads_report_unmapped_fraction)} unit="%" caption="in the workflow’s alignment table; shown, not judged" /> : null}
    </div>}
  </Section>;
}
function CheckTile({ b, id, label, v, k, caption }) {
  if (!F.isNum(v)) return null;
  const r = ((b.status || {}).rules || []).find(q => q.id === id), st = r && r.outcome !== 'not_evaluated' ? ruleState(r) : null;
  return <Tile label={k ? withTag(b, k, label) : label} glyph={st} state={st === 'warn' || st === 'fail' ? st : null} value={F.pctNum(v)} unit="%" caption={caption} />;
}
// the worst result among the checks whose evidence lives in one section
function sectionResult(b, sec) {
  let worst = null; const order = { pass: 0, warn: 1, fail: 2 };
  for (const r of (b.status || {}).rules || []) { if (SECTION_OF[r.id] !== sec || r.role === 'info' || r.outcome === 'not_evaluated') continue; const k = ruleState(r); if (worst == null || order[k] > order[worst]) worst = k; }
  return worst;
}

function MappingSection({ b }) {
  const s = S.summaryOf(b); if (!F.isNum(s.sr_mapped_genome)) return null;
  const conf = s.sr_mapped_confidently_genome || 0, parts = [['Exonic', s.sr_mapped_exonic, '--fc-fate-usable'], ['Intronic', s.sr_mapped_intronic, '--fc-cat-2'], ['Intergenic', s.sr_mapped_intergenic, '--fc-cat-4'],
    ['Mapped, not confidently', Math.max(0, s.sr_mapped_genome - conf), '--fc-cat-6'], ['Not mapped', Math.max(0, 1 - s.sr_mapped_genome), '--fc-fate-no-barcode']].filter(p => F.isNum(p[1]));
  const rows = [{ label: 'Reads', parts: parts.map(p => ({ value: p[1], color: tk(p[2]) })) }];
  const q = [['Barcode', s.sr_q30_barcode], ['UMI', s.sr_q30_umi], ['RNA read', s.sr_q30_rna_read]].filter(x => F.isNum(x[1])), qRows = q.map(x => ({ label: x[0], value: x[1], text: F.pct1(x[1]), note: '' }));
  return <Section id="mapping" title="Mapping" result={sectionResult(b, 'mapping')} reading={F.pct(s.sr_mapped_transcriptome) + ' of the reads map confidently to the transcriptome; ' + F.pct1(s.sr_antisense) + ' map antisense to a gene.'}>
    <div class="fc-tiles fc-section__tiles">
      <CheckTile b={b} id="sr_mapping" label="Confidently to the transcriptome" v={s.sr_mapped_transcriptome} caption="warn below 30 %" />
      <CheckTile b={b} id="sr_antisense_fraction" label="Antisense to a gene" v={s.sr_antisense} caption="warn above 10 %: reads not oriented" />
      <Tile label="Mapped to the genome" value={F.pctNum(s.sr_mapped_genome)} unit="%" />
      {F.isNum(s.reads_chrm_fraction) ? <CheckTile b={b} id="sc_reads_chrm_fraction" label="On the mitochondrial genome" v={s.reads_chrm_fraction} caption="of the aligned reads, from the BAM index" /> : null}
    </div>
    <Grid>
      <ChartPanel title="Where the reads map" reading="Shares of all reads, as Space Ranger reports them." height={22 + 34} lines={1}
        spec={{ form: 'hstack', norm: true, rowH: 22, rows, readout: k => { const i = k % 1000; return [[F.pct1(parts[i][1]), F.lower(parts[i][0])], 'of ' + F.count(s.sr_reads) + ' reads']; } }}
        legend={parts.map(p => ({ name: p[0], color: tk(p[2]) }))} table={() => ({ head: ['Where', 'Share of the reads'], rows: parts.map(p => [p[0], F.pct1(p[1])]).concat([['Confidently to the transcriptome', F.pct1(s.sr_mapped_transcriptome)], ['Antisense to a gene', F.pct1(s.sr_antisense)]]) })} />
      {qRows.length ? <ChartPanel title="Bases at quality 30 or better" reading="As Space Ranger sees the three parts of a read. Nanopore reads score lower here than short reads: context, not a verdict." height={qRows.length * 26 + 10} lines={1}
        spec={{ form: 'hbars', max: 1, rows: qRows, readout: i => [[qRows[i].text, 'of the bases of the ' + F.lower(qRows[i].label)]] }} table={() => ({ head: ['Part of the read', 'Bases at Q30 or better'], rows: q.map(x => [x[0], F.pct1(x[1])]) })} /> : null}
    </Grid>
  </Section>;
}

function RankSection({ b }) {
  const c = S.cellOf(b).curves || {}, k = c.barcodes_knee; if (!k) return null;
  const a = S.analysisOf(b), u = S.unitWords(a.assay), s = S.summaryOf(b), name = Object.keys(k.y)[0], off = c.barcodes_knee_off_tissue, what = k.y_unit === 'umis' ? 'UMIs' : 'Reads';
  const ser = [{ name: name === 'under_tissue' ? 'Under tissue' : S.cap(u.plural), x: k.x, y: k.y[name], color: tk('--fc-this-run') }];
  if (off) ser.push({ name: 'Off tissue', x: off.x, y: off.y.off_tissue, color: tk('--fc-cat-4'), dash: [4, 3] });
  const maxX = Math.max(k.x[k.x.length - 1], off ? off.x[off.x.length - 1] : 0), thr = s.barcodes_cell_threshold;
  const reading = name === 'under_tissue' ? 'The ' + u.plural + ' under tissue hold from ' + F.int(k.y[name][k.y[name].length - 1]) + ' to ' + F.int(k.y[name][0]) + ' ' + what + '.' + (off ? ' ' + S.cap(u.plural) + ' off tissue are drawn dashed: molecules that landed on bare slide.' : '')
    : F.int(k.x[k.x.length - 1]) + ' ' + u.plural + ' have at least one read; the richest has ' + F.int(k.y[name][0]) + '.' + (F.isNum(thr) ? ' A barcode with more than ' + F.int(thr) + ' reads counts as a cell.' : '');
  return <Section id="rank" title="Barcode rank" result={sectionResult(b, 'rank')} reading={reading}>
    <ChartPanel title={what + ' per ' + u.name + ', ranked'} height={240} lines={1} legend={ser.length > 1 ? ser.map(q => ({ name: q.name, kind: 'line', color: q.color, dash: !!q.dash })) : null}
      reading={'Every ' + u.name + ' sorted from the richest to the poorest; both axes are logarithmic.'}
      spec={{ form: 'xy', x: { domain: [1, Math.max(10, maxX)], log: true, kind: 'count', title: S.cap(u.plural) + ', richest first (rank, logarithmic)', fmt: v => F.int(v) }, y: { log: true, title: what + ' per ' + u.name + ', logarithmic' }, series: ser,
        limits: F.isNum(thr) ? [{ y: thr, label: 'cell threshold', kind: 'target' }] : undefined, cursor: { x: k.x }, readout: i => [['Rank ' + F.int(k.x[i]), ''], [F.int(k.y[name][i]), F.lower(what)]] }}
      table={() => ({ head: ['Rank', ser[0].name + ', ' + what], rows: k.x.map((x, i) => [F.int(x), F.int(k.y[name][i])]) })} />
  </Section>;
}

function HistPanel({ b, name, title, what, u, log, pct, c, of }) {
  const h = S.histOf(b, name); if (!h) return null;
  const [a, z] = span(h), e = h.edges, cx = centers(e, !pct), n = total(h) + h.under + h.over, s = S.summaryOf(b);
  const prefix = name.replace('_hist', '').replace('mito_fraction', 'mito_fraction'), med = pct ? s.mito_fraction_median : s[prefix + '_median'], mean = pct ? null : s[prefix + '_mean'];
  const fx = v => pct ? Math.round(v * 100) + ' %' : v >= 100 ? F.int(v) : String(+v.toPrecision(2)), dom = pct ? [0, Math.min(1, e[Math.min(e.length - 1, z + 2)])] : [Math.max(1, e[a]), e[Math.min(e.length - 1, z + 1)]];
  const y = Array.from(h.n), yc = c.logY ? { log: true, title: 'Number of ' + u.plural + ', logarithmic' } : { title: 'Number of ' + u.plural };
  // "11,492 8 µm bins" runs two numbers together: beside a count the unit is written "bins of 8 µm"
  const many = u.plural.replace(/^(\d+ µm) (\w+)$/, '$2 of $1');
  const reading = (F.isNum(med) ? 'Half of ' + (of || 'the ' + u.plural) + ' have ' + (pct ? F.pct1(med) + ' or less' : fmtCount(med) + ' ' + what + ' or fewer') + (F.isNum(mean) ? '; the mean is ' + fmtCount(mean) : '') + '.' : '')
    + (h.under ? ' ' + (h.under === 1 ? 'One ' + u.name + ' has' : F.int(h.under) + ' of the ' + u.plural + ' have') + ' none and ' + (h.under === 1 ? 'is' : 'are') + ' left of the axis.' : '');
  return <ChartPanel title={title} reading={reading} height={190} lines={1} scope={c.logY ? 'log scale' : ''}
    spec={{ form: 'xy', x: pct ? { domain: dom, kind: 'linear', title: 'Mitochondrial share of the ' + u.name + '’s UMIs', fmt: fx } : { domain: dom, log: true, kind: 'count', title: S.cap(what) + ' per ' + u.name + ', logarithmic', fmt: v => F.int(v) }, y: yc,
      series: [{ name: S.cap(u.plural), type: 'hist', edges: e, y }], vlines: F.isNum(med) && med >= dom[0] ? [{ x: Math.max(med, dom[0]), label: 'median ' + (pct ? F.pct1(med) : fmtCount(med)), strong: true }] : [], cursor: { x: cx },
      readout: i => [[pct ? fx(e[i]) + ' to ' + fx(e[i + 1]) : fmtCount(e[i]) + ' to ' + fmtCount(e[i + 1]), pct ? '' : what], [F.int(h.n[i]), many], [F.pct1(h.n[i] / (n || 1)), 'of them']] }}
    table={() => ({ head: [pct ? 'Mitochondrial share' : S.cap(what) + ' per ' + u.name, S.cap(u.plural), 'Share'], rows: (h.under ? [['none', F.int(h.under), F.pct1(h.under / n)]] : []).concat(Array.from(h.n).map((v, i) => [pct ? fx(e[i]) + ' to ' + fx(e[i + 1]) : fmtCount(e[i]) + ' to ' + fmtCount(e[i + 1]), F.int(v), F.pct1(v / (n || 1))]).filter((_, i) => i >= a && i <= z)) })} />;
}
const fmtCount = v => v >= 100 ? F.int(v) : String(+(+v).toPrecision(3));
function UnitsSection({ b, c }) {
  const a = S.analysisOf(b), u0 = S.unitWords(a.assay), u = { name: u0.name, plural: u0.plural }, s = S.summaryOf(b), under = a.pipeline === 'spaceranger';
  const has = n => !!S.histOf(b, n); if (!has('barcodes_umis_hist') && !has('barcodes_genes_hist') && !has('barcodes_reads_hist') && !has('mito_fraction_hist')) return null;
  // The median mitochondrial share of the ONT workflow is the median of the workflow's own table, which holds only the
  // units it kept (those with at least "matrix_min_genes" genes): the tile and the chart say of how many, and of which.
  const mh = S.histOf(b, 'mito_fraction_hist'), mitoN = mh ? total(mh) + mh.under + mh.over : null, withData = F.isNum(s.barcodes_detected) ? s.barcodes_detected : s.barcodes_cells;
  const minGenes = +(((S.tableOf(b, 'workflow_parameters') || { rows: [] }).rows.find(r => r[0] === 'matrix_min_genes') || [])[1]) || null;
  const kept = !under && mitoN != null && F.isNum(withData) && mitoN < withData, mitoOf = kept ? 'of the ' + F.int(mitoN) + ' ' + u.plural + ' the workflow kept' + (minGenes > 1 ? ' (' + F.int(minGenes) + ' genes or more)' : '') + ', of ' + F.int(withData) + ' with a molecule · ' : '';
  const reading = (under ? S.cap(u.plural) + ' under tissue. ' : S.cap(u.plural) + ' with at least one molecule. ') + (F.isNum(s.barcodes_umis_median) && F.isNum(s.barcodes_genes_median) ? 'The median ' + u.name + ' holds ' + fmtCount(s.barcodes_umis_median) + ' UMIs of ' + fmtCount(s.barcodes_genes_median) + ' genes.' : '');
  return <Section id="units" title={'UMIs and genes per ' + u.name} result={sectionResult(b, 'units')} reading={reading}
    tools={<Seg name="an-y" label="Scale" value={c.logY ? 'log' : 'lin'} onChange={v => setParams({ y: v === 'lin' ? '' : v })} opts={[{ v: 'lin', label: 'Linear' }, { v: 'log', label: 'Log' }]} />}>
    <Grid>
      <HistPanel b={b} c={c} u={u} name="barcodes_umis_hist" title={'UMIs per ' + u.name} what="UMIs" />
      <HistPanel b={b} c={c} u={u} name="barcodes_genes_hist" title={'Genes per ' + u.name} what="genes" />
      <HistPanel b={b} c={c} u={u} name="barcodes_reads_hist" title={'Reads per ' + u.name} what="reads" />
      <HistPanel b={b} c={c} u={u} name="mito_fraction_hist" title={'Mitochondrial share per ' + u.name} what="share" pct of={kept ? 'the ' + F.int(mitoN) + ' ' + u.plural + ' the workflow kept' + (minGenes > 1 ? ' (' + F.int(minGenes) + ' genes or more)' : '') : null} />
    </Grid>
    <div class="fc-tiles fc-section__tiles">
      <CheckTile b={b} id="sc_mito_fraction_per_barcode" label={'Mitochondrial share, median ' + u.name} v={s.mito_fraction_median} caption={mitoOf + 'warn above 25 %; never fails'} />
      {F.isNum(s.mito_umi_fraction) ? <Tile label="Mitochondrial share of all UMIs" value={F.pctNum(s.mito_umi_fraction)} unit="%" caption={'every ' + u.name + ' with a molecule'} /> : null}
      {a.pipeline !== 'spaceranger' ? <CheckTile b={b} id="sc_reads_chrm_fraction" label="Reads on the mitochondrial genome" v={s.reads_chrm_fraction} caption="of the aligned reads, from the BAM index" /> : null}
      {F.isNum(s.barcodes_genes_detected) ? <Tile label="Genes detected" value={F.int(s.barcodes_genes_detected)} caption={F.isNum(s.barcodes_transcripts_detected) ? F.int(s.barcodes_transcripts_detected) + ' isoforms' : null} /> : null}
    </div>
  </Section>;
}

const SEG_WORD = { cells: 'Cells', reads_in_cells: 'Reads in cells', umis_in_cells: 'UMIs in cells', mean_reads_per_cell: 'Reads per cell, mean', median_genes_per_cell: 'Genes per cell, median', median_umis_per_cell: 'UMIs per cell, median',
  median_cell_area_um2: 'Cell area, median', median_nucleus_area_um2: 'Nucleus area, median' };
function BinsSection({ b }) {
  const t = S.tableOf(b, 'spatial_bin_levels'), seg = S.tableOf(b, 'spatial_segmentation'), s = S.summaryOf(b); if (!t && !seg) return null;
  const segV = (k, v) => /in_cells$/.test(k) ? [F.pctNum(v), '%'] : /area/.test(k) ? [S.num1(v), 'µm²'] : k === 'cells' ? [F.int(v), ''] : [S.num1(v), ''];
  return <Section id="bins" title="Bin sizes and cells" result={sectionResult(b, 'bins')} reading={t ? 'The same four numbers at each resolution Space Ranger reports. The page uses the 8 µm bins.' : null}>
    {t ? <div class="fc-tablewrap"><table class="fc-table fc-table--comfortable"><caption class="fc-sr">Numbers under tissue at each bin size</caption>
      <thead><tr><th scope="col">Bin size</th><th scope="col" class="fc-num">Under tissue</th><th scope="col" class="fc-num">Share of the area</th><th scope="col" class="fc-num">Reads per bin, mean</th><th scope="col" class="fc-num">UMIs per bin, mean</th><th scope="col" class="fc-num">Genes per bin, mean</th><th scope="col" class="fc-num">Genes detected</th></tr></thead>
      <tbody>{t.rows.map(r => <tr><th scope="row">{r[0]} µm {r[0] === 2 ? 'squares' : 'bins'}</th><td class="fc-num">{F.int(r[1])}</td><td class="fc-num">{F.pct1(r[2])}</td><td class="fc-num">{S.num1(r[3])}</td><td class="fc-num">{S.num1(r[4])}</td><td class="fc-num">{S.num1(r[5])}</td><td class="fc-num">{F.int(r[6])}</td></tr>)}</tbody></table></div> : null}
    {seg ? <div class="fc-tiles fc-section__tiles">{seg.rows.map(r => r[0] === 'umis_in_cells' ? <CheckTile b={b} id="sp_segmentation_summary" label={SEG_WORD[r[0]]} v={r[1]} caption="of the molecules under tissue; warn below 50 %" />
      : <Tile label={SEG_WORD[r[0]] || r[0]} value={segV(r[0], r[1])[0]} unit={segV(r[0], r[1])[1]} />)}</div> : null}
    {F.isNum(s.spatial_umis_per_mm2) ? <p class="fc-section__note">{F.count(s.spatial_umis_per_mm2)} UMIs and {F.count(s.spatial_reads_per_mm2)} reads per mm² of tissue: comparable between sections of different size.</p> : null}
  </Section>;
}

function PairSection({ wf, sr }) {
  if (!wf || !sr) return null;
  const a = S.summaryOf(wf), z = S.summaryOf(sr), u = S.unitWords(S.analysisOf(sr).assay), ratio = F.isNum(a.barcodes_unique_molecules) && z.barcodes_unique_molecules ? a.barcodes_unique_molecules / z.barcodes_unique_molecules : null;
  const rule = ((wf.status || {}).rules || []).find(r => r.id === 'sp_wf_sr_umi_ratio');
  const rows = [['Unique molecules', F.int(a.barcodes_unique_molecules), F.int(z.barcodes_unique_molecules)], ['Saturation', S.pct(a.saturation_value), S.pct(z.sr_saturation)],
    ['Reads matched to a barcode', F.isNum(a.reads_valid_barcode) ? F.int(a.reads_valid_barcode) : F.DASH, F.isNum(z.sr_valid_barcodes) && F.isNum(z.sr_reads) ? F.int(z.sr_valid_barcodes * z.sr_reads) : F.DASH],
    ['UMIs per ' + u.name, F.isNum(a.barcodes_umis_median) ? fmtCount(a.barcodes_umis_median) + ' median, every ' + u.name + ' with data' : F.DASH, F.isNum(z.barcodes_umis_median) ? fmtCount(z.barcodes_umis_median) + ' median, under tissue' : F.DASH],
    ['Genes per ' + u.name, F.isNum(a.barcodes_genes_median) ? fmtCount(a.barcodes_genes_median) + ' median' : F.DASH, F.isNum(z.barcodes_genes_median) ? fmtCount(z.barcodes_genes_median) + ' median' : F.DASH],
    ['Reads on the mitochondrial genome', S.pct(a.reads_chrm_fraction), S.pct(z.reads_chrm_fraction)]];
  return <Section id="pair" title="Two pipelines compared" result={rule && rule.outcome !== 'not_evaluated' ? ruleState(rule) : null}
    reading={ratio != null ? 'The same reads counted two ways: the ONT workflow finds ' + ratio.toFixed(2) + ' times the molecules Space Ranger finds. ' + (ratio > 2 || ratio < 0.5 ? 'A gap this large means one of the two pipelines is losing reads.' : 'Full-length long reads usually give somewhat more.') : 'The same reads counted two ways.'}>
    <div class="fc-tablewrap"><table class="fc-table fc-table--comfortable fc-table--wrap"><caption class="fc-sr">The ONT workflow and Space Ranger side by side</caption>
      <thead><tr><th scope="col">Number</th><th scope="col" class="fc-num">ONT workflow</th><th scope="col" class="fc-num">Space Ranger</th></tr></thead>
      <tbody>{rows.map(r => <tr><th scope="row">{r[0]}</th><td class="fc-num">{r[1]}</td><td class="fc-num">{r[2]}</td></tr>)}</tbody></table></div>
    <p class="fc-section__note">The two saturation figures are not interchangeable: each pipeline keeps a different set of reads and counts molecules its own way.</p>
  </Section>;
}

const PARAM_WORD = { kit: 'Kit', full_length_only: 'Full-length reads only', min_read_qual: 'Minimum read quality', barcode_max_ed: 'Barcode: most edits', barcode_min_ed_diff: 'Barcode: least gap to the next match', barcode_min_quality: 'Barcode: minimum quality',
  gene_assigns_minqv: 'Gene assignment: minimum mapping quality', matrix_min_genes: 'Matrix: minimum genes per unit', matrix_min_cells: 'Matrix: minimum units per gene', matrix_max_mito: 'Matrix: most mitochondrial %', mito_prefix: 'Mitochondrial gene prefix',
  expected_cells: 'Expected cells', estimate_cell_count: 'Estimate the cell count', ref_genome_dir: 'Reference' };
function WorkflowSection({ b }) {
  const s = S.summaryOf(b), top = S.tableOf(b, 'workflow_execution'), par = S.tableOf(b, 'workflow_parameters'), sample = S.tableOf(b, 'spatial_sample'); if (!top && !par && !sample && !F.isNum(s.workflow_tasks)) return null;
  const mem = v => F.isNum(v) ? F.bytes(v) : F.DASH;
  return <Section id="workflow" title="Workflow" result={sectionResult(b, 'workflow')} reading={F.isNum(s.workflow_tasks) ? (s.workflow_failed_tasks ? F.plural(s.workflow_failed_tasks, 'task') + ' of ' + F.int(s.workflow_tasks) + ' failed.' : 'All ' + F.int(s.workflow_tasks) + ' tasks finished') + (s.workflow_failed_tasks ? '' : F.isNum(s.workflow_wall_hours) ? ' in ' + F.hours(s.workflow_wall_hours) + '.' : '.') : null}>
    {F.isNum(s.workflow_tasks) ? <div class="fc-tiles fc-section__tiles">
      <CheckTile0 b={b} id="sc_wf_execution" label="Tasks" value={F.int(s.workflow_tasks)} caption={s.workflow_failed_tasks ? F.plural(s.workflow_failed_tasks, 'failed task') : s.workflow_cached_tasks ? F.int(s.workflow_cached_tasks) + ' taken from an earlier run' : 'none failed'} />
      <Tile label="Wall time" value={F.hours(s.workflow_wall_hours)} caption="from the first task to the last" />
      <Tile label="Task time" value={F.hours(s.workflow_task_hours)} caption="all tasks added up" />
      <Tile label="Peak memory" value={mem(s.workflow_peak_rss_bytes)} caption="of the largest task" />
    </div> : null}
    {top ? <div class="fc-tablewrap"><table class="fc-table fc-table--comfortable"><caption class="fc-sr">The five steps that took the longest</caption>
      <thead><tr><th scope="col">Step, slowest first</th><th scope="col" class="fc-num">Tasks</th><th scope="col" class="fc-num">Time</th><th scope="col" class="fc-num">Peak memory</th></tr></thead>
      <tbody>{top.rows.map(r => <tr><th scope="row" class="fc-mono">{r[0]}</th><td class="fc-num">{F.int(r[1])}</td><td class="fc-num">{F.hours(r[2])}</td><td class="fc-num">{mem(r[3])}</td></tr>)}</tbody></table></div> : null}
    {par ? <div class="fc-tablewrap fc-anparams"><table class="fc-table fc-table--comfortable fc-table--wrap"><caption class="fc-sr">Settings that change what the numbers mean</caption>
      <thead><tr><th scope="col">Setting{S.fromReport(b, 'workflow_parameters') ? <> <FromReport /></> : null}</th><th scope="col">Value</th></tr></thead>
      <tbody>{par.rows.map(r => <tr><th scope="row" title={r[0]}>{PARAM_WORD[r[0]] || r[0].replace(/[-_]/g, ' ')}</th><td class="fc-mono fc-break">{/^true$/i.test(String(r[1])) ? 'yes' : /^false$/i.test(String(r[1])) ? 'no' : String(r[1])}</td></tr>)}</tbody></table></div> : null}
    {sample ? <div class="fc-tablewrap fc-anparams"><table class="fc-table fc-table--comfortable fc-table--wrap"><caption class="fc-sr">Sample table of the pipeline</caption>
      <thead><tr><th scope="col">Sample table <FromReport /></th><th scope="col">Value</th></tr></thead>
      <tbody>{sample.rows.map(r => <tr><th scope="row">{r[0]}</th><td class="fc-break">{String(r[1])}</td></tr>)}</tbody></table></div> : null}
  </Section>;
}
function CheckTile0({ b, id, label, value, caption }) {
  const r = ((b.status || {}).rules || []).find(q => q.id === id), st = r && r.outcome !== 'not_evaluated' ? ruleState(r) : null;
  return <Tile label={label} glyph={st} state={st === 'warn' || st === 'fail' ? st : null} value={value} caption={caption} />;
}

const STATE_WORD = { ok: 'read', skipped: 'not read', missing: 'not found', error: 'could not be read', partially_read: 'read in part', truncated: 'cut off' };
function SourcesSection({ b }) {
  const inputs = (b.provenance || {}).inputs || [], sec = S.section(b) || {}, src = sec.value_sources || {}, from = Object.keys(src).filter(k => src[k].kind === 'html_table' || src[k].kind === 'embedded_json');
  const notes = (sec.notes || []).filter(n => n.code !== 'pipeline_alarm'), a = S.analysisOf(b);
  return <Section id="sources" title="Data sources" reading={inputs.length + ' files were read for this page.' + (from.length ? ' ' + F.plural(from.length, 'value') + ' came from a report’s HTML because no data file holds ' + (from.length === 1 ? 'it' : 'them') + '.' : ' Every number comes from a data file.')}>
    {notes.map(n => <p class={'fc-section__note fc-notice' + (sec.state === 'partial' ? ' fc-partial' : '')}>{n.message}</p>)}
    <div class="fc-tablewrap"><table class="fc-table fc-table--comfortable fc-table--wrap"><caption class="fc-sr">Files read for this analysis</caption>
      <thead><tr><th scope="col">File</th><th scope="col">What it is</th><th scope="col" class="fc-num">Size</th><th scope="col">Written</th><th scope="col">State</th></tr></thead>
      <tbody>{inputs.map(i => <tr><th scope="row" title={i.path}><span class="fc-mono fc-break">{(i.path || '').split('/').slice(-2).join('/')}</span>{(i.notes || []).map(n => <span class="fc-table__sub">{n.message}</span>)}</th>
        <td>{S.fileWord(i.role)}</td><td class="fc-num fc-nowrap">{F.bytes(i.size)}</td><td class="fc-nowrap">{i.mtime ? <When iso={i.mtime} /> : F.DASH}</td><td>{STATE_WORD[i.state] || i.state}{F.isNum(i.records) ? ' · ' + F.int(i.records) + (i.records === 1 ? ' record' : ' records') : ''}</td></tr>)}</tbody></table></div>
    {from.length ? <p class="fc-section__note">Read from a report’s HTML, because no data file holds {from.length === 1 ? 'it' : 'them'}: {[...new Set(from.map(k => S.valueWord(k)))].join(', ')}.</p> : null}
    {a.path ? <p class="fc-section__note">Folder: <span class="fc-mono fc-break">{a.path}</span></p> : null}
  </Section>;
}

// ---------- rail ----------
function Rail({ items, current }) {
  return <nav class="fc-rail fc-runrail" aria-label="Sections of this analysis">
    {items.map(it => it.group ? <div class="fc-rail__group fc-caps">{it.group}</div>
      : <a class={'fc-rail__item' + (current === it.id ? ' is-on' : '')} href={makeHash(route.value.path, { ...route.value.q, at: it.id === 'summary' ? '' : it.id })} aria-current={current === it.id ? 'location' : undefined}
          onClick={e => { e.preventDefault(); setParams({ at: it.id === 'summary' ? '' : it.id }); (it.id === 'summary' ? window.scrollTo(0, 0) : document.getElementById('s-' + it.id)?.scrollIntoView({ behavior: 'instant', block: 'start' })); }}>
        <span>{it.label}</span><span class="fc-rail__mark">{it.result ? <Glyph v={it.result} /> : it.mark || ''}</span></a>)}
  </nav>;
}

// ---------- the page ----------
function crumbs(label) {
  if (isReport) return <Crumbs>Analysis report {'›'} <b>{label}</b></Crumbs>;
  return <Crumbs><a href="#/home">Library</a> {'›'} <a href="#/analyses">Analyses</a> {'›'} <b>{label}</b></Crumbs>;
}
export function Analysis({ id }) {
  const ix = index.value, rows = (ix && ix.analyses) || [], byId = (ix && ix.analysisById) || new Map();
  // an address may name a sample (its slide and area) instead of one analysis: the newest analysis of that sample
  const asked = byId.get(id) || rows.find(r => r.group === id) || null;
  const fam = asked ? S.family(rows, byId, asked.id) : [];
  const want = param('p', ''), shown = fam.find(r => r.pipeline === want) || asked, other = fam.find(r => r !== shown) || null;
  const b = useBundle(shown ? shown.id : null), bo = useBundle(other ? other.id : null), state = shown ? bundleState(shown.id) : 'absent', err = shown ? bundleError(shown.id) : null;
  const [current, setCurrent] = useState('summary'), at = param('at', '');
  useEffect(() => { if (!b || !at) return; const t = setTimeout(() => { const el = document.getElementById('s-' + at); if (el) el.scrollIntoView({ behavior: 'instant', block: 'start' }); }, 80); return () => clearTimeout(t); }, [b && b.bundle_id]);
  useEffect(() => {
    if (!b) return; const secs = [...document.querySelectorAll('[data-sec]')];
    const spy = () => { let cur = 'summary'; for (const s of secs) if (s.getBoundingClientRect().top < 180) cur = s.getAttribute('data-sec'); setCurrent(cur); };
    addEventListener('scroll', spy, { passive: true }); return () => removeEventListener('scroll', spy);
  }, [b && b.bundle_id]);
  if (!asked) return <div class="fc-page">{crumbs(id)}<div class="fc-empty fc-page__empty"><b>This analysis is not in the library.</b><span>No analysis has the id {id}. {isReport ? null : <a href="#/analyses">All analyses</a>}</span></div></div>;
  const lead = fam.find(r => r.pipeline === 'spaceranger') || asked, label = lead.label || asked.label || asked.id, u = S.unitWords(asked.assay);
  const links = S.linksOf(fam), same = rows.filter(r => r.group && r.group === asked.group && !fam.includes(r));
  const head = <header class="fc-runhead">
    <div class="fc-runhead__row">
      <h1 class="fc-runhead__id"><Mid text={label} tail={14} /></h1>
      <span class="fc-runhead__chips">
        <Chip>{S.ASSAY_WORD[asked.assay] || 'Assay not known'}{lead.chemistry && !/visium/i.test(S.ASSAY_WORD[asked.assay] || '') ? ' · ' + lead.chemistry : ''}</Chip>
        {lead.chemistry ? <Chip title="Chemistry, as the pipeline recorded it">{lead.chemistry}</Chip> : asked.kit ? <Chip title="Kit, as given to the workflow">{asked.kit}</Chip> : null}
        {lead.slide ? <Chip title="Slide serial number and capture area">slide {lead.slide}</Chip> : null}
        {(lead.reference || asked.reference) ? <Chip title="Reference">{lead.reference || asked.reference}</Chip> : null}
        {fam.map(r => <Chip title={r.variant ? 'Output folder: ' + r.variant : undefined}>{S.PIPELINE_WORD[r.pipeline] || r.pipeline}{r.version ? ' ' + r.version : ''}</Chip>)}
      </span>
      {b ? <span class="fc-runhead__actions"><SavePdf /></span> : null}
    </div>
    <p class="fc-runhead__line fc-anruns"><span>{links.length ? 'Reads from ' + F.plural(links.length, 'flow cell run') + (links.length > 1 ? ', merged before the analysis' : '') + ':' : 'Flow cell runs:'}</span> <RunChips links={links} /></p>
    {same.length && !isReport ? <p class="fc-runhead__line">Other analyses of this sample: {same.map((r, k) => <>{k ? ', ' : ''}<a href={S.hrefAnalysis(r.id)}>{(S.PIPELINE_WORD[r.pipeline] || r.pipeline) + (r.variant ? ' · ' + r.variant : '')}</a></>)}</p> : null}
  </header>;
  if (!b) return <div class="fc-page fc-run fc-analysis">{crumbs(label)}{head}
    <section class="fc-verdict"><div class="fc-verdict__head"><Badge lg v={shown.verdict} /><div><p class="fc-verdict__headline">{shown.headline}</p></div></div>
      {(shown.reasons || []).length ? <ul class="fc-reasons">{shown.reasons.map((t, k) => <li><Glyph v={shown.reason_levels[k] === 'fail' ? 'fail' : 'warn'} /><span>{t}</span></li>)}</ul> : null}</section>
    {state === 'failed' ? <div class={'fc-banner fc-runfail' + (err && err.damaged ? ' fc-banner--strong' : '')} role="alert"><span>{err && err.damaged ? 'The data of this analysis is damaged: its checksum does not match, so nothing of it is drawn.' : 'The data file of this analysis could not be read (' + (err ? err.message : 'unknown reason') + '). The summary above comes from the library index.'}</span>
      {isReport ? null : <button type="button" class="fc-btn" onClick={() => retryBundle(shown.id)}>Try again</button>}</div> : <p class="fc-loading">Loading the map and the charts of this analysis …</p>}</div>;
  const c = { logY: param('y', 'lin') === 'log' }, a = S.analysisOf(b), sr = a.pipeline === 'spaceranger', res = s => sectionResult(b, s), cell = S.cellOf(b);
  const wfB = sr ? bo : b, srB = sr ? b : bo, spatial = a.assay === 'visium' || a.assay === 'visium_hd';
  const items = [{ id: 'summary', label: 'Summary', result: b.status && ['pass', 'warn', 'fail'].includes(b.status.verdict) ? b.status.verdict : null }];
  if (spatial) items.push({ group: 'Tissue' }, { id: 'map', label: 'Spatial map', result: res('map') });
  if (S.tableOf(b, 'spatial_bin_levels') || S.tableOf(b, 'spatial_segmentation')) items.push({ id: 'bins', label: 'Bin sizes and cells', result: res('bins') });
  items.push({ group: 'Reads' }, { id: 'fate', label: 'Read fate', result: res('fate') });
  if (sr && F.isNum(S.summaryOf(b).sr_mapped_genome)) items.push({ id: 'mapping', label: 'Mapping', result: res('mapping') });
  items.push({ group: 'Per ' + u.name.replace(/^\d+ µm /, '') });      // the heading is set in capitals: a "µ" would read as "M"
  if ((cell.curves || {}).barcodes_knee) items.push({ id: 'rank', label: 'Barcode rank' });
  items.push({ id: 'units', label: 'UMIs and genes per ' + u.name.replace('8 µm ', ''), result: res('units') }, { id: 'saturation', label: 'Saturation' });
  items.push({ group: 'Details' });
  if (wfB && srB) items.push({ id: 'pair', label: 'Two pipelines compared', result: res('pair') });
  items.push({ id: 'workflow', label: 'Workflow', result: res('workflow') }, { id: 'sources', label: 'Data sources', mark: String(((b.provenance || {}).inputs || []).length) });
  const tabWord = r => S.PIPELINE_TAB[r.pipeline] || r.pipeline;
  return <div class="fc-page fc-run fc-analysis">
    <PaperHead kind="Analysis report" bundles={[b]} file={'FlowCellar analysis ' + label} />
    {crumbs(label)}
    {head}
    <div class="fc-runtop"><Verdict b={b} runs={links.length} /><KeyNumbers b={b} other={bo} /></div>
    <PaperScope items={[['Pipeline', tabWord(shown) + (other ? ' (the numbers of ' + tabWord(other) + ' are on its own tab and not in this PDF)' : '')], S.histOf(b, 'barcodes_umis_hist') || S.histOf(b, 'barcodes_genes_hist') ? ['Scale of the counts per ' + u.name, c.logY ? 'logarithmic' : 'linear'] : null]} />
    {fam.length > 1 ? <div class="fc-controls fc-runcontrols fc-ancontrols" role="group" aria-label="Which pipeline's numbers are shown">
      <Seg name="an-p" label="Pipeline" value={shown.pipeline} onChange={v => setParams({ p: v === asked.pipeline ? '' : v, at: '' })} opts={fam.slice().sort((x, y) => (x.pipeline === 'spaceranger' ? 1 : 0) - (y.pipeline === 'spaceranger' ? 1 : 0)).map(r => ({ v: r.pipeline, label: S.PIPELINE_TAB[r.pipeline] || r.pipeline }))} />
      <span class="fc-ancontrols__verdicts">{fam.slice().sort((x, y) => (x.pipeline === 'spaceranger' ? 1 : 0) - (y.pipeline === 'spaceranger' ? 1 : 0)).map(r => <span class="fc-ancontrols__v">{S.PIPELINE_TAB[r.pipeline] || r.pipeline}: <Badge v={r.verdict} /></span>)}</span>
      <span class="fc-ancontrols__note">{sr ? 'Space Ranger’s numbers: tissue, bins, mapping. ' : 'The ONT workflow’s numbers: full-length reads, genes and isoforms. '}The two are made from the same reads; each has a verdict of its own.</span>
    </div> : null}
    <div class="fc-runbody">
      <Rail items={items} current={current} />
      <main class="fc-runmain">
        {spatial ? <MapSection b={b} label={label} /> : null}
        <BinsSection b={b} />
        <FateSection b={b} />
        {sr ? <MappingSection b={b} /> : null}
        <RankSection b={b} />
        <UnitsSection b={b} c={c} />
        <SaturationSection b={b} other={bo} />
        <PairSection wf={wfB} sr={srB} />
        <WorkflowSection b={b} />
        <SourcesSection b={b} />
      </main>
    </div>
  </div>;
}
