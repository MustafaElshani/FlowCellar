// Read data: what linking basecalled reads and the pipeline's result files adds to a run. The numbers come from
// every read of the linked files (exact); that is the normal case and carries no mark. When a section rests on part
// of its data (not every sample, a file read in part, files still waiting) it says so in a line under its title,
// in a tag on every chart and in the captions of its tiles. A stage nobody counted is "not known", never zero.
import { setParams, hrefRun } from '../core/state.js';
import { index, isReport } from '../core/data.js';
import { Tile, Mid } from '../shell/parts.jsx';
import { ChartPanel, ScaleLegend, Empty } from '../charts/ChartPanel.jsx';
import { tk } from '../charts/theme.js';
import * as F from '../core/format.js';
import * as D from '../core/derive.js';
import { Section, Disclosure, AllNumbers, countNumbers } from './section.jsx';
import { sectionResult, scopeTag, ruleById, partialInfo, readSetWords, sampleName, chosenWords, waitingOnly, variantNote, sharedTables, ofReads } from './common.js';
import { AlignmentOpening, AlignmentMore, Calibration, SamplesChart } from './alignment.jsx';
import { Grid } from '../charts/Grid.jsx';

const slot = (s, i) => { const m = /(\d+)$/.exec((s && s.barcode) || ''); return ((m ? +m[1] - 1 : i) % 12) + 1; };
export const sampleColor = (s, i) => tk('--fc-cat-' + slot(s, i));
const centersLog = e => Array.from({ length: e.length - 1 }, (_, i) => Math.sqrt(e[i] * e[i + 1]));
const centers = e => Array.from({ length: e.length - 1 }, (_, i) => (e[i] + e[i + 1]) / 2);
const lower = s => s ? s.charAt(0).toLowerCase() + s.slice(1) : s;
const tagOf = (st, p) => [st.tag, p ? p.tag : null].filter(Boolean).join(' · ');
const Partial = ({ p }) => p ? <p class="fc-section__note fc-notice fc-partial" role="status">{p.text}</p> : null;
// a section whose files are linked but not read, and that holds no number yet: the notice and the command, nothing else
const Waiting = ({ w }) => <p class="fc-section__note fc-notice fc-partial fc-waiting" role="status">Partial: {w.text} The next scan reads them, or now: <span class="fc-mono fc-break">{w.command}</span></p>;
const Variant = ({ text }) => text ? <p class="fc-section__note fc-notice fc-variant" role="status">{text}</p> : null;
// a section that has numbers, but none for the samples chosen in the Samples control
const NoneFor = ({ c, what }) => <Empty title={'No ' + what + ' for ' + chosenWords(c)}>The files read so far hold no {what} of {c.sKeys.length > 1 ? 'these samples' : 'this sample'}: other sections of the page may. <button type="button" class="fc-linkbtn" onClick={() => setParams({ s: '' })}>Show all samples</button></Empty>;

export function ReadDataCard({ b }) {
  const y = D.summary(b, 'yield'), bar = b.identity.assay.barcoded, g = b.identity.assay.group;
  const items = [['Read fate', F.int(y.yield_reads) + ' reads sequenced → basecalled → ' + (bar ? 'given to a sample → ' : '') + (g === 'cdna' ? 'full-length → ' : '') + 'mapped'],
    bar ? ['Samples', 'reads and bases per barcode'] : null, ['Read length and quality', 'read N50 and Q-scores from every basecalled read, in place of the estimates above'], ['Alignment', 'summary of read alignment, base breakdown, depth along the reference, identity'],
    ['Reads over run time', 'output, read length, quality and speed hour by hour, from the reads themselves'],
    g === 'cdna' || g === 'drna' ? ['Transcripts', 'poly(A) tails, coverage along transcripts, transcripts and genes found'] : null,
    g === 'cdna' ? ['Pychopper', 'the full-length share and everything else pychopper recorded, per sample'] : null,
    g === 'cdna' || g === 'drna' ? ['RNA composition', 'reads by transcript type, tail length by type and how the reads spread over genes, from the annotated read tables'] : null].filter(Boolean);
  return <section class="fc-section" id="s-read_data" data-sec="read_data" aria-labelledby="h-read_data">
    <div class="fc-section__head"><h2 class="fc-section__title" id="h-read_data">Read data</h2><span class="fc-badge fc-badge--unknown">Not available yet</span></div>
    <div class="fc-empty fc-readcard">
      <b>These sections appear when basecalled reads are linked to this run</b>
      <dl class="fc-readcard__grid">{items.map(([t, d]) => <div><dt>{t}</dt><dd>{d}</dd></div>)}</dl>
      <span>How to get them: basecall the run, then run <span class="fc-mono">flowcellar scan</span> again on the folder that holds the basecalled files. FlowCellar matches files to this run by the IDs inside them, not by folder names. Until then this is a state, not a warning: it never lowers the verdict.</span>
    </div>
  </section>;
}

export function ReadFate({ b, c, evRule, evSec }) {
  const fu = D.funnel(b); if (!fu || !fu.length) return null;
  const fs = D.funnelSamples(b), sm = D.summary(b, 'funnel'), seq = fu[0].n, known = fu.filter(s => !s.missing), last = known[known.length - 1], missing = fu.filter(s => s.missing), p = partialInfo(b, 'funnel');
  const rows = fu.map(s => s.missing ? { label: s.label, unknown: true } : { label: s.label, value: s.n, text: F.count(s.n), note: s.id === 'sequenced' ? 'reads' : F.pct1(s.n / seq) + ' of sequenced' });
  const reading = sm.funnel_usable_frac != null && last.id !== 'sequenced' ? 'Of every 100 reads the flow cell produced, ' + Math.round(sm.funnel_usable_frac * 100) + ' ended up ' + lower(last.label) + '.'
    : last.id === 'sequenced' ? 'Read fate is known only up to sequencing: no count of the reads that reached a sample is linked.' : 'Reads are followed from the flow cell to “' + lower(last.label) + '”; the usable share is not assessed while the read data is partial.';
  const stages = ['classified', 'full_length', 'genome_mapped', 'quantified'];
  const samples = c.samples.filter(s => fs[s.key] && fs[s.key].classified != null), st = scopeTag(c, []);
  const endOf = f => { for (let i = stages.length - 1; i > 0; i--) if (f[stages[i]] != null) return stages[i]; return null; }, end = samples.length ? endOf(fs[samples[0].key]) : null;
  const names = end ? [lower(D.FUNNEL_LABELS[end]), 'given to the sample, not ' + lower(D.FUNNEL_LABELS[end])] : [];
  const srows = end ? samples.map(s => { const f = fs[s.key], q = f[end] != null ? f[end] : 0;
    return { label: s.name, short: c.short(s.key, 22), strong: c.sKeys.includes(s.key), parts: [{ value: q, color: tk('--fc-fate-usable') }, { value: Math.max(0, f.classified - q), color: tk('--fc-fate-not-full-length') }], f }; }) : [];
  return <Section id="read_fate" title="Read fate" result={sectionResult(b, 'read_fate')} reading={reading} x={{ evRule, evSec }}>
    <Partial p={p} />
    <Grid>
      <ChartPanel id="p-read-fate" title="Where every read went" reading={'Each bar starts from the same baseline: the ' + F.count(seq) + ' reads the flow cell produced.'} tag={tagOf(st, p)} height={rows.length * 28 + 10}
        spec={{ form: 'hbars', rowH: 28, rows, readout: i => fu[i].missing ? [fu[i].label + ': not known', fu[i].source || ''] : [[F.int(fu[i].n), 'reads'], [F.pct1(fu[i].n / seq), 'of the reads sequenced'], fu[i].source || ''] }}
        table={() => ({ head: ['Stage', 'Reads', 'Share of sequenced', 'Where the number comes from'], rows: fu.map(s => [s.label, s.missing ? 'not known' : F.int(s.n), s.missing ? F.DASH : F.pct1(s.n / seq), s.source || '']) })} />
      {srows.length ? <ChartPanel id="p-read-fate-samples" title="The same, per sample" reading={'Reads given to each sample, and how many of them were ' + names[0] + '.'} tag={tagOf(scopeTag(c, ['samples']), partialInfo(b, 'funnel.samples'))} height={srows.length * 24 + 30}
        spec={{ form: 'hstack', rows: srows, fmt: v => v >= 1e6 ? v / 1e6 + ' M' : F.int(v), readout: k => { const idx = [Math.floor(k / 1000), k % 1000], r = srows[idx[0]]; const tot = r.parts.reduce((a, q) => a + q.value, 0); return [[r.short || r.label, ''], [F.count(r.parts[idx[1]].value), names[idx[1]]], [F.pct1(r.parts[idx[1]].value / (tot || 1)), 'of the reads given to this sample']]; } }}
        legend={[{ name: D.FUNNEL_LABELS[end], color: tk('--fc-fate-usable') }, { name: 'Given to the sample, not ' + lower(D.FUNNEL_LABELS[end]), color: tk('--fc-fate-not-full-length') }]}
        table={() => ({ head: ['Sample'].concat(stages.filter(k => samples.some(s => fs[s.key][k] != null)).map(k => D.FUNNEL_LABELS[k])), rows: samples.map(s => [s.name].concat(stages.filter(k => samples.some(q => fs[q.key][k] != null)).map(k => fs[s.key][k] == null ? 'not known' : F.int(fs[s.key][k])))) })} />
        : <Empty title="Per sample: not known">No per-sample count of reads is linked to this run.</Empty>}
    </Grid>
    {missing.length && !p ? <p class="fc-section__note fc-notice">{missing.length === 1 ? 'One stage is drawn as an outline' : missing.length + ' stages are drawn as outlines'}, never as zero: {missing.map(s => <><b>{lower(s.label)}</b> ({(s.source || 'not known').replace(/^not known: /, '')})</>).reduce((a, x, i) => a.concat(i ? ['; '] : [], [x]), [])}.</p> : null}
    {fu.some(s => s.id === 'genome_mapped' && !s.missing) && fu.some(s => s.id === 'quantified' && !s.missing) ? <p class="fc-section__note">“Mapped to the genome” and “Counted towards a transcript” are two branches from the full-length reads, not one step after the other: a read can be counted towards a transcript without a genome alignment, so either number may be the larger.</p> : null}
    <Disclosure id="n-fate" c={c} title={'All numbers (' + countNumbers(b, ['funnel']) + ')'}><AllNumbers b={b} mods={['funnel']} /></Disclosure>
  </Section>;
}

// what stands under "Reads" and "Bases" in the samples table: the reads that were counted, in one word
const POP_SUB = { basecalled_pass: 'passed', basecalled_all: 'basecalled', classified: 'given to the sample', full_length_trimmed: 'full-length', aligned_primary: 'aligned' };
export function SamplesSection({ b, c, evRule, evSec }) {
  if (!c.named.length) return null;
  const fs = D.funnelSamples(b), cd = D.cdna(b), ix = index.value, sy = D.sampleYield(b), many = new Set(c.named.map(s => s.project).filter(Boolean)).size > 1, p = partialInfo(b, 'len') || partialInfo(b, 'yield.basecalled');
  const totalPass = Object.values(sy).reduce((t, o) => t + ((o.pass || {}).reads || 0), 0) || 1, withData = new Set(c.samples.map(s => s.key)), hasCdna = Object.keys(cd).length > 0, hasFate = Object.keys(fs).length > 0;
  const max = Math.max(1, ...c.named.map(s => ((sy[s.key] || {}).pass || {}).reads || 0));
  const pick = (e, k) => { if (e.target.closest('a')) return; const cur = c.sKeys; let next; if (e.ctrlKey || e.metaKey) next = cur.includes(k) ? cur.filter(x => x !== k) : cur.concat(k).slice(0, 8); else next = cur.length === 1 && cur[0] === k ? [] : [k]; setParams({ s: next.join(',') }); };
  const also = s => isReport || !ix ? [] : (ix.raw.samples || []).filter(o => o.sample_uid && o.sample_uid === s.sample_uid && o.run !== b.identity.run_key);
  const pop = (D.mod(b, 'len') || {}).population, un = sy.unclassified && sy.unclassified.pass, cols = 7 + (many ? 1 : 0) + (hasCdna ? 1 : 0) + (hasFate ? 1 : 0);
  const usableOf = f => f ? (f.quantified != null ? f.quantified : f.genome_mapped != null ? f.genome_mapped : f.full_length) : null;
  return <Section id="samples" title="Samples" reading={F.plural(c.named.length, 'sample') + (many ? ' of ' + new Set(c.named.map(s => s.project).filter(Boolean)).size + ' projects' : '') + ' on this flow cell' + (c.samples.length < c.named.length ? '; ' + c.samples.length + ' have read data so far.' : '.') + ' Counts are ' + readSetWords(pop, ['pass']) + '.'} x={{ evRule, evSec }}>
    <Partial p={p} />
    <div class="fc-tablewrap"><table class="fc-table fc-table--comfortable fc-samples"><caption class="fc-sr">Samples of this run</caption>
      <thead><tr><th scope="col">Sample</th><th scope="col">Barcode</th>{many ? <th scope="col">Project</th> : null}<th scope="col">Reads<span class="fc-table__sub">{POP_SUB[pop && pop.stage] || 'read'}</span></th><th scope="col" class="fc-num">Bases<span class="fc-table__sub">{POP_SUB[pop && pop.stage] || 'read'}</span></th>
        <th scope="col" class="fc-num">Read N50</th><th scope="col" class="fc-num">Median quality</th>{hasCdna ? <th scope="col" class="fc-num">Full-length<span class="fc-table__sub">of the reads classified</span></th> : null}{hasFate ? <th scope="col" class="fc-num">Usable reads</th> : null}<th scope="col">Also on</th></tr></thead>
      <tbody>{c.named.map(s => { const y = (sy[s.key] || {}).pass, has = withData.has(s.key) && y, H = has ? D.readHists(b, [s.key], ['pass']) : null, n50 = H && H.len ? D.n50(H.len.edges, H.len.w) : null, qm = H && H.qual ? D.quantile(H.qual.edges, H.qual.n, 0.5) : null;
        const on = c.sKeys.includes(s.key), others = also(s), cdn = cd[s.key], u = usableOf(fs[s.key]), can = withData.has(s.key), none = <span class="fc-muted">{can ? 'not read yet' : 'no files yet'}</span>;
        return <tr class={on ? 'is-marked' : ''} onClick={can ? e => pick(e, s.key) : undefined} aria-selected={on}>
          <th scope="row">{can ? <button type="button" class="fc-linkbtn fc-samples__name" aria-pressed={on} title={s.alias && s.alias !== s.name ? s.name + ' (written ' + s.alias + ' on this flow cell)' : s.name}><Mid text={s.name} tail={12} /></button> : <span class="fc-samples__name fc-samples__name--off" title={s.name}><Mid text={s.name} tail={12} /></span>}</th>
          <td>{s.barcode ? s.barcode.replace('barcode', '') : F.DASH}</td>{many ? <td>{s.project || F.DASH}</td> : null}
          <td class="fc-nowrap">{has ? <><span class="fc-cellbar" style={{ width: Math.max(1, Math.round(y.reads / max * 70)) + 'px' }} />{F.count(y.reads)} <span class="fc-muted">{F.pct1(y.reads / totalPass)}</span></> : none}</td>
          <td class="fc-num">{has ? F.bases(y.bases) : F.DASH}</td><td class="fc-num">{n50 ? F.bases(n50) : F.DASH}</td><td class="fc-num">{qm ? F.q(qm) : F.DASH}</td>
          {hasCdna ? <td class="fc-num">{cdn && cdn.pass_reads ? F.pct1(((cdn.primers_found || 0) + (cdn.rescued_reads || 0)) / cdn.pass_reads) : <span class="fc-muted">not known</span>}</td> : null}
          {hasFate ? <td class="fc-num">{u != null ? F.count(u) : <span class="fc-muted">not known</span>}</td> : null}
          <td>{others.length ? others.map(o => <a class="fc-mono" href={hrefRun(o.run)}>{o.run.split('_')[0]}</a>).reduce((a, x, i) => a.concat(i ? [', '] : [], [x]), []) : F.DASH}</td></tr>; })}
        <tr class={'is-dim' + (c.sKeys.includes('unclassified') ? ' is-marked' : '')} onClick={un ? e => pick(e, 'unclassified') : undefined}><th scope="row">{un ? <button type="button" class="fc-linkbtn" aria-pressed={c.sKeys.includes('unclassified')}>No barcode found</button> : 'No barcode found'}</th><td>{F.DASH}</td>{many ? <td>{F.DASH}</td> : null}
          <td>{un ? <>{F.count(un.reads)} <span class="fc-muted">{F.pct1(un.reads / totalPass)}</span></> : <span class="fc-muted">not known</span>}</td><td class="fc-num">{un ? F.bases(un.bases) : F.DASH}</td>
          {Array.from({ length: cols - (many ? 5 : 4) - 1 }, () => <td class="fc-num">{F.DASH}</td>)}<td>{F.DASH}</td></tr>
      </tbody></table></div>
    <SamplesChart b={b} c={c} />
    <p class="fc-section__note fc-screen">Choose a sample to show only that sample in every chart below; Ctrl-click adds further samples (up to 8): distributions are then drawn one over the other, counts are added up, and the scope line of each chart says which. The Samples control above does the same. {c.sKeys.length ? <button type="button" class="fc-linkbtn" onClick={() => setParams({ s: '' })}>Show all samples added up</button> : null}</p>
  </Section>;
}

// reads, bases and N50 kept by a custom filter, from the length x quality grid (cuts snap to the stored bin edges)
function customCut(grid, minq, minlen) {
  const ix = grid.xe.findIndex(e => e >= minlen - 1e-9), jy = grid.ye.findIndex(e => e >= minq - 1e-9), x0 = ix < 0 ? grid.nx : Math.min(ix, grid.nx), y0 = jy < 0 ? grid.ny : Math.min(jy, grid.ny);
  let all = 0, kept = 0, keptBases = 0; const len = new Float64Array(grid.nx), w = new Float64Array(grid.nx), qual = new Float64Array(grid.ny);
  for (let i = 0; i < grid.nx; i++) for (let j = 0; j < grid.ny; j++) { const n = grid.grid[i][j]; if (!n) continue; all += n; if (i >= x0 && j >= y0) { const bases = grid.w ? grid.w[i][j] : n * Math.sqrt(grid.xe[i] * grid.xe[i + 1]); kept += n; keptBases += bases; len[i] += n; w[i] += bases; qual[j] += n; } }
  return { all, kept, keptBases, len, w, qual, exactBases: !!grid.w, snapLen: grid.xe[Math.min(x0, grid.nx)], snapQ: grid.ye[Math.min(y0, grid.ny)] };
}

export function LengthQuality({ b, c, evRule, evSec }) {
  const mod = D.mod(b, 'len'); if (!mod) return null;
  const wait = waitingOnly(b, ['len', 'qual'], 'read length, read quality and their charts');
  if (wait) return <Section id="length_quality" title="Read length and quality" result={sectionResult(b, 'length_quality')} x={{ evRule, evSec }}><Waiting w={wait} /></Section>;
  const keys = c.keys, H = D.readHists(b, keys, c.filters);
  if (!H.len) return c.sKeys.length && (mod.cells || []).length ? <Section id="length_quality" title="Read length and quality" result={sectionResult(b, 'length_quality')} x={{ evRule, evSec }}><Partial p={partialInfo(b, 'len')} /><NoneFor c={c} what="reads" /></Section> : null;
  const grid = D.lenQualGrid(b, keys, c.filters), pop = mod.population || {}, set = readSetWords(pop, c.filters), p = partialInfo(b, 'len'), pq = partialInfo(b, 'qual') || p;
  const cut = c.custom && grid && (c.minq || c.minlen) ? customCut(grid, c.minq, c.minlen) : null, byB = c.by === 'bases', noGrid = c.custom && !grid;
  const n50 = cut ? D.n50(grid.xe, cut.w) : D.n50(H.len.edges, H.len.w), med = cut ? D.quantile(grid.xe, cut.len, 0.5, true) : D.quantile(H.len.edges, H.len.n, 0.5, true);
  const qmed = cut ? D.quantile(grid.ye, cut.qual, 0.5) : H.qual ? D.quantile(H.qual.edges, H.qual.n, 0.5) : null;
  // the two shares follow the custom filter like every other number of the section (its cuts snap to the bins of the grid)
  const total = D.sum(H.len.n), below = (e, n, edge) => { let t = 0; for (let i = 0; i < n.length; i++) if (e[i + 1] <= edge) t += n[i]; return t / (D.sum(n) || 1); }, from = (e, n, edge) => { let t = 0; for (let i = 0; i < n.length; i++) if (e[i] >= edge) t += n[i]; return t / (D.sum(n) || 1); };
  const short = cut ? below(grid.xe, cut.len, 200) : below(H.len.edges, H.len.n, 200), q20 = cut ? from(grid.ye, cut.qual, 20) : H.qual ? from(H.qual.edges, H.qual.n, 20) : null;
  const overlay = c.sKeys.length > 1, ycfg = title => c.logY ? { log: true, title: title + ', logarithmic' } : { title };
  const unit = v => byB ? (v >= 1e9 ? 1e9 : v >= 1e6 ? 1e6 : 1) : (v >= 1e6 ? 1e6 : 1), unitWord = (u, what) => what + ' per bin' + (u === 1e9 ? ', Gb' : u === 1e6 ? (byB ? ', Mb' : ', millions') : '');
  const info = k => ({ s: c.named.find(q => q.key === k), i: c.named.findIndex(q => q.key === k), name: sampleName(c, k) });
  // ---- length ----
  let lenSpec, lenLegend = null, lenTable;
  const st = scopeTag(c, ['samples', 'by', 'y', 'reads', 'custom', 'overlay']);
  if (overlay) { // each sample as an outline scaled to its own total
    const ser = c.sKeys.map(k => { const { s, i, name } = info(k), h = D.readHists(b, [k], c.filters), a = byB ? h.len.w : h.len.n, t = D.sum(a) || 1; return { name, type: 'outline', edges: h.len.edges, y: Array.from(a, v => v / t * 100), color: sampleColor(s, i < 0 ? 11 : i) }; });
    lenSpec = { form: 'xy', x: { domain: D.spanQ(H.len.edges, H.len.n), log: true, kind: 'length', title: 'Read length, logarithmic' }, y: ycfg('Share of each sample’s ' + (byB ? 'bases' : 'reads') + ' per bin, %'), series: ser, cursor: { x: centersLog(H.len.edges) },
      readout: i => [[F.int(H.len.edges[i]) + ' to ' + F.int(H.len.edges[i + 1]) + ' bases', '']].concat(ser.map(q => [q.y[i].toFixed(2) + ' %', c.shorten(q.name, 18)])) };
    lenLegend = ser.map(q => ({ name: q.name, label: q.name, kind: 'line', color: q.color }));
    lenTable = () => ({ head: ['Read length, bases'].concat(ser.map(q => q.name + ', %')), rows: Array.from(H.len.n, (_, i) => [F.int(H.len.edges[i]) + ' to ' + F.int(H.len.edges[i + 1])].concat(ser.map(q => q.y[i].toFixed(2)))).filter((_, i) => H.len.n[i] > 0) });
  } else {
    const e = cut ? grid.xe : H.len.edges, nn = cut ? cut.len : H.len.n, ww = cut ? cut.w : H.len.w, a = byB ? ww : nn, tot = D.sum(a) || 1, u = unit(Math.max(...a));
    lenSpec = { form: 'xy', x: { domain: D.spanQ(H.len.edges, H.len.n), log: true, kind: 'length', title: 'Read length, logarithmic' }, y: ycfg(unitWord(u, byB ? 'Bases' : 'Reads')), series: [{ name: cut ? 'Kept by the custom filter' : 'Reads', type: 'hist', edges: e, y: Array.from(a, v => v / u) }],
      vlines: [n50 ? { x: n50, label: 'N50 ' + F.bases(n50), strong: true } : null, !cut && med ? { x: med, label: 'median ' + F.bases(med), side: 'left' } : null].filter(Boolean), cursor: { x: centersLog(e) },
      readout: i => [[F.int(e[i]) + ' to ' + F.int(e[i + 1]) + ' bases', ''], [F.int(nn[i]), 'reads'], [F.bases(ww[i]), ''], [F.pct1(a[i] / tot), 'of the ' + (byB ? 'bases' : 'reads')]] };
    lenTable = () => ({ head: ['Read length, bases', 'Reads', 'Bases'], rows: Array.from(nn, (v, i) => [F.int(e[i]) + ' to ' + F.int(e[i + 1]), F.int(v), F.int(ww[i])]).filter((_, i) => nn[i] > 0) });
  }
  // ---- quality ----
  let qSpec = null, qLegend = null, qTable = null;
  if (H.qual) {
    const qdom = D.spanQ(H.qual.edges, H.qual.n, 0.0002, 4);
    if (overlay) { const ser = c.sKeys.map(k => { const { s, i, name } = info(k), h = D.readHists(b, [k], c.filters), t = D.sum(h.qual.n) || 1; return { name, type: 'outline', edges: h.qual.edges, y: Array.from(h.qual.n, v => v / t * 100), color: sampleColor(s, i < 0 ? 11 : i) }; });
      qSpec = { form: 'xy', x: { domain: qdom, kind: 'linear', title: 'Mean read quality, Q', fmt: v => 'Q ' + v }, y: ycfg('Share of each sample’s reads per bin, %'), series: ser, cursor: { x: centers(H.qual.edges) }, readout: i => [['Q ' + H.qual.edges[i].toFixed(2) + ' to ' + H.qual.edges[i + 1].toFixed(2), '']].concat(ser.map(q => [q.y[i].toFixed(2) + ' %', c.shorten(q.name, 18)])) };
      qLegend = ser.map(q => ({ name: q.name, label: q.name, kind: 'line', color: q.color }));
      qTable = () => ({ head: ['Mean read quality'].concat(ser.map(q => q.name + ', %')), rows: Array.from(H.qual.n, (_, i) => ['Q ' + H.qual.edges[i].toFixed(2) + ' to ' + H.qual.edges[i + 1].toFixed(2)].concat(ser.map(q => q.y[i].toFixed(2)))).filter((_, i) => H.qual.n[i] > 0) });
    } else { const e = cut ? grid.ye : H.qual.edges, a = cut ? cut.qual : H.qual.n, u = unit(Math.max(...a)) === 1e9 ? 1e6 : (Math.max(...a) >= 1e6 ? 1e6 : 1);
      qSpec = { form: 'xy', x: { domain: qdom, kind: 'linear', title: 'Mean read quality, Q', fmt: v => 'Q ' + v }, y: ycfg('Reads per bin' + (u === 1e6 ? ', millions' : '')), series: [{ name: 'Reads', type: 'hist', edges: e, y: Array.from(a, v => v / u) }],
        vlines: [qmed ? { x: qmed, label: 'median ' + F.q(qmed), strong: true } : null, pop.min_qscore ? { x: pop.min_qscore, label: 'pass: Q ≥ ' + pop.min_qscore, side: 'left' } : null].filter(Boolean), cursor: { x: centers(e) }, readout: i => [['Q ' + e[i] + ' to ' + e[i + 1], ''], [F.int(a[i]), 'reads']] };
      qTable = () => ({ head: ['Mean read quality', 'Reads'], rows: Array.from(a, (v, i) => ['Q ' + e[i] + ' to ' + e[i + 1], F.int(v)]).filter((_, i) => a[i] > 0) }); }
  }
  // ---- length x quality ----
  let heat = null, heatTable = null;
  if (grid) { const I = [], J = [], T = []; let mx = 0; for (let i = 0; i < grid.nx; i++) for (let j = 0; j < grid.ny; j++) if (grid.grid[i][j] > mx) mx = grid.grid[i][j];
    for (let i = 0; i < grid.nx; i++) for (let j = 0; j < grid.ny; j++) { const n = grid.grid[i][j]; if (n > 0) { I.push(i); J.push(j); T.push(Math.sqrt(n / mx)); } }
    const xn = new Float64Array(grid.nx); for (let i = 0; i < grid.nx; i++) xn[i] = D.sum(grid.grid[i]);
    const yn = new Float64Array(grid.ny); for (let i = 0; i < grid.nx; i++) for (let j = 0; j < grid.ny; j++) yn[j] += grid.grid[i][j];
    const colMed = Array.from({ length: grid.nx }, (_, i) => xn[i] > 0 ? D.quantile(grid.ye, grid.grid[i], 0.5) : null);
    heat = { form: 'xy', x: { domain: D.span(grid.xe, xn, 1), log: true, kind: 'length', title: 'Read length, logarithmic' }, y: { domain: D.span(grid.ye, yn, 1), title: 'Mean read quality, Q', fmt: v => String(v) }, heat: { xe: grid.xe, ye: grid.ye, i: I, j: J, t: T }, series: [],
      cursor: { x: centersLog(grid.xe) }, readout: i => [[F.int(grid.xe[i]) + ' to ' + F.int(grid.xe[i + 1]) + ' bases', ''], [F.int(xn[i]), 'reads of this length']].concat(colMed[i] != null ? [[F.q(colMed[i]), 'their median quality']] : []),
      vlines: cut && c.minlen ? [{ x: cut.snapLen, label: '≥ ' + F.int(cut.snapLen) + ' b', strong: true }] : [], limits: cut && c.minq ? [{ y: cut.snapQ, kind: 'target', label: 'Q ≥ ' + cut.snapQ }] : [] };
    heatTable = () => ({ head: ['Read length, bases', 'Mean read quality', 'Reads'].concat(grid.w ? ['Bases'] : []), rows: I.map((i, k) => [F.int(grid.xe[i]) + ' to ' + F.int(grid.xe[i + 1]), 'Q ' + grid.ye[J[k]] + ' to ' + grid.ye[J[k] + 1], F.int(grid.grid[i][J[k]])].concat(grid.w ? [F.int(grid.w[i][J[k]])] : [])) }); }
  const reading = 'Half of the bases are in reads of ' + F.bases(n50) + ' or longer' + (qmed ? '; typical read quality is ' + F.q(qmed) : '') + '. These are ' + set + (c.sKeys.length ? ', of ' + chosenWords(c) + (overlay ? ' added up' : '') : '') + (cut ? ', cut by the custom filter' : '') + '.';
  const cap = text => p ? text + ' · ' + p.word : text, reads = cut ? cut.kept : H.reads || total, basesT = cut ? cut.keptBases : H.bases || D.sum(H.len.w);
  return <Section id="length_quality" title="Read length and quality" result={sectionResult(b, 'length_quality')} reading={reading} x={{ evRule, evSec }}>
    <Partial p={p} />
    <Variant text={variantNote(b, ['len', 'qual'])} />
    {cut ? <div class="fc-evidence"><span>Custom filter{c.minq ? ', Q ≥ ' + cut.snapQ : ''}{c.minlen ? ', length ≥ ' + F.int(cut.snapLen) + ' bases' : ''} (cuts snap to the stored bin edges): <b>{F.int(cut.kept)}</b> of {F.int(cut.all)} reads kept ({F.pct1(cut.kept / (cut.all || 1))}). The numbers below differ from the ones the verdict used.</span><button type="button" class="fc-btn" onClick={() => setParams({ reads: '', minq: '', minlen: '' })}>Remove the filter</button></div>
      : noGrid ? <div class="fc-evidence"><span>The Custom filter needs the length against quality grid, which is stored per sample only for runs of up to 12 samples: the numbers below are not filtered.</span><button type="button" class="fc-btn" onClick={() => setParams({ reads: '', minq: '', minlen: '' })}>Remove the filter</button></div> : null}
    <div class="fc-tiles fc-section__tiles">
      <Tile label="Read N50" value={F.basesParts(n50)[0]} unit={F.basesParts(n50)[1]} caption={cap(cut ? 'under the custom filter' : set)} title="len_n50" />
      <Tile label="Median read length" value={F.basesParts(med)[0]} unit={F.basesParts(med)[1]} caption={cap(cut ? 'under the custom filter' : set)} title="len_median" />
      <Tile label="Reads" value={F.countParts(reads)[0]} unit={F.countParts(reads)[1]} caption={cap(F.bases(basesT) + (keys.length < c.keys.length || c.sKeys.length ? ' · ' + F.plural(c.sKeys.length, 'sample') : ''))} />
      {qmed ? <Tile label="Median read quality" value={F.q(qmed)} caption={pq ? set + ' · ' + pq.word : set} title="qual_median" /> : <Tile label="Median read quality" state="unknown" value="not known" caption="the files carry no quality per read" />}
      {q20 != null ? <Tile label="Reads at Q 20 or better" value={F.pctNum(q20)} unit="%" caption={(cut ? 'of the reads kept by the custom filter' : 'of the reads') + (pq ? ' · ' + pq.word : '')} /> : null}
      <Tile label="Reads under 200 bases" value={F.pctNum(short)} unit="%" caption={cap(cut ? 'of the reads kept by the custom filter' : 'of the reads')} />
    </div>
    <Grid wide={heat ? 3 : 2}>
      <ChartPanel id="p-len" title="Read length" reading={overlay ? 'Each sample is drawn as an outline scaled to its own total, so shapes compare.' : 'Half of the bases are in reads of ' + F.bases(n50) + ' or longer.'} scope={st.scope} tag={tagOf(st, p)} height={220} spec={lenSpec} legend={lenLegend} legendToggle={overlay} table={lenTable} />
      {qSpec ? <ChartPanel id="p-qual" title="Read quality" reading={'Typical read quality is ' + F.q(qmed) + '.'} scope={scopeTag(c, ['samples', 'y', 'reads', 'custom', 'overlay']).scope} tag={tagOf(scopeTag(c, ['samples', 'y', 'reads', 'custom', 'overlay']), pq)} height={220} spec={qSpec} legend={qLegend} legendToggle={overlay} table={qTable} /> : null}
      {heat ? <ChartPanel id="p-lenqual" title="Length against quality" reading={'Where the reads sit; the Custom filter of the Reads control cuts this grid.' + (grid.outside ? ' ' + F.int(grid.outside) + ' reads fall outside its axes and are not drawn.' : '')} scope={scopeTag(c, ['samples', 'reads', 'custom']).scope} tag={tagOf(scopeTag(c, ['samples', 'reads', 'custom']), pq)} height={220} spec={heat} table={heatTable} footer={<div class="fc-chart__foot"><ScaleLegend lo="few reads" hi="most reads" none={null} /></div>} /> : null}
    </Grid>
    <Disclosure id="n-len" c={c} title={'All numbers (' + countNumbers(b, ['yield.basecalled', 'len', 'qual']) + ')'}><AllNumbers b={b} mods={['yield.basecalled', 'len', 'qual']} /></Disclosure>
  </Section>;
}

export function Alignment({ b, c, evRule, evSec }) {
  const mod = D.mod(b, 'aln'); if (!mod) return null;
  const H = D.readHists(b, c.keys, c.filters), pop = mod.population || {}, r = ruleById(b, 'aln_mapped_frac'), p = partialInfo(b, 'aln'), set = readSetWords(pop, ['pass']);
  // the alignment files wait to be read: say that, and nothing that would read as a measurement
  const wait = waitingOnly(b, ['aln'], 'the alignment numbers and charts');
  if (wait) return <Section id="alignment" title="Alignment" result={sectionResult(b, 'alignment')} x={{ evRule, evSec }}><Waiting w={wait} /></Section>;
  if (!H.ident || !H.mapq || !H.alnReads) return <Section id="alignment" title="Alignment" result={sectionResult(b, 'alignment')} x={{ evRule, evSec }}><Partial p={p} />
    {c.sKeys.length ? <NoneFor c={c} what="alignment numbers" /> : <Empty title="No alignment numbers">{(mod.cells || []).length ? 'The alignment files that were read hold no read of the read set chosen in the Reads control.' : 'No alignment file of this run has been read: ' + (((mod.notes || [])[0] || {}).message || 'the files that are linked hold no read of this run.')}</Empty>}</Section>;
  const idTotal = D.sum(H.ident.n) + H.ident.over, idQ = D.sum(H.ident.n) > 0 ? D.quantile(H.ident.edges, H.ident.n, Math.min(1, 0.5 * idTotal / D.sum(H.ident.n))) : null;
  const perfectHalf = H.ident.over >= idTotal / 2, ident = perfectHalf ? 1 : idQ != null ? 1 - Math.pow(10, -idQ / 10) : null, identText = ident == null ? 'not known' : perfectHalf ? '100 %' : (ident * 100).toFixed(1) + ' %';
  const q2p = q => (100 * (1 - Math.pow(10, -q / 10))), idom = D.spanQ(H.ident.edges, H.ident.n, 0.001, 3), mq = H.mapq, top = mq.n.length - 1, mqTop = mq.n[top] / (D.sum(mq.n) || 1), mapped = H.mapped / (H.alnReads || 1);
  const ticks = [10, 13, 20, 23, 30, 40].filter(v => v >= idom[0] && v <= idom[1]).map(v => ({ v, label: +q2p(v).toFixed(v >= 30 ? 2 : 1) + ' %', major: true }));
  const u = Math.max(...H.ident.n) >= 1e6 ? 1e6 : 1, um = Math.max(...mq.n) >= 1e6 ? 1e6 : 1, per = k => k === 1e6 ? ', millions' : '';
  const iSpec = { form: 'xy', x: { domain: idom, kind: 'linear', ticks, title: 'Identity with the genome, spaced by error rate' }, y: c.logY ? { log: true, title: 'Reads per bin' + per(u) + ', logarithmic' } : { title: 'Reads per bin' + per(u) }, series: [{ name: 'Reads', type: 'hist', edges: H.ident.edges, y: Array.from(H.ident.n, v => v / u) }],
    vlines: idQ != null && !perfectHalf ? [{ x: idQ, label: 'median ' + identText, strong: true }] : [], cursor: { x: centers(H.ident.edges) }, readout: i => [[q2p(H.ident.edges[i]).toFixed(2) + ' to ' + q2p(H.ident.edges[i + 1]).toFixed(2) + ' % identity', ''], [F.int(H.ident.n[i]), 'reads']] };
  const mSpec = { form: 'xy', x: { domain: [0, mq.n.length], kind: 'linear', title: 'Mapping quality', fmt: v => String(v) }, y: c.logY ? { log: true, title: 'Reads per value' + per(um) + ', logarithmic' } : { title: 'Reads per value' + per(um) }, series: [{ name: 'Reads', type: 'hist', edges: mq.edges, y: Array.from(mq.n, v => v / um) }], cursor: { x: centers(mq.edges) }, readout: i => [['Mapping quality ' + i, ''], [F.int(mq.n[i]), 'reads']] };
  const st0 = scopeTag(c, ['samples', 'y', 'reads']), st = { scope: [set, st0.scope].filter(Boolean).join(' · '), tag: st0.tag }, cap = text => p ? text + ' · ' + p.word : text;
  return <Section id="alignment" title="Alignment" result={sectionResult(b, 'alignment')} reading={F.pct1(mapped) + ' of ' + ofReads(pop) + (c.sKeys.length ? ' of ' + chosenWords(c) : '') + ' map to the genome; typical identity is ' + identText + '.' + (r && r.outcome === 'not_evaluated' ? ' ' + r.message : '')} x={{ evRule, evSec }}>
    <Partial p={p} />
    <Variant text={variantNote(b, ['aln'])} />
    <div class="fc-tiles fc-section__tiles">
      <Tile label="Mapped reads" value={F.pctNum(mapped)} unit="%" caption={cap(F.count(H.mapped) + ' of ' + F.count(H.alnReads) + ' ' + set)} title="aln_mapped_frac" />
      <Tile label="Median identity" value={ident != null ? (perfectHalf ? '100' : (ident * 100).toFixed(1)) : F.DASH} unit="%" caption={cap(set)} title="aln_identity_median" />
      <Tile label="Reads without a single difference" value={F.pctNum(H.perfect / (H.mapped || 1))} unit="%" caption={cap('of mapped reads')} />
      <Tile label={'Mapping quality ' + top} value={F.pctNum(mqTop)} unit="%" caption={cap('of mapped reads: placed without doubt')} />
      <Tile label="Extra alignment records" value={F.count(H.secondary + H.supplementary)} caption={cap(F.count(H.secondary) + ' secondary · ' + F.count(H.supplementary) + ' supplementary; counted, not read as reads')} />
      {H.columns ? <Tile label="Differences per 100 aligned bases" value={(H.edits / H.columns * 100).toFixed(2)} caption={cap('mismatches, insertions and deletions')} /> : null}
    </div>
    <AlignmentOpening b={b} c={c} p={p} more={<>
      <ChartPanel id="p-ident" title="Identity with the genome" reading={'Half of the mapped reads match the genome at ' + identText + ' or better; ' + F.pct(H.ident.over / (idTotal || 1)) + ' have no difference at all and are not drawn.'} scope={st.scope} tag={tagOf(st, p)} height={220} spec={iSpec}
        table={() => ({ head: ['Identity, %', 'Reads'], rows: Array.from(H.ident.n, (v, i) => [q2p(H.ident.edges[i]).toFixed(2) + ' to ' + q2p(H.ident.edges[i + 1]).toFixed(2), F.int(v)]).filter((_, i) => H.ident.n[i] > 0).concat([['100 (no difference)', F.int(H.ident.over)]]) })} />
      <Calibration b={b} c={c} p={p} /></>} />
    <Disclosure id="aln-more" c={c} charts title="More charts: mapping quality, soft clipping, accuracy by read length, the length and identity grids">
      <AlignmentMore b={b} c={c} p={p} first={<ChartPanel id="p-mapq" title="Mapping quality" reading={F.pct(mqTop) + ' of mapped reads have mapping quality ' + top + ': the aligner had no doubt where they belong.'} scope={st.scope} tag={tagOf(st, p)} height={240} spec={mSpec}
        table={() => ({ head: ['Mapping quality', 'Reads'], rows: Array.from(mq.n, (v, i) => [String(i), F.int(v)]).filter((_, i) => mq.n[i] > 0) })} />} /></Disclosure>
    <Disclosure id="n-aln" c={c} title={'All numbers (' + countNumbers(b, ['aln']) + ')'}><AllNumbers b={b} mods={['aln']} /></Disclosure>
  </Section>;
}

// The classes of the cDNA step from the section "cdna". A bundle that holds the whole pychopper table (cdna.pychopper)
// shows them in the Pychopper section, with the check; this section is drawn only for a bundle without that table
// (one made before the table was read whole). Names and denominators are those of the Pychopper section.
const CDNA = [['primers_found', 'Full-length, one cDNA', '--fc-cdna-full-length'], ['rescued_reads', 'Full-length, rescued', '--fc-cdna-rescued'], ['unusable', 'Unusable', '--fc-cdna-unusable'], ['qc_fail', 'Dropped: low quality', '--fc-cdna-quality-fail']];
export function CdnaLibrary({ b, c, evRule, evSec }) {
  const cd = D.cdna(b), sm = D.summary(b, 'cdna'); if (!D.mod(b, 'cdna') || D.mod(b, 'cdna.pychopper')) return null;
  const samples = c.named.filter(s => cd[s.key]), r = ruleById(b, 'cdna_full_length_frac'), p = partialInfo(b, 'cdna');
  if (!samples.length) return null;
  const tot = k => samples.reduce((a, s) => a + (cd[s.key][k] || 0), 0), input = tot('input') || 1, classified = tot('pass_reads') || 1, pieces = tot('primers_found') + tot('rescue_segments');
  const shared = sharedTables(b, ['cdna'], c, 'pychopper'), capS = t => shared ? t + ' · ' + shared.caption : t;
  const rows = samples.map(s => ({ label: s.name + (shared ? shared.mark(s.key) : ''), short: c.short(s.key, 22) + (shared ? shared.mark(s.key) : ''), strong: c.sKeys.includes(s.key), parts: CDNA.map(([k, , token]) => ({ value: cd[s.key][k] || 0, color: tk(token) })) }));
  const full = sm.cdna_full_length_frac != null ? sm.cdna_full_length_frac : (tot('primers_found') + tot('rescued_reads')) / classified;
  return <Section id="cdna" title="cDNA library" result={sectionResult(b, 'cdna')} reading={F.pct1(full) + ' of the ' + F.count(classified) + ' reads pychopper classified' + (shared ? ', merged from ' + shared.both + ' of the sample,' : '') + ' are full-length cDNA (one cDNA with both primers, or rescued from a read that held several). ' + (r && r.outcome === 'not_evaluated' ? r.message : '')} x={{ evRule, evSec }}>
    {shared ? <p class="fc-section__note fc-notice fc-shared" role="status">{shared.text}</p> : null}
    <Partial p={p} />
    <div class="fc-tiles fc-section__tiles">
      <Tile label="Full-length reads" value={F.pct1(full).replace(/\s*%$/, '')} unit="%" caption={capS(F.count(tot('primers_found') + tot('rescued_reads')) + ' of ' + F.count(classified) + ' reads classified, with the rescued ones')} title="cdna_full_length_frac" />
      <Tile label="Full-length, rescued" value={F.pctNum(tot('rescued_reads') / classified)} unit="%" caption={capS(F.int(tot('rescued_reads')) + ' of ' + F.count(classified) + ' reads classified, cut into ' + F.int(tot('rescue_segments')) + ' pieces')} />
      <Tile label="Unusable" value={F.pctNum(tot('unusable') / classified)} unit="%" caption={capS(F.count(tot('unusable')) + ' of ' + F.count(classified) + ' reads classified')} />
      <Tile label="Dropped: low quality" value={F.pctNum(tot('qc_fail') / input)} unit="%" caption={capS(F.int(tot('qc_fail')) + ' of ' + F.count(input) + ' reads given, before classifying')} />
      {samples.every(s => cd[s.key].umi_detected != null) && pieces ? <Tile label="With a UMI" value={F.pctNum(tot('umi_detected') / pieces)} unit="%" caption={capS(F.count(tot('umi_detected')) + ' of ' + F.count(pieces) + ' pieces kept')} /> : null}
      <Tile label="Samples" value={String(samples.length)} caption="one pychopper statistics table each" />
    </div>
    <ChartPanel id="p-cdna" title="What pychopper made of the reads" reading="Each bar is every read of one sample that was given to pychopper." tag={tagOf(scopeTag(c, ['samples']), p)} height={rows.length * 24 + 30}
      spec={{ form: 'hstack', norm: true, rows, readout: k => { const ri = Math.floor(k / 1000), pi = k % 1000, row = rows[ri], t = row.parts.reduce((a, q) => a + q.value, 0); return [[row.short || row.label, ''], [F.pct1(row.parts[pi].value / (t || 1)), F.lower(CDNA[pi][1])], [F.int(row.parts[pi].value), 'reads'], 'of ' + F.count(t) + ' reads given to pychopper']; } }}
      legend={CDNA.map(([, name, token]) => ({ name, color: tk(token) }))}
      table={() => ({ head: ['Sample'].concat(CDNA.map(x => x[1])), rows: samples.map(s => [s.name].concat(CDNA.map(([k]) => F.int(cd[s.key][k])))) })} />
    <Disclosure id="n-cdna" c={c} title={'All numbers (' + countNumbers(b, ['cdna']) + ')'}><AllNumbers b={b} mods={['cdna']} /></Disclosure>
  </Section>;
}
