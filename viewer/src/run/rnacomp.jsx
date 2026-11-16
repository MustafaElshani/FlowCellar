// RNA composition: what kinds of transcript the reads come from, how long their poly(A) tails are by transcript type,
// and how the reads spread over genes (sections rna.biotype and rna.genes of a bundle). The numbers are counted from
// each sample's annotated read table: one row per read that the pipeline gave to a transcript.
import { Tile } from '../shell/parts.jsx';
import { ChartPanel } from '../charts/ChartPanel.jsx';
import { tk } from '../charts/theme.js';
import * as F from '../core/format.js';
import * as D from '../core/derive.js';
import { Section, Disclosure, AllNumbers, countNumbers } from './section.jsx';
import { scopeTag, sharedTables } from './common.js';
import { sampleColor } from './readdata.jsx';
import { Grid } from '../charts/Grid.jsx';

export const BIOTYPE = 'rna.biotype', GENES = 'rna.genes';
const TAIL = ['reads', 'tail_estimated', 'tail_zero', 'tail_not_found', 'tail_no_value', 'tail_sum'];
const share = (a, b) => F.isNum(a) && b > 0 ? a / b : null;
// a share in words: one decimal, two for a share under one in a thousand (so that a few reads never read as "0.0 %")
const pc = f => !F.isNum(f) ? F.DASH : f > 0 && f < 0.00005 ? 'under 0.01 %' : f > 0 && f < 0.001 ? (f * 100).toFixed(2) + ' %' : F.pct1(f);
const pcNum = f => !F.isNum(f) ? F.DASH : f > 0 && f < 0.001 ? (f * 100).toFixed(2) : F.pctNum(f);
const word = t => t === 'other' ? 'Other types' : t.replace(/_/g, ' ');
const classColor = (t, i, line) => t === 'other' ? tk(line ? '--fc-chart-none-line' : '--fc-chart-none') : tk('--fc-cat-' + (i % 12 + 1));
const NO_TAIL = [['With a tail estimate', '--fc-fate-usable', 'tail_estimated'], ['Primer found, no tail length', '--fc-cat-3', 'tail_zero'], ['Primer not found', '--fc-cat-5', 'tail_not_found'], ['No value in the table', '--fc-chart-none', 'tail_no_value']];
const nt = v => v == null ? F.DASH : String(Math.round(v));

// cells of the chosen samples added up, class by class: the numbers, the tail histogram, the reads of every biotype
function addUp(b, cells, classes) {
  const out = {}, types = new Map();
  for (const cls of classes.concat(['all'])) out[cls] = { sc: Object.fromEntries(TAIL.map(k => [k, 0])), n: null, edges: null, over: 0 };
  for (const cell of cells) {
    const cls = cell.at[1], h = cell.hists && cell.hists.rna_biotype_polya_hist, t = cell.counts && cell.counts.rna_biotype_table;
    for (const o of [out[cls], out.all]) { if (!o) continue; for (const k of TAIL) o.sc[k] += cell.scalars[k] || 0;
      if (h) { const edges = D.edgesOf(b, h.axis), n = D.expand(h, edges.length - 1); if (!o.n) { o.n = n; o.edges = edges; } else D.addInto(o.n, n); o.over += h.over || 0; } }
    if (t) t.cats.forEach((name, i) => { const e = types.get(name) || { name, cls, reads: 0, tail: 0 }; e.reads += t.n[i]; e.tail += t.w ? t.w[i] : 0; types.set(name, e); });
    if (t && t.other) { const e = types.get('(rarer types)') || { name: '(rarer types)', cls, reads: 0, tail: 0 }; e.reads += t.other; types.set('(rarer types)', e); }
  }
  return { by: out, types: [...types.values()].sort((a, x) => x.reads - a.reads) };
}
const histMedian = o => o && o.n ? D.quantile(o.edges, o.n, Math.min(1, 0.5 * (D.sum(o.n) + o.over) / (D.sum(o.n) || 1))) : null;
const histQ = (o, q) => o && o.n ? D.quantile(o.edges, o.n, Math.min(1, q * (D.sum(o.n) + o.over) / (D.sum(o.n) || 1))) : null;

export function RnaComposition({ b, c, evRule, evSec }) {
  const m = D.mod(b, BIOTYPE), g = D.mod(b, GENES); if (!m) return null;
  const waitingNote = (m.notes || []).find(n => n.code === 'waiting'), st = scopeTag(c, ['samples']);
  const have = new Set((m.cells || []).map(x => x.at[0])), all = c.named.filter(s => have.has(s.key));
  const partial = m.state === 'partial' ? (m.partial_reason === 'waiting' ? (waitingNote ? waitingNote.message : 'Some annotated read tables of this run have not been read yet.') + ' The numbers here will change when the next scan has read them.'
    : m.partial_reason === 'partly_read' ? 'A table was not read to its end; counts are a lower bound.' : all.length + ' of ' + c.named.length + ' samples have an annotated read table linked to this run; the numbers here cover those samples only.') : null;
  if (!all.length) return <Section id="rna_composition" title="RNA composition" x={{ evRule, evSec }}>
    <p class="fc-section__note fc-notice fc-partial" role="status">{partial ? 'Partial: ' + partial : 'The annotated read tables of this run hold no read.'}</p></Section>;
  const classes = (m.params && m.params.classes) || [...new Set(m.cells.map(x => x.at[1]))];
  const chosen = all.filter(s => c.sKeys.includes(s.key)), sel = chosen.length ? chosen : all, keys = new Set(sel.map(s => s.key)), one = sel.length === 1;
  // a table that covers the sample on two flow cells: said at the top, in the sentence, under every count and on every row
  const shared = sharedTables(b, [BIOTYPE, GENES], c, 'table'), sharedSel = !!shared && sel.some(q => shared.keys.has(q.key)), capS = t => sharedSel ? t + ' · ' + shared.caption : t;
  const mark = k => shared ? shared.mark(k) : '', mid = (name, n = 26) => c.shorten(name, n), noteS = shared ? shared.legend : undefined;
  const A = addUp(b, m.cells.filter(x => keys.has(x.at[0])), classes), tot = A.by.all.sc, reads = tot.reads || 1;
  const perSample = {}; for (const s of all) perSample[s.key] = addUp(b, m.cells.filter(x => x.at[0] === s.key), classes);
  const tag = [st.tag, m.state === 'partial' ? 'PARTIAL · ' + (m.partial_reason === 'waiting' ? 'files waiting' : m.partial_reason === 'partly_read' ? 'read in part' : all.length + ' of ' + c.named.length + ' samples') : null].filter(Boolean).join(' · ');
  const strong = k => chosen.length > 0 && chosen.length < all.length && c.sKeys.includes(k);
  const sumLabel = chosen.length ? sel.length + ' chosen' : 'All samples', height = rows => rows.length * 22 + 34;
  const hasTail = tot.tail_estimated > 0, typeShare = t => share((A.types.find(x => x.name === t) || {}).reads || 0, reads);
  const rrna = ((m.params && m.params.rrna_biotypes) || ['rRNA', 'rRNA_pseudogene', 'Mt_rRNA']).reduce((a, t) => a + ((A.types.find(x => x.name === t) || {}).reads || 0), 0) / reads;
  // exact quartiles are stored per sample and for all samples together; a choice of several samples reads them from the 5 nt histogram
  const stats = {}; for (const r of (((m.run_level || {}).tables || {}).rna_polya_by_biotype || { rows: [] }).rows) stats[r[0] + '|' + r[1]] = r;
  const exactKey = one ? sel[0].key : sel.length === all.length ? '*' : null;
  const tailOf = cls => { const r = exactKey ? stats[exactKey + '|' + cls] : null, o = A.by[cls];
    return r ? { median: r[4], q25: r[5], q75: r[6], mean: r[7], exact: true } : { median: histMedian(o), q25: histQ(o, 0.25), q75: histQ(o, 0.75), mean: share(o.sc.tail_sum, o.sc.tail_estimated), exact: false }; };
  const overall = tailOf('all');

  // ---- 1. reads by transcript type, per sample ----
  const typeRows = (sel.length > 1 ? [{ label: sumLabel, strong: true, parts: classes.map((t, i) => ({ value: A.by[t].sc.reads, color: classColor(t, i) })) }] : [])
    .concat(all.map(s => ({ label: s.name + mark(s.key), short: c.short(s.key, 22) + mark(s.key), strong: strong(s.key), parts: classes.map((t, i) => ({ value: perSample[s.key].by[t].sc.reads, color: classColor(t, i) })) })));
  const top = A.types[0], second = A.types[1];
  // ---- 2. every transcript type, ranked ----
  const listed = A.types.slice(0, 14), restReads = A.types.slice(14).reduce((a, t) => a + t.reads, 0);
  const rankRows = listed.map(t => ({ label: word(t.name), value: t.reads, text: F.count(t.reads), note: pc(t.reads / reads), color: classColor(t.cls, classes.indexOf(t.cls)) }))
    .concat(restReads ? [{ label: F.plural(A.types.length - 14, 'rarer type'), value: restReads, text: F.count(restReads), note: pc(restReads / reads), color: tk('--fc-chart-none') }] : []);
  // ---- 3. poly(A) tail length by transcript type ----
  let tailSpec = null, tailLegend = null, tailTable = null, tailReading = '';
  if (hasTail) {
    const base = A.by.all, edges = base.edges, dom = D.spanQ(edges, base.n, 0.005, 2), drawn = classes.filter(t => A.by[t].n && A.by[t].sc.tail_estimated >= 50);
    // seven lines are not told apart by colour alone: every second one is dashed, and the first (most reads) is thicker
    const ser = drawn.map((t, k) => { const o = A.by[t], n = D.sum(o.n) + o.over || 1; return { name: word(t), type: 'outline', edges, y: Array.from(o.n, v => v / n * 100), color: classColor(t, classes.indexOf(t), true), width: t === classes[0] ? 2 : 1.5, dash: k % 2 ? [5, 3] : null }; });
    const med = drawn.map(t => ({ t, m: tailOf(t).median })).filter(o => o.m != null).sort((a, x) => a.m - x.m), lo = med[0], hi = med[med.length - 1];
    tailSpec = { form: 'xy', x: { domain: [0, dom[1]], kind: 'linear', title: 'Poly(A) tail length, nt' }, y: { title: 'Share of each type’s reads with a tail, %' }, series: ser,
      vlines: overall.median != null ? [{ x: overall.median, label: 'median ' + nt(overall.median) + ' nt', strong: true }] : [], cursor: { x: Array.from({ length: edges.length - 1 }, (_, i) => (edges[i] + edges[i + 1]) / 2) },
      readout: i => { const by = ser.slice().sort((a, q) => q.y[i] - a.y[i]), shown = by.slice(0, 4);      // the types with the largest share at this length first
        return [[edges[i] + ' to ' + (edges[i + 1] - 1) + ' nt', '']].concat(shown.map(q => [q.y[i].toFixed(1) + ' %', mid(q.name, 18)])).concat(by.length > shown.length ? ['+ ' + (by.length - shown.length) + ' more in the table'] : []); } };
    tailLegend = ser.map(q => ({ name: q.name, kind: 'line', color: q.color, dash: !!q.dash }));
    tailReading = 'Half of the tails are ' + nt(overall.median) + ' nt or shorter' + (lo && hi && lo.t !== hi.t ? '; by type the median runs from ' + nt(lo.m) + ' nt (' + word(lo.t) + ') to ' + nt(hi.m) + ' nt (' + word(hi.t) + ')' : '') + '.';
    tailTable = () => ({ head: ['Transcript type', 'Reads', 'With a tail estimate', 'Median tail, nt', 'Middle half, nt', 'Mean tail, nt'],
      rows: ['all'].concat(classes).filter(t => A.by[t].sc.reads > 0).map(t => { const o = A.by[t], q = tailOf(t); return [t === 'all' ? 'All types' : word(t), F.int(o.sc.reads), F.int(o.sc.tail_estimated) + ' (' + F.pct1(share(o.sc.tail_estimated, o.sc.reads)) + ')', nt(q.median), q.q25 != null ? nt(q.q25) + ' to ' + nt(q.q75) : F.DASH, q.mean != null ? q.mean.toFixed(1) : F.DASH]; }),
      note: overall.exact ? 'Medians and quartiles are exact (read from 1 nt steps when the table was read).' : 'Medians and quartiles of several samples chosen together are read from the stored histogram (5 nt steps).' });
  }
  // ---- 4. tail estimate found, per sample ----
  const estRows = (sel.length > 1 ? [{ label: sumLabel, strong: true, parts: NO_TAIL.map(([, token, k]) => ({ value: tot[k], color: tk(token) })) }] : [])
    .concat(all.map(s => ({ label: s.name + mark(s.key), short: c.short(s.key, 22) + mark(s.key), strong: strong(s.key), parts: NO_TAIL.map(([, token, k]) => ({ value: perSample[s.key].by.all.sc[k], color: tk(token) })) })));
  const estBy = all.map(s => share(perSample[s.key].by.all.sc.tail_estimated, perSample[s.key].by.all.sc.reads)).filter(F.isNum);
  const readStack = (rows, names, unit, of) => k => { const r = rows[Math.floor(k / 1000)], i = k % 1000, t = r.parts.reduce((a, q) => a + q.value, 0);
    return [[r.short || r.label, ''], [pc(r.parts[i].value / (t || 1)), F.lower(names[i])], [F.int(r.parts[i].value), unit], 'of ' + F.count(t) + ' ' + of]; };
  // ---- 5. genes ----
  const gCells = {}; for (const x of (g && g.cells) || []) gCells[x.at[0]] = x;
  const det = {}; for (const r of ((((g || {}).run_level || {}).tables || {}).rna_genes_detected || { rows: [] }).rows) det[r[0]] = { reads: r[1], ge1: r[2], ge10: r[3], ge100: r[4], tx: r[5], top1: r[6], top10: r[7], top100: r[8] };
  const topRows = (((g || {}).run_level || {}).tables || {}).rna_top_gene_table, gAll = all.filter(s => gCells[s.key] && det[s.key]), gSel = gAll.filter(s => keys.has(s.key));
  const mito = gSel.length ? share(gSel.reduce((a, s) => a + (gCells[s.key].scalars.reads_mito_genes || 0), 0), gSel.reduce((a, s) => a + (gCells[s.key].scalars.reads_with_gene || 0), 0)) : null;
  const ge10 = gSel.map(s => det[s.key].ge10), t10 = gSel.map(s => share(det[s.key].top10, gCells[s.key].scalars.reads_with_gene)).filter(F.isNum);
  let rankSpec = null, rankLegend = null;
  const withCurve = gAll.filter(s => gCells[s.key].curves && gCells[s.key].curves.rna_gene_rank_curve);
  if (withCurve.length) {
    const few = withCurve.length <= 8, some = !few && withCurve.some(s => c.sKeys.includes(s.key)), lit = s => few || c.sKeys.includes(s.key), maxRank = Math.max(...withCurve.map(s => { const x = gCells[s.key].curves.rna_gene_rank_curve.x; return x[x.length - 1]; }));
    const ser = withCurve.map(s => { const v = gCells[s.key].curves.rna_gene_rank_curve, t = v.y.reads_cum[v.y.reads_cum.length - 1] || 1, i = c.named.findIndex(q => q.key === s.key);
      return { name: s.name, type: 'line', x: v.x, y: v.y.reads_cum.map(n => n / t * 100), v, lit: lit(s), color: lit(s) ? sampleColor(s, i) : some ? tk('--fc-chart-none-line') : tk('--fc-this-run'), width: lit(s) ? 1.75 : 1, alpha: few || some ? 1 : 0.45 }; })
      .sort((a, q) => (a.lit > q.lit) - (a.lit < q.lit));
    const ref = ser[ser.length - 1];
    rankSpec = { form: 'xy', x: { domain: [1, Math.max(10, maxRank)], log: true, kind: 'count', title: 'Genes, most reads first (rank, logarithmic)', fmt: v => F.int(v) }, y: { domain: [0, 100], title: 'Share of the sample’s reads in these genes, %' }, series: ser, cursor: { x: ref.x },
      readout: i => { const list = (few ? ser : ser.filter(q => q.lit)).filter(q => q.v.x[i] === ref.x[i]), shown = list.slice(0, 4);
        return [[ref.x[i] === 1 ? 'The gene with most reads' : 'The ' + F.int(ref.x[i]) + ' genes with most reads', '']].concat(shown.map(q => [q.y[i].toFixed(1) + ' %', mid(q.name, 14)])).concat(list.length > shown.length ? ['+ ' + (list.length - shown.length) + ' more in the table'] : []).concat(few ? [] : [ser.length + ' samples drawn in one colour: they are not told apart here, the table names each']); } };
    rankLegend = few ? ser.map(q => ({ name: q.name, label: q.name, kind: 'line', color: q.color })) : some ? ser.filter(q => q.lit).slice(0, 8).map(q => ({ name: q.name, label: q.name, kind: 'line', color: q.color })).concat([{ name: 'other', label: 'Other samples', kind: 'line', color: tk('--fc-chart-none-line') }]) : null;
  }
  const showTop = gSel[0], topOf = s => topRows ? topRows.rows.filter(r => r[0] === s.key) : [], topList = showTop ? topOf(showTop) : [];
  const topBars = topList.map(r => ({ label: r[2], value: r[3], text: F.count(r[3]), note: pc(share(r[3], gCells[showTop.key].scalars.reads_with_gene)) }));
  const rangeOf = (a, f) => a.length > 1 && f(Math.min(...a)) !== f(Math.max(...a)) ? f(Math.min(...a)) + ' to ' + f(Math.max(...a)) : a.length ? f(a[0]) : F.DASH;
  const reading = 'Of the ' + F.count(tot.reads) + ' reads given to a transcript' + (one ? ' in ' + sel[0].name : sel.length < all.length ? ' (' + F.plural(sel.length, 'sample') + ' chosen)' : all.length > 1 ? ' (' + all.length + ' samples)' : '')
    + (sharedSel ? ', merged from ' + shared.both + ' of the sample' : '') + ', '
    + F.pct1(typeShare('protein_coding')) + ' come from protein-coding transcripts and ' + pc(rrna) + ' from ribosomal RNA' + (hasTail ? '; ' + F.pct1(share(tot.tail_estimated, tot.reads)) + ' carry a poly(A) tail estimate, half of them ' + nt(overall.median) + ' nt or shorter.' : '. These tables hold no tail length.');
  return <Section id="rna_composition" title="RNA composition" reading={reading} x={{ evRule, evSec }}>
    {shared ? <p class="fc-section__note fc-notice fc-shared" role="status">{shared.text}</p> : null}
    {partial ? <p class="fc-section__note fc-notice fc-partial" role="status">Partial: {partial}</p> : null}
    <div class="fc-tiles fc-section__tiles">
      <Tile label="Reads given to a transcript" value={F.countParts(tot.reads)[0]} unit={F.countParts(tot.reads)[1]} caption={capS((one ? mid(sel[0].name, 24) : F.plural(sel.length, 'sample')) + ' · rows of the annotated read tables')} title="rna_biotype_reads" />
      <Tile label="Protein coding" value={F.pctNum(typeShare('protein_coding'))} unit="%" caption="of the reads given to a transcript" title="rna_protein_coding_frac" />
      <Tile label="Ribosomal RNA" value={pcNum(rrna)} unit="%" caption="rRNA and mitochondrial rRNA transcripts" title="rna_rrna_frac" />
      {mito != null ? <Tile label="Mitochondrial genes" value={pcNum(mito)} unit="%" caption="reads of genes named MT-…" title="rna_mito_frac" /> : null}
      {hasTail ? <Tile label="With a tail estimate" value={F.pctNum(share(tot.tail_estimated, tot.reads))} unit="%" caption={capS(F.count(tot.tail_estimated) + ' of ' + F.count(tot.reads) + ' reads given to a transcript')} title="rna_biotype_polya_frac" />
        : <Tile label="With a tail estimate" state="unknown" value="not known" caption="the tables hold no tail length" />}
      {hasTail ? <Tile label="Median poly(A) tail, annotated reads" value={nt(overall.median)} unit="nt" caption={'reads given to a transcript, with a tail estimate' + (overall.q25 != null ? ' · middle half ' + nt(overall.q25) + ' to ' + nt(overall.q75) + ' nt' : '')} title="rna_biotype_polya_median_nt" /> : null}
      {ge10.length ? <Tile label="Genes with 10 reads or more" value={rangeOf(ge10, F.int)} caption={gSel.length > 1 ? 'per sample, lowest to highest' : F.int(det[gSel[0].key].ge1) + ' genes with any read'} title="rna_genes_detected" /> : null}
      {t10.length ? <Tile label="Reads in the ten top genes" value={rangeOf(t10, F.pctNum)} unit="%" caption={gSel.length > 1 ? 'per sample, lowest to highest' : 'of the reads that name a gene'} title="rna_top_gene_share" /> : null}
    </div>
    {/* one bar per sample */}
    <Grid>
      <ChartPanel id="p-rna-types" title="Reads by transcript type" reading={top ? F.pct1(top.reads / reads) + ' of the reads are ' + word(top.name) + (second ? ', ' + F.pct1(second.reads / reads) + ' ' + word(second.name) : '') + '; ' + F.plural(A.types.length, 'type') + ' seen.' : ''} scope={st.scope} tag={tag} height={height(typeRows)}
        spec={{ form: 'hstack', norm: true, rowH: 22, rows: typeRows, readout: readStack(typeRows, classes.map(word), 'reads', 'reads given to a transcript') }} legend={classes.map((t, i) => ({ name: word(t), color: classColor(t, i), kind: t === 'other' ? 'none' : undefined }))}
        how={'Each bar is one sample. A read counts for the type of the transcript it was given to, as the annotation names it. Types that hold less than one read in a thousand of the run are drawn together as “Other types”; the ranked list beside this chart names them.' + (shared ? ' ' + shared.legend + '.' : '')}
        table={() => ({ head: ['Sample', 'Reads'].concat(classes.map(word)), rows: (sel.length > 1 ? [[sumLabel, F.int(tot.reads)].concat(classes.map(t => F.int(A.by[t].sc.reads)))] : []).concat(all.map(s => [s.name + mark(s.key), F.int(perSample[s.key].by.all.sc.reads)].concat(classes.map(t => F.int(perSample[s.key].by[t].sc.reads))))), note: noteS })} />
      {hasTail ? <ChartPanel id="p-rna-est" title="Tail estimate found, per sample" reading={'The basecaller estimated a tail for ' + F.pct1(share(tot.tail_estimated, tot.reads)) + ' of the reads' + (estBy.length > 1 ? ' (' + F.pct1(Math.min(...estBy)) + ' to ' + F.pct1(Math.max(...estBy)) + ' across ' + (chosen.length ? 'all ' + all.length : 'the') + ' samples)' : '') + '.'} scope={st.scope} tag={tag} height={height(estRows)}
        spec={{ form: 'hstack', norm: true, rowH: 22, rows: estRows, readout: readStack(estRows, NO_TAIL.map(x => x[0]), 'reads', 'reads given to a transcript') }} legend={NO_TAIL.filter(([, , k]) => tot[k] > 0).map(([name, token, k]) => ({ name, color: tk(token), kind: k === 'tail_no_value' ? 'none' : undefined }))}
        table={() => ({ head: ['Sample', 'Reads'].concat(NO_TAIL.map(x => x[0])), rows: (sel.length > 1 ? [[sumLabel, F.int(tot.reads)].concat(NO_TAIL.map(x => F.int(tot[x[2]])))] : []).concat(all.map(s => { const q = perSample[s.key].by.all.sc; return [s.name + mark(s.key), F.int(q.reads)].concat(NO_TAIL.map(x => F.int(q[x[2]]))); })), note: noteS })} /> : null}
    </Grid>
    {/* the samples chosen, added up; and the genes of each sample */}
    <Grid>
      <ChartPanel id="p-rna-rank" title="Every transcript type, most reads first" reading={(one ? sel[0].name : chosen.length ? F.plural(sel.length, 'sample') + ' chosen' : 'All samples') + ': ' + F.plural(A.types.length, 'transcript type') + '; the bars use the colours of the chart beside.'} scope={st.scope} tag={tag} height={rankRows.length * 22 + 10} lines={1}
        spec={{ form: 'hbars', rowH: 22, log: c.logY, rows: rankRows, readout: i => [[rankRows[i].label, ''], [F.int(rankRows[i].value), 'reads'], [rankRows[i].note, 'of the reads given to a transcript']] }}
        table={() => ({ head: ['Transcript type', 'Reads', 'Share', 'With a tail estimate'], rows: A.types.map(t => [word(t.name), F.int(t.reads), pc(t.reads / reads), hasTail ? F.int(t.tail) : F.DASH]) })} />
      {hasTail ? <ChartPanel id="p-rna-tail" title="Poly(A) tail length by transcript type" reading={tailReading} scope={st.scope} tag={tag} height={240} spec={tailSpec} legend={tailLegend} legendToggle table={tailTable}
        how="Each line is one transcript type, scaled to its own reads, so shapes compare; every second line is dashed, and the legend switches a line off. Mitochondrial transcripts have short tails; a type whose tails are much shorter than usual may be degraded. Types with fewer than 50 tails are left out of the picture and are in the table." /> : null}
      {rankSpec ? <ChartPanel id="p-rna-genes" title="How the reads spread over genes" reading={(t10.length ? 'The ten genes with most reads take ' + rangeOf(t10, F.pct1) + ' of a sample’s reads; ' + rangeOf(ge10, F.int) + ' genes have 10 reads or more.' : '') + (withCurve.length > 8 && !withCurve.some(s => c.sKeys.includes(s.key)) ? ' ' + withCurve.length + ' samples are drawn in one colour and are not told apart here: choose samples to colour them, or see the table.' : '')} tag={tag} height={240} spec={rankSpec} legend={rankLegend} legendToggle={withCurve.length <= 8}
        how="Genes are put in order of their reads. A line that rises early means few genes take most of the reads: a library of low complexity. Each line is one sample and ends at 100 % at its last gene."
        table={() => ({ head: ['Sample', 'Reads naming a gene', 'Genes with a read', 'with 10 or more', 'with 100 or more', 'Transcripts with a read', 'Top gene', 'Top 10 genes', 'Top 100 genes', 'Mitochondrial genes'],
          note: noteS, rows: gAll.map(s => { const d = det[s.key], q = gCells[s.key].scalars, w = q.reads_with_gene; return [s.name + mark(s.key), F.int(w), F.int(d.ge1), F.int(d.ge10), F.int(d.ge100), d.tx == null ? F.DASH : F.int(d.tx), F.pct1(share(d.top1, w)), F.pct1(share(d.top10, w)), F.pct1(share(d.top100, w)), F.pct1(share(q.reads_mito_genes, w))]; }) })} /> : null}
      {topBars.length ? <ChartPanel id="p-rna-topgenes" title={'The genes with the most reads in ' + showTop.name} reading={'Together ' + F.pct1(share(det[showTop.key].top10, gCells[showTop.key].scalars.reads_with_gene)) + ' of this sample’s reads that name a gene.' + (gSel.length > 1 ? ' This is the first sample: choose one in the Samples control to see its list; the table holds every sample.' : '')} tag={tag} height={topBars.length * 22 + 10} lines={1}
        spec={{ form: 'hbars', rowH: 22, rows: topBars, readout: i => [[topBars[i].label, ''], [F.int(topBars[i].value), 'reads'], [topBars[i].note, 'of the reads that name a gene']] }}
        table={() => ({ head: ['Sample', 'Rank', 'Gene', 'Reads', 'Share'], rows: gAll.flatMap(s => topOf(s).map(r => [s.name + mark(s.key), String(r[1]), r[2], F.int(r[3]), F.pct1(share(r[3], gCells[s.key].scalars.reads_with_gene))])), note: noteS })} /> : null}
    </Grid>
    <Disclosure id="n-rna" c={c} title={'All numbers (' + countNumbers(b, [BIOTYPE, GENES]) + ')'}><AllNumbers b={b} mods={[BIOTYPE, GENES]} /></Disclosure>
  </Section>;
}
