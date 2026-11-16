// What a page shows on paper only, and the control that takes it there (core/print.js does the work).
//
//   SavePdf     the control of a page: "Save as PDF", and on a page with closed groups of further charts "PDF with
//               every chart". Both open the print dialog of the browser.
//   PaperHead   the strip at the head of the first sheet: the mark, what these sheets are, when the numbers were
//               computed and by which version of FlowCellar, when the PDF was made.
//   PaperScope  the selection the sheets show, in words: on paper there is no controls bar to read it from.
//
// PaperHead and PaperScope draw nothing on screen: they are in the page only while it is laid out for paper.
import { useEffect } from 'preact/hooks';
import brand from '../gen/brand.js';
import { printing, printAll, savePdf, setPrintTitle } from '../core/print.js';
import { index, isReport } from '../core/data.js';
import * as F from '../core/format.js';

export function SavePdf({ charts }) {
  return <span class="fc-savepdf">
    <button type="button" class="fc-btn" data-fk="save-pdf" title="Opens the print dialog of the browser with a layout made for paper. Choose “Save as PDF” (in Firefox “Save to PDF”) as the printer there."
      onClick={() => savePdf(false)}>Save as PDF</button>
    {charts ? <button type="button" class="fc-btn fc-btn--ghost" data-fk="save-pdf-all" title="The same, and the closed groups of further charts on this page (“More charts”, “more tracks”) are opened for the PDF. Tables of numbers that are closed stay out."
      onClick={() => savePdf(true)}>PDF with every chart</button> : null}
  </span>;
}

// the lock-up in the colours it has on a light ground
const PaperLogo = () => <svg class="fc-paper__logo" viewBox={brand.viewBox} role="img" aria-label="FlowCellar" focusable="false">
  {brand.parts.map((p, k) => <path fill={brand.paper.parts[k] || p.fill} d={p.d} />)}
  <g transform={brand.wordTransform}>{brand.word.map((p, k) => <path fill={brand.paper.word[k] || p.fill} d={p.d} />)}</g>
</svg>;

// When and by which version the numbers of some runs or analyses were computed: one moment and one version, or the
// span of them. bundles: the data of what the page shows.
export function computedWords(bundles) {
  const times = [], versions = [];
  for (const b of bundles) { const p = (b && b.provenance) || {}, t = p.updated || p.created, v = (p.generator || {}).version; if (t) times.push(t); if (v && !versions.includes(v)) versions.push(v); }
  times.sort();
  const first = times[0], last = times[times.length - 1];
  const when = !first ? 'at a time that was not recorded' : F.date(first) === F.date(last) ? F.dateTime(last) + ' ' + F.zone(last) : 'between ' + F.date(first) + ' and ' + F.date(last);
  return 'Numbers computed ' + when + ' with FlowCellar ' + (versions.length ? versions.join(' and ') : '(version not recorded)');
}

// kind: "Run report", "Analysis report", "Comparison of 3 runs", "Flow cell history". file: the name offered for the PDF.
export function PaperHead({ kind, bundles, about, file }) {
  useEffect(() => { setPrintTitle(file); return () => setPrintTitle(''); }, [file]);
  if (!printing.value) return null;
  const ix = index.value, lib = !isReport && ix && ix.library && ix.library.name, now = new Date().toISOString();
  return <div class="fc-paper fc-paper__head">
    <PaperLogo />
    <div class="fc-paper__what"><b class="fc-paper__kind">{kind}</b>
      <span class="fc-paper__made">{computedWords(bundles)}{about ? ' · ' + about : ''}{lib ? ' · library ' + lib : ''} {'·'} PDF made {F.dateTime(now)} {F.zone(now)}</span></div>
  </div>;
}

// items: [[name of the control, what it is set to], ...]. groups: the page has groups that can be closed (a run page):
// the line then says what was left out because it was closed (core/print.js counts the groups once the page is drawn).
export function PaperScope({ items, groups }) {
  if (!printing.value) return null;
  return <div class="fc-paper fc-paper__scope">
    <p><b>This PDF shows</b> {items.filter(Boolean).map((it, k) => <>{k ? ' · ' : ''}<span class="fc-paper__item">{it[0]}: <b>{it[1]}</b></span></>)}</p>
    {groups ? <p class="fc-paper__left" data-every={printAll.value ? '1' : '0'} /> : null}
  </div>;
}
