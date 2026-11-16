// Small checks of the viewer's plain functions, without a browser: how a number is written from its id, and how a
// time is shown on the clock of the computer the page is read on.
//
//   node verify/units.mjs            runs itself once per time zone (TZ) and reports
//
// Exit code 0 when everything holds.
import * as esbuild from 'esbuild';
import { spawnSync } from 'node:child_process';
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath, pathToFileURL } from 'node:url';

const here = path.dirname(fileURLToPath(import.meta.url)), viewerDir = path.resolve(here, '..'), repo = path.resolve(viewerDir, '..');
const ZONES = { 'UTC': ['17 Apr 2025, 12:45', 'UTC', '30 Jun 2025'], 'Europe/London': ['17 Apr 2025, 13:45', 'UTC+1', '1 Jul 2025'], 'Asia/Tokyo': ['17 Apr 2025, 21:45', 'UTC+9', '1 Jul 2025'],
  'America/Los_Angeles': ['17 Apr 2025, 05:45', 'UTC−7', '30 Jun 2025'], 'Asia/Kolkata': ['17 Apr 2025, 18:15', 'UTC+5:30', '1 Jul 2025'] };

if (!process.env.FC_UNITS_CHILD) {
  let failed = 0;
  for (const tz of Object.keys(ZONES)) {
    const r = spawnSync(process.execPath, [fileURLToPath(import.meta.url)], { env: { ...process.env, TZ: tz, FC_UNITS_CHILD: '1' }, encoding: 'utf8' });
    process.stdout.write(r.stdout); process.stderr.write(r.stderr); if (r.status !== 0) failed++;
  }
  console.log(failed ? 'FAILED in ' + failed + ' time zone(s)' : 'units: all checks passed in ' + Object.keys(ZONES).length + ' time zones');
  process.exit(failed ? 1 : 0);
}

const tz = process.env.TZ, out = path.join(viewerDir, '.cache', 'units.' + tz.replace(/\W/g, '_') + '.mjs');
await esbuild.build({ absWorkingDir: viewerDir, stdin: { contents: "export * as F from './src/core/format.js'; export { fmtValue } from './src/core/values.js'; export { fitColumns, fitWords } from './src/hub/columns.js'; export { columnsFor } from './src/charts/columns.js'; export * as S from './src/analysis/sc.js'; export { stepValue } from './src/core/derive_align.js'; export { Chart } from './src/charts/kit.js';", resolveDir: viewerDir, loader: 'js' },
  bundle: true, format: 'esm', outfile: out, logLevel: 'warning', charset: 'utf8' });
// the chart kit reads its colours from the page: here a token answers with its own name, so a test can see which one was used
globalThis.document = { documentElement: {} }; globalThis.getComputedStyle = () => ({ getPropertyValue: name => name });
const { F, fmtValue, fitColumns, fitWords, columnsFor, S, stepValue, Chart } = await import(pathToFileURL(out).href + '?t=' + Date.now());
let bad = 0; const eq = (name, got, want) => { if (got !== want) { bad++; console.log('FAIL [' + tz + '] ' + name + ': got ' + JSON.stringify(got) + ', want ' + JSON.stringify(want)); } };

// ---- times: this computer's clock, the zone never left unsaid ----
const [local, zone, lateDay] = ZONES[tz], iso = '2025-04-17T12:45:00Z';
eq('dateTime', F.dateTime(iso), local);
eq('zone', F.zone(iso), zone);
eq('a time from a folder name is shown as written', F.dateTime(iso, true), '17 Apr 2025, 12:45');
eq('a run started late in the evening (UTC) is on the local day', F.date('2025-06-30T23:30:00Z'), lateDay);
eq('monthKey follows the local day', F.monthKey('2025-06-30T23:30:00Z'), lateDay.slice(lateDay.indexOf(' ') + 1));
const title = F.timeTitle(iso);
if (!/17 Apr 2025, 12:45 UTC/.test(title) || (zone !== 'UTC' && !title.includes(zone))) { bad++; console.log('FAIL [' + tz + '] timeTitle: ' + title); }
if (!/folder/.test(F.timeTitle(iso, true))) { bad++; console.log('FAIL [' + tz + '] timeTitle of a floating time: ' + F.timeTitle(iso, true)); }
eq('no time', F.dateTime(null), F.DASH);
eq('winter zone in London is UTC', tz === 'Europe/London' ? F.zone('2026-11-13T19:47:00Z') : 'UTC', 'UTC');

// ---- numbers: the unit follows the end of the id ----
const cat = JSON.parse(fs.readFileSync(path.join(repo, 'src', 'flowcellar', 'catalogue', 'metrics_nanopore.json'), 'utf8'));
const ids = cat.metrics.concat(cat.derived || []).map(m => m.id);
const hourIds = ids.filter(id => /_h$/.test(id) && !/_per_h$/.test(id));
if (hourIds.length < 4) { bad++; console.log('FAIL: the catalogue names fewer than four numbers in hours: ' + hourIds.join(', ')); }
for (const id of hourIds) { const t = fmtValue(id, 81.17); if (!/^81\.2\sh$/.test(t)) { bad++; console.log('FAIL ' + id + ' = 81.17 is written ' + JSON.stringify(t) + ', not in hours'); } }
eq('fc_duty_half_life_h', fmtValue('fc_duty_half_life_h', 81.17).replace(/\s/g, ' '), '81.2 h');
eq('fc_sequencing_duty stays a share', fmtValue('fc_sequencing_duty', 0.5).replace(/\s/g, ' '), '50 %');
eq('occupancy', fmtValue('fc_occupancy_first_hour', 0.963).replace(/\s/g, ' '), '96 %');
eq('seconds', fmtValue('dyn_duration_median_s', 2.5).replace(/\s/g, ' '), '2.5 s');
eq('degrees', fmtValue('fc_temperature_max_dev_c', 0.44), '0.4 °C');
eq('a poly(A) tail is nucleotides', fmtValue('rna_biotype_polya_median_nt', 124), '124 nt');
eq('a poly(A) tail is nucleotides, whatever the id ends in', fmtValue('rna_polya_median', 123.5), '123.5 nt');
eq('a share that is not zero is never written as 0.0', fmtValue('cdna_rescued_frac', 0.000315).replace(/\s/g, ' '), '0.03 %');
eq('a share of zero stays zero', fmtValue('cdna_rescued_frac', 0).replace(/\s/g, ' '), '0.0 %');
eq('voltage drift', fmtValue('fc_bias_voltage_drift_mv_per_h', -0.31), '-0.31 mV per hour');
eq('output rate', fmtValue('yield_rate_peak_gb_per_h', 4.2).replace(/\s/g, ' '), '4.2 Gb per hour');
eq('Mb per pore', fmtValue('yield_per_pore_mb_48h', 12.34).replace(/\s/g, ' '), '12.3 Mb');
eq('N50', fmtValue('len_est_n50', 2752).replace(/\s/g, ' '), '2.8 kb');
eq('the stage that counts as usable is named as on the chart', fmtValue('funnel_usable_stage', 'quantified'), 'counted towards a transcript');
eq('missing', fmtValue('len_n50', null), F.DASH);

// ---- the Runs table: columns grow with the text, and whole phrases only ----
const one = fitColumns(1300, 1), two = fitColumns(1300, 2);
if (!(two.cols.length < one.cols.length)) { bad++; console.log('FAIL: text twice as large must hide columns at the same width (' + one.cols.length + ' -> ' + two.cols.length + ')'); }
for (const c of two.cols) if (two.widths[c.id] < 2 * c.w - 1) { bad++; console.log('FAIL: column ' + c.id + ' is ' + two.widths[c.id] + ' px at twice the text size, designed ' + c.w); }
eq('every column on request', fitColumns(600, 1, true).hidden, 0);
eq('whole phrases, and how many more', fitWords(['re-used · use 2', 'no history', 'reads partly counted'], 150), 're-used · use 2 · +2 more');

// ---- one check per defect once found in a plain function ----
const nb = t => String(t).replace(/\s/g, ' ');
// a share is never written as more than it is
eq('25,140 of 54 M reads is not 0.0 %', nb(F.pct1(25140 / 54e6)), '0.05 %');
eq('a few reads in millions', nb(F.pct1(22 / 54e6)), 'under 0.01 %');
eq('99.995 % is not 100.0 %', nb(F.pct1(0.99995)), 'over 99.99 %');
eq('99.96 % keeps its second decimal', nb(F.pct(0.9996)), '99.96 %');
eq('the whole is the whole', nb(F.pct1(1)), '100.0 %');
eq('nothing is nothing', nb(F.pct1(0)), '0.0 %');
eq('a label inside a sentence keeps a name in capitals', F.lower('UMI found and kept'), 'UMI found and kept');
eq('a label inside a sentence starts small', F.lower('Given to Space Ranger'), 'given to Space Ranger');
// short names of the samples of one run stay different, whatever room they get
{ const names = ['lib_101_plate_ctrl_day1', 'lib_101_plate_trt_day1', 'lib_102_plate_ctrl_day1', 'lib_102_plate_trt_day1', 'lib_103_plate_ctrl_day1', 'lib_103_plate_trt_day1'], short = F.shortNamer(names);
  for (const max of [8, 10, 12, 16, 20, 22]) { const got = names.map(n => short(n, max)); if (new Set(got).size !== names.length) { bad++; console.log('FAIL short names at ' + max + ' characters are not distinct: ' + got.join(' | ')); } }
  eq('the words all samples share are left out', short(names[0], 12), '101_ctrl');
  eq('a name that fits is written whole', short(names[0], 40), names[0]);
  const odd = F.shortNamer(['alpha_ctrl', 'alpha_treated_with_a_long_name', 'beta']); eq('names of different shapes: only the shared start is dropped, else the name is cut', odd('alpha_treated_with_a_long_name', 40), 'alpha_treated_with_a_long_name'); }
// rows of charts: equal rows, never three and one
eq('4 charts, room for 3', columnsFor(4, 3), 2); eq('4 charts, room for 4', columnsFor(4, 4), 4); eq('6 charts, room for 4', columnsFor(6, 4), 3); eq('6 charts, room for 5', columnsFor(6, 5), 3);
eq('10 charts, room for 5', columnsFor(10, 5), 5); eq('10 charts, room for 3', columnsFor(10, 3), 3); eq('2 charts, room for 7', columnsFor(2, 7), 2); eq('5 charts, room for 4', columnsFor(5, 4), 3); eq('1 chart', columnsFor(1, 5), 1);
// analyses: a probed share is "about", a count read off the map is a whole number, a cell without a value says of what
eq('a probed share', S.aboutShare(0.555), 'about 56 %'); eq('a probed share of everything', S.aboutShare(1), 'about 100 %'); eq('no share', S.aboutShare(null), '');
eq('genes on a spot are whole', S.channelText('genes', 37.6, 'spot'), 'about 38 genes'); eq('few genes on a spot', S.channelText('genes', 5.0, 'spot'), '5 genes'); eq('one gene', S.channelText('genes', 1, 'spot'), '1 gene');
eq('no data says of what', S.channelText('genes', NaN, 'spot'), 'genes: no data'); eq('no data of UMIs', S.channelText('umis', NaN, 'spot'), 'UMIs: no data');
{ const ch = { transform: 'linear', scale: 1, offset: 0, nodata: true }, r = S.channelRange('umis', ch, { hi: 100, top: 180 }), r2 = S.channelRange('umis', ch, { hi: 100, top: 100 });
  if (!/or more/.test(r.unit) || /or more/.test(r2.unit)) { bad++; console.log('FAIL the end of a colour scale below the largest cell must say "or more": ' + JSON.stringify([r, r2])); } }
// whole-number qualities: the value of the step that holds the rank
eq('median of whole-number qualities', stepValue([17, 18, 19, 20], [10, 60, 30], 0.5), 18); eq('upper quartile', stepValue([17, 18, 19, 20], [10, 60, 30], 0.75), 19);
// a line over heat cells is drawn on a casing in the colour of the chart surface (it vanished on the densest cells)
{ const calls = [], g = new Proxy({}, { get: (t, k) => k in t ? t[k] : (...a) => { calls.push([k, t.strokeStyle, t.lineWidth, ...a]); }, set: (t, k, v) => { t[k] = v; return true; } });
  const L = { X: v => v, Y: v => v, py: 0, ph: 100 }, q = { x: [0, 1, 2], y: [1, 2, 3], color: 'ink', width: 2 };
  Chart.prototype.paintSeries.call({}, g, q, L, true);
  const strokes = calls.filter(c => c[0] === 'stroke');
  if (!(strokes.length === 2 && strokes[0][1] === '--fc-chart-surface' && strokes[0][2] === 5 && strokes[1][1] === 'ink' && strokes[1][2] === 2)) { bad++; console.log('FAIL a line over heat cells needs a casing under it: ' + JSON.stringify(strokes)); }
  calls.length = 0; Chart.prototype.paintSeries.call({}, g, q, L, false);
  if (calls.filter(c => c[0] === 'stroke').length !== 1) { bad++; console.log('FAIL a line on a plain chart has no casing'); } }
// a histogram whose bars are narrower than a few pixels is one filled outline, not one rectangle per bar
{ const calls = [], g = new Proxy({}, { get: (t, k) => k in t ? t[k] : (...a) => { calls.push(k); }, set: (t, k, v) => { t[k] = v; return true; } });
  const n = 400, e = Array.from({ length: n + 1 }, (_, i) => i), L = { X: v => v * 0.8, Y: v => 100 - v, py: 0, ph: 100 };
  Chart.prototype.paintSeries.call({}, g, { type: 'hist', edges: e, y: e.slice(0, n).map(i => 10 + (i % 7)) }, L, false);
  if (calls.includes('fillRect') || calls.filter(c => c === 'fill').length !== 1) { bad++; console.log('FAIL dense histogram: ' + calls.filter(c => c === 'fillRect').length + ' rectangles, ' + calls.filter(c => c === 'fill').length + ' fills'); }
  calls.length = 0; Chart.prototype.paintSeries.call({}, g, { type: 'hist', edges: [0, 1, 2, 3], y: [1, 2, 3] }, { X: v => v * 40, Y: v => 100 - v, py: 0, ph: 100 }, false);
  if (calls.filter(c => c === 'fillRect').length !== 3) { bad++; console.log('FAIL wide bars stay separate bars'); } }
if (bad) process.exit(1);
