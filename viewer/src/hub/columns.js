// The columns of the Runs table and how they fit a width. Plain functions (no markup): checked in verify/units.mjs.
// hide: the order in which columns leave when the table is too narrow (1 first). Columns without it always stay.
export const COLS = [
  { id: 'tick', w: 28 },
  { id: 'health', label: 'Run health', w: 150, max: 186, sort: 'health' },
  { id: 'data', label: 'Read data', w: 124, sort: 'data', hide: 6 },
  { id: 'fc', label: 'Flow cell', sub: 'run name', w: 152, max: 300, sort: 'fc' },
  { id: 'started', label: 'Started', sub: 'ran · assay', w: 106, sort: 'started', hide: 7 },
  { id: 'spark', label: 'Output', sub: 'over time', w: 92, hide: 4 },
  { id: 'out', label: 'Output at 48 h', sub: 'judged unit · other', w: 112, sort: 'out', num: true, hide: 8 },
  { id: 'pores', label: 'Pores', sub: '1st scan', w: 56, hide: 3, sort: 'pores', num: true },
  { id: 'occ', label: 'Occup.', sub: '1st hour', w: 56, hide: 2, sort: 'occ', num: true },
  { id: 'passed', label: 'Bases', sub: 'read', w: 68, hide: 1, sort: 'passed', num: true, title: 'Bases in the reads FlowCellar read for this run. Which reads those are (passed, full-length, …) stands under Read data.' },
  { id: 'n50', label: 'Read', sub: 'N50', w: 56, hide: 1, sort: 'n50', num: true },
  { id: 'why', label: 'Differs from typical', w: 132, max: 252, want: 228, hide: 5 },
  { id: 'peek', w: 30 },
];
const LAST_LEVEL = 8;
const HIDE_NAME = { 1: 'Bases read, Read N50', 2: 'Occupancy', 3: 'Pores', 4: 'Output over time', 5: 'Differs from typical', 6: 'Read data', 7: 'Started', 8: 'Output at 48 h' };
export const MIN_TABLE_WIDTH = COLS.filter(c => !c.hide).reduce((a, c) => a + c.w, 0);
// Which columns fit: hide in the stated order until every column that is left is whole. Column widths are stated
// for text of the designed size and grow with `scale` (larger text in the browser), so that a number is never
// squeezed into a column made for smaller letters. `all`: every column, whatever the width (the table then
// scrolls sideways).
export function fitColumns(width, scale = 1, all = false) {
  let level = 0; const W = c => Math.round(c.w * scale), px = v => Math.round(v * scale);
  const need = (l, min) => COLS.filter(c => !(c.hide && c.hide <= l)).reduce((a, c) => a + (min ? W(c) : px(c.want || c.w)), 0);
  // the first columns go so that the words column can show its longest phrase; after that the least width decides
  if (!all) while (level < LAST_LEVEL && need(level) > width) { if (level >= 3 && need(level, true) <= width) break; level++; }
  const cols = COLS.filter(c => !(c.hide && c.hide <= level));
  let spare = Math.max(0, width - cols.reduce((a, c) => a + W(c), 0));
  const widths = {}; for (const c of cols) widths[c.id] = W(c);
  // spare width: first the words column up to its longest phrase, then the modifiers, then the run name
  const give = (id, upTo) => { const c = cols.find(x => x.id === id); if (!c) return; const add = Math.max(0, Math.min(spare, px(upTo) - widths[id])); widths[id] += add; spare -= add; };
  give('why', 228); give('health', 186); give('fc', 300); give('why', 252);
  // what a wide window still has to give goes to every column in proportion (the tick box and the peek button stay as they are): the name of a run
  // stays near its numbers instead of one column taking a gap of half a window
  if (spare > 0) {
    const flex = cols.filter(c => c.label), sum = flex.reduce((a, c) => a + widths[c.id], 0); let left = spare;
    for (const c of flex) { const add = Math.floor(spare * widths[c.id] / sum); widths[c.id] += add; left -= add; }
    widths.fc += left;
  }
  const wouldHide = []; { let l = 0; while (l < LAST_LEVEL && need(l) > width) { if (l >= 3 && need(l, true) <= width) break; l++; } for (let k = 1; k <= l; k++) if (COLS.some(c => c.hide === k)) wouldHide.push(k); }
  const hiddenLevels = []; for (let l = 1; l <= level; l++) if (COLS.some(c => c.hide === l)) hiddenLevels.push(l);
  return { cols, widths, hidden: COLS.length - cols.length, hiddenNames: hiddenLevels.map(l => HIDE_NAME[l]), total: cols.reduce((a, c) => a + widths[c.id], 0),
    min: Math.round(MIN_TABLE_WIDTH * scale), wouldHide: COLS.filter(c => c.hide && wouldHide.includes(c.hide)).length, wouldHideNames: wouldHide.map(l => HIDE_NAME[l]) };
}
// As many whole phrases as the cell is wide, and "+n more" for the rest (the cell's tooltip lists them all): a
// phrase is never cut in the middle.
export function fitWords(words, px) {
  if (!words.length) return '';
  const max = Math.max(8, Math.floor((px - 12) / 6.0)); let s = '', n = 0;
  for (const w of words) {
    const next = s ? s + ' · ' + w : w, rest = words.length - n - 1;
    if (n > 0 && next.length + (rest ? 10 : 0) > max) break;
    s = next; n++;
  }
  return n < words.length ? s + ' · +' + (words.length - n) + ' more' : s;
}
