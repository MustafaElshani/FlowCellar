// Number, unit and date formats.
// Comma for thousands and point for decimals in every locale; a no-break space before a unit.
const NB = ' ';
const group = s => s.replace(/\B(?=(\d{3})+(?!\d))/g, ',');
export const DASH = '–';
export const isNum = v => typeof v === 'number' && isFinite(v);

export function int(v) { return isNum(v) ? group(String(Math.round(v))) : DASH; }
// three significant digits, one decimal below 100 ("208", "65.5", "2.5")
export function sig(v) {
  if (!isNum(v)) return DASH;
  const a = Math.abs(v);
  if (a >= 100) return group(String(Math.round(v)));
  if (a >= 10) return (Math.round(v * 10) / 10).toFixed(1);
  if (a >= 1) return (Math.round(v * 10) / 10).toFixed(1);
  if (a === 0) return '0';
  return String(+v.toPrecision(2));
}
export function basesParts(v) {
  if (!isNum(v)) return [DASH, ''];
  const a = Math.abs(v);
  if (a >= 1e9) return [sig(v / 1e9), 'Gb'];
  if (a >= 1e6) return [sig(v / 1e6), 'Mb'];
  if (a >= 1e3) return [sig(v / 1e3), 'kb'];
  return [String(Math.round(v)), 'b'];
}
export function bases(v) { const [n, u] = basesParts(v); return u ? n + NB + u : n; }
export function countParts(v) {
  if (!isNum(v)) return [DASH, ''];
  if (Math.abs(v) >= 1e6) return [sig(v / 1e6), 'M'];
  return [int(v), ''];
}
export function count(v) { const [n, u] = countParts(v); return u ? n + NB + u : n; }
// A share is never written as more than it is: one that is not zero is never "0.0", and one that is not the whole
// is never "100.0". Near the two ends it gets a second decimal, and beyond that the words "under" and "over".
const nearEnd = (frac, p) => frac > 0 && p < 0.05 ? (p >= 0.005 ? p.toFixed(2) : 'under 0.01') : frac < 1 && p >= 99.95 ? (p < 99.995 ? p.toFixed(2) : 'over 99.99') : null;
// share: whole numbers from 10 to 99, one decimal outside
export function pctNum(frac) {
  if (!isNum(frac)) return DASH;
  const p = frac * 100, end = nearEnd(frac, p);
  if (end) return end;
  if (p >= 10 && p < 99.5) return String(Math.round(p));
  return (Math.round(p * 10) / 10).toFixed(1);
}
export function pct(frac) { return isNum(frac) ? pctNum(frac) + NB + '%' : DASH; }
// share with one decimal
export function pct1(frac) { if (!isNum(frac)) return DASH; const p = frac * 100; return (nearEnd(frac, p) || p.toFixed(1)) + NB + '%'; }
export function hours(h) {
  if (!isNum(h)) return DASH;
  if (h < 1) return Math.max(1, Math.round(h * 60)) + NB + 'min';
  if (Math.abs(h - Math.round(h)) < 0.05) return Math.round(h) + NB + 'h';
  return h.toFixed(1) + NB + 'h';
}
export const q = v => isNum(v) ? 'Q' + NB + v.toFixed(1) : DASH;
export const fixed = (v, d = 1) => isNum(v) ? v.toFixed(d) : DASH;

const MONTHS = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'];
// Times. MinKNOW's files give UTC. Every time is shown on the clock of the computer the page is read on: the
// clock on the wall of whoever reads it, and for a lab that sits beside its sequencer also the clock in the name
// of the run folder. The zone is said beside the time or in its tooltip (zone, timeTitle), never left unsaid.
// `floating`: the time has no zone at all (it was read from a folder name) and is shown exactly as written, so
// that a run with raw data only and a run with reports show the same clock.
const two = n => String(n).padStart(2, '0');
function clock(iso, floating) {
  const t = new Date(iso); if (isNaN(t.getTime())) return null;
  return floating ? { y: t.getUTCFullYear(), mo: t.getUTCMonth(), d: t.getUTCDate(), h: t.getUTCHours(), mi: t.getUTCMinutes() }
    : { y: t.getFullYear(), mo: t.getMonth(), d: t.getDate(), h: t.getHours(), mi: t.getMinutes() };
}
export function date(iso, floating) { const c = iso ? clock(iso, floating) : null; return c ? c.d + ' ' + MONTHS[c.mo] + ' ' + c.y : DASH; }
export function dateTime(iso, floating) { const c = iso ? clock(iso, floating) : null; return c ? c.d + ' ' + MONTHS[c.mo] + ' ' + c.y + ', ' + two(c.h) + ':' + two(c.mi) : DASH; }
export function monthKey(iso, floating) { const c = iso ? clock(iso, floating) : null; return c ? MONTHS[c.mo] + ' ' + c.y : DASH; }
// "UTC", "UTC+1", "UTC−5:30": this computer's zone at that moment (summer time included)
export function zone(iso) {
  const t = new Date(iso || Date.now()), off = isNaN(t.getTime()) ? 0 : -t.getTimezoneOffset(); if (!off) return 'UTC';
  const a = Math.abs(off); return 'UTC' + (off > 0 ? '+' : '−') + Math.floor(a / 60) + (a % 60 ? ':' + two(a % 60) : '');
}
// the tooltip of a time: the same moment in UTC, and which clock the page shows
export function timeTitle(iso, floating) {
  if (!iso) return undefined; const t = new Date(iso); if (isNaN(t.getTime())) return undefined;
  if (floating) return 'From the name of the run folder: the sequencer’s own clock (its time zone is not recorded)';
  const utc = t.getUTCDate() + ' ' + MONTHS[t.getUTCMonth()] + ' ' + t.getUTCFullYear() + ', ' + two(t.getUTCHours()) + ':' + two(t.getUTCMinutes()) + ' UTC';
  return zone(iso) === 'UTC' ? utc + ' (this computer’s clock)' : utc + ' · shown on this computer’s clock, ' + zone(iso);
}
// "2026-11" -> "Nov 2026"; short: "Nov" (the year is said once, under January)
export function monthLabel(ym, short) { const m = /^(\d{4})-(\d{2})$/.exec(ym || ''); if (!m) return ym || DASH; return MONTHS[+m[2] - 1] + (short ? '' : ' ' + m[1]); }
// "4 Nov 2026" or "Jul 2025 to Nov 2026" for a span of dates
export function period(a, b) { if (!a) return DASH; if (!b || date(a) === date(b)) return date(a); const ma = monthKey(a), mb = monthKey(b); return ma === mb ? ma : ma + ' to ' + mb; }
export function bytes(n) { return !isNum(n) ? DASH : n >= 1e12 ? sig(n / 1e12) + NB + 'TB' : n >= 1e9 ? sig(n / 1e9) + NB + 'GB' : n >= 1e6 ? sig(n / 1e6) + NB + 'MB' : n >= 1e3 ? sig(n / 1e3) + NB + 'kB' : n + NB + 'B'; }
// estimated bases of a whole library: "4.1 Tb", "812 Gb"
export function basesBig(v) { if (!isNum(v)) return DASH; return Math.abs(v) >= 1e12 ? sig(v / 1e12) + NB + 'Tb' : bases(v); }

// axis tick helpers
export function niceTicks(lo, hi, n) {
  const span = hi - lo || 1, raw = span / Math.max(1, n), p = Math.pow(10, Math.floor(Math.log10(raw))), f = raw / p;
  const step = (f > 5 ? 10 : f > 2 ? 5 : f > 1 ? 2 : 1) * p, out = [];
  for (let v = Math.ceil(lo / step - 1e-9) * step; v <= hi + step * 1e-9; v += step) out.push(+v.toPrecision(12));
  return out;
}
export function tickLabel(v) {
  const a = Math.abs(v);
  if (a >= 1000) return group(String(+v.toPrecision(6)));
  return String(+v.toPrecision(6));
}
// middle truncation: keep the first `head` and the last `tail` characters
export function midCut(s, max = 24, tail = 10) {
  if (!s || s.length <= max) return s || '';
  return s.slice(0, Math.max(3, max - tail - 1)) + '…' + s.slice(-tail);
}
export const plural = (n, one, many) => n + ' ' + (n === 1 ? one : (many || one + 's'));

// Short names for the samples of one run. The names of a run mostly differ in a few words ("lib_101_plate_ctrl_day1",
// "lib_101_plate_trt_day1"): a cut in the middle keeps what they share and drops what tells them apart. Here the
// words every name has in the same place are left out, so the short names stay different: "101_ctrl", "101_trt".
// -> a function (name, max) that gives the name itself when it fits in max characters, else its short form.
export function shortNamer(names) {
  const list = [...new Set((names || []).filter(Boolean))], short = new Map();
  if (list.length > 1) {
    const parts = list.map(n => n.split(/([_\-. ]+)/)), words = parts.map(p => p.filter((_, i) => i % 2 === 0));
    const same = words.every(w => w.length === words[0].length);
    let made;
    if (same) { const keep = words[0].map((_, i) => words.some(w => w[i] !== words[0][i])); made = parts.map(p => { let out = ''; p.forEach((t, i) => { if (i % 2 === 0 && keep[i / 2]) out += (out ? (p[i - 1] && p[i - 1][0]) || '_' : '') + t; }); return out; }); }
    else { // names of different shapes: only the words at the start and at the end that all of them share
      let a = 0, z = 0; const min = Math.min(...words.map(w => w.length));
      while (a < min - 1 && words.every(w => w[a] === words[0][a])) a++;
      while (z < min - 1 - a && words.every(w => w[w.length - 1 - z] === words[0][words[0].length - 1 - z])) z++;
      made = parts.map(p => p.slice(2 * a, p.length - 2 * z).join('')); }
    if (made.every(Boolean) && new Set(made).size === list.length) list.forEach((n, i) => short.set(n, made[i]));
  }
  return (name, max = 24) => { if (!name || name.length <= max) return name || ''; const s = short.get(name) || name; return s.length <= max ? s : max < 8 ? s.slice(0, max) : midCut(s, max, Math.min(10, Math.floor(max / 2) - 1)); };
}
// a label inside a sentence: its first letter small, unless the label starts with a name written in capitals ("UMI found")
export const lower = s => s && /^[A-Z][a-z]/.test(s) ? s.charAt(0).toLowerCase() + s.slice(1) : s;
