// The search grammar, shared in the real tool with `flowcellar ls --where`:
//
//   query := term (space term)*          all terms must match
//   term  := ["-"] field ":" value | word | "quoted words"
//   value := item ("," item)*            any of the items
//   item  := text | number[unit] | cmp number[unit] | low ".." high | date
//
// Terms are evaluated over the row objects of the index; an unknown field or a value that cannot be read makes
// the term inactive (it is shown struck through) and never empties the table.
const UNIT = { bases: { gb: 1e9, mb: 1e6, kb: 1e3, b: 1 }, count: { m: 1e6, k: 1e3 }, pct: { '%': 0.01 }, hours: { h: 1, d: 24 }, num: { k: 1e3, m: 1e6 } };
export const fold = s => String(s == null ? '' : s).normalize('NFD').replace(/[̀-ͯ]/g, '').toLowerCase();

export function tokenize(q) {
  const out = [], re = /(-?)(?:([a-z_][a-z0-9_]*):)?(?:"([^"]*)"?|(\S+))/gi; let m;
  while ((m = re.exec(q || ''))) { const v = m[3] != null ? m[3] : m[4]; if (v == null || (v === '' && !m[2])) continue; out.push({ neg: m[1] === '-', field: m[2] ? m[2].toLowerCase() : null, value: v, raw: m[0] }); }
  return out;
}
function parseNum(s, unit) {
  const m = /^(<=|>=|<|>)?\s*(-?\d+(?:\.\d+)?)\s*([a-z%]*)$/i.exec(s); if (!m) return null;
  const table = UNIT[unit] || UNIT.num, u = m[3].toLowerCase(); let mul = u ? table[u] : (unit === 'pct' ? 0.01 : 1); if (u && mul == null) return null;
  const dec = (m[2].split('.')[1] || '').length;
  return { cmp: m[1] || null, v: +m[2] * mul, half: 0.5 * Math.pow(10, -dec) * mul };
}
function numTest(item, unit) {
  const r = item.split('..');
  if (r.length === 2) { const a = parseNum(r[0], unit), b = parseNum(r[1], unit); if (!a || !b) return null; return x => x >= a.v && x <= b.v; }
  const p = parseNum(item, unit); if (!p) return null;
  switch (p.cmp) { case '<': return x => x < p.v; case '<=': return x => x <= p.v; case '>': return x => x > p.v; case '>=': return x => x >= p.v; }
  return x => x >= p.v - p.half && x < p.v + p.half;       // a bare number matches what rounds to it
}
function dayRange(s) {
  const m = /^(\d{4})(?:-(\d{1,2}))?(?:-(\d{1,2}))?$/.exec(s); if (!m) return null;
  const y = +m[1], mo = m[2] ? +m[2] - 1 : 0, d = m[3] ? +m[3] : 1, a = Date.UTC(y, mo, d);
  if (mo < 0 || mo > 11 || d < 1 || new Date(a).getUTCMonth() !== mo) return null;        // 2026-13 or 2026-02-31 is not a date
  const b = m[3] ? Date.UTC(y, mo, d + 1) : m[2] ? Date.UTC(y, mo + 1, 1) : Date.UTC(y + 1, 0, 1);
  return [a, b];
}
function lev(a, b) { const d = Array.from({ length: a.length + 1 }, (_, i) => [i]); for (let j = 1; j <= b.length; j++) d[0][j] = j;
  for (let i = 1; i <= a.length; i++) for (let j = 1; j <= b.length; j++) d[i][j] = Math.min(d[i - 1][j] + 1, d[i][j - 1] + 1, d[i - 1][j - 1] + (a[i - 1] === b[j - 1] ? 0 : 1)); return d[a.length][b.length]; }

const gb = k => r => r[k] == null ? null : r[k] * 1e9;
// The same tables as flowcellar/library/query.py: the two are held equal by the vectors in tests/vectors/query.
const STATUS_WORDS = { warn: ['warn', 'warning'], pass: ['pass'], fail: ['fail'], incomplete: ['incomplete'], running: ['running'], unknown: ['unknown', 'not_assessed'] };
const PART_WORDS = { ...STATUS_WORDS, not_available: ['not_available', 'none', 'not', 'no', 'na'] };
const mods = r => r.modifiers || [];
const hasReport = r => r.detail_t1 != null && r.detail_t1 !== 'none';
// project: a value that is exactly the id or the name of a project of some run means that project and no other
// (ticking "alpha" must not also give "alpha_repeat"); any other value is looked for inside the names.
const projectsOf = r => [r.project].concat(r.project_ids || [], r.project_names || []).filter(v => v != null && v !== '').map(fold);
// field -> how to read it from a row
export const FIELDS = {
  status: { type: 'enum', get: r => r.verdict, words: STATUS_WORDS },
  health: { type: 'enum', get: r => r.health, words: PART_WORDS },
  data: { type: 'enum', get: r => r.data, words: PART_WORDS },
  fc: { type: 'text', get: r => r.fc, match: 'prefix' },
  id: { type: 'text', get: r => r.id, match: 'prefix' },
  project: { type: 'project', get: r => [r.project || ''].concat(r.project_ids || []).join(' ') },
  root: { type: 'text', get: r => (r.roots || []).join(' ') },
  sample: { type: 'text', get: r => (r.samples || []).join(' ') },
  name: { type: 'text', get: r => r.run_name },
  exp: { type: 'text', get: r => r.experiment },
  assay: { type: 'text', get: r => (r.assay_group || '') + ' ' + (r.assay || '') },
  kit: { type: 'text', get: r => r.kit },
  fctype: { type: 'text', get: r => r.product },
  device: { type: 'text', get: r => r.device_type },
  pos: { type: 'text', get: r => r.position },
  minknow: { type: 'text', get: r => r.minknow, match: 'prefix' },
  detail: { type: 'text', get: r => r.detail_t1 },
  check: { type: 'text', get: r => (r.reason_checks || []).join(' ') },
  yield: { type: 'num', unit: 'bases', get: gb('gb') },
  yield24: { type: 'num', unit: 'bases', get: gb('gb_24h') },
  yield48: { type: 'num', unit: 'bases', get: gb('gb_48h') },
  yield72: { type: 'num', unit: 'bases', get: gb('gb_72h') },
  passed: { type: 'num', unit: 'bases', get: gb('t2_passed_gb') },
  reads: { type: 'num', unit: 'count', get: r => r.reads },
  pores: { type: 'num', unit: 'num', get: r => r.pores },
  hours: { type: 'num', unit: 'hours', get: r => r.hours },
  n50: { type: 'num', unit: 'bases', get: r => r.t2_n50 != null ? r.t2_n50 : r.n50_est },
  q: { type: 'num', unit: 'num', get: r => r.t2_q_med },
  occ: { type: 'num', unit: 'pct', get: r => r.occ_first },
  usable: { type: 'num', unit: 'pct', get: r => r.t2_usable },
  samples: { type: 'num', unit: 'num', get: r => r.t2_samples },
  use: { type: 'num', unit: 'num', get: r => r.use },
  is: { type: 'flag', flags: {
    reused: r => (r.use || 0) > 1, 'first-use': r => r.use === 1 && hasReport(r), restart: r => !!r.restart,
    // counting: checks of this run still wait for the read count. A run whose read table is cut off is in quick
    // detail for good and is not "counting".
    counting: r => r.detail_t1 === 'quick' && (r.pending || 0) > 0 && !mods(r).includes('reads_partly_counted') && !mods(r).includes('reads_not_counted'), quick: r => r.detail_t1 === 'quick',
    cut: r => mods(r).includes('reads_partly_counted') || mods(r).includes('reads_not_counted'), conflict: r => (r.t2_conflicts || 0) > 0, waiting: r => (r.t2_files_waiting || 0) > 0,
    full: r => r.detail_t1 === 'full', early: r => mods(r).includes('baseline_provisional'), raw: r => r.cause === 'raw_only',
    running: r => r.run_state === 'running' || r.verdict === 'running', offline: r => r.source_state === 'offline' || r.source_state === 'missing', stale: r => mods(r).includes('judged_before_older_runs') } },
  has: { type: 'flag', flags: { reads: r => r.read_data != null && r.read_data !== 'absent', history: r => !mods(r).includes('no_history') && hasReport(r), report: hasReport, project: r => (r.project_ids || []).length > 0 || !!r.project } },
};
const ALIASES = { verdict: 'status', flowcell: 'fc', flow_cell: 'fc', run: 'name', experiment: 'exp', product: 'fctype', position: 'pos', read_data: 'data', occupancy: 'occ' };
const DATE_FIELDS = new Set(['after', 'before', 'on', 'year', 'last']);
// what a word without a field looks through
const HAY = ['id', 'fc', 'run_name', 'experiment', 'project', 'kit', 'product', 'assay', 'assay_group', 'device_type', 'position', 'folder', 'headline'];
export const hayOf = r => fold(HAY.map(k => r[k]).concat(r.samples || [], r.project_ids || [], r.roots || []).filter(Boolean).join(' '));
const startedMs = r => { if (r.startedMs != null) return r.startedMs || null; const t = Date.parse(r.started); return t === t ? t : null; };

// the words a verdict is shown with, for a hint when one is typed as plain text
const VERDICT_TERMS = { warning: 'warn', warnings: 'warn', warn: 'warn', pass: 'pass', passed: 'pass', fail: 'fail', failed: 'fail', incomplete: 'incomplete', running: 'running' };
export function verdictHint(q) {
  const toks = tokenize(q), hit = toks.find(t => !t.field && !t.neg && VERDICT_TERMS[fold(t.value)]); if (!hit) return null;
  const term = 'status:' + VERDICT_TERMS[fold(hit.value)];
  return { term, query: toks.map(t => t === hit ? term : t.raw).join(' ') };
}

export function makeEngine(rows, now) {
  const columns = new Set(); for (const r of rows) for (const k of Object.keys(r)) columns.add(k);
  const projectNames = new Set(); for (const r of rows) for (const p of projectsOf(r)) projectNames.add(p);
  function test(t) {
    if (!t.field) { const w = fold(t.value); return { fn: r => (r.hay != null ? r.hay : hayOf(r)).includes(w) }; }
    let f = ALIASES[t.field] || t.field, flip = false;
    if (f === 'no') { f = 'has'; flip = true; }
    if (DATE_FIELDS.has(f)) {
      let lo = -Infinity, hi = Infinity;
      if (f === 'last') { const m = /^(\d+)([dwm])$/i.exec(t.value); if (!m) return { error: 'expected a span such as 14d, 8w or 6m' }; lo = now - +m[1] * { d: 1, w: 7, m: 30.44 }[m[2].toLowerCase()] * 864e5; }
      else { const r = dayRange(t.value); if (!r) return { error: 'expected a date such as 2026, 2026-11 or 2026-11-04' }; if (f === 'after') lo = r[0]; else if (f === 'before') hi = r[0]; else { lo = r[0]; hi = r[1]; } }
      return { field: f, fn: r => { const ms = startedMs(r); return ms != null && ms >= lo && ms < hi; } };
    }
    let s = FIELDS[f];
    if (!s && columns.has(f)) s = { type: 'auto', get: r => r[f] };          // any run column of the index, by its id
    if (!s) { let best = null, bd = 3; for (const k of Object.keys(FIELDS).concat(Object.keys(ALIASES), [...DATE_FIELDS], 'no')) { const d = lev(f, k); if (d < bd) { bd = d; best = k; } } return { error: `'${t.field}' is not a field.` + (best ? ` Did you mean ${best}?` : '') }; }
    const items = t.value.split(',').filter(Boolean); if (!items.length) return { error: 'no value' };
    if (s.type === 'auto') { const tests = items.map(it => numTest(it, 'num')), want = items.map(fold);
      return { field: f, fn: r => { const v = s.get(r); if (v == null) return false; if (typeof v === 'boolean') return want.includes(String(v));
        if (typeof v === 'number') return tests.some(fn => fn && fn(v)); const text = fold(Array.isArray(v) ? v.join(' ') : v); return want.some(w => text.includes(w)); } }; }
    if (s.type === 'num') { const tests = items.map(it => numTest(it, s.unit)); if (tests.some(x => !x)) return { error: 'could not read the number' }; return { field: f, fn: r => { const x = s.get(r); return typeof x === 'number' && x === x && tests.some(fn => fn(x)); } }; }
    if (s.type === 'enum') {
      // a value that is no verdict is an error, as for is: and has: (a misspelt "wran" must not read as "no warnings")
      const want = items.map(fold), known = Object.keys(s.words).concat(...Object.values(s.words)), bad = want.find(x => !known.includes(x));
      if (bad != null) { let best = null, bd = 3; for (const k of known) { const d = lev(bad, k); if (d < bd) { bd = d; best = k; } } return { error: `'${bad}' is not a value of ${f}; one of ${Object.keys(s.words).join(', ')}.` + (best ? ` Did you mean ${best}?` : '') }; }
      return { field: f, fn: r => { const v = s.get(r), w = (s.words && s.words[v]) || [fold(v)]; return want.some(x => w.includes(x) || fold(v) === x); } };
    }
    if (s.type === 'project') { const want = items.map(fold); return { field: f, fn: r => { const mine = projectsOf(r), text = mine.join(' '); return want.some(w => projectNames.has(w) ? mine.includes(w) : text.includes(w)); } }; }
    if (s.type === 'text') { const want = items.map(fold); return { field: f, fn: r => { const v = fold(s.get(r)); return want.some(w => s.match === 'prefix' ? v.startsWith(w) : v.includes(w)); } }; }
    if (s.type === 'flag') { const fns = items.map(it => s.flags[fold(it)]); if (fns.some(x => !x)) return { error: 'unknown value; one of ' + Object.keys(s.flags).join(', ') }; return { field: f, fn: r => fns.some(fn => fn(r)) !== flip }; }
    return { error: 'unsupported field' };
  }
  const cache = new Map();
  function run(q, skipField) {
    const terms = tokenize(q).map(t => { let c = cache.get(t.raw); if (!c) { c = test(t); cache.set(t.raw, c); } return { raw: t.raw, field: c.field || (ALIASES[t.field] || t.field), value: t.value, neg: t.neg, ok: !!c.fn, fn: c.fn, error: c.error }; });
    const live = terms.filter(t => t.ok && (!skipField || t.field !== skipField));
    const out = rows.filter(r => live.every(t => t.fn(r) !== t.neg));
    return { terms, rows: out, count: out.length };
  }
  // counts per value of a field, under every active term except those on that same field; a row with several
  // values (the projects of a run) counts under each
  function facet(q, field, get) { const res = run(q, field), counts = new Map(); for (const r of res.rows) { const v = get(r); for (const x of Array.isArray(v) ? v : [v]) counts.set(x, (counts.get(x) || 0) + 1); } return counts; }
  return { run, facet };
}
// add or remove "field:value" in a query string; several values of one field become "field:a,b"
export function toggleTerm(q, field, value) {
  const toks = tokenize(q), v = /[\s,]/.test(value) ? '"' + value + '"' : value; let found = false; const out = [];
  for (const t of toks) {
    if (t.field === field && !t.neg) {
      found = true; const items = t.value.split(',').filter(Boolean), i = items.findIndex(x => fold(x) === fold(value));
      if (i >= 0) items.splice(i, 1); else items.push(value);
      if (items.length) out.push(field + ':' + items.map(x => /\s/.test(x) ? '"' + x + '"' : x).join(','));
    } else out.push(t.raw);
  }
  if (!found) out.push(field + ':' + v);
  return out.join(' ');
}
export const hasTerm = (q, field, value) => tokenize(q).some(t => t.field === field && !t.neg && t.value.split(',').some(x => fold(x) === fold(value)));
export const removeTerm = (q, raw) => tokenize(q).filter(t => t.raw !== raw).map(t => t.raw).join(' ');
