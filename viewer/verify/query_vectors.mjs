// The search grammar of the Runs table against the shared vectors (tests/vectors/query/*.json), which the Python
// side (flowcellar ls --where) runs too: the two implementations must give the same runs for every query.
//
//   node verify/query_vectors.mjs [DIR]          default DIR: <repo>/tests/vectors/query
//
// A vector: {"name", "function": "query.filter", "args": {"rows", "queries", "now"}, "expect": [{"ids", "inactive"}]}.
// "inactive" are the terms that could not be read: the hub shows them struck through and never applies them.
import * as esbuild from 'esbuild';
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath, pathToFileURL } from 'node:url';

const here = path.dirname(fileURLToPath(import.meta.url)), viewerDir = path.resolve(here, '..'), repo = path.resolve(viewerDir, '..');
const dir = path.resolve(process.argv[2] || path.join(repo, 'tests', 'vectors', 'query'));
// the viewer's module, bundled for Node exactly as it is bundled for the page
const out = path.join(viewerDir, '.cache', 'query.test.mjs');
await esbuild.build({ absWorkingDir: viewerDir, entryPoints: ['src/hub/query.js'], bundle: true, format: 'esm', outfile: out, logLevel: 'warning', charset: 'utf8' });
const { makeEngine } = await import(pathToFileURL(out).href + '?t=' + Date.now());

const files = fs.existsSync(dir) ? fs.readdirSync(dir).filter(f => f.endsWith('.json')).sort() : [];
let queries = 0, failed = 0;
for (const f of files) {
  const v = JSON.parse(fs.readFileSync(path.join(dir, f), 'utf8'));
  if (v.function !== 'query.filter') { console.log('skip ' + f + ': function ' + v.function); continue; }
  const rows = v.args.rows, eng = makeEngine(rows, Date.parse(v.args.now));
  v.args.queries.forEach((q, k) => {
    const res = eng.run(q), got = { ids: res.rows.map(r => r.id), inactive: res.terms.filter(t => !t.ok).map(t => t.raw) }, want = v.expect[k];
    queries++;
    if (JSON.stringify(got) !== JSON.stringify(want)) { failed++; console.log('FAIL ' + f + ' ' + JSON.stringify(q) + '\n  got  ' + JSON.stringify(got) + '\n  want ' + JSON.stringify(want)); }
  });
}
console.log(files.length + ' vector files, ' + queries + ' queries, ' + failed + ' differ');
process.exit(failed || !queries ? 1 : 0);
