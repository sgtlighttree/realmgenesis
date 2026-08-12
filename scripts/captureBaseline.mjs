// Captures a worldDigest baseline into tmp/ for before/after refactor comparison.
// Deliberately NOT committed as a fixture: Math.sin/cos/pow are
// implementation-defined, so a committed bit-exact baseline drifts across V8
// versions. Same engine, same session, same machine is the only honest use.
// Usage: node scripts/captureBaseline.mjs before
//        node scripts/captureBaseline.mjs after && node scripts/captureBaseline.mjs compare
// The recorded gitSha is stamped from HEAD, which does NOT capture working-tree
// state — a `before` baseline captured after editing (uncommitted changes) is
// indistinguishable from one captured correctly, and the comparison then proves
// nothing. Whenever the result actually matters, capture `before` from a
// pristine `git worktree add --detach <pre-change-sha>` (symlinking
// node_modules in), not from the working tree.
import { mkdirSync, writeFileSync, readFileSync } from 'node:fs';
import { execSync } from 'node:child_process';
import { createServer } from 'vite';

const label = process.argv[2];
if (!['before', 'after', 'compare'].includes(label)) {
  console.error('usage: node scripts/captureBaseline.mjs before|after|compare');
  process.exit(1);
}

if (label === 'compare') {
  const a = JSON.parse(readFileSync('tmp/baseline-before.json', 'utf8'));
  const b = JSON.parse(readFileSync('tmp/baseline-after.json', 'utf8'));
  console.log(`before: ${a.capturedAt} @ ${a.gitSha}`);
  console.log(`after:  ${b.capturedAt} @ ${b.gitSha}`);
  const keys = [...new Set([...Object.keys(a.digest), ...Object.keys(b.digest)])].sort();
  const diffs = keys.filter(k => a.digest[k] !== b.digest[k]);
  if (diffs.length === 0) { console.log('IDENTICAL — 0 fields differ'); process.exit(0); }
  console.error(`DIFFERS in ${diffs.length} field(s):\n  ${diffs.join('\n  ')}`);
  process.exit(1);
}

const server = await createServer({ server: { middlewareMode: true } });
const { generateWorld } = await server.ssrLoadModule('/utils/worldGen.ts');
const { digestWorld } = await server.ssrLoadModule('/tests/helpers/worldDigest.ts');
const { makeParams } = await server.ssrLoadModule('/tests/helpers.ts');

const worlds = {};
for (const [name, params] of [['small', makeParams()], ['medium', makeParams({ points: 5000, erosionIterations: 3 })]]) {
  const w = await generateWorld(params);
  const d = digestWorld(w);
  for (const [k, v] of Object.entries(d)) worlds[`${name}.${k}`] = v;
}
await server.close();

mkdirSync('tmp', { recursive: true });
writeFileSync(`tmp/baseline-${label}.json`, JSON.stringify({
  capturedAt: new Date().toISOString(),
  gitSha: execSync('git rev-parse HEAD').toString().trim(),
  node: process.version,
  digest: worlds,
}, null, 2));
console.log(`wrote tmp/baseline-${label}.json (${Object.keys(worlds).length} fields)`);
