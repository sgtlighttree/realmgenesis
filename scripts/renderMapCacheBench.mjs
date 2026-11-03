// Runs the F4 Map2D geometry-cache instrument (`scripts/perf/mapCacheBench.ts`).
//
// A script and not a test on purpose: at the default 200k cells it generates a
// world for ~20s before measuring. `npm test` must never pick it up.
//
// Usage:
//   node scripts/renderMapCacheBench.mjs [--points=200000] [--reps=3] [--out=FILE]
//
// Every timing line that builds paths carries a `hash=` of the path commands
// emitted. Run before and after a change and diff: the ms should move, the
// hashes must not.
import { createServer } from 'vite';
import { readFileSync } from 'fs';

const flags = Object.fromEntries(
  process.argv.slice(2).filter(a => a.startsWith('--')).map(a => {
    const [k, v] = a.replace(/^--/, '').split('=');
    return [k, v ?? true];
  }),
);
const out = String(flags.out ?? 'tmp/mapCacheBench.txt');
process.env.F4_OUT = out;
process.env.F4_POINTS = String(flags.points ?? 200000);
process.env.F4_REPS = String(flags.reps ?? 3);

const server = await createServer({ server: { middlewareMode: true }, appType: 'custom', logLevel: 'error' });
const { run } = await server.ssrLoadModule('/scripts/perf/mapCacheBench.ts');
await run();
await server.close();
console.log(readFileSync(out, 'utf8'));
