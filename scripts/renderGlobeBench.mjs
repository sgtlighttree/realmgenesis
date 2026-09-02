// Runs the F4 globe-render instrument (`scripts/perf/globeBench.ts`).
//
// It is a script and not a test on purpose: it generates a 30k world and can
// measure for minutes. `npm test` must never pick it up.
//
// Usage:
//   node scripts/renderGlobeBench.mjs [--points=30000] [--frames=120] [--out=FILE]
//
// Writes a plain-text report and echoes it. Run it once before a change and
// once after, and diff the two files — the per-tenant ms/frame lines are the
// ones that tell you which overlay to attack.
import { createServer } from 'vite';
import { readFileSync } from 'fs';

const flags = Object.fromEntries(
  process.argv.slice(2).filter(a => a.startsWith('--')).map(a => {
    const [k, v] = a.replace(/^--/, '').split('=');
    return [k, v ?? true];
  }),
);
const out = String(flags.out ?? 'tmp/f4bench.txt');
process.env.F4_OUT = out;
process.env.F4_POINTS = String(flags.points ?? 30000);
process.env.F4_FRAMES = String(flags.frames ?? 120);

const server = await createServer({ server: { middlewareMode: true }, appType: 'custom', logLevel: 'error' });
const { run } = await server.ssrLoadModule('/scripts/perf/globeBench.ts');
await run();
await server.close();
console.log(readFileSync(out, 'utf8'));
