// F4 Map2D geometry-cache instrument.
//
// Deliberately NOT a `*.test.ts`: at 200k cells it generates for ~20s and then
// measures. Run it on purpose through `scripts/renderMapCacheBench.mjs`.
//
// Measures the Phase-1 main-thread costs Map2D pays on mount
// (`buildMapGeometryCache`) and on every zoom-settle (`buildBoundaryPaths`),
// broken down by stage so the hot spot is named rather than guessed.
//
// Correctness signal: Path2D is stubbed with a recorder that folds every
// moveTo/lineTo/closePath/addPath into an FNV hash. The hashes of the mount
// cache and of each settle tolerance are printed beside the timings — an
// optimisation that changes any of them changed the drawn geometry.
import { appendFileSync, writeFileSync } from 'node:fs';

import { generateWorld } from '../../utils/worldGen';
import { makeParams } from '../../tests/helpers';
import { WorldData } from '../../types';
import { buildProjection } from '../../utils/projections';
import {
  chainSegments, computeCoastlineSegments, computeFactionBorderSegments,
  computeLakeOutlineSegments,
} from '../../utils/boundaries';

const OUT = process.env.F4_OUT ?? '/tmp/mapCacheBench.txt';
const POINTS = Number(process.env.F4_POINTS ?? 200000);
const REPS = Number(process.env.F4_REPS ?? 3);
const say = (m: string) => { appendFileSync(OUT, m + '\n'); };

// --- Path2D recorder stub ---------------------------------------------------
// One global running hash across every path built, reset per measurement, so a
// whole cache build reduces to a single comparable number.
let H = 0x811c9dc5;
const mix = (v: number) => {
  // Round to 1e-6 px: well below anything visible, well above float noise
  // from reordered-but-equivalent arithmetic.
  const q = Math.round(v * 1e6) | 0;
  H = Math.imul(H ^ q, 0x01000193) >>> 0;
};
class Path2DStub {
  moveTo(x: number, y: number) { mix(1); mix(x); mix(y); }
  lineTo(x: number, y: number) { mix(2); mix(x); mix(y); }
  closePath() { mix(3); }
  addPath() { mix(4); }
  rect() { mix(5); }
  arc() { mix(6); }
}
(globalThis as unknown as { Path2D: unknown }).Path2D = Path2DStub;

const hex = (h: number) => (h >>> 0).toString(16).padStart(8, '0');

const time = <T>(fn: () => T): { ms: number; out: T } => {
  const t0 = performance.now();
  const out = fn();
  return { ms: performance.now() - t0, out };
};

const median = (xs: number[]) => [...xs].sort((a, b) => a - b)[Math.floor(xs.length / 2)];

export async function run(): Promise<void> {
  writeFileSync(OUT, '');
  // Imported after the stub is installed so any module-level Path2D use sees it.
  const { buildMapGeometryCache, buildBoundaryPaths } = await import('../../utils/mapCache');
  const { buildCellQuadtree } = await import('../../utils/mapPick');

  say(`mapCacheBench — ${POINTS} cells, ${REPS} reps (median ms)`);
  const tg = performance.now();
  const world: WorldData = await generateWorld(
    makeParams({ points: POINTS, plates: 12, numFactions: 8, erosionIterations: 3 }),
  );
  say(`generateWorld: ${(performance.now() - tg).toFixed(0)} ms (not part of the cache cost)`);

  const W = 1600;
  const Hh = 800;
  const projection = buildProjection('mercator', W, Hh);
  const n = world.cells.length;

  // Whole mount build.
  const mount: number[] = [];
  let mountHash = 0;
  for (let r = 0; r < REPS; r++) {
    H = 0x811c9dc5;
    mount.push(time(() => buildMapGeometryCache(world, projection, W, Hh)).ms);
    mountHash = H;
  }
  say(`buildMapGeometryCache: ${median(mount).toFixed(1)}   hash=${hex(mountHash)}`);

  // Stage breakdown of the boundary half.
  const segCoast = time(() => computeCoastlineSegments(world));
  const segBorder = time(() => computeFactionBorderSegments(world));
  const segLake = time(() => computeLakeOutlineSegments(world));
  say(`  computeCoastlineSegments:     ${segCoast.ms.toFixed(1)}  (${segCoast.out.length} segs)`);
  say(`  computeFactionBorderSegments: ${segBorder.ms.toFixed(1)}  (${segBorder.out.length} segs)`);
  say(`  computeLakeOutlineSegments:   ${segLake.ms.toFixed(1)}  (${segLake.out.length} segs)`);
  const chCoast = time(() => chainSegments(segCoast.out, n));
  const chBorder = time(() => chainSegments(segBorder.out, n));
  const chLake = time(() => chainSegments(segLake.out, n));
  const longest = (cs: unknown[][]) => cs.reduce((m, c) => Math.max(m, c.length), 0);
  say(`  chainSegments coast:  ${chCoast.ms.toFixed(1)}  (${chCoast.out.length} chains, longest ${longest(chCoast.out)})`);
  say(`  chainSegments border: ${chBorder.ms.toFixed(1)}  (${chBorder.out.length} chains, longest ${longest(chBorder.out)})`);
  say(`  chainSegments lake:   ${chLake.ms.toFixed(1)}  (${chLake.out.length} chains, longest ${longest(chLake.out)})`);

  // Settle LOD rebuilds at the tolerances a zoom actually produces.
  for (const scale of [1, 2, 4, 10]) {
    const tol = Math.max(0.05, 0.5 / Math.max(scale, 1));
    const ms: number[] = [];
    let hsh = 0;
    for (let r = 0; r < REPS; r++) {
      H = 0x811c9dc5;
      ms.push(time(() => buildBoundaryPaths(world, n, projection, tol)).ms);
      hsh = H;
    }
    say(`buildBoundaryPaths scale=${scale} tol=${tol}: ${median(ms).toFixed(1)}   hash=${hex(hsh)}`);
  }

  const qt: number[] = [];
  for (let r = 0; r < REPS; r++) qt.push(time(() => buildCellQuadtree(world, projection, W, Hh)).ms);
  say(`buildCellQuadtree: ${median(qt).toFixed(1)}`);
}
