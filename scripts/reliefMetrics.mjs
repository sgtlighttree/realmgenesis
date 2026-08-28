// Agent-facing terrain-relief metrics — the instrument that root-caused D10.
//
// An aggregate depth histogram CANNOT tell "smooth swells" from "bumpy": the sea
// bed already spanned ~8,000 m against land's ~2,700 m while still reading flat.
// This tool reports relief at TWO scales, so the two are separable:
//
//   s1 = mean |Δheight| to immediate neighbours        → TEXTURE
//   s2 = mean |Δheight| to neighbours-of-neighbours    → SWELL
//
// Both in normalized height units (what the renderer displaces), split by ocean
// vs land. The land/ocean s1 ratio is the headline number: 1.0 means the sea bed
// carries the same texture as land, higher means the sea bed is comparatively
// smooth. It also reports land-cell count and coastal-ocean depth, because a
// change that improves ocean texture by resizing the continents or by welding a
// deposition rim to every coastline is not an improvement.
//
// Companion to scripts/queryWorld.mjs. Loads the real engine through vite
// ssrLoadModule — no extra deps, no browser.
//
// Usage:
//   node scripts/reliefMetrics.mjs [--seeds=a,b,c] [--points=N] [--param=value ...]
//
// Any WorldParams key can be overridden as a numeric flag, e.g.:
//   node scripts/reliefMetrics.mjs --seafloorRelief=2.0
//   node scripts/reliefMetrics.mjs --points=80000 --erosionIterations=0
//
// Defaults: seeds=realmgenesis,s2,s17  points=20000
import { createServer } from 'vite';

const flags = Object.fromEntries(
  process.argv.slice(2).filter(a => a.startsWith('--')).map(a => {
    const [k, v] = a.replace(/^--/, '').split('=');
    return [k, v ?? true];
  }),
);

const seeds = String(flags.seeds ?? 'realmgenesis,s2,s17').split(',');
const points = Number(flags.points ?? 20000);
const overrides = {};
for (const [k, v] of Object.entries(flags)) {
  if (k === 'seeds' || k === 'points') continue;
  const n = Number(v);
  overrides[k] = Number.isFinite(n) ? n : v;
}

const server = await createServer({ server: { middlewareMode: true } });
const { generateWorld, isLakeCell } = await server.ssrLoadModule('/utils/worldGen.ts');
const { makeParams } = await server.ssrLoadModule('/tests/helpers.ts');

const mean = a => a.reduce((x, y) => x + y, 0) / (a.length || 1);
const pct = (a, p) => { const s = [...a].sort((x, y) => x - y); return s[Math.floor(p * (s.length - 1))]; };

// Mean |Δheight| from each cell to the ring `hops` steps away, staying inside
// `set` so ocean relief never picks up the land/ocean step at the coastline.
const relief = (cells, set, hops) => {
  const ids = new Set(set.map(c => c.id));
  const out = [];
  for (const c of set) {
    let seen = new Set([c.id]);
    let frontier = [c.id];
    for (let h = 0; h < hops; h++) {
      const next = [];
      for (const id of frontier) {
        for (const n of cells[id].neighbors ?? []) {
          if (ids.has(n) && !seen.has(n)) { seen.add(n); next.push(n); }
        }
      }
      frontier = next;
    }
    if (frontier.length) out.push(mean(frontier.map(n => Math.abs(cells[n].height - c.height))));
  }
  return out;
};

const agg = { o1: [], o2: [], l1: [], l2: [], land: [], rim: [] };

console.log(`points=${points}  overrides=${JSON.stringify(overrides)}`);
for (const seed of seeds) {
  const world = await generateWorld(makeParams({ seed, points, ...overrides }));
  const sea = world.params.seaLevel;
  const cells = world.cells;
  const ocean = cells.filter(c => c.height < sea && !isLakeCell(c));
  const land = cells.filter(c => c.height >= sea && !isLakeCell(c));
  const coastal = ocean.filter(c => (c.neighbors ?? []).some(n => cells[n].height >= sea));

  const o1 = mean(relief(cells, ocean, 1));
  const o2 = mean(relief(cells, ocean, 2));
  const l1 = mean(relief(cells, land, 1));
  const l2 = mean(relief(cells, land, 2));
  const rim = pct(coastal.map(c => sea - c.height), 0.05);

  agg.o1.push(o1); agg.o2.push(o2); agg.l1.push(l1); agg.l2.push(l2);
  agg.land.push(land.length); agg.rim.push(rim);

  console.log(
    `  ${seed.padEnd(14)} ocean s1 ${o1.toFixed(4)} s2 ${o2.toFixed(4)} | ` +
    `land s1 ${l1.toFixed(4)} s2 ${l2.toFixed(4)} | ratio ${(l1 / o1).toFixed(2)}x | ` +
    `land ${land.length}`,
  );
}

const O1 = mean(agg.o1), L1 = mean(agg.l1);
console.log('');
console.log(`AVG  ocean s1 ${O1.toFixed(4)} s2 ${mean(agg.o2).toFixed(4)} | land s1 ${L1.toFixed(4)} s2 ${mean(agg.l2).toFixed(4)}`);
console.log(`AVG  TEXTURE RATIO land/ocean ${(L1 / O1).toFixed(2)}x   (1.0 = sea bed as textured as land)`);
console.log(`AVG  land cells ${Math.round(mean(agg.land))}   coastal-ocean p5 depth ${mean(agg.rim).toFixed(4)}`);
console.log(`     ocean s2/s1 ${(mean(agg.o2) / O1).toFixed(2)}  land s2/s1 ${(mean(agg.l2) / L1).toFixed(2)}   (high = swells, low = texture)`);

await server.close();
process.exit(0);
