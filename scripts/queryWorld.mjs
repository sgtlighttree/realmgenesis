// Agent-facing world/cell query tool — data-only, no browser, no clicking.
//
// Loads the real TS engine through vite ssrLoadModule (same trick as
// captureBaseline.mjs — no extra deps), generates a world deterministically,
// and answers structured questions about cells, elevations, gradients, biomes.
// Built so an agent can debug terrain/climate without driving the 3D UI, which
// Playwright cannot introspect.
//
// Usage:
//   node scripts/queryWorld.mjs cell <id> [--seed=S] [--points=N]
//   node scripts/queryWorld.mjs hypsometry [--seed=S] [--points=N] [--datum=M]
//   node scripts/queryWorld.mjs gradients [--top=N] [--seed=S] [--points=N]
//   node scripts/queryWorld.mjs biomes [--seed=S] [--points=N]
//   node scripts/queryWorld.mjs near <x> <y> <z> [--seed=S]   (nearest cell to a unit vector)
//
// Flags default: seed=realmgenesis, points=20000, datum=9000.
import { createServer } from 'vite';

const argv = process.argv.slice(2);
const cmd = argv[0];
const positional = argv.slice(1).filter((a) => !a.startsWith('--'));
const flags = Object.fromEntries(
  argv.filter((a) => a.startsWith('--')).map((a) => {
    const [k, v] = a.replace(/^--/, '').split('=');
    return [k, v ?? true];
  }),
);

if (!cmd || cmd === 'help') {
  console.log('commands: cell <id> | hypsometry | gradients | biomes | near <x> <y> <z>');
  console.log('flags: --seed=S --points=N --datum=M --top=N');
  process.exit(cmd ? 0 : 1);
}

const seed = flags.seed ?? 'realmgenesis';
const points = Number(flags.points ?? 20000);
const datum = Number(flags.datum ?? 9000);

const server = await createServer({ server: { middlewareMode: true } });
const { generateWorld, isLakeCell } = await server.ssrLoadModule('/utils/worldGen.ts');
const { makeParams } = await server.ssrLoadModule('/tests/helpers.ts');
const { elevationMetres, formatElevation } = await server.ssrLoadModule('/utils/datum.ts');

const world = await generateWorld(makeParams({ seed, points, maxElevationM: datum }));
const sea = world.params.seaLevel;
const cells = world.cells;
const isLand = (c) => c.height >= sea && !isLakeCell(c);
const elevM = (c) => elevationMetres(c.height, sea, datum);

const factionName = (id) =>
  id === undefined ? null : world.civData?.factions?.find((f) => f.id === id)?.name ?? `#${id}`;

const dumpCell = (c) => ({
  id: c.id,
  biome: c.biome,
  height: Number(c.height.toFixed(4)),
  elevationM: Math.round(elevM(c)),
  elevationLabel: formatElevation(c.height, sea, datum),
  temperatureC: Number(c.temperature.toFixed(1)),
  moisture: Number(c.moisture.toFixed(3)),
  population: c.population ?? 0,
  faction: factionName(c.regionId),
  isLand: isLand(c),
});

function run() {
  if (cmd === 'cell') {
    const id = Number(positional[0]);
    const c = cells[id];
    if (!c) return { error: `no cell ${id} (world has ${cells.length})` };
    return {
      meta: { seed, points: cells.length, seaLevel: sea, datum },
      cell: dumpCell(c),
      neighbors: c.neighbors.map((n) => {
        const nc = cells[n];
        return { ...dumpCell(nc), elevDeltaM: Math.round(elevM(nc) - elevM(c)) };
      }),
    };
  }

  if (cmd === 'hypsometry') {
    // The realism test: is land height spread evenly (linear datum wrong) or
    // concentrated low like Earth (~70% under 1 km)? Reports the land-height
    // distribution AND the metre distribution. --curve=P applies a concave
    // hypsometric power curve elev = datum * frac^P (frac = above-sea fraction),
    // to see if a curve fixes the perception without touching the height field.
    const curve = Number(flags.curve ?? 1);
    const toM = (h) => datum * Math.pow(Math.max(0, (h - sea) / (1 - sea)), curve);
    const land = cells.filter(isLand).map((c) => c.height);
    land.sort((a, b) => a - b);
    const n = land.length;
    const q = (p) => land[Math.floor(p * (n - 1))];
    const bands = [0, 500, 1000, 2000, 3000, 4000, 6000, 9000];
    const metreHist = {};
    for (let i = 0; i < bands.length - 1; i++) {
      const lo = bands[i], hi = bands[i + 1];
      const frac = land.filter((h) => toM(h) >= lo && toM(h) < hi).length / n;
      metreHist[`${lo}-${hi}m`] = `${(frac * 100).toFixed(1)}%`;
    }
    const mean = Math.round(land.reduce((s, h) => s + toM(h), 0) / n);
    // Ocean-side check: is the symmetric curve making abyssal plains too shallow?
    const ocean = cells.filter((c) => c.height < sea && !isLakeCell(c)).map((c) => elevationMetres(c.height, sea, datum));
    ocean.sort((a, b) => a - b); // most negative first
    const on = ocean.length;
    const oq = (p) => Math.round(ocean[Math.floor(p * (on - 1))]);
    const oceanDepth = { cells: on, meanM: Math.round(ocean.reduce((s, m) => s + m, 0) / on), deepest: oq(0), p25: oq(0.25), median: oq(0.5), p75_shallow: oq(0.75), earthRef: 'mean ~-3700m, abyssal plains -3000..-6000m' };
    return {
      oceanDepthActualDatum: oceanDepth,
      meta: { seed, landCells: n, seaLevel: sea, datum, curve },
      normalizedHeightQuantiles: {
        min: q(0).toFixed(3), p10: q(0.1).toFixed(3), p25: q(0.25).toFixed(3),
        median: q(0.5).toFixed(3), p75: q(0.75).toFixed(3), p90: q(0.9).toFixed(3),
        max: q(1).toFixed(3),
      },
      metresWithCurve: {
        curve,
        meanM: mean,
        median: Math.round(toM(q(0.5))),
        p90: Math.round(toM(q(0.9))),
        distribution: metreHist,
      },
      earthReference: '~70% of real land is below 1000m; mean ~840m',
    };
  }

  if (cmd === 'gradients') {
    const top = Number(flags.top ?? 15);
    const pairs = [];
    for (const c of cells) {
      if (!isLand(c)) continue;
      for (const n of c.neighbors) {
        const nc = cells[n];
        if (n <= c.id) continue; // each undirected pair once
        const dM = Math.abs(elevM(nc) - elevM(c));
        pairs.push({ a: c.id, b: nc.id, deltaM: Math.round(dM), aBiome: c.biome, bBiome: nc.biome, aElevM: Math.round(elevM(c)), bElevM: Math.round(elevM(nc)) });
      }
    }
    pairs.sort((x, y) => y.deltaM - x.deltaM);
    return { meta: { seed, points: cells.length, datum }, steepestAdjacentPairs: pairs.slice(0, top) };
  }

  if (cmd === 'biomes') {
    const hist = {};
    let land = 0;
    for (const c of cells) {
      if (!isLand(c)) continue;
      land++;
      hist[c.biome] = (hist[c.biome] ?? 0) + 1;
    }
    const out = {};
    for (const [k, v] of Object.entries(hist).sort((a, b) => b[1] - a[1])) out[k] = `${((v / land) * 100).toFixed(1)}%`;
    return { meta: { seed, landCells: land }, landBiomeShare: out };
  }

  if (cmd === 'climate') {
    // Why does one biome dominate? Reports land moisture + temperature
    // distributions and, for steppe specifically, which determineBiome branch
    // produced it (arid-temperate vs semiarid-cool).
    const land = cells.filter(isLand);
    const mo = land.map((c) => c.moisture).sort((a, b) => a - b);
    const te = land.map((c) => c.temperature).sort((a, b) => a - b);
    const n = land.length;
    const q = (arr, p) => arr[Math.floor(p * (n - 1))];
    const dryLt15 = land.filter((c) => c.moisture < 0.15).length / n;
    const dry15to40 = land.filter((c) => c.moisture >= 0.15 && c.moisture < 0.4).length / n;
    // steppe branches (see determineBiome): A = moisture<0.15 & 10<temp<=25;
    // B = moisture 0.15-0.4 & temp<=10.
    const steppe = land.filter((c) => c.biome === 'Steppe');
    const branchA = steppe.filter((c) => c.moisture < 0.15).length;
    const branchB = steppe.filter((c) => c.moisture >= 0.15 && c.moisture < 0.4).length;
    return {
      meta: { seed, landCells: n },
      moistureQuantiles: { min: q(mo, 0).toFixed(3), p25: q(mo, 0.25).toFixed(3), median: q(mo, 0.5).toFixed(3), p75: q(mo, 0.75).toFixed(3), max: q(mo, 1).toFixed(3) },
      tempQuantilesC: { min: q(te, 0).toFixed(1), p25: q(te, 0.25).toFixed(1), median: q(te, 0.5).toFixed(1), p75: q(te, 0.75).toFixed(1), max: q(te, 1).toFixed(1) },
      dryLandShare: { 'moisture<0.15': `${(dryLt15 * 100).toFixed(1)}%`, 'moisture0.15-0.4': `${(dry15to40 * 100).toFixed(1)}%` },
      steppeBreakdown: { total: steppe.length, aridTemperate_moLt15: branchA, semiaridCool_mo15to40: branchB },
      earthReference: 'arid+semiarid ~33% of real land',
    };
  }

  if (cmd === 'near') {
    const [x, y, z] = positional.map(Number);
    let best = null, bestD = Infinity;
    for (const c of cells) {
      const d = (c.center.x - x) ** 2 + (c.center.y - y) ** 2 + (c.center.z - z) ** 2;
      if (d < bestD) { bestD = d; best = c; }
    }
    return { query: { x, y, z }, nearest: dumpCell(best) };
  }

  return { error: `unknown command "${cmd}"` };
}

console.log(JSON.stringify(run(), null, 2));
await server.close();
process.exit(0);
