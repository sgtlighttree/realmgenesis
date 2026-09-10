import { describe, it, expect, beforeAll } from 'vitest';
import { generateWorld } from '../utils/worldGen';
import { computeRoutes, ARTERIAL_TIERS } from '../utils/routes';
import { isWaterCell } from '../utils/pathfinding';
import { makeParams } from './helpers';
import type { Settlement, WorldData, WorldParams } from '../types';

// Routes are computed lazily in the app (App.tsx), not by generateWorld, so
// tests compute them explicitly here — exactly what the lazy effect does.
async function build(seed = 'route-test', overrides: Partial<WorldParams> = {}) {
  const world = await generateWorld(makeParams({ seed, civSeed: seed, ...overrides }));
  world.routes = computeRoutes(world, world.params);
  return world;
}

// Reconstruct land-connected components the same way computeRoutes does, so the
// forest invariant can be checked without exposing routes.ts internals.
function landComponents(world: WorldData): number[] {
  const seaLevel = world.params.seaLevel;
  const comp = new Array<number>(world.cells.length).fill(-1);
  let next = 0;
  for (const c of world.cells) {
    if (comp[c.id] !== -1 || isWaterCell(c, seaLevel)) continue;
    const id = next++;
    const stack = [c.id];
    comp[c.id] = id;
    while (stack.length) {
      const cur = stack.pop()!;
      for (const nId of world.cells[cur].neighbors) {
        if (comp[nId] === -1 && !isWaterCell(world.cells[nId], seaLevel)) { comp[nId] = id; stack.push(nId); }
      }
    }
  }
  return comp;
}

// cellId -> { factionId, group key } for every arterial-tier settlement
// (capital|city|town) — mirrors computeRoutes' gatherTowns() gating, since
// arterial MST/trunk edges only ever join that subset.
function townIndex(world: WorldData, comp: number[]) {
  const map = new Map<number, { factionId: number; group: string }>();
  (world.civData?.settlements ?? []).forEach(s => {
    if (!ARTERIAL_TIERS.has(s.tier)) return;
    map.set(s.cellId, { factionId: s.factionId, group: `${s.factionId}:${comp[s.cellId]}` });
  });
  return map;
}

// cellId -> Settlement, for tier lookups in the arterial-gating tests.
function settlementByCell(world: WorldData): Map<number, Settlement> {
  const map = new Map<number, Settlement>();
  (world.civData?.settlements ?? []).forEach(s => map.set(s.cellId, s));
  return map;
}

// Scale + density overrides so the test map actually produces every settlement
// tier. Density is bumped above the DEFAULT_PARAMS 1.0 because a faction's quota
// is `round(landCells * density * 0.02)`; `points` is bumped off makeParams()'s
// 300 for two reasons (S36):
//   1. At 300 cells a faction can hold only a capital plus one other settlement,
//      so there are no village/hamlet feeders for Task 8's pass to feed.
//   2. 300 cells span the sphere so coarsely that a single faction's two
//      settlements sit 1.26 rad apart — past routes.ts's MAX_ROAD_ANGLE of 1.2 —
//      so computeRoutes deliberately SKIPS the MST edge and the group is no
//      longer a spanning tree. That is the cap working as designed, not a bug;
//      the fix is a test map at a scale the cap was tuned for.
// Measured at 1000/3.0: 6 capital / 4 city / 3 town / 4 village / 2 hamlet,
// 8 arterial / 5 local / 9 sea, zero spanning-tree violations. generateWorld
// costs the same as at 300 (~2.5s — fixed pipeline overhead dominates).
// Kept as a named override rather than changing makeParams()'s baseline, which
// other suites depend on staying untouched.
const DENSE_OVERRIDES: Partial<WorldParams> = { points: 1000, settlementDensity: 3.0 };

describe('computeRoutes — roads', () => {
  // Single shared build for every non-determinism assertion in this describe
  // block — keeps world-gen count from growing as tier-gating/feeder
  // assertions are added (routes.test.ts already ran 7 full generateWorld
  // calls before this task; see CLAUDE.md's per-file CPU budget note).
  let w: WorldData;
  beforeAll(async () => {
    w = await build('route-test', DENSE_OVERRIDES);
  });

  it('is deterministic: same seed → identical routes', async () => {
    const b = await build('route-test', DENSE_OVERRIDES);
    expect(JSON.stringify(w.routes)).toEqual(JSON.stringify(b.routes));
  });

  it('produces at least some road routes with valid town endpoints', () => {
    const roads = (w.routes ?? []).filter(r => r.kind === 'road');
    expect(roads.length).toBeGreaterThan(0);
    for (const r of roads) {
      expect(r.path.length).toBeGreaterThanOrEqual(2);
      expect(r.fromCellId).toBeGreaterThanOrEqual(0);
      expect(r.toCellId).toBeGreaterThanOrEqual(0);
    }
  });

  it('purity: recomputing routes does not mutate terrain or civ fields', () => {
    const before = w.cells.map(c => `${c.height}|${c.biome}|${c.regionId}|${c.provinceId}`).join(',');
    computeRoutes(w, w.params);
    const after = w.cells.map(c => `${c.height}|${c.biome}|${c.regionId}|${c.provinceId}`).join(',');
    expect(after).toEqual(before);
  });

  it('forest invariant: intra-group road edges are acyclic and connect each group', () => {
    const comp = landComponents(w);
    const towns = townIndex(w, comp);

    // Group sizes (towns per (faction, land-component)).
    const groupSize = new Map<string, number>();
    for (const { group } of towns.values()) groupSize.set(group, (groupSize.get(group) ?? 0) + 1);

    // Union-find over intra-group road edges (same-faction endpoints => MST edge,
    // not a cross-faction capital trunk).
    const parent = new Map<number, number>();
    const find = (x: number): number => { while (parent.get(x)! !== x) { parent.set(x, parent.get(parent.get(x)!)!); x = parent.get(x)!; } return x; };
    for (const cellId of towns.keys()) parent.set(cellId, cellId);
    const intraEdges = new Map<string, number>();

    for (const r of (w.routes ?? [])) {
      if (r.kind !== 'road') continue;
      const a = towns.get(r.fromCellId), b = towns.get(r.toCellId);
      if (!a || !b || a.group !== b.group) continue; // trunk / non-town edge
      const ra = find(r.fromCellId), rb = find(r.toCellId);
      expect(ra, 'road MST subgraph must be acyclic (no edge unites already-connected towns)').not.toEqual(rb);
      parent.set(ra, rb);
      intraEdges.set(a.group, (intraEdges.get(a.group) ?? 0) + 1);
    }

    // Each multi-town group is a spanning tree: edges == size - 1.
    for (const [group, size] of groupSize) {
      if (size < 2) continue;
      expect(intraEdges.get(group) ?? 0, `group ${group} (size ${size}) should be a spanning tree`).toEqual(size - 1);
    }
  });

  it('gates arterial roads to top settlement tiers: no village/hamlet cell is an arterial road endpoint', () => {
    const byCell = settlementByCell(w);
    let checked = 0;
    for (const r of (w.routes ?? [])) {
      if (r.kind !== 'road' || r.class !== 'arterial') continue;
      for (const cellId of [r.fromCellId, r.toCellId]) {
        const s = byCell.get(cellId);
        expect(s, `road endpoint ${cellId} should be a settlement`).toBeDefined();
        expect(ARTERIAL_TIERS.has(s!.tier), `road endpoint ${cellId} has tier ${s!.tier}`).toBe(true);
        checked++;
      }
    }
    expect(checked).toBeGreaterThan(0);
  });

  it('every emitted route has a valid class field', () => {
    let checked = 0;
    for (const r of (w.routes ?? [])) {
      expect(['arterial', 'local']).toContain(r.class);
      checked++;
    }
    expect(checked).toBeGreaterThan(0);
  });

  it('every MST/trunk/sea route is tagged class "arterial"', () => {
    let checked = 0;
    const byCell = settlementByCell(w);
    for (const r of (w.routes ?? [])) {
      if (r.kind === 'searoute') { expect(r.class).toEqual('arterial'); checked++; continue; }
      // Roads with a top-tier settlement at BOTH ends are arterial (MST/trunk);
      // roads with a village/hamlet at either end are local feeders.
      const from = byCell.get(r.fromCellId), to = byCell.get(r.toCellId);
      if (from && to && ARTERIAL_TIERS.has(from.tier) && ARTERIAL_TIERS.has(to.tier)) {
        expect(r.class).toEqual('arterial');
        checked++;
      }
    }
    expect(checked).toBeGreaterThan(0);
  });

  it('produces capped local feeder roads from villages/hamlets to their nearest hub', () => {
    const byCell = settlementByCell(w);
    const locals = (w.routes ?? []).filter(r => r.class === 'local');
    expect(locals.length).toBeGreaterThan(0);
    let checked = 0;
    for (const r of locals) {
      expect(r.kind).toEqual('road');
      const from = byCell.get(r.fromCellId);
      const to = byCell.get(r.toCellId);
      expect(from, `local route "from" ${r.fromCellId} should be a settlement`).toBeDefined();
      expect(to, `local route "to" ${r.toCellId} should be a settlement`).toBeDefined();
      expect(ARTERIAL_TIERS.has(from!.tier), `local route "from" ${r.fromCellId} should be village/hamlet`).toBe(false);
      expect(ARTERIAL_TIERS.has(to!.tier), `local route "to" ${r.toCellId} should be a hub tier`).toBe(true);
      expect(from!.factionId).toEqual(to!.factionId);
      checked++;
    }
    expect(checked).toBeGreaterThan(0);
  });
});

// Rescanned for D7 Task 1c (plateJitter/boundaryRoughness defaults raised to
// 1.5): 'sea-test' no longer produces distinct-landmass coastal towns under
// the new terrain (0 sea routes at the new defaults) — 'islands' does.
const SEA_SEED = 'islands';

describe('computeRoutes — sea', () => {
  it('produces sea routes between distinct coastal towns, structurally valid', async () => {
    const w = await build(SEA_SEED);
    const sea = (w.routes ?? []).filter(r => r.kind === 'searoute');
    expect(sea.length).toBeGreaterThan(0);
    for (const r of sea) {
      expect(r.path.length).toBeGreaterThanOrEqual(2);
      expect(r.fromCellId).not.toEqual(r.toCellId);
    }
  });

  it('sea-route endpoints are coastal cells (a water neighbor exists)', async () => {
    const w = await build(SEA_SEED);
    const seaLevel = w.params.seaLevel;
    const isCoastal = (cellId: number) =>
      w.cells[cellId].neighbors.some(n => isWaterCell(w.cells[n], seaLevel));
    for (const r of (w.routes ?? []).filter(r => r.kind === 'searoute')) {
      expect(isCoastal(r.fromCellId), `from ${r.fromCellId} should be coastal`).toBe(true);
      expect(isCoastal(r.toCellId), `to ${r.toCellId} should be coastal`).toBe(true);
    }
  });
});
