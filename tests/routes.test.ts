import { describe, it, expect } from 'vitest';
import { generateWorld } from '../utils/worldGen';
import { computeRoutes } from '../utils/routes';
import { isWaterCell } from '../utils/pathfinding';
import { makeParams } from './helpers';
import type { WorldData } from '../types';

// Routes are computed lazily in the app (App.tsx), not by generateWorld, so
// tests compute them explicitly here — exactly what the lazy effect does.
async function build(seed = 'route-test') {
  const world = await generateWorld(makeParams({ seed, civSeed: seed }));
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

// cellId -> { factionId, group key } for every town.
function townIndex(world: WorldData, comp: number[]) {
  const map = new Map<number, { factionId: number; group: string }>();
  (world.civData?.factions ?? []).forEach((f, fi) => {
    for (const prov of f.provinces) {
      for (const t of prov.towns) {
        map.set(t.cellId, { factionId: fi, group: `${fi}:${comp[t.cellId]}` });
      }
    }
  });
  return map;
}

describe('computeRoutes — roads', () => {
  it('is deterministic: same seed → identical routes', async () => {
    const a = await build();
    const b = await build();
    expect(JSON.stringify(a.routes)).toEqual(JSON.stringify(b.routes));
  });

  it('produces at least some road routes with valid town endpoints', async () => {
    const w = await build();
    const roads = (w.routes ?? []).filter(r => r.kind === 'road');
    expect(roads.length).toBeGreaterThan(0);
    for (const r of roads) {
      expect(r.path.length).toBeGreaterThanOrEqual(2);
      expect(r.fromCellId).toBeGreaterThanOrEqual(0);
      expect(r.toCellId).toBeGreaterThanOrEqual(0);
    }
  });

  it('purity: recomputing routes does not mutate terrain or civ fields', async () => {
    const w = await build();
    const before = w.cells.map(c => `${c.height}|${c.biome}|${c.regionId}|${c.provinceId}`).join(',');
    computeRoutes(w, w.params);
    const after = w.cells.map(c => `${c.height}|${c.biome}|${c.regionId}|${c.provinceId}`).join(',');
    expect(after).toEqual(before);
  });

  it('forest invariant: intra-group road edges are acyclic and connect each group', async () => {
    const w = await build();
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
});

describe('computeRoutes — sea', () => {
  it('produces sea routes between distinct coastal towns, structurally valid', async () => {
    const w = await build('sea-test');
    const sea = (w.routes ?? []).filter(r => r.kind === 'searoute');
    expect(sea.length).toBeGreaterThan(0);
    for (const r of sea) {
      expect(r.path.length).toBeGreaterThanOrEqual(2);
      expect(r.fromCellId).not.toEqual(r.toCellId);
    }
  });

  it('sea-route endpoints are coastal cells (a water neighbor exists)', async () => {
    const w = await build('sea-test');
    const seaLevel = w.params.seaLevel;
    const isCoastal = (cellId: number) =>
      w.cells[cellId].neighbors.some(n => isWaterCell(w.cells[n], seaLevel));
    for (const r of (w.routes ?? []).filter(r => r.kind === 'searoute')) {
      expect(isCoastal(r.fromCellId), `from ${r.fromCellId} should be coastal`).toBe(true);
      expect(isCoastal(r.toCellId), `to ${r.toCellId} should be coastal`).toBe(true);
    }
  });
});
