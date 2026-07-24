import { isWaterCell, landTerrainStepCost, aStar } from './pathfinding';
import type { Cell, Point, RouteData, WorldData, WorldParams } from '../types';

const MAX_ROAD_EXPANSIONS = 20000; // per-edge A* cap
const MAX_ROAD_ANGLE = 1.2;        // skip MST edges longer than this (radians)

function angle(a: Point, b: Point): number {
  const d = a.x * b.x + a.y * b.y + a.z * b.z;
  return Math.acos(Math.max(-1, Math.min(1, d)));
}

// BFS-label every land cell with a component id (-1 for water).
function landComponentIds(world: WorldData, params: WorldParams): number[] {
  const comp = new Array<number>(world.cells.length).fill(-1);
  let next = 0;
  for (const c of world.cells) {
    if (comp[c.id] !== -1 || isWaterCell(c, params.seaLevel)) continue;
    const id = next++;
    const stack = [c.id];
    comp[c.id] = id;
    while (stack.length) {
      const cur = stack.pop()!;
      for (const nId of world.cells[cur].neighbors) {
        const n = world.cells[nId];
        if (comp[nId] === -1 && !isWaterCell(n, params.seaLevel)) { comp[nId] = id; stack.push(nId); }
      }
    }
  }
  return comp;
}

interface Town { cellId: number; factionId: number; population: number; isCapital: boolean; }

function gatherTowns(world: WorldData): Town[] {
  const towns: Town[] = [];
  (world.civData?.factions ?? []).forEach((f, fi) => {
    for (const prov of f.provinces) {
      for (const t of prov.towns) {
        towns.push({ cellId: t.cellId, factionId: fi, population: t.population, isCapital: t.isCapital });
      }
    }
  });
  return towns;
}

// Kruskal MST over the given town indices, weighted by great-circle distance.
// Deterministic: edges sorted by (weight, minCellId, maxCellId).
function mstEdges(towns: Town[], idx: number[], cells: Cell[]): Array<[number, number]> {
  const edges: Array<{ a: number; b: number; w: number }> = [];
  for (let i = 0; i < idx.length; i++) {
    for (let j = i + 1; j < idx.length; j++) {
      const a = towns[idx[i]].cellId, b = towns[idx[j]].cellId;
      edges.push({ a, b, w: angle(cells[a].center, cells[b].center) });
    }
  }
  edges.sort((e1, e2) =>
    e1.w - e2.w ||
    Math.min(e1.a, e1.b) - Math.min(e2.a, e2.b) ||
    Math.max(e1.a, e1.b) - Math.max(e2.a, e2.b));
  const parent = new Map<number, number>();
  const find = (x: number): number => { while (parent.get(x)! !== x) { parent.set(x, parent.get(parent.get(x)!)!); x = parent.get(x)!; } return x; };
  for (const e of edges) { parent.set(e.a, e.a); parent.set(e.b, e.b); }
  const out: Array<[number, number]> = [];
  for (const e of edges) {
    const ra = find(e.a), rb = find(e.b);
    if (ra !== rb) { parent.set(ra, rb); out.push([e.a, e.b]); }
  }
  return out;
}

function landHeuristic(cells: Cell[], goalId: number, stepsPerRadian: number): (id: number) => number {
  const goal = cells[goalId].center;
  return (id: number) => angle(cells[id].center, goal) * stepsPerRadian;
}

export function computeRoutes(world: WorldData, params: WorldParams): RouteData[] {
  const cells = world.cells;
  const comp = landComponentIds(world, params);
  const towns = gatherTowns(world);
  const routes: RouteData[] = [];

  // Rough lower bound on cells crossed per radian, for the A* heuristic.
  const stepsPerRadian = Math.max(1, Math.sqrt(cells.length / (4 * Math.PI)));

  const landStep = (fromId: number, toId: number): number =>
    isWaterCell(cells[toId], params.seaLevel) ? Infinity : landTerrainStepCost(cells[fromId], cells[toId]);

  // Roads: one MST per (faction, land-component) group.
  const groups = new Map<string, number[]>();
  towns.forEach((t, i) => {
    if (comp[t.cellId] < 0) return; // town on water (shouldn't happen) — skip
    const key = `${t.factionId}:${comp[t.cellId]}`;
    if (!groups.has(key)) groups.set(key, []);
    groups.get(key)!.push(i);
  });

  for (const idx of groups.values()) {
    if (idx.length < 2) continue;
    for (const [a, b] of mstEdges(towns, idx, cells)) {
      if (angle(cells[a].center, cells[b].center) > MAX_ROAD_ANGLE) continue;
      const path = aStar(cells, a, b, landStep, landHeuristic(cells, b, stepsPerRadian), MAX_ROAD_EXPANSIONS);
      if (path && path.length >= 2) {
        routes.push({ path: path.map(id => cells[id].center), kind: 'road', fromCellId: a, toCellId: b });
      }
    }
  }

  // Capital trunk roads between bordering factions in the same land component.
  const capitals = towns.filter(t => t.isCapital);
  for (let i = 0; i < capitals.length; i++) {
    for (let j = i + 1; j < capitals.length; j++) {
      const A = capitals[i], B = capitals[j];
      if (A.factionId === B.factionId) continue;
      if (comp[A.cellId] < 0 || comp[A.cellId] !== comp[B.cellId]) continue;
      if (angle(cells[A.cellId].center, cells[B.cellId].center) > MAX_ROAD_ANGLE) continue;
      const path = aStar(cells, A.cellId, B.cellId, landStep, landHeuristic(cells, B.cellId, stepsPerRadian), MAX_ROAD_EXPANSIONS);
      if (path && path.length >= 2) {
        routes.push({ path: path.map(id => cells[id].center), kind: 'road', fromCellId: A.cellId, toCellId: B.cellId });
      }
    }
  }

  return routes;
}
