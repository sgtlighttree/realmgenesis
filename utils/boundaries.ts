import { BiomeType, Cell, LakeData, Point, WorldData } from '../types';
import { Point3 } from './geo';

const toPoint3 = (p: Point): Point3 => [p.x, p.y, p.z];

/**
 * Shared Voronoi edge geometry — the boundary segments between neighbouring
 * cells that differ by some predicate.
 *
 * Extracted from `exportVector.ts` (A3 Task 8) so the SVG export, the PNG
 * export and Map2D all draw the SAME coastline. It lives in its own module
 * rather than in `shading.ts`: this scan is generic boundary geometry serving
 * coastlines and faction borders alike, while shading.ts is about light and
 * contours.
 *
 * NOTE: this duplicates Map2D's `getFactionBorders` adjacency scan on purpose —
 * `utils/` must never import from `components/`.
 */

export const isLandCell = (cell: Cell, seaLevel: number): boolean =>
  cell.height >= seaLevel && cell.biome !== BiomeType.LAKE && cell.biome !== BiomeType.SALT_LAKE;

// Shared Voronoi edges between neighbor cells whose `differs` predicate is
// true — the same shared-vertex adjacency scan as Map2D's getFactionBorders
// (copied rather than imported: utils/ must never import from components/),
// generalized so both coastlines and faction borders can reuse it.
export const computeBoundarySegments = (
  world: WorldData,
  differs: (a: Cell, b: Cell) => boolean,
): Array<[Point3, Point3]> => {
  const segments: Array<[Point3, Point3]> = [];
  const threshold = 0.000001;
  world.cells.forEach((cellA) => {
    cellA.neighbors.forEach((nId) => {
      const cellB = world.cells[nId];
      if (!cellB || cellA.id >= cellB.id) return;
      if (!differs(cellA, cellB)) return;

      const shared: Point3[] = [];
      for (const vA of cellA.vertices) {
        for (const vB of cellB.vertices) {
          const distSq = (vA.x - vB.x) ** 2 + (vA.y - vB.y) ** 2 + (vA.z - vB.z) ** 2;
          if (distSq < threshold) {
            shared.push(toPoint3(vA));
            break;
          }
        }
        if (shared.length === 2) break;
      }
      if (shared.length === 2) segments.push([shared[0], shared[1]]);
    });
  });
  return segments;
};

export const computeCoastlineSegments = (world: WorldData): Array<[Point3, Point3]> => {
  const seaLevel = world.params.seaLevel;
  return computeBoundarySegments(world, (a, b) => isLandCell(a, seaLevel) !== isLandCell(b, seaLevel));
};

export const computeFactionBorderSegments = (world: WorldData): Array<[Point3, Point3]> => {
  if (!world.civData) return [];
  return computeBoundarySegments(
    world,
    (a, b) => a.regionId !== b.regionId && !(a.regionId === undefined && b.regionId === undefined),
  );
};

// Endpoints come from a distance-threshold vertex match (see
// computeBoundarySegments: shared vertices are found within sqrt(1e-6) ~= 1e-3
// of each other), and the *stored* float for "the same" conceptual vertex can
// differ depending on which cell's vertex array contributed it — a T-junction
// vertex shared by 3+ cells uses cellA's copy for two of its edges but a
// different cell's copy for the third (see computeBoundarySegments: it always
// pushes `vA`, the first cell's own vertex, never `vB`). Observed drift on a
// generated world tops out around 4e-3, two orders of magnitude below the
// shortest real Voronoi edge (~0.05 on a unit sphere at typical point counts),
// so a single fixed rounding grid is the wrong tool here: whatever grid size
// is chosen, points can land on opposite sides of a cell boundary purely by
// chance and fail to merge (confirmed empirically — coarsening the grid from
// 1e-3 to 5e-3 made MORE endpoints split, not fewer, because the exact same
// point pairs straddled the newly picked bucket lines differently).
// Instead, weld endpoints by an actual neighbor-search: bucket points in a
// grid sized to the weld threshold, then when placing a new point, look for an
// existing representative in any of the 27 neighboring cells (covers every
// cell the point could be adjacent to, regardless of where it falls within its
// own bucket) within WELD_THRESHOLD. This merges "same vertex, different
// floats" reliably without any risk of merging two genuinely distinct nearby
// vertices (threshold is far below true minimum vertex spacing).
const WELD_THRESHOLD = 0.01;
const WELD_CELL = WELD_THRESHOLD;

const cellIndex = (v: number): number => Math.floor(v / WELD_CELL);

class VertexWelder {
  private grid = new Map<string, number[]>(); // cell key -> representative ids
  private reps: Point3[] = [];

  private cellKey(ix: number, iy: number, iz: number): string {
    return `${ix},${iy},${iz}`;
  }

  /** Returns a stable integer id for `p`, merging it with any existing
   * representative within WELD_THRESHOLD. */
  idFor(p: Point3): number {
    const ix = cellIndex(p[0]);
    const iy = cellIndex(p[1]);
    const iz = cellIndex(p[2]);
    for (let dx = -1; dx <= 1; dx++) {
      for (let dy = -1; dy <= 1; dy++) {
        for (let dz = -1; dz <= 1; dz++) {
          const cands = this.grid.get(this.cellKey(ix + dx, iy + dy, iz + dz));
          if (!cands) continue;
          for (const id of cands) {
            const r = this.reps[id];
            const distSq = (r[0] - p[0]) ** 2 + (r[1] - p[1]) ** 2 + (r[2] - p[2]) ** 2;
            if (distSq < WELD_THRESHOLD * WELD_THRESHOLD) return id;
          }
        }
      }
    }
    const id = this.reps.length;
    this.reps.push(p);
    const key = this.cellKey(ix, iy, iz);
    const arr = this.grid.get(key);
    if (arr) arr.push(id); else this.grid.set(key, [id]);
    return id;
  }

  repOf(id: number): Point3 {
    return this.reps[id];
  }
}

/**
 * Assemble disjoint 2-point boundary edges into contiguous polylines. A chain
 * that returns to its start is closed (last point === first). At a branch point
 * (endpoint shared by >2 edges) the walk continues deterministically (lowest
 * unused edge index) and the remaining edges seed their own chains, so no edge
 * is ever dropped. Edge count is conserved: sum(chain.length - 1) === segments.length.
 */
export const chainSegments = (segments: Array<[Point3, Point3]>): Point3[][] => {
  const welder = new VertexWelder();
  const ids: Array<[number, number]> = segments.map(([a, b]) => [welder.idFor(a), welder.idFor(b)]);

  // welded id -> list of segment indices touching it
  const byId = new Map<number, number[]>();
  const push = (id: number, i: number) => {
    const arr = byId.get(id);
    if (arr) arr.push(i); else byId.set(id, [i]);
  };
  ids.forEach(([a, b], i) => { push(a, i); push(b, i); });

  const used = new Uint8Array(segments.length);
  const chains: Point3[][] = [];

  const nextFrom = (id: number): number => {
    const cands = byId.get(id);
    if (!cands) return -1;
    for (const i of cands) if (!used[i]) return i; // lowest unused index
    return -1;
  };

  for (let start = 0; start < segments.length; start++) {
    if (used[start]) continue;
    used[start] = 1;
    const [a, b] = segments[start];
    const [aId, bId] = ids[start];
    const chain: Point3[] = [a, b];
    const chainIds: number[] = [aId, bId];
    let closed = false;
    // Extend forward from the tail.
    let tailId = bId;
    for (;;) {
      const ni = nextFrom(tailId);
      if (ni === -1) break;
      used[ni] = 1;
      const [naId, nbId] = ids[ni];
      const [na, nb] = segments[ni];
      const nextId = naId === tailId ? nbId : naId;
      const next = naId === tailId ? nb : na;
      chain.push(next);
      chainIds.push(nextId);
      tailId = nextId;
      if (tailId === chainIds[0]) { closed = true; break; } // closed ring
    }
    // Segments can arrive in any order (e.g. reversed / out-of-order relative
    // to the walk direction), so also extend backward from the head — without
    // this a segment placed "before" its true predecessor in the input array
    // never gets attached and the chain fragments.
    if (!closed) {
      let headId = chainIds[0];
      for (;;) {
        const pi = nextFrom(headId);
        if (pi === -1) break;
        used[pi] = 1;
        const [paId, pbId] = ids[pi];
        const [pa, pb] = segments[pi];
        const prevId = paId === headId ? pbId : paId;
        const prev = paId === headId ? pb : pa;
        chain.unshift(prev);
        chainIds.unshift(prevId);
        headId = prevId;
        if (headId === chainIds[chainIds.length - 1]) break; // closed ring
      }
    }
    // Snap the two end representatives together when the chain closed, so
    // consumers checking `chain[0] === chain[last]` see exact equality rather
    // than the sub-threshold float drift the welder absorbed.
    if (chainIds[0] === chainIds[chainIds.length - 1]) {
      chain[chain.length - 1] = welder.repOf(chainIds[0]);
      chain[0] = welder.repOf(chainIds[0]);
    }
    chains.push(chain);
  }
  return chains;
};

/**
 * Boundary edges around lake cells — a lake cell adjacent to any cell not in the
 * SAME lake. Lakes carry only member cell ids (LakeData.cellIds); this derives
 * their outline for the vector map. Reuses the shared adjacency scan.
 */
export const computeLakeOutlineSegments = (world: WorldData): Array<[Point3, Point3]> => {
  const lakes = world.lakes;
  if (!lakes || lakes.length === 0) return [];
  const lakeOf = new Int32Array(world.cells.length).fill(-1);
  lakes.forEach((lake: LakeData, li: number) => {
    for (const id of lake.cellIds) lakeOf[id] = li;
  });
  return computeBoundarySegments(world, (a, b) => lakeOf[a.id] !== lakeOf[b.id]
    && (lakeOf[a.id] !== -1 || lakeOf[b.id] !== -1));
};
