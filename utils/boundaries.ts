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
  const threshold = 0.000001; // sq-distance match threshold for "same vertex" between cellA/cellB (~1e-3 raw)
  // A shared Voronoi edge has exactly 2 DISTINCT endpoints. A near-degenerate
  // cell can carry two vertices in its own `vertices` array that are almost
  // coincident (observed drift ~7e-4 raw distance on a 3000-point generated
  // world). Without a distinctness check, the scan below can match against
  // BOTH near-duplicate points and accept them as "the two shared vertices" —
  // emitting a spurious near-zero-length segment while never reaching the
  // true, genuinely distinct second shared vertex. That both fabricates a
  // bogus edge and drops the real one, corrupting the boundary graph
  // (confirmed on a generated world: one endpoint of the drop went from
  // degree 2 to degree 3, the other from degree 2 to degree 1 — exactly the
  // two odd-degree vertices you'd get from swapping a real edge for a
  // duplicate self-adjacent one).
  //
  // DISTINCT_SQ rejects a second candidate within this squared-distance of
  // the first accepted match and keeps scanning cellA's vertices for the
  // real second vertex instead. It must scale with point density, the same
  // way weldThresholdFor (below) scales the chaining weld threshold: unit
  // -sphere spacing shrinks as ~sqrt(4*pi/N), and average real edge length
  // shrinks right along with it — a FIXED epsilon that's safe at 3000 points
  // thins out at higher N (a fixed ~2e-3 raw epsilon sits ~33x below the
  // *average* spacing at 3000 points but only ~4x below it at 200k, eating
  // into the margin against genuinely short — not just average-length — real
  // edges, which do occur off the average). DISTINCT_FRACTION * spacing(N)
  // keeps the same ~33x margin below average spacing at every N (this is a
  // margin against the *average*, not a hard worst-case bound: an
  // unusually short real edge below the epsilon could in principle still be
  // misclassified, same caveat as the weld threshold below), while comfortably
  // clearing the ~7e-4 duplicate-vertex artifact, assumed to scale down with
  // mesh density along with everything else geometric (0.031 chosen to
  // reproduce the previously-verified fixed value, ~2e-3, at 3000 points):
  //   N=3000    -> spacing ~0.0647, epsilon ~0.0020 (~33x below avg spacing, ~2.9x the ~7e-4 drift)
  //   N=30000   -> spacing ~0.0205, epsilon ~0.00064 (~33x below avg spacing)
  //   N=200000  -> spacing ~0.0079, epsilon ~0.00024 (~33x below avg spacing)
  const DISTINCT_FRACTION = 0.031;
  const distinctEps = DISTINCT_FRACTION * Math.sqrt((4 * Math.PI) / Math.max(world.cells.length, 1));
  const DISTINCT_SQ = distinctEps * distinctEps;
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
            const p = toPoint3(vA);
            if (shared.length === 1) {
              const [fx, fy, fz] = shared[0];
              const distToFirstSq = (p[0] - fx) ** 2 + (p[1] - fy) ** 2 + (p[2] - fz) ** 2;
              // Near-duplicate of the already-accepted vertex — not a
              // genuinely distinct second endpoint. Keep scanning cellA's
              // remaining vertices instead of accepting this one.
              if (distToFirstSq < DISTINCT_SQ) break;
            }
            shared.push(p);
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
// generated world tops out around 7e-4, orders of magnitude below the
// shortest real Voronoi edge, so a single fixed rounding grid is the wrong
// tool here: whatever grid size is chosen, points can land on opposite sides
// of a cell boundary purely by chance and fail to merge (confirmed
// empirically — coarsening the grid from 1e-3 to 5e-3 made MORE endpoints
// split, not fewer, because the exact same point pairs straddled the newly
// picked bucket lines differently).
// Instead, weld endpoints by an actual neighbor-search: bucket points in a
// grid sized to the weld threshold, then when placing a new point, look for an
// existing representative in any of the 27 neighboring cells (covers every
// cell the point could be adjacent to, regardless of where it falls within its
// own bucket) within the threshold. This merges "same vertex, different
// floats" reliably without any risk of merging two genuinely distinct nearby
// vertices, AS LONG AS the threshold stays well below true minimum vertex
// spacing — which shrinks as the world gets denser (spacing on a unit sphere
// scales ~sqrt(4*pi / cellCount)). A threshold hand-tuned for one point count
// is not safe at another: 0.01 was fine at 3000 points (spacing ~0.065, an
// ~6.5x margin) but would leave almost no margin at 200k (spacing ~0.008).
// So derive the threshold from cell count instead of hardcoding it.
//   N=3000    -> spacing ~0.0647, threshold ~0.0097 (~14x the ~7e-4 drift, ~6.7x below spacing)
//   N=30000   -> spacing ~0.0205, threshold ~0.0031 (~4.4x the drift, ~6.7x below spacing)
//   N=200000  -> spacing ~0.0079, threshold ~0.0012 (~1.7x the drift, ~6.7x below spacing)
const WELD_FRACTION = 0.15; // threshold = WELD_FRACTION * spacing; keeps ~6.7x margin below true spacing at every N
const weldThresholdFor = (cellCount: number): number =>
  WELD_FRACTION * Math.sqrt((4 * Math.PI) / Math.max(cellCount, 1));

class VertexWelder {
  private grid = new Map<string, number[]>(); // cell key -> representative ids
  private reps: Point3[] = [];
  private readonly threshold: number;
  private readonly thresholdSq: number;
  private readonly cellSize: number;

  constructor(threshold: number) {
    this.threshold = threshold;
    this.thresholdSq = threshold * threshold;
    this.cellSize = threshold;
  }

  private cellIndex(v: number): number {
    return Math.floor(v / this.cellSize);
  }

  private cellKey(ix: number, iy: number, iz: number): string {
    return `${ix},${iy},${iz}`;
  }

  /** Returns a stable integer id for `p`, merging it with any existing
   * representative within this welder's threshold. */
  idFor(p: Point3): number {
    const ix = this.cellIndex(p[0]);
    const iy = this.cellIndex(p[1]);
    const iz = this.cellIndex(p[2]);
    for (let dx = -1; dx <= 1; dx++) {
      for (let dy = -1; dy <= 1; dy++) {
        for (let dz = -1; dz <= 1; dz++) {
          const cands = this.grid.get(this.cellKey(ix + dx, iy + dy, iz + dz));
          if (!cands) continue;
          for (const id of cands) {
            const r = this.reps[id];
            const distSq = (r[0] - p[0]) ** 2 + (r[1] - p[1]) ** 2 + (r[2] - p[2]) ** 2;
            if (distSq < this.thresholdSq) return id;
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
 *
 * `cellCount` scales the vertex-welding threshold (see weldThresholdFor above)
 * and is REQUIRED, not optional with a fallback — a silent default here would
 * be a silent behavior-changing default for any caller that forgot to pass
 * it, exactly the kind of bug this whole module exists to catch. Callers on
 * real generated geometry must pass `world.cells.length`; synthetic/test
 * callers should still pass an explicit count (any value is fine as long as
 * the resulting weld threshold stays well below the synthetic points'
 * spacing — see tests/boundaries.test.ts for examples).
 */
export const chainSegments = (
  segments: Array<[Point3, Point3]>,
  cellCount: number,
): Point3[][] => {
  const welder = new VertexWelder(weldThresholdFor(cellCount));
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
