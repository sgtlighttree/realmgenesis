import { BiomeType, Cell, Point, WorldData } from '../types';
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
