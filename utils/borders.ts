import { Cell, Point } from '../types';

// Faction border edges extracted from the cell graph (F2 borders tenant).
//
// Shape mirrors `ContourSegment` (utils/shading.ts) deliberately: both are
// precomputed edge lists that the ScreenOverlay draws at a per-segment radius,
// so keeping them structurally identical keeps the two tenants readable side by
// side.
export interface BorderSegment {
  a: Point;
  b: Point;
  // Height of the HIGHER of the two cells sharing this edge. The globe mesh
  // renders a cell boundary as a vertical step, so riding the taller cell's
  // radius crowns that step. `displayRadius` is affine and monotonic in height,
  // so max-then-map equals map-then-max — storing the scalar is equivalent to
  // storing max(displayRadius(hA), displayRadius(hB)).
  height: number;
}

// Two shared vertices are "the same" corner when they land within this squared
// distance on the unit sphere. Same tolerance the 3D FactionBorders used.
const SHARED_EPS_SQ = 0.000001;

/**
 * Every edge where the region on one side differs from the region on the other,
 * including faction-vs-unclaimed. Cells with no region on BOTH sides never
 * match, because `undefined !== undefined` is false.
 *
 * O(cells x neighbors x vertices^2) — the nested vertex scan makes this far too
 * costly for a per-frame redraw. Memoize on world identity and close over the
 * result, exactly as WorldViewer does for contour segments.
 */
export const computeBorderSegments = (cells: Cell[]): BorderSegment[] => {
  const segments: BorderSegment[] = [];
  for (const cellA of cells) {
    for (const nId of cellA.neighbors) {
      const cellB = cells[nId];
      if (!cellB || cellA.id >= cellB.id) continue; // each pair once
      if (cellA.regionId === cellB.regionId) continue;
      if (!cellA.vertices || !cellB.vertices) continue;

      // The shared edge is the pair of corners both cells carry.
      let v1: Point | null = null;
      let v2: Point | null = null;
      for (const vA of cellA.vertices) {
        for (const vB of cellB.vertices) {
          const dx = vA.x - vB.x;
          const dy = vA.y - vB.y;
          const dz = vA.z - vB.z;
          if (dx * dx + dy * dy + dz * dz < SHARED_EPS_SQ) {
            if (v1) v2 = vA;
            else v1 = vA;
            break;
          }
        }
        if (v2) break;
      }

      if (v1 && v2) {
        segments.push({ a: v1, b: v2, height: Math.max(cellA.height, cellB.height) });
      }
    }
  }
  return segments;
};
