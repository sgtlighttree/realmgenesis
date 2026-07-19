import { Cell } from '../types';
import { Point3, toLonLat } from './geo';

// d3.GeoPath-compatible callable; typed structurally so this module stays
// free of a d3 import (only the canvas-bound path generator is ever passed,
// and only ever with LineString geometry).
type GeoPathLike = (object: { type: 'LineString'; coordinates: [number, number][] }) => unknown;

// Hillshade multiplier clamp band. Centered on 1.0 so flat land and water both
// read neutral; the small overshoot above 1.0 is a mild highlight allowance for
// slopes facing the light. Exported so the test band-check can't drift.
export const SHADE_MIN = 0.6;
export const SHADE_MAX = 1.15;

// Fixed world-space light direction. Classic NW-ish cartographic convention
// (upper-left) so relief reads the way a reader expects on a map.
const LIGHT: Point3 = (() => {
  const v: Point3 = [-0.6, 0.5, -0.6];
  const len = Math.hypot(v[0], v[1], v[2]);
  return [v[0] / len, v[1] / len, v[2] / len];
})();

// Slope exaggeration and shading strength. Heights are 0-1 and neighbor
// differences are small, so the tangential gradient is scaled up hard to make
// relief legible; both are visual constants, not physical values.
const EXAGGERATION = 1.6;
const STRENGTH = 1.6;

const dot3 = (a: Point3, b: Point3): number => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];

// Per-cell Lambert-style hillshade factor in [SHADE_MIN, SHADE_MAX].
//
// The shade is the RELIEF component only, decoupled from where the cell sits on
// the globe: we subtract the radial baseline dot(unitNormal, light) from the
// perturbed dot(surfaceNormal, light). Flat terrain has no perturbation, so it
// lands at exactly 1.0 — same neutral value as water — instead of swinging
// across a globe-wide day/night terminator. This keeps the pass usable as an
// overlay in any view (the 3D globe already has its own directional lighting).
export const computeShadeMap = (cells: Cell[], seaLevel: number): Float32Array => {
  const shade = new Float32Array(cells.length);

  for (let i = 0; i < cells.length; i++) {
    const cell = cells[i];
    // Water carries no relief shading.
    if (cell.height < seaLevel) {
      shade[cell.id] = 1.0;
      continue;
    }

    const c = cell.center;
    const nLen = Math.hypot(c.x, c.y, c.z) || 1;
    const u: Point3 = [c.x / nLen, c.y / nLen, c.z / nLen]; // unit (radial) normal

    // Average tangential uphill gradient across neighbors.
    let gx = 0, gy = 0, gz = 0;
    let count = 0;
    for (const nId of cell.neighbors) {
      const nb = cells[nId];
      if (!nb) continue;
      let dx = nb.center.x - c.x;
      let dy = nb.center.y - c.y;
      let dz = nb.center.z - c.z;
      // Project the direction onto the tangent plane at the cell.
      const radial = dx * u[0] + dy * u[1] + dz * u[2];
      dx -= radial * u[0];
      dy -= radial * u[1];
      dz -= radial * u[2];
      const run = Math.hypot(dx, dy, dz);
      if (run < 1e-9) continue;
      const slope = (nb.height - cell.height) / run; // rise over angular run
      gx += (dx / run) * slope;
      gy += (dy / run) * slope;
      gz += (dz / run) * slope;
      count++;
    }
    if (count > 0) {
      gx /= count; gy /= count; gz /= count;
    }

    // Tilt the normal downhill (negative gradient) and renormalize.
    let nx = u[0] - EXAGGERATION * gx;
    let ny = u[1] - EXAGGERATION * gy;
    let nz = u[2] - EXAGGERATION * gz;
    const len = Math.hypot(nx, ny, nz) || 1;
    nx /= len; ny /= len; nz /= len;

    const relief = dot3([nx, ny, nz], LIGHT) - dot3(u, LIGHT);
    let s = 1 + STRENGTH * relief;
    if (s < SHADE_MIN) s = SHADE_MIN;
    else if (s > SHADE_MAX) s = SHADE_MAX;
    shade[cell.id] = s;
  }

  return shade;
};

const SHARED_EPS_SQ = 0.000001;

// Shared Voronoi edge between two neighboring cells, found the same way
// FactionBorders does: two vertices within epsilon of each other.
const sharedEdge = (a: Cell, b: Cell): [Point3, Point3] | null => {
  const shared: Point3[] = [];
  for (const vA of a.vertices) {
    for (const vB of b.vertices) {
      const distSq = (vA.x - vB.x) ** 2 + (vA.y - vB.y) ** 2 + (vA.z - vB.z) ** 2;
      if (distSq < SHARED_EPS_SQ) {
        shared.push([vA.x, vA.y, vA.z]);
        break;
      }
    }
    if (shared.length === 2) break;
  }
  return shared.length === 2 ? [shared[0], shared[1]] : null;
};

// Elevation isolines as chunky cell-edge segments. For each pair of neighboring
// LAND cells whose heights fall in different contour bands, emit their shared
// Voronoi edge. Bands are seaLevel + k*interval (k>=1), so any band difference
// between two land cells implies a crossed contour level above seaLevel.
// Draws contour segments through a projection-bound d3 path generator with the
// shared subdued stroke used by every 2D pipeline (Map2D + PNG export).
export const drawContourPaths = (
  ctx: CanvasRenderingContext2D,
  pathGenerator: GeoPathLike,
  segments: Array<[Point3, Point3]>,
  lineWidth: number,
): void => {
  if (segments.length === 0) return;
  ctx.save();
  ctx.strokeStyle = 'rgba(58, 42, 26, 0.4)';
  ctx.lineWidth = lineWidth;
  ctx.lineJoin = 'round';
  ctx.lineCap = 'round';
  segments.forEach(([a, b]) => {
    ctx.beginPath();
    pathGenerator({
      type: 'LineString',
      coordinates: [toLonLat(a), toLonLat(b)],
    });
    ctx.stroke();
  });
  ctx.restore();
};

export const computeContourSegments = (
  cells: Cell[],
  seaLevel: number,
  interval: number,
): Array<[Point3, Point3]> => {
  const segments: Array<[Point3, Point3]> = [];
  if (interval <= 0) return segments;

  for (const a of cells) {
    if (a.height < seaLevel) continue;
    const bandA = Math.floor((a.height - seaLevel) / interval);
    for (const nId of a.neighbors) {
      const b = cells[nId];
      if (!b || a.id >= b.id) continue; // process each pair once
      if (b.height < seaLevel) continue; // both cells must be land
      const bandB = Math.floor((b.height - seaLevel) / interval);
      if (bandA === bandB) continue; // no contour level between them
      const edge = sharedEdge(a, b);
      if (edge) segments.push(edge);
    }
  }

  return segments;
};
