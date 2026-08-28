import { Cell } from '../types';
import { Point3, toLonLat } from './geo';
import { formatElevation } from './datum';

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
    // D10: water DOES carry relief shading. It used to short-circuit to 1.0, so
    // the sea bed rendered as a flat colour ramp in Map2D and every export no
    // matter what the bathymetry said — half of the "sea bed looks flat" report.
    //
    // Two rules keep it honest:
    //   * A water cell averages its gradient over WATER neighbours only. Mixing
    //     in a coastal land neighbour puts the full land/ocean height step into
    //     the gradient, which would draw a hard shaded rim around every
    //     coastline instead of sea-floor relief. Land cells are UNCHANGED — they
    //     still average over all neighbours, so existing maps do not shift.
    //   * Same STRENGTH and clamp band as land, no separate water constant.
    //     Ocean and land now carry equal measured texture at the default
    //     `seafloorRelief` (0.99x, see D10), so one strength gives both the same
    //     visual weight. Flat sea floor still lands at exactly 1.0 on its own,
    //     because the radial baseline is subtracted below.
    //
    // Lakes are unaffected: their height is >= seaLevel, so they take the land
    // path exactly as before.
    const isWater = cell.height < seaLevel;

    const c = cell.center;
    const nLen = Math.hypot(c.x, c.y, c.z) || 1;
    const u: Point3 = [c.x / nLen, c.y / nLen, c.z / nLen]; // unit (radial) normal

    // Average tangential uphill gradient across neighbors.
    let gx = 0, gy = 0, gz = 0;
    let count = 0;
    for (const nId of cell.neighbors) {
      const nb = cells[nId];
      if (!nb) continue;
      if (isWater && nb.height >= seaLevel) continue; // water shades off water only
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

// One elevation isoline: the shared Voronoi edge between two land cells that
// fall in different contour bands.
export interface ContourSegment {
  a: Point3;
  b: Point3;
  // Height of the HIGHER of the two cells. The globe mesh renders a cell
  // boundary as a vertical step, so drawing the isoline at the taller cell's
  // radius crowns that step instead of cutting into it or floating over it.
  // Unit-sphere consumers (Map2D, PNG export) ignore this.
  height: number;
  // Index contour: every CONTOUR_INDEX_EVERY-th level, drawn thick and bright
  // so elevation can be counted by eye (standard topographic convention).
  index: boolean;
  // Normalized height of the contour LINE itself (seaLevel + level*interval) —
  // not either cell's height. This is what an elevation label must read out.
  elevation: number;
}

// Every Nth contour level is an index contour — the standard cartographic 5th,
// which is only meaningful because the interval now adapts to give ~20 levels.
// (It was briefly 2, to salvage a hardcoded 0.1 interval that yielded four
// possible levels on a default world. The interval was the real bug.)
export const CONTOUR_INDEX_EVERY = 5;

// Contour levels to aim for between sea level and the highest land.
export const CONTOUR_TARGET_LEVELS = 20;

// Preferred intervals, in normalized height. Snapping to one of these keeps
// label values readable (2%, 5%) and spacing roughly comparable between worlds,
// while still adapting: a fixed interval cannot serve both a flat world and an
// alpine one. Default params put sea level at 0.55, so land spans only ~0.45 —
// the old fixed 0.1 gave four levels, and most terrain saw one or two lines.
const NICE_INTERVALS = [0.002, 0.005, 0.01, 0.02, 0.025, 0.05, 0.1];

/**
 * Interval for `computeContourSegments`, derived from the world's actual
 * relief. Returns 0 when there is no land above sea level, which
 * computeContourSegments treats as "no contours".
 *
 * O(cells) — memoize per world alongside the segments themselves.
 */
export const contourInterval = (
  cells: Cell[],
  seaLevel: number,
  targetLevels: number = CONTOUR_TARGET_LEVELS,
): number => {
  let max = seaLevel;
  for (const c of cells) if (c.height > max) max = c.height;
  const range = max - seaLevel;
  if (range <= 0 || targetLevels <= 0) return 0;

  const ideal = range / targetLevels;
  let best = NICE_INTERVALS[0];
  let bestErr = Infinity;
  for (const step of NICE_INTERVALS) {
    const err = Math.abs(Math.log(step / ideal)); // ratio distance, not absolute
    if (err < bestErr) { bestErr = err; best = step; }
  }
  return best;
};

const CONTOUR_INK = '255, 244, 224'; // warm off-white, legible on land and ice
const CONTOUR_ALPHA = 0.38;
const CONTOUR_INDEX_ALPHA = 0.75;
const CONTOUR_INDEX_WIDTH = 2; // multiplier on the caller's base line width

export const contourStroke = (index: boolean): string =>
  `rgba(${CONTOUR_INK}, ${index ? CONTOUR_INDEX_ALPHA : CONTOUR_ALPHA})`;

/**
 * Elevation readout for an index contour, in metres against the D8a
 * presentation datum (`utils/datum.ts`) — the same units the Inspector shows.
 *
 * NOTE: contour elevation labels are DORMANT — deliberately not rendered
 * (HANDOFF S19e). This is the documented single change point for when they are
 * revived, kept in metres so it never re-earns the old %-vs-km inconsistency.
 */
export const contourLabel = (
  elevation: number,
  seaLevel: number,
  maxElevationM?: number,
): string => formatElevation(elevation, seaLevel, maxElevationM);

// Draws contour segments through a projection-bound d3 path generator. Two
// passes so index contours land on top of the intermediates they cross.
export const drawContourPaths = (
  ctx: CanvasRenderingContext2D,
  pathGenerator: GeoPathLike,
  segments: ContourSegment[],
  lineWidth: number,
): void => {
  if (segments.length === 0) return;
  ctx.save();
  ctx.lineJoin = 'round';
  ctx.lineCap = 'round';
  for (const indexPass of [false, true]) {
    ctx.strokeStyle = contourStroke(indexPass);
    ctx.lineWidth = indexPass ? lineWidth * CONTOUR_INDEX_WIDTH : lineWidth;
    for (const seg of segments) {
      if (seg.index !== indexPass) continue;
      ctx.beginPath();
      pathGenerator({
        type: 'LineString',
        coordinates: [toLonLat(seg.a), toLonLat(seg.b)],
      });
      ctx.stroke();
    }
  }
  ctx.restore();
};

// For each pair of neighboring LAND cells whose heights fall in different
// contour bands, emit their shared Voronoi edge. Bands are seaLevel + k*interval
// (k>=1), so any band difference between two land cells implies a crossed level.
export const computeContourSegments = (
  cells: Cell[],
  seaLevel: number,
  interval: number,
): ContourSegment[] => {
  const segments: ContourSegment[] = [];
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
      if (!edge) continue;
      // Lowest level crossed between the two bands identifies the isoline.
      const level = Math.min(bandA, bandB) + 1;
      segments.push({
        a: edge[0],
        b: edge[1],
        height: Math.max(a.height, b.height),
        index: level % CONTOUR_INDEX_EVERY === 0,
        elevation: seaLevel + level * interval,
      });
    }
  }

  return segments;
};
