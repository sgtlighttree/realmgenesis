// F2 screen-space overlay tenants — pure Canvas2D draw functions dispatched by
// ScreenOverlay. Each receives the projected cell screen coords, the world, and
// a projector for arbitrary local-frame points.

import { WorldData } from '../../types';
import { ProjectedCells } from '../../utils/screenProject';
import { displayRadius } from '../../utils/displayRadius';
import { nearestCellWalk, nearestCellBrute } from '../../utils/nearestCell';
import { LocalProjector } from './ScreenOverlay';

// --- currents draw constants (single source; also consumed by Map2D, Task 6) ---
export const CURRENT_SPEED_MIN = 0.04; // skip near-still cells
export const CURRENT_ARC = 0.05;       // max arrow arc length (unit-sphere)
export const CURRENT_SPEED_REF = 0.6;  // speed at which arrows reach full length
export const CURRENT_STRIDE = 1;       // draw every Nth qualifying ocean cell
export const SST_CLAMP = 6;            // °C anomaly mapped to full cold/warm
const COLD: [number, number, number] = [0x3b, 0x6f, 0xb0];
const WARM: [number, number, number] = [0xc0, 0x44, 0x2e];

// Warm/cold tint from SST anomaly (°C), clamped to ±SST_CLAMP.
export function currentTint(sst: number): string {
  const t = Math.max(0, Math.min(1, (sst + SST_CLAMP) / (2 * SST_CLAMP)));
  const r = Math.round(COLD[0] + (WARM[0] - COLD[0]) * t);
  const g = Math.round(COLD[1] + (WARM[1] - COLD[1]) * t);
  const b = Math.round(COLD[2] + (WARM[2] - COLD[2]) * t);
  return `rgb(${r},${g},${b})`;
}

export function drawCurrentsTenant(
  ctx: CanvasRenderingContext2D,
  proj: ProjectedCells,
  world: WorldData,
  project: LocalProjector,
  smooth = false,
): void {
  const cur = world.currents;
  if (!cur) return;
  const cells = world.cells;
  const sea = world.params.seaLevel;
  const tip: [number, number] = [0, 0];
  ctx.lineWidth = 1.2;
  for (let i = 0; i < proj.n; i++) {
    if (!proj.visible[i]) continue;
    if (cells[i].height >= sea) continue; // ocean cells only
    if (CURRENT_STRIDE > 1 && i % CURRENT_STRIDE !== 0) continue;
    const sp = Math.hypot(cur.vx[i], cur.vy[i], cur.vz[i]);
    if (sp < CURRENT_SPEED_MIN) continue;
    // arrow tip = center + unit(vel) * (arc scaled by speed, capped)
    const mag = CURRENT_ARC * Math.min(1, sp / CURRENT_SPEED_REF);
    const k = mag / sp;
    const c = cells[i].center;
    // Lift the tip to the cell's rendered radius so it sits on the terrain
    // surface exactly like the arrow base (see LocalProjector radius contract).
    const r = displayRadius(cells[i].height, smooth);
    if (!project((c.x + cur.vx[i] * k) * r, (c.y + cur.vy[i] * k) * r, (c.z + cur.vz[i] * k) * r, tip)) continue;
    const color = currentTint(cur.sst[i]);
    ctx.strokeStyle = color;
    ctx.beginPath();
    ctx.moveTo(proj.x[i], proj.y[i]);
    ctx.lineTo(tip[0], tip[1]);
    ctx.stroke();
    ctx.fillStyle = color;
    ctx.beginPath();
    ctx.arc(tip[0], tip[1], 1.4, 0, Math.PI * 2);
    ctx.fill();
  }
}

const D2R = Math.PI / 180;
const GRAT_SEG = 96; // samples per line (smooth circles)

// 10° graticule drawn in screen space: parallels −80..80°, meridians 0..350°,
// each sampled on the sphere, projected via the horizon-culling projector.
// The polyline breaks wherever a sample crosses behind the limb, so lines never
// draw across the globe silhouette (the win over the old always-visible 3D grid).
//
// Radius: the grid is not cell-bound, so it needs a deliberate radius (see the
// LocalProjector contract). On the SMOOTH globe every cell is r=1, so the grid
// sits on the unit sphere (zero parallax). On the RAISED globe it DRAPES: each
// sample is projected at the terrain radius at its lat/lon (nearest cell height,
// RAW — see below) so the line rides relief and meets terrain at coastlines,
// with zero parallax at any zoom (the Google-Earth read). The nearest cell is
// found by a greedy Voronoi hill-climb (utils/nearestCell) seeded by the previous
// sample — 1–3 hops along a polyline — so the whole grid costs one O(n) brute
// seed + ~5100 short walks per (gated) redraw.
//
// The height is RAW, never clamped to sea level. S18 shipped a
// `max(height, seaLevel)` clamp, reasoning the grid should ride the water
// surface; there is no water surface — the mesh renders ocean cells at their
// true seafloor radius (`displayRadius(cell.height, smooth)`, WorldViewer
// refill) and merely colours them blue. The clamp therefore floated the grid
// above every ocean cell by `(seaLevel − height)·0.05` — most of the globe, and
// the residual parallax Matt reported after S18. Clamping again reintroduces it.
export function drawGraticuleTenant(
  ctx: CanvasRenderingContext2D,
  _proj: ProjectedCells,
  world: WorldData,
  project: LocalProjector,
  smooth = false,
): void {
  ctx.strokeStyle = 'rgba(255,255,255,0.28)';
  ctx.lineWidth = 1;
  const pt: [number, number] = [0, 0];
  const cells = world.cells;
  // Running nearest-cell id carried across the whole draw; -1 = seed via brute
  // scan on the first sample, then hill-climb from the previous sample.
  let startId = -1;

  // Projects a UNIT direction, applying the drape radius on the raised globe.
  const projSample = (x: number, y: number, z: number): boolean => {
    if (smooth || cells.length === 0) return project(x, y, z, pt); // unit sphere
    startId = startId < 0 ? nearestCellBrute(cells, x, y, z) : nearestCellWalk(cells, x, y, z, startId);
    const r = displayRadius(cells[startId].height, false);
    return project(x * r, y * r, z * r, pt);
  };

  // parallels (constant latitude)
  for (let lat = -80; lat <= 80; lat += 10) {
    const la = lat * D2R;
    const cy = Math.sin(la);
    const cr = Math.cos(la);
    let drawing = false;
    ctx.beginPath();
    for (let s = 0; s <= GRAT_SEG; s++) {
      const lon = (s / GRAT_SEG) * Math.PI * 2;
      if (projSample(cr * Math.cos(lon), cy, cr * Math.sin(lon))) {
        if (drawing) ctx.lineTo(pt[0], pt[1]);
        else { ctx.moveTo(pt[0], pt[1]); drawing = true; }
      } else {
        drawing = false;
      }
    }
    ctx.stroke();
  }

  // meridians (constant longitude, pole to pole)
  for (let lon = 0; lon < 360; lon += 10) {
    const lo = lon * D2R;
    const cl = Math.cos(lo);
    const sl = Math.sin(lo);
    let drawing = false;
    ctx.beginPath();
    for (let s = 0; s <= GRAT_SEG; s++) {
      const lat = (s / GRAT_SEG) * Math.PI - Math.PI / 2;
      const cla = Math.cos(lat);
      if (projSample(cla * cl, Math.sin(lat), cla * sl)) {
        if (drawing) ctx.lineTo(pt[0], pt[1]);
        else { ctx.moveTo(pt[0], pt[1]); drawing = true; }
      } else {
        drawing = false;
      }
    }
    ctx.stroke();
  }
}

// --- routes (C3 roads + sea trade routes), migrated off 3D LineSegments ---

export const ROUTE_LIFT = 0.008; // sits just above the surface, over rivers
const ROAD_COLOR = '#c8a25a';
const SEAROUTE_COLOR = '#5eb8c8';
const SEAROUTE_DASH = [4, 3];
const KINDS = ['road', 'searoute'] as const;

// Roads and sea routes drawn in screen space. Each route is a polyline over its
// path of cell centers; the polyline BREAKS wherever a point falls past the
// horizon, so a route never draws a chord across the globe silhouette — the win
// over the old 3D LineSegments, which also drew at a fixed r = 1.008 and so sank
// into any terrain above that (mountains reach 1.05).
//
// Radius (LocalProjector contract): on the SMOOTH globe every cell is r=1, so
// routes sit on the unit sphere plus ROUTE_LIFT. On the RAISED globe each point
// DRAPES at its own cell's terrain radius, so a road climbs a mountain range
// instead of tunnelling through it. Routes are already cell-bound — paths are
// built from cell centers in utils/routes.ts and RouteData carries the parallel
// cellIds — so this needs no nearestCell walk, unlike the graticule.
//
// The height is RAW, never clamped to sea level: sea routes run over ocean, and
// the mesh renders ocean cells at their true seafloor radius. Clamping would
// float them, which is the S18 graticule bug (see drawGraticuleTenant).
//
// Fallback: a route without usable cellIds draws flat at the sea-level radius.
//
// No curve smoothing: cell centers are dense enough at any cell count, and
// resampling a spline would have to re-derive the horizon breaks.
export function drawRoutesTenant(
  ctx: CanvasRenderingContext2D,
  _proj: ProjectedCells,
  world: WorldData,
  project: LocalProjector,
  smooth = false,
): void {
  const routes = world.routes;
  if (!routes || routes.length === 0) return;
  const cells = world.cells;
  const pt: [number, number] = [0, 0];
  // Flat radius: used on the smooth globe, and as the fallback for a route
  // whose cellIds are missing or out of step with its path.
  const flat = displayRadius(smooth ? 0 : world.params.seaLevel, smooth, ROUTE_LIFT);

  ctx.lineWidth = 1.5;
  for (const kind of KINDS) {
    ctx.strokeStyle = kind === 'road' ? ROAD_COLOR : SEAROUTE_COLOR;
    ctx.setLineDash(kind === 'road' ? [] : SEAROUTE_DASH);
    for (const route of routes) {
      if (route.kind !== kind || route.path.length < 2) continue;
      const ids = route.cellIds;
      const drape = !smooth && !!ids && ids.length === route.path.length;
      let drawing = false;
      ctx.beginPath();
      for (let k = 0; k < route.path.length; k++) {
        const p = route.path[k];
        const cell = drape ? cells[ids[k]] : undefined;
        const rad = cell ? displayRadius(cell.height, false, ROUTE_LIFT) : flat;
        if (project(p.x * rad, p.y * rad, p.z * rad, pt)) {
          if (drawing) ctx.lineTo(pt[0], pt[1]);
          else { ctx.moveTo(pt[0], pt[1]); drawing = true; }
        } else {
          drawing = false;
        }
      }
      ctx.stroke();
    }
  }
  // The 2D context is shared across tenants — hand it back undashed.
  ctx.setLineDash([]);
}
