// F2 screen-space overlay tenants — pure Canvas2D draw functions dispatched by
// ScreenOverlay. Each receives the projected cell screen coords, the world, and
// a projector for arbitrary local-frame points.

import { WorldData, Point, LabelVisibility } from '../../types';
import { ProjectedCells } from '../../utils/screenProject';
import { displayRadius } from '../../utils/displayRadius';
import { nearestCellWalk, nearestCellBrute } from '../../utils/nearestCell';
import { ContourSegment, contourStroke } from '../../utils/shading';
import { BorderSegment } from '../../utils/borders';
import { LabelAnchor } from '../../utils/labelAnchors';
import { MapLabel, drawMapLabels } from '../../utils/labels';
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

// --- contour lines (A4 isolines), migrated off 3D LineSegments ---

export const CONTOUR_LIFT = 0.002; // clears the terrain step it crowns
const CONTOUR_BASE_WIDTH = 1;

// Elevation isolines drawn in screen space. Unlike the other tenants this one
// takes its segments as an argument: computeContourSegments is an O(cells x
// neighbors) sweep, far too costly to redo on every redraw, so WorldViewer
// memoizes it on world identity and closes over the result.
//
// Radius (LocalProjector contract): each segment carries the height of the
// TALLER of its two cells. The globe mesh renders a cell boundary as a vertical
// step, so riding the taller cell's radius crowns that step. This replaces a
// fixed r = 1.053 that floated every isoline above all terrain — the worst
// parallax offender of the remaining overlays.
//
// Two passes so index contours draw over the intermediates they cross, matching
// drawContourPaths in utils/shading.ts (one style across globe, 2D, and export).
//
// NO ELEVATION LABELS. They shipped briefly and were pulled: anchors were picked
// from segment midpoints in array order and thinned greedily per redraw, so a
// different subset survived every frame and the readouts wandered across the
// terrain as the globe moved. Reviving them needs anchors chosen ONCE in world
// space (fixed segments per level), projected each frame and merely hidden on
// collision — never re-selected. `ContourSegment.elevation` and `contourLabel()`
// are kept as the seam for that, and for ROADMAP D8 (World Datum), which would
// make the readout metres instead of a percentage.
export function drawContoursTenant(
  ctx: CanvasRenderingContext2D,
  segments: ContourSegment[],
  project: LocalProjector,
  smooth: boolean,
): void {
  if (segments.length === 0) return;
  const p1: [number, number] = [0, 0];
  const p2: [number, number] = [0, 0];
  ctx.lineCap = 'round';

  for (const indexPass of [false, true]) {
    ctx.strokeStyle = contourStroke(indexPass);
    ctx.lineWidth = indexPass ? CONTOUR_BASE_WIDTH * 2 : CONTOUR_BASE_WIDTH;
    ctx.beginPath();
    for (const seg of segments) {
      if (seg.index !== indexPass) continue;
      const r = displayRadius(seg.height, smooth, CONTOUR_LIFT);
      // Both ends must be on the near hemisphere; a segment is one cell edge,
      // far too short to be worth clipping against the limb.
      if (!project(seg.a[0] * r, seg.a[1] * r, seg.a[2] * r, p1)) continue;
      if (!project(seg.b[0] * r, seg.b[1] * r, seg.b[2] * r, p2)) continue;
      ctx.moveTo(p1[0], p1[1]);
      ctx.lineTo(p2[0], p2[1]);
    }
    ctx.stroke();
  }
}


// --- faction borders, migrated off 3D LineSegments ---

export const BORDER_LIFT = 0.002; // clears the terrain step it crowns
const BORDER_COLOR = 'rgba(255,255,255,0.85)';
const BORDER_WIDTH = 1.4;

// Faction borders drawn in screen space.
//
// NOT a parallax fix. The 3D `FactionBorders` this replaces already used
// `max(displayRadius(hA), displayRadius(hB)) + 0.002` with RAW heights keyed off
// `smoothGlobe`, so its radius was correct — do not "fix" something here that
// was never broken. The migration buys three other things:
//
//  1. Horizon culling by the analytic limb test, instead of leaning on
//     `depthTest` against the mesh, which bleeds at grazing angles and across
//     the S19f plate overhang.
//  2. Real stroke width. `LineBasicMaterial.linewidth` is a no-op in WebGL on
//     every platform, so the old borders were locked at 1px whatever it said.
//  3. One paint order with the other overlays instead of a 3D object
//     interleaved among 2D tenants.
//
// Like contours, this takes precomputed segments: the extraction is
// O(cells x neighbors x vertices^2) (utils/borders.ts) and is memoized on world
// identity by WorldViewer, not recomputed per redraw.
//
// Radius (LocalProjector contract): each segment rides the TALLER of its two
// cells, so the line crowns the vertical step the mesh draws at a cell boundary.
// The height is RAW, never clamped to sea level — see drawGraticuleTenant for
// what that clamp cost in S18.
export function drawBordersTenant(
  ctx: CanvasRenderingContext2D,
  segments: BorderSegment[],
  project: LocalProjector,
  smooth: boolean,
): void {
  if (segments.length === 0) return;
  const p1: [number, number] = [0, 0];
  const p2: [number, number] = [0, 0];
  ctx.strokeStyle = BORDER_COLOR;
  ctx.lineWidth = BORDER_WIDTH;
  ctx.lineCap = 'round';
  ctx.beginPath();
  for (const seg of segments) {
    const r = displayRadius(seg.height, smooth, BORDER_LIFT);
    // Both ends must be on the near hemisphere; a segment is one cell edge,
    // far too short to be worth clipping against the limb.
    if (!project(seg.a.x * r, seg.a.y * r, seg.a.z * r, p1)) continue;
    if (!project(seg.b.x * r, seg.b.y * r, seg.b.z * r, p2)) continue;
    ctx.moveTo(p1[0], p1[1]);
    ctx.lineTo(p2[0], p2[1]);
  }
  ctx.stroke();
}

// --- rivers, migrated off 3D LineSegments ---

// Offset, like ROUTE_LIFT/CONTOUR_LIFT/BORDER_LIFT — NOT the old `RIVER_LIFT`
// that lived in `RiverLines`, which was 1.005, an ABSOLUTE radius. Smooth-globe
// radius here is `1 + RIVER_LIFT` (displayRadius(0, true, RIVER_LIFT)), i.e.
// 1.005 — the same number as before. Do not write `1 + 1.005`.
export const RIVER_LIFT = 0.005;
// Matches the old LineBasicMaterial in RiverLines, INCLUDING its alpha: that
// material was `opacity={0.8} transparent`, so a solid stroke here would make
// every river ~25% more prominent than before the migration.
const RIVER_COLOR = 'rgba(56,189,248,0.8)'; // #38bdf8 @ 0.8
const RIVER_WIDTH = 1.2;

// River polylines drawn in screen space. `polylines` is precomputed by
// `computeRiverPolylines` (utils/riverPaths.ts) and closed over by the tenant,
// same deal as contours/borders: CatmullRom smoothing over ~1741 paths is far
// too costly to redo on every redraw, so WorldViewer memoizes it on
// `world.rivers` alone (rivers are stable across paint strokes).
//
// THE TRAP this tenant exists to avoid: `world.rivers` points are NOT unit
// directions. `getRenderPoint` (utils/worldGen.ts ~L310) pre-scales each river
// point at generation time to `r = 1 + cell.height·0.05 + 0.005` — the same
// baked radius the terrain mesh itself uses for that cell, plus a small lift.
// Every other tenant receives unit directions and applies `displayRadius`
// itself; a river point has ALREADY been scaled. Applying `displayRadius` to it
// again would double-scale it and float every river off the globe. So on the
// RAISED globe this projects each point AS-IS — no radius math — and only on
// the SMOOTH globe does it normalize to unit and re-lift by RIVER_LIFT, exactly
// mirroring what the old `RiverLines` did when `smoothGlobe` was true.
//
// Horizon breaks exactly like drawRoutesTenant: `moveTo` on re-entry after a
// culled point, never `lineTo` across the gap.
export function drawRiversTenant(
  ctx: CanvasRenderingContext2D,
  polylines: Point[][],
  project: LocalProjector,
  smooth: boolean,
): void {
  if (polylines.length === 0) return;
  const pt: [number, number] = [0, 0];
  ctx.strokeStyle = RIVER_COLOR;
  ctx.lineWidth = RIVER_WIDTH;
  ctx.lineCap = 'round';
  ctx.beginPath();
  for (const path of polylines) {
    if (path.length < 2) continue;
    let drawing = false;
    for (const p of path) {
      let ok: boolean;
      if (smooth) {
        const len = Math.hypot(p.x, p.y, p.z) || 1;
        const r = (1 + RIVER_LIFT) / len;
        ok = project(p.x * r, p.y * r, p.z * r, pt);
      } else {
        // Already baked at the correct relief radius — project as-is.
        ok = project(p.x, p.y, p.z, pt);
      }
      if (ok) {
        if (drawing) ctx.lineTo(pt[0], pt[1]);
        else { ctx.moveTo(pt[0], pt[1]); drawing = true; }
      } else {
        drawing = false;
      }
    }
  }
  ctx.stroke();
}

// --- ruler arc (A5 great-circle measurement), migrated off 3D LineSegments ---

const RULER_COLOR = '#fbbf24';
// Deliberate FIXED radii, not draped: the ruler measures a great circle and
// must clear the highest peaks (relief tops out at 1.05), so it floats above
// terrain rather than riding it — the one case where a tenant does NOT track
// the mesh (see the LocalProjector radius contract). Mirrors the old 3D
// RulerArc's RULER_RADIUS / 1.01 constants exactly.
export const RULER_RADIUS_RAISED = 1.062;
export const RULER_RADIUS_SMOOTH = 1.01;
const RULER_DOT_PX = 3.5; // endpoint marker radius, screen pixels

// The measurement arc is a polyline over unit-sphere great-circle samples plus a
// dot at each endpoint. The polyline BREAKS at the horizon like every other
// tenant, so it never draws a chord across the far side. `points` are unit
// vectors (A5 computes them on the unit sphere); the tenant lifts them to the
// fixed ruler radius itself.
export function drawRulerTenant(
  ctx: CanvasRenderingContext2D,
  points: Point[],
  project: LocalProjector,
  smooth: boolean,
): void {
  if (points.length < 2) return;
  const r = smooth ? RULER_RADIUS_SMOOTH : RULER_RADIUS_RAISED;
  const pt: [number, number] = [0, 0];

  ctx.strokeStyle = RULER_COLOR;
  ctx.lineWidth = 1.4;
  ctx.beginPath();
  let drawing = false;
  for (const p of points) {
    if (project(p.x * r, p.y * r, p.z * r, pt)) {
      if (drawing) ctx.lineTo(pt[0], pt[1]);
      else { ctx.moveTo(pt[0], pt[1]); drawing = true; }
    } else {
      drawing = false;
    }
  }
  ctx.stroke();

  // Endpoint dots, drawn only when the endpoint is on the near hemisphere.
  ctx.fillStyle = RULER_COLOR;
  for (const end of [points[0], points[points.length - 1]]) {
    if (project(end.x * r, end.y * r, end.z * r, pt)) {
      ctx.beginPath();
      ctx.arc(pt[0], pt[1], RULER_DOT_PX, 0, Math.PI * 2);
      ctx.fill();
    }
  }
}

// --- point labels (capitals/provinces/towns/geography/markers), migrated off
// canvas-texture sprites (F2 S22, last tenant) ---

// Above ROUTE_LIFT (0.008) so labels clear the overlays they sit on.
export const LABEL_LIFT = 0.01;

// `drawMapLabels`'s own LOD compares `scale` against `ZOOM_THRESHOLDS`, whose
// max entry is 2.0. This tenant applies its OWN LOD first (the camDist cutoffs
// below, kept verbatim from the old 3D `PointLabels` — see plan §3, "do not
// unify the zoom LOD"), so `scale` must be pushed past every threshold to be a
// no-op. Infinity is the most self-evident way to say "LOD off" at the call site.
const LABEL_SCALE_LOD_OFF = Infinity;

// Chosen to land close to the old 3D sprites' on-screen size at the plan's
// measurement baseline (fov 45, 1000px-tall viewport, camDist 2.5 — the
// "default" row of the §2 parallax table). At that baseline the old capital
// sprite (world-space height 0.042, canvas padding baked in) rendered its
// glyph at roughly 13px on screen; `LABEL_CONFIG.capital.baseFontSize` is 11,
// so 11 * 1.2 ≈ 13.2 lands on the same target. `drawMapLabels` decouples
// `fontScale` from `scale` (verified in plan §3), so this is a fixed on-screen
// size — unlike the old sprites, whose size shrank with camera distance the
// way any 3D billboard does. Nobody asked for that zoom-linked shrink to
// survive the migration; the camDist LOD below keeps *which* labels show, not
// *how big* they render.
const LABEL_FONT_SCALE = 1.2;

// Geographic LOD tiers, transcribed verbatim from the old 3D `PointLabels`
// (GEO_SMALL_KINDS / GEO_MID_KINDS), which itself already mirrors the 2D
// ZOOM_THRESHOLDS relationship inversely (camDist UP = zoomed OUT, so a
// smaller camDist cutoff means the kind disappears sooner while zooming out).
const LOD_SMALL_GEO_KINDS = new Set(['island', 'lake']);
const LOD_MID_GEO_KINDS = new Set(['sea', 'range', 'desert', 'forest']);

// True when `label` should be considered for this frame at `camDist`, before
// declutter/visibility (drawMapLabels handles per-kind visibility and overlap
// itself). Faction labels are excluded unconditionally: they stay 3D curved
// meshes (CurvedFactionLabel/CountryLabels, untouched by F2 — plan §1), so this
// tenant must never emit one, on pain of a duplicate flat label under the
// curved one.
function passesGlobeLod(label: MapLabel, camDist: number): boolean {
  if (label.kind === 'faction') return false;
  if (label.kind === 'town') return camDist <= 2;
  if (label.kind === 'province') return camDist <= 3;
  if (label.kind === 'capital') return camDist <= 4;
  if (LOD_SMALL_GEO_KINDS.has(label.kind)) return camDist <= 2.5;
  if (LOD_MID_GEO_KINDS.has(label.kind)) return camDist <= 3.5;
  return true;
}

// Point labels drawn in screen space, reusing `drawMapLabels` (utils/labels.ts)
// for style + declutter instead of duplicating them a third time — see plan
// §3. `anchors` is precomputed by `computeLabelAnchors` (utils/labelAnchors.ts)
// and closed over by WorldViewer like contours/borders/rivers: a per-redraw
// nearest-cell scan is far too costly at 200k cells (plan §2).
//
// Radius (LocalProjector contract): each anchor carries the height of the cell
// nearest the label's position, RAW, never clamped to sea level — an ocean
// label rides the seafloor radius like every other tenant (see
// drawGraticuleTenant for what that clamp cost in S18). `displayRadius` turns
// that height into the drape radius per `smooth`, exactly mirroring
// ContourSegment.height / BorderSegment.height.
//
// `camDist` cannot be read inside a tenant (no camera access), so WorldViewer
// captures it (`camera.position.length()`, same as the old `PointLabels`) and
// passes it in — a closed-over argument, the established pattern here (see
// drawContoursTenant/drawBordersTenant taking precomputed segments).
//
// Declutter is a greedy first-wins pass over `anchors` in the order given —
// this function does not sort. `collectLabels` (utils/labels.ts) already
// emits labels ascending by `priority`, and `computeLabelAnchors` preserves
// input order verbatim, so the array WorldViewer closes over is already
// priority-sorted by construction. Passing an unsorted array here would let
// the wrong label win a collision silently; there is no runtime guard for it.
export function drawLabelsTenant(
  ctx: CanvasRenderingContext2D,
  anchors: LabelAnchor[],
  project: LocalProjector,
  smooth: boolean,
  camDist: number,
  visibility: LabelVisibility,
): void {
  if (anchors.length === 0) return;

  const visibleAnchors = anchors.filter((a) => passesGlobeLod(a.label, camDist));
  if (visibleAnchors.length === 0) return;

  // `drawMapLabels` only sees `label.position`; it has no notion of `height`.
  // Keyed by the position OBJECT (each label owns exactly one, reused verbatim
  // through `anchors`), so the radius lookup below survives the trip through
  // `drawMapLabels`'s generic `project(position)` callback.
  const heightByPosition = new Map<MapLabel['position'], number>();
  for (const a of visibleAnchors) heightByPosition.set(a.label.position, a.height);

  const pt: [number, number] = [0, 0];
  const projectLabel = (position: MapLabel['position']): [number, number] | null => {
    const height = heightByPosition.get(position) ?? 0;
    const r = displayRadius(height, smooth, LABEL_LIFT);
    if (!project(position.x * r, position.y * r, position.z * r, pt)) return null;
    return [pt[0], pt[1]];
  };

  drawMapLabels(
    ctx,
    visibleAnchors.map((a) => a.label),
    projectLabel,
    LABEL_SCALE_LOD_OFF,
    visibility,
    LABEL_FONT_SCALE,
  );
}
