import * as d3 from 'd3';

import { MarkerData, Point, RouteData, WorldData, LabelVisibility } from '../../types';
import { toLonLat } from '../../utils/geo';
import { drawMapLabels, MapLabel } from '../../utils/labels';
import { OverlayInk } from '../../utils/mapStyle/overlayInk';
import { LabelTheme } from '../../utils/mapStyle/labelTheme';
import { ContourSegment, drawContourPaths } from '../../utils/shading';
import { drawCurrents2D } from '../overlays/currents2D';

// Labels are precomputed by the host via `collectLabels` (`utils/labels.ts`)
// and passed in through `mapLabels` — this file never recomputes them itself.

/**
 * Per-frame overlay draw sequence for the Canvas2D vector layer that sits
 * above the WebGL cell-fill surface (F3 Phase 2). This module is the ONE
 * place that sequence lives — `VectorOverlay.tsx` calls `paintVectorOverlay`
 * after applying the composed transform; nothing here touches `ctx`'s
 * transform itself.
 *
 * Task 6 must route `Map2D.tsx`'s offscreen (non-Dymaxion) render path
 * through this same helper instead of keeping its own inline copy of this
 * sequence (currently ~Map2D.tsx:756-894) — see the module-level functions
 * exported below (`drawRiverPaths`, `drawMarkerPins`, `strokeBoundaryPaths`,
 * `drawRoutePaths`), which Task 6 should import here rather than leaving
 * Map2D's private duplicates in place.
 */

/** Toggles for the optional overlay passes — mirrors Map2D's show* props. */
export interface VectorOverlayToggles {
  showGrid: boolean;
  showRivers: boolean;
  showRoutes: boolean;
  showCurrents: boolean;
}

/**
 * Cached, already-projected (CSS-px, un-flipped) boundary geometry from
 * `MapGeometryCache` (`utils/mapCache.ts`). These are resolution-independent
 * Path2Ds — stroking one at any `strokeScale` yields a crisp line at that
 * zoom. The host (Task 6) swaps this prop's Path2Ds on zoom-settle for
 * versions rebuilt at a lower simplification tolerance (`tol = 0.5 / scale`
 * instead of the base-fit `tol = 0.5`).
 *
 * IMPORTANT: this is settle-refined GEOMETRY (vector-tile LOD), not
 * settle-refined PIXELS — the overlay itself is resolution-independent at
 * EVERY zoom (it is redrawn through the live composed transform every
 * frame); only the simplification budget behind these Path2Ds tracks zoom,
 * to keep vertex counts sane at deep zoom without ever going blurry.
 */
export interface VectorOverlayBoundaryPaths {
  coast: Path2D;
  borders: Path2D;
  lakes: Path2D;
}

/** Everything `paintVectorOverlay` needs to run the draw sequence once. */
export interface PaintVectorOverlayParams {
  world: WorldData;
  /** Same `d3.GeoProjection` the boundary/contour/river/route passes project through. */
  projection: d3.GeoProjection;
  size: { width: number; height: number };
  toggles: VectorOverlayToggles;
  boundaryPaths: VectorOverlayBoundaryPaths;
  /** Precomputed by the host via `computeContourSegments` — this file never derives contours itself. */
  contourSegments: ContourSegment[];
  /** Precomputed by the host via `collectLabels` — see the note on `_collectLabels` above. */
  mapLabels: MapLabel[];
  labelVisibility: LabelVisibility;
  labelTheme: LabelTheme;
  overlayInk: OverlayInk;
  /**
   * Coastline/lake-outline stroke colour. `OverlayInk` (rivers/borders/roads/
   * grid/contours/currents) has no `coast` field — coastline colour lives on
   * each map style's `StylePalette.coast` (`utils/mapStyle/types.ts`) instead,
   * because today's cell-fill pass draws no separate coastline stroke; the
   * boundary IS the land/sea fill-colour edge. This overlay adds an explicit
   * stroke (the WebGL fill layer needs one), so the host passes
   * `style.palette.coast` through here.
   */
  coastColor: string;
  /**
   * Divisor for every overlay line width, replacing `qualityDpr` from
   * Map2D's offscreen path. The vector transform this overlay is drawn
   * through already carries `dpr` (see `VectorOverlay.redraw`), so dividing
   * widths by `scale` alone keeps stroke weight constant on screen at every
   * zoom level, instead of Map2D's raster-DPR-only compensation.
   */
  strokeScale: number;
}

// ---------------------------------------------------------------------------
// Relocated draw helpers (module-private in Map2D.tsx today — not exported
// there, so they cannot be imported without editing Map2D.tsx, which this
// task is scoped to avoid). Exported here, byte-faithful to Map2D's current
// logic, so Task 6 can delete Map2D's copies and import these instead rather
// than the two files drifting.
// ---------------------------------------------------------------------------

/**
 * Rivers, in projected space. Ported unchanged from `Map2D.tsx`'s
 * module-private `drawRiverPaths` (components/Map2D.tsx:43-71).
 */
export const drawRiverPaths = (
  ctx: CanvasRenderingContext2D,
  rivers: Point[][],
  project: (p: [number, number]) => [number, number] | null,
  lineWidth: number,
  ink: OverlayInk,
): void => {
  ctx.strokeStyle = ink.river;
  ctx.lineWidth = lineWidth;
  ctx.globalAlpha = 0.8;
  for (const path of rivers) {
    if (path.length < 2) continue;
    ctx.beginPath();
    let lastLon: number | null = null;
    path.forEach((p, i) => {
      const lon = Math.atan2(p.z, p.x) * (180 / Math.PI);
      const lat = Math.asin(Math.max(-1, Math.min(1, p.y))) * (180 / Math.PI);
      // Antimeridian crossing: break the subpath rather than draw across the map.
      const isJump = lastLon !== null && Math.abs(lon - lastLon) > 180;
      const pt = project([lon, lat]);
      if (pt) {
        if (i === 0 || isJump) ctx.moveTo(pt[0], pt[1]);
        else ctx.lineTo(pt[0], pt[1]);
      }
      lastLon = lon;
    });
    ctx.stroke();
  }
  ctx.globalAlpha = 1;
};

/**
 * Roads/sea-routes, in projected space. Ported unchanged from `Map2D.tsx`'s
 * inline routes block (components/Map2D.tsx:772-797), just parameterised.
 */
export const drawRoutePaths = (
  ctx: CanvasRenderingContext2D,
  routes: RouteData[],
  projection: d3.GeoProjection,
  roadLineWidth: number,
  seaRouteLineWidth: number,
  dash: [number, number],
): void => {
  ctx.globalAlpha = 0.9;
  routes.forEach(route => {
    if (route.path.length < 2) return;
    ctx.strokeStyle = route.kind === 'road' ? '#c8a25a' : '#5eb8c8';
    ctx.lineWidth = route.kind === 'road' ? roadLineWidth : seaRouteLineWidth;
    ctx.setLineDash(route.kind === 'searoute' ? dash : []);
    ctx.beginPath();
    let lastLon: number | null = null;
    route.path.forEach((p, i) => {
      const lon = Math.atan2(p.z, p.x) * (180 / Math.PI);
      const lat = Math.asin(Math.max(-1, Math.min(1, p.y))) * (180 / Math.PI);
      const isJump = lastLon !== null && Math.abs(lon - lastLon) > 180;
      const pt = projection([lon, lat]);
      if (pt) {
        if (i === 0 || isJump) ctx.moveTo(pt[0], pt[1]);
        else ctx.lineTo(pt[0], pt[1]);
      }
      lastLon = lon;
    });
    ctx.stroke();
  });
  ctx.setLineDash([]);
  ctx.globalAlpha = 1.0;
};

/**
 * Small amber diamond per marker. Ported unchanged from `Map2D.tsx`'s
 * module-private `drawMarkerPins` (components/Map2D.tsx:108-131).
 */
export const drawMarkerPins = (
  ctx: CanvasRenderingContext2D,
  markers: MarkerData[],
  project: (position: { x: number; y: number; z: number }) => [number, number] | null,
  halfSize: number,
): void => {
  if (markers.length === 0) return;
  ctx.save();
  ctx.fillStyle = '#f59e0b';
  ctx.strokeStyle = 'rgba(28, 18, 7, 0.9)';
  ctx.lineWidth = Math.max(0.75, halfSize * 0.3);
  for (const marker of markers) {
    const projected = project(marker.position);
    if (!projected) continue;
    const [x, y] = projected;
    ctx.beginPath();
    ctx.moveTo(x, y - halfSize);
    ctx.lineTo(x + halfSize, y);
    ctx.lineTo(x, y + halfSize);
    ctx.lineTo(x - halfSize, y);
    ctx.closePath();
    ctx.fill();
    ctx.stroke();
  }
  ctx.restore();
};

/**
 * Strokes the three cached boundary Path2Ds (coast, lakes, then borders).
 * These paths are already in projected CSS-px space (built by
 * `buildMapGeometryCache` in `utils/mapCache.ts`), so this is a plain
 * `ctx.stroke(path)` per boundary — no live d3 re-projection, which is what
 * makes the boundary pass cheap enough to run every frame.
 *
 * Borders get Map2D's existing double-stroke treatment (a wide casing stroke
 * under a thinner, possibly-dashed border stroke — see the module-private
 * `drawFactionBorders` in components/Map2D.tsx:74-98) so a border reads
 * clearly over any fill colour underneath it.
 */
export const strokeBoundaryPaths = (
  ctx: CanvasRenderingContext2D,
  paths: VectorOverlayBoundaryPaths,
  coastColor: string,
  coastLineWidth: number,
  lakeLineWidth: number,
  borderLineWidth: number,
  ink: OverlayInk,
): void => {
  ctx.save();
  ctx.lineJoin = 'round';
  ctx.lineCap = 'round';

  ctx.setLineDash([]);
  ctx.globalAlpha = 1;
  ctx.strokeStyle = coastColor;
  ctx.lineWidth = coastLineWidth;
  ctx.stroke(paths.coast);

  ctx.lineWidth = lakeLineWidth;
  ctx.stroke(paths.lakes);

  ctx.globalAlpha = 0.95;
  ctx.strokeStyle = ink.borderCasing;
  ctx.lineWidth = borderLineWidth * ink.borderWidthScale * 2.5;
  ctx.setLineDash([]);
  ctx.stroke(paths.borders);
  ctx.strokeStyle = ink.border;
  ctx.lineWidth = borderLineWidth * ink.borderWidthScale;
  ctx.setLineDash(ink.borderDash);
  ctx.stroke(paths.borders);
  ctx.setLineDash([]);

  ctx.restore();
};

// ---------------------------------------------------------------------------
// The sequence itself.
// ---------------------------------------------------------------------------

/**
 * Paints the full overlay draw sequence into `ctx`, which the caller
 * (`VectorOverlay.redraw`) has ALREADY set to the composed transform
 * (mirror ∘ dpr ∘ pan/zoom). This function never touches `ctx`'s transform —
 * it only issues draw calls, in the same order Map2D's non-Dymaxion main
 * render path uses today (components/Map2D.tsx:756-894), with the boundary
 * pass (coast/lakes/borders) prepended as the resolution-independent,
 * LOD-swappable replacement for the coloured-fill coastline edge and for
 * Map2D's live-recomputed `drawFactionBorders` call — see
 * `VectorOverlayBoundaryPaths` above for why.
 */
export const paintVectorOverlay = (
  ctx: CanvasRenderingContext2D,
  params: PaintVectorOverlayParams,
): void => {
  const {
    world, projection, size, toggles, boundaryPaths, contourSegments,
    mapLabels, labelVisibility, labelTheme, overlayInk, coastColor, strokeScale,
  } = params;
  const s = strokeScale > 0 ? strokeScale : 1;

  const pathGenerator = d3.geoPath(projection, ctx);

  // 1. Boundaries — coast, lakes, borders (cached Path2Ds, LOD-swappable).
  strokeBoundaryPaths(
    ctx,
    boundaryPaths,
    coastColor,
    Math.max(1, 1.5 / s),
    Math.max(0.75, 1.2 / s),
    Math.max(1, 2 / s),
    overlayInk,
  );

  // 2. Contours.
  drawContourPaths(ctx, pathGenerator, contourSegments, Math.max(0.75, 1.5 / s));

  // 3. Graticule / grid.
  if (toggles.showGrid) {
    ctx.strokeStyle = 'rgba(255,255,255,0.15)';
    ctx.lineWidth = Math.max(0.5, 1 / s);
    ctx.beginPath();
    pathGenerator(d3.geoGraticule().step([10, 10])());
    ctx.stroke();
  }

  // 4. Rivers.
  if (toggles.showRivers && world.rivers) {
    drawRiverPaths(ctx, world.rivers, projection, 1.5 / s, overlayInk);
  }

  // 5. Routes (roads / sea routes).
  if (toggles.showRoutes && world.routes) {
    drawRoutePaths(ctx, world.routes, projection, 1.4 / s, 1.2 / s, [5 / s, 4 / s]);
  }

  // 6. Ocean currents.
  if (toggles.showCurrents && world.currents) {
    drawCurrents2D(ctx, world, projection, s);
  }

  // markerProject applies its own manual mirror (size.width - x), so both the
  // marker and label passes below must run under an UN-mirrored transform,
  // same as Map2D's marker/label passes (components/Map2D.tsx:864-894) — the
  // ambient ctx transform here still carries the caller's mirror (applied
  // before `paintVectorOverlay` was called), which would double it up.
  const markerProject = (position: { x: number; y: number; z: number }): [number, number] | null => {
    const projected = projection(toLonLat([position.x, position.y, position.z]));
    return projected ? [size.width - projected[0], projected[1]] : null;
  };

  // 7. Markers.
  if (labelVisibility.markers && world.markers) {
    ctx.save();
    unmirror(ctx, size);
    drawMarkerPins(ctx, world.markers, markerProject, Math.max(1.5, 2.5 / s));
    ctx.restore();
  }

  // 8. Labels.
  if (mapLabels.length > 0) {
    ctx.save();
    unmirror(ctx, size);
    drawMapLabels(ctx, mapLabels, markerProject, s, labelVisibility, { theme: labelTheme });
    ctx.restore();
  }
};

// Cancels the mirror the caller applied before invoking `paintVectorOverlay`
// (`translate(size.width,0); scale(-1,1)`), leaving the dpr/pan/zoom part of
// the ambient transform untouched. Re-applying the SAME translate+scale a
// second time is the cancellation: mirroring twice about the same axis is
// the identity (mirror(mirror(x)) = width - (width - x) = x), so this needs
// no access to the caller's raw dpr/scale/offset values.
const unmirror = (ctx: CanvasRenderingContext2D, size: { width: number; height: number }): void => {
  ctx.translate(size.width, 0);
  ctx.scale(-1, 1);
};
