import * as d3 from 'd3';
import type { GeoProjection } from 'd3-geo';

import { GeoJsonFeature, WorldData } from '../types';
import { Point3, toLonLat } from './geo';
import {
  chainSegments, computeCoastlineSegments, computeFactionBorderSegments,
  computeLakeOutlineSegments,
} from './boundaries';
import { simplifyPolyline } from './simplify';

export interface MapGeometryCache {
  cellPaths: Path2D[];
  // Flat [x0,y0,x1,y1,...] CSS-px, un-flipped. This is a PARITY buffer, not a
  // renderable one: a cell whose ring crosses the antimeridian is clipped by
  // d3.geoPath into 2+ disjoint closed subpaths (see the capture note below),
  // and cellVerts concatenates every subpath's points in emission order with
  // no boundary marker between them. Feeding it through a naive
  // moveTo(first)+lineTo(rest) draws a spurious line straight across the map
  // for those cells. Only `cellPaths` (built from the SAME capture, but
  // preserving each subpath's own moveTo/closePath) is safe to render as-is;
  // a future consumer that wants to redraw straight from cellVerts must first
  // add a subpath-boundary index of its own.
  cellVerts: Float32Array;
  cellOffsets: Uint32Array; // cell i's verts are [offsets[i], offsets[i+1])
  cellSubStart: Uint32Array;   // length n+1: cell i's sub-ring entries are [cellSubStart[i], cellSubStart[i+1])
  cellSubOffsets: Uint32Array; // per sub-ring: start vertex index into cellVerts (vertex units)
  coast: Path2D;
  borders: Path2D;
  lakes: Path2D;
}

// Capture what a d3.geoPath trace of one feature actually emits, split into
// subpaths (one per moveTo, closed if geoPath called closePath on it). Using
// the real d3.geoPath algorithm — not a per-vertex projection loop — is
// required for correctness, not just to match the parity test: a Voronoi
// mesh over a full sphere always has cells straddling the antimeridian and
// (for unbounded cylindrical projections like Mercator) cells enclosing a
// pole. d3's clip stream splits the former into disjoint closed pieces along
// the map edge and can legitimately emit NOTHING for the latter — Map2D
// already draws cells this way today (see `d3.geoPath(projection, ctx)` at
// components/Map2D.tsx:486/719/967), so the cache must reproduce it exactly
// rather than "fixing" it into a different (uncached) visual result.
interface CapturedSubpath {
  points: [number, number][];
  closed: boolean;
}

const captureFeatureSubpaths = (
  feature: GeoJsonFeature,
  geoPathGen: d3.GeoPath<unknown, d3.GeoPermissibleObjects>,
): CapturedSubpath[] => {
  const subpaths: CapturedSubpath[] = [];
  let current: CapturedSubpath | null = null;
  const ctx: Partial<CanvasRenderingContext2D> = {
    beginPath() {},
    moveTo(x: number, y: number) {
      current = { points: [[x, y]], closed: false };
      subpaths.push(current);
    },
    lineTo(x: number, y: number) {
      current?.points.push([x, y]);
    },
    closePath() {
      if (current) current.closed = true;
    },
    arc() {},
  };
  // Re-target the same generator's context per call rather than constructing
  // a new d3.geoPath per cell — the whole point of this cache is to pay this
  // setup cost once, not once per cell.
  geoPathGen.context(ctx as CanvasRenderingContext2D);
  geoPathGen(feature as unknown as d3.GeoPermissibleObjects);
  return subpaths;
};

// Split a polyline into wrap-free runs wherever consecutive points' raw
// longitude jumps by more than 180 degrees (the antimeridian-crossing
// heuristic already used for rivers — components/Map2D.tsx:104-116 —
// reused here for chained coastline/border/lake polylines instead of
// routing them through d3.geoPath, since these are welded chains of
// arbitrary length rather than single GeoJSON features).
const splitOnAntimeridian = (chain: Point3[]): [number, number][][] => {
  const runs: [number, number][][] = [];
  let run: [number, number][] = [];
  let lastLon: number | null = null;
  for (const v of chain) {
    const [lon, lat] = toLonLat(v);
    if (lastLon !== null && Math.abs(lon - lastLon) > 180) {
      if (run.length > 0) runs.push(run);
      run = [];
    }
    run.push([lon, lat]);
    lastLon = lon;
  }
  if (run.length > 0) runs.push(run);
  return runs;
};

// Project a chained boundary into one Path2D of simplified polylines.
// `cellCount` is threaded through to chainSegments, which uses it to scale
// the vertex-welding threshold to point density (see utils/boundaries.ts) —
// there is no argument-less form any more.
const chainedToPath = (
  segments: Array<[Point3, Point3]>,
  cellCount: number,
  projection: GeoProjection,
  tolPx: number,
): Path2D => {
  const path = new Path2D();
  for (const chain of chainSegments(segments, cellCount)) {
    // Split BEFORE projecting/simplifying: simplifyPolyline measures
    // perpendicular distance in projected space, and a run that still
    // straddles the wrap point would have that distance computed against a
    // segment spanning the entire map width, corrupting which points it
    // keeps or drops.
    for (const run of splitOnAntimeridian(chain)) {
      let pts: [number, number][] = [];
      for (const lonLat of run) {
        const p = projection(lonLat);
        if (p && Number.isFinite(p[0])) pts.push([p[0], p[1]]);
      }
      if (tolPx > 0 && pts.length > 2) pts = simplifyPolyline(pts, tolPx);
      if (pts.length < 2) continue;
      path.moveTo(pts[0][0], pts[0][1]);
      for (let k = 1; k < pts.length; k++) path.lineTo(pts[k][0], pts[k][1]);
    }
  }
  return path;
};

// Project all cell rings ONCE to CSS-px (un-flipped) + build a Path2D per cell.
// Coordinates are filled through Map2D's flipped/DPR context, so the mirror and
// DPR are applied at draw, never baked here (see the geometry-cache spec §2).
// Cell rings are projected but NOT simplified (only boundary polylines are,
// via simplifyPolyline) — the parity contract requires cached cell vertices to
// match a fresh d3.geoPath trace to sub-pixel precision.
export const buildMapGeometryCache = (
  world: WorldData,
  projection: GeoProjection,
  width: number,
  height: number,
  opts?: { simplifyTolerancePx?: number },
): MapGeometryCache => {
  void width; void height; // projection is expected to already be fit to [width,height]
  const tolPx = opts?.simplifyTolerancePx ?? 0.5;
  const n = world.cells.length;
  const features = world.geoJson?.features ?? [];

  // One geoPath generator, reused (via .context()) across every cell — never
  // mutate the caller's `projection` itself (no .precision()/.clipExtent()),
  // only this cache-local generator's context.
  const geoPathGen = d3.geoPath(projection as unknown as d3.GeoProjection);

  const perCell: CapturedSubpath[][] = new Array(n);
  const cellOffsets = new Uint32Array(n + 1);
  const cellSubStart = new Uint32Array(n + 1);
  let totalVerts = 0;
  let totalSubs = 0;
  for (let i = 0; i < n; i++) {
    const feature = features[i];
    const subpaths = feature ? captureFeatureSubpaths(feature, geoPathGen) : [];
    perCell[i] = subpaths;
    cellOffsets[i] = totalVerts;
    cellSubStart[i] = totalSubs;
    totalSubs += subpaths.length;
    for (const sp of subpaths) totalVerts += sp.points.length;
  }
  cellOffsets[n] = totalVerts;
  cellSubStart[n] = totalSubs;

  const cellVerts = new Float32Array(totalVerts * 2);
  const cellSubOffsets = new Uint32Array(totalSubs);
  const cellPaths: Path2D[] = new Array(n);
  for (let i = 0; i < n; i++) {
    const subpaths = perCell[i];
    const path = new Path2D();
    let cursor = cellOffsets[i];
    let subCursor = cellSubStart[i];
    for (const sp of subpaths) {
      cellSubOffsets[subCursor++] = cursor; // cursor is the vertex index at sub-ring start
      sp.points.forEach(([x, y], k) => {
        cellVerts[cursor * 2] = x;
        cellVerts[cursor * 2 + 1] = y;
        cursor++;
        if (k === 0) path.moveTo(x, y);
        else path.lineTo(x, y);
      });
      if (sp.closed) path.closePath();
    }
    cellPaths[i] = path;
  }

  return {
    cellPaths, cellVerts, cellOffsets, cellSubStart, cellSubOffsets,
    coast: chainedToPath(computeCoastlineSegments(world), n, projection, tolPx),
    borders: chainedToPath(computeFactionBorderSegments(world), n, projection, tolPx),
    lakes: chainedToPath(computeLakeOutlineSegments(world), n, projection, tolPx),
  };
};
