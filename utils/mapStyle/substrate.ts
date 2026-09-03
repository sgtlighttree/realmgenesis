import { Point3 } from '../geo';
import { PlacedGlyph } from './types';

/** A d3 GeoJSON feature, typed structurally so this module needs no d3 import. */
export type GeoFeatureLike = { type: string; geometry: unknown };

export interface HatchSpec {
  color: string;
  spacingPx: number;
  widthPx: number;
  angleDeg: number;
  /**
   * 0-1, default 1. EXPLICIT, never baked into `color` as an `rgba()` string,
   * for the same reason `fillFeature` takes one: SVG 1.1 has no `rgba()` syntax.
   */
  opacity?: number;
}

export interface GrainSpec {
  seed: string;
  opacity: number;
  /** Feature size multiplier; 1 = fine paper tooth. */
  scale: number;
}

/**
 * The drawing surface a style pass writes to. Two implementations —
 * Canvas2D and SVG — so every parchment pass is authored ONCE.
 *
 * Deliberately narrow: passes may only do things both substrates can express
 * natively. Anything raster-only (a blur, a composite mode) would silently
 * degrade the SVG export, so it does not belong here.
 */
export interface Substrate {
  fillRect(x: number, y: number, w: number, h: number, fill: string): void;
  /**
   * `opacity` is an EXPLICIT parameter, never baked into `fill` as an `rgba()`
   * string. Canvas2D accepts `rgba()`; SVG 1.1 does not — it needs a separate
   * `fill-opacity`, and an `rgba()` fill renders inconsistently in Illustrator
   * and Inkscape. Colour strings stay opaque so both substrates agree.
   */
  fillFeature(feature: GeoFeatureLike, fill: string, opacity?: number): void;
  strokeFeature(feature: GeoFeatureLike, stroke: string, width: number): void;
  strokeSegments(segments: Array<[Point3, Point3]>, stroke: string, width: number): void;
  hatchRect(x: number, y: number, w: number, h: number, spec: HatchSpec): void;
  /**
   * Hatch a SET of features as one composite region.
   *
   * Plural, not singular, for cost. A world has 13,000-17,000 ocean cells at the
   * default 20k points, and Map2D re-renders on every viewport change. Hatching
   * per cell means one clip plus a full-canvas line sweep EACH — roughly 400
   * clipped lines per cell. Building one composite path and sweeping once turns
   * that into a single operation.
   *
   * Full-bleed `hatchRect` is not an option either: under the `bare` fill policy
   * the land is unpainted parchment, so a full-canvas hatch would sit on top of
   * it.
   */
  hatchFeatures(features: GeoFeatureLike[], spec: HatchSpec): void;
  /**
   * Hatch a SET of cells addressed by INDEX, reusing the geometry cache's
   * prebuilt per-cell Path2Ds instead of reprojecting each feature. The cached
   * equivalent of `hatchFeatures` for the hot on-screen path: F3 Phase 1 found
   * the ocean hatch's per-redraw reprojection of ~15k cells was the dominant
   * pan/zoom cost. SVG has no reprojection cost, so its impl may delegate to the
   * per-feature path.
   */
  hatchCells(indices: number[] | Uint32Array, spec: HatchSpec): void;
  grain(spec: GrainSpec): void;
  /**
   * Run `draw` clipped to the PROJECTED SPHERE — the map's own outline, which
   * is what `d3.geoPath({ type: 'Sphere' })` traces for whatever projection is
   * in force.
   *
   * Needed because a projection does not necessarily fill its canvas.
   * `d3.geoMercator` clipped at +/-85 degrees is square, so fitting it into a
   * wide viewport leaves a margin down each side. Anything painted by raw
   * canvas coordinates — paper tone, grain — spills into that margin, where
   * there is no map, and the parchment then reads as extra land off the edges
   * of the world. Reported by Matt against the Mercator view.
   *
   * A no-op in effect for equirectangular, where the sphere IS the canvas.
   */
  withSphereClip(draw: () => void): void;
  drawGlyph(g: PlacedGlyph, ink: string, width: number): void;
  /**
   * Fill a set of cells by INDEX from a prebuilt geometry cache, each with its
   * own colour (colours indexed by cell id). The batched hot-loop primitive for
   * cached rendering: Canvas2D fills cached Path2Ds (no reprojection); SVG emits
   * one <path> per cell as fillFeature does (no aggregation — batching helps only
   * the raster/GPU backends). See the F3 Phase 1 spec.
   */
  fillCells(indices: number[] | Uint32Array, colors: string[]): void;
}
