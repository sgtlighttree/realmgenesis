import { Point3 } from '../geo';
import { PlacedGlyph } from './types';

/** A d3 GeoJSON feature, typed structurally so this module needs no d3 import. */
export type GeoFeatureLike = { type: string; geometry: unknown };

export interface HatchSpec {
  color: string;
  spacingPx: number;
  widthPx: number;
  angleDeg: number;
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
  /** Hatch ONE feature. Full-bleed hatching would cover bare-paper land. */
  hatchFeature(feature: GeoFeatureLike, spec: HatchSpec): void;
  grain(spec: GrainSpec): void;
  drawGlyph(g: PlacedGlyph, ink: string, width: number): void;
}
