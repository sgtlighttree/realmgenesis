import { ViewMode, WorldData } from '../../types';
import { ColorContext } from '../colors';
import { Point3 } from '../geo';
import { Substrate } from './substrate';

/** How a style paints the land fill for a given view mode. See spec §4. */
export type FillPolicy =
  | 'bare'         // no fill; glyphs, coastline and hillshade carry the terrain
  | 'categorical'  // bare paper plus a muted categorical fill on top
  | 'ramp';        // keep the mode's own continuous fill; suppress glyphs

export type MapStyleId = 'default' | 'parchment';

export type GlyphKind = 'mountain' | 'hill' | 'forest' | 'conifer' | 'dune' | 'marsh' | 'ice';

/** A glyph resolved to output pixels. Both substrates just draw this. */
export interface PlacedGlyph {
  x: number;
  y: number;
  kind: GlyphKind;
  scale: number;   // pixels, tip to base
  seedRot: number; // radians, deterministic per cell
  cellId: number;
}

export interface StylePalette {
  paper: string;
  ink: string;
  inkLight: string;
  sea: string;
  seaHatch: string;
  coast: string;
  /**
   * Hillshade tints. These exist as palette entries rather than hardcoded
   * black/white because on bare parchment a neutral grey shadow reads as dirt
   * on the paper rather than as relief — a drawn map shades in its own ink.
   */
  shadow: string;
  highlight: string;
  /**
   * Permanent ice — polar caps and frozen sea. A named entry because ice is the
   * one terrain a bare-paper style cannot leave to glyphs alone: unpainted, an
   * ice cap is indistinguishable from temperate land, and 6.8% of land in a
   * default world is ICE_CAP. Period maps render polar ice pale and stippled,
   * not blank.
   */
  ice: string;
}

/**
 * Everything a style pass needs. Assembled once per render by the caller
 * (Map2D, the PNG export, the SVG export) and handed to every pass.
 */
export interface StyleRenderContext {
  world: WorldData;
  viewMode: ViewMode;
  widthPx: number;
  heightPx: number;
  /** Pre-computed by placeGlyphs; empty when the fill policy suppresses glyphs. */
  glyphs: PlacedGlyph[];
  shadeMap: Float32Array | null;
  /**
   * Multiplier on line weights, hatch density and glyph size, on top of the
   * widthPx scaling. 1 suits a view showing the WHOLE world; the globe shows one
   * hemisphere at a time, so weights tuned for a flat map read as heavy blobs
   * there and it bakes at a lower value.
   */
  lineScale?: number;
  colorCtx: ColorContext;
  coastlines: Array<[Point3, Point3]>;
}

export type StylePass = (ctx: StyleRenderContext, sub: Substrate) => void;

export interface MapStyle {
  id: MapStyleId;
  name: string;
  palette: StylePalette;
  /** Per-view-mode land fill rule. See spec §4. */
  fillPolicy: (mode: ViewMode) => FillPolicy;
  /**
   * Ordered draw passes. EMPTY means "this style draws nothing" — that is the
   * single test for whether a render path should run its own legacy fill loop
   * instead. Do not compare against `id`; one invariant, one place.
   */
  passes: StylePass[];
}
