import { ViewMode, WorldData } from '../../types';
import { ColorContext } from '../colors';
import { LabelTheme } from './labelTheme';
import { OverlayInk } from './overlayInk';
import { Point3 } from '../geo';
import { Substrate } from './substrate';

/** How a style paints the land fill for a given view mode. See spec §4. */
export type FillPolicy =
  | 'bare'         // no fill; glyphs, coastline and hillshade carry the terrain
  | 'categorical'  // bare paper plus a muted categorical fill on top
  | 'ramp';        // keep the mode's own continuous fill; suppress glyphs

export type MapStyleId = 'default' | 'parchment' | 'blueprint';

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
  /**
   * What fills the viewport OUTSIDE the map — the desk the paper lies on. A
   * projection rarely fills its canvas, and raw black in that margin reads as a
   * rendering failure rather than a choice.
   *
   * 2D only, screen space, and painted by EVERY style including `default` —
   * see `paintDesk` in `desk.ts`, which is called from Map2D's display effect
   * rather than from a pass, and never by an export.
   */
  desk: string;
  /** Shadow hugging the map's edge, so the paper rests rather than floats. */
  deskShadow: string;
  /**
   * Optional vignette colour at the corners of the viewport, blending out of
   * `desk` at the centre. A lit surface rather than a flat swatch. Omit for a
   * plain fill.
   */
  deskEdge?: string;
  /**
   * Optional repeating image tiled over `desk` — wood, leather, cloth. A URL
   * (or a data URI) for an asset in `public/`. It must be TILEABLE: it is drawn
   * as a `repeat` pattern in device pixels.
   *
   * Load it through `useDeskTexture`, never with a bare `new Image()`: an
   * undecoded image draws nothing and reports nothing, which presents as "the
   * texture did nothing" — the same trap as the webfonts.
   */
  deskTexture?: string;
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
  /**
   * `[startDeg, endDeg]` latitude band over which the ocean hatch fades to
   * nothing, or undefined for no fade. Set by the GLOBE BAKE only.
   *
   * Equirectangular content converges at the poles, and no texture coordinate
   * scheme changes that. A fixed-frequency screen pattern is what turns the
   * convergence into a visible defect: the ocean hatch wound round the
   * singularity as a spiral rosette. Fading it out leaves plain paper there,
   * which reads as paper. Coastlines and fills still pinch, which is
   * geographically honest.
   *
   * A flat map has no singularity, so `Map2D` and the exports leave this unset
   * and keep hatching to the pole.
   */
  polarHatchFadeDeg?: [number, number];
  /**
   * True when the output is wrapped end to end — the globe bake, whose left and
   * right edges meet at the antimeridian.
   *
   * A hatch is drawn in output pixels, so unless its pattern repeats a whole
   * number of times across the width, its phase differs either side of the join
   * and the diagonals visibly jog down the seam. A flat map has no join, so
   * `Map2D` and the exports leave this unset.
   */
  wrapsHorizontally?: boolean;
  /**
   * DEBUG. Ocean hatch angle in degrees; 45 when unset.
   *
   * Set it to 0 to get horizontal lines. A diagonal hatch hides a vertical
   * offset or a skew inside its own slope — horizontal lines are a ruled grid,
   * and any misregistration between the baked texture and the cell mesh reads
   * off them at a glance. That is how the flipY inversion should have been
   * caught, instead of surviving three sessions of UV fixes.
   */
  hatchAngleDeg?: number;
  colorCtx: ColorContext;
  coastlines: Array<[Point3, Point3]>;
}

export type StylePass = (ctx: StyleRenderContext, sub: Substrate) => void;

export interface MapStyle {
  id: MapStyleId;
  name: string;
  palette: StylePalette;
  /**
   * Lettering. Separate from `palette` because labels are drawn by
   * `utils/labels.ts`, outside the pass/substrate machinery entirely.
   */
  labelTheme: LabelTheme;
  /**
   * Line overlays — rivers, borders, roads, graticule, contours. Drawn outside
   * the pass machinery by three different renderers, so like `labelTheme` this
   * is handed to them rather than reached for.
   */
  overlayInk: OverlayInk;
  /** Per-view-mode land fill rule. See spec §4. */
  fillPolicy: (mode: ViewMode) => FillPolicy;
  /**
   * Ordered draw passes. EMPTY means "this style draws nothing" — that is the
   * single test for whether a render path should run its own legacy fill loop
   * instead. Do not compare against `id`; one invariant, one place.
   */
  passes: StylePass[];
}
