import { ViewMode } from '../../types';

/** How a style paints the land fill for a given view mode. See spec §4. */
export type FillPolicy =
  | 'bare'         // no fill; glyphs, coastline and hillshade carry the terrain
  | 'categorical'  // bare paper plus a muted categorical fill on top
  | 'ramp';        // keep the mode's own continuous fill; suppress glyphs

export type MapStyleId = 'default' | 'parchment';

export type GlyphKind = 'mountain' | 'hill' | 'forest' | 'conifer' | 'dune' | 'marsh';

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
}

export interface MapStyle {
  id: MapStyleId;
  name: string;
  palette: StylePalette;
  /** Per-view-mode land fill rule. See spec §4. */
  fillPolicy: (mode: ViewMode) => FillPolicy;
}
