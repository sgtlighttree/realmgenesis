/**
 * How a map style inks the LINE overlays — rivers, borders, roads, sea routes,
 * the graticule, contours, the ruler and the Dymaxion cage.
 *
 * These are not painted by style passes. Rivers are drawn by an R3F screen
 * overlay tenant on the globe, by `Map2D` on the flat map, and again by each
 * export path — so before this existed the river colour was hardcoded in FOUR
 * files across FIVE call sites, and a style could only ever reach the fills
 * underneath them. That is why parchment read as "nothing happened": the paper
 * was right and everything drawn on top of it was still Tailwind sky-400.
 *
 * **Every value must be an OPAQUE colour.** The SVG export consumes these, and
 * SVG 1.1 has no `rgba()` syntax — it renders inconsistently in Illustrator and
 * Inkscape. Transparency goes through a separate opacity attribute, exactly as
 * `Substrate.fillFeature` documents.
 */
export interface OverlayInk {
  river: string;
  /** Multiplier on the base river line width. */
  riverWidthScale: number;
  border: string;
  /**
   * Wider line laid under the border so it stays readable over any fill. Dark
   * on a dark map; on paper it is the PAPER tone, which knocks the border clear
   * of the ocean hatch instead of adding a black outline to it.
   */
  borderCasing: string;
  borderWidthScale: number;
  /** Canvas dash pattern in px; empty is solid. */
  borderDash: number[];
  road: string;
  seaRoute: string;
  graticule: string;
  ruler: string;
  /** Dymaxion unfold cage. */
  cage: string;
  /**
   * Contour lines. `null` keeps the height-ramped stroke `contourStroke`
   * computes, which is the default look and carries elevation information a
   * flat ink would throw away.
   */
  contour: string | null;
  /**
   * Ocean current arrows. `null` keeps the sea-surface-temperature tint ramp —
   * again, information a flat ink would discard.
   */
  currents: string | null;
}

/** The app's dark theme. Reproduces the pre-style look exactly. */
export const DEFAULT_OVERLAY_INK: OverlayInk = {
  river: 'rgba(56,189,248,0.8)',
  riverWidthScale: 1,
  border: 'rgba(255,255,255,0.85)',
  borderCasing: 'rgba(2, 6, 23, 0.9)',
  borderWidthScale: 1,
  borderDash: [],
  road: '#c8a25a',
  seaRoute: '#5eb8c8',
  graticule: 'rgba(255,255,255,0.28)',
  ruler: '#fbbf24',
  cage: 'rgba(251,191,36,0.95)',
  contour: null,
  currents: null,
};

/**
 * One ink, several strengths — which is how an engraved map actually works.
 *
 * There is no colour coding here on purpose. A period map had a single plate of
 * ink (sometimes a second for water), so hierarchy comes from WEIGHT and from
 * whether a line is broken, not from hue. Rivers run slightly blue and thinner
 * than the coast; political borders are DASHED, because that is how a border —
 * a claim rather than a thing you can see — was distinguished from a coastline,
 * which is a fact.
 */
export const PARCHMENT_OVERLAY_INK: OverlayInk = {
  river: '#5f7480',
  // Thinner than the coastline. A river as heavy as the coast reads as a strait.
  riverWidthScale: 0.7,
  border: '#4a3c24',
  borderCasing: '#efe2c0',
  borderWidthScale: 0.9,
  borderDash: [7, 5],
  road: '#8a6a3a',
  seaRoute: '#7c8f8a',
  // inkLight: present enough to read as a ruled grid, faint enough to sit under
  // the map rather than over it.
  graticule: '#9a8355',
  ruler: '#8c3b1c',
  cage: '#8c3b1c',
  // Flat ink for both, unlike the default: the height ramp and the SST ramp are
  // modern data-viz gradients and they fight a single-plate engraving.
  contour: '#8a7550',
  currents: '#7c8f8a',
};

/**
 * Blueprint: white line work on cyanotype blue, with a red-pencil annotation.
 *
 * A drawing office marks up in one colour, so the ruler and the unfold cage —
 * both of them measurement furniture rather than map content — take the pencil,
 * and everything else stays in the drawing's own white-to-blue range.
 *
 * The graticule is deliberately the strongest it is in any style: a grid is
 * native to a technical drawing, not an optional overlay on top of one.
 */
export const BLUEPRINT_OVERLAY_INK: OverlayInk = {
  river: '#7fc4e8',
  riverWidthScale: 0.75,
  border: '#dce9f4',
  borderCasing: '#1a5183',
  borderWidthScale: 0.9,
  borderDash: [6, 4],
  road: '#9fc0da',
  seaRoute: '#5f93b8',
  graticule: '#4c81ab',
  ruler: '#ff9e6b',
  cage: '#ff9e6b',
  contour: '#5e93bd',
  currents: '#5f93b8',
};

/**
 * Ink & wash: one grey-black plate, thinner than parchment's throughout.
 *
 * A wash drawing carries its elevation in the wash, so its lines do less work
 * and run lighter. Borders are dashed for the same reason parchment's are.
 */
export const INKWASH_OVERLAY_INK: OverlayInk = {
  river: '#5c6a70',
  riverWidthScale: 0.65,
  border: '#3a3a36',
  borderCasing: '#f4f1ea',
  borderWidthScale: 0.85,
  borderDash: [6, 4],
  road: '#7a7368',
  seaRoute: '#8d9a98',
  graticule: '#b4b1a8',
  ruler: '#8c3b1c',
  cage: '#8c3b1c',
  contour: '#6f6f68',
  currents: '#8d9a98',
};

/**
 * Boardgame: heavy black rules, and a border that is SOLID rather than dashed.
 *
 * This is the one style where that is right. Everywhere else a dashed border
 * says "a claim, not a fact"; on a board the territory outline IS the fact —
 * it is the line you move pieces across — so it gets the same weight and
 * treatment as the coast.
 */
export const BOARDGAME_OVERLAY_INK: OverlayInk = {
  river: '#4a7f9e',
  riverWidthScale: 0.8,
  border: '#141414',
  borderCasing: '#f2e9d5',
  borderWidthScale: 1.6,
  borderDash: [],
  road: '#8a5a2a',
  seaRoute: '#3f6b85',
  graticule: '#c3b795',
  ruler: '#b3401c',
  cage: '#b3401c',
  contour: '#a89670',
  currents: '#3f6b85',
};
