/**
 * How a map style letters its labels.
 *
 * Separate from `StylePalette` because labels are drawn by `utils/labels.ts`,
 * which knows nothing about passes or substrates — it is handed a theme rather
 * than reaching into the style registry, so there is no import cycle.
 */
export interface LabelTheme {
  /** Canvas font stack for territorial and settlement names. */
  displayFamily: string;
  /** Canvas font stack for geographic names. */
  geoFamily: string;
  /** Default text fill. */
  ink: string;
  /** Water names — seas, oceans, lakes. Cartographically these run apart. */
  waterInk: string;
  /** User markers. */
  markerInk: string;
  /** Halo drawn behind the text so it knocks a hole in whatever it sits on. */
  halo: string;
  /** Halo width as a fraction of font size. */
  haloScale: number;
  /**
   * Tracking on uppercase territorial names, in em. Real maps letterspace a
   * country name across its territory; 0 keeps the browser default.
   */
  trackingEm: number;
  /**
   * Faces to await before drawing. Canvas2D does NOT wait for webfonts — it
   * silently falls back and draws anyway — so these are loaded explicitly.
   * Format is the CSS shorthand `document.fonts.load` expects.
   */
  faces: string[];
}

/** The app's own dark-theme lettering. Unchanged from before styles existed. */
export const DEFAULT_LABEL_THEME: LabelTheme = {
  displayFamily: 'Inter, ui-sans-serif, system-ui, sans-serif',
  geoFamily: 'Inter, ui-sans-serif, system-ui, sans-serif',
  ink: '#f8fafc',
  waterInk: '#dbeafe',
  markerInk: '#fde68a',
  halo: 'rgba(2, 6, 23, 0.95)',
  haloScale: 0.28,
  trackingEm: 0,
  faces: [],
};

/**
 * Parchment lettering: engraved Roman caps for territory, a 17th-century
 * Oxford face for geography.
 *
 * This follows the actual cartographic rule rather than picking two fonts that
 * look old — roman for land, italic for water — so the hierarchy reads even
 * before you can make out the names. Cinzel is cut from Roman inscriptional
 * capitals, which is what map engravers imitated; IM Fell English is a
 * digitisation of the types Oxford was printing with in the 1600s.
 *
 * The halo is PAPER, not a dark outline. A dark halo is right on a dark map and
 * looks like soot on a light one; paper makes the label knock a clean hole in
 * the ocean hatch, which is how a real map keeps names legible over engraving.
 */
export const PARCHMENT_LABEL_THEME: LabelTheme = {
  displayFamily: "'Cinzel', 'Times New Roman', Georgia, serif",
  geoFamily: "'IM Fell English', 'Times New Roman', Georgia, serif",
  ink: '#3b2f1c',
  // Faded to sit behind the land names without a colour change — a period map
  // has one ink, and hierarchy comes from weight and spacing, not from hue.
  waterInk: '#6a6f66',
  markerInk: '#7a3b1c',
  halo: '#efe2c0',
  // Wider than the dark theme's: it is knocking out hatching, not adding a
  // contrast edge, so it needs to actually clear the lines.
  haloScale: 0.38,
  trackingEm: 0.18,
  faces: [
    '400 16px Cinzel',
    '700 16px Cinzel',
    '400 16px "IM Fell English"',
    'italic 400 16px "IM Fell English"',
  ],
};

/**
 * Await a theme's faces so Canvas2D draws with them.
 *
 * **This is not optional and its absence is invisible.** `ctx.font` accepts a
 * family that has not loaded, falls back to the next stack entry, and draws
 * without error — so a map rendered too early looks exactly like an unstyled
 * one and nothing anywhere reports a problem. Whatever draws labels must wait
 * on this and redraw when it settles.
 *
 * Resolves (rather than rejecting) on failure: a missing font should degrade to
 * the fallback stack, never blank the map.
 */
export const ensureLabelFonts = async (theme: LabelTheme): Promise<void> => {
  if (typeof document === 'undefined' || !document.fonts || theme.faces.length === 0) return;
  await Promise.all(theme.faces.map(f => document.fonts.load(f).catch(() => [])));
};

/**
 * Blueprint lettering: a draughtsman's hand, which on a drawing board means a
 * monospaced technical face and nothing else. One family for everything —
 * a drawing office does not switch typefaces between a mountain range and a
 * country, it changes size and spacing.
 *
 * `faces` is EMPTY on purpose. The stack is the system's own monospace, so
 * there is nothing to await and nothing to embed: an exported SVG will letter
 * itself in whatever monospace the reader's machine has, which is the correct
 * trade for a look that does not depend on one specific face.
 */
export const BLUEPRINT_LABEL_THEME: LabelTheme = {
  displayFamily: "ui-monospace, 'SF Mono', SFMono-Regular, Menlo, Consolas, monospace",
  geoFamily: "ui-monospace, 'SF Mono', SFMono-Regular, Menlo, Consolas, monospace",
  ink: '#e6f0f8',
  waterInk: '#9dc2de',
  markerInk: '#ff9e6b',
  // The paper tone, so a name knocks a clean hole in the hatch and the
  // graticule rather than gaining a dark outline over them.
  halo: '#1a5183',
  haloScale: 0.36,
  trackingEm: 0.16,
  faces: [],
};

/**
 * Ink & wash lettering: IM Fell English throughout, no display face.
 *
 * Parchment spends its Roman caps on territory; this look has no engraved
 * plate to imitate, so it letters everything in the one text face and leans on
 * the roman/italic split alone — upright for land, italic for water.
 *
 * Both faces are already loaded and already in `EMBEDDABLE_FACES`, so this
 * theme costs the bundle nothing and exports carry it.
 */
export const INKWASH_LABEL_THEME: LabelTheme = {
  displayFamily: "'IM Fell English', 'Times New Roman', Georgia, serif",
  geoFamily: "'IM Fell English', 'Times New Roman', Georgia, serif",
  ink: '#22221f',
  waterInk: '#6f7472',
  markerInk: '#8c3b1c',
  halo: '#f4f1ea',
  haloScale: 0.34,
  trackingEm: 0.14,
  faces: [
    '400 16px "IM Fell English"',
    'italic 400 16px "IM Fell English"',
  ],
};

/**
 * Boardgame lettering: the system sans, set heavy and wide.
 *
 * A printed board is read from a metre away across a table, so its names are
 * short, uppercase and generously tracked. `faces` is empty for the same reason
 * as Blueprint's — the stack is the system's own.
 */
export const BOARDGAME_LABEL_THEME: LabelTheme = {
  displayFamily: 'Inter, ui-sans-serif, system-ui, sans-serif',
  geoFamily: 'Inter, ui-sans-serif, system-ui, sans-serif',
  ink: '#1a1a18',
  waterInk: '#2f4a5c',
  markerInk: '#8c2f1c',
  halo: '#f2e9d5',
  // Wide: these sit on saturated territory fills, and a thin halo leaves the
  // letterforms fighting the colour underneath.
  haloScale: 0.42,
  trackingEm: 0.2,
  faces: [],
};
