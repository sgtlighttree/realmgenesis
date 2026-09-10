import { DUNE, FOREST, MOUNTAIN, UnitPoint, UnitSubpath } from './mapStyle/glyphPaths';

/**
 * Small symbols drawn beside a map label so a settlement tier or a natural
 * feature is identifiable without reading the name.
 *
 * Same convention as `mapStyle/glyphPaths.ts`, deliberately: open polylines in a
 * unit box, stroked and never filled. `range`/`forest`/`desert` REFERENCE that
 * module's shapes rather than restating them — A3's relief pass and a range
 * label must not draw two different mountains.
 *
 * The difference from the relief glyphs is normalization. A relief glyph is
 * placed at an authored scale; a label glyph has to sit next to text of a known
 * size, so every shape here is fitted into a 1x1 box (long axis touching, aspect
 * preserved) at module load. `drawMapLabels` then multiplies by `fontSize`, and
 * grows its declutter rectangle by `labelGlyphMetrics` — a glyph drawn outside
 * that rectangle would overlap neighbouring labels invisibly to the declutterer.
 */

// --- Authored shapes -------------------------------------------------------
//
// The settlement ladder reads by COUNT and silhouette, not by size: every shape
// is normalized to the same box, so a "smaller" hamlet house would render at
// exactly the village's size and the two would be indistinguishable. Three
// towers -> two towers -> one house -> a bare ring, with the capital's star
// standing outside the series entirely.

/** Five-pointed star, outer r = 1, inner r = 0.45, first point up. */
const CAPITAL: UnitSubpath[] = [
  [
    [0, -1], [0.264, -0.364], [0.951, -0.309], [0.428, 0.139], [0.588, 0.809],
    [0, 0.45], [-0.588, 0.809], [-0.428, 0.139], [-0.951, -0.309], [-0.264, -0.364],
    [0, -1],
  ],
];

const CITY: UnitSubpath[] = [
  [[-1, 0.6], [1, 0.6]],
  [[-0.9, 0.6], [-0.9, -0.2], [-0.45, -0.2], [-0.45, 0.6]],
  [[-0.9, -0.2], [-0.675, -0.55], [-0.45, -0.2]],
  [[-0.25, 0.6], [-0.25, -0.55], [0.25, -0.55], [0.25, 0.6]],
  [[-0.25, -0.55], [0, -0.95], [0.25, -0.55]],
  [[0.45, 0.6], [0.45, -0.2], [0.9, -0.2], [0.9, 0.6]],
  [[0.45, -0.2], [0.675, -0.55], [0.9, -0.2]],
];

const TOWN: UnitSubpath[] = [
  [[-0.8, 0.6], [0.8, 0.6]],
  [[-0.7, 0.6], [-0.7, -0.35], [-0.2, -0.35], [-0.2, 0.6]],
  [[-0.7, -0.35], [-0.45, -0.75], [-0.2, -0.35]],
  [[0.2, 0.6], [0.2, -0.05], [0.7, -0.05], [0.7, 0.6]],
  [[0.2, -0.05], [0.45, -0.4], [0.7, -0.05]],
];

const VILLAGE: UnitSubpath[] = [
  [[-0.5, 0.6], [-0.5, -0.1], [0.5, -0.1], [0.5, 0.6], [-0.5, 0.6]],
  [[-0.72, -0.1], [0, -0.7], [0.72, -0.1]],
];

/** A bare ring: the end of the ladder, no built form left to draw. */
const HAMLET: UnitSubpath[] = [
  [
    [0.45, 0], [0.318, 0.318], [0, 0.45], [-0.318, 0.318], [-0.45, 0],
    [-0.318, -0.318], [0, -0.45], [0.318, -0.318], [0.45, 0],
  ],
];

/** Open water: two crests. Serves `sea` and `ocean`, which differ only in size. */
const WAVES: UnitSubpath[] = [
  [[-0.9, -0.25], [-0.45, -0.6], [0, -0.25], [0.45, -0.6], [0.9, -0.25]],
  [[-0.9, 0.35], [-0.45, 0], [0, 0.35], [0.45, 0], [0.9, 0.35]],
];

/** Enclosed water, so a lake reads as bounded where a sea runs off the edge. */
const LAKE: UnitSubpath[] = [
  [
    [0.9, 0], [0.64, 0.42], [0, 0.6], [-0.64, 0.42], [-0.9, 0],
    [-0.64, -0.42], [0, -0.6], [0.64, -0.42], [0.9, 0],
  ],
  [[-0.4, 0.05], [-0.13, -0.15], [0.13, 0.05], [0.4, -0.15]],
];

/** The inverse of LAKE: land bounded by water rather than water bounded by land. */
const ISLAND: UnitSubpath[] = [
  [[-1, 0.3], [-0.62, 0.3]],
  [[-0.62, 0.3], [-0.32, -0.3], [0.05, -0.55], [0.42, -0.2], [0.68, 0.3]],
  [[0.68, 0.3], [1, 0.3]],
  [[-0.6, 0.72], [-0.2, 0.5], [0.2, 0.72], [0.6, 0.5]],
];

export const LABEL_GLYPH_KINDS = [
  'capital', 'city', 'town', 'village', 'hamlet',
  'mountain', 'trees', 'dune', 'waves', 'lake', 'island',
] as const;

export type LabelGlyphKind = (typeof LABEL_GLYPH_KINDS)[number];

const RAW_SHAPES: Record<LabelGlyphKind, UnitSubpath[]> = {
  capital: CAPITAL,
  city: CITY,
  town: TOWN,
  village: VILLAGE,
  hamlet: HAMLET,
  // Borrowed from the A3 relief set — see the module header.
  mountain: MOUNTAIN,
  trees: FOREST,
  dune: DUNE,
  waves: WAVES,
  lake: LAKE,
  island: ISLAND,
};

// --- Normalization ---------------------------------------------------------

export interface LabelGlyphMetrics {
  /** Width as a fraction of the caller's `size`. At most 1. */
  width: number;
  /** Height as a fraction of the caller's `size`. At most 1. */
  height: number;
}

interface NormalizedGlyph {
  subpaths: UnitSubpath[];
  metrics: LabelGlyphMetrics;
}

/**
 * Fit a shape's bounding box to unit HEIGHT, centred on the origin, preserving
 * aspect. Width follows and may exceed 1 — MOUNTAIN is roughly 2:1, so it comes
 * out twice as wide as it is tall.
 *
 * Height, not the long axis: the glyph stacks above the name, so height is what
 * sets its visual weight beside the text and width is nearly free (the label's
 * box takes `max(textWidth, glyphWidth)`, and names are usually the wider of
 * the two). Normalizing the LONG axis instead would render the wide kinds —
 * mountain, dune, waves — at half the height of the square ones, so a hamlet's
 * ring would stand twice as tall as a mountain range's peak.
 */
const normalize = (subpaths: UnitSubpath[]): NormalizedGlyph => {
  let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
  for (const sub of subpaths) {
    for (const [x, y] of sub) {
      if (x < minX) minX = x;
      if (x > maxX) maxX = x;
      if (y < minY) minY = y;
      if (y > maxY) maxY = y;
    }
  }
  const w = maxX - minX;
  const h = maxY - minY;
  const scale = 1 / h;
  const cx = (minX + maxX) / 2;
  const cy = (minY + maxY) / 2;

  return {
    subpaths: subpaths.map(sub =>
      sub.map(([x, y]) => [(x - cx) * scale, (y - cy) * scale] as UnitPoint)),
    metrics: { width: w * scale, height: h * scale },
  };
};

const NORMALIZED = Object.fromEntries(
  LABEL_GLYPH_KINDS.map(kind => [kind, normalize(RAW_SHAPES[kind])]),
) as Record<LabelGlyphKind, NormalizedGlyph>;

/** Shape in the unit box, ready to scale by `size` and translate to an anchor. */
export const normalizedLabelGlyph = (kind: LabelGlyphKind): UnitSubpath[] =>
  NORMALIZED[kind].subpaths;

/** Extent as a fraction of `size`, for the caller's layout and declutter maths. */
export const labelGlyphMetrics = (kind: LabelGlyphKind): LabelGlyphMetrics =>
  NORMALIZED[kind].metrics;

// --- Drawing ---------------------------------------------------------------

/**
 * Stroke the glyph centred on (cx, cy), fitted to `size`.
 *
 * Walks the points directly instead of building a `Path2D` from path data: this
 * runs per label per frame on the globe overlay, and the string round-trip would
 * allocate on every one of them. `labelGlyphPathData` keeps the string form for
 * SVG export, where allocation is irrelevant and markup is the output.
 *
 * Stroke style, width, alpha and halo are the CALLER's — the glyph inherits the
 * label's ink so the two read as one mark.
 */
export const strokeLabelGlyph = (
  ctx: CanvasRenderingContext2D,
  kind: LabelGlyphKind,
  cx: number,
  cy: number,
  size: number,
): void => {
  ctx.beginPath();
  for (const sub of NORMALIZED[kind].subpaths) {
    ctx.moveTo(cx + sub[0][0] * size, cy + sub[0][1] * size);
    for (let i = 1; i < sub.length; i++) {
      ctx.lineTo(cx + sub[i][0] * size, cy + sub[i][1] * size);
    }
  }
  ctx.stroke();
};

const round = (n: number): string => (Math.round(n * 100) / 100).toString();

/** SVG `d` for the same glyph, so the vector export cannot drift from the raster. */
export const labelGlyphPathData = (
  kind: LabelGlyphKind,
  cx: number,
  cy: number,
  size: number,
): string => {
  const parts: string[] = [];
  for (const sub of NORMALIZED[kind].subpaths) {
    const pts = sub.map(([x, y]) => `${round(cx + x * size)} ${round(cy + y * size)}`);
    parts.push(`M${pts[0]}` + pts.slice(1).map(p => `L${p}`).join(''));
  }
  return parts.join('');
};
