import { WorldData, Cell, LabelVisibility, GeoFeature } from '../types';
import { normalizeVec, Point3 } from './geo';
import { DEFAULT_LABEL_THEME, LabelTheme } from './mapStyle/labelTheme';
import { LabelGlyphKind, labelGlyphMetrics, strokeLabelGlyph } from './labelGlyphs';

export type LabelKind =
  | 'faction' | 'capital' | 'town' | 'province'
  | 'city' | 'village' | 'hamlet'
  | 'range' | 'desert' | 'forest' | 'sea' | 'ocean' | 'island' | 'lake'
  | 'marker';

export interface MapLabel {
  kind: LabelKind;
  name: string;
  position: { x: number; y: number; z: number };
  priority: number;
  factionId?: number;
}

export interface LabelStyleConfig {
  fontWeight: number;
  baseFontSize: number;
  alpha: number;
  uppercase: boolean;
  visibilityKey: keyof LabelVisibility;
  italic?: boolean; // cartographic feel for geographic (esp. water) labels
  fill?: string; // override the default light fill (water kinds run slightly blued)
  /**
   * Which of the theme's two faces letters this kind.
   *
   * `display` is territory and settlement — the political layer. `geo` is the
   * physical layer. Real maps set those in different faces, and the split is
   * what carries the hierarchy before you can read any of the names.
   */
  face?: 'display' | 'geo';
  /** Water names take the theme's `waterInk` and run apart from land names. */
  water?: boolean;
  /**
   * Symbol drawn to the left of the name, so the kind is readable before the
   * name is. Optional: territory (`faction`, `province`) has no single point to
   * symbolise, and `marker` already carries its own C4 pin.
   *
   * Lives on this table rather than in a side map because `exportVector.ts`
   * mirrors LABEL_CONFIG to build the SVG label pass — a second table is a
   * second thing to forget.
   */
  glyph?: LabelGlyphKind;
}

// Exported so vector exporters (SVG/GeoJSON) can mirror the exact per-kind
// styling used by the canvas label pass instead of duplicating the table.
export const LABEL_CONFIG: Record<LabelKind, LabelStyleConfig> = {
  faction: { fontWeight: 700, baseFontSize: 14, alpha: 1.0, uppercase: true, visibilityKey: 'factions', face: 'display' },
  capital: { fontWeight: 700, baseFontSize: 11, alpha: 0.95, uppercase: false, visibilityKey: 'capitals', face: 'display', glyph: 'capital' },
  province: { fontWeight: 400, baseFontSize: 10, alpha: 0.85, uppercase: false, visibilityKey: 'provinces', face: 'display' },
  town: { fontWeight: 400, baseFontSize: 9, alpha: 0.75, uppercase: false, visibilityKey: 'towns', face: 'display', glyph: 'town' },
  city: { fontWeight: 600, baseFontSize: 10, alpha: 0.85, uppercase: false, visibilityKey: 'cities', face: 'display', glyph: 'city' },
  village: { fontWeight: 400, baseFontSize: 8, alpha: 0.7, uppercase: false, visibilityKey: 'villages', face: 'display', glyph: 'village' },
  hamlet: { fontWeight: 400, baseFontSize: 7, alpha: 0.65, uppercase: false, visibilityKey: 'hamlets', face: 'display', glyph: 'hamlet' },
  // Geographic features (B3) — one shared visibility toggle. Water kinds italic.
  ocean: { fontWeight: 400, baseFontSize: 13, alpha: 0.9, uppercase: false, visibilityKey: 'geography', italic: true, face: 'geo', water: true, glyph: 'waves' },
  sea: { fontWeight: 400, baseFontSize: 10, alpha: 0.85, uppercase: false, visibilityKey: 'geography', italic: true, face: 'geo', water: true, glyph: 'waves' },
  lake: { fontWeight: 400, baseFontSize: 9, alpha: 0.8, uppercase: false, visibilityKey: 'geography', italic: true, face: 'geo', water: true, glyph: 'lake' },
  // Land geography stays ROMAN. Italic is the water convention; setting the
  // mountains in it too throws away the one distinction it buys.
  range: { fontWeight: 400, baseFontSize: 10, alpha: 0.85, uppercase: false, visibilityKey: 'geography', face: 'geo', glyph: 'mountain' },
  desert: { fontWeight: 400, baseFontSize: 10, alpha: 0.85, uppercase: false, visibilityKey: 'geography', face: 'geo', glyph: 'dune' },
  forest: { fontWeight: 400, baseFontSize: 10, alpha: 0.85, uppercase: false, visibilityKey: 'geography', face: 'geo', glyph: 'trees' },
  island: { fontWeight: 400, baseFontSize: 9, alpha: 0.8, uppercase: false, visibilityKey: 'geography', face: 'geo', glyph: 'island' },
  marker: { fontWeight: 600, baseFontSize: 10, alpha: 0.95, uppercase: false, visibilityKey: 'markers', face: 'display' },
};

/**
 * Map2D's pan/zoom clamp — the full range `scale` can take. Exported so the
 * thresholds below can be tested against the zoom the user can actually reach:
 * a cutoff outside this range is dead code that silently hides a whole tier.
 * The globe had exactly that bug (hamlets gated at camDist 1.15 against a
 * minDistance of 1.2) — see `components/markerLod.ts`.
 */
export const MAP2D_MIN_SCALE = 0.6;
export const MAP2D_MAX_SCALE = 6;

export const ZOOM_THRESHOLDS: Record<LabelKind, number> = {
  faction: 0,
  // At MAP2D_MIN_SCALE so a capital is never missing from a fully zoomed-out
  // map, matching the globe where the capital cutoff is the camera's maximum.
  capital: MAP2D_MIN_SCALE,
  province: 1.2,
  town: 2.0,
  city: 0.9,
  village: 2.6,
  hamlet: 3.2,
  // Oceans always; mid-scale geography from ~1.0; small features only when zoomed in.
  ocean: 0,
  sea: 1.0,
  range: 1.0,
  desert: 1.0,
  forest: 1.0,
  island: 1.8,
  lake: 1.8,
  marker: 0.7,
};

interface LabelRect {
  x: number;
  y: number;
  w: number;
  h: number;
}

const rectsOverlap = (a: LabelRect, b: LabelRect): boolean =>
  a.x < b.x + b.w && a.x + a.w > b.x && a.y < b.y + b.h && a.y + a.h > b.y;

const computeCentroid = (cells: Cell[]): Point3 => {
  let sumX = 0, sumY = 0, sumZ = 0;
  for (const cell of cells) {
    sumX += cell.center.x;
    sumY += cell.center.y;
    sumZ += cell.center.z;
  }
  return normalizeVec([sumX, sumY, sumZ]);
};

// Geographic features (B3) carry their own anchor; priorities interleave with
// the civ labels so the shared declutter pass ranks them sensibly. Oceans sit
// just below faction names; mid kinds around provinces; islands/lakes lowest.
const GEO_PRIORITY: Record<string, number> = {
  ocean: 0.5,
  sea: 2.5,
  range: 2.5,
  desert: 2.5,
  forest: 2.5,
  island: 3.5,
  lake: 3.5,
};

// Cap on how many features of one kind may produce labels, as a fraction of
// total cells. Feature COUNT scales with cell count — a 200k-cell world detects
// ~355 lakes where a 5k one finds a dozen — so a fixed cap would starve small
// worlds and a fixed fraction floods large ones. This keeps the map's label
// budget roughly constant while letting big worlds surface their big features.
const GEO_LABEL_CAP: Record<string, number> = {
  ocean: 8,
  sea: 12,
  range: 16,
  desert: 12,
  forest: 16,
  island: 12,
  lake: 16,
};

/**
 * Keeps only the largest features of each kind. `GeoFeature.size` (member cell
 * count) was already computed by the detection pass and previously ignored, so
 * a three-cell pond competed with an inland sea for the same label slot — the
 * visible cause of label flooding at high cell counts.
 *
 * Ranking is by size within kind, so relative importance is preserved no matter
 * what the absolute sizes are at a given resolution.
 */
const rankGeoFeatures = (features: GeoFeature[]): GeoFeature[] => {
  const byKind = new Map<string, GeoFeature[]>();
  for (const f of features) {
    const list = byKind.get(f.kind);
    if (list) list.push(f);
    else byKind.set(f.kind, [f]);
  }
  const kept: GeoFeature[] = [];
  for (const [kind, list] of byKind) {
    list.sort((a, b) => b.size - a.size);
    kept.push(...list.slice(0, GEO_LABEL_CAP[kind] ?? 12));
  }
  return kept;
};

export const collectLabels = (world: WorldData): MapLabel[] => {
  const labels: MapLabel[] = [];
  const seaLevel = world.params.seaLevel;

  // Geographic labels are terrain-derived — emitted even when civData is absent.
  for (const feature of rankGeoFeatures(world.features ?? [])) {
    labels.push({
      kind: feature.kind,
      name: feature.name,
      position: { x: feature.anchor.x, y: feature.anchor.y, z: feature.anchor.z },
      priority: GEO_PRIORITY[feature.kind] ?? 3,
    });
  }

  // User-placed markers are sphere-anchored, not civ/terrain-derived — emitted
  // independent of civData, same as geo features above.
  for (const marker of world.markers ?? []) {
    labels.push({
      kind: 'marker',
      name: marker.name,
      position: { x: marker.position.x, y: marker.position.y, z: marker.position.z },
      priority: 1.5,
    });
  }

  if (!world.civData) {
    labels.sort((a, b) => a.priority - b.priority);
    return labels;
  }

  const cellsByFaction = new Map<number, Cell[]>();
  const cellsByProvince = new Map<string, Cell[]>();

  for (const cell of world.cells) {
    if (cell.regionId !== undefined) {
      let bucket = cellsByFaction.get(cell.regionId);
      if (!bucket) { bucket = []; cellsByFaction.set(cell.regionId, bucket); }
      bucket.push(cell);

      if (cell.provinceId !== undefined) {
        const key = `${cell.regionId}-${cell.provinceId}`;
        let provBucket = cellsByProvince.get(key);
        if (!provBucket) { provBucket = []; cellsByProvince.set(key, provBucket); }
        provBucket.push(cell);
      }
    }
  }

  for (const faction of world.civData.factions) {
    const factionCells = cellsByFaction.get(faction.id) || [];
    const landCells = factionCells.filter(c => c.height >= seaLevel);
    const labelCells = landCells.length > 0 ? landCells : factionCells;

    if (labelCells.length > 0) {
      const [nx, ny, nz] = computeCentroid(labelCells);
      labels.push({
        kind: 'faction',
        name: faction.name,
        position: { x: nx, y: ny, z: nz },
        priority: 0,
        factionId: faction.id,
      });
    }

    for (const province of faction.provinces) {
      const key = `${faction.id}-${province.id}`;
      const provinceCells = cellsByProvince.get(key) || [];
      const provLandCells = provinceCells.filter(c => c.height >= seaLevel);
      const provLabelCells = provLandCells.length > 0 ? provLandCells : provinceCells;

      if (provLabelCells.length > 0) {
        const [nx, ny, nz] = computeCentroid(provLabelCells);
        labels.push({
          kind: 'province',
          name: province.name,
          position: { x: nx, y: ny, z: nz },
          priority: 2,
          factionId: faction.id,
        });
      }

      for (const town of province.towns) {
        const cell = world.cells[town.cellId];
        if (!cell) continue;
        const [nx, ny, nz] = normalizeVec([cell.center.x, cell.center.y, cell.center.z]);
        labels.push({
          kind: town.tier === 'capital' ? 'capital' : town.tier, // tiers share kind names
          name: town.name,
          position: { x: nx, y: ny, z: nz },
          priority: town.tier === 'capital' ? 1 : town.tier === 'city' ? 2 : town.tier === 'town' ? 3 : town.tier === 'village' ? 4 : 5,
          factionId: faction.id,
        });
      }
    }
  }

  labels.sort((a, b) => a.priority - b.priority);
  return labels;
};

/**
 * Trailing options, not more positional arguments — the same call this file
 * already learned with `getCellColor`'s ColorContext. Both fields are optional,
 * so an unstyled caller passes nothing.
 */
export interface DrawLabelsOptions {
  /** Override for high-resolution export canvases. */
  fontScale?: number;
  /** Lettering; the app's dark-theme default when omitted. */
  theme?: LabelTheme;
  /**
   * Smallest font size, in the ctx's CURRENT units. Defaults to 8.
   *
   * Exists for callers drawing under a transform that already carries the zoom:
   * the floor is meant to keep text legible in SCREEN pixels, so under a scaled
   * transform it has to be divided by that scale like `fontScale` is, or it
   * clamps every label back up as you zoom in and the sizing stops being
   * constant. Pass `8 / scale` alongside `fontScale: 1 / scale`.
   */
  minFontSize?: number;
  /**
   * Draw the glyph on GEOGRAPHIC labels (range/desert/forest/water/island).
   * Callers pass `!style.paintsRelief`: A3's relief styles already stamp
   * mountains, trees and dunes across the land, and a range label's glyph would
   * land a second mountain next to one of theirs. Settlement glyphs are never
   * suppressed — no style paints towns. Defaults to true.
   */
  geoGlyphs?: boolean;
}

/**
 * Glyph size as a fraction of the label's font size, and the gap between glyph
 * and name. Tied to `fontSize` rather than to the caller's `scale`, which is a
 * sentinel on the globe path (`LABEL_SCALE_LOD_OFF`).
 *
 * TUNED, not picked. A glyph makes its label's declutter rectangle taller, and
 * the greedy declutterer pays for that in labels dropped. Measured on a dense
 * world (40k points, settlementDensity 2.0, 314 candidate labels) against the
 * same run with glyphs stripped: at the original 0.95 the fully-zoomed-out map
 * lost 26.7% of its labels. At 0.7 that worst case is 15.6%, and 6-8% at
 * working zoom. The gap is nearly free — `max(2, ...)` floors it at every font
 * size the app actually uses — so size is the only real lever here.
 */
export const GLYPH_SIZE_FACTOR = 0.7;
export const GLYPH_GAP_EM = 0.15;

/**
 * Where the glyph sits relative to the name, as offsets from the label anchor.
 * Both are pure HEIGHT arithmetic, deliberately: `exportVector.ts` builds the
 * SVG label pass with no DOM and therefore no `measureText`, so any layout that
 * needed the text WIDTH — a glyph to the left of the name, say — could not be
 * reproduced there and the vector map would drift from the raster. Stacking the
 * glyph above the name keeps one layout for both substrates, and symbol-above-
 * name is the ordinary convention for point features anyway.
 */
export const glyphLayout = (
  fontSize: number,
  glyph: LabelGlyphKind,
): { glyphSize: number; glyphDy: number; textDy: number; width: number; height: number } => {
  const glyphSize = fontSize * GLYPH_SIZE_FACTOR;
  const metrics = labelGlyphMetrics(glyph);
  const glyphH = metrics.height * glyphSize;
  const gap = Math.max(2, fontSize * GLYPH_GAP_EM);
  const textHeight = fontSize * 1.2;
  const unionH = glyphH + gap + textHeight;
  return {
    glyphSize,
    glyphDy: -unionH / 2 + glyphH / 2,
    textDy: unionH / 2 - textHeight / 2,
    width: metrics.width * glyphSize,
    height: unionH,
  };
};

export const drawMapLabels = (
  ctx: CanvasRenderingContext2D,
  labels: MapLabel[],
  project: (position: { x: number; y: number; z: number }) => [number, number] | null,
  scale: number,
  visibility: LabelVisibility,
  opts: DrawLabelsOptions = {},
): void => {
  if (labels.length === 0) return;

  const theme = opts.theme ?? DEFAULT_LABEL_THEME;
  const effectiveFontScale = opts.fontScale ?? Math.min(scale, 2);
  const geoGlyphs = opts.geoGlyphs ?? true;
  const minFontSize = opts.minFontSize ?? 8;
  const placed: LabelRect[] = [];

  ctx.save();

  for (const label of labels) {
    const config = LABEL_CONFIG[label.kind];
    if (!visibility[config.visibilityKey]) continue;
    if (scale < ZOOM_THRESHOLDS[label.kind]) continue;

    const projected = project(label.position);
    if (!projected) continue;

    const fontSize = Math.max(minFontSize, config.baseFontSize * effectiveFontScale);
    const text = config.uppercase ? label.name.toUpperCase() : label.name;
    const padding = Math.max(4, fontSize * 0.3);

    const glyph = config.glyph && (geoGlyphs || config.face !== 'geo')
      ? config.glyph : undefined;

    const style = config.italic ? 'italic ' : '';
    const family = config.face === 'geo' ? theme.geoFamily : theme.displayFamily;
    ctx.font = `${style}${config.fontWeight} ${fontSize}px ${family}`;
    // Tracking on the letterspaced territorial caps only. `letterSpacing` is
    // set BEFORE measureText, or the declutter rectangle is narrower than the
    // text it is supposed to bound and names overlap. Unsupported in older
    // engines, where it is simply ignored — the labels render untracked rather
    // than not at all.
    const tracking = config.uppercase && theme.trackingEm !== 0
      ? `${(theme.trackingEm * fontSize).toFixed(2)}px`
      : '0px';
    ctx.letterSpacing = tracking;
    const textWidth = ctx.measureText(text).width;
    const textHeight = fontSize * 1.2;

    // The glyph is laid out INSIDE the declutter rectangle, not beside it. A
    // rectangle sized from text metrics alone would let the glyph overlap the
    // neighbouring label while the declutterer saw no collision — it never
    // measured the glyph, so the overlap is invisible to it.
    const layout = glyph ? glyphLayout(fontSize, glyph) : null;
    const totalWidth = Math.max(textWidth, layout ? layout.width : 0);
    const totalHeight = layout ? layout.height : textHeight;
    const textY = projected[1] + (layout ? layout.textDy : 0);

    const rect: LabelRect = {
      x: projected[0] - totalWidth / 2 - padding,
      y: projected[1] - totalHeight / 2 - padding,
      w: totalWidth + padding * 2,
      h: totalHeight + padding * 2,
    };

    if (placed.some(r => rectsOverlap(rect, r))) continue;
    placed.push(rect);

    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';
    ctx.lineJoin = 'round';
    ctx.lineWidth = Math.max(3, fontSize * theme.haloScale);
    ctx.strokeStyle = theme.halo;
    ctx.fillStyle = config.fill
      ?? (label.kind === 'marker' ? theme.markerInk : config.water ? theme.waterInk : theme.ink);
    ctx.globalAlpha = config.alpha;
    if (layout && glyph) {
      const glyphY = projected[1] + layout.glyphDy;
      // These are OPEN polylines, so the cap is visible geometry, not a detail:
      // the SVG export emits stroke-linecap="round" and the default `butt`
      // would leave the raster's ends chopped where the vector's are rounded.
      ctx.lineCap = 'round';
      // Halo first at the text's own halo width, then the mark itself in the
      // label's ink, so glyph and name read as one object over any fill.
      strokeLabelGlyph(ctx, glyph, projected[0], glyphY, layout.glyphSize);
      ctx.strokeStyle = ctx.fillStyle;
      ctx.lineWidth = Math.max(1, fontSize * 0.09);
      strokeLabelGlyph(ctx, glyph, projected[0], glyphY, layout.glyphSize);
      // Restored for the text below. `save`/`restore` wraps the whole loop, not
      // each label, so anything left set here leaks into every later label too.
      ctx.lineCap = 'butt';
      ctx.strokeStyle = theme.halo;
      ctx.lineWidth = Math.max(3, fontSize * theme.haloScale);
    }

    ctx.strokeText(text, projected[0], textY);
    ctx.fillText(text, projected[0], textY);
  }

  ctx.restore();
};
