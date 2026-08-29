import * as d3 from 'd3';

import { BiomeType, Cell } from '../../types';
import { toLonLat } from '../geo';
import { GlyphKind, PlacedGlyph } from './types';

export interface GlyphOptions {
  seaLevel: number;
  /** World seed, so glyph variation is stable across re-renders. */
  seed?: string;
  /** Minimum centre-to-centre distance, in output pixels, at widthPx = 1024. */
  minSpacingPx?: number;
  /** Glyph height in output pixels, at widthPx = 1024. */
  sizePx?: number;
}

const REFERENCE_WIDTH = 1024;
const DEFAULT_SPACING_PX = 22;
const DEFAULT_SIZE_PX = 16;

/** Height above sea level, as a 0-1 fraction of the land band. */
const landFrac = (cell: Cell, seaLevel: number): number =>
  (cell.height - seaLevel) / Math.max(1e-6, 1 - seaLevel);

/**
 * Which glyph a land cell earns, or null for open ground. Relief wins over
 * vegetation: a forested mountain reads as a mountain, the way a drawn map
 * would show it.
 */
const glyphFor = (cell: Cell, seaLevel: number): GlyphKind | null => {
  const f = landFrac(cell, seaLevel);
  // Ice outranks relief: a glaciated peak reads as ice on a drawn map, and
  // without this an ice cap carries mountain glyphs and no sign of being ice.
  if (cell.biome === BiomeType.ICE_CAP) return 'ice';
  if (f > 0.45) return 'mountain';
  if (f > 0.22) return 'hill';
  switch (cell.biome) {
    case BiomeType.BOREAL_FOREST:
      return 'conifer';
    case BiomeType.TEMPERATE_FOREST:
    case BiomeType.TEMPERATE_RAINFOREST:
    case BiomeType.TROPICAL_RAINFOREST:
      return 'forest';
    case BiomeType.HOT_DESERT:
      return 'dune';
    default:
      return null;
  }
};

/**
 * Prominence orders the greedy thinning pass: the most significant feature in
 * a crowded neighbourhood is the one that survives. Relief sorts by elevation;
 * vegetation sorts below all relief so mountains never lose to a tree.
 */
const prominence = (cell: Cell, kind: GlyphKind, seaLevel: number): number => {
  const f = landFrac(cell, seaLevel);
  if (kind === 'ice') return 2000 + f * 100; // never thinned away by relief
  if (kind === 'mountain' || kind === 'hill') return 1000 + f * 1000;
  return f * 100;
};

/** Deterministic 0-1 hash from a cell id and the world seed. */
const hash01 = (cellId: number, seed: string): number => {
  let h = 2166136261 >>> 0;
  const str = `${seed}:${cellId}`;
  for (let i = 0; i < str.length; i++) {
    h ^= str.charCodeAt(i);
    h = Math.imul(h, 16777619) >>> 0;
  }
  return (h >>> 8) / 0x1000000;
};

/**
 * Resolve every land cell that earns a glyph into output-pixel coordinates,
 * then thin greedily so no two glyphs collide.
 *
 * Substrate-independent BY DESIGN (spec §3.3): glyph size and collision are
 * screen-space decisions, and Mercator vs Mollweide distort area wildly at high
 * latitude. Both substrates draw the returned list, so this logic exists once.
 *
 * Spacing and size scale with `widthPx`, so an 8192px export is the same map as
 * a 1024px one at higher resolution — not a denser map.
 */
export const placeGlyphs = (
  cells: Cell[],
  projection: d3.GeoProjection,
  widthPx: number,
  opts: GlyphOptions,
): PlacedGlyph[] => {
  const seed = opts.seed ?? '';
  const k = widthPx / REFERENCE_WIDTH;
  const minSpacing = (opts.minSpacingPx ?? DEFAULT_SPACING_PX) * k;
  const size = (opts.sizePx ?? DEFAULT_SIZE_PX) * k;
  const heightPx = widthPx / 2;

  interface Candidate extends PlacedGlyph { rank: number }
  const candidates: Candidate[] = [];

  for (const cell of cells) {
    if (cell.height < opts.seaLevel) continue;
    if (cell.biome === BiomeType.LAKE || cell.biome === BiomeType.SALT_LAKE) continue;
    const kind = glyphFor(cell, opts.seaLevel);
    if (!kind) continue;

    const projected = projection(toLonLat([cell.center.x, cell.center.y, cell.center.z]));
    if (!projected) continue;
    const [x, y] = projected;
    if (!Number.isFinite(x) || !Number.isFinite(y)) continue;
    // Clip to the output box. Projections wrap or blow up past the antimeridian
    // and the poles; an off-canvas glyph is wasted work and can smear.
    if (x < 0 || y < 0 || x > widthPx || y > heightPx * 2) continue;

    const r = hash01(cell.id, seed);
    candidates.push({
      x, y, kind, cellId: cell.id,
      scale: size * (0.85 + r * 0.3),
      seedRot: (r - 0.5) * 0.35,
      rank: prominence(cell, kind, opts.seaLevel),
    });
  }

  // Most prominent first, id as a tiebreak so the order is total and stable.
  candidates.sort((a, b) => (b.rank - a.rank) || (a.cellId - b.cellId));

  // Greedy thinning against a uniform grid: only cells within one spacing can
  // collide, so this stays linear instead of quadratic on large worlds.
  const cellSize = Math.max(1, minSpacing);
  const grid = new Map<string, PlacedGlyph[]>();
  const key = (gx: number, gy: number) => `${gx},${gy}`;
  const accepted: PlacedGlyph[] = [];
  const minSq = minSpacing * minSpacing;

  for (const c of candidates) {
    const gx = Math.floor(c.x / cellSize);
    const gy = Math.floor(c.y / cellSize);
    let collides = false;
    for (let dx = -1; dx <= 1 && !collides; dx++) {
      for (let dy = -1; dy <= 1 && !collides; dy++) {
        for (const other of grid.get(key(gx + dx, gy + dy)) ?? []) {
          const ddx = other.x - c.x;
          const ddy = other.y - c.y;
          if (ddx * ddx + ddy * ddy < minSq) { collides = true; break; }
        }
      }
    }
    if (collides) continue;
    const glyph: PlacedGlyph = {
      x: c.x, y: c.y, kind: c.kind, scale: c.scale, seedRot: c.seedRot, cellId: c.cellId,
    };
    accepted.push(glyph);
    const bucket = grid.get(key(gx, gy));
    if (bucket) bucket.push(glyph); else grid.set(key(gx, gy), [glyph]);
  }

  return accepted;
};
