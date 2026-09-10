import { BiomeType, Cell, Settlement, SettlementTier, WorldData, WorldParams } from '../types';
import { RNG } from './rng';
import { getCultureNameStyle, isLakeCell, makeStyleNameGenCache } from './civNaming';

export const TIER_ORDER: SettlementTier[] = ['capital', 'city', 'town', 'village', 'hamlet'];
export const TIER_MULTIPLIER: Record<SettlementTier, number> = {
  capital: 8,
  city: 5,
  town: 3,
  village: 1.5,
  hamlet: 1,
};

// Marker geometry per tier — the SINGLE source for both the 3D globe
// (WorldViewer's CityMarkers) and the GLB export, so the two cannot drift.
// Colour is a hex number: R3F's `color` prop and THREE materials both take it.
export const TIER_MARKER_STYLE: Record<SettlementTier, {
  radius: number; height: number; segments: number; color: number;
}> = {
  capital: { radius: 0.008, height: 0.08, segments: 6, color: 0xef4444 },
  city: { radius: 0.0065, height: 0.06, segments: 6, color: 0xfb923c },
  town: { radius: 0.005, height: 0.04, segments: 5, color: 0xffffff },
  village: { radius: 0.004, height: 0.025, segments: 4, color: 0xcbd5e1 },
  hamlet: { radius: 0.003, height: 0.018, segments: 4, color: 0x94a3b8 },
};

/**
 * Group settlement CELLS by tier. `cell.isTown` is set on EVERY settlement
 * including hamlets, so bucketing on it lumps a hamlet in with a city — every
 * marker/export consumer must group through this instead.
 */
export function settlementCellsByTier(world: WorldData): Record<SettlementTier, Cell[]> {
  const g: Record<SettlementTier, Cell[]> = { capital: [], city: [], town: [], village: [], hamlet: [] };
  (world.civData?.settlements ?? []).forEach(s => {
    const c = world.cells[s.cellId];
    if (c) g[s.tier].push(c);
  });
  return g;
}

export function settlementScore(cell: Cell, world: WorldData, params: WorldParams): number {
  if (cell.height < params.seaLevel || isLakeCell(cell)) return 0;
  let s = 0;
  switch (cell.biome) {
    case BiomeType.TROPICAL_RAINFOREST: s = 0.4; break;
    case BiomeType.TROPICAL_SAVANNA: s = 0.7; break;
    case BiomeType.HOT_DESERT: s = 0.1; break;
    case BiomeType.COLD_DESERT: s = 0.1; break;
    case BiomeType.TEMPERATE_FOREST: s = 0.9; break;
    case BiomeType.TEMPERATE_RAINFOREST: s = 0.8; break;
    case BiomeType.MEDITERRANEAN: s = 1.0; break;
    case BiomeType.STEPPE: s = 0.5; break;
    case BiomeType.GRASSLAND: s = 0.9; break;
    case BiomeType.BOREAL_FOREST: s = 0.4; break;
    case BiomeType.TUNDRA: s = 0.2; break;
    case BiomeType.ICE_CAP: s = 0.0; break;
    case BiomeType.VOLCANIC: s = 0.1; break;
    case BiomeType.BEACH: s = 0.6; break;
    default: s = 0.5;
  }
  if ((cell.flux || 0) > 0.5) s += 0.3;
  if ((cell.flux || 0) > 2.0) s += 0.2;
  let coast = false;
  for (const nId of cell.neighbors) {
    const n = world.cells[nId];
    if (n.height < params.seaLevel || isLakeCell(n)) { coast = true; break; }
  }
  if (coast) s += 0.3;
  if (cell.height > 0.6) s -= (cell.height - 0.6) * 2;
  return Math.max(0, s);
}

// Per-faction quota-then-score settlement placement with tiers. Not yet wired
// into the pipeline (recalculateCivs) — a later task does that. Populates
// world.civData.settlements, sets cell.isTown on settlement cells, and
// assigns tier/isPort/name/factionId. Does NOT set population — that's a
// later task, done in recalculateProvinces.
/**
 * The resolution `K = 0.02` was tuned at, and the constant itself. Kept only to
 * derive `SETTLEMENTS_PER_STERADIAN` below, so the calibration is visible rather
 * than a magic 318.
 */
const REFERENCE_POINTS = 200_000;
const REFERENCE_K = 0.02;

/**
 * Settlements per steradian of faction land.
 *
 * The quota used to be `landCells * density * K`, which made settlement count a
 * function of the MESH RESOLUTION rather than of the planet: the same world had
 * ~30 settlements at 5,000 cells and ~600 at 200,000 (measured in-app at 5k:
 * 1,553 land cells, 6 factions, 30 settlements total). At the resolution Matt
 * iterates on, the planet was empty and Settlement Density had almost nothing
 * to scale — the slider moved a total of roughly 15 to 60 across its whole
 * range, which is why it read as inert.
 *
 * Area is a property of the planet; cell count is a property of the
 * tessellation. This constant is derived so the new formula is IDENTICAL to the
 * old one at `REFERENCE_POINTS` and only diverges below it, which is precisely
 * the intent — high-resolution behaviour is unchanged.
 */
export const SETTLEMENTS_PER_STERADIAN = (REFERENCE_POINTS * REFERENCE_K) / (4 * Math.PI);

/**
 * Ceiling on settlements as a fraction of a faction's land CELLS, per unit of
 * density. A low-resolution world cannot host an area-derived quota — there are
 * not enough cells — so something has to bind, and it must not be a flat cap:
 * a density-independent ceiling would make the slider go dead the moment it
 * binds, exactly the trap the settlementDensity-0 finding documented. Scaling
 * the cap with density keeps the control live at every resolution.
 *
 * 0.20 is Matt's pick (S39): ~310 settlements on a 5k world against 30 before.
 */
export const MAX_SETTLEMENT_CELL_FRACTION = 0.20;

/**
 * How many settlements a faction gets. Pure, so it can be tested across
 * resolutions without generating a world.
 *
 * Floors at 1 so `density = 0` remains the deliberate capitals-only mode
 * (S36 finding 4), and never exceeds the land cells available.
 */
export function settlementQuota(
  landCellCount: number,
  totalCells: number,
  density: number,
): number {
  const avgCellArea = (4 * Math.PI) / totalCells; // steradians per cell
  const areaTarget = landCellCount * avgCellArea * SETTLEMENTS_PER_STERADIAN * density;
  // The cap FLOORS rather than rounds: rounding a cap can land above it, which
  // would make the ceiling a suggestion instead of a bound.
  const cellCap = Math.floor(landCellCount * MAX_SETTLEMENT_CELL_FRACTION * density);
  return Math.max(1, Math.min(Math.round(areaTarget), cellCap, landCellCount));
}

export function placeSettlements(world: WorldData, params: WorldParams): void {
  if (!world.civData) return;
  const rng = new RNG(params.civSeed + '_settle');
  const getNameGen = makeStyleNameGenCache(params.civSeed + '_settlenames_');
  const settlements: Settlement[] = [];
  const avgCellArea = (4 * Math.PI) / world.cells.length; // steradians per cell
  let nextId = 0;

  for (const faction of world.civData.factions) {
    const land = world.cells.filter(c =>
      c.regionId === faction.id && c.height >= params.seaLevel && !isLakeCell(c));
    if (land.length === 0) continue;

    // `?? 1.0` guards the browser-save load path (Controls.handleLoadBrowser ->
    // onLoadWorld) which passes stored params straight through without
    // withParamDefaults. An undefined density makes quota NaN, and NaN defeats
    // BOTH loop guards below (`>= NaN` and `< NaN` are false), so every ranked
    // land cell would be accepted as a settlement.
    const density = params.settlementDensity ?? 1.0;
    const quota = settlementQuota(land.length, world.cells.length, density);
    // Quota-derived spacing (Azgaar mechanism): disc-pack `quota` points into the
    // faction's area. chord² threshold, tier-scaled per candidate below.
    const factionArea = land.length * avgCellArea;
    const baseSpacingRad = Math.sqrt(factionArea / quota);
    const baseChordSq = (2 * Math.sin(baseSpacingRad / 2)) ** 2;

    // Rank by score; RNG side-stream breaks ties deterministically. Scores are
    // precomputed rather than called from the comparator — settlementScore
    // re-walks the cell's neighbour list, and a comparator evaluates it
    // O(n log n) times instead of O(n). Same values, so the order is identical;
    // this mirrors the `tiebreak` map immediately above.
    const tiebreak = new Map<number, number>();
    const score = new Map<number, number>();
    for (const c of land) {
      tiebreak.set(c.id, rng.next());
      score.set(c.id, settlementScore(c, world, params));
    }
    const ranked = [...land].sort((a, b) =>
      score.get(b.id)! - score.get(a.id)!
      || tiebreak.get(a.id)! - tiebreak.get(b.id)!);

    // Force the capital in first, then greedily accept with spacing.
    const accepted: Cell[] = [];
    const capCell = world.cells[faction.capitalId];
    accepted.push(capCell);
    const tierSpacingMult = (rankFrac: number) => rankFrac < 0.3 ? 1.4 : rankFrac < 0.6 ? 1.0 : 0.6;
    for (const c of ranked) {
      if (accepted.length >= quota) break;
      if (c.id === capCell.id) continue;
      const rankFrac = accepted.length / quota;
      const minChordSq = baseChordSq * tierSpacingMult(rankFrac);
      let ok = true;
      for (const a of accepted) {
        const d = (c.center.x - a.center.x) ** 2 + (c.center.y - a.center.y) ** 2 + (c.center.z - a.center.z) ** 2;
        if (d < minChordSq) { ok = false; break; }
      }
      if (ok) accepted.push(c);
    }

    // Tier by rank position among accepted. The capital is forced in at index 0
    // and always tiers 'capital', so rank fractions are taken over the NON-capital
    // candidates only — otherwise the first real candidate sits at i=1 and clears
    // the city cutoff only when accepted.length > 8.3, killing the whole city tier
    // on small worlds (measured: zero cities at 5k points).
    const nonCapCount = accepted.length - 1;
    accepted.forEach((c, i) => {
      const rank = i - 1; // -1 for the capital itself; unused on that branch
      const tier: SettlementTier = c.id === capCell.id ? 'capital'
        : rank < nonCapCount * 0.12 ? 'city'
        : rank < nonCapCount * 0.35 ? 'town'
        : rank < nonCapCount * 0.70 ? 'village' : 'hamlet';
      let isPort = false;
      for (const nId of c.neighbors) {
        const n = world.cells[nId];
        if (n.height < params.seaLevel && !isLakeCell(n)) { isPort = true; break; }
      }
      c.isTown = true;
      const nameGen = getNameGen(getCultureNameStyle(world, params, c.id));
      settlements.push({
        id: nextId++, cellId: c.id, name: nameGen.town(), population: 0,
        tier, factionId: faction.id, isCapital: tier === 'capital', isPort,
      });
    });
  }
  world.civData.settlements = settlements;
}
