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
export function placeSettlements(world: WorldData, params: WorldParams): void {
  if (!world.civData) return;
  const rng = new RNG(params.civSeed + '_settle');
  const getNameGen = makeStyleNameGenCache(params.civSeed + '_settlenames_');
  const settlements: Settlement[] = [];
  const K = 0.02; // tuned in Task-4 perf pass; ~1500 total at 200k
  const avgCellArea = (4 * Math.PI) / world.cells.length; // steradians per cell
  let nextId = 0;

  for (const faction of world.civData.factions) {
    const land = world.cells.filter(c =>
      c.regionId === faction.id && c.height >= params.seaLevel && !isLakeCell(c));
    if (land.length === 0) continue;

    const quota = Math.max(1, Math.round(land.length * params.settlementDensity * K));
    // Quota-derived spacing (Azgaar mechanism): disc-pack `quota` points into the
    // faction's area. chord² threshold, tier-scaled per candidate below.
    const factionArea = land.length * avgCellArea;
    const baseSpacingRad = Math.sqrt(factionArea / quota);
    const baseChordSq = (2 * Math.sin(baseSpacingRad / 2)) ** 2;

    // Rank by score; RNG side-stream breaks ties deterministically.
    const tiebreak = new Map<number, number>();
    for (const c of land) tiebreak.set(c.id, rng.next());
    const ranked = [...land].sort((a, b) =>
      settlementScore(b, world, params) - settlementScore(a, world, params)
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

    // Tier by rank position among accepted (capital forced to 'capital').
    accepted.forEach((c, i) => {
      const tier: SettlementTier = c.id === capCell.id ? 'capital'
        : i < accepted.length * 0.12 ? 'city'
        : i < accepted.length * 0.35 ? 'town'
        : i < accepted.length * 0.70 ? 'village' : 'hamlet';
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
