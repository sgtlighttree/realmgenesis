import { BiomeType, Cell, WorldData, WorldParams } from '../types';
import { isLakeCell } from './civNaming';

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
