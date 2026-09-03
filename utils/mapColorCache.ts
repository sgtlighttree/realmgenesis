import { WorldData, ViewMode } from '../types';
import { ColorContext, getCellColor } from './colors';
import { seasonalTemperatureDelta } from './seasons';

// Per-cell fill colour, computed ONCE per (viewMode, season, colormaps, seaLevel,
// shadeMap). The fill passes read cache[i] instead of calling getCellColor +
// seasonalTemperatureDelta + allocating a THREE.Color every redraw — the other
// half of the 2.2s stall. Strings (#rrggbb) because Canvas2D fillStyle needs them.
export const buildCellColorCache = (
  world: WorldData,
  viewMode: ViewMode,
  colorCtx: ColorContext,
  shadeMap: Float32Array | null,
): string[] => {
  const out: string[] = new Array(world.cells.length);
  for (let i = 0; i < world.cells.length; i++) {
    const cell = world.cells[i];
    const color = getCellColor(cell, viewMode, {
      ...colorCtx,
      seasonalDelta: seasonalTemperatureDelta(cell, world.params),
    });
    if (shadeMap) color.multiplyScalar(shadeMap[cell.id]);
    out[i] = `#${color.getHexString()}`;
  }
  return out;
};
