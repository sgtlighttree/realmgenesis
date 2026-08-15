// F2 screen-space overlay tenants — pure Canvas2D draw functions dispatched by
// ScreenOverlay. Each receives the projected cell screen coords + the world.
// Filled in by later tasks: currents (Task 4), graticule (Task 5).

import { WorldData } from '../../types';
import { ProjectedCells } from '../../utils/screenProject';

export function drawCurrentsTenant(
  _ctx: CanvasRenderingContext2D,
  _proj: ProjectedCells,
  _world: WorldData,
): void {
  // implemented in Task 4
}

export function drawGraticuleTenant(
  _ctx: CanvasRenderingContext2D,
  _proj: ProjectedCells,
  _world: WorldData,
): void {
  // implemented in Task 5
}
