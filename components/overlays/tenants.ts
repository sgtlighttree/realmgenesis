// F2 screen-space overlay tenants — pure Canvas2D draw functions dispatched by
// ScreenOverlay. Each receives the projected cell screen coords, the world, and
// a projector for arbitrary local-frame points.

import { WorldData } from '../../types';
import { ProjectedCells } from '../../utils/screenProject';
import { LocalProjector } from './ScreenOverlay';

// --- currents draw constants (single source; also consumed by Map2D, Task 6) ---
export const CURRENT_SPEED_MIN = 0.04; // skip near-still cells
export const CURRENT_ARC = 0.05;       // max arrow arc length (unit-sphere)
export const CURRENT_SPEED_REF = 0.6;  // speed at which arrows reach full length
export const CURRENT_STRIDE = 1;       // draw every Nth qualifying ocean cell
export const SST_CLAMP = 6;            // °C anomaly mapped to full cold/warm
const COLD: [number, number, number] = [0x3b, 0x6f, 0xb0];
const WARM: [number, number, number] = [0xc0, 0x44, 0x2e];

// Warm/cold tint from SST anomaly (°C), clamped to ±SST_CLAMP.
export function currentTint(sst: number): string {
  const t = Math.max(0, Math.min(1, (sst + SST_CLAMP) / (2 * SST_CLAMP)));
  const r = Math.round(COLD[0] + (WARM[0] - COLD[0]) * t);
  const g = Math.round(COLD[1] + (WARM[1] - COLD[1]) * t);
  const b = Math.round(COLD[2] + (WARM[2] - COLD[2]) * t);
  return `rgb(${r},${g},${b})`;
}

export function drawCurrentsTenant(
  ctx: CanvasRenderingContext2D,
  proj: ProjectedCells,
  world: WorldData,
  project: LocalProjector,
): void {
  const cur = world.currents;
  if (!cur) return;
  const cells = world.cells;
  const sea = world.params.seaLevel;
  const tip: [number, number] = [0, 0];
  ctx.lineWidth = 1.2;
  for (let i = 0; i < proj.n; i++) {
    if (!proj.visible[i]) continue;
    if (cells[i].height >= sea) continue; // ocean cells only
    if (CURRENT_STRIDE > 1 && i % CURRENT_STRIDE !== 0) continue;
    const sp = Math.hypot(cur.vx[i], cur.vy[i], cur.vz[i]);
    if (sp < CURRENT_SPEED_MIN) continue;
    // arrow tip = center + unit(vel) * (arc scaled by speed, capped)
    const mag = CURRENT_ARC * Math.min(1, sp / CURRENT_SPEED_REF);
    const k = mag / sp;
    const c = cells[i].center;
    if (!project(c.x + cur.vx[i] * k, c.y + cur.vy[i] * k, c.z + cur.vz[i] * k, tip)) continue;
    const color = currentTint(cur.sst[i]);
    ctx.strokeStyle = color;
    ctx.beginPath();
    ctx.moveTo(proj.x[i], proj.y[i]);
    ctx.lineTo(tip[0], tip[1]);
    ctx.stroke();
    ctx.fillStyle = color;
    ctx.beginPath();
    ctx.arc(tip[0], tip[1], 1.4, 0, Math.PI * 2);
    ctx.fill();
  }
}

export function drawGraticuleTenant(
  _ctx: CanvasRenderingContext2D,
  _proj: ProjectedCells,
  _world: WorldData,
  _project: LocalProjector,
): void {
  // implemented in Task 5
}
