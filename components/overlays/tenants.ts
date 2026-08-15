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
    // Lift the tip to the cell's rendered radius so it sits on the terrain
    // surface exactly like the arrow base (see LocalProjector radius contract).
    const r = 1 + cells[i].height * 0.05;
    if (!project((c.x + cur.vx[i] * k) * r, (c.y + cur.vy[i] * k) * r, (c.z + cur.vz[i] * k) * r, tip)) continue;
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

const D2R = Math.PI / 180;
const GRAT_SEG = 96; // samples per line (smooth circles)

// 10° graticule drawn in screen space: parallels −80..80°, meridians 0..350°,
// each sampled on the sphere, projected via the horizon-culling projector.
// The polyline breaks wherever a sample crosses behind the limb, so lines never
// draw across the globe silhouette (the win over the old always-visible 3D grid).
//
// Radius: the grid is not cell-bound, so it needs a deliberate radius (see the
// LocalProjector contract). We put it at the SEA-LEVEL rendered radius
// (1 + seaLevel·0.05) — coastlines render right there, so the grid is
// parallax-free at the edge the eye actually locks onto, and land peaks occlude
// it (the Google-Earth read). Locking it above max relief (~1.055) instead would
// float a visible halo off the ocean limb wherever there's no mountain.
export function drawGraticuleTenant(
  ctx: CanvasRenderingContext2D,
  _proj: ProjectedCells,
  world: WorldData,
  project: LocalProjector,
): void {
  ctx.strokeStyle = 'rgba(255,255,255,0.28)';
  ctx.lineWidth = 1;
  const pt: [number, number] = [0, 0];
  const R = 1 + world.params.seaLevel * 0.05;

  // parallels (constant latitude)
  for (let lat = -80; lat <= 80; lat += 10) {
    const la = lat * D2R;
    const cy = Math.sin(la) * R;
    const cr = Math.cos(la) * R;
    let drawing = false;
    ctx.beginPath();
    for (let s = 0; s <= GRAT_SEG; s++) {
      const lon = (s / GRAT_SEG) * Math.PI * 2;
      if (project(cr * Math.cos(lon), cy, cr * Math.sin(lon), pt)) {
        if (drawing) ctx.lineTo(pt[0], pt[1]);
        else { ctx.moveTo(pt[0], pt[1]); drawing = true; }
      } else {
        drawing = false;
      }
    }
    ctx.stroke();
  }

  // meridians (constant longitude, pole to pole)
  for (let lon = 0; lon < 360; lon += 10) {
    const lo = lon * D2R;
    const cl = Math.cos(lo);
    const sl = Math.sin(lo);
    let drawing = false;
    ctx.beginPath();
    for (let s = 0; s <= GRAT_SEG; s++) {
      const lat = (s / GRAT_SEG) * Math.PI - Math.PI / 2;
      const cla = Math.cos(lat) * R;
      if (project(cla * cl, Math.sin(lat) * R, cla * sl, pt)) {
        if (drawing) ctx.lineTo(pt[0], pt[1]);
        else { ctx.moveTo(pt[0], pt[1]); drawing = true; }
      } else {
        drawing = false;
      }
    }
    ctx.stroke();
  }
}
