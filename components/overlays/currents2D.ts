import { WorldData } from '../../types';
import { CURRENT_SPEED_MIN, CURRENT_ARC, CURRENT_SPEED_REF, CURRENT_STRIDE, currentTint } from './tenants';

// F2 currents on the 2D map. Mirrors Map2D's river/route drawing idiom: convert
// each ocean cell's velocity into a short lon/lat arrow and draw it through the
// caller's d3 projection. Standalone (no Map2D internals) so it drops into any
// projection pass. Same tint/speed constants as the globe tenant (single source).

const R2D = 180 / Math.PI;
const toLon = (x: number, z: number) => Math.atan2(z, x) * R2D;
const toLat = (y: number) => Math.asin(Math.max(-1, Math.min(1, y))) * R2D;

export function drawCurrents2D(
  ctx: CanvasRenderingContext2D,
  world: WorldData,
  projection: (lonlat: [number, number]) => [number, number] | null,
  dpr: number,
): void {
  const cur = world.currents;
  if (!cur) return;
  const cells = world.cells;
  const sea = world.params.seaLevel;
  ctx.lineWidth = Math.max(0.6, 1.1 / dpr);
  const headR = Math.max(1, 1.3 / dpr);

  for (let i = 0; i < cells.length; i++) {
    if (cells[i].height >= sea) continue; // ocean only
    if (CURRENT_STRIDE > 1 && i % CURRENT_STRIDE !== 0) continue;
    const vx = cur.vx[i], vy = cur.vy[i], vz = cur.vz[i];
    const sp = Math.hypot(vx, vy, vz);
    if (sp < CURRENT_SPEED_MIN) continue;

    const c = cells[i].center;
    const k = (CURRENT_ARC * Math.min(1, sp / CURRENT_SPEED_REF)) / sp;
    const baseLon = toLon(c.x, c.z), baseLat = toLat(c.y);
    const tipLon = toLon(c.x + vx * k, c.z + vz * k), tipLat = toLat(c.y + vy * k);
    if (Math.abs(tipLon - baseLon) > 180) continue; // antimeridian guard

    const a = projection([baseLon, baseLat]);
    const b = projection([tipLon, tipLat]);
    if (!a || !b) continue;

    const color = currentTint(cur.sst[i]);
    ctx.strokeStyle = color;
    ctx.beginPath();
    ctx.moveTo(a[0], a[1]);
    ctx.lineTo(b[0], b[1]);
    ctx.stroke();
    ctx.fillStyle = color;
    ctx.beginPath();
    ctx.arc(b[0], b[1], headR, 0, Math.PI * 2);
    ctx.fill();
  }
}
