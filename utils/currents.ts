import { Cell, WorldParams } from '../types';
import { annualMeanLatTemp } from './seasons';

// D2 ocean currents — a fixed-pass, fixed-order relaxation (no Poisson solve,
// no RNG) so byte-identical reruns hold, matching the 8-pass moisture loop's
// construction. See docs/superpowers/specs/2026-08-15-d2-ocean-currents-design.md.

export interface OceanCurrentField {
  vx: Float32Array;
  vy: Float32Array;
  vz: Float32Array;
}

// Tuned in the render-verify step (plan Task 5). Keep as named constants.
const DRAG = 0.6;        // wind-stress -> initial current speed
const CORIOLIS_K = 0.9;  // deflection magnitude scale
const MIX = 0.3;         // advective neighbour blend per pass
const N_VEL = 12;        // velocity relaxation passes
const MAX_SPEED = 1.0;   // speed ceiling
const N_HEAT = 12;       // heat-advection passes
const HEAT_MIX = 0.5;    // upstream blend per pass

// Coupling constants consumed by worldGen (centralized here for tuning).
export const COAST_K = 0.6; // land coastal-moderation weight
export const EVAP_K = 0.02; // evaporation boost per +°C of warm anomaly

const projectTangent = (
  vx: number, vy: number, vz: number,
  nx: number, ny: number, nz: number,
) => {
  const d = vx * nx + vy * ny + vz * nz;
  return { x: vx - d * nx, y: vy - d * ny, z: vz - d * nz };
};

export function computeOceanCurrents(
  cells: Cell[],
  windVectors: { x: number; y: number; z: number }[],
  seaLevel: number,
  currentStrength: number,
): OceanCurrentField {
  const n = cells.length;
  const vx = new Float32Array(n);
  const vy = new Float32Array(n);
  const vz = new Float32Array(n);
  const isOcean = (i: number) => cells[i].height < seaLevel;

  // Init from wind stress (tangent), ocean only.
  for (let i = 0; i < n; i++) {
    if (!isOcean(i)) continue;
    const c = cells[i].center;
    const t = projectTangent(windVectors[i].x, windVectors[i].y, windVectors[i].z, c.x, c.y, c.z);
    vx[i] = t.x * DRAG; vy[i] = t.y * DRAG; vz[i] = t.z * DRAG;
  }

  const sx = new Float32Array(n), sy = new Float32Array(n), sz = new Float32Array(n);
  for (let pass = 0; pass < N_VEL; pass++) {
    for (let i = 0; i < n; i++) {
      if (!isOcean(i)) { sx[i] = 0; sy[i] = 0; sz[i] = 0; continue; }
      const c = cells[i].center;
      // (a) Coriolis: rotate about local normal by theta ~ currentStrength*sin(lat).
      // sin(lat) === c.y on the unit sphere. Rodrigues about normal k=c: for a
      // tangent v, v_rot = v*cosT + (k x v)*sinT (the (k·v)k term vanishes since
      // v is tangent). At lat≈0, theta→0: equatorial flow is pure wind-driven
      // zonal flow with no deflection — physically correct, NOT a bug.
      const lat = c.y;
      const theta = CORIOLIS_K * currentStrength * lat;
      const cosT = Math.cos(theta), sinT = Math.sin(theta);
      let ax = vx[i], ay = vy[i], az = vz[i];
      const kx = c.x, ky = c.y, kz = c.z;
      const crx = ky * az - kz * ay;
      const cry = kz * ax - kx * az;
      const crz = kx * ay - ky * ax;
      ax = ax * cosT + crx * sinT;
      ay = ay * cosT + cry * sinT;
      az = az * cosT + crz * sinT;
      // (b) advective smoothing toward mean ocean-neighbour velocity.
      let mx = 0, my = 0, mz = 0, cnt = 0;
      for (const nb of cells[i].neighbors) {
        if (!isOcean(nb)) continue;
        mx += vx[nb]; my += vy[nb]; mz += vz[nb]; cnt++;
      }
      if (cnt > 0) {
        ax = ax * (1 - MIX) + (mx / cnt) * MIX;
        ay = ay * (1 - MIX) + (my / cnt) * MIX;
        az = az * (1 - MIX) + (mz / cnt) * MIX;
      }
      // (c) boundary tangency: remove the NET into-land component (component
      // along the summed unit land-normal). Per-edge tangency is geometrically
      // unsatisfiable when land surrounds a cell on multiple sides; the net
      // normal enforces "no net flow into the coast" and is satisfiable, so a
      // concave cell relaxes toward along-shore / near-zero flow.
      let lx = 0, ly = 0, lz = 0;
      for (const nb of cells[i].neighbors) {
        if (isOcean(nb)) continue;
        let dx = cells[nb].center.x - c.x, dy = cells[nb].center.y - c.y, dz = cells[nb].center.z - c.z;
        const dp = dx * c.x + dy * c.y + dz * c.z; // project edge dir to tangent
        dx -= dp * c.x; dy -= dp * c.y; dz -= dp * c.z;
        const dl = Math.hypot(dx, dy, dz);
        if (dl < 1e-9) continue;
        lx += dx / dl; ly += dy / dl; lz += dz / dl;
      }
      const ll = Math.hypot(lx, ly, lz);
      if (ll > 1e-9) {
        lx /= ll; ly /= ll; lz /= ll;
        const into = ax * lx + ay * ly + az * lz;
        if (into > 0) { ax -= into * lx; ay -= into * ly; az -= into * lz; }
      }
      // reproject to tangent + clamp speed
      const t = projectTangent(ax, ay, az, c.x, c.y, c.z);
      const sp = Math.hypot(t.x, t.y, t.z);
      const scale = sp > MAX_SPEED ? MAX_SPEED / sp : 1;
      sx[i] = t.x * scale; sy[i] = t.y * scale; sz[i] = t.z * scale;
    }
    vx.set(sx); vy.set(sy); vz.set(sz);
  }
  return { vx, vy, vz };
}

export function computeSstAnomaly(
  cells: Cell[],
  field: OceanCurrentField,
  params: WorldParams,
  seaLevel: number,
): Float32Array {
  const n = cells.length;
  const isOcean = (i: number) => cells[i].height < seaLevel;
  const base = new Float32Array(n);
  const T = new Float32Array(n);
  for (let i = 0; i < n; i++) {
    if (!isOcean(i)) continue;
    const phi = Math.asin(Math.max(-1, Math.min(1, cells[i].center.y)));
    base[i] = annualMeanLatTemp(phi, params);
    T[i] = base[i];
  }
  const scratch = new Float32Array(n);
  for (let pass = 0; pass < N_HEAT; pass++) {
    for (let i = 0; i < n; i++) {
      if (!isOcean(i)) { scratch[i] = 0; continue; }
      const c = cells[i].center;
      let incoming = 0, cnt = 0;
      for (const nb of cells[i].neighbors) {
        if (!isOcean(nb)) continue;
        // edgeDir from nb toward i; if the current at nb flows toward i, nb is upstream.
        const dx = c.x - cells[nb].center.x;
        const dy = c.y - cells[nb].center.y;
        const dz = c.z - cells[nb].center.z;
        const dot = dx * field.vx[nb] + dy * field.vy[nb] + dz * field.vz[nb];
        if (dot > 0) { incoming += T[nb]; cnt++; }
      }
      if (cnt === 0) scratch[i] = T[i] * 0.95 + base[i] * 0.05;
      else scratch[i] = T[i] * (1 - HEAT_MIX) + (incoming / cnt) * HEAT_MIX;
    }
    T.set(scratch);
  }
  const anomaly = new Float32Array(n);
  for (let i = 0; i < n; i++) if (isOcean(i)) anomaly[i] = T[i] - base[i];
  return anomaly;
}
