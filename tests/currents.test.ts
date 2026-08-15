import { describe, it, expect } from 'vitest';
import { generateWorld } from '../utils/worldGen';
import { computeOceanCurrents, computeSstAnomaly } from '../utils/currents';
import { makeParams } from './helpers';
import type { Cell } from '../types';

// Build a lightweight world once for structural fixtures.
const buildCells = async (): Promise<Cell[]> =>
  (await generateWorld(makeParams({ seed: 'currents_test' }))).cells;

// A wind field that pushes every cell 'eastward' (tangent), for a deterministic probe.
const eastwardWind = (cells: Cell[]) =>
  cells.map(c => {
    const len = Math.hypot(c.center.x, c.center.z) || 1;
    return { x: -c.center.z / len, y: 0, z: c.center.x / len };
  });

describe('ocean currents', () => {
  it('is deterministic: same inputs -> identical velocity and anomaly', async () => {
    const cells = await buildCells();
    const wind = eastwardWind(cells);
    const p = makeParams({ seed: 'currents_test' });
    const a = computeOceanCurrents(cells, wind, p.seaLevel, 1.0);
    const b = computeOceanCurrents(cells, wind, p.seaLevel, 1.0);
    expect(Array.from(a.vx)).toEqual(Array.from(b.vx));
    expect(Array.from(a.vy)).toEqual(Array.from(b.vy));
    expect(Array.from(a.vz)).toEqual(Array.from(b.vz));
    const sa = computeSstAnomaly(cells, a, p, p.seaLevel);
    const sb = computeSstAnomaly(cells, b, p, p.seaLevel);
    expect(Array.from(sa)).toEqual(Array.from(sb));
  });

  it('leaves land cells with zero velocity', async () => {
    const cells = await buildCells();
    const p = makeParams({ seed: 'currents_test' });
    const field = computeOceanCurrents(cells, eastwardWind(cells), p.seaLevel, 1.0);
    cells.forEach((c, i) => {
      if (c.height >= p.seaLevel) {
        expect(field.vx[i]).toBe(0);
        expect(field.vy[i]).toBe(0);
        expect(field.vz[i]).toBe(0);
      }
    });
  });

  it('boundary tangency: coastal ocean flow has no NET into-land component', async () => {
    const cells = await buildCells();
    const p = makeParams({ seed: 'currents_test' });
    const field = computeOceanCurrents(cells, eastwardWind(cells), p.seaLevel, 1.0);
    let checked = 0;
    cells.forEach((c, i) => {
      if (c.height >= p.seaLevel) return;
      // Summed unit land-normal (tangent-projected), matching the solver invariant.
      let lx = 0, ly = 0, lz = 0, hasLand = false;
      for (const n of c.neighbors) {
        if (cells[n].height < p.seaLevel) continue; // only land neighbours
        let dx = cells[n].center.x - c.center.x;
        let dy = cells[n].center.y - c.center.y;
        let dz = cells[n].center.z - c.center.z;
        const dp = dx * c.center.x + dy * c.center.y + dz * c.center.z;
        dx -= dp * c.center.x; dy -= dp * c.center.y; dz -= dp * c.center.z;
        const dl = Math.hypot(dx, dy, dz);
        if (dl < 1e-9) continue;
        lx += dx / dl; ly += dy / dl; lz += dz / dl;
        hasLand = true;
      }
      if (!hasLand) return;
      const ll = Math.hypot(lx, ly, lz);
      if (ll < 1e-9) return;
      lx /= ll; ly /= ll; lz /= ll;
      const into = field.vx[i] * lx + field.vy[i] * ly + field.vz[i] * lz;
      const speed = Math.hypot(field.vx[i], field.vy[i], field.vz[i]);
      // Net into-land flow is removed as the last step of each pass -> ~0.
      if (speed > 1e-6) expect(into).toBeLessThan(speed * 0.05 + 1e-6);
      checked++;
    });
    expect(checked).toBeGreaterThan(0); // fixture actually has coastlines
  });

  it('poleward current warms high latitudes: SST anomaly correlates with poleward flow', async () => {
    const cells = await buildCells();
    const p = makeParams({ seed: 'currents_test' });
    const field = computeOceanCurrents(cells, eastwardWind(cells), p.seaLevel, 1.0);
    const anomaly = computeSstAnomaly(cells, field, p, p.seaLevel);
    // At a high-latitude ocean cell whose current flows poleward (|y| increasing),
    // anomaly should be >= 0; equatorward flow <= 0. Aggregate over the sphere.
    let polewardWarm = 0, polewardCold = 0;
    cells.forEach((c, i) => {
      if (c.height >= p.seaLevel || Math.abs(c.center.y) < 0.3) return;
      const poleward = Math.sign(c.center.y) * field.vy[i]; // >0 means flowing toward its pole
      if (poleward > 1e-4) { anomaly[i] >= 0 ? polewardWarm++ : polewardCold++; }
    });
    expect(polewardWarm).toBeGreaterThanOrEqual(polewardCold);
  });
});
