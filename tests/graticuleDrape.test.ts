import { describe, it, expect } from 'vitest';
import { drawGraticuleTenant } from '../components/overlays/tenants';
import { ProjectedCells } from '../utils/screenProject';
import { makeFakeCtx, makeProjector, makeRingWorld } from './helpers/overlayCanvas';

const emptyProj: ProjectedCells = {
  x: new Float32Array(0), y: new Float32Array(0), visible: new Uint8Array(0), n: 0,
};

describe('graticule drape radius', () => {
  it('drapes ocean cells at their RAW seafloor radius, not at sea level', () => {
    // Every cell is deep ocean (height 0) under a seaLevel of 0.5. The terrain
    // mesh renders those at displayRadius(0, false) = 1.0. S18 clamped the grid
    // to max(height, seaLevel) = 0.5 → r = 1.025, floating it above every ocean
    // cell. That gap was the residual parallax; this asserts it is gone.
    const world = makeRingWorld(64, [0], 0.5);
    const p = makeProjector();
    drawGraticuleTenant(makeFakeCtx(), emptyProj, world, p.project, false);

    expect(p.radii.length).toBeGreaterThan(0);
    for (const r of p.radii) expect(r).toBeCloseTo(1.0, 6);
  });

  it('drapes land cells at 1 + height * 0.05', () => {
    const world = makeRingWorld(64, [1.0], 0.5);
    const p = makeProjector();
    drawGraticuleTenant(makeFakeCtx(), emptyProj, world, p.project, false);

    for (const r of p.radii) expect(r).toBeCloseTo(1.05, 6);
  });

  it('projects on the unit sphere when the globe is smooth', () => {
    const world = makeRingWorld(64, [1.0], 0.5);
    const p = makeProjector();
    drawGraticuleTenant(makeFakeCtx(), emptyProj, world, p.project, true);

    for (const r of p.radii) expect(r).toBeCloseTo(1.0, 6);
  });
});
