import { describe, it, expect } from 'vitest';

import { WorldData } from '../types';
import { buildGlobeUVs } from '../utils/mapStyle/bakeTexture';

// Unit vectors, so the UV maths is exercised with real geometry.
const pt = (lonDeg: number, latDeg: number) => {
  const lon = (lonDeg * Math.PI) / 180;
  const lat = (latDeg * Math.PI) / 180;
  return { x: Math.cos(lat) * Math.cos(lon), y: Math.sin(lat), z: Math.cos(lat) * Math.sin(lon) };
};

const worldWith = (cells: { center: ReturnType<typeof pt>; vertices: ReturnType<typeof pt>[] }[]) =>
  ({ cells } as unknown as WorldData);

describe('buildGlobeUVs', () => {
  it('emits two floats per vertex, three vertices per triangle', () => {
    const world = worldWith([
      { center: pt(0, 0), vertices: [pt(-1, -1), pt(1, -1), pt(0, 1)] },
    ]);
    // 3 vertices in the fan -> 3 triangles -> 9 vertices -> 18 floats.
    expect(buildGlobeUVs(world).length).toBe(18);
  });

  it('maps longitude to u across the full range', () => {
    const world = worldWith([
      { center: pt(0, 0), vertices: [pt(0, 0), pt(0, 1), pt(1, 0)] },
    ]);
    const uv = buildGlobeUVs(world);
    // lon 0 sits at the centre of the texture.
    expect(uv[0]).toBeCloseTo(0.5, 3);
  });

  it('maps latitude to v with north at the top', () => {
    const north = worldWith([{ center: pt(0, 80), vertices: [pt(0, 80), pt(1, 80), pt(0, 81)] }]);
    const south = worldWith([{ center: pt(0, -80), vertices: [pt(0, -80), pt(1, -80), pt(0, -81)] }]);
    // v grows downward, so the northern cell must sit nearer 0.
    expect(buildGlobeUVs(north)[1]).toBeLessThan(buildGlobeUVs(south)[1]);
  });

  // THE seam bug this exists to prevent: a cell straddling lon +/-180 has
  // vertices at u ~ 0.99 and u ~ 0.01. Interpolated as-is, that triangle samples
  // the ENTIRE texture backwards and draws a bright smear down the seam.
  it('keeps antimeridian triangles narrow instead of wrapping the whole texture', () => {
    const world = worldWith([
      { center: pt(179.5, 0), vertices: [pt(179, 0), pt(-179, 0), pt(180, 1)] },
    ]);
    const uv = buildGlobeUVs(world);
    for (let t = 0; t < uv.length; t += 6) {
      const us = [uv[t], uv[t + 2], uv[t + 4]];
      // Every triangle stays local. Without the fix, spans approach 1.0.
      expect(Math.max(...us) - Math.min(...us)).toBeLessThan(0.5);
    }
  });

  it('pushes wrapped coordinates past 1.0 rather than clamping them', () => {
    // Values > 1 are correct ONLY because the texture uses RepeatWrapping.
    // Clamping instead would pile vertices onto the texture edge.
    const world = worldWith([
      { center: pt(179.5, 0), vertices: [pt(179, 0), pt(-179, 0), pt(180, 1)] },
    ]);
    const uv = buildGlobeUVs(world);
    const us = [];
    for (let i = 0; i < uv.length; i += 2) us.push(uv[i]);
    expect(us.some(u => u > 1)).toBe(true);
  });

  it('leaves ordinary triangles untouched', () => {
    const world = worldWith([
      { center: pt(10, 10), vertices: [pt(9, 9), pt(11, 9), pt(10, 11)] },
    ]);
    const uv = buildGlobeUVs(world);
    for (let i = 0; i < uv.length; i += 2) {
      expect(uv[i]).toBeGreaterThanOrEqual(0);
      expect(uv[i]).toBeLessThanOrEqual(1);
    }
  });
});
