import { describe, it, expect } from 'vitest';

import { WorldData } from '../types';
import { buildGlobeDirs } from '../utils/mapStyle/bakeTexture';

// Unit vectors, so the direction maths is exercised with real geometry.
const pt = (lonDeg: number, latDeg: number) => {
  const lon = (lonDeg * Math.PI) / 180;
  const lat = (latDeg * Math.PI) / 180;
  return { x: Math.cos(lat) * Math.cos(lon), y: Math.sin(lat), z: Math.cos(lat) * Math.sin(lon) };
};

const worldWith = (cells: { center: ReturnType<typeof pt>; vertices: ReturnType<typeof pt>[] }[]) =>
  ({ cells } as unknown as WorldData);

// This file used to test `buildGlobeUVs`, which assigned an equirectangular UV
// per vertex and needed two guards to survive linear interpolation: an
// antimeridian seam wrap and a polar collapse. The globe now derives its
// texture coordinate in the FRAGMENT shader, so neither invariant exists any
// more and both cases are gone with the code they covered. What is left to
// test is the buffer this replaced them with.
describe('buildGlobeDirs', () => {
  it('emits three floats per vertex, three vertices per triangle', () => {
    const world = worldWith([
      { center: pt(0, 0), vertices: [pt(-1, -1), pt(1, -1), pt(0, 1)] },
    ]);
    // 3 vertices in the fan -> 3 triangles -> 9 vertices -> 27 floats.
    expect(buildGlobeDirs(world).length).toBe(27);
  });

  it('lays the fan out as (centre, v1, v2), matching the position buffer', () => {
    const c = pt(30, 20);
    const a = pt(29, 19);
    const b = pt(31, 19);
    const world = worldWith([{ center: c, vertices: [a, b] }]);
    const d = buildGlobeDirs(world);
    // First triangle: centre, vertices[0], vertices[1].
    expect(d[0]).toBeCloseTo(c.x, 6);
    expect(d[1]).toBeCloseTo(c.y, 6);
    expect(d[2]).toBeCloseTo(c.z, 6);
    expect(d[3]).toBeCloseTo(a.x, 6);
    expect(d[6]).toBeCloseTo(b.x, 6);
    // Second triangle wraps back to vertices[0].
    expect(d[9]).toBeCloseTo(c.x, 6);
    expect(d[12]).toBeCloseTo(b.x, 6);
    expect(d[15]).toBeCloseTo(a.x, 6);
  });

  // The shader reads longitude and latitude off this vector with atan2 and
  // asin, both of which need it on the unit sphere. It renormalizes after
  // interpolation, but a non-unit input would still bias asin at the vertices.
  it('normalizes every direction onto the unit sphere', () => {
    const scale = (p: ReturnType<typeof pt>, k: number) => ({ x: p.x * k, y: p.y * k, z: p.z * k });
    const world = worldWith([
      { center: scale(pt(0, 45), 1.4), vertices: [scale(pt(1, 44), 0.7), scale(pt(-1, 44), 2.1)] },
    ]);
    const d = buildGlobeDirs(world);
    for (let i = 0; i < d.length; i += 3) {
      expect(Math.hypot(d[i], d[i + 1], d[i + 2])).toBeCloseTo(1, 6);
    }
  });

  // Directions must come from the undisplaced sphere. The refill loop scales
  // each vertex by displayRadius(cell.height), so if height ever leaked in
  // here, two neighbours at different heights would sample different texture
  // content either side of their shared edge.
  it('ignores cell height entirely', () => {
    const cells = [{ center: pt(10, 10), vertices: [pt(9, 9), pt(11, 9), pt(10, 11)] }];
    const flat = buildGlobeDirs(worldWith(cells));
    const tall = buildGlobeDirs(
      { cells: cells.map(c => ({ ...c, height: 0.9 })) } as unknown as WorldData,
    );
    expect(Array.from(tall)).toEqual(Array.from(flat));
  });
});
