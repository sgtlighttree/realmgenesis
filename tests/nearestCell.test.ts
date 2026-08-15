import { describe, it, expect } from 'vitest';
import { nearestCellWalk, nearestCellBrute } from "../utils/nearestCell";
import type { Cell } from '../types';

const PHI = (1 + Math.sqrt(5)) / 2;

function makeIcosahedron(): { center: { x: number; y: number; z: number }; neighbors: number[] }[] {
  const raw: number[][] = [];
  for (const s of [-1, 1]) {
    raw.push([0, s, s * PHI], [0, -s, s * PHI]);
  }
  for (const s of [-1, 1]) {
    raw.push([s, s * PHI, 0], [-s, s * PHI, 0]);
  }
  for (const s of [-1, 1]) {
    raw.push([s * PHI, 0, s], [s * PHI, 0, -s]);
  }

  const verts = raw.map(([x, y, z]) => {
    const len = Math.hypot(x, y, z);
    return { x: x / len, y: y / len, z: z / len };
  });

  const dots: number[][] = verts.map((a, i) => verts.map((b, j) => a.x * b.x + a.y * b.y + a.z * b.z));

  const neighbors: number[][] = verts.map((_, i) =>
    verts
      .map((_, j) => j)
      .filter((j) => j !== i)
      .sort((a, b) => dots[i][b] - dots[i][a])
      .slice(0, 5)
  );

  return verts.map((center, i) => ({ center, neighbors: neighbors[i] }));
}

function randomUnit(): [number, number, number] {
  const theta = Math.random() * 2 * Math.PI;
  const z = Math.random() * 2 - 1;
  const r = Math.sqrt(1 - z * z);
  return [r * Math.cos(theta), r * Math.sin(theta), z];
}

describe('nearestCell', () => {
  const graph = makeIcosahedron();
  const cells = graph as unknown as Cell[];

  it('walk returns the true argmax-dot from several startIds including the far side', () => {
    const [dx, dy, dz] = [0.5, 0.3, 0.8];

    const trueBest = nearestCellBrute(cells, dx, dy, dz);
    for (let start = 0; start < cells.length; start++) {
      expect(nearestCellWalk(cells, dx, dy, dz, start)).toBe(trueBest);
    }
  });

  it('walk === brute for several random directions', () => {
    for (let t = 0; t < 200; t++) {
      const [dx, dy, dz] = randomUnit();
      const brute = nearestCellBrute(cells, dx, dy, dz);
      for (const start of [0, 3, 5, 7, 11]) {
        expect(nearestCellWalk(cells, dx, dy, dz, start)).toBe(brute);
      }
    }
  });

  it('clamps out-of-range startIds and terminates', () => {
    const [dx, dy, dz] = [-0.9, 0.2, -0.3];
    const brute = nearestCellBrute(cells, dx, dy, dz);
    expect(nearestCellWalk(cells, dx, dy, dz, -100)).toBe(brute);
    expect(nearestCellWalk(cells, dx, dy, dz, 999)).toBe(brute);
    expect(nearestCellWalk(cells, dx, dy, dz, NaN)).toBe(brute);
  });

  it('handles an empty cell array', () => {
    expect(nearestCellWalk([], 1, 0, 0, 0)).toBe(-1);
    expect(nearestCellBrute([], 1, 0, 0)).toBe(-1);
  });
});
