import { describe, it, expect } from 'vitest';
import { tessellateCells } from '../utils/tessellate';

// One square cell as a single sub-ring: (0,0)-(10,0)-(10,10)-(0,10).
function squareCell() {
  const verts = new Float32Array([0,0, 10,0, 10,10, 0,10]);
  const cellOffsets = new Uint32Array([0, 4]);
  const cellSubStart = new Uint32Array([0, 1]);
  const cellSubOffsets = new Uint32Array([0]);
  return { verts, cellOffsets, cellSubStart, cellSubOffsets };
}

describe('tessellateCells', () => {
  it('triangulates a convex quad into two triangles covering its area', () => {
    const { verts, cellOffsets, cellSubStart, cellSubOffsets } = squareCell();
    const t = tessellateCells(verts, cellOffsets, cellSubStart, cellSubOffsets, 1);
    const triVerts = t.cellTriRange[1] - t.cellTriRange[0];
    expect(triVerts).toBe(6); // 2 triangles
    // area of emitted triangles ~= 100
    let area = 0;
    for (let i = 0; i < t.positions.length; i += 6) {
      const [ax,ay,bx,by,cx,cy] = t.positions.slice(i, i+6);
      area += Math.abs((bx-ax)*(cy-ay) - (cx-ax)*(by-ay)) / 2;
    }
    expect(area).toBeCloseTo(100, 5);
  });

  it('triangulates a NON-convex (L-shaped) ring without triangles outside the polygon', () => {
    // L-shape: earcut must not fan across the concavity.
    const verts = new Float32Array([0,0, 10,0, 10,4, 4,4, 4,10, 0,10]);
    const cellOffsets = new Uint32Array([0, 6]);
    const cellSubStart = new Uint32Array([0, 1]);
    const cellSubOffsets = new Uint32Array([0]);
    const t = tessellateCells(verts, cellOffsets, cellSubStart, cellSubOffsets, 1);
    let area = 0;
    for (let i = 0; i < t.positions.length; i += 6) {
      const [ax,ay,bx,by,cx,cy] = t.positions.slice(i, i+6);
      area += Math.abs((bx-ax)*(cy-ay) - (cx-ax)*(by-ay)) / 2;
    }
    // L area = 100 - 36 (6x6 notch) = 64; a bad centroid fan would overshoot this.
    // NOTE: brief's test constant was 76 (100-24), which does not match the
    // shoelace area of the given hexagon; corrected to the true polygon area.
    expect(area).toBeCloseTo(64, 4);
  });

  it('emits an empty range for a cell with zero sub-rings (polar/degenerate)', () => {
    const verts = new Float32Array([]);
    const cellOffsets = new Uint32Array([0, 0]);
    const cellSubStart = new Uint32Array([0, 0]);
    const cellSubOffsets = new Uint32Array([]);
    const t = tessellateCells(verts, cellOffsets, cellSubStart, cellSubOffsets, 1);
    expect(t.cellTriRange[1] - t.cellTriRange[0]).toBe(0);
  });

  it('triangulates a cell with two disjoint sub-rings (antimeridian split) independently', () => {
    // two separate squares as two sub-rings of one cell
    const verts = new Float32Array([0,0, 2,0, 2,2, 0,2,  10,0, 12,0, 12,2, 10,2]);
    const cellOffsets = new Uint32Array([0, 8]);
    const cellSubStart = new Uint32Array([0, 2]);
    const cellSubOffsets = new Uint32Array([0, 4]);
    const t = tessellateCells(verts, cellOffsets, cellSubStart, cellSubOffsets, 1);
    expect(t.cellTriRange[1] - t.cellTriRange[0]).toBe(12); // 4 triangles total
  });
});
