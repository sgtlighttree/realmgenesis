import { describe, it, expect } from 'vitest';
import { buildFillColorBuffer, writeCellColors } from '../utils/mapFillColorBuffer';

describe('mapFillColorBuffer', () => {
  // cell 0 -> 3 verts (1 tri), cell 1 -> 6 verts (2 tris)
  const cellTriRange = new Uint32Array([0, 3, 9]);

  it('fills each cell vertex range with the cell colour as normalized RGB', () => {
    const hex = ['#ff0000', '#0080ff'];
    const buf = buildFillColorBuffer(hex, cellTriRange);
    expect(buf.length).toBe(9 * 3);
    // cell 0 vertices -> red
    expect(buf[0]).toBeCloseTo(1); expect(buf[1]).toBeCloseTo(0); expect(buf[2]).toBeCloseTo(0);
    // cell 1 first vertex (index 3) -> (0, 128/255, 1)
    expect(buf[3 * 3]).toBeCloseTo(0);
    expect(buf[3 * 3 + 1]).toBeCloseTo(128 / 255, 3);
    expect(buf[3 * 3 + 2]).toBeCloseTo(1);
  });

  it('writeCellColors rewrites only the named cell, leaving others untouched', () => {
    const hex = ['#ff0000', '#0080ff'];
    const buf = buildFillColorBuffer(hex, cellTriRange);
    const hex2 = ['#00ff00', '#0080ff'];
    writeCellColors(buf, hex2, cellTriRange, [0]);
    expect(buf[0]).toBeCloseTo(0); expect(buf[1]).toBeCloseTo(1); expect(buf[2]).toBeCloseTo(0);
    // cell 1 unchanged
    expect(buf[3 * 3 + 2]).toBeCloseTo(1);
  });
});
