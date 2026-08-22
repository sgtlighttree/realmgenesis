import { describe, it, expect } from 'vitest';
import { readFileSync } from 'node:fs';
import { resolve } from 'node:path';
import { buildDymaxionNet } from '../utils/dymaxion';

// Ground-truth icosphere UV + vertex data, dumped live from Blender 5.1's
// DEFAULT icosphere (1 subdivision, default UV unwrap) via the MCP bridge on
// 2026-08-22. See docs/dymaxion.md for the extraction script. This guards
// buildBlenderNet('blender') against ever drifting from Blender's real unwrap —
// which is the whole point of the Blender layout (drop the exported PNG straight
// onto the icosphere as a UV texture, no tweaking).
const gt = JSON.parse(
  readFileSync(resolve(__dirname, 'fixtures/dymaxionBlenderUV.json'), 'utf8'),
);

describe('Blender Dymaxion net matches Blender ground truth', () => {
  const net = buildDymaxionNet('blender');

  it('has 20 faces', () => {
    expect(net.faces).toHaveLength(20);
  });

  it('reproduces every face UV to within 1e-3', () => {
    let maxUV = 0;
    for (let i = 0; i < 20; i++) {
      for (let k = 0; k < 3; k++) {
        const nu = net.faces[i].vertices[k];
        const gu = gt.faces[i].loops[k].uv;
        maxUV = Math.max(maxUV, Math.abs(nu[0] - gu[0]), Math.abs(nu[1] - gu[1]));
      }
    }
    expect(maxUV).toBeLessThan(1e-3);
  });

  it('reproduces every face 3D vertex (Blender z-up -> engine y-up) to within 1e-3', () => {
    let max3D = 0;
    for (let i = 0; i < 20; i++) {
      for (let k = 0; k < 3; k++) {
        const n3 = net.faces[i].vertices3D[k];
        const gc = gt.faces[i].loops[k].co; // z-up [x,y,z]
        const gy = [gc[0], gc[2], gc[1]];   // -> y-up [x,z,y]
        max3D = Math.max(max3D, Math.abs(n3[0] - gy[0]), Math.abs(n3[1] - gy[1]), Math.abs(n3[2] - gy[2]));
      }
    }
    expect(max3D).toBeLessThan(1e-3);
  });
});
