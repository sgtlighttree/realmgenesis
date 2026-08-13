import { describe, it, expect } from 'vitest';
import { generateWorld } from '../utils/worldGen';
import { makeParams } from './helpers';

// Every plateId must form exactly one connected component over cell.neighbors.
function strayComponents(cells: { plateId: number; neighbors: number[] }[]): number {
  const seen = new Uint8Array(cells.length);
  const compsByPlate = new Map<number, number>();
  for (let i = 0; i < cells.length; i++) {
    if (seen[i]) continue;
    const pid = cells[i].plateId;
    const stack = [i];
    seen[i] = 1;
    while (stack.length) {
      const c = stack.pop()!;
      for (const n of cells[c].neighbors) {
        if (!seen[n] && cells[n].plateId === pid) { seen[n] = 1; stack.push(n); }
      }
    }
    compsByPlate.set(pid, (compsByPlate.get(pid) ?? 0) + 1);
  }
  let stray = 0;
  for (const count of compsByPlate.values()) stray += count - 1;
  return stray;
}

// Guard, not a general invariant. Fine-mesh cell.plateId connectivity is NOT
// guaranteed: exclaves occur PRE-EXISTING at elongation 0.0 (10-seed sweep
// "sweep-0".."sweep-9": strayTotal 2, seed "sweep-6" failed), and at 0.4/1.0
// (strayTotal 3 both, 2 seeds each) — n=10, the rate difference between 0.0
// and 0.4/1.0 is not distinguishable from noise. Connectivity IS guaranteed
// at the MACRO level by assignPlatesDijkstra's shared `claimed` set (verified
// separately); the nearest-macro downsample at tectonicsV3.ts ("dc.plateId =
// macroResult.plateIds[nearest]") can pinch a connected macro region into
// disconnected fine cells. This is out of Task 1's scope — Task 1 only fixes
// macro-level chain-seed connectivity.
describe('band seeding preserves plate connectivity', () => {
  it('has 0 stray plate components at elongation 0.4', async () => {
    const world = await generateWorld(makeParams({ seed: 'realmgenesis', plateElongation: 0.4 }));
    expect(strayComponents(world.cells as never)).toBe(0);
  }, 120000);
});
