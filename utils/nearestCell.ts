import { Cell } from '../types';

export function nearestCellWalk(cells: Cell[], dx: number, dy: number, dz: number, startId: number): number {
  if (cells.length === 0) return -1;

  let current = Number.isFinite(startId)
    ? Math.max(0, Math.min(cells.length - 1, Math.trunc(startId)))
    : 0;

  for (let i = 0; i < cells.length; i++) {
    let best = current;
    let bestDot = dx * cells[current].center.x + dy * cells[current].center.y + dz * cells[current].center.z;

    for (const n of cells[current].neighbors) {
      if (n < 0 || n >= cells.length) continue;
      const c = cells[n].center;
      const d = dx * c.x + dy * c.y + dz * c.z;
      if (d > bestDot) {
        bestDot = d;
        best = n;
      }
    }

    if (best === current) break;
    current = best;
  }

  return current;
}

export function nearestCellBrute(cells: Cell[], dx: number, dy: number, dz: number): number {
  let best = -1;
  let bestDot = -Infinity;

  for (let i = 0; i < cells.length; i++) {
    const c = cells[i].center;
    const d = dx * c.x + dy * c.y + dz * c.z;
    if (d > bestDot) {
      bestDot = d;
      best = i;
    }
  }

  return best;
}
