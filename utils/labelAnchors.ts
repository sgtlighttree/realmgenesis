import { Cell } from '../types';
import { MapLabel } from './labels';
import { nearestCellBrute, nearestCellWalk } from './nearestCell';

// F2 point-labels migration (S22). `MapLabel.position` is a UNIT direction
// (collectLabels runs normalizeVec), so the overlay tenant needs a deliberate
// per-label radius, exactly like the graticule and contour/border segments —
// see drawLabelsTenant in components/overlays/tenants.ts.
//
// The naive approach — call nearestCellBrute/nearestCellWalk per label inside
// the tenant, every redraw — is wrong for two reasons: it re-walks on every
// frame the overlay redraws, and labels are sorted by PRIORITY, so consecutive
// labels jump across the sphere and a hill-climb walk restarts cold each time
// (at 200k cells that is ~20M distance tests per redraw). Instead this computes
// one height per label ONCE, walking in a SPATIALLY sorted order so consecutive
// lookups are geometric neighbors and the walk actually gets to hill-climb.
export interface LabelAnchor {
  label: MapLabel;
  height: number;
}

export function computeLabelAnchors(labels: MapLabel[], cells: Cell[]): LabelAnchor[] {
  const n = labels.length;
  if (n === 0 || cells.length === 0) {
    return labels.map((label) => ({ label, height: 0 }));
  }

  // Index array sorted by direction (y then x is enough — see plan §2), NOT
  // by priority. This is the whole point: it decouples anchor lookup order
  // from label draw/declutter order so the walk stays warm.
  const order = labels.map((_, i) => i);
  order.sort((a, b) => {
    const pa = labels[a].position;
    const pb = labels[b].position;
    if (pa.y !== pb.y) return pa.y - pb.y;
    return pa.x - pb.x;
  });

  const heights = new Array<number>(n);
  let cellId = -1; // -1 = seed via brute scan on the first (sorted) label
  for (const i of order) {
    const p = labels[i].position;
    cellId = cellId < 0
      ? nearestCellBrute(cells, p.x, p.y, p.z)
      : nearestCellWalk(cells, p.x, p.y, p.z, cellId);
    heights[i] = cells[cellId].height;
  }

  // Scatter back into INPUT order — callers (the declutter pass, tests) index
  // this array positionally against the original `labels` array.
  return labels.map((label, i) => ({ label, height: heights[i] }));
}
