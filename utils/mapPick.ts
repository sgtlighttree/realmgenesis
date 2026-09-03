import { quadtree, Quadtree } from 'd3-quadtree';
import type { GeoProjection } from 'd3-geo';

import { WorldData } from '../types';
import { toLonLat } from './geo';

// Index each cell id by its projected CSS-px centre (un-flipped, matching the
// coordinate space getCellIdAtMapPoint queries after size.width - mapX). Replaces
// the O(cells) geodesic scan; nearest-projected-centre differs from geodesic only
// near the antimeridian/poles, which the parity test samples deliberately.
export const buildCellQuadtree = (
  world: WorldData,
  projection: GeoProjection,
  width: number,
  height: number,
): Quadtree<number> => {
  void width; void height; // reserved: projection already fits [width,height]
  const projected: [number, number][] = new Array(world.cells.length);
  for (let i = 0; i < world.cells.length; i++) {
    const c = world.cells[i].center;
    const p = projection(toLonLat([c.x, c.y, c.z]));
    projected[i] = p && Number.isFinite(p[0]) ? [p[0], p[1]] : [NaN, NaN];
  }
  const qt = quadtree<number>()
    .x((id) => projected[id][0])
    .y((id) => projected[id][1]);
  // Add only cells that projected to a finite point.
  for (let i = 0; i < world.cells.length; i++) {
    if (Number.isFinite(projected[i][0])) qt.add(i);
  }
  return qt;
};

export const findCellIdAtPoint = (
  qt: Quadtree<number>,
  xUnflipped: number,
  y: number,
): number | null => {
  const id = qt.find(xUnflipped, y);
  return id === undefined ? null : id;
};
