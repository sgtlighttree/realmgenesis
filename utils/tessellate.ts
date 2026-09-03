import earcut from 'earcut';

export interface CellTessellation {
  positions: Float32Array;
  cellTriRange: Uint32Array;
}

// Earcut each cell's sub-rings independently and concatenate the resulting
// triangle-list positions. Cells are convex on the sphere but NOT after
// projection (Mercator polar stretch) or antimeridian clipping (a sub-ring
// bounded partly by the map edge), so a centroid fan would emit triangles
// outside the polygon exactly at cell corners — earcut avoids that. A cell
// whose clip produced no sub-rings (a pole under a cylindrical projection)
// contributes zero triangles.
export const tessellateCells = (
  cellVerts: Float32Array,
  cellOffsets: Uint32Array,
  cellSubStart: Uint32Array,
  cellSubOffsets: Uint32Array,
  cellCount: number,
): CellTessellation => {
  const out: number[] = [];
  const cellTriRange = new Uint32Array(cellCount + 1);
  let vertCount = 0;
  for (let i = 0; i < cellCount; i++) {
    cellTriRange[i] = vertCount;
    const subLo = cellSubStart[i];
    const subHi = cellSubStart[i + 1];
    for (let s = subLo; s < subHi; s++) {
      const ringStart = cellSubOffsets[s];             // vertex index
      const ringEnd = s + 1 < subHi ? cellSubOffsets[s + 1] : cellOffsets[i + 1];
      const count = ringEnd - ringStart;
      if (count < 3) continue;                          // degenerate sub-ring
      const coords: number[] = new Array(count * 2);
      for (let k = 0; k < count; k++) {
        coords[k * 2] = cellVerts[(ringStart + k) * 2];
        coords[k * 2 + 1] = cellVerts[(ringStart + k) * 2 + 1];
      }
      const tris = earcut(coords, undefined, 2);        // no holes; indices into `coords`
      for (const idx of tris) {
        out.push(coords[idx * 2], coords[idx * 2 + 1]);
        vertCount++;
      }
    }
  }
  cellTriRange[cellCount] = vertCount;
  return { positions: new Float32Array(out), cellTriRange };
};
