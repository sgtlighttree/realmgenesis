// Shared globe display radius (F2). The 3D globe renders each cell at
// r = 1 + height·0.05 (per-cell relief). "Smooth globe" mode collapses every
// cell to a single unit sphere (r = 1) so screen-space overlays — the graticule
// and ocean-current arrows — can never parallax against terrain at any zoom:
// terrain and overlay then share one sphere by construction. `offset` lifts an
// overlay just above the surface (rivers, routes, contours, selection rings).
export function displayRadius(height: number, smooth: boolean, offset = 0): number {
  return (smooth ? 1 : 1 + height * 0.05) + offset;
}

/**
 * How far each cell's rim is pushed out from its own centre on the 3D globe.
 *
 * Cell "outlines" there are open SEAMS, not lines: neighbours sit at different
 * radii, so they do not share an edge and the background shows through the gap
 * as dark dashes. Widening every cell about its own centre makes neighbours
 * overlap, so the taller one overhangs the shorter — which is what a cliff looks
 * like. Only the rim moves; the centre and the radius are untouched, so the
 * drape invariant (overlay radius == cell radius) is unaffected.
 *
 * Lives here rather than in `WorldViewer` because anything that rebuilds the
 * globe mesh needs it. `scripts/preview/globe.html` omitted it once and the
 * resulting seams were read as a rendering defect that did not exist.
 */
export const CELL_OVERHANG = 1.03;
