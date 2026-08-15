// Shared globe display radius (F2). The 3D globe renders each cell at
// r = 1 + height·0.05 (per-cell relief). "Smooth globe" mode collapses every
// cell to a single unit sphere (r = 1) so screen-space overlays — the graticule
// and ocean-current arrows — can never parallax against terrain at any zoom:
// terrain and overlay then share one sphere by construction. `offset` lifts an
// overlay just above the surface (rivers, routes, contours, selection rings).
export function displayRadius(height: number, smooth: boolean, offset = 0): number {
  return (smooth ? 1 : 1 + height * 0.05) + offset;
}
