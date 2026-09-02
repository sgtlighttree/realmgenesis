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
 *
 * **Why 1.10 and not 1.03.** Measured, Boardgame at 20k, one seed, camera at the
 * equator, over the ocean pixels of the globe disc — Boardgame because its flat
 * pale sea has no hatch or grain to hide a seam, so it is the worst case:
 *
 *     overhang   core (r<0.80R)   limb (0.80-0.93R)   darkest ocean px
 *     1.03            0.283%            0.316%              17  (background)
 *     1.06            0.072%            0.062%              17  (background)
 *     1.10            0.000%            0.000%             168  (no gap left)
 *
 * At 1.03 the darkest ocean pixel is 17 — the dark background showing straight
 * through an open seam. At 1.10 nothing in the ocean is even 40% below the sea
 * tone, in the core AND at the limb, so the seams are closed rather than merely
 * thinned. The flat texture bake is clean at every value, which is what makes
 * this a mesh defect and not a style-pass one.
 *
 * **The counter-cost is smaller than expected.** Raising it was assumed to fatten
 * cliff overlap. It does not visibly: land at 1.10 is CLEANER, because the dark
 * dashes over land read as terrain detail but were the same seams. Coastline ink,
 * biome patches and relief glyphs are unchanged.
 *
 * A fixed constant is still the wrong shape long term — the seam is a radial step
 * proportional to relief, the closure is lateral and proportional to cell radius,
 * and cell radius falls as sqrt(1/N). So this weakens as cell count rises and a
 * count-derived overhang is the principled fix. 1.10 clears it at the counts in
 * use today.
 *
 * `showCellEdges` deliberately drops the overhang to 1 so the seams reopen and
 * the tessellation shows; that path is independent of this value.
 */
export const CELL_OVERHANG = 1.10;
