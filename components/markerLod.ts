/**
 * Camera-distance cutoffs for the globe's small settlement LABELS, and the
 * camera's own range.
 *
 * Originally these gated 3D marker pins as well, which is why they were split
 * out of `WorldViewer.tsx` — so `overlays/tenants.ts` and the tests could read
 * them without importing R3F. Only capitals draw a pin now (S40), so they
 * govern labels alone; the file keeps its name and location because
 * `passesGlobeLod` is still the consumer and the reachability invariant below
 * is the reason it exists.
 *
 * A cutoff is only meaningful INSIDE the camera's reachable range. The hamlet
 * cutoff was originally 1.15 against a `minDistance` of 1.2, so no camera
 * position could ever satisfy it and hamlets were dead on the globe by
 * construction. `tests/globeLodReachable.test.ts` pins that.
 */

/** OrbitControls `minDistance` — the closest the camera can get to the globe. */
export const GLOBE_MIN_CAM_DIST = 1.2;
/** OrbitControls `maxDistance` — fully zoomed out. */
export const GLOBE_MAX_CAM_DIST = 6;

export const VILLAGE_LABEL_MAX_DIST = 1.9;
export const HAMLET_LABEL_MAX_DIST = 1.5;
