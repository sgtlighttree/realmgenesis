/**
 * Camera-distance cutoffs for the globe's per-tier settlement markers, and the
 * camera's own range.
 *
 * Split out of `WorldViewer.tsx` so `overlays/tenants.ts` and the test suite can
 * read them without importing R3F. That matters here specifically: these
 * constants and the label cutoffs in `passesGlobeLod` must stay equal — a
 * settlement's dot and its name have to appear at the same zoom, or a dot
 * floats with no label.
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

export const VILLAGE_MARKER_MAX_DIST = 1.9;
export const HAMLET_MARKER_MAX_DIST = 1.5;
