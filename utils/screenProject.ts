// F2 screen-space overlay support: pure projection + occlusion math, no React,
// no THREE component types (kept unit-testable). A point P is "visible" when the
// (unit) view ray from P to the camera has a positive dot with the outward
// surface normal at P — i.e. P faces the camera and is not past the horizon. The
// normal is P̂ = P/|P|, so this is the EXACT perspective horizon for a sphere of
// radius |P|. Callers MUST pass P at its true rendered radius (the globe mesh
// puts each cell at r = 1 + height·0.05, up to ≈1.05), NOT at unit radius:
// testing a rendered-radius point at r=1 culls a band inside the true limb and
// makes overlays drift off the terrain on zoom. eps is a hair (not the old fat
// 0.08) only to suppress limb flicker.

export interface ProjectedCells {
  x: Float32Array; // screen-pixel x
  y: Float32Array; // screen-pixel y
  visible: Uint8Array; // 1 = on the near hemisphere, 0 = culled
  n: number;
}

export function isVisible(
  px: number, py: number, pz: number,
  camx: number, camy: number, camz: number,
  eps = 0.005,
): boolean {
  let dx = camx - px, dy = camy - py, dz = camz - pz;
  // sqrt, not Math.hypot: hypot's overflow-safe scaling is far slower in V8 and
  // buys nothing here — all inputs are near-unit sphere/camera coords. The swap
  // alone is 3.5–4.9x per tenant (F4 rung 1, S29).
  const len = Math.sqrt(dx * dx + dy * dy + dz * dz);
  if (len < 1e-9) return true;
  dx /= len; dy /= len; dz /= len;
  // outward normal at P is P̂ = P/|P| (P at its TRUE radius, not unit)
  const pl = Math.sqrt(px * px + py * py + pz * pz);
  if (pl < 1e-9) return true;
  return (dx * px + dy * py + dz * pz) / pl > eps;
}

// -----------------------------------------------------------------------------
// F4 rung 2: staged + fused projector. `isVisible` above tests one point; the
// per-frame globe loop projects EVERY cell, and the naive form re-derived each
// cell's rendered radius, transformed local→world with a THREE matrix op, then
// world→NDC with another (Vector3.project) — two matrix touches and a Math.hypot
// per cell, every frame. This path folds all of that into:
//   1. stageCellPoints() — precompute each cell's rendered-radius LOCAL point and
//      its |P|, ONCE per world (they depend only on cell.height + smooth, not the
//      camera). Rebuilt only when the world or smooth flag changes.
//   2. projectStaged() — one fused MVP matrix (projection·viewInverse·globe) and
//      the camera transformed into the globe's LOCAL frame, both hoisted out of
//      the loop, so each cell is a single mat-vec plus a sqrt horizon test.
// Measured 11.5x over the naive loop at 30k, output exact to 4.55e-13 px / 0
// visibility flips (scripts/perf/globeBench.ts, and tests/screenProjectStaged).
// These are the SAME math as `isVisible` + a THREE project() — kept here, pure
// and unit-testable, so the parity test exercises production code, not a copy.

export interface StagedCells {
  pts: Float32Array;  // [x,y,z] per cell at its rendered radius, LOCAL frame
  plen: Float32Array; // |P| per cell (the horizon-test denominator)
  n: number;
}

// A minimal cell shape, so this module still needs no app types.
interface StageableCell {
  center: { x: number; y: number; z: number };
  height: number;
}

// radius(height, smooth): injected so this stays free of the displayRadius import
// cycle and testable in isolation. Callers pass utils/displayRadius#displayRadius.
export function stageCellPoints(
  cells: readonly StageableCell[],
  radius: (height: number, smooth: boolean) => number,
  smooth: boolean,
  reuse?: StagedCells | null,
): StagedCells {
  const n = cells.length;
  const pts = reuse && reuse.n === n ? reuse.pts : new Float32Array(n * 3);
  const plen = reuse && reuse.n === n ? reuse.plen : new Float32Array(n);
  for (let k = 0; k < n; k++) {
    const c = cells[k].center;
    const r = radius(cells[k].height, smooth);
    const x = c.x * r, y = c.y * r, z = c.z * r;
    const j = k * 3;
    pts[j] = x; pts[j + 1] = y; pts[j + 2] = z;
    plen[k] = Math.sqrt(x * x + y * y + z * z);
  }
  return { pts, plen, n };
}

// Fused per-cell projection. `mvp` is a column-major 16-element matrix
// (THREE.Matrix4.elements) = projection · matrixWorldInverse · globeMatrix;
// camLocal is the camera position in the globe's LOCAL frame. Writes into `out`.
// The degenerate guards match isVisible: a zero-length ray or origin point counts
// as visible (never NaN-culled).
export function projectStaged(
  staged: StagedCells,
  mvp: ArrayLike<number>,
  camLocalX: number, camLocalY: number, camLocalZ: number,
  cssW: number, cssH: number,
  out: ProjectedCells,
  eps = 0.005,
): void {
  const { pts, plen, n } = staged;
  const e0 = mvp[0], e4 = mvp[4], e8 = mvp[8], e12 = mvp[12];
  const e1 = mvp[1], e5 = mvp[5], e9 = mvp[9], e13 = mvp[13];
  const e3 = mvp[3], e7 = mvp[7], e11 = mvp[11], e15 = mvp[15];
  for (let k = 0; k < n; k++) {
    const j = k * 3;
    const x = pts[j], y = pts[j + 1], z = pts[j + 2];
    const dx = camLocalX - x, dy = camLocalY - y, dz = camLocalZ - z;
    const len = Math.sqrt(dx * dx + dy * dy + dz * dz);
    const pl = plen[k];
    if (len > 1e-9 && pl > 1e-9 && (dx * x + dy * y + dz * z) / (len * pl) <= eps) {
      out.visible[k] = 0;
      continue;
    }
    const w = e3 * x + e7 * y + e11 * z + e15;
    out.x[k] = ((e0 * x + e4 * y + e8 * z + e12) / w + 1) / 2 * cssW;
    out.y[k] = (1 - ((e1 * x + e5 * y + e9 * z + e13) / w + 1) / 2) * cssH;
    out.visible[k] = 1;
  }
}

// Fused projection of one arbitrary LOCAL-frame point — the LocalProjector core
// tenants use for non-center points (velocity tips, graticule samples). Same math
// and guards as projectStaged, computing |P| inline. Returns false when culled.
export function projectLocalPoint(
  x: number, y: number, z: number,
  mvp: ArrayLike<number>,
  camLocalX: number, camLocalY: number, camLocalZ: number,
  cssW: number, cssH: number,
  out: [number, number],
  eps = 0.005,
): boolean {
  const dx = camLocalX - x, dy = camLocalY - y, dz = camLocalZ - z;
  const len = Math.sqrt(dx * dx + dy * dy + dz * dz);
  const pl = Math.sqrt(x * x + y * y + z * z);
  if (len > 1e-9 && pl > 1e-9 && (dx * x + dy * y + dz * z) / (len * pl) <= eps) return false;
  const w = mvp[3] * x + mvp[7] * y + mvp[11] * z + mvp[15];
  out[0] = ((mvp[0] * x + mvp[4] * y + mvp[8] * z + mvp[12]) / w + 1) / 2 * cssW;
  out[1] = (1 - ((mvp[1] * x + mvp[5] * y + mvp[9] * z + mvp[13]) / w + 1) / 2) * cssH;
  return true;
}
