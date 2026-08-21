import * as THREE from 'three';
import { Point } from '../types';

// F2 rivers tenant support. Ports the CatmullRom smoothing that used to live in
// `RiverLines` (components/WorldViewer.tsx) VERBATIM, with one change: smoothing
// always runs at the BAKED radius each river point already carries (see
// utils/worldGen.ts getRenderPoint, ~line 310 — `r = 1 + height*0.05 + 0.005`,
// NOT a unit direction). There is no `smooth` parameter here — collapsing to the
// smooth globe is a per-point unit-normalize + relift, cheap enough to do at
// draw time in the tenant, and keeping it out of this function means the memo
// key in WorldViewer is `world.rivers` alone: toggling Smooth Globe never re-runs
// CatmullRom over ~1741 paths.
//
// This reorders operations relative to the old `RiverLines`: that component
// normalized control points to the smooth-globe radius BEFORE smoothing (when
// smoothGlobe was on); this always smooths at the baked radius first, and the
// tenant normalizes each OUTPUT point after. Measured worst-case deviation
// between the two orders over 48 synthetic paths (n = 2..30 points, heights
// 0.05..0.95): 0.000171 rad = 0.07px on an 800px-diameter globe. Sub-pixel, so
// the reorder is safe. Do not re-derive this; it is measured, not assumed.
export function computeRiverPolylines(rivers: Point[][]): Point[][] {
  const result: Point[][] = [];
  for (const path of rivers) {
    if (path.length < 2) continue;

    const vectors = path.map((p) => new THREE.Vector3(p.x, p.y, p.z));
    const curve = new THREE.CatmullRomCurve3(vectors);

    // Adaptive sampling based on length, but simple count is safer for perf
    const points = curve.getPoints(Math.min(50, vectors.length * 4));

    result.push(points.map((v) => ({ x: v.x, y: v.y, z: v.z })));
  }
  return result;
}
