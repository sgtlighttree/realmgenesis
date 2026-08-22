// F2 Dymaxion cage geometry. The cage is an icosahedron WIREFRAME floated as a
// reference frame over the 3D globe (Export tab -> Dymaxion -> Show Overlay). It
// used to be a physical 3D object (`DymaxionOverlay` in WorldViewer); this holds
// the geometry so the ScreenOverlay `dymaxion` tenant can draw it in Canvas2D.
//
// The 12 unit vertices and 30 edges are the EXACT output of
// `THREE.IcosahedronGeometry(1, 0)` + `THREE.EdgesGeometry` (dumped from three
// 5.1, deduped), so the migrated cage sits in the same base orientation as the
// old one — same `lon/lat/roll` settings, same cage on screen.

import * as THREE from 'three';

import { Point } from '../types';
import { DymaxionSettings } from '../types';

// 12 unit-sphere icosahedron vertices (THREE.IcosahedronGeometry orientation).
export const ICOSA_VERTS: Point[] = [
  { x: -0.8506507874, y: 0, z: 0.5257310867 },
  { x: 0, y: 0.5257310867, z: 0.8506507874 },
  { x: -0.5257310867, y: 0.8506507874, z: 0 },
  { x: 0.5257310867, y: 0.8506507874, z: 0 },
  { x: 0, y: 0.5257310867, z: -0.8506507874 },
  { x: -0.8506507874, y: 0, z: -0.5257310867 },
  { x: 0.8506507874, y: 0, z: 0.5257310867 },
  { x: 0, y: -0.5257310867, z: 0.8506507874 },
  { x: -0.5257310867, y: -0.8506507874, z: 0 },
  { x: 0, y: -0.5257310867, z: -0.8506507874 },
  { x: 0.8506507874, y: 0, z: -0.5257310867 },
  { x: 0.5257310867, y: -0.8506507874, z: 0 },
];

// 30 undirected edges as index pairs into ICOSA_VERTS.
export const ICOSA_EDGES: [number, number][] = [
  [1, 2], [2, 3], [2, 4], [0, 2], [2, 5], [1, 3], [0, 1], [0, 5], [4, 5], [3, 4],
  [7, 11], [8, 11], [9, 11], [6, 11], [10, 11], [1, 6], [1, 7], [6, 7], [0, 7],
  [0, 8], [7, 8], [5, 8], [5, 9], [8, 9], [4, 9], [4, 10], [9, 10], [3, 10],
  [3, 6], [6, 10],
];

// Rotate the cage by the settings euler and return its 30 edges as pairs of
// unit-vector endpoints (LOCAL frame — ScreenOverlay applies the globe spin).
// The euler matches the old `<Group rotation={new THREE.Euler(lat, -lon, roll,
// 'YXZ')}>`: same convention, so the cage tracks the sliders identically.
export function cageEdges(settings: DymaxionSettings): [Point, Point][] {
  const euler = new THREE.Euler(
    THREE.MathUtils.degToRad(settings.lat),
    -THREE.MathUtils.degToRad(settings.lon),
    THREE.MathUtils.degToRad(settings.roll),
    'YXZ',
  );
  const v = new THREE.Vector3();
  const rotated: Point[] = ICOSA_VERTS.map((p) => {
    v.set(p.x, p.y, p.z).applyEuler(euler);
    return { x: v.x, y: v.y, z: v.z };
  });
  return ICOSA_EDGES.map(([i, j]) => [rotated[i], rotated[j]] as [Point, Point]);
}
