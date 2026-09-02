import { describe, it, expect } from 'vitest';
import * as THREE from 'three';

import {
  ProjectedCells, isVisible, stageCellPoints, projectStaged, projectLocalPoint,
} from '../utils/screenProject';
import { displayRadius } from '../utils/displayRadius';

// F4 rung 2 parity guard. The staged+fused projector (stageCellPoints +
// projectStaged / projectLocalPoint) must produce byte-comparable output to the
// naive path it replaced in ScreenOverlay: transform each LOCAL point → world
// with the globe matrix, horizon-test with isVisible against the world-space
// camera, then world → NDC with THREE's Vector3.project. This test IS the
// discriminating check called out at the S29 commitment boundary — the earlier
// tenant tests never run the projection loop, so they cannot see a fused-matrix
// or staging regression. It exercises the production functions, not a copy.

interface FakeCell {
  center: { x: number; y: number; z: number };
  height: number;
}

// A deterministic scatter of unit-ish cell centers with varied heights.
function makeCells(n: number): FakeCell[] {
  const cells: FakeCell[] = [];
  // Fibonacci sphere so points cover both hemispheres (front + back of globe).
  const golden = Math.PI * (3 - Math.sqrt(5));
  for (let i = 0; i < n; i++) {
    const y = 1 - (i / (n - 1)) * 2;
    const rad = Math.sqrt(Math.max(0, 1 - y * y));
    const theta = golden * i;
    cells.push({
      center: { x: Math.cos(theta) * rad, y, z: Math.sin(theta) * rad },
      height: (Math.sin(i * 1.37) * 0.5 + 0.5), // 0..1
    });
  }
  return cells;
}

// The naive reference — exactly what ScreenOverlay did before rung 2.
function naiveProject(
  cells: FakeCell[], smooth: boolean, gm: THREE.Matrix4 | null,
  camera: THREE.Camera, cssW: number, cssH: number,
): ProjectedCells {
  const out: ProjectedCells = {
    x: new Float32Array(cells.length), y: new Float32Array(cells.length),
    visible: new Uint8Array(cells.length), n: cells.length,
  };
  const camPos = new THREE.Vector3();
  camera.getWorldPosition(camPos);
  const v = new THREE.Vector3();
  for (let i = 0; i < cells.length; i++) {
    const c = cells[i].center;
    const r = displayRadius(cells[i].height, smooth);
    v.set(c.x * r, c.y * r, c.z * r);
    if (gm) v.applyMatrix4(gm);
    if (!isVisible(v.x, v.y, v.z, camPos.x, camPos.y, camPos.z)) { out.visible[i] = 0; continue; }
    v.project(camera);
    out.x[i] = (v.x + 1) / 2 * cssW;
    out.y[i] = (1 - (v.y + 1) / 2) * cssH;
    out.visible[i] = 1;
  }
  return out;
}

// The production fused path, mirroring ScreenOverlay's useFrame math.
function fusedProject(
  cells: FakeCell[], smooth: boolean, gm: THREE.Matrix4 | null,
  camera: THREE.Camera, cssW: number, cssH: number,
): ProjectedCells {
  const out: ProjectedCells = {
    x: new Float32Array(cells.length), y: new Float32Array(cells.length),
    visible: new Uint8Array(cells.length), n: cells.length,
  };
  const staged = stageCellPoints(cells, displayRadius, smooth);
  const mvp = new THREE.Matrix4().multiplyMatrices(camera.projectionMatrix, camera.matrixWorldInverse);
  if (gm) mvp.multiply(gm);
  const camPos = new THREE.Vector3();
  camera.getWorldPosition(camPos);
  const camLocal = camPos.clone();
  if (gm) camLocal.applyMatrix4(new THREE.Matrix4().copy(gm).invert());
  projectStaged(staged, mvp.elements, camLocal.x, camLocal.y, camLocal.z, cssW, cssH, out);
  return out;
}

function makeCamera(cssW: number, cssH: number, camWorldPos: THREE.Vector3, rotY = 0): THREE.PerspectiveCamera {
  const cam = new THREE.PerspectiveCamera(50, cssW / cssH, 0.1, 100);
  cam.position.copy(camWorldPos);
  cam.rotateY(rotY);
  cam.lookAt(0, 0, 0);
  cam.updateMatrixWorld(true); // also refreshes matrixWorldInverse (THREE.Camera)
  return cam;
}

function assertParity(a: ProjectedCells, b: ProjectedCells): { maxErr: number; flips: number; tested: number } {
  let maxErr = 0, flips = 0, tested = 0;
  for (let i = 0; i < a.n; i++) {
    if (a.visible[i] !== b.visible[i]) { flips++; continue; }
    if (!a.visible[i]) continue;
    tested++;
    maxErr = Math.max(maxErr, Math.abs(a.x[i] - b.x[i]), Math.abs(a.y[i] - b.y[i]));
  }
  return { maxErr, flips, tested };
}

describe('staged+fused projector parity with the naive path', () => {
  const cssW = 1280, cssH = 720;
  const cells = makeCells(2000);

  const cases: Array<[string, boolean, THREE.Matrix4 | null, THREE.Vector3]> = [
    ['relief, no globe transform', false, null, new THREE.Vector3(0, 0, 2.5)],
    ['smooth, no globe transform', true, null, new THREE.Vector3(0, 0, 2.5)],
    ['relief, spun globe', false, new THREE.Matrix4().makeRotationY(1.1).multiply(new THREE.Matrix4().makeRotationX(0.4)), new THREE.Vector3(1.2, 0.8, 2.0)],
    ['relief, zoomed in close', false, new THREE.Matrix4().makeRotationY(2.3), new THREE.Vector3(0, 0, 1.4)],
    ['smooth, spun + off-axis cam', true, new THREE.Matrix4().makeRotationY(-0.7).multiply(new THREE.Matrix4().makeRotationZ(0.9)), new THREE.Vector3(-1.5, 1.0, 1.8)],
  ];

  for (const [name, smooth, gm, camPos] of cases) {
    it(`matches for: ${name}`, () => {
      const camera = makeCamera(cssW, cssH, camPos);
      const naive = naiveProject(cells, smooth, gm, camera, cssW, cssH);
      const fused = fusedProject(cells, smooth, gm, camera, cssW, cssH);
      const { maxErr, flips, tested } = assertParity(naive, fused);
      expect(tested).toBeGreaterThan(200); // a real chunk of the near hemisphere
      expect(flips).toBe(0);
      expect(maxErr).toBeLessThan(1e-3); // Float32 staging vs Float64 naive: sub-pixel
    });
  }

  it('projectLocalPoint matches projectStaged for cell centers', () => {
    const camera = makeCamera(cssW, cssH, new THREE.Vector3(0.9, 0.5, 2.1));
    const gm = new THREE.Matrix4().makeRotationY(0.6);
    const staged = stageCellPoints(cells, displayRadius, false);
    const mvp = new THREE.Matrix4().multiplyMatrices(camera.projectionMatrix, camera.matrixWorldInverse).multiply(gm);
    const camLocal = new THREE.Vector3();
    camera.getWorldPosition(camLocal).applyMatrix4(new THREE.Matrix4().copy(gm).invert());

    const out: ProjectedCells = {
      x: new Float32Array(cells.length), y: new Float32Array(cells.length),
      visible: new Uint8Array(cells.length), n: cells.length,
    };
    projectStaged(staged, mvp.elements, camLocal.x, camLocal.y, camLocal.z, cssW, cssH, out);

    const pt: [number, number] = [0, 0];
    let checked = 0;
    for (let i = 0; i < cells.length; i++) {
      const j = i * 3;
      const vis = projectLocalPoint(
        staged.pts[j], staged.pts[j + 1], staged.pts[j + 2],
        mvp.elements, camLocal.x, camLocal.y, camLocal.z, cssW, cssH, pt,
      );
      expect(vis ? 1 : 0).toBe(out.visible[i]);
      if (!vis) continue;
      checked++;
      expect(Math.abs(pt[0] - out.x[i])).toBeLessThan(1e-4);
      expect(Math.abs(pt[1] - out.y[i])).toBeLessThan(1e-4);
    }
    expect(checked).toBeGreaterThan(200);
  });

  it('stageCellPoints reuses buffers when the cell count is unchanged', () => {
    const s1 = stageCellPoints(cells, displayRadius, false);
    const s2 = stageCellPoints(cells, displayRadius, false, s1);
    expect(s2.pts).toBe(s1.pts);
    expect(s2.plen).toBe(s1.plen);
    const s3 = stageCellPoints(cells.slice(0, 100), displayRadius, false, s1);
    expect(s3.pts).not.toBe(s1.pts); // count changed → fresh allocation
  });
});
