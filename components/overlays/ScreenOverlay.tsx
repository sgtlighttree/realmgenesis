import { useEffect, useMemo, useRef } from 'react';
import { useThree, useFrame } from '@react-three/fiber';
import * as THREE from 'three';

import { WorldData } from '../../types';
import { isVisible, ProjectedCells } from '../../utils/screenProject';

// F2 screen-space overlay layer. Mounted INSIDE the R3F <Canvas>; it owns a
// sibling 2D <canvas> stacked over the WebGL canvas (pointer-events:none), and
// each frame projects visible cells to screen space and dispatches to tenant
// draw callbacks. This is the foundation that replaces physical 3D overlay
// objects (graticule/borders/…) — see spec §2.2/§3.2.

// Projects a LOCAL-frame point (globe's own coords) to screen pixels, applying
// the globe's live world matrix + camera. Writes px into `out`; returns false
// when the point is culled by the horizon. Tenants use it to project points
// other than cell centers (velocity tips, graticule samples).
export type LocalProjector = (x: number, y: number, z: number, out: [number, number]) => boolean;

export interface OverlayTenant {
  id: string;
  visible: boolean;
  draw(ctx: CanvasRenderingContext2D, proj: ProjectedCells, world: WorldData, project: LocalProjector): void;
}

// The globe mesh is auto-rotated via a parent group; we read its live world
// matrix (found by name) so projection tracks the spin, not just the camera.
const GLOBE_NAME = 'globe-mesh';

const IDENT = new THREE.Matrix4().elements;

function matrixKey(a: THREE.Matrix4, b: Float32Array | number[]): string {
  const ea = a.elements;
  let s = '';
  for (let i = 0; i < 16; i++) s += ((ea[i] * 1000) | 0) + ',' + ((b[i] * 1000) | 0) + ';';
  return s;
}

export const ScreenOverlay: React.FC<{ world: WorldData; tenants: OverlayTenant[] }> = ({ world, tenants }) => {
  const { gl, camera, scene } = useThree();
  const canvasRef = useRef<HTMLCanvasElement | null>(null);
  const ctxRef = useRef<CanvasRenderingContext2D | null>(null);
  const projRef = useRef<ProjectedCells | null>(null);
  const globeRef = useRef<THREE.Object3D | null>(null);
  const lastKey = useRef<string>('');

  const scratch = useMemo(() => new THREE.Vector3(), []);
  const projScratch = useMemo(() => new THREE.Vector3(), []);
  const camPos = useMemo(() => new THREE.Vector3(), []);

  // Create the sibling 2D canvas over the WebGL canvas.
  useEffect(() => {
    const parent = gl.domElement.parentElement;
    if (!parent) return;
    const c = document.createElement('canvas');
    Object.assign(c.style, {
      position: 'absolute', inset: '0', width: '100%', height: '100%',
      pointerEvents: 'none', zIndex: '1',
    } as CSSStyleDeclaration);
    parent.appendChild(c);
    canvasRef.current = c;
    ctxRef.current = c.getContext('2d');
    lastKey.current = '';
    return () => {
      if (c.parentElement) c.parentElement.removeChild(c);
      canvasRef.current = null;
      ctxRef.current = null;
    };
  }, [gl]);

  // (Re)allocate the projection buffers when the cell count changes.
  const nCells = world.cells.length;
  useEffect(() => {
    projRef.current = {
      x: new Float32Array(nCells), y: new Float32Array(nCells),
      visible: new Uint8Array(nCells), n: nCells,
    };
    globeRef.current = null; // world changed → re-find the globe
    lastKey.current = '';
  }, [nCells, world]);

  useFrame(() => {
    const ctx = ctxRef.current;
    const canvas = canvasRef.current;
    const proj = projRef.current;
    if (!ctx || !canvas || !proj) return;

    const cssW = gl.domElement.clientWidth;
    const cssH = gl.domElement.clientHeight;
    if (cssW === 0 || cssH === 0) return;
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    const bw = Math.round(cssW * dpr);
    const bh = Math.round(cssH * dpr);
    if (canvas.width !== bw || canvas.height !== bh) {
      canvas.width = bw; canvas.height = bh; lastKey.current = '';
    }

    if (!globeRef.current) globeRef.current = scene.getObjectByName(GLOBE_NAME) ?? null;
    const globe = globeRef.current;
    const gm = globe ? globe.matrixWorld : null;

    // Redraw only when the camera, globe rotation, or active tenants changed.
    // Keyed on camera.matrixWorld (pose) + globe matrix; assumes a fixed
    // projection (dolly-zoom via OrbitControls). If FOV/zoom is ever animated,
    // add camera.projectionMatrix to the key so it doesn't go stale.
    const active = tenants.filter((t) => t.visible);
    const key = active.map((t) => t.id).join(',') + '|' + matrixKey(camera.matrixWorld, gm ? gm.elements : IDENT);
    if (key === lastKey.current) return;
    lastKey.current = key;

    ctx.setTransform(dpr, 0, 0, dpr, 0, 0); // draw in CSS pixels
    ctx.clearRect(0, 0, cssW, cssH);
    if (active.length === 0) return;

    camera.getWorldPosition(camPos);
    for (let i = 0; i < nCells; i++) {
      const c = world.cells[i].center;
      scratch.set(c.x, c.y, c.z);
      if (gm) scratch.applyMatrix4(gm); // local → world (tracks spin)
      if (!isVisible(scratch.x, scratch.y, scratch.z, camPos.x, camPos.y, camPos.z)) {
        proj.visible[i] = 0;
        continue;
      }
      scratch.project(camera); // world → NDC
      proj.x[i] = (scratch.x + 1) / 2 * cssW;
      proj.y[i] = (1 - (scratch.y + 1) / 2) * cssH;
      proj.visible[i] = 1;
    }

    const project: LocalProjector = (x, y, z, out) => {
      projScratch.set(x, y, z);
      if (gm) projScratch.applyMatrix4(gm);
      if (!isVisible(projScratch.x, projScratch.y, projScratch.z, camPos.x, camPos.y, camPos.z)) return false;
      projScratch.project(camera);
      out[0] = (projScratch.x + 1) / 2 * cssW;
      out[1] = (1 - (projScratch.y + 1) / 2) * cssH;
      return true;
    };

    for (const t of active) t.draw(ctx, proj, world, project);
  });

  return null;
};
