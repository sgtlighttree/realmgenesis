import { useEffect, useMemo, useRef } from 'react';
import { useThree, useFrame } from '@react-three/fiber';
import * as THREE from 'three';

import { WorldData } from '../../types';
import {
  ProjectedCells, StagedCells, stageCellPoints, projectStaged, projectLocalPoint,
} from '../../utils/screenProject';
import { displayRadius } from '../../utils/displayRadius';

// F2 screen-space overlay layer. Mounted INSIDE the R3F <Canvas>; it owns a
// sibling 2D <canvas> stacked over the WebGL canvas (pointer-events:none), and
// each frame projects visible cells to screen space and dispatches to tenant
// draw callbacks. This is the foundation that replaces physical 3D overlay
// objects (graticule/borders/…) — see spec §2.2/§3.2.

// Projects a LOCAL-frame point (globe's own coords) to screen pixels, applying
// the globe's live world matrix + camera. Writes px into `out`; returns false
// when the point is culled by the horizon. Tenants use it to project points
// other than cell centers (velocity tips, graticule samples).
//
// RADIUS CONTRACT (load-bearing — F2 parallax/occlusion fix): the horizon test
// uses the point's OWN |P| as the sphere radius, so callers MUST pass points at
// their intended RENDERED radius, not unit radius. The terrain mesh puts each
// cell at r = 1 + cell.height·0.05 (WorldViewer.tsx ~L898). A point passed at
// r=1 is culled inside the true limb and drifts off the terrain on zoom. Cell
// centers are pre-scaled here in ScreenOverlay; every other point a tenant
// projects must be scaled by the tenant (cell radius for cell-bound points, or a
// deliberate fixed radius e.g. sea level for the graticule). New tenants
// migrated onto this seam (roads/routes, contours, borders, rivers, labels) MUST
// honor this or they silently reintroduce the bug.
export type LocalProjector = (x: number, y: number, z: number, out: [number, number]) => boolean;

export interface OverlayTenant {
  id: string;
  visible: boolean;
  // `proj` cell coords are pre-scaled to each cell's rendered radius. Any point
  // a tenant projects itself via `project` must ALSO be at its rendered radius —
  // see the LocalProjector radius contract above. `smooth` is the smooth-globe
  // flag: when true the globe is a unit sphere (r=1), so tenants must project
  // their own points at r=1 too (use displayRadius).
  draw(ctx: CanvasRenderingContext2D, proj: ProjectedCells, world: WorldData, project: LocalProjector, smooth: boolean): void;
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

export const ScreenOverlay: React.FC<{ world: WorldData; tenants: OverlayTenant[]; smoothGlobe?: boolean }> = ({ world, tenants, smoothGlobe = false }) => {
  const { gl, camera, scene } = useThree();
  const canvasRef = useRef<HTMLCanvasElement | null>(null);
  const ctxRef = useRef<CanvasRenderingContext2D | null>(null);
  const projRef = useRef<ProjectedCells | null>(null);
  const stagedRef = useRef<StagedCells | null>(null);
  const staleRef = useRef<boolean>(true);
  const globeRef = useRef<THREE.Object3D | null>(null);
  const lastKey = useRef<string>('');

  const camPos = useMemo(() => new THREE.Vector3(), []);
  const camLocal = useMemo(() => new THREE.Vector3(), []);
  const mvp = useMemo(() => new THREE.Matrix4(), []);
  const gmInv = useMemo(() => new THREE.Matrix4(), []);

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

  // (Re)allocate the projection buffers when the cell count changes, and mark the
  // staged point cache stale. Any height mutation (paint, undo) ends in
  // setWorld({...world}) — a fresh `world` ref — so this fires before the next
  // frame and the cache can never silently decouple from the terrain radius (the
  // historical parallax failure mode). smoothGlobe toggles are handled below.
  const nCells = world.cells.length;
  useEffect(() => {
    projRef.current = {
      x: new Float32Array(nCells), y: new Float32Array(nCells),
      visible: new Uint8Array(nCells), n: nCells,
    };
    globeRef.current = null; // world changed → re-find the globe
    staleRef.current = true;
    lastKey.current = '';
  }, [nCells, world]);

  // Smooth toggle changes every cell's rendered radius (displayRadius) without a
  // new `world`, so it must invalidate the staged cache AND force a redraw.
  useEffect(() => {
    staleRef.current = true;
    lastKey.current = '';
  }, [smoothGlobe]);

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

    // Three recomputes world matrices inside render(), which runs AFTER every
    // useFrame callback. Reading matrixWorld here without forcing an update
    // therefore yields the PREVIOUS frame's transform while WebGL is about to
    // draw the current one — a permanent one-frame offset between the overlay
    // and the terrain, invisible when still and worse the further you zoom in,
    // because the error is angular.
    //
    // Camera.updateMatrixWorld() also recomputes matrixWorldInverse (three's
    // Camera overrides it to do exactly that), which is what `.project()` reads
    // — so this one call covers both, and the renderer recomputes the same
    // values from the same source moments later rather than fighting us.
    //
    // This is only half the fix: the globe's idle spin must also have advanced
    // before this callback runs, which is why it lives in <GlobeSpin/> mounted
    // ahead of this component (see WorldViewer).
    camera.updateMatrixWorld();
    if (globe) globe.updateWorldMatrix(true, false);

    const gm = globe ? globe.matrixWorld : null;

    // Redraw only when the camera, globe rotation, or active tenants changed.
    // Keyed on camera.matrixWorld (pose) + globe matrix; assumes a fixed
    // projection (dolly-zoom via OrbitControls). If FOV/zoom is ever animated,
    // add camera.projectionMatrix to the key so it doesn't go stale.
    const active = tenants.filter((t) => t.visible);
    const key = (smoothGlobe ? 's|' : 'r|') + active.map((t) => t.id).join(',') + '|' + matrixKey(camera.matrixWorld, gm ? gm.elements : IDENT);
    if (key === lastKey.current) return;
    lastKey.current = key;

    ctx.setTransform(dpr, 0, 0, dpr, 0, 0); // draw in CSS pixels
    ctx.clearRect(0, 0, cssW, cssH);
    if (active.length === 0) return;

    // Rebuild the staged rendered-radius points only when the world or smooth
    // flag changed (marked by the effects above). They depend on cell.height +
    // smooth, never the camera, so this is per-world, not per-frame (~0.5 ms).
    if (staleRef.current || !stagedRef.current || stagedRef.current.n !== nCells) {
      stagedRef.current = stageCellPoints(world.cells, displayRadius, smoothGlobe, stagedRef.current);
      staleRef.current = false;
    }
    const staged = stagedRef.current;

    // Fold the whole projection into one MVP matrix + the camera expressed in the
    // globe's LOCAL frame, both hoisted out of the per-cell loop. The naive form
    // transformed each cell local→world then world→NDC (two matrix ops per cell);
    // this is one. mvp = projection · matrixWorldInverse · globe. Reads the SAME
    // freshly-updated matrices as the renderer (camera.updateMatrixWorld and
    // globe.updateWorldMatrix above) — identical frame, only faster math.
    mvp.multiplyMatrices(camera.projectionMatrix, camera.matrixWorldInverse);
    if (gm) mvp.multiply(gm);
    camera.getWorldPosition(camPos);
    if (gm) {
      gmInv.copy(gm).invert();
      camLocal.copy(camPos).applyMatrix4(gmInv);
    } else {
      camLocal.copy(camPos);
    }

    projectStaged(staged, mvp.elements, camLocal.x, camLocal.y, camLocal.z, cssW, cssH, proj);

    // LocalProjector for tenant-projected points (velocity tips, graticule
    // samples): same fused mvp + camLocal, points already in the LOCAL frame.
    const project: LocalProjector = (x, y, z, out) =>
      projectLocalPoint(x, y, z, mvp.elements, camLocal.x, camLocal.y, camLocal.z, cssW, cssH, out);

    for (const t of active) t.draw(ctx, proj, world, project, smoothGlobe);
  });

  return null;
};
