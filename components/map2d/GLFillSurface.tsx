import React, { useEffect, useRef, useState } from 'react';
import * as THREE from 'three';

export interface GLFillHandle {
  setTransform: (scale: number, offsetX: number, offsetY: number) => void; // used in Task 6
  setColors: (colors: Float32Array, range?: { start: number; count: number }) => void; // Task 7
  redraw: () => void;
  isContextLost: () => boolean;
}

export interface GLFillSurfaceProps {
  positions: Float32Array;
  colors: Float32Array;
  width: number;
  height: number;
  dpr: number;
  /**
   * Fired only after GLFillSurface has exhausted its internal fresh-canvas
   * retries — i.e. WebGL is genuinely unavailable and the host should latch its
   * blit fallback. A single construction failure does NOT fire this; it is
   * recovered internally by remounting on a new canvas element.
   */
  onContextLost?: () => void;
  handleRef?: React.MutableRefObject<GLFillHandle | null>;
}

type FillMesh = THREE.Mesh<THREE.BufferGeometry, THREE.MeshBasicMaterial>;

// Expand tessellator output (2-component [x,y] pairs, CSS px) into a
// 3-component [x,y,0] position buffer for BufferGeometry — the ortho camera
// lives in the z=0 plane, so this is a one-time upload-time expansion rather
// than a custom shader.
const expandPositions = (flat: Float32Array): Float32Array => {
  const count = flat.length / 2;
  const out = new Float32Array(count * 3);
  for (let i = 0; i < count; i++) {
    out[i * 3] = flat[i * 2];
    out[i * 3 + 1] = flat[i * 2 + 1];
    out[i * 3 + 2] = 0;
  }
  return out;
};

// Pure helper (no refs/closures) so it's callable identically from both the
// geometry-rebuild effect's initial upload and the imperative handle's
// setColors — a single colour-upload code path either way.
const uploadColorAttribute = (mesh: FillMesh | null, nextColors: Float32Array): void => {
  if (!mesh) return;
  const colorAttribute = mesh.geometry.getAttribute('color') as THREE.BufferAttribute;
  (colorAttribute.array as Float32Array).set(nextColors);
  colorAttribute.needsUpdate = true;
};

export const GLFillSurface: React.FC<GLFillSurfaceProps> = ({
  positions,
  colors,
  width,
  height,
  dpr,
  onContextLost,
  handleRef,
}) => {
  // A canvas element can back exactly one WebGLRenderer for its lifetime: once a
  // renderer is constructed (and later force-lost/disposed) on a canvas, calling
  // getContext on that SAME element returns null → THREE throws "reading
  // 'precision'". React StrictMode's mount→cleanup→mount double-invoke (dev) and
  // the 3D→2D switch both re-run the mount effect against the same canvas, so the
  // second construction throws and the whole GL path silently drops to blit —
  // the bug behind "not crisp / distorted tri-points still there" (the user was
  // on the blit fallback the entire time). Fix: on a construction failure,
  // discard the poisoned element by bumping `canvasGen` (React mounts a brand-new
  // <canvas>), and let the mount effect re-run against the clean one. A fresh
  // canvas is never poisoned, so 2–3 attempts suffice; beyond that WebGL is
  // genuinely unavailable → signal the host to use blit.
  const MAX_CANVAS_ATTEMPTS = 3;
  const [canvasGen, setCanvasGen] = useState(0);
  const constructAttemptsRef = useRef(0);
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const rendererRef = useRef<THREE.WebGLRenderer | null>(null);
  const sceneRef = useRef<THREE.Scene | null>(null);
  const cameraRef = useRef<THREE.OrthographicCamera | null>(null);
  const meshRef = useRef<FillMesh | null>(null);
  const contextLostRef = useRef<boolean>(false);
  const widthRef = useRef<number>(width);
  const colorsRef = useRef<Float32Array>(colors);
  const onContextLostRef = useRef<(() => void) | undefined>(onContextLost);

  useEffect(() => {
    widthRef.current = width;
  }, [width]);

  useEffect(() => {
    colorsRef.current = colors;
  }, [colors]);

  useEffect(() => {
    onContextLostRef.current = onContextLost;
  }, [onContextLost]);

  // Renderer, scene, camera, the webglcontextlost listener, and the imperative
  // handle. Keyed on `canvasGen` (not `[]`): the ONLY thing that re-runs this is
  // a construction failure bumping canvasGen to swap in a fresh <canvas> — a
  // WebGL context is a scarce resource (browsers cap ~16 live), so `positions`/
  // size changes must NOT recreate it; those are handled in place below.
  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return undefined;

    // Guard renderer construction: THREE throws (null .precision in
    // WebGLCapabilities) when getContext returns null — which happens when this
    // canvas element was already used by a prior renderer (StrictMode remount /
    // 3D→2D teardown). Recover by discarding the poisoned element: bump canvasGen
    // so React mounts a clean <canvas> and this effect re-runs against it. Only
    // after MAX_CANVAS_ATTEMPTS fresh canvases still fail is WebGL genuinely
    // unavailable → tell the host to use the blit fallback.
    let renderer: THREE.WebGLRenderer;
    try {
      renderer = new THREE.WebGLRenderer({ canvas, antialias: true, alpha: true });
    } catch {
      if (constructAttemptsRef.current < MAX_CANVAS_ATTEMPTS) {
        constructAttemptsRef.current += 1;
        setCanvasGen((g) => g + 1);
        return undefined;
      }
      contextLostRef.current = true;
      onContextLostRef.current?.();
      return undefined;
    }
    // Clean construction — reset the attempt budget so a later independent loss
    // gets its own fresh round of retries.
    constructAttemptsRef.current = 0;
    // Colour parity with the Canvas2D blit fallback: the per-vertex colours are
    // the SAME #rrggbb bytes the blit draws verbatim (already sRGB-encoded, from
    // mapColorCache). THREE r152+ defaults outputColorSpace to sRGB, which would
    // sRGB-encode them a SECOND time → washed-out / gamma-off fills. Setting the
    // output space to linear-sRGB writes the vertex bytes straight to the frame-
    // buffer unchanged, matching the blit exactly.
    renderer.outputColorSpace = THREE.LinearSRGBColorSpace;
    renderer.setPixelRatio(dpr);
    renderer.setSize(width, height, false);

    // Robust ortho depth range for a z=0 plane; the earlier near=-1/far=1 at
    // camera z=1 put the mesh on the far clip boundary.
    const camera = new THREE.OrthographicCamera(0, width, 0, height, 0.1, 1000);
    camera.position.z = 10;

    const scene = new THREE.Scene();

    rendererRef.current = renderer;
    sceneRef.current = scene;
    cameraRef.current = camera;
    contextLostRef.current = false;

    const handleContextLost = (event: Event): void => {
      event.preventDefault();
      contextLostRef.current = true;
      onContextLostRef.current?.();
    };
    canvas.addEventListener('webglcontextlost', handleContextLost);

    const redraw = (): void => {
      if (!rendererRef.current || !sceneRef.current || !cameraRef.current) return;
      rendererRef.current.render(sceneRef.current, cameraRef.current);
    };

    const setTransform = (scale: number, offsetX: number, offsetY: number): void => {
      const currentMesh = meshRef.current;
      if (!currentMesh) return;
      // Compose pan/zoom on top of the X mirror. Screen x = scale*(width - wx)
      // + offsetX (matching VectorOverlay's transform), and worldX =
      // -scale*wx + position.x, so position.x = scale*width + offsetX — NOT
      // width + offsetX, which only coincides at scale 1 and slides every cell
      // off-screen when zoomed.
      currentMesh.scale.set(-scale, scale, 1);
      currentMesh.position.set(scale * widthRef.current + offsetX, offsetY, 0);
      redraw();
    };

    const setColors = (nextColors: Float32Array, _range?: { start: number; count: number }): void => {
      // Full re-upload for this task; true partial GPU upload lands in Task 7.
      uploadColorAttribute(meshRef.current, nextColors);
      redraw();
    };

    const isContextLost = (): boolean => contextLostRef.current;

    if (handleRef) {
      handleRef.current = { setTransform, setColors, redraw, isContextLost };
    }

    return () => {
      // Remove the listener BEFORE forceContextLoss so the loss it triggers
      // doesn't call back into an unmounting parent. dispose() alone does not
      // synchronously release the underlying context; forceContextLoss() does,
      // which matters under the browser's ~16-live-context cap as style/
      // projection toggles mount and unmount this surface repeatedly.
      canvas.removeEventListener('webglcontextlost', handleContextLost);
      renderer.forceContextLoss();
      renderer.dispose();
      rendererRef.current = null;
      sceneRef.current = null;
      cameraRef.current = null;
      if (handleRef) {
        handleRef.current = null;
      }
    };
    // width/height/dpr are read here only as the renderer/camera's initial
    // values — later changes are applied in place by the resize effect below,
    // not by rerunning this effect. `handleRef` is a stable ref object. The sole
    // trigger for a rebuild is `canvasGen` (a fresh canvas after a construction
    // failure); a normal size/positions change must not recreate the context.
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [canvasGen]);

  // Geometry/material/mesh. Keyed on `positions` AND `canvasGen`: a canvasGen
  // bump rebuilds the scene above, so the mesh must be re-added to the NEW scene
  // even when `positions` is unchanged — otherwise the fresh renderer draws an
  // empty scene.
  useEffect(() => {
    const scene = sceneRef.current;
    if (!scene) return undefined;

    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute('position', new THREE.Float32BufferAttribute(expandPositions(positions), 3));
    geometry.setAttribute('color', new THREE.Float32BufferAttribute(colorsRef.current, 3));

    // DoubleSide: this is a flat 2D fill and the X mirror (scale.x = -1) flips
    // triangle winding, so front-face culling would drop every cell. Winding is
    // meaningless for a flat fill — render both sides.
    const material = new THREE.MeshBasicMaterial({ vertexColors: true, side: THREE.DoubleSide });

    const mesh = new THREE.Mesh(geometry, material);
    // Mirror the 2D Canvas2D map's X flip (translate(width,0); scale(-1,1))
    // as a mesh transform; positions themselves stay un-flipped.
    mesh.scale.x = -1;
    mesh.position.x = widthRef.current;

    const previousMesh = meshRef.current;
    if (previousMesh) {
      scene.remove(previousMesh);
    }
    scene.add(mesh);
    meshRef.current = mesh;

    if (rendererRef.current && cameraRef.current) {
      rendererRef.current.render(scene, cameraRef.current);
    }

    return () => {
      scene.remove(mesh);
      geometry.dispose();
      material.dispose();
      if (meshRef.current === mesh) {
        meshRef.current = null;
      }
    };
  }, [positions, canvasGen]);

  // Resize / dpr changes: update the renderer + camera + mirror position
  // in place, without tearing down GPU resources.
  useEffect(() => {
    const renderer = rendererRef.current;
    const camera = cameraRef.current;
    const mesh = meshRef.current;
    if (!renderer || !camera) return;
    renderer.setPixelRatio(dpr);
    renderer.setSize(width, height, false);
    camera.right = width;
    camera.bottom = height;
    camera.updateProjectionMatrix();
    if (mesh) {
      mesh.position.x = width;
    }
    if (sceneRef.current) {
      renderer.render(sceneRef.current, camera);
    }
  }, [width, height, dpr]);

  // key={canvasGen}: a construction failure bumps canvasGen, so React discards
  // the poisoned canvas element and mounts a brand-new one for the retry.
  return <canvas key={canvasGen} ref={canvasRef} style={{ width: '100%', height: '100%', display: 'block' }} />;
};
