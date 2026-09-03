import React, { useEffect, useRef } from 'react';
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

  // Mount-only: renderer, scene, camera, the webglcontextlost listener, and
  // the imperative handle. This is a scarce resource (browsers cap ~16 live
  // WebGL contexts) so it must NOT be recreated on every `positions`/size
  // change — only on unmount. Geometry/material/mesh live in the effect
  // below and get added into this persisted scene.
  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return undefined;

    const renderer = new THREE.WebGLRenderer({ canvas, antialias: true, alpha: true });
    renderer.setPixelRatio(dpr);
    renderer.setSize(width, height, false);

    const camera = new THREE.OrthographicCamera(0, width, 0, height, -1, 1);
    camera.position.z = 1;

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
      // Compose pan/zoom on top of the fixed X mirror.
      currentMesh.scale.set(-scale, scale, 1);
      currentMesh.position.set(widthRef.current + offsetX, offsetY, 0);
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
      canvas.removeEventListener('webglcontextlost', handleContextLost);
      renderer.dispose();
      rendererRef.current = null;
      sceneRef.current = null;
      cameraRef.current = null;
      if (handleRef) {
        handleRef.current = null;
      }
    };
    // Mount-only by design: width/height/dpr are read here only as the
    // renderer/camera's initial values — later changes are applied in place
    // by the resize effect below via setSize/camera bounds, not by rerunning
    // this effect. `handleRef` is a stable ref object supplied by the
    // parent. Re-running this effect would tear down and recreate the WebGL
    // context, which is the exact bug this split fixes.
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  // Geometry/material/mesh only, keyed on `positions`. Rebuilds and swaps
  // the mesh into the persisted scene without touching the renderer/camera.
  useEffect(() => {
    const scene = sceneRef.current;
    if (!scene) return undefined;

    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute('position', new THREE.Float32BufferAttribute(expandPositions(positions), 3));
    geometry.setAttribute('color', new THREE.Float32BufferAttribute(colorsRef.current, 3));

    const material = new THREE.MeshBasicMaterial({ vertexColors: true });

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
  }, [positions]);

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

  return <canvas ref={canvasRef} style={{ width: '100%', height: '100%', display: 'block' }} />;
};
