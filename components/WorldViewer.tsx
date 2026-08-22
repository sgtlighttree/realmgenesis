import React, { useMemo, useRef, useEffect, useLayoutEffect, useState, useCallback } from 'react';
import { Canvas, useFrame, useThree } from '@react-three/fiber';
import { OrbitControls, Stars, Line } from '@react-three/drei';
import * as THREE from 'three';
import { WorldData, ViewMode, Cell, Point, InspectMode, DymaxionSettings, EditMode, LabelVisibility, DEFAULT_LABEL_VISIBILITY, MarkerData } from '../types';
import { getCellColor } from '../utils/colors';
import { seasonalTemperatureDelta } from '../utils/seasons';
import { displayRadius } from '../utils/displayRadius';
import { computeShadeMap, computeContourSegments, contourInterval } from '../utils/shading';
import { computeBorderSegments } from '../utils/borders';
import { collectLabels } from '../utils/labels';
import { computeLabelAnchors } from '../utils/labelAnchors';
import { ScreenOverlay, OverlayTenant } from './overlays/ScreenOverlay';
import { drawCurrentsTenant, drawGraticuleTenant, drawRoutesTenant, drawContoursTenant, drawBordersTenant, drawRiversTenant, drawLabelsTenant, drawRulerTenant, drawDymaxionTenant } from './overlays/tenants';
import { cageEdges } from '../utils/dymaxionCage';
import { computeRiverPolylines } from '../utils/riverPaths';

const Mesh = 'mesh' as any;
const Group = 'group' as any;
const AmbientLight = 'ambientLight' as any;
const PointLight = 'pointLight' as any;
const DirectionalLight = 'directionalLight' as any;
const MeshStandardMaterial = 'meshStandardMaterial' as any;
const InstancedMesh = 'instancedMesh' as any;
const CylinderGeometry = 'cylinderGeometry' as any;
const MeshBasicMaterial = 'meshBasicMaterial' as any;
const LineSegments = 'lineSegments' as any;
const LineBasicMaterial = 'lineBasicMaterial' as any;
const IcosahedronGeometry = 'icosahedronGeometry' as any;
type R3FIntrinsic = React.FC<{ children?: React.ReactNode } & Record<string, unknown>>;
const OctahedronGeometry = 'octahedronGeometry' as unknown as R3FIntrinsic;

// Idle globe spin rate, radians per second.
const SPIN_RATE = 0.05;

// The idle spin lives in its own component for one reason: ORDERING. R3F runs
// useFrame callbacks in subscription order, and React subscribes a child's
// effects before its parent's. While this lived in the parent that also renders
// <ScreenOverlay/>, the overlay's callback ran FIRST — projecting a rotation the
// spin had not yet advanced, one frame behind what WebGL then drew. Mounted as a
// sibling AHEAD of <ScreenOverlay/>, it subscribes first and runs first.
//
// The other half of that fix is in ScreenOverlay, which must also force the world
// matrices current (Three only recomputes them inside render()). Both are needed;
// either alone still leaves the overlay one frame behind the terrain.
const GlobeSpin: React.FC<{ target: React.RefObject<THREE.Group | null>; paused: boolean }> = ({ target, paused }) => {
  useFrame((_state, delta) => {
    if (paused || !target.current) return;
    target.current.rotation.y += delta * SPIN_RATE;
  });
  return null;
};

// Plate widening used to close the seams between cells of differing height.
// 1 = untouched (seams visible, the pre-existing look).
const CELL_OVERHANG = 1.03;

const CityMarkers: React.FC<{ world: WorldData; viewMode: ViewMode; smoothGlobe?: boolean }> = ({ world, viewMode, smoothGlobe = false }) => {
    const capitalsRef = useRef<THREE.InstancedMesh>(null);
    const townsRef = useRef<THREE.InstancedMesh>(null);
    const dummy = useMemo(() => new THREE.Object3D(), []);

    const { capitals, towns } = useMemo(() => {
        const c: Cell[] = [];
        const t: Cell[] = [];
        world.cells.forEach(cell => {
            if (cell.isCapital) c.push(cell);
            else if (cell.isTown) t.push(cell);
        });
        return { capitals: c, towns: t };
    }, [world]);

    useEffect(() => {
        const yAxis = new THREE.Vector3(0, 1, 0);
        const pos = new THREE.Vector3();
        const up = new THREE.Vector3();
        const updateMesh = (mesh: THREE.InstancedMesh, data: Cell[], heightBase: number) => {
            if (!mesh) return;
            data.forEach((cell, i) => {
                const h = displayRadius(cell.height, smoothGlobe);
                pos.set(cell.center.x * h, cell.center.y * h, cell.center.z * h);
                up.copy(pos).normalize();
                dummy.position.copy(pos).addScaledVector(up, heightBase * 0.5);
                dummy.quaternion.setFromUnitVectors(yAxis, up);
                dummy.scale.set(1, 1, 1);
                dummy.updateMatrix();
                mesh.setMatrixAt(i, dummy.matrix);
            });
            mesh.instanceMatrix.needsUpdate = true;
        };
        if (capitalsRef.current) updateMesh(capitalsRef.current, capitals, 0.08);
        if (townsRef.current) updateMesh(townsRef.current, towns, 0.04);
    }, [world, capitals, towns, dummy, smoothGlobe]);

    return (
        <>
            {capitals.length > 0 && (
                <InstancedMesh ref={capitalsRef} args={[undefined, undefined, capitals.length]} visible={true}>
                    <CylinderGeometry args={[0.008, 0.008, 0.08, 6]} />
                    <MeshBasicMaterial color="#ef4444" toneMapped={false} />
                </InstancedMesh>
            )}
            {towns.length > 0 && (
                <InstancedMesh ref={townsRef} args={[undefined, undefined, towns.length]} visible={viewMode === 'political'}>
                    <CylinderGeometry args={[0.005, 0.005, 0.04, 5]} />
                    <MeshBasicMaterial color="#ffffff" toneMapped={false} />
                </InstancedMesh>
            )}
        </>
    );
};

// User-placed POI pins (C4). Markers are sphere-anchored (MarkerData.position),
// not cell-derived, so they're positioned directly rather than via cell.center
// like CityMarkers above. Offset (1.055) sits just under the point-label offset
// (1.09+) so the pin renders under its name sprite without z-fighting.
const MarkerPins: React.FC<{ markers: MarkerData[]; visible: boolean }> = ({ markers, visible }) => {
    const meshRef = useRef<THREE.InstancedMesh>(null);
    const dummy = useMemo(() => new THREE.Object3D(), []);

    useEffect(() => {
        const mesh = meshRef.current;
        if (!mesh) return;
        markers.forEach((marker, i) => {
            const len = Math.hypot(marker.position.x, marker.position.y, marker.position.z) || 1;
            dummy.position.set(
                (marker.position.x / len) * 1.055,
                (marker.position.y / len) * 1.055,
                (marker.position.z / len) * 1.055,
            );
            dummy.rotation.set(0, 0, 0);
            dummy.scale.set(1, 1, 1);
            dummy.updateMatrix();
            mesh.setMatrixAt(i, dummy.matrix);
        });
        mesh.instanceMatrix.needsUpdate = true;
    }, [markers, dummy]);

    if (markers.length === 0) return null;

    return (
        <InstancedMesh ref={meshRef} args={[undefined, undefined, markers.length]} visible={visible}>
            <OctahedronGeometry args={[0.016, 0]} />
            <MeshBasicMaterial color="#f59e0b" toneMapped={false} />
        </InstancedMesh>
    );
};

const CurvedFactionLabel: React.FC<{ name: string; position: THREE.Vector3 }> = ({ name, position }) => {
    const { texture, scale } = useMemo(() => {
        const canvas = document.createElement('canvas');
        const ctx = canvas.getContext('2d');
        const pixelRatio = Math.min(2, window.devicePixelRatio || 1);
        const fontSize = 34 * pixelRatio;
        const paddingX = 16 * pixelRatio;
        const paddingY = 10 * pixelRatio;
        const label = name.toUpperCase();

        if (!ctx) {
            const fallbackTexture = new THREE.CanvasTexture(canvas);
            return { texture: fallbackTexture, scale: [0.2, 0.06, 1] as [number, number, number] };
        }

        ctx.font = `800 ${fontSize}px Inter, ui-sans-serif, system-ui, sans-serif`;
        const textWidth = Math.ceil(ctx.measureText(label).width);
        canvas.width = Math.max(1, textWidth + paddingX * 2);
        canvas.height = Math.max(1, fontSize + paddingY * 2);

        ctx.clearRect(0, 0, canvas.width, canvas.height);
        ctx.font = `800 ${fontSize}px Inter, ui-sans-serif, system-ui, sans-serif`;
        ctx.textAlign = 'center';
        ctx.textBaseline = 'middle';
        ctx.lineJoin = 'round';
        ctx.lineWidth = Math.max(5, 7 * pixelRatio);
        ctx.strokeStyle = '#020617';
        ctx.fillStyle = '#f8fafc';
        ctx.strokeText(label, canvas.width / 2, canvas.height / 2);
        ctx.fillText(label, canvas.width / 2, canvas.height / 2);

        const labelTexture = new THREE.CanvasTexture(canvas);
        labelTexture.minFilter = THREE.LinearFilter;
        labelTexture.magFilter = THREE.LinearFilter;
        labelTexture.generateMipmaps = false;
        labelTexture.needsUpdate = true;

        const height = 0.075;
        const width = height * (canvas.width / canvas.height);
        return { texture: labelTexture, scale: [width, height, 1] as [number, number, number] };
    }, [name]);

    const geometry = useMemo(() => {
        const radius = Math.max(1.08, position.length());
        const center = position.clone().normalize();
        const worldUp = new THREE.Vector3(0, 1, 0);
        let right = new THREE.Vector3().crossVectors(worldUp, center);
        if (right.lengthSq() < 1e-6) {
            right = new THREE.Vector3().crossVectors(new THREE.Vector3(1, 0, 0), center);
        }
        right.normalize();
        const up = new THREE.Vector3().crossVectors(center, right).normalize();
        const width = scale[0];
        const height = scale[1];
        const cols = 28;
        const rows = 6;
        const positions: number[] = [];
        const uvs: number[] = [];
        const indices: number[] = [];

        for (let yStep = 0; yStep <= rows; yStep++) {
            const y = (yStep / rows - 0.5) * height;
            for (let xStep = 0; xStep <= cols; xStep++) {
                const x = (xStep / cols - 0.5) * width;
                const dir = center.clone()
                    .multiplyScalar(radius)
                    .add(right.clone().multiplyScalar(x))
                    .add(up.clone().multiplyScalar(y))
                    .normalize();
                positions.push(dir.x * radius, dir.y * radius, dir.z * radius);
                uvs.push(xStep / cols, yStep / rows);
            }
        }

        for (let yStep = 0; yStep < rows; yStep++) {
            for (let xStep = 0; xStep < cols; xStep++) {
                const a = yStep * (cols + 1) + xStep;
                const b = a + 1;
                const c = a + cols + 1;
                const d = c + 1;
                indices.push(a, b, d, a, d, c);
            }
        }

        const geo = new THREE.BufferGeometry();
        geo.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
        geo.setAttribute('uv', new THREE.Float32BufferAttribute(uvs, 2));
        geo.setIndex(indices);
        return geo;
        // Value-based deps: label positions are recomputed as fresh Vector3s on
        // every world change, but the patch only needs rebuilding when the
        // centroid actually moves
        // eslint-disable-next-line react-hooks/exhaustive-deps
    }, [position.x, position.y, position.z, scale]);

    useEffect(() => () => {
        texture.dispose();
    }, [texture]);

    useEffect(() => () => {
        geometry.dispose();
    }, [geometry]);

    return (
        <Mesh geometry={geometry} renderOrder={7}>
            <MeshBasicMaterial
                map={texture}
                transparent
                depthTest
                depthWrite={false}
                toneMapped={false}
                side={THREE.FrontSide}
                polygonOffset
                polygonOffsetFactor={-1}
                polygonOffsetUnits={-1}
            />
        </Mesh>
    );
};

const CountryLabels: React.FC<{ world: WorldData; visible: boolean }> = ({ world, visible }) => {
    const labels = useMemo(() => {
        if (!world.civData || !visible) return [];
        return world.civData.factions.map(f => {
            let sumX = 0, sumY = 0, sumZ = 0, count = 0;
            const factionCells = world.cells.filter(cell => cell.regionId === f.id);
            const labelCells = factionCells.some(cell => cell.height >= world.params.seaLevel)
                ? factionCells.filter(cell => cell.height >= world.params.seaLevel)
                : factionCells;

            for (const cell of labelCells) {
                sumX += cell.center.x; sumY += cell.center.y; sumZ += cell.center.z;
                count++;
            }
            if (count === 0) return null;
            const position = new THREE.Vector3(sumX / count, sumY / count, sumZ / count).normalize().multiplyScalar(1.09);
            return { id: f.id, name: f.name, position };
        }).filter(Boolean) as {id: number, name: string, position: THREE.Vector3}[];
    }, [world, visible]);

    if (!visible) return null;

    return (
        <Group>
            {labels.map(l => (
                <CurvedFactionLabel
                    key={l.id}
                    name={l.name}
                    position={l.position}
                />
            ))}
        </Group>
    );
};

// The Dymaxion reference cage migrated to the ScreenOverlay `dymaxion` tenant
// (F2 — last named overlay). Geometry lives in utils/dymaxionCage.ts.

const BrushRing: React.FC<{ center: [number, number, number]; radius: number }> = ({ center, radius }) => {
  const geometry = useMemo(() => {
    const normal = new THREE.Vector3(...center).normalize();
    let tangent = new THREE.Vector3(1, 0, 0);
    if (Math.abs(normal.dot(tangent)) > 0.9) tangent = new THREE.Vector3(0, 1, 0);
    const u = tangent.clone().cross(normal).normalize().multiplyScalar(radius);
    const v = normal.clone().cross(u).normalize().multiplyScalar(radius);
    const N = 48;
    const pts: number[] = [];
    for (let i = 0; i < N; i++) {
      const a0 = (i / N) * Math.PI * 2;
      const a1 = ((i + 1) / N) * Math.PI * 2;
      pts.push(
        center[0] + Math.cos(a0) * u.x + Math.sin(a0) * v.x,
        center[1] + Math.cos(a0) * u.y + Math.sin(a0) * v.y,
        center[2] + Math.cos(a0) * u.z + Math.sin(a0) * v.z,
        center[0] + Math.cos(a1) * u.x + Math.sin(a1) * v.x,
        center[1] + Math.cos(a1) * u.y + Math.sin(a1) * v.y,
        center[2] + Math.cos(a1) * u.z + Math.sin(a1) * v.z,
      );
    }
    const geo = new THREE.BufferGeometry();
    geo.setAttribute('position', new THREE.Float32BufferAttribute(pts, 3));
    return geo;
  }, [center, radius]);
  useEffect(() => () => { geometry.dispose(); }, [geometry]);
  return (
    <LineSegments geometry={geometry} renderOrder={10}>
      <LineBasicMaterial color="#ffffff" opacity={0.8} transparent depthTest={false} />
    </LineSegments>
  );
};

const CellSelectionOverlay: React.FC<{ cell: Cell; smoothGlobe?: boolean }> = ({ cell, smoothGlobe = false }) => {
  const { fillGeometry, outlinePoints } = useMemo(() => {
    const hMult = displayRadius(cell.height, smoothGlobe, 0.012);
    const center = new THREE.Vector3(cell.center.x * hMult, cell.center.y * hMult, cell.center.z * hMult);
    const positions: number[] = [];
    const outlinePoints = cell.vertices.map(v => new THREE.Vector3(v.x * hMult, v.y * hMult, v.z * hMult));

    for (let i = 0; i < cell.vertices.length; i++) {
      const next = (i + 1) % cell.vertices.length;
      const v1 = cell.vertices[i];
      const v2 = cell.vertices[next];
      positions.push(
        center.x, center.y, center.z,
        v1.x * hMult, v1.y * hMult, v1.z * hMult,
        v2.x * hMult, v2.y * hMult, v2.z * hMult,
      );
    }

    if (outlinePoints.length > 0) outlinePoints.push(outlinePoints[0].clone());

    const fillGeometry = new THREE.BufferGeometry();
    fillGeometry.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));

    return { fillGeometry, outlinePoints };
    // cell.height mutates in place during terrain painting
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [cell, cell.height, smoothGlobe]);

  useEffect(() => () => { fillGeometry.dispose(); }, [fillGeometry]);

  return (
    <Group>
      <Mesh geometry={fillGeometry} renderOrder={8}>
        <MeshBasicMaterial
          color="#fef08a"
          opacity={0.28}
          transparent
          depthWrite={false}
          side={THREE.DoubleSide}
        />
      </Mesh>
      <Line
        points={outlinePoints}
        color="#facc15"
        lineWidth={3.5}
        transparent
        opacity={1}
        depthWrite={false}
        renderOrder={9}
      />
    </Group>
  );
};

const CellHighlightOutline: React.FC<{ cell: Cell; smoothGlobe?: boolean }> = ({ cell, smoothGlobe = false }) => {
  const geometry = useMemo(() => {
    const hm = displayRadius(cell.height, smoothGlobe, 0.004);
    const verts = cell.vertices;
    const pts: number[] = [];
    for (let i = 0; i < verts.length; i++) {
      const a = verts[i], b = verts[(i + 1) % verts.length];
      pts.push(a.x * hm, a.y * hm, a.z * hm, b.x * hm, b.y * hm, b.z * hm);
    }
    const geo = new THREE.BufferGeometry();
    geo.setAttribute('position', new THREE.Float32BufferAttribute(pts, 3));
    return geo;
    // cell.height mutates in place during terrain painting
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [cell, cell.height, smoothGlobe]);
  useEffect(() => () => { geometry.dispose(); }, [geometry]);
  return (
    <LineSegments geometry={geometry} renderOrder={9}>
      <LineBasicMaterial color="#ffff00" opacity={0.9} transparent depthTest={false} />
    </LineSegments>
  );
};

const WorldMesh: React.FC<{
  world: WorldData,
  viewMode: ViewMode,
  onHover: (cell: Cell | null) => void,
  paused: boolean,
  showGrid: boolean,
  smoothGlobe: boolean,
  showRivers: boolean,
  showRoutes: boolean,
  showHillshade: boolean,
  showContours: boolean,
  showCurrents: boolean,
  showCellEdges: boolean,
  inspectMode: InspectMode;
  onInspect: (cellId: number | null) => void;
  dymaxionSettings: DymaxionSettings;
  editMode: EditMode;
  onPaint: (cellId: number, phase: 'start' | 'stroke' | 'end', isRightClick?: boolean) => void;
  factionColors?: Map<number, string>;
  cultureColors?: Map<number, string>;
  religionColors?: Map<number, string>;
  brushSize: number;
  selectedCellId?: number | null;
  labelVisibility: LabelVisibility;
  rulerArc?: Point[] | null;
}> = ({ world, viewMode, onHover, paused, showGrid, smoothGlobe, showRivers, showRoutes, showHillshade, showContours, showCurrents, showCellEdges, inspectMode, onInspect, dymaxionSettings, editMode, onPaint, factionColors, cultureColors, religionColors, brushSize, selectedCellId = null, labelVisibility, rulerArc = null }) => {
  const spinRef = useRef<THREE.Group>(null);
  const meshRef = useRef<THREE.Mesh>(null);
  const lastUpdate = useRef<number>(0);
  const lastPaintedCell = useRef<number | null>(null);
  const paintPointerActive = useRef(false);
  const pointerStart = useRef<{ x: number; y: number } | null>(null);
  const mapLabels = useMemo(() => collectLabels(world), [world]);
  // Per-label terrain height for the F2 labels tenant (plan §2). Memoized on
  // [mapLabels, world.cells] — NOT smoothGlobe, since height is
  // radius-independent (displayRadius applies smooth/raised at draw time).
  const labelAnchors = useMemo(
    () => computeLabelAnchors(mapLabels, world.cells),
    [mapLabels, world.cells],
  );
  const raycaster = useMemo(() => new THREE.Raycaster(), []);
  const pointer = useMemo(() => new THREE.Vector2(), []);
  const { camera, gl } = useThree();
  const [brushCenter, setBrushCenter] = useState<[number, number, number] | null>(null);
  const [highlightCellId, setHighlightCellId] = useState<number | null>(null);
  

  // The world geometry is allocated once per world structure (cells array
  // identity is stable across paint strokes — App mutates cells in place and
  // shallow-copies WorldData) and refilled in place on data/view changes.
  // No normal attribute: the unlit basic material ignores normals and the
  // standard material uses flatShading, which derives face normals in-shader.
  const geometry = useMemo(() => {
    let triCount = 0;
    for (const cell of world.cells) triCount += cell.vertices.length;
    const geo = new THREE.BufferGeometry();
    const posAttr = new THREE.BufferAttribute(new Float32Array(triCount * 9), 3);
    const colAttr = new THREE.BufferAttribute(new Float32Array(triCount * 9), 3);
    posAttr.setUsage(THREE.DynamicDrawUsage);
    colAttr.setUsage(THREE.DynamicDrawUsage);
    geo.setAttribute('position', posAttr);
    geo.setAttribute('color', colAttr);
    // The globe always fits inside r = 1.05 + margin; a fixed bounding sphere
    // avoids an O(vertices) recomputation on every refill and keeps raycasts valid
    geo.boundingSphere = new THREE.Sphere(new THREE.Vector3(0, 0, 0), 1.1);
    return geo;
  }, [world.cells]);

  useEffect(() => () => { geometry.dispose(); }, [geometry]);

  // Per-cell hillshade relief factor. Keyed on world identity (heights mutate in
  // place on paint, WorldData shallow-copied) so it refreshes exactly like colors.
  const shadeMap = useMemo(
    () => computeShadeMap(world.cells, world.params.seaLevel),
    [world],
  );

  // Refill runs synchronously before the next painted frame so a fresh
  // geometry is never displayed with empty buffers
  useLayoutEffect(() => {
    const posAttr = geometry.getAttribute('position') as THREE.BufferAttribute;
    const colAttr = geometry.getAttribute('color') as THREE.BufferAttribute;
    const pos = posAttr.array as Float32Array;
    const col = colAttr.array as Float32Array;
    // Cell "outlines" are open seams, not lines: neighbours at different radii
    // do not share an edge, so the inner sphere shows through the gap. Widening
    // each plate about its own centre makes neighbours overlap, so the taller
    // one overhangs the shorter — which is what a cliff looks like. Only the
    // rim moves; the centre and hMult are untouched, so the drape invariant
    // (overlay radius == cell radius) is unaffected.
    const inflate = showCellEdges ? 1 : CELL_OVERHANG;
    let o = 0;
    for (const cell of world.cells) {
      const c = getCellColor(cell, viewMode, world.params.seaLevel, factionColors, cultureColors, religionColors, seasonalTemperatureDelta(cell, world.params));
      // Multiply in relief shading only when toggled; the off path stays a
      // straight color copy so rendering is unchanged.
      if (showHillshade) c.multiplyScalar(shadeMap[cell.id]);
      const hMult = displayRadius(cell.height, smoothGlobe);
      const uc = cell.center;
      const cx = uc.x * hMult; const cy = uc.y * hMult; const cz = uc.z * hMult;
      for (let i = 0; i < cell.vertices.length; i++) {
        const v1 = cell.vertices[i]; const v2 = cell.vertices[(i + 1) % cell.vertices.length];
        pos[o] = cx; pos[o + 1] = cy; pos[o + 2] = cz;
        pos[o + 3] = (uc.x + (v1.x - uc.x) * inflate) * hMult;
        pos[o + 4] = (uc.y + (v1.y - uc.y) * inflate) * hMult;
        pos[o + 5] = (uc.z + (v1.z - uc.z) * inflate) * hMult;
        pos[o + 6] = (uc.x + (v2.x - uc.x) * inflate) * hMult;
        pos[o + 7] = (uc.y + (v2.y - uc.y) * inflate) * hMult;
        pos[o + 8] = (uc.z + (v2.z - uc.z) * inflate) * hMult;
        col[o] = c.r; col[o + 1] = c.g; col[o + 2] = c.b;
        col[o + 3] = c.r; col[o + 4] = c.g; col[o + 5] = c.b;
        col[o + 6] = c.r; col[o + 7] = c.g; col[o + 8] = c.b;
        o += 9;
      }
    }
    posAttr.needsUpdate = true;
    colAttr.needsUpdate = true;
  }, [geometry, world, viewMode, factionColors, cultureColors, religionColors, showHillshade, shadeMap, smoothGlobe, showCellEdges]);

  const faceMap = useMemo(() => {
     const map: number[] = [];
     world.cells.forEach(cell => {
       for(let i=0; i<cell.vertices.length; i++) map.push(cell.id);
     });
     return map;
  }, [world.cells]);

  const getTriangleIndex = useCallback((e: any) => {
      if (e.faceIndex !== undefined && e.faceIndex !== null) return e.faceIndex;
      if (e.face && e.face.a !== undefined && e.face.a !== null) return Math.floor(e.face.a / 3);
      return null;
  }, []);

  const isPaintMode = editMode !== 'off' && editMode !== 'world-edit';
  const tracksPointerMove = !isPaintMode && inspectMode === 'hover';
  const handlesSelectionClick = !isPaintMode && (inspectMode === 'click' || editMode === 'world-edit');

  const pickNativeCellId = useCallback((e: PointerEvent) => {
      const mesh = meshRef.current;
      if (!mesh) return null;

      const canvas = gl.domElement;
      const rect = canvas.getBoundingClientRect();
      pointer.set(
          ((e.clientX - rect.left) / rect.width) * 2 - 1,
          -((e.clientY - rect.top) / rect.height) * 2 + 1,
      );
      camera.updateMatrixWorld();
      raycaster.setFromCamera(pointer, camera);

      const hit = raycaster.intersectObject(mesh, false)[0];
      if (!hit || hit.faceIndex === undefined || hit.faceIndex === null) return null;

      return faceMap[hit.faceIndex] ?? null;
  }, [gl.domElement, pointer, camera, raycaster, faceMap]);

  const showBrushCell = useCallback((cellId: number) => {
      const cell = world.cells[cellId];
      if (!cell) return;
      const hMult = displayRadius(cell.height, smoothGlobe, 0.005);
      setBrushCenter([cell.center.x * hMult, cell.center.y * hMult, cell.center.z * hMult]);
      setHighlightCellId(cellId);
  }, [world.cells, smoothGlobe]);

  useEffect(() => {
      if (tracksPointerMove) return;
      setBrushCenter(null);
      setHighlightCellId(null);
      onHover(null);
  }, [tracksPointerMove, onHover]);

  useEffect(() => {
      if (!handlesSelectionClick || dymaxionSettings.mode === 'overlay') return;

      const canvas = gl.domElement;

      const handleNativePointerDown = (e: PointerEvent) => {
          if (e.button !== 0) return;
          pointerStart.current = { x: e.clientX, y: e.clientY };
      };

      const handleNativePointerUp = (e: PointerEvent) => {
          const start = pointerStart.current;
          pointerStart.current = null;
          if (!start || e.button !== 0) return;

          const dx = e.clientX - start.x;
          const dy = e.clientY - start.y;
          if (Math.hypot(dx, dy) > 4) return;

          const cellId = pickNativeCellId(e);
          if (cellId !== null) onInspect(cellId);
      };

      canvas.addEventListener('pointerdown', handleNativePointerDown);
      canvas.addEventListener('pointerup', handleNativePointerUp);

      return () => {
          canvas.removeEventListener('pointerdown', handleNativePointerDown);
          canvas.removeEventListener('pointerup', handleNativePointerUp);
      };
  }, [handlesSelectionClick, dymaxionSettings.mode, gl.domElement, pickNativeCellId, onInspect]);

  useEffect(() => {
      if (!isPaintMode || dymaxionSettings.mode === 'overlay') {
          paintPointerActive.current = false;
          lastPaintedCell.current = null;
          setBrushCenter(null);
          setHighlightCellId(null);
          return;
      }

      const canvas = gl.domElement;

      const finishPaint = () => {
          if (paintPointerActive.current && lastPaintedCell.current !== null) {
              onPaint(lastPaintedCell.current, 'end');
          }
          paintPointerActive.current = false;
          lastPaintedCell.current = null;
          setBrushCenter(null);
          setHighlightCellId(null);
      };

      const handleNativePointerDown = (e: PointerEvent) => {
          const isRight = e.button === 2;
          if (e.button !== 0 && !(isRight && editMode === 'terrain-flatten')) return;

          const cellId = pickNativeCellId(e);
          if (cellId === null) return;

          e.preventDefault();
          showBrushCell(cellId);

          if (isRight) {
              onPaint(cellId, 'start', true);
              return;
          }

          paintPointerActive.current = true;
          lastPaintedCell.current = cellId;
          onPaint(cellId, 'start');
          onPaint(cellId, 'stroke');
      };

      const handleNativePointerMove = (e: PointerEvent) => {
          if (!paintPointerActive.current) return;
          const now = Date.now();
          if (now - lastUpdate.current < 50) return;
          lastUpdate.current = now;

          const cellId = pickNativeCellId(e);
          if (cellId === null) return;

          showBrushCell(cellId);
          lastPaintedCell.current = cellId;
          onPaint(cellId, 'stroke');
      };

      const handleContextMenu = (e: MouseEvent) => {
          e.preventDefault();
      };

      canvas.addEventListener('pointerdown', handleNativePointerDown);
      window.addEventListener('pointermove', handleNativePointerMove);
      window.addEventListener('pointerup', finishPaint);
      canvas.addEventListener('contextmenu', handleContextMenu);

      return () => {
          canvas.removeEventListener('pointerdown', handleNativePointerDown);
          window.removeEventListener('pointermove', handleNativePointerMove);
          window.removeEventListener('pointerup', finishPaint);
          canvas.removeEventListener('contextmenu', handleContextMenu);
      };
  }, [isPaintMode, editMode, dymaxionSettings.mode, gl.domElement, pickNativeCellId, showBrushCell, onPaint]);

  const handlePointerMove = useCallback((e: any) => {
      if (dymaxionSettings.mode === 'overlay') return;
      const now = Date.now();
      if (now - lastUpdate.current < 50) return;
      lastUpdate.current = now;

      const triIndex = getTriangleIndex(e);
      if (triIndex === null) { setBrushCenter(null); return; }
      const cellId = faceMap[triIndex];
      if (cellId === undefined) return;

      const cell = world.cells[cellId];
      if (inspectMode === 'hover') onHover(cell);
  }, [inspectMode, faceMap, world.cells, onHover, getTriangleIndex, dymaxionSettings.mode]);

  const handlePointerOut = useCallback(() => {
      setBrushCenter(null);
      setHighlightCellId(null);
      if (inspectMode === 'hover' && !isPaintMode) onHover(null);
  }, [inspectMode, isPaintMode, onHover]);

  // Isolines are an O(cells x neighbors) sweep, so they are memoized here and
  // closed over by the tenant rather than recomputed on every overlay redraw.
  // Keyed on world identity: heights mutate in place on paint and WorldData is
  // shallow-copied, so isolines must recompute per stroke.
  const contourSegments = useMemo(
    () => (showContours ? computeContourSegments(world.cells, world.params.seaLevel, contourInterval(world.cells, world.params.seaLevel)) : []),
    [world, showContours],
  );

  // Faction border edges. Same deal as contour segments: the extraction is
  // O(cells x neighbors x vertices^2), so it is memoized on world identity and
  // closed over by the tenant. Regions mutate in place on a civ edit and
  // WorldData is shallow-copied, so borders recompute per edit.
  const borderSegments = useMemo(
    () => (labelVisibility.borders && world.civData ? computeBorderSegments(world.cells) : []),
    [world, labelVisibility.borders],
  );

  // River polylines. Keyed on world.rivers (stable across paint strokes), NOT
  // world identity — the CatmullRom smoothing in computeRiverPolylines must not
  // re-run per paint stroke, same contract the old RiverLines kept.
  const riverPolylines = useMemo(
    () => (showRivers && world.rivers ? computeRiverPolylines(world.rivers) : []),
    [world.rivers, showRivers],
  );

  // F2 screen-space overlay tenants. Painting order is array order: rivers sit
  // under everything (ROUTE_LIFT's comment notes routes sit "over rivers"), then
  // contours (terrain annotation), then currents, then routes, then borders,
  // then the graticule, with LABELS LAST (plan §5) so they sit on top of every
  // other overlay.
  // Rotated cage edges for the dymaxion tenant. Recomputed only when the cage
  // orientation changes (the tenant's id also encodes lon/lat/roll so the
  // ScreenOverlay redraw gate fires — see the note in the tenant array below).
  const dymaxionEdges = useMemo(() => cageEdges(dymaxionSettings), [dymaxionSettings]);

  const overlayTenants = useMemo<OverlayTenant[]>(() => [
    { id: 'rivers', visible: showRivers && riverPolylines.length > 0,
      draw: (ctx, _proj, _world, project, smooth) => drawRiversTenant(ctx, riverPolylines, project, smooth) },
    { id: 'contours', visible: showContours && contourSegments.length > 0,
      draw: (ctx, _proj, _world, project, smooth) => drawContoursTenant(ctx, contourSegments, project, smooth) },
    { id: 'currents', visible: showCurrents && !!world.currents, draw: drawCurrentsTenant },
    // Routes above the current field so dashed sea routes read over the arrows.
    { id: 'routes', visible: showRoutes && !!world.routes, draw: drawRoutesTenant },
    // Borders above the routes they cut across, below the reference grid.
    { id: 'borders', visible: labelVisibility.borders && borderSegments.length > 0,
      draw: (ctx, _proj, _world, project, smooth) => drawBordersTenant(ctx, borderSegments, project, smooth) },
    { id: 'graticule', visible: showGrid, draw: drawGraticuleTenant },
    // Dymaxion reference cage. Fixed radius (clears peaks), not draped; see
    // drawDymaxionTenant. lon/lat/roll are folded into the id so the redraw
    // gate fires when the cage is rotated with the globe paused — the same
    // many-state-in-the-draw-closure trap the labels tenant hit (S22).
    { id: `dymaxion:${dymaxionSettings.lon},${dymaxionSettings.lat},${dymaxionSettings.roll}`,
      visible: dymaxionSettings.showOverlay && dymaxionEdges.length > 0,
      draw: (ctx, _proj, _world, project, smooth) => drawDymaxionTenant(ctx, dymaxionEdges, project, smooth) },
    // Labels last, on top of everything (plan §5). camDist is read live at
    // draw time (camera is a stable reference from useThree, not a memo dep)
    // so it never goes stale between redraws — same as the old PointLabels.
    //
    // ScreenOverlay's redraw gate is keyed on `active tenant ids + camera/globe
    // matrix` (see ScreenOverlay.tsx) — it has no idea what a tenant's `draw`
    // closure reads. Every OTHER visibility toggle in this array (borders,
    // routes, contours...) is one boolean mapped straight to `visible`, so
    // toggling it changes which ids are in `active` and the key changes with
    // it. `drawLabelsTenant` instead multiplexes FIVE independent toggles
    // (capitals/towns/provinces/geography/markers) through `labelVisibility`
    // inside the draw body, which the key can't see — so a toggle silently
    // no-ops until the next unrelated redraw (e.g. an orbit move). Folding the
    // flags into the id is the id-based idiom this seam already uses, applied
    // to a many-flags tenant instead of a one-flag tenant.
    { id: `labels:${[labelVisibility.capitals, labelVisibility.towns, labelVisibility.provinces, labelVisibility.geography, labelVisibility.markers].map(Number).join('')}`,
      visible: labelAnchors.length > 0,
      draw: (ctx, _proj, _world, project, smooth) =>
        drawLabelsTenant(ctx, labelAnchors, project, smooth, camera.position.length(), labelVisibility) },
    // Ruler on top of everything — it is an active measurement affordance. Fixed
    // radius (clears peaks), not draped; see drawRulerTenant.
    { id: 'ruler', visible: !!rulerArc && rulerArc.length > 1,
      draw: (ctx, _proj, _world, project, smooth) => drawRulerTenant(ctx, rulerArc!, project, smooth) },
  ], [showRivers, riverPolylines, showContours, contourSegments, showCurrents, world.currents, showRoutes, world.routes, labelVisibility, borderSegments, showGrid, dymaxionSettings.showOverlay, dymaxionSettings.lon, dymaxionSettings.lat, dymaxionSettings.roll, dymaxionEdges, labelAnchors, camera, rulerArc]);

  return (
    <Group>
        <GlobeSpin target={spinRef} paused={paused} />
        <Group ref={spinRef}>
            <Mesh
            ref={meshRef}
            name="globe-mesh"
            geometry={geometry}
            onPointerMove={tracksPointerMove ? handlePointerMove : undefined}
            onPointerOut={tracksPointerMove ? handlePointerOut : undefined}
            >
                {viewMode === 'political' ? (
                    <MeshBasicMaterial vertexColors toneMapped={false} side={THREE.FrontSide} />
                ) : (
                    <MeshStandardMaterial vertexColors roughness={0.8} metalness={0.1} flatShading side={THREE.FrontSide} />
                )}
                <CityMarkers world={world} viewMode={viewMode} smoothGlobe={smoothGlobe} />
                <MarkerPins markers={world.markers ?? []} visible={labelVisibility.markers} />
                <React.Suspense fallback={null}>
                    {/* CurvedFactionLabel stays 3D (curved textured meshes) — F2 does not
                        flatten faction names to Canvas2D (plan §1). */}
                    <CountryLabels world={world} visible={labelVisibility.factions} />
                </React.Suspense>
                {/* Faction borders migrated to ScreenOverlay (F2 borders tenant). */}
                {/* Rivers migrated to ScreenOverlay (F2 rivers tenant). */}
                {/* Roads & sea routes migrated to ScreenOverlay (F2 routes tenant). */}
                {/* Contour lines migrated to ScreenOverlay (F2 contours tenant). */}
                {/* Lat/long grid migrated to ScreenOverlay (F2 graticule tenant). */}
                {/* Point labels (capitals/provinces/towns/geography/markers) migrated
                    to ScreenOverlay (F2 labels tenant, S22 — last tenant). */}
                {showGrid && <TiltAxisLine radius={1.35} />}
                {/* Cell highlight outline */}
                {highlightCellId !== null && world.cells[highlightCellId] && (
                    <CellHighlightOutline cell={world.cells[highlightCellId]} smoothGlobe={smoothGlobe} />
                )}
                {selectedCellId !== null && (() => {
                    const selectedCell = world.cells[selectedCellId];
                    return selectedCell ? <CellSelectionOverlay cell={selectedCell} smoothGlobe={smoothGlobe} /> : null;
                })()}
                {/* Ruler arc migrated to ScreenOverlay (F2 ruler tenant). */}
            </Mesh>
            {/* Brush size ring */}
            {isPaintMode && brushCenter && (
                <BrushRing center={brushCenter} radius={Math.max(0.018, brushSize * 0.038)} />
            )}
            {/* Dymaxion cage migrated to ScreenOverlay (F2 dymaxion tenant). */}
            <Mesh scale={[0.99, 0.99, 0.99]}>
                <IcosahedronGeometry args={[1, 16]} />
                <MeshBasicMaterial color="#000000" side={THREE.FrontSide} />
            </Mesh>
        </Group>
        <ScreenOverlay world={world} tenants={overlayTenants} smoothGlobe={smoothGlobe} />
    </Group>
  );
};

const WorldViewer: React.FC<{ world: WorldData | null; viewMode: ViewMode; showGrid?: boolean; smoothGlobe?: boolean; showRivers?: boolean; showRoutes?: boolean; showHillshade?: boolean; showContours?: boolean; showCurrents?: boolean; showCellEdges?: boolean; labelVisibility?: LabelVisibility; inspectMode: InspectMode; onInspect: (cellId: number | null) => void; selectedCellId?: number | null; dymaxionSettings: DymaxionSettings; onDymaxionChange: React.Dispatch<React.SetStateAction<DymaxionSettings>>; editMode: EditMode; onPaint: (cellId: number, phase: 'start' | 'stroke' | 'end', isRightClick?: boolean) => void; factionColors?: Map<number, string>; cultureColors?: Map<number, string>; religionColors?: Map<number, string>; brushSize?: number; rulerArc?: Point[] | null; overlayClassName?: string; paused?: boolean; onPausedChange?: (v: boolean) => void; showPauseControl?: boolean; }> = ({ world, viewMode, showGrid = false, smoothGlobe = false, showRivers = true, showRoutes = false, showHillshade = false, showContours = false, showCurrents = false, showCellEdges = false, labelVisibility = DEFAULT_LABEL_VISIBILITY, inspectMode, onInspect, selectedCellId = null, dymaxionSettings, onDymaxionChange, editMode, onPaint, factionColors, cultureColors, religionColors, brushSize = 1, rulerArc = null, overlayClassName = 'absolute top-4 right-4 z-overlay flex gap-2', paused: pausedProp, onPausedChange, showPauseControl = true }) => {
  const [hoveredCell, setHoveredCell] = useState<Cell | null>(null);
  // Rotation pause is controlled-OPTIONAL, the same contract as a native input:
  // pass `paused` + `onPausedChange` to own it from outside (the shell lifts it
  // so the control can live in the top strip), or pass neither and this keeps
  // its own state (classic App). Not a `bare`-style personality flag — the host
  // chooses where the STATE lives, and the rendering follows from that.
  const [internalPaused, setInternalPaused] = useState(false);
  const isPauseControlled = pausedProp !== undefined;
  const paused = isPauseControlled ? pausedProp : internalPaused;
  const setPaused = useCallback((v: boolean) => {
    if (!isPauseControlled) setInternalPaused(v);
    onPausedChange?.(v);
  }, [isPauseControlled, onPausedChange]);
  const [isSpaceHeld, setIsSpaceHeld] = useState(false);
  const dragRef = useRef<{ active: boolean; lastX: number; lastY: number }>({ active: false, lastX: 0, lastY: 0 });
  const overlayMode = dymaxionSettings.mode === 'overlay';
  const isPaintModeActive = editMode !== 'off' && editMode !== 'world-edit';
  const orbitMouseButtons = useMemo(() => ({
    LEFT: isPaintModeActive && !isSpaceHeld ? -1 : THREE.MOUSE.ROTATE,
    MIDDLE: THREE.MOUSE.ROTATE,
    RIGHT: -1,
  }), [isPaintModeActive, isSpaceHeld]);

  useEffect(() => {
    const onDown = (e: KeyboardEvent) => { if (e.code === 'Space' && !e.repeat) { e.preventDefault(); setIsSpaceHeld(true); } };
    const onUp = (e: KeyboardEvent) => { if (e.code === 'Space') setIsSpaceHeld(false); };
    window.addEventListener('keydown', onDown);
    window.addEventListener('keyup', onUp);
    return () => { window.removeEventListener('keydown', onDown); window.removeEventListener('keyup', onUp); };
  }, []);

  useEffect(() => {
    if (overlayMode) setPaused(true);
  }, [overlayMode, setPaused]);

  const wrapAngle = useCallback((v: number) => {
    let r = ((v + 180) % 360 + 360) % 360 - 180;
    if (r === -180) r = 180;
    return r;
  }, []);

  const clampLat = useCallback((v: number) => Math.max(-90, Math.min(90, v)), []);

  const handleOverlayPointerDown = useCallback((e: any) => {
    if (!overlayMode || !dymaxionSettings.showOverlay) return;
    dragRef.current = { active: true, lastX: e.clientX, lastY: e.clientY };
  }, [overlayMode, dymaxionSettings.showOverlay]);

  const handleOverlayPointerMove = useCallback((e: any) => {
    if (!dragRef.current.active) return;
    const dx = e.clientX - dragRef.current.lastX;
    const dy = e.clientY - dragRef.current.lastY;
    dragRef.current.lastX = e.clientX;
    dragRef.current.lastY = e.clientY;
    const sensitivity = 0.25;
    onDymaxionChange((prev) => ({
      ...prev,
      // Grab-and-drag feel, derived from the cage euler Euler(lat, -lon, roll):
      // longitude drives a -lon rotation about the up axis, so dragging right
      // (+dx) must DECREASE lon for the front to move right; dragging up (dy<0)
      // must pitch the front up, which is lat += dy. Both were inverted before.
      lon: e.shiftKey ? prev.lon : wrapAngle(prev.lon - dx * sensitivity),
      lat: e.shiftKey ? prev.lat : clampLat(prev.lat + dy * sensitivity),
      roll: e.shiftKey ? wrapAngle(prev.roll - dx * sensitivity) : prev.roll,
    }));
  }, [onDymaxionChange, clampLat, wrapAngle]);

  const handleOverlayPointerUp = useCallback(() => {
    dragRef.current.active = false;
  }, []);

  useEffect(() => {
    if (inspectMode !== 'hover') return;
    if (hoveredCell) onInspect(hoveredCell.id);
    else onInspect(null);
  }, [hoveredCell, inspectMode, onInspect]);

  return (
    <div className="w-full h-full bg-black relative group">
      <Canvas
        camera={{ position: [0, 0, 2.5], fov: 45 }}
        onPointerDown={handleOverlayPointerDown}
        onPointerMove={handleOverlayPointerMove}
        onPointerUp={handleOverlayPointerUp}
        onPointerLeave={handleOverlayPointerUp}
      >
        <AmbientLight intensity={0.5} />
        <PointLight position={[10, 10, 10]} intensity={1.5} />
        <DirectionalLight position={[-5, 5, 2]} intensity={0.5} />
        <Stars radius={100} depth={50} count={5000} factor={4} saturation={0} fade speed={1} />
        {world && (
          <Group rotation={[0,0, (world.params.axialTilt || 0) * Math.PI / 180]}>
             <WorldMesh
               world={world}
               viewMode={viewMode}
               onHover={setHoveredCell}
               paused={paused}
               showGrid={showGrid}
               smoothGlobe={smoothGlobe}
               showRivers={showRivers}
               showRoutes={showRoutes}
               showHillshade={showHillshade}
               showContours={showContours}
               showCurrents={showCurrents}
               showCellEdges={showCellEdges}
               labelVisibility={labelVisibility}
               inspectMode={inspectMode}
               onInspect={onInspect}
               selectedCellId={selectedCellId}
               dymaxionSettings={dymaxionSettings}
               editMode={editMode}
               onPaint={onPaint}
               factionColors={factionColors}
               cultureColors={cultureColors}
               religionColors={religionColors}
               brushSize={brushSize}
               rulerArc={rulerArc}
             />
          </Group>
        )}
        <OrbitControls
          enablePan={false}
          minDistance={1.2}
          maxDistance={6}
          enableRotate={!overlayMode}
          mouseButtons={orbitMouseButtons}
        />
      </Canvas>
      {!world && <div className="absolute inset-0 flex items-center justify-center text-ink-strong/50">Forging World...</div>}
      
      {showPauseControl && (
        <div className={overlayClassName}>
           <button
             onClick={() => { setPaused(!paused); }}
             disabled={overlayMode}
             aria-label={paused ? 'Resume globe rotation' : 'Pause globe rotation'}
             className={`bg-surface-raised/80 text-ink-strong p-2 backdrop-blur border border-white/10 shadow-lg ${overlayMode ? 'opacity-40 cursor-not-allowed' : 'hover:bg-surface-hover'}`}
           >
             {paused ? "▶" : "⏸"}
           </button>
        </div>
      )}
    </div>
  );
};

export default WorldViewer;

// The planet's rotation axis (pole-to-pole), rendered inside the tilted world
// group so it visibly leans by `axialTilt`. Along local Y, so it stays fixed as
// the globe spins around it. depthTest ON → the globe occludes the back half
// (reads as a real 3D axis, not an overlay pasted on top); the two poles poke
// out past the sphere so the tilt stays legible while navigating.
const TiltAxisLine: React.FC<{ radius?: number }> = ({ radius = 1.35 }) => {
  const geometry = useMemo(() => {
      const geo = new THREE.BufferGeometry();
      geo.setAttribute('position', new THREE.Float32BufferAttribute([0, -radius, 0, 0, radius, 0], 3));
      return geo;
  }, [radius]);
  useEffect(() => () => { geometry.dispose(); }, [geometry]);
  return (
      <LineSegments geometry={geometry}>
          <LineBasicMaterial color="#ffd166" opacity={0.85} transparent depthTest={true} />
      </LineSegments>
  );
};
