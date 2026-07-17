import React, { useMemo, useRef, useEffect, useLayoutEffect, useState, useCallback } from 'react';
import { Canvas, useFrame, useThree } from '@react-three/fiber';
import { OrbitControls, Stars, Line } from '@react-three/drei';
import * as THREE from 'three';
import { WorldData, ViewMode, Cell, Point, InspectMode, DymaxionSettings, EditMode } from '../types';
import { getCellColor } from '../utils/colors';

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

const CityMarkers: React.FC<{ world: WorldData; viewMode: ViewMode }> = ({ world, viewMode }) => {
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
                const h = 1 + (cell.height * 0.05);
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
    }, [world, capitals, towns, dummy]);

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

const RiverLines: React.FC<{ world: WorldData; visible: boolean }> = ({ world, visible }) => {
    // Keyed on world.rivers (stable across paint strokes), not world identity,
    // so painting never re-runs the CatmullRom smoothing or reallocates buffers
    const rivers = world.rivers;
    const geometry = useMemo(() => {
        if (!rivers || !visible) return null;

        const positions: number[] = [];

        // Batch all river segments into a single LineSegments geometry for performance
        // Rendering thousands of individual <Line> components causes massive overhead/freezes
        rivers.forEach(path => {
            if (path.length < 2) return;
            
            // Create Curve for smoothing
            const vectors = path.map(p => new THREE.Vector3(p.x, p.y, p.z));
            const curve = new THREE.CatmullRomCurve3(vectors);
            
            // Adaptive sampling based on length, but simple count is safer for perf
            const points = curve.getPoints(Math.min(50, vectors.length * 4));
            
            for (let i = 0; i < points.length - 1; i++) {
                positions.push(points[i].x, points[i].y, points[i].z);
                positions.push(points[i+1].x, points[i+1].y, points[i+1].z);
            }
        });

        if (positions.length === 0) return null;

        const geo = new THREE.BufferGeometry();
        geo.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
        return geo;
    }, [rivers, visible]);

    useEffect(() => () => { geometry?.dispose(); }, [geometry]);

    if (!visible || !geometry) return null;

    return (
        <LineSegments geometry={geometry}>
            <LineBasicMaterial color="#38bdf8" opacity={0.8} transparent linewidth={1.5} />
        </LineSegments>
    );
}

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

const FactionBorders: React.FC<{ world: WorldData; visible: boolean }> = ({ world, visible }) => {
  const geometry = useMemo(() => {
      if (!visible || !world.civData) return null;

      const positions: number[] = [];
      const threshold = 0.000001; 

      // Iterate unique pairs of neighbors to find borders
      world.cells.forEach(cellA => {
          cellA.neighbors.forEach(nId => {
              const cellB = world.cells[nId];
              if (!cellB || cellA.id >= cellB.id) return; // Process pair once
              
              const rA = cellA.regionId;
              const rB = cellB.regionId;
              
              // Draw border if regions are different
              // This includes border between Faction A and Faction B
              // AND border between Faction A and Unclaimed (International Waters)
              if (rA !== rB) {
                  // Find shared vertices between cellA and cellB to define the edge
                  const shared: Point[] = [];
                  for (const vA of cellA.vertices) {
                      for (const vB of cellB.vertices) {
                          const distSq = (vA.x - vB.x)**2 + (vA.y - vB.y)**2 + (vA.z - vB.z)**2;
                          if (distSq < threshold) {
                              shared.push(vA);
                              break; 
                          }
                      }
                      if (shared.length === 2) break;
                  }
                  
                  if (shared.length === 2) {
                      const hA = 1 + (cellA.height * 0.05);
                      const hB = 1 + (cellB.height * 0.05);
                      // Slight offset to prevent z-fighting with mesh
                      const h = Math.max(hA, hB) + 0.002; 
                      
                      const p1 = shared[0];
                      const p2 = shared[1];
                      
                      positions.push(p1.x * h, p1.y * h, p1.z * h);
                      positions.push(p2.x * h, p2.y * h, p2.z * h);
                  }
              }
          });
      });
      
      if (positions.length === 0) return null;
      
      const geo = new THREE.BufferGeometry();
      geo.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
      return geo;
  }, [world, visible]);

  useEffect(() => () => { geometry?.dispose(); }, [geometry]);

  if (!geometry) return null;

  return (
    <LineSegments geometry={geometry}>
      <LineBasicMaterial color="white" linewidth={1} opacity={0.8} transparent depthTest={true} />
    </LineSegments>
  );
};

const DymaxionOverlay: React.FC<{ settings: DymaxionSettings }> = ({ settings }) => {
  const { faceGeometry, edgeGeometry } = useMemo(() => {
    const faceGeometry = new THREE.IcosahedronGeometry(1.12, 0);
    const edgeGeometry = new THREE.EdgesGeometry(faceGeometry);
    return { faceGeometry, edgeGeometry };
  }, []);

  useEffect(() => () => { faceGeometry.dispose(); edgeGeometry.dispose(); }, [faceGeometry, edgeGeometry]);

  const rotation = useMemo(() => {
    const lon = THREE.MathUtils.degToRad(settings.lon);
    const lat = THREE.MathUtils.degToRad(settings.lat);
    const roll = THREE.MathUtils.degToRad(settings.roll);
    return new THREE.Euler(lat, -lon, roll, 'YXZ');
  }, [settings.lon, settings.lat, settings.roll]);

  return (
    <Group rotation={rotation}>
      <Mesh geometry={faceGeometry} renderOrder={5}>
        <MeshBasicMaterial
          color="#fbbf24"
          opacity={0.18}
          transparent
          depthWrite={false}
          depthTest={false}
          side={THREE.DoubleSide}
          polygonOffset
          polygonOffsetFactor={-2}
          polygonOffsetUnits={-2}
        />
      </Mesh>
      <LineSegments geometry={edgeGeometry} renderOrder={6}>
        <LineBasicMaterial color="#fbbf24" linewidth={1} opacity={0.95} transparent depthTest={false} />
      </LineSegments>
    </Group>
  );
};

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

const CellSelectionOverlay: React.FC<{ cell: Cell }> = ({ cell }) => {
  const { fillGeometry, outlinePoints } = useMemo(() => {
    const hMult = 1 + cell.height * 0.05 + 0.012;
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
  }, [cell, cell.height]);

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

const CellHighlightOutline: React.FC<{ cell: Cell }> = ({ cell }) => {
  const geometry = useMemo(() => {
    const hm = 1 + cell.height * 0.05 + 0.004;
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
  }, [cell, cell.height]);
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
  showRivers: boolean,
  inspectMode: InspectMode;
  onInspect: (cellId: number | null) => void;
  dymaxionSettings: DymaxionSettings;
  editMode: EditMode;
  onPaint: (cellId: number, phase: 'start' | 'stroke' | 'end', isRightClick?: boolean) => void;
  factionColors?: Map<number, string>;
  brushSize: number;
  selectedCellId?: number | null;
  showFactionOverlay: boolean;
}> = ({ world, viewMode, onHover, paused, showGrid, showRivers, inspectMode, onInspect, dymaxionSettings, editMode, onPaint, factionColors, brushSize, selectedCellId = null, showFactionOverlay }) => {
  const spinRef = useRef<THREE.Group>(null);
  const meshRef = useRef<THREE.Mesh>(null);
  const lastUpdate = useRef<number>(0);
  const lastPaintedCell = useRef<number | null>(null);
  const paintPointerActive = useRef(false);
  const pointerStart = useRef<{ x: number; y: number } | null>(null);
  const raycaster = useMemo(() => new THREE.Raycaster(), []);
  const pointer = useMemo(() => new THREE.Vector2(), []);
  const { camera, gl } = useThree();
  const [brushCenter, setBrushCenter] = useState<[number, number, number] | null>(null);
  const [highlightCellId, setHighlightCellId] = useState<number | null>(null);
  
  useFrame((state, delta) => {
    if (!paused) {
        if (spinRef.current) spinRef.current.rotation.y += delta * 0.05;
    }
  });

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

  // Refill runs synchronously before the next painted frame so a fresh
  // geometry is never displayed with empty buffers
  useLayoutEffect(() => {
    const posAttr = geometry.getAttribute('position') as THREE.BufferAttribute;
    const colAttr = geometry.getAttribute('color') as THREE.BufferAttribute;
    const pos = posAttr.array as Float32Array;
    const col = colAttr.array as Float32Array;
    let o = 0;
    for (const cell of world.cells) {
      const c = getCellColor(cell, viewMode, world.params.seaLevel, factionColors);
      const hMult = 1 + (cell.height * 0.05);
      const cx = cell.center.x * hMult; const cy = cell.center.y * hMult; const cz = cell.center.z * hMult;
      for (let i = 0; i < cell.vertices.length; i++) {
        const v1 = cell.vertices[i]; const v2 = cell.vertices[(i + 1) % cell.vertices.length];
        pos[o] = cx; pos[o + 1] = cy; pos[o + 2] = cz;
        pos[o + 3] = v1.x * hMult; pos[o + 4] = v1.y * hMult; pos[o + 5] = v1.z * hMult;
        pos[o + 6] = v2.x * hMult; pos[o + 7] = v2.y * hMult; pos[o + 8] = v2.z * hMult;
        col[o] = c.r; col[o + 1] = c.g; col[o + 2] = c.b;
        col[o + 3] = c.r; col[o + 4] = c.g; col[o + 5] = c.b;
        col[o + 6] = c.r; col[o + 7] = c.g; col[o + 8] = c.b;
        o += 9;
      }
    }
    posAttr.needsUpdate = true;
    colAttr.needsUpdate = true;
  }, [geometry, world, viewMode, factionColors]);

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
      const hMult = 1 + cell.height * 0.05 + 0.005;
      setBrushCenter([cell.center.x * hMult, cell.center.y * hMult, cell.center.z * hMult]);
      setHighlightCellId(cellId);
  }, [world.cells]);

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

  return (
    <Group>
        <Group ref={spinRef}>
            <Mesh
            ref={meshRef}
            geometry={geometry}
            onPointerMove={tracksPointerMove ? handlePointerMove : undefined}
            onPointerOut={tracksPointerMove ? handlePointerOut : undefined}
            >
                {viewMode === 'political' ? (
                    <MeshBasicMaterial vertexColors toneMapped={false} side={THREE.FrontSide} />
                ) : (
                    <MeshStandardMaterial vertexColors roughness={0.8} metalness={0.1} flatShading side={THREE.FrontSide} />
                )}
                <CityMarkers world={world} viewMode={viewMode} />
                <React.Suspense fallback={null}>
                    <CountryLabels world={world} visible={showFactionOverlay} />
                </React.Suspense>
                <FactionBorders world={world} visible={showFactionOverlay} />
                <RiverLines world={world} visible={showRivers} />
                {showGrid && <LatLongGrid radius={1.06} />}
                {/* Cell highlight outline */}
                {highlightCellId !== null && world.cells[highlightCellId] && (
                    <CellHighlightOutline cell={world.cells[highlightCellId]} />
                )}
                {selectedCellId !== null && (() => {
                    const selectedCell = world.cells[selectedCellId];
                    return selectedCell ? <CellSelectionOverlay cell={selectedCell} /> : null;
                })()}
            </Mesh>
            {/* Brush size ring */}
            {isPaintMode && brushCenter && (
                <BrushRing center={brushCenter} radius={Math.max(0.018, brushSize * 0.038)} />
            )}
            {dymaxionSettings.showOverlay && <DymaxionOverlay settings={dymaxionSettings} />}
            <Mesh scale={[0.99, 0.99, 0.99]}>
                <IcosahedronGeometry args={[1, 16]} />
                <MeshBasicMaterial color="#000000" side={THREE.FrontSide} />
            </Mesh>
        </Group>
    </Group>
  );
};

const WorldViewer: React.FC<{ world: WorldData | null; viewMode: ViewMode; showGrid?: boolean; showRivers?: boolean; showFactionOverlay?: boolean; inspectMode: InspectMode; onInspect: (cellId: number | null) => void; selectedCellId?: number | null; dymaxionSettings: DymaxionSettings; onDymaxionChange: React.Dispatch<React.SetStateAction<DymaxionSettings>>; editMode: EditMode; onPaint: (cellId: number, phase: 'start' | 'stroke' | 'end', isRightClick?: boolean) => void; factionColors?: Map<number, string>; brushSize?: number; }> = ({ world, viewMode, showGrid = false, showRivers = true, showFactionOverlay = true, inspectMode, onInspect, selectedCellId = null, dymaxionSettings, onDymaxionChange, editMode, onPaint, factionColors, brushSize = 1 }) => {
  const [hoveredCell, setHoveredCell] = useState<Cell | null>(null);
  const [paused, setPaused] = useState(false);
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
  }, [overlayMode]);

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
      lon: e.shiftKey ? prev.lon : wrapAngle(prev.lon + dx * sensitivity),
      lat: e.shiftKey ? prev.lat : clampLat(prev.lat + dy * sensitivity),
      roll: e.shiftKey ? wrapAngle(prev.roll + dx * sensitivity) : prev.roll,
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
               showRivers={showRivers}
               showFactionOverlay={showFactionOverlay}
               inspectMode={inspectMode}
               onInspect={onInspect}
               selectedCellId={selectedCellId}
               dymaxionSettings={dymaxionSettings}
               editMode={editMode}
               onPaint={onPaint}
               factionColors={factionColors}
               brushSize={brushSize}
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
      {!world && <div className="absolute inset-0 flex items-center justify-center text-white/50">Forging World...</div>}
      
      <div className="absolute top-4 right-4 z-10 flex gap-2">
         <button
           onClick={() => { setPaused(!paused); }}
           disabled={overlayMode}
           className={`bg-gray-800/80 text-white p-2 backdrop-blur border border-white/10 shadow-lg ${overlayMode ? 'opacity-40 cursor-not-allowed' : 'hover:bg-gray-700'}`}
         >
           {paused ? "▶" : "⏸"}
         </button>
      </div>
    </div>
  );
};

export default WorldViewer;

const LatLongGrid: React.FC<{ radius: number }> = ({ radius }) => {
  const geometry = useMemo(() => {
      const segments = 64; const positions: number[] = [];
      for (let i = 1; i < 18; i++) { 
          const lat = (i * 10 - 90) * (Math.PI / 180); const r = Math.cos(lat) * radius; const y = Math.sin(lat) * radius;
          for (let j = 0; j <= segments; j++) {
              const lon = (j / segments) * Math.PI * 2; const nextLon = ((j + 1) / segments) * Math.PI * 2;
              positions.push(Math.cos(lon) * r, y, Math.sin(lon) * r, Math.cos(nextLon) * r, y, Math.sin(nextLon) * r);
          }
      }
      for (let i = 0; i < 36; i++) { 
          const lon = (i * 10) * (Math.PI / 180); const cosLon = Math.cos(lon); const sinLon = Math.sin(lon);
          for (let j = 0; j <= segments; j++) {
              const lat = (j / segments) * Math.PI - Math.PI/2; const nextLat = ((j + 1) / segments) * Math.PI - Math.PI/2;
              positions.push(Math.cos(lat) * cosLon * radius, Math.sin(lat) * radius, Math.cos(lat) * sinLon * radius, Math.cos(nextLat) * cosLon * radius, Math.sin(nextLat) * radius, Math.cos(nextLat) * sinLon * radius);
          }
      }
      const geo = new THREE.BufferGeometry();
      geo.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
      return geo;
  }, [radius]);
  useEffect(() => () => { geometry.dispose(); }, [geometry]);
  return (
      <LineSegments geometry={geometry}>
          <LineBasicMaterial color="#ffffff" opacity={0.15} transparent depthTest={true} />
      </LineSegments>
  );
};
