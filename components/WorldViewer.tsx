import React, { useMemo, useRef, useEffect, useState, useCallback } from 'react';
import { Canvas, useFrame } from '@react-three/fiber';
import { OrbitControls, Stars, Text, Line } from '@react-three/drei';
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
        const updateMesh = (mesh: THREE.InstancedMesh, data: Cell[], heightBase: number) => {
            if (!mesh) return;
            data.forEach((cell, i) => {
                const h = 1 + (cell.height * 0.05);
                const pos = new THREE.Vector3(cell.center.x * h, cell.center.y * h, cell.center.z * h);
                const offsetPos = pos.clone().add(pos.clone().normalize().multiplyScalar(heightBase * 0.5));
                dummy.position.copy(offsetPos);
                dummy.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), pos.normalize());
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
    const geometry = useMemo(() => {
        if (!world.rivers || !visible) return null;
        
        const positions: number[] = [];
        
        // Batch all river segments into a single LineSegments geometry for performance
        // Rendering thousands of individual <Line> components causes massive overhead/freezes
        world.rivers.forEach(path => {
            if (path.length < 2) return;
            
            // Create Curve for smoothing
            const vectors = path.map(p => new THREE.Vector3(p.x, p.y, p.z));
            const curve = new THREE.CatmullRomCurve3(vectors);
            
            // Adaptive sampling based on length, but simple count is safer for perf
            const points = curve.getPoints(Math.min(50, vectors.length * 4));
            
            for (let i = 0; i < points.length - 1; i++) {
                // eslint-disable-next-line @typescript-eslint/no-explicit-any
        // codacy-disable-next-line
                positions.push(points[i].x, points[i].y, points[i].z);
                positions.push(points[i+1].x, points[i+1].y, points[i+1].z);
            }
        });

        if (positions.length === 0) return null;

        const geo = new THREE.BufferGeometry();
        geo.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
        return geo;
    }, [world, visible]);

    if (!visible || !geometry) return null;

    return (
        <LineSegments geometry={geometry}>
            <LineBasicMaterial color="#38bdf8" opacity={0.8} transparent linewidth={1.5} />
        </LineSegments>
    );
}

const CountryLabels: React.FC<{ world: WorldData; viewMode: ViewMode }> = ({ world, viewMode }) => {
    const labels = useMemo(() => {
        if (!world.civData || viewMode !== 'political') return [];
        return world.civData.factions.map(f => {
            let sumX = 0, sumY = 0, sumZ = 0, count = 0;
            for (const cell of world.cells) {
                if (cell.regionId === f.id) {
                    sumX += cell.center.x; sumY += cell.center.y; sumZ += cell.center.z;
                    count++;
                }
            }
            if (count === 0) return null;
            const avg = new THREE.Vector3(sumX/count, sumY/count, sumZ/count).normalize().multiplyScalar(1.08);
            return { id: f.id, name: f.name, position: avg };
        }).filter(Boolean) as {id: number, name: string, position: THREE.Vector3}[];
    }, [world, viewMode]);

    if (viewMode !== 'political') return null;

    return (
        <Group>
            {labels.map(l => (
                <Text
                    key={l.id}
                    position={l.position}
                    color="white"
                    anchorX="center" anchorY="middle"
                    fontSize={0.07}
                    outlineWidth={0.008} outlineColor="#000000"
                    quaternion={new THREE.Quaternion().setFromUnitVectors(new THREE.Vector3(0,0,1), l.position.clone().normalize())}
                >
                    {l.name}
                </Text>
            ))}
        </Group>
    );
};

const FactionBorders: React.FC<{ world: WorldData; viewMode: ViewMode }> = ({ world, viewMode }) => {
  const geometry = useMemo(() => {
      // Only show for political view
      if (viewMode !== 'political') return null;

      const positions: number[] = [];
      const threshold = 0.000001; 

      // Iterate unique pairs of neighbors to find borders
      world.cells.forEach(cellA => {
          cellA.neighbors.forEach(nId => {
              // eslint-disable-next-line @typescript-eslint/no-explicit-any
        // codacy-disable-next-line
              const cellB = world.cells[nId];
              if (cellA.id >= cellB.id) return; // Process pair once
              
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
  }, [world, viewMode]);

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
  return (
    <LineSegments geometry={geometry} renderOrder={10}>
      <LineBasicMaterial color="#ffffff" opacity={0.8} transparent depthTest={false} />
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
}> = ({ world, viewMode, onHover, paused, showGrid, showRivers, inspectMode, onInspect, dymaxionSettings, editMode, onPaint, factionColors, brushSize }) => {
  const spinRef = useRef<THREE.Group>(null);
  const lastUpdate = useRef<number>(0);
  const lastPaintedCell = useRef<number | null>(null);
  const [brushCenter, setBrushCenter] = useState<[number, number, number] | null>(null);
  const [highlightCellId, setHighlightCellId] = useState<number | null>(null);
  
  useFrame((state, delta) => {
    if (!paused) {
        if (spinRef.current) spinRef.current.rotation.y += delta * 0.05;
    }
  });

  const geometry = useMemo(() => {
    const positions: number[] = []; const colors: number[] = [];
    world.cells.forEach(cell => {
      const c = getCellColor(cell, viewMode, world.params.seaLevel, factionColors);
      const hMult = 1 + (cell.height * 0.05); 
      const cx = cell.center.x * hMult; const cy = cell.center.y * hMult; const cz = cell.center.z * hMult;
      for (let i = 0; i < cell.vertices.length; i++) {
        const next = (i + 1) % cell.vertices.length; const v1 = cell.vertices[i]; const v2 = cell.vertices[next];
        positions.push(cx, cy, cz, v1.x * hMult, v1.y * hMult, v1.z * hMult, v2.x * hMult, v2.y * hMult, v2.z * hMult);
        colors.push(c.r, c.g, c.b, c.r, c.g, c.b, c.r, c.g, c.b);
      }
    });
    const geo = new THREE.BufferGeometry();
    geo.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
    geo.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3));
    geo.computeVertexNormals();
    return geo;
  }, [world, viewMode]);

  const faceMap = useMemo(() => {
     const map: number[] = []; 
     world.cells.forEach(cell => { 
       for(let i=0; i<cell.vertices.length; i++) map.push(cell.id); 
     }); 
     return map;
  }, [world]);

  const getTriangleIndex = useCallback((e: any) => {
      if (e.faceIndex !== undefined && e.faceIndex !== null) return e.faceIndex;
      if (e.face && e.face.a !== undefined && e.face.a !== null) return Math.floor(e.face.a / 3);
      return null;
  }, []);

  const isPaintMode = editMode !== 'off' && editMode !== 'world-edit';

  const handlePointerMove = useCallback((e: any) => {
      if (dymaxionSettings.mode === 'overlay') return;
      const now = Date.now();
      if (now - lastUpdate.current < 50) return;
      lastUpdate.current = now;

      const triIndex = getTriangleIndex(e);
      if (triIndex === null) { if (!isPaintMode) setBrushCenter(null); return; }
      const cellId = faceMap[triIndex];
      if (cellId === undefined) return;

      const cell = world.cells[cellId];
      const hMult = 1 + cell.height * 0.05 + 0.005;
      setBrushCenter([cell.center.x * hMult, cell.center.y * hMult, cell.center.z * hMult]);
      setHighlightCellId(cellId);

      if (isPaintMode) {
          lastPaintedCell.current = cellId;
          onPaint(cellId, 'stroke');
          return;
      }
      if (inspectMode === 'hover') onHover(cell);
  }, [isPaintMode, inspectMode, faceMap, world.cells, onHover, onPaint, getTriangleIndex, dymaxionSettings.mode]);

  const handlePointerDown = useCallback((e: any) => {
      if (dymaxionSettings.mode === 'overlay') return;
      const triIndex = getTriangleIndex(e);
      if (triIndex === null) return;
      const cellId = faceMap[triIndex];
      if (cellId === undefined) return;

      if (isPaintMode) {
          lastPaintedCell.current = cellId;
          const isRight = e.button === 2 || e.buttons === 2;
          onPaint(cellId, 'start', isRight);
          if (!isRight) onPaint(cellId, 'stroke');
          return;
      }
      if (inspectMode === 'click') onInspect(cellId);
  }, [isPaintMode, inspectMode, faceMap, onInspect, onPaint, getTriangleIndex, dymaxionSettings.mode]);

  const handlePointerUp = useCallback(() => {
      if (isPaintMode && lastPaintedCell.current !== null) {
          onPaint(lastPaintedCell.current, 'end');
          lastPaintedCell.current = null;
      }
  }, [isPaintMode, onPaint]);

  const handleClick = useCallback((e: any) => {
      if (dymaxionSettings.mode === 'overlay') return;
      if (isPaintMode) return; // paint handled by pointerDown/Up
      if (inspectMode !== 'click') return;
      const triIndex = getTriangleIndex(e);
      if (triIndex !== null) {
          const cellId = faceMap[triIndex];
          onInspect(cellId !== undefined ? cellId : null);
      }
  }, [isPaintMode, inspectMode, faceMap, onInspect, getTriangleIndex, dymaxionSettings.mode]);

  return (
    <Group>
        <Group ref={spinRef}>
            <Mesh
            geometry={geometry}
            onPointerMove={handlePointerMove}
            onPointerOut={() => { setBrushCenter(null); setHighlightCellId(null); if (inspectMode === 'hover' && !isPaintMode) onHover(null); }}
            onPointerDown={handlePointerDown}
            onPointerUp={isPaintMode ? handlePointerUp : undefined}
            onClick={handleClick}
            onContextMenu={(e: any) => { e.nativeEvent?.preventDefault?.(); }}
            >
                <MeshStandardMaterial vertexColors roughness={0.8} metalness={0.1} flatShading side={THREE.FrontSide} />
                <CityMarkers world={world} viewMode={viewMode} />
                <React.Suspense fallback={null}>
                    <CountryLabels world={world} viewMode={viewMode} />
                </React.Suspense>
                <FactionBorders world={world} viewMode={viewMode} />
                <RiverLines world={world} visible={showRivers} />
                {showGrid && <LatLongGrid radius={1.06} />}
                {/* Cell highlight outline */}
                {highlightCellId !== null && (() => {
                    const hCell = world.cells[highlightCellId];
                    if (!hCell) return null;
                    const hm = 1 + hCell.height * 0.05 + 0.004;
                    const verts = hCell.vertices;
                    const pts: number[] = [];
                    for (let i = 0; i < verts.length; i++) {
                        const a = verts[i], b = verts[(i + 1) % verts.length];
                        pts.push(a.x*hm, a.y*hm, a.z*hm, b.x*hm, b.y*hm, b.z*hm);
                    }
                    const geo = new THREE.BufferGeometry();
                    geo.setAttribute('position', new THREE.Float32BufferAttribute(pts, 3));
                    return (
                        <LineSegments geometry={geo} renderOrder={9}>
                            <LineBasicMaterial color="#ffff00" opacity={0.9} transparent depthTest={false} />
                        </LineSegments>
                    );
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

const WorldViewer: React.FC<{ world: WorldData | null; viewMode: ViewMode; showGrid?: boolean; showRivers?: boolean; inspectMode: InspectMode; onInspect: (cellId: number | null) => void; dymaxionSettings: DymaxionSettings; onDymaxionChange: React.Dispatch<React.SetStateAction<DymaxionSettings>>; editMode: EditMode; onPaint: (cellId: number, phase: 'start' | 'stroke' | 'end', isRightClick?: boolean) => void; factionColors?: Map<number, string>; brushSize?: number; }> = ({ world, viewMode, showGrid = false, showRivers = true, inspectMode, onInspect, dymaxionSettings, onDymaxionChange, editMode, onPaint, factionColors, brushSize = 1 }) => {
  const [hoveredCell, setHoveredCell] = useState<Cell | null>(null);
  const [paused, setPaused] = useState(false);
  const [isSpaceHeld, setIsSpaceHeld] = useState(false);
  const dragRef = useRef<{ active: boolean; lastX: number; lastY: number }>({ active: false, lastX: 0, lastY: 0 });
  const overlayMode = dymaxionSettings.mode === 'overlay';
  const isPaintModeActive = editMode !== 'off' && editMode !== 'world-edit';

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
        onPointerMissed={() => {
          if (inspectMode === 'click') onInspect(null);
        }}
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
               inspectMode={inspectMode}
               onInspect={onInspect}
               dymaxionSettings={dymaxionSettings}
               editMode={editMode}
               onPaint={onPaint}
               factionColors={factionColors}
               brushSize={brushSize}
             />
          </Group>
        )}
        <OrbitControls enablePan={false} minDistance={1.2} maxDistance={6} enableRotate={!overlayMode && (!isPaintModeActive || isSpaceHeld)} />
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
  return (
      <LineSegments geometry={geometry}>
          <LineBasicMaterial color="#ffffff" opacity={0.15} transparent depthTest={true} />
      </LineSegments>
  );
};
