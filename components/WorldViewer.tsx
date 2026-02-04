import React, { useMemo, useRef, useEffect, useState, useCallback } from 'react';
import { Canvas, useFrame, useThree } from '@react-three/fiber';
import { OrbitControls, Stars, Text } from '@react-three/drei';
import * as THREE from 'three';
import { WorldData, ViewMode, Cell } from '../types';
import { getCellColor } from '../utils/colors';
import { MousePointer2, EyeOff, ChevronDown, ChevronUp } from 'lucide-react';

const Mesh = 'mesh' as any;
const Group = 'group' as any;
const AmbientLight = 'ambientLight' as any;
const PointLight = 'pointLight' as any;
const DirectionalLight = 'directionalLight' as any;
const MeshStandardMaterial = 'meshStandardMaterial' as any;
const InstancedMesh = 'instancedMesh' as any;

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
                    <cylinderGeometry args={[0.008, 0.008, 0.08, 6]} />
                    <meshBasicMaterial color="#ef4444" toneMapped={false} />
                </InstancedMesh>
            )}
            {towns.length > 0 && (
                <InstancedMesh ref={townsRef} args={[undefined, undefined, towns.length]} visible={viewMode === 'political'}>
                    <cylinderGeometry args={[0.005, 0.005, 0.04, 5]} />
                    <meshBasicMaterial color="#ffffff" toneMapped={false} />
                </InstancedMesh>
            )}
        </>
    );
};

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

const WorldMesh: React.FC<{ 
  world: WorldData, 
  viewMode: ViewMode, 
  onHover: (cell: Cell | null) => void, 
  paused: boolean, 
  showGrid: boolean,
  hudActive: boolean 
}> = ({ world, viewMode, onHover, paused, showGrid, hudActive }) => {
  const meshRef = useRef<THREE.Mesh>(null);
  const coreRef = useRef<THREE.Mesh>(null);
  const lastUpdate = useRef<number>(0);
  
  useFrame((state, delta) => {
    if (!paused) {
        if (meshRef.current) meshRef.current.rotation.y += delta * 0.05;
        if (coreRef.current) coreRef.current.rotation.y += delta * 0.05;
    }
  });

  const geometry = useMemo(() => {
    const positions: number[] = []; const colors: number[] = [];
    world.cells.forEach(cell => {
      const c = getCellColor(cell, viewMode, world.params.seaLevel);
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

  const handlePointerMove = useCallback((e: any) => {
      if (!hudActive) return;
      
      // Throttle lookup to every 100ms for smoother manipulation performance
      const now = Date.now();
      if (now - lastUpdate.current < 80) return; 
      lastUpdate.current = now;

      if (e.faceIndex !== undefined) {
          const cellId = faceMap[e.faceIndex];
          if (cellId !== undefined) onHover(world.cells[cellId]);
      } else { onHover(null); }
  }, [hudActive, faceMap, world.cells, onHover]);

  return (
    <Group>
        <Mesh 
          ref={meshRef} 
          geometry={geometry} 
          onPointerMove={hudActive ? handlePointerMove : undefined} 
          onPointerOut={hudActive ? () => onHover(null) : undefined}
        >
          <MeshStandardMaterial vertexColors roughness={0.8} metalness={0.1} flatShading side={THREE.FrontSide} />
          <CityMarkers world={world} viewMode={viewMode} />
          <CountryLabels world={world} viewMode={viewMode} />
          {showGrid && <LatLongGrid radius={1.06} />}
        </Mesh>
        <Mesh ref={coreRef} scale={[0.99, 0.99, 0.99]}>
            <icosahedronGeometry args={[1, 16]} />
            <meshBasicMaterial color="#000000" side={THREE.FrontSide} />
        </Mesh>
    </Group>
  );
};

const WorldViewer: React.FC<{ world: WorldData | null; viewMode: ViewMode; showGrid?: boolean }> = ({ world, viewMode, showGrid = false }) => {
  const [hoveredCell, setHoveredCell] = useState<Cell | null>(null);
  const [paused, setPaused] = useState(false);
  const [hudActive, setHudActive] = useState(true);
  const [hudCollapsed, setHudCollapsed] = useState(false);

  const factionMap = useMemo(() => {
    if (!world?.civData) return new Map();
    return new Map(world.civData.factions.map(f => [f.id, f]));
  }, [world?.civData]);

  const locationName = useMemo(() => {
      if (!hoveredCell || !world?.civData || !hudActive) return null;
      const { regionId, provinceId } = hoveredCell;
      if (regionId === undefined) return "Unclaimed Land";
      
      const faction = factionMap.get(regionId);
      if (!faction) return `Region ${regionId}`;
      const province = provinceId !== undefined && faction.provinces[provinceId];
      return (
          <div className="flex flex-col">
            <span className="font-bold text-blue-200">{faction.name}</span>
            {province && <span className="text-xs text-blue-100/70">{province.name}</span>}
          </div>
      );
  }, [hoveredCell, factionMap, hudActive]);

  return (
    <div className="w-full h-full bg-black relative group">
      <Canvas camera={{ position: [0, 0, 2.5], fov: 45 }}>
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
               hudActive={hudActive} 
             />
          </Group>
        )}
        <OrbitControls enablePan={false} minDistance={1.2} maxDistance={6} />
      </Canvas>
      {!world && <div className="absolute inset-0 flex items-center justify-center text-white/50">Forging World...</div>}
      
      {/* HUD - TOP CENTER */}
      <div className="absolute top-6 left-1/2 -translate-x-1/2 flex flex-col items-center gap-2 pointer-events-none z-10">
          <div className={`bg-black/80 backdrop-blur text-white rounded shadow-xl border border-white/20 transition-all duration-300 pointer-events-auto ${hudCollapsed ? 'w-10 overflow-hidden' : 'min-w-[200px]'}`}>
              <div className="flex items-center justify-between p-2 border-b border-white/10">
                  {!hudCollapsed && <span className="text-[10px] font-bold text-gray-400 uppercase tracking-tighter">Inspector</span>}
                  <div className="flex items-center gap-2">
                      <button 
                        onClick={() => setHudActive(!hudActive)}
                        className={`p-1 rounded transition-colors ${hudActive ? 'text-blue-400 hover:bg-blue-900/40' : 'text-gray-600 hover:bg-gray-800'}`}
                        title={hudActive ? "Disable Hover Info (Performance Boost)" : "Enable Hover Info"}
                      >
                        {hudActive ? <MousePointer2 size={12}/> : <EyeOff size={12}/>}
                      </button>
                      <button 
                        onClick={() => setHudCollapsed(!hudCollapsed)}
                        className="p-1 rounded text-gray-400 hover:bg-gray-800"
                      >
                        {hudCollapsed ? <ChevronDown size={12}/> : <ChevronUp size={12}/>}
                      </button>
                  </div>
              </div>

              {!hudCollapsed && hudActive && hoveredCell && (
                <div className="p-2 text-xs">
                    <div className="font-bold flex justify-between gap-4 mb-2 border-b border-white/10 pb-1 items-start">
                        <div className="flex flex-col">
                          <span>Cell {hoveredCell.id}</span>
                          {locationName && <div className="mt-1">{locationName}</div>}
                        </div>
                        <span style={{color: '#aaa'}}>{hoveredCell.biome}</span>
                    </div>
                    <div className="grid grid-cols-2 gap-x-4 gap-y-1">
                        <div className="text-gray-400">Temp: <span className="text-white">{hoveredCell.temperature.toFixed(1)}°C</span></div>
                        <div className="text-gray-400">Rain: <span className="text-white">{(hoveredCell.moisture*100).toFixed(0)}%</span></div>
                        <div className="text-gray-400">Elev: <span className="text-white">{(hoveredCell.height*100).toFixed(0)}%</span></div>
                        <div className="text-gray-400">Pop: <span className="text-white">{hoveredCell.population?.toLocaleString()}</span></div>
                    </div>
                </div>
              )}
              {!hudCollapsed && !hudActive && (
                <div className="p-4 text-[10px] text-gray-500 text-center italic">
                  Hover Info Disabled
                </div>
              )}
              {!hudCollapsed && hudActive && !hoveredCell && (
                <div className="p-4 text-[10px] text-gray-500 text-center italic">
                  Hover over a cell...
                </div>
              )}
          </div>
      </div>
      
      <div className="absolute top-4 right-4 z-10 flex gap-2">
         <button onClick={() => setPaused(!paused)} className="bg-gray-800/80 text-white p-2 rounded hover:bg-gray-700 backdrop-blur border border-white/10 shadow-lg">{paused ? "▶" : "⏸"}</button>
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
      <lineSegments geometry={geometry}>
          <lineBasicMaterial color="#ffffff" opacity={0.15} transparent depthTest={true} />
      </lineSegments>
  );
};
