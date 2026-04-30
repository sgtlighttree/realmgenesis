import * as THREE from 'three';
import { GLTFExporter } from 'three/examples/jsm/exporters/GLTFExporter.js';
import { WorldData, ViewMode, Cell } from '../types';
import { getCellColor } from './colors';

const buildWorldMesh = (world: WorldData, viewMode: ViewMode): THREE.Mesh => {
  const positions: number[] = [];
  const colors: number[] = [];

  world.cells.forEach(cell => {
    const c = getCellColor(cell, viewMode, world.params.seaLevel);
    const hMult = 1 + cell.height * 0.05;
    const cx = cell.center.x * hMult;
    const cy = cell.center.y * hMult;
    const cz = cell.center.z * hMult;
    for (let i = 0; i < cell.vertices.length; i++) {
      const v1 = cell.vertices[i];
      const v2 = cell.vertices[(i + 1) % cell.vertices.length];
      positions.push(cx, cy, cz, v1.x * hMult, v1.y * hMult, v1.z * hMult, v2.x * hMult, v2.y * hMult, v2.z * hMult);
      colors.push(c.r, c.g, c.b, c.r, c.g, c.b, c.r, c.g, c.b);
    }
  });

  const geo = new THREE.BufferGeometry();
  geo.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
  geo.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3));
  geo.computeVertexNormals();

  const mat = new THREE.MeshStandardMaterial({ vertexColors: true, roughness: 0.8, metalness: 0.1 });
  const mesh = new THREE.Mesh(geo, mat);
  mesh.name = 'World';
  return mesh;
};

const buildRivers = (world: WorldData): THREE.LineSegments | null => {
  if (!world.rivers || world.rivers.length === 0) return null;
  const positions: number[] = [];
  world.rivers.forEach(path => {
    for (let i = 0; i < path.length - 1; i++) {
      positions.push(path[i].x, path[i].y, path[i].z);
      positions.push(path[i + 1].x, path[i + 1].y, path[i + 1].z);
    }
  });
  if (positions.length === 0) return null;
  const geo = new THREE.BufferGeometry();
  geo.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
  const lines = new THREE.LineSegments(geo, new THREE.LineBasicMaterial({ color: 0x38bdf8 }));
  lines.name = 'Rivers';
  return lines;
};

const buildMarkerMesh = (
  cells: Cell[],
  rTop: number, rBot: number, height: number, segments: number,
  color: number, name: string,
): THREE.Mesh | null => {
  if (cells.length === 0) return null;

  // Expand indexed cylinder template to non-indexed so each instance
  // can be baked in with its own transform without index offsets.
  const template = new THREE.CylinderGeometry(rTop, rBot, height, segments);
  const nonIndexed = template.toNonIndexed();
  template.dispose();

  const posAttr = nonIndexed.getAttribute('position') as THREE.BufferAttribute;
  const normAttr = nonIndexed.getAttribute('normal') as THREE.BufferAttribute;
  const vertsPerMarker = posAttr.count;

  const allPos  = new Float32Array(cells.length * vertsPerMarker * 3);
  const allNorm = new Float32Array(cells.length * vertsPerMarker * 3);

  const yAxis = new THREE.Vector3(0, 1, 0);
  const quat  = new THREE.Quaternion();
  const m4    = new THREE.Matrix4();
  const m3    = new THREE.Matrix3();
  const vp    = new THREE.Vector3();
  const vn    = new THREE.Vector3();

  cells.forEach((cell, ci) => {
    const hMult = 1 + cell.height * 0.05;
    const surfPos = new THREE.Vector3(
      cell.center.x * hMult, cell.center.y * hMult, cell.center.z * hMult
    );
    const up = surfPos.clone().normalize();
    quat.setFromUnitVectors(yAxis, up);
    // Place base on sphere surface, extend outward by full height
    const center = surfPos.clone().add(up.multiplyScalar(height * 0.5));
    m4.compose(center, quat, new THREE.Vector3(1, 1, 1));
    m3.getNormalMatrix(m4);

    const base = ci * vertsPerMarker * 3;
    for (let v = 0; v < vertsPerMarker; v++) {
      vp.set(posAttr.getX(v), posAttr.getY(v), posAttr.getZ(v)).applyMatrix4(m4);
      vn.set(normAttr.getX(v), normAttr.getY(v), normAttr.getZ(v)).applyMatrix3(m3).normalize();
      allPos[base + v * 3]     = vp.x;
      allPos[base + v * 3 + 1] = vp.y;
      allPos[base + v * 3 + 2] = vp.z;
      allNorm[base + v * 3]     = vn.x;
      allNorm[base + v * 3 + 1] = vn.y;
      allNorm[base + v * 3 + 2] = vn.z;
    }
  });

  nonIndexed.dispose();

  const geo = new THREE.BufferGeometry();
  geo.setAttribute('position', new THREE.BufferAttribute(allPos, 3));
  geo.setAttribute('normal',   new THREE.BufferAttribute(allNorm, 3));
  const mesh = new THREE.Mesh(geo, new THREE.MeshBasicMaterial({ color }));
  mesh.name = name;
  return mesh;
};

export const exportGLB = (world: WorldData, viewMode: ViewMode) => {
  const scene = new THREE.Scene();

  scene.add(buildWorldMesh(world, viewMode));

  const rivers = buildRivers(world);
  if (rivers) scene.add(rivers);

  if (world.civData) {
    const capitals: Cell[] = [];
    const towns: Cell[] = [];
    world.cells.forEach(cell => {
      if (cell.isCapital) capitals.push(cell);
      else if (cell.isTown) towns.push(cell);
    });
    // Capitals: 6-sided, height 0.08  |  Towns: 5-sided, height 0.04
    const capMesh  = buildMarkerMesh(capitals, 0.008, 0.008, 0.08, 6, 0xef4444, 'Capitals');
    const townMesh = buildMarkerMesh(towns,    0.005, 0.005, 0.04, 5, 0xffffff, 'Towns');
    if (capMesh)  scene.add(capMesh);
    if (townMesh) scene.add(townMesh);
  }

  const exporter = new GLTFExporter();
  const mapName = world.params.mapName || 'map';
  exporter.parse(
    scene,
    (result) => {
      const blob = new Blob([result as ArrayBuffer], { type: 'model/gltf-binary' });
      const url = URL.createObjectURL(blob);
      const link = document.createElement('a');
      link.download = `realmgenesis_${mapName}_${world.params.seed}_${viewMode}.glb`;
      link.href = url;
      link.click();
      URL.revokeObjectURL(url);
    },
    (error) => { console.error('GLB export error:', error); },
    { binary: true },
  );
};
