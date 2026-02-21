import React, { useEffect, useMemo, useRef } from 'react';
import * as d3 from 'd3';
import { WorldData, ViewMode, DymaxionSettings } from '../types';
import { getCellColor } from '../utils/colors';

const PHI = (1 + Math.sqrt(5)) / 2;
const RAW_VERTS = [
  [-1, PHI, 0],
  [1, PHI, 0],
  [-1, -PHI, 0],
  [1, -PHI, 0],
  [0, -1, PHI],
  [0, 1, PHI],
  [0, -1, -PHI],
  [0, 1, -PHI],
  [PHI, 0, -1],
  [PHI, 0, 1],
  [-PHI, 0, -1],
  [-PHI, 0, 1],
];

const FACES = [
  [0, 11, 5],
  [0, 5, 1],
  [0, 1, 7],
  [0, 7, 10],
  [0, 10, 11],
  [1, 5, 9],
  [5, 11, 4],
  [11, 10, 2],
  [10, 7, 6],
  [7, 1, 8],
  [3, 9, 4],
  [3, 4, 2],
  [3, 2, 6],
  [3, 6, 8],
  [3, 8, 9],
  [4, 9, 5],
  [2, 4, 11],
  [6, 2, 10],
  [8, 6, 7],
  [9, 8, 1],
];

const normalize = (v: number[]) => {
  const len = Math.hypot(v[0], v[1], v[2]) || 1;
  return [v[0] / len, v[1] / len, v[2] / len];
};

const cartesianToLonLat = (v: number[]) => {
  const lon = Math.atan2(v[2], v[0]) * (180 / Math.PI);
  const lat = Math.asin(Math.max(-1, Math.min(1, v[1]))) * (180 / Math.PI);
  return [lon, lat] as [number, number];
};

const buildEdges = () => {
  const verts = RAW_VERTS.map(normalize);
  const edgeSet = new Set<string>();
  const edges: [number, number][] = [];
  FACES.forEach((face) => {
    for (let i = 0; i < 3; i++) {
      // eslint-disable-next-line security/detect-object-injection, nosemgrep
        // codacy-disable-next-line
      const a = face[i];
      const b = face[(i + 1) % 3];
      const key = a < b ? `${a}-${b}` : `${b}-${a}`;
      if (!edgeSet.has(key)) {
        edgeSet.add(key);
        edges.push([a, b]);
      }
    }
  });
  const lonLatVerts = verts.map(cartesianToLonLat);
  // eslint-disable-next-line security/detect-object-injection, nosemgrep
        // codacy-disable-next-line
  return edges.map(([a, b]) => [lonLatVerts[a], lonLatVerts[b]] as [[number, number], [number, number]]);
};

type Props = {
  world: WorldData | null;
  viewMode: ViewMode;
  settings: DymaxionSettings;
  onChange?: React.Dispatch<React.SetStateAction<DymaxionSettings>>;
  width?: number;
  height?: number;
};

const clamp = (v: number, min: number, max: number) => Math.max(min, Math.min(max, v));
const wrapAngle = (v: number) => {
  let out = v;
  while (out > 180) out -= 360;
  while (out < -180) out += 360;
  return out;
};

const DymaxionPreview2D: React.FC<Props> = ({ world, viewMode, settings, onChange, width = 240, height = 120 }) => {
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const edges = useMemo(buildEdges, []);
  const dragging = useRef(false);
  const lastPos = useRef({ x: 0, y: 0 });

  useEffect(() => {
    if (!world) return;
    const canvas = canvasRef.current;
    if (!canvas) return;
    const dpr = Math.min(2, window.devicePixelRatio || 1);
    canvas.width = Math.floor(width * dpr);
    canvas.height = Math.floor(height * dpr);
      canvas.style.width = '100%';
      canvas.style.height = 'auto';
    const ctx = canvas.getContext('2d');
    if (!ctx) return;

    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);

    const projection = d3.geoMercator().fitSize([width, height], { type: 'Sphere' } as any);
    const path = d3.geoPath(projection, ctx);

    ctx.fillStyle = viewMode === 'satellite' || viewMode === 'biome' ? '#050505' : '#000000';
    ctx.fillRect(0, 0, width, height);

    ctx.save();
    ctx.translate(width, 0);
    ctx.scale(-1, 1);

    world.cells.forEach((cell, i) => {
      // eslint-disable-next-line security/detect-object-injection, nosemgrep
        // codacy-disable-next-line
      const feature = world.geoJson?.features?.[i];
      if (!feature) return;
      const color = getCellColor(cell, viewMode, world.params.seaLevel);
      ctx.beginPath();
      path(feature);
      ctx.fillStyle = '#' + color.getHexString();
      ctx.fill();
    });

    ctx.restore();

    const rotate = d3.geoRotation([settings.lon, settings.lat, settings.roll]);
    ctx.save();
    ctx.strokeStyle = 'rgba(255,255,255,0.9)';
    ctx.lineWidth = 1;

    edges.forEach(([a, b]) => {
      const ia = rotate(a);
      const ib = rotate(b);
      const interp = d3.geoInterpolate(ia, ib);
      ctx.beginPath();
      for (let i = 0; i <= 30; i++) {
        const p = interp(i / 30);
        const projected = projection(p);
        if (!projected) continue;
        const x = width - projected[0];
        const y = projected[1];
        if (i === 0) ctx.moveTo(x, y);
        else ctx.lineTo(x, y);
      }
      ctx.stroke();
    });
    ctx.restore();
  }, [world, viewMode, settings, edges, width, height]);

  const onPointerDown = (event: React.PointerEvent<HTMLCanvasElement>) => {
    dragging.current = true;
    lastPos.current = { x: event.clientX, y: event.clientY };
    event.currentTarget.setPointerCapture(event.pointerId);
  };

  const onPointerMove = (event: React.PointerEvent<HTMLCanvasElement>) => {
    if (!dragging.current || !onChange) return;
    const dx = event.clientX - lastPos.current.x;
    const dy = event.clientY - lastPos.current.y;
    lastPos.current = { x: event.clientX, y: event.clientY };
    onChange((prev) => {
      if (event.shiftKey) {
        return { ...prev, roll: wrapAngle(prev.roll + dx * 0.3) };
      }
      return {
        ...prev,
        lon: wrapAngle(prev.lon + dx * 0.3),
        lat: clamp(prev.lat - dy * 0.3, -90, 90)
      };
    });
  };

  const onPointerUp = (event: React.PointerEvent<HTMLCanvasElement>) => {
    dragging.current = false;
    event.currentTarget.releasePointerCapture(event.pointerId);
  };

  return (
    <div className="rounded border border-gray-800 bg-black/60 p-2">
      <div className="text-[10px] text-gray-400 mb-2">2D Mercator Preview (drag to rotate, Shift+drag to roll)</div>
      <canvas
        ref={canvasRef}
        onPointerDown={onPointerDown}
        onPointerMove={onPointerMove}
        onPointerUp={onPointerUp}
        className="block w-full rounded border border-gray-800"
      />
    </div>
  );
};

export default DymaxionPreview2D;
