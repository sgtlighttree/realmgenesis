import React, { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import * as d3 from 'd3';
import { WorldData, ViewMode, InspectMode, DymaxionSettings, EditMode } from '../types';
import { getCellColor } from '../utils/colors';
import { buildDymaxionNet } from '../utils/dymaxion';
import { insideTri, barycentric, normalizeVec, toLonLat, lonLatToPoint3, Point2, Point3 } from '../utils/geo';

type Size = { width: number; height: number };
type FactionLabel = { id: number; name: string; position: Point3 };
type FactionOverlayData = { labels: FactionLabel[]; borders: Array<[Point3, Point3]> };

const clamp = (v: number, min: number, max: number) => Math.max(min, Math.min(max, v));
const INTERACTION_DPR = 1;
const MAX_SHARP_DPR = 3;
const MAX_SHARP_SCALE = 2.5;
const SETTLE_MS = 200;

const getDymaxionNetTransform = (layout: DymaxionSettings['layout'], canvasWidth: number, canvasHeight: number) => {
  const net = buildDymaxionNet(layout);
  const faces = net.faces;
  let minX = Infinity;
  let minY = Infinity;
  let maxX = -Infinity;
  let maxY = -Infinity;
  faces.forEach((face) => {
    face.vertices.forEach((v) => {
      minX = Math.min(minX, v[0]);
      minY = Math.min(minY, v[1]);
      maxX = Math.max(maxX, v[0]);
      maxY = Math.max(maxY, v[1]);
    });
  });

  const pad = 8;
  const netWidth = Math.max(1e-6, maxX - minX);
  const netHeight = Math.max(1e-6, maxY - minY);
  const scale = Math.min((canvasWidth - pad * 2) / netWidth, (canvasHeight - pad * 2) / netHeight);
  const offsetX = (canvasWidth - netWidth * scale) / 2 - minX * scale;
  const offsetY = (canvasHeight - netHeight * scale) / 2 - minY * scale;

  return { faces, scale, offsetX, offsetY };
};

const dot3 = (a: Point3, b: Point3) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const sub3 = (a: Point3, b: Point3): Point3 => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const cross3 = (a: Point3, b: Point3): Point3 => [
  a[1] * b[2] - a[2] * b[1],
  a[2] * b[0] - a[0] * b[2],
  a[0] * b[1] - a[1] * b[0],
];

const barycentric3D = (p: Point3, a: Point3, b: Point3, c: Point3) => {
  const v0 = sub3(b, a);
  const v1 = sub3(c, a);
  const v2 = sub3(p, a);
  const d00 = dot3(v0, v0);
  const d01 = dot3(v0, v1);
  const d11 = dot3(v1, v1);
  const d20 = dot3(v2, v0);
  const d21 = dot3(v2, v1);
  const denom = d00 * d11 - d01 * d01;
  if (!denom) return null;
  const v = (d11 * d20 - d01 * d21) / denom;
  const w = (d00 * d21 - d01 * d20) / denom;
  return [1 - v - w, v, w] as [number, number, number];
};

const pointInsideSphericalFace = (p: Point3, vertices: Point3[]) => {
  const centroid = normalizeVec([
    vertices[0][0] + vertices[1][0] + vertices[2][0],
    vertices[0][1] + vertices[1][1] + vertices[2][1],
    vertices[0][2] + vertices[1][2] + vertices[2][2],
  ]);
  for (let i = 0; i < 3; i++) {
    const a = vertices[i];
    const b = vertices[(i + 1) % 3];
    let edgeNormal = normalizeVec(cross3(a, b));
    if (dot3(edgeNormal, centroid) < 0) {
      edgeNormal = [-edgeNormal[0], -edgeNormal[1], -edgeNormal[2]];
    }
    if (dot3(edgeNormal, p) < -1e-7) return false;
  }
  return true;
};

const projectDymaxionPoint = (
  position: Point3,
  layout: DymaxionSettings['layout'],
  canvasWidth: number,
  canvasHeight: number,
  lon: number,
  lat: number,
  roll: number,
): Point2 | null => {
  const { faces, scale, offsetX, offsetY } = getDymaxionNetTransform(layout, canvasWidth, canvasHeight);
  const rotate = d3.geoRotation([lon, lat, roll]);
  const [sourceLon, sourceLat] = toLonLat(position);
  const inverted = rotate.invert([-sourceLon, sourceLat]);
  if (!inverted) return null;
  const p3 = lonLatToPoint3(inverted as Point2);

  let selectedFace = faces[0];
  let selectedInside = false;
  let bestScore = -Infinity;

  for (const face of faces) {
    const faceCenter = normalizeVec([
      face.vertices3D[0][0] + face.vertices3D[1][0] + face.vertices3D[2][0],
      face.vertices3D[0][1] + face.vertices3D[1][1] + face.vertices3D[2][1],
      face.vertices3D[0][2] + face.vertices3D[1][2] + face.vertices3D[2][2],
    ]);
    const score = dot3(p3, faceCenter);
    const inside = pointInsideSphericalFace(p3, face.vertices3D);
    if (inside) {
      selectedFace = face;
      selectedInside = true;
      break;
    }
    if (!selectedInside && score > bestScore) {
      bestScore = score;
      selectedFace = face;
    }
  }

  const weights = barycentric3D(
    p3,
    selectedFace.vertices3D[0],
    selectedFace.vertices3D[1],
    selectedFace.vertices3D[2],
  );
  if (!weights) return null;

  const clamped = selectedInside ? weights : weights.map(w => Math.max(0, w)) as [number, number, number];
  const weightSum = clamped[0] + clamped[1] + clamped[2] || 1;
  const u = clamped[0] / weightSum;
  const v = clamped[1] / weightSum;
  const w = clamped[2] / weightSum;
  const x = selectedFace.vertices[0][0] * u + selectedFace.vertices[1][0] * v + selectedFace.vertices[2][0] * w;
  const y = selectedFace.vertices[0][1] * u + selectedFace.vertices[1][1] * v + selectedFace.vertices[2][1] * w;

  return [x * scale + offsetX, y * scale + offsetY];
};

const getNearestCellId = (world: WorldData, lon: number, lat: number) => {
  const lonRad = lon * (Math.PI / 180);
  const latRad = lat * (Math.PI / 180);
  const cosLat = Math.cos(latRad);
  const x = cosLat * Math.cos(lonRad);
  const y = Math.sin(latRad);
  const z = cosLat * Math.sin(lonRad);

  let bestId: number | null = null;
  let bestDot = -Infinity;
  for (const cell of world.cells) {
    const d = cell.center.x * x + cell.center.y * y + cell.center.z * z;
    if (d > bestDot) {
      bestDot = d;
      bestId = cell.id;
    }
  }

  return bestId;
};

const getFactionOverlayData = (world: WorldData | null, visible: boolean): FactionOverlayData => {
  if (!world?.civData || !visible) return { labels: [], borders: [] };

  const labels = world.civData.factions.map((faction) => {
    const factionCells = world.cells.filter(cell => cell.regionId === faction.id);
    const landCells = factionCells.filter(cell => cell.height >= world.params.seaLevel);
    const labelCells = landCells.length > 0 ? landCells : factionCells;
    if (labelCells.length === 0) return null;

    let sumX = 0;
    let sumY = 0;
    let sumZ = 0;
    labelCells.forEach((cell) => {
      sumX += cell.center.x;
      sumY += cell.center.y;
      sumZ += cell.center.z;
    });

    return {
      id: faction.id,
      name: faction.name,
      position: normalizeVec([sumX, sumY, sumZ]),
    };
  }).filter(Boolean) as FactionLabel[];

  const borders: Array<[Point3, Point3]> = [];
  const threshold = 0.000001;
  world.cells.forEach((cellA) => {
    cellA.neighbors.forEach((nId) => {
      const cellB = world.cells[nId];
      if (!cellB || cellA.id >= cellB.id) return;
      if (cellA.regionId === cellB.regionId) return;
      if (cellA.regionId === undefined && cellB.regionId === undefined) return;

      const shared: Point3[] = [];
      for (const vA of cellA.vertices) {
        for (const vB of cellB.vertices) {
          const distSq = (vA.x - vB.x) ** 2 + (vA.y - vB.y) ** 2 + (vA.z - vB.z) ** 2;
          if (distSq < threshold) {
            shared.push([vA.x, vA.y, vA.z]);
            break;
          }
        }
        if (shared.length === 2) break;
      }
      if (shared.length === 2) borders.push([shared[0], shared[1]]);
    });
  });

  return { labels, borders };
};

const drawFactionBorders = (ctx: CanvasRenderingContext2D, pathGenerator: d3.GeoPath, borders: Array<[Point3, Point3]>, lineWidth: number) => {
  if (borders.length === 0) return;

  const drawPass = () => {
    borders.forEach(([a, b]) => {
      ctx.beginPath();
      pathGenerator({
        type: 'LineString',
        coordinates: [toLonLat(a), toLonLat(b)],
      } as any);
      ctx.stroke();
    });
  };

  ctx.save();
  ctx.lineJoin = 'round';
  ctx.lineCap = 'round';
  ctx.globalAlpha = 0.95;
  ctx.strokeStyle = 'rgba(2, 6, 23, 0.9)';
  ctx.lineWidth = lineWidth * 2.5;
  drawPass();
  ctx.strokeStyle = 'rgba(248, 250, 252, 0.95)';
  ctx.lineWidth = lineWidth;
  drawPass();
  ctx.restore();
};

const drawFactionLabels = (
  ctx: CanvasRenderingContext2D,
  labels: FactionLabel[],
  project: (position: Point3) => Point2 | null,
  fontSize: number,
  maxWidth: number,
) => {
  if (labels.length === 0) return;

  ctx.save();
  ctx.font = `700 ${fontSize}px Inter, ui-sans-serif, system-ui, sans-serif`;
  ctx.textAlign = 'center';
  ctx.textBaseline = 'middle';
  ctx.lineJoin = 'round';
  ctx.lineWidth = Math.max(3, fontSize * 0.28);
  ctx.strokeStyle = 'rgba(2, 6, 23, 0.95)';
  ctx.fillStyle = '#f8fafc';
  labels.forEach((label) => {
    const projected = project(label.position);
    if (!projected) return;
    ctx.strokeText(label.name, projected[0], projected[1], maxWidth);
    ctx.fillText(label.name, projected[0], projected[1], maxWidth);
  });
  ctx.restore();
};

const rasterizeDymaxionSource = ({
  sourceData,
  sourceWidth,
  sourceHeight,
  outputData,
  outputWidth,
  outputHeight,
  canvasWidth,
  canvasHeight,
  renderDpr,
  layout,
  lon,
  lat,
  roll,
}: {
  sourceData: Uint8ClampedArray;
  sourceWidth: number;
  sourceHeight: number;
  outputData: Uint8ClampedArray;
  outputWidth: number;
  outputHeight: number;
  canvasWidth: number;
  canvasHeight: number;
  renderDpr: number;
  layout: DymaxionSettings['layout'];
  lon: number;
  lat: number;
  roll: number;
}) => {
  const { faces, scale, offsetX, offsetY } = getDymaxionNetTransform(layout, canvasWidth, canvasHeight);
  const rotate = d3.geoRotation([lon, lat, roll]);

  faces.forEach((face) => {
    const verts = face.vertices.map((v) => [v[0] * scale + offsetX, v[1] * scale + offsetY]) as Point2[];
    const [a, b, c] = verts;
    const minBX = Math.max(0, Math.floor(Math.min(a[0], b[0], c[0])));
    const maxBX = Math.min(canvasWidth - 1, Math.ceil(Math.max(a[0], b[0], c[0])));
    const minBY = Math.max(0, Math.floor(Math.min(a[1], b[1], c[1])));
    const maxBY = Math.min(canvasHeight - 1, Math.ceil(Math.max(a[1], b[1], c[1])));

    const startOY = Math.floor(minBY * renderDpr);
    const endOY = Math.min(outputHeight - 1, Math.ceil(maxBY * renderDpr));
    const startOX = Math.floor(minBX * renderDpr);
    const endOX = Math.min(outputWidth - 1, Math.ceil(maxBX * renderDpr));

    for (let oy = startOY; oy <= endOY; oy++) {
      for (let ox = startOX; ox <= endOX; ox++) {
        const x = ox / renderDpr;
        const y = oy / renderDpr;
        if (!insideTri([x, y], a, b, c)) continue;

        const netPoint: Point2 = [(x - offsetX) / scale, (y - offsetY) / scale];
        const weights = barycentric(netPoint, face.vertices[0], face.vertices[1], face.vertices[2]);
        if (!weights) continue;
        const [u, v, w] = weights;
        const v0 = face.vertices3D[0];
        const v1 = face.vertices3D[1];
        const v2 = face.vertices3D[2];
        const p3 = normalizeVec([
          u * v0[0] + v * v1[0] + w * v2[0],
          u * v0[1] + v * v1[1] + w * v2[1],
          u * v0[2] + v * v1[2] + w * v2[2],
        ]);
        const [rotatedLon, rotatedLat] = rotate(toLonLat(p3));
        const srcX = Math.min(sourceWidth - 1, Math.max(0, Math.floor((rotatedLon + 180) / 360 * sourceWidth)));
        const srcY = Math.min(sourceHeight - 1, Math.max(0, Math.floor((90 - rotatedLat) / 180 * sourceHeight)));
        const srcIdx = (srcY * sourceWidth + srcX) * 4;
        const outIdx = (oy * outputWidth + ox) * 4;
        if (outIdx >= 0 && outIdx < outputData.length - 3) {
          outputData[outIdx] = sourceData[srcIdx];
          outputData[outIdx + 1] = sourceData[srcIdx + 1];
          outputData[outIdx + 2] = sourceData[srcIdx + 2];
          outputData[outIdx + 3] = 255;
        }
      }
    }
  });
};

const Map2D: React.FC<{
  world: WorldData | null;
  viewMode: ViewMode;
  inspectMode: InspectMode;
  onInspect: (cellId: number | null) => void;
  highlightCellId?: number | null;
  projectionType?: 'mercator' | 'dymaxion';
  dymaxionSettings?: DymaxionSettings;
  showGrid?: boolean;
  showRivers?: boolean;
  showFactionOverlay?: boolean;
  editMode?: EditMode;
  onPaint?: (cellId: number, phase: 'start' | 'stroke' | 'end', isRightClick?: boolean) => void;
  factionColors?: Map<number, string>;
  brushSize?: number;
}> = ({ world, viewMode, inspectMode, onInspect, highlightCellId = null, projectionType = 'mercator', dymaxionSettings, showGrid = false, showRivers = true, showFactionOverlay = true, editMode = 'off', onPaint, factionColors, brushSize = 1 }) => {
  const containerRef = useRef<HTMLDivElement>(null);
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const offscreenRef = useRef<HTMLCanvasElement | null>(null);
  const pickCanvasRef = useRef<HTMLCanvasElement | null>(null);
  const pickCtxRef = useRef<CanvasRenderingContext2D | null>(null);
  const [size, setSize] = useState<Size>({ width: 0, height: 0 });
  const [scale, setScale] = useState(1);
  const [offset, setOffset] = useState({ x: 0, y: 0 });
  const dragging = useRef(false);
  const dragDistance = useRef(0);
  const lastPos = useRef({ x: 0, y: 0 });
  const isPaintingRef = useRef(false);
  const lastPaintCell = useRef<number | null>(null);
  const [isSpaceHeld, setIsSpaceHeld] = useState(false);
  const [mousePos, setMousePos] = useState<{ x: number; y: number } | null>(null);
  const dymaxionLayout = dymaxionSettings?.layout || 'classic';
  const dymaxionLon = dymaxionSettings?.lon ?? 0;
  const dymaxionLat = dymaxionSettings?.lat ?? 0;
  const dymaxionRoll = dymaxionSettings?.roll ?? 0;
  const factionOverlay = useMemo(
    () => getFactionOverlayData(world, showFactionOverlay),
    [world, showFactionOverlay],
  );

  useEffect(() => {
    const onDown = (e: KeyboardEvent) => { if (e.code === 'Space' && !e.repeat) { e.preventDefault(); setIsSpaceHeld(true); } };
    const onUp = (e: KeyboardEvent) => { if (e.code === 'Space') setIsSpaceHeld(false); };
    window.addEventListener('keydown', onDown);
    window.addEventListener('keyup', onUp);
    return () => { window.removeEventListener('keydown', onDown); window.removeEventListener('keyup', onUp); };
  }, []);
  const [isInteracting, setIsInteracting] = useState(false);
  const [qualityDpr, setQualityDpr] = useState(1);
  const rafId = useRef<number | null>(null);
  const pendingOffset = useRef<{ x: number; y: number } | null>(null);
  const settleTimer = useRef<number | null>(null);
  const wheelTimer = useRef<number | null>(null);
  const [renderCount, setRenderCount] = useState(0);
  const inspectEnabled = inspectMode === 'click';

  useEffect(() => {
    if (!containerRef.current) return;
    const ro = new ResizeObserver(entries => {
      const entry = entries[0];
      setSize({
        width: Math.floor(entry.contentRect.width),
        height: Math.floor(entry.contentRect.height)
      });
    });
    ro.observe(containerRef.current);
    return () => { ro.disconnect(); };
  }, []);

  useEffect(() => {
    if (!size.width || !size.height) return;
    setScale(1);
    setOffset({ x: 0, y: 0 });
  }, [size.width, size.height, world?.params.seed]);

  useEffect(() => {
    if (settleTimer.current) {
      window.clearTimeout(settleTimer.current);
      settleTimer.current = null;
    }
    if (isInteracting) {
      setQualityDpr(INTERACTION_DPR);
      return;
    }
    settleTimer.current = window.setTimeout(() => {
      const baseDpr = Math.min(2, window.devicePixelRatio || 1);
      const scaled = baseDpr * Math.min(scale, MAX_SHARP_SCALE);
      const target = Math.min(MAX_SHARP_DPR, scaled);
      setQualityDpr(target);
      settleTimer.current = null;
    }, SETTLE_MS);
  }, [isInteracting, scale]);

  const projection = useMemo(() => {
    if (!size.width || !size.height) return null;
    if (projectionType === 'dymaxion') return null;
    return d3.geoMercator().fitSize([size.width, size.height], { type: 'Sphere' } as d3.GeoPermissibleObjects);
  }, [size.width, size.height, projectionType]);

  useEffect(() => {
    if (!world || !size.width || !size.height) return;
    const offscreen = offscreenRef.current ?? document.createElement('canvas');
    offscreenRef.current = offscreen;

    const renderDpr = qualityDpr;
    const ctx = offscreen.getContext('2d');
    if (!ctx) return;

    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.clearRect(0, 0, offscreen.width, offscreen.height);
    ctx.setTransform(renderDpr, 0, 0, renderDpr, 0, 0);
    ctx.translate(size.width, 0);
    ctx.scale(-1, 1);

    if (projectionType === 'dymaxion') {
      const srcWidth = Math.max(1, Math.floor(size.width * renderDpr));
      const srcHeight = Math.max(1, Math.round(srcWidth / 2));
      const source = document.createElement('canvas');
      source.width = srcWidth;
      source.height = srcHeight;
      const srcCtx = source.getContext('2d');
      if (!srcCtx) return;
      const projection = d3.geoEquirectangular().fitSize([srcWidth, srcHeight], { type: 'Sphere' } as d3.GeoPermissibleObjects);
      const pathGenerator = d3.geoPath(projection, srcCtx);
      srcCtx.fillStyle = viewMode === 'satellite' || viewMode === 'biome' ? '#050505' : '#000000';
      srcCtx.fillRect(0, 0, srcWidth, srcHeight);
      srcCtx.save();
      srcCtx.translate(srcWidth, 0);
      srcCtx.scale(-1, 1);
      for (let i = 0; i < world.cells.length; i++) {
        const feature = world.geoJson?.features?.[i];
        if (!feature || !feature.geometry) continue;
        const color = getCellColor(world.cells[i], viewMode, world.params.seaLevel, factionColors);
        const hexColor = '#' + color.getHexString();
        srcCtx.beginPath();
        pathGenerator(feature);
        srcCtx.fillStyle = hexColor;
        srcCtx.strokeStyle = hexColor;
        srcCtx.lineWidth = 1;
        srcCtx.fill();
        srcCtx.stroke();
      }

      if (highlightCellId !== null) {
        const feature = world.geoJson?.features?.[highlightCellId];
        if (feature && feature.geometry) {
          srcCtx.save();
          srcCtx.strokeStyle = '#f9a8a8';
          srcCtx.lineWidth = Math.max(2, 3 * renderDpr);
          srcCtx.lineJoin = 'round';
          srcCtx.lineCap = 'round';
          srcCtx.beginPath();
          pathGenerator(feature);
          srcCtx.stroke();
          srcCtx.restore();
        }
      }

      // Draw Grid on source equirectangular canvas
      if (showGrid) {
        srcCtx.strokeStyle = 'rgba(255,255,255,0.15)';
        srcCtx.lineWidth = 1;
        srcCtx.beginPath();
        pathGenerator(d3.geoGraticule().step([10, 10])());
        srcCtx.stroke();
      }

      // Draw Rivers on source equirectangular canvas
      if (showRivers && world.rivers) {
        srcCtx.strokeStyle = '#38bdf8';
        srcCtx.lineWidth = Math.max(0.5, 1.5 / renderDpr);
        srcCtx.globalAlpha = 0.8;
        world.rivers.forEach(path => {
          if (path.length < 2) return;
          srcCtx.beginPath();
          let lastLon: number | null = null;
          path.forEach((p, i) => {
            const lon = Math.atan2(p.z, p.x) * (180 / Math.PI);
            const lat = Math.asin(Math.max(-1, Math.min(1, p.y))) * (180 / Math.PI);
            
            // Detect antimeridian crossing
            const isJump = lastLon !== null && Math.abs(lon - lastLon) > 180;
            
            const pt = projection([lon, lat]);
            if (pt) {
              if (i === 0 || isJump) srcCtx.moveTo(pt[0], pt[1]);
              else srcCtx.lineTo(pt[0], pt[1]);
            }
            lastLon = lon;
          });
          srcCtx.stroke();
        });
        srcCtx.globalAlpha = 1.0;
      }

      drawFactionBorders(
        srcCtx,
        pathGenerator,
        factionOverlay.borders,
        Math.max(1.5, 2 * renderDpr),
      );

      srcCtx.restore();

      const srcImage = srcCtx.getImageData(0, 0, srcWidth, srcHeight);
      const srcData = srcImage.data;

      const canvasWidth = size.width;
      const canvasHeight = size.height;
      const outWidth = Math.floor(canvasWidth * renderDpr);
      const outHeight = Math.floor(canvasHeight * renderDpr);

      // Create a temporary canvas for heavy rendering to avoid clearing the screen too early
      const tempCanvas = document.createElement('canvas');
      tempCanvas.width = outWidth;
      tempCanvas.height = outHeight;
      const tCtx = tempCanvas.getContext('2d');
      if (!tCtx) return;

      tCtx.fillStyle = viewMode === 'satellite' || viewMode === 'biome' ? '#050505' : '#000000';
      tCtx.fillRect(0, 0, outWidth, outHeight);

      const output = tCtx.getImageData(0, 0, outWidth, outHeight);
      const outData = output.data;

      rasterizeDymaxionSource({
        sourceData: srcData,
        sourceWidth: srcWidth,
        sourceHeight: srcHeight,
        outputData: outData,
        outputWidth: outWidth,
        outputHeight: outHeight,
        canvasWidth,
        canvasHeight,
        renderDpr,
        layout: dymaxionLayout,
        lon: dymaxionLon,
        lat: dymaxionLat,
        roll: dymaxionRoll,
      });

      tCtx.putImageData(output, 0, 0);
      drawFactionLabels(
        tCtx,
        factionOverlay.labels,
        (position) => {
          const projected = projectDymaxionPoint(
            position,
            dymaxionLayout,
            canvasWidth,
            canvasHeight,
            dymaxionLon,
            dymaxionLat,
            dymaxionRoll,
          );
          return projected ? [projected[0] * renderDpr, projected[1] * renderDpr] : null;
        },
        Math.max(12, 12 * renderDpr),
        Math.max(120, 120 * renderDpr),
      );
      
      // Update the main offscreen canvas only once the heavy rendering is complete
      offscreen.width = outWidth;
      offscreen.height = outHeight;
      const finalCtx = offscreen.getContext('2d');
      if (finalCtx) {
        finalCtx.drawImage(tempCanvas, 0, 0);
        setRenderCount(c => c + 1);
      }
      return;
    }

    if (!projection) return;
    offscreen.width = Math.max(1, Math.floor(size.width * renderDpr));
    offscreen.height = Math.max(1, Math.floor(size.height * renderDpr));
    // Reset ctx state after resize
    ctx.setTransform(renderDpr, 0, 0, renderDpr, 0, 0);
    ctx.translate(size.width, 0);
    ctx.scale(-1, 1);
    
    const pathGenerator = d3.geoPath(projection, ctx);

    for (let i = 0; i < world.cells.length; i++) {
        const feature = world.geoJson?.features?.[i];
      if (!feature || !feature.geometry) continue;
        const color = getCellColor(world.cells[i], viewMode, world.params.seaLevel, factionColors);
        const hexColor = '#' + color.getHexString();
      ctx.beginPath();
      pathGenerator(feature);
      ctx.fillStyle = hexColor;
      ctx.strokeStyle = hexColor;
      ctx.lineWidth = 1;
      ctx.fill();
      ctx.stroke();
    }

    // Draw Grid
    if (showGrid) {
      ctx.strokeStyle = 'rgba(255,255,255,0.15)';
      ctx.lineWidth = 1;
      ctx.beginPath();
      pathGenerator(d3.geoGraticule().step([10, 10])());
      ctx.stroke();
    }

    // Draw Rivers
    if (showRivers && world.rivers) {
      ctx.strokeStyle = '#38bdf8';
      ctx.lineWidth = 1.5 / qualityDpr;
      ctx.globalAlpha = 0.8;
      world.rivers.forEach(path => {
        if (path.length < 2) return;
        ctx.beginPath();
        let lastLon: number | null = null;
        path.forEach((p, i) => {
          const lon = Math.atan2(p.z, p.x) * (180 / Math.PI);
          const lat = Math.asin(Math.max(-1, Math.min(1, p.y))) * (180 / Math.PI);
          
          // Detect antimeridian crossing
          const isJump = lastLon !== null && Math.abs(lon - lastLon) > 180;
          
          const pt = projection([lon, lat]);
          if (pt) {
            if (i === 0 || isJump) ctx.moveTo(pt[0], pt[1]);
            else ctx.lineTo(pt[0], pt[1]);
          }
          lastLon = lon;
        });
        ctx.stroke();
      });
      ctx.globalAlpha = 1.0;
    }

    drawFactionBorders(
      ctx,
      pathGenerator,
      factionOverlay.borders,
      Math.max(1, 2 / qualityDpr),
    );

    if (viewMode === 'political') {
      ctx.save();
      ctx.globalAlpha = 0.25;
      ctx.strokeStyle = '#ffffff';
      ctx.lineWidth = 0.5;
      for (let i = 0; i < world.cells.length; i++) {
        const feature = world.geoJson?.features?.[i];
        if (!feature || !feature.geometry) continue;
        ctx.beginPath();
        pathGenerator(feature);
        ctx.stroke();
      }
      ctx.restore();
    }

    if (highlightCellId !== null) {
      const feature = world.geoJson?.features?.[highlightCellId];
      if (feature && feature.geometry) {
        ctx.save();
        ctx.strokeStyle = '#f9a8a8';
        ctx.lineWidth = 3 / qualityDpr;
        ctx.lineJoin = 'round';
        ctx.lineCap = 'round';
        ctx.beginPath();
        pathGenerator(feature);
        ctx.stroke();
        ctx.restore();
      }
    }

    if (factionOverlay.labels.length > 0) {
      ctx.save();
      ctx.setTransform(renderDpr, 0, 0, renderDpr, 0, 0);
      drawFactionLabels(
        ctx,
        factionOverlay.labels,
        (position) => {
          const projected = projection(toLonLat(position));
          return projected ? [size.width - projected[0], projected[1]] : null;
        },
        12,
        120,
      );
      ctx.restore();
    }
  }, [
    projection,
    size.width,
    size.height,
    world,
    viewMode,
    qualityDpr,
    highlightCellId,
    projectionType,
    dymaxionLayout,
    dymaxionLon,
    dymaxionLat,
    dymaxionRoll,
    showGrid,
    showRivers,
    factionOverlay,
    factionColors,
  ]);

  useEffect(() => {
    if (!world || !size.width || !size.height || projectionType !== 'dymaxion') {
      pickCtxRef.current = null;
      return;
    }
    const pickCanvas = pickCanvasRef.current ?? document.createElement('canvas');
    pickCanvasRef.current = pickCanvas;
    pickCanvas.width = size.width;
    pickCanvas.height = size.height;
    const ctx = pickCanvas.getContext('2d', { willReadFrequently: true });
    if (!ctx) return;
    pickCtxRef.current = ctx;

    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.clearRect(0, 0, size.width, size.height);

    const srcWidth = Math.max(1, Math.floor(size.width));
    const srcHeight = Math.max(1, Math.round(srcWidth / 2));
    const source = document.createElement('canvas');
    source.width = srcWidth;
    source.height = srcHeight;
    const srcCtx = source.getContext('2d', { willReadFrequently: true });
    if (!srcCtx) return;

    const sourceProjection = d3.geoEquirectangular().fitSize(
      [srcWidth, srcHeight],
      { type: 'Sphere' } as d3.GeoPermissibleObjects,
    );
    const pathGenerator = d3.geoPath(sourceProjection, srcCtx);

    srcCtx.fillStyle = 'rgb(0,0,0)';
    srcCtx.fillRect(0, 0, srcWidth, srcHeight);
    srcCtx.save();
    srcCtx.translate(srcWidth, 0);
    srcCtx.scale(-1, 1);

    for (let i = 0; i < world.cells.length; i++) {
      const feature = world.geoJson?.features?.[i];
      if (!feature || !feature.geometry) continue;
      const id = i + 1;
      const r = id & 255;
      const g = (id >> 8) & 255;
      const b = (id >> 16) & 255;
      const pickColor = `rgb(${r},${g},${b})`;
      srcCtx.beginPath();
      pathGenerator(feature);
      srcCtx.fillStyle = pickColor;
      srcCtx.strokeStyle = pickColor;
      srcCtx.lineWidth = 1;
      srcCtx.fill();
      srcCtx.stroke();
    }

    srcCtx.restore();

    ctx.fillStyle = 'rgb(0,0,0)';
    ctx.fillRect(0, 0, size.width, size.height);
    const srcImage = srcCtx.getImageData(0, 0, srcWidth, srcHeight);
    const output = ctx.getImageData(0, 0, size.width, size.height);
    rasterizeDymaxionSource({
      sourceData: srcImage.data,
      sourceWidth: srcWidth,
      sourceHeight: srcHeight,
      outputData: output.data,
      outputWidth: size.width,
      outputHeight: size.height,
      canvasWidth: size.width,
      canvasHeight: size.height,
      renderDpr: 1,
      layout: dymaxionLayout,
      lon: dymaxionLon,
      lat: dymaxionLat,
      roll: dymaxionRoll,
    });
    ctx.putImageData(output, 0, 0);
    // The pick buffer encodes cell IDs from the world's structure (cells +
    // geoJson are stable references across paint strokes), so keying on
    // world.cells instead of world identity skips a full-canvas per-pixel
    // reprojection on every stroke event while painting
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [
    size.width,
    size.height,
    world?.cells,
    projectionType,
    dymaxionLayout,
    dymaxionLon,
    dymaxionLat,
    dymaxionRoll,
  ]);

  useEffect(() => {
    const canvas = canvasRef.current;
    const offscreen = offscreenRef.current;
    if (!canvas || !offscreen || !size.width || !size.height) return;
    const ctx = canvas.getContext('2d');
    if (!ctx) return;

    const displayDpr = qualityDpr;
    canvas.width = Math.max(1, Math.floor(size.width * displayDpr));
    canvas.height = Math.max(1, Math.floor(size.height * displayDpr));

    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    ctx.setTransform(displayDpr * scale, 0, 0, displayDpr * scale, displayDpr * offset.x, displayDpr * offset.y);
    ctx.drawImage(offscreen, 0, 0, size.width, size.height);
  }, [size.width, size.height, scale, offset.x, offset.y, qualityDpr, viewMode, world?.params.seed, world, projectionType, renderCount]);

  const scaleRef = useRef(scale);
  const offsetRef = useRef(offset);

  useEffect(() => {
    scaleRef.current = scale;
  }, [scale]);

  useEffect(() => {
    offsetRef.current = offset;
  }, [offset]);

  const handleWheel = useCallback((event: WheelEvent) => {
    event.preventDefault();
    setIsInteracting(true);
    if (wheelTimer.current) window.clearTimeout(wheelTimer.current);
    wheelTimer.current = window.setTimeout(() => {
      setIsInteracting(false);
      wheelTimer.current = null;
    }, 180);
    const rect = canvasRef.current?.getBoundingClientRect();
    if (!rect) return;
    const mx = event.clientX - rect.left;
    const my = event.clientY - rect.top;

    const prevScale = scaleRef.current;
    const zoomFactor = event.deltaY < 0 ? 1.1 : 0.9;
    const nextScale = clamp(prevScale * zoomFactor, 0.6, 6);

    const prevOffset = offsetRef.current;
    const worldX = (mx - prevOffset.x) / prevScale;
    const worldY = (my - prevOffset.y) / prevScale;

    const nextOffsetX = mx - worldX * nextScale;
    const nextOffsetY = my - worldY * nextScale;

    setScale(nextScale);
    setOffset({ x: nextOffsetX, y: nextOffsetY });
  }, []);

  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    const listener = (event: WheelEvent) => { handleWheel(event); };
    canvas.addEventListener('wheel', listener, { passive: false });
    return () => { canvas.removeEventListener('wheel', listener); };
  }, [handleWheel]);

  const isPaintMode = editMode !== 'off' && editMode !== 'world-edit';

  const getPickBufferCellId = useCallback((mapX: number, mapY: number) => {
    if (!pickCtxRef.current || !world) return null;
    const data = pickCtxRef.current.getImageData(Math.floor(mapX), Math.floor(mapY), 1, 1).data;
    const id = data[0] + (data[1] << 8) + (data[2] << 16);
    if (id <= 0 || id > world.cells.length) return null;
    return id - 1;
  }, [world]);

  const getCellIdAtMapPoint = useCallback((mapX: number, mapY: number) => {
    if (!world) return null;

    if (projectionType === 'dymaxion') {
      return getPickBufferCellId(mapX, mapY);
    }

    if (!projection) return null;
    const lonLat = projection.invert?.([size.width - mapX, mapY]);
    if (!lonLat) return null;
    const [lon, lat] = lonLat;
    if (!Number.isFinite(lon) || !Number.isFinite(lat)) return null;
    return getNearestCellId(world, lon, lat);
  }, [getPickBufferCellId, projection, projectionType, size.width, world]);

  const paintPickAt = useCallback((clientX: number, clientY: number, phase: 'start' | 'stroke' | 'end') => {
    if (!canvasRef.current || !onPaint) return;
    const rect = canvasRef.current.getBoundingClientRect();
    const mapX = (clientX - rect.left - offset.x) / scale;
    const mapY = (clientY - rect.top - offset.y) / scale;
    if (mapX < 0 || mapY < 0 || mapX >= size.width || mapY >= size.height) return;
    const cellId = getCellIdAtMapPoint(mapX, mapY);
    if (cellId === null) return;
    lastPaintCell.current = cellId;
    onPaint(cellId, phase);
  }, [getCellIdAtMapPoint, offset.x, offset.y, scale, size.width, size.height, onPaint]);

  const handleMouseDown = (e: React.MouseEvent<HTMLCanvasElement>) => {
    if (e.button === 1) e.preventDefault();
    dragging.current = true;
    dragDistance.current = 0;
    lastPos.current = { x: e.clientX, y: e.clientY };
    setIsInteracting(true);

    if (isPaintMode && !isSpaceHeld && e.button !== 1) {
      const isRight = e.button === 2;
      if (isRight) {
        // Right-click: send start with isRightClick=true, don't begin a stroke
        if (!canvasRef.current) return;
        const rect = canvasRef.current.getBoundingClientRect();
        const mapX = (e.clientX - rect.left - offset.x) / scale;
        const mapY = (e.clientY - rect.top - offset.y) / scale;
        if (mapX >= 0 && mapY >= 0 && mapX < size.width && mapY < size.height) {
          const cellId = getCellIdAtMapPoint(mapX, mapY);
          if (cellId !== null && onPaint) onPaint(cellId, 'start', true);
        }
        return;
      }
      isPaintingRef.current = true;
      paintPickAt(e.clientX, e.clientY, 'start');
      paintPickAt(e.clientX, e.clientY, 'stroke');
    }
  };

  const handleMouseMove = (e: React.MouseEvent<HTMLCanvasElement>) => {
    // Track mouse for brush circle overlay
    if (isPaintMode) {
      const rect = canvasRef.current?.getBoundingClientRect();
      if (rect) setMousePos({ x: e.clientX - rect.left, y: e.clientY - rect.top });
    }

    if (!dragging.current) return;
    const dx = e.clientX - lastPos.current.x;
    const dy = e.clientY - lastPos.current.y;
    lastPos.current = { x: e.clientX, y: e.clientY };
    dragDistance.current += Math.abs(dx) + Math.abs(dy);

    if (isPaintMode && isPaintingRef.current && !isSpaceHeld) {
      paintPickAt(e.clientX, e.clientY, 'stroke');
      return; // skip pan during paint stroke
    }

    const next = { x: offset.x + dx, y: offset.y + dy };
    pendingOffset.current = next;
    if (rafId.current === null) {
      rafId.current = requestAnimationFrame(() => {
        rafId.current = null;
        if (pendingOffset.current) {
          setOffset(pendingOffset.current);
          pendingOffset.current = null;
        }
      });
    }
  };

  const pickAt = (clientX: number, clientY: number, clearOnMiss = true) => {
    if (!canvasRef.current) return;
    const rect = canvasRef.current.getBoundingClientRect();
    const mapX = (clientX - rect.left - offset.x) / scale;
    const mapY = (clientY - rect.top - offset.y) / scale;
    if (mapX < 0 || mapY < 0 || mapX >= size.width || mapY >= size.height) {
      if (clearOnMiss) onInspect(null);
      return;
    }
    const cellId = getCellIdAtMapPoint(mapX, mapY);
    if (cellId === null) {
      if (clearOnMiss) onInspect(null);
      return;
    }
    onInspect(cellId);
  };

  const handleHover = (e: React.MouseEvent<HTMLCanvasElement>) => {
    if (inspectMode !== 'hover' || dragging.current) return;
    pickAt(e.clientX, e.clientY);
  };

  const endDrag = () => {
    if (isPaintMode && isPaintingRef.current) {
      if (lastPaintCell.current !== null && onPaint) {
        onPaint(lastPaintCell.current, 'end');
      }
      isPaintingRef.current = false;
      lastPaintCell.current = null;
    }
    dragging.current = false;
    setIsInteracting(false);
  };

  const handleMouseUp = (e: React.MouseEvent<HTMLCanvasElement>) => {
    endDrag();
    if (inspectEnabled && dragDistance.current < 6) {
      pickAt(e.clientX, e.clientY, false);
    }
  };

  return (
    <div ref={containerRef} className="w-full h-full bg-black relative">
      <canvas
        ref={canvasRef}
        className="w-full h-full cursor-grab active:cursor-grabbing"
        onMouseDown={handleMouseDown}
        onMouseMove={(e) => { handleMouseMove(e); handleHover(e); }}
        onMouseUp={handleMouseUp}
        onMouseLeave={() => { endDrag(); setMousePos(null); if (inspectMode === 'hover') onInspect(null); }}
        onAuxClick={e => { if (e.button === 1) e.preventDefault(); }}
        onContextMenu={e => e.preventDefault()}
      />
      {/* Brush size circle overlay */}
      {isPaintMode && mousePos && world && (() => {
        const avgCellPx = Math.sqrt((size.width * size.height) / world.cells.length / Math.PI);
        const brushPx = Math.max(8, avgCellPx * (brushSize + 0.5) * scale);
        return (
          <div
            className="pointer-events-none absolute border-2 border-white/70 rounded-full"
            style={{
              left: mousePos.x - brushPx,
              top: mousePos.y - brushPx,
              width: brushPx * 2,
              height: brushPx * 2,
            }}
          />
        );
      })()}
      {!world && (
        <div className="absolute inset-0 flex items-center justify-center text-white/50">
          Forging World...
        </div>
      )}
      {world && (
        <div className="absolute bottom-4 right-4 text-[10px] bg-black/60 border border-white/10 px-2 py-1 text-gray-300">
          {projectionType === 'dymaxion' ? '2D Dymaxion' : '2D Mercator'} • Scroll to zoom, drag to pan
        </div>
      )}
    </div>
  );
};

export default Map2D;
