import React, { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import * as d3 from 'd3';
import { WorldData, ViewMode, InspectMode, DymaxionSettings } from '../types';
import { getCellColor } from '../utils/colors';
import { buildDymaxionNet } from '../utils/dymaxion';

type Size = { width: number; height: number };

const clamp = (v: number, min: number, max: number) => Math.max(min, Math.min(max, v));
const INTERACTION_DPR = 1;
const MAX_SHARP_DPR = 3;
const MAX_SHARP_SCALE = 2.5;
const SETTLE_MS = 200;

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
}> = ({ world, viewMode, inspectMode, onInspect, highlightCellId = null, projectionType = 'mercator', dymaxionSettings, showGrid = false, showRivers = true }) => {
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
        // eslint-disable-next-line security/detect-object-injection
        const feature = world.geoJson?.features?.[i];
        if (!feature || !feature.geometry) continue;
        // eslint-disable-next-line security/detect-object-injection
        const color = getCellColor(world.cells[i], viewMode, world.params.seaLevel);
        srcCtx.beginPath();
        pathGenerator(feature);
        srcCtx.fillStyle = '#' + color.getHexString();
        srcCtx.fill();
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

      srcCtx.restore();

      const srcImage = srcCtx.getImageData(0, 0, srcWidth, srcHeight);
      const srcData = srcImage.data;

      const canvasWidth = size.width;
      const canvasHeight = size.height;
      const outWidth = Math.floor(canvasWidth * renderDpr);
      const outHeight = Math.floor(canvasHeight * renderDpr);

      const net = buildDymaxionNet(dymaxionSettings?.layout || 'classic');
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

      const rotate = dymaxionSettings ? d3.geoRotation([dymaxionSettings.lon, dymaxionSettings.lat, dymaxionSettings.roll]) : null;

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

      const insideTri = (p: [number, number], a: [number, number], b: [number, number], c: [number, number]) => {
        const v0 = [c[0] - a[0], c[1] - a[1]];
        const v1 = [b[0] - a[0], b[1] - a[1]];
        const v2 = [p[0] - a[0], p[1] - a[1]];
        const dot00 = v0[0] * v0[0] + v0[1] * v0[1];
        const dot01 = v0[0] * v1[0] + v0[1] * v1[1];
        const dot02 = v0[0] * v2[0] + v0[1] * v2[1];
        const dot11 = v1[0] * v1[0] + v1[1] * v1[1];
        const dot12 = v1[0] * v2[0] + v1[1] * v2[1];
        const invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
        const u = (dot11 * dot02 - dot01 * dot12) * invDenom;
        const v = (dot00 * dot12 - dot01 * dot02) * invDenom;
        return u >= -1e-6 && v >= -1e-6 && u + v <= 1 + 1e-6;
      };

      const barycentric = (p: [number, number], a: [number, number], b: [number, number], c: [number, number]) => {
        const v0 = [b[0] - a[0], b[1] - a[1]];
        const v1 = [c[0] - a[0], c[1] - a[1]];
        const v2 = [p[0] - a[0], p[1] - a[1]];
        const d00 = v0[0] * v0[0] + v0[1] * v0[1];
        const d01 = v0[0] * v1[0] + v0[1] * v1[1];
        const d11 = v1[0] * v1[0] + v1[1] * v1[1];
        const d20 = v2[0] * v0[0] + v2[1] * v0[1];
        const d21 = v2[0] * v1[0] + v2[1] * v1[1];
        const denom = d00 * d11 - d01 * d01;
        if (!denom) return null;
        const v = (d11 * d20 - d01 * d21) / denom;
        const w = (d00 * d21 - d01 * d20) / denom;
        const u = 1 - v - w;
        return [u, v, w] as [number, number, number];
      };

      const normalizeVec = (v: [number, number, number]) => {
        const len = Math.hypot(v[0], v[1], v[2]) || 1;
        return [v[0] / len, v[1] / len, v[2] / len] as [number, number, number];
      };

      const toLonLat = (v: [number, number, number]) => {
        const lon = Math.atan2(v[2], v[0]) * (180 / Math.PI);
        const lat = Math.asin(Math.max(-1, Math.min(1, v[1]))) * (180 / Math.PI);
        return [lon, lat] as [number, number];
      };

      faces.forEach((face) => {
        const verts = face.vertices.map((v) => [v[0] * scale + offsetX, v[1] * scale + offsetY]) as [number, number][];
        const [a, b, c] = verts;
        const minBX = Math.max(0, Math.floor(Math.min(a[0], b[0], c[0])));
        const maxBX = Math.min(canvasWidth - 1, Math.ceil(Math.max(a[0], b[0], c[0])));
        const minBY = Math.max(0, Math.floor(Math.min(a[1], b[1], c[1])));
        const maxBY = Math.min(canvasHeight - 1, Math.ceil(Math.max(a[1], b[1], c[1])));

        const startOY = Math.floor(minBY * renderDpr);
        const endOY = Math.min(outHeight - 1, Math.ceil(maxBY * renderDpr));
        const startOX = Math.floor(minBX * renderDpr);
        const endOX = Math.min(outWidth - 1, Math.ceil(maxBX * renderDpr));

        for (let oy = startOY; oy <= endOY; oy++) {
          for (let ox = startOX; ox <= endOX; ox++) {
            const x = ox / renderDpr;
            const y = oy / renderDpr;
            const p: [number, number] = [x, y];
            if (!insideTri(p, a, b, c)) continue;
            
            const netPoint: [number, number] = [(x - offsetX) / scale, (y - offsetY) / scale];
            const weights = barycentric(netPoint, face.vertices[0], face.vertices[1], face.vertices[2]);
            if (!weights) continue;
            const [u, v, w] = weights;
            const v0 = face.vertices3D[0];
            const v1 = face.vertices3D[1];
            const v2 = face.vertices3D[2];
            const p3 = normalizeVec([
              u * v0[0] + v * v1[0] + w * v2[0],
              u * v0[1] + v * v1[1] + w * v2[1],
              u * v0[2] + v * v1[2] + w * v2[2]
            ]);
            const lonLat = toLonLat(p3);
            const rotated = rotate ? rotate(lonLat) : lonLat;
            const lon = rotated[0];
            const lat = rotated[1];
            const srcX = Math.min(srcWidth - 1, Math.max(0, Math.floor((lon + 180) / 360 * srcWidth)));
            const srcY = Math.min(srcHeight - 1, Math.max(0, Math.floor((90 - lat) / 180 * srcHeight)));
            const srcIdx = (srcY * srcWidth + srcX) * 4;
            const outIdx = (oy * outWidth + ox) * 4;
            if (outIdx >= 0 && outIdx < outData.length - 3) {
              // eslint-disable-next-line security/detect-object-injection
              outData[outIdx] = srcData[srcIdx];
              outData[outIdx + 1] = srcData[srcIdx + 1];
              outData[outIdx + 2] = srcData[srcIdx + 2];
              outData[outIdx + 3] = 255;
            }
          }
        }
      });

      tCtx.putImageData(output, 0, 0);
      
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
      // eslint-disable-next-line security/detect-object-injection
        const feature = world.geoJson?.features?.[i];
      if (!feature || !feature.geometry) continue;
      // eslint-disable-next-line security/detect-object-injection
        const color = getCellColor(world.cells[i], viewMode, world.params.seaLevel);
      ctx.beginPath();
      pathGenerator(feature);
      ctx.fillStyle = '#' + color.getHexString();
      ctx.fill();
    }

    // Draw Grid (Mercator)
    if (showGrid) {
      ctx.strokeStyle = 'rgba(255,255,255,0.15)';
      ctx.lineWidth = 1;
      ctx.beginPath();
      pathGenerator(d3.geoGraticule().step([10, 10])());
      ctx.stroke();
    }

    // Draw Rivers (Mercator)
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

    if (viewMode === 'political') {
      ctx.save();
      ctx.globalAlpha = 0.25;
      ctx.strokeStyle = '#ffffff';
      ctx.lineWidth = 0.5;
      for (let i = 0; i < world.cells.length; i++) {
        // eslint-disable-next-line security/detect-object-injection
        const feature = world.geoJson?.features?.[i];
        if (!feature || !feature.geometry) continue;
        ctx.beginPath();
        pathGenerator(feature);
        ctx.stroke();
      }
      ctx.restore();
    }

    if (highlightCellId !== null) {
      // eslint-disable-next-line security/detect-object-injection
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
  }, [
    projection,
    size.width,
    size.height,
    world,
    viewMode,
    qualityDpr,
    highlightCellId,
    projectionType,
    dymaxionSettings?.layout,
    dymaxionSettings?.lon,
    dymaxionSettings?.lat,
    dymaxionSettings?.roll,
    showGrid,
    showRivers
  ]);

  useEffect(() => {
    if (projectionType === 'dymaxion') return;
    if (!world || !projection || !size.width || !size.height) return;
    const pickCanvas = pickCanvasRef.current ?? document.createElement('canvas');
    pickCanvasRef.current = pickCanvas;
    pickCanvas.width = size.width;
    pickCanvas.height = size.height;
    const ctx = pickCanvas.getContext('2d', { willReadFrequently: true });
    if (!ctx) return;
    pickCtxRef.current = ctx;

    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.clearRect(0, 0, size.width, size.height);
    ctx.translate(size.width, 0);
    ctx.scale(-1, 1);

    const pathGenerator = d3.geoPath(projection, ctx);
    for (let i = 0; i < world.cells.length; i++) {
      // eslint-disable-next-line security/detect-object-injection
        const feature = world.geoJson?.features?.[i];
      if (!feature || !feature.geometry) continue;
      const id = i + 1;
      const r = id & 255;
      const g = (id >> 8) & 255;
      const b = (id >> 16) & 255;
      ctx.beginPath();
      pathGenerator(feature);
      ctx.fillStyle = `rgb(${r},${g},${b})`;
      ctx.fill();
    }
  }, [projection, size.width, size.height, world, projectionType]);

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

  const handleMouseDown = (e: React.MouseEvent<HTMLCanvasElement>) => {
    dragging.current = true;
    dragDistance.current = 0;
    lastPos.current = { x: e.clientX, y: e.clientY };
    setIsInteracting(true);
  };

  const handleMouseMove = (e: React.MouseEvent<HTMLCanvasElement>) => {
    if (!dragging.current) return;
    const dx = e.clientX - lastPos.current.x;
    const dy = e.clientY - lastPos.current.y;
    lastPos.current = { x: e.clientX, y: e.clientY };
    dragDistance.current += Math.abs(dx) + Math.abs(dy);
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

  const pickAt = (clientX: number, clientY: number) => {
    if (!pickCtxRef.current || !canvasRef.current) return;
    const rect = canvasRef.current.getBoundingClientRect();
    const mapX = (clientX - rect.left - offset.x) / scale;
    const mapY = (clientY - rect.top - offset.y) / scale;
    if (mapX < 0 || mapY < 0 || mapX >= size.width || mapY >= size.height) {
      onInspect(null);
      return;
    }
    const data = pickCtxRef.current.getImageData(Math.floor(mapX), Math.floor(mapY), 1, 1).data;
    const id = data[0] + (data[1] << 8) + (data[2] << 16);
    onInspect(id === 0 ? null : id - 1);
  };

  const handleHover = (e: React.MouseEvent<HTMLCanvasElement>) => {
    if (inspectMode !== 'hover' || dragging.current) return;
    pickAt(e.clientX, e.clientY);
  };

  const handleClick = (e: React.MouseEvent<HTMLCanvasElement>) => {
    if (inspectMode !== 'click') return;
    pickAt(e.clientX, e.clientY);
  };

  const endDrag = () => {
    dragging.current = false;
    setIsInteracting(false);
  };

  const handleMouseUp = (e: React.MouseEvent<HTMLCanvasElement>) => {
    endDrag();
    if (inspectEnabled && dragDistance.current < 6) {
      pickAt(e.clientX, e.clientY);
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
        onClick={handleClick}
        onMouseLeave={() => { endDrag(); if (inspectMode === 'hover') onInspect(null); }}
      />
      {!world && (
        <div className="absolute inset-0 flex items-center justify-center text-white/50">
          Forging World...
        </div>
      )}
      {world && (
        <div className="absolute bottom-4 right-4 text-[10px] bg-black/60 border border-white/10 rounded px-2 py-1 text-gray-300">
          2D Mercator • Scroll to zoom, drag to pan
        </div>
      )}
    </div>
  );
};

export default Map2D;
