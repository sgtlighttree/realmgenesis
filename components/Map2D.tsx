import React, { useEffect, useMemo, useRef, useState } from 'react';
import * as d3 from 'd3';
import { WorldData, ViewMode, InspectMode } from '../types';
import { getCellColor } from '../utils/colors';

type Size = { width: number; height: number };

const clamp = (v: number, min: number, max: number) => Math.max(min, Math.min(max, v));
const INTERACTION_DPR = 1;
const MAX_SHARP_DPR = 3;
const MAX_SHARP_SCALE = 2.5;
const SETTLE_MS = 200;

const Map2D: React.FC<{ world: WorldData | null; viewMode: ViewMode; inspectMode: InspectMode; onInspect: (cellId: number | null) => void; }> = ({ world, viewMode, inspectMode, onInspect }) => {
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
  const hoverRaf = useRef<number | null>(null);
  const pendingHover = useRef<{ x: number; y: number } | null>(null);

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
    return () => ro.disconnect();
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
    return d3.geoMercator().fitSize([size.width, size.height], { type: 'Sphere' } as any);
  }, [size.width, size.height]);

  useEffect(() => {
    if (!world || !projection || !size.width || !size.height) return;
    const offscreen = offscreenRef.current ?? document.createElement('canvas');
    offscreenRef.current = offscreen;

    offscreen.width = Math.max(1, Math.floor(size.width * qualityDpr));
    offscreen.height = Math.max(1, Math.floor(size.height * qualityDpr));
    const ctx = offscreen.getContext('2d');
    if (!ctx) return;

    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.clearRect(0, 0, offscreen.width, offscreen.height);
    ctx.setTransform(qualityDpr, 0, 0, qualityDpr, 0, 0);
    ctx.translate(size.width, 0);
    ctx.scale(-1, 1);

    const pathGenerator = d3.geoPath(projection, ctx);

    for (let i = 0; i < world.cells.length; i++) {
      const feature = world.geoJson?.features?.[i];
      if (!feature || !feature.geometry) continue;
      const color = getCellColor(world.cells[i], viewMode, world.params.seaLevel);
      ctx.beginPath();
      pathGenerator(feature);
      ctx.fillStyle = '#' + color.getHexString();
      ctx.fill();
    }

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
  }, [projection, size.width, size.height, world, viewMode, qualityDpr]);

  useEffect(() => {
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
  }, [projection, size.width, size.height, world]);

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
  }, [size.width, size.height, scale, offset.x, offset.y, qualityDpr, viewMode, world?.params.seed]);

  const handleWheel = (e: React.WheelEvent<HTMLCanvasElement>) => {
    e.preventDefault();
    setIsInteracting(true);
    if (wheelTimer.current) window.clearTimeout(wheelTimer.current);
    wheelTimer.current = window.setTimeout(() => {
      setIsInteracting(false);
      wheelTimer.current = null;
    }, 180);
    const rect = e.currentTarget.getBoundingClientRect();
    const mx = e.clientX - rect.left;
    const my = e.clientY - rect.top;

    const prevScale = scale;
    const zoomFactor = e.deltaY < 0 ? 1.1 : 0.9;
    const nextScale = clamp(prevScale * zoomFactor, 0.6, 6);

    const worldX = (mx - offset.x) / prevScale;
    const worldY = (my - offset.y) / prevScale;

    const nextOffsetX = mx - worldX * nextScale;
    const nextOffsetY = my - worldY * nextScale;

    setScale(nextScale);
    setOffset({ x: nextOffsetX, y: nextOffsetY });
  };

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

  const endDrag = () => {
    dragging.current = false;
    setIsInteracting(false);
  };

  const handleMouseUp = (e: React.MouseEvent<HTMLCanvasElement>) => {
    endDrag();
    if (inspectMode === 'click' && dragDistance.current < 6) {
      pickAt(e.clientX, e.clientY);
    }
  };

  const handleClick = (e: React.MouseEvent<HTMLCanvasElement>) => {
    if (inspectMode !== 'click') return;
    if (dragDistance.current < 6) {
      pickAt(e.clientX, e.clientY);
    }
  };

  const handleHover = (e: React.MouseEvent<HTMLCanvasElement>) => {
    if (inspectMode !== 'hover' || dragging.current) return;
    pendingHover.current = { x: e.clientX, y: e.clientY };
    if (hoverRaf.current !== null) return;
    hoverRaf.current = requestAnimationFrame(() => {
      hoverRaf.current = null;
      if (pendingHover.current) {
        pickAt(pendingHover.current.x, pendingHover.current.y);
        pendingHover.current = null;
      }
    });
  };

  return (
    <div ref={containerRef} className="w-full h-full bg-black relative">
      <canvas
        ref={canvasRef}
        className="w-full h-full cursor-grab active:cursor-grabbing"
        onWheel={handleWheel}
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
