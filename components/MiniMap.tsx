import React, { useEffect, useRef, useState } from 'react';
import * as d3 from 'd3';
import { WorldData, ViewMode } from '../types';
import { buildFactionColorMap, buildCultureColorMap, buildReligionColorMap, getCellColor } from '../utils/colors';
import { ChevronDown, ChevronUp } from 'lucide-react';

interface MiniMapProps {
  world: WorldData | null;
  viewMode: ViewMode;
}

/**
 * The projection canvas alone, no chrome and no collapse — the shell's Read
 * rail supplies both via `Panel`. NOTE: the redraw is debounced but not
 * otherwise gated here, so hosts must UNMOUNT this to stop the work; hiding it
 * with CSS would keep redrawing 5k paths on every paint stroke.
 */
export const MiniMapCanvas: React.FC<MiniMapProps> = ({ world, viewMode }) => {
  const canvasRef = useRef<HTMLCanvasElement>(null);

  useEffect(() => {
    if (!world || !canvasRef.current) return;
    // Debounced: paint strokes replace `world` ~20x/sec; redrawing 5k paths
    // per stroke event is wasted work for an overview map
    const timer = window.setTimeout(() => {
      const canvas = canvasRef.current;
      if (!canvas) return;
      const ctx = canvas.getContext('2d');
      if (!ctx) return;
      const width = canvas.width; const height = canvas.height;
      ctx.fillStyle = '#111'; ctx.fillRect(0, 0, width, height);
      ctx.save();
      ctx.translate(width, 0); ctx.scale(-1, 1);
      const projection = d3.geoEquirectangular().fitSize([width, height], { type: "Sphere" } as d3.GeoPermissibleObjects);
      const pathGenerator = d3.geoPath(projection, ctx);
      const factionColors = buildFactionColorMap(world.civData);
      const cultureColors = buildCultureColorMap(world.cultures);
      const religionColors = buildReligionColorMap(world.religions);
      world.cells.forEach((cell, i) => {
          if (!world.geoJson || !world.geoJson.features[i]) { return; }
          const feature = world.geoJson.features[i];
          if (!feature.geometry || feature.geometry.coordinates.length === 0) return;
          const color = getCellColor(cell, viewMode, world.params.seaLevel, factionColors, cultureColors, religionColors);
          ctx.beginPath(); pathGenerator(feature);
          ctx.fillStyle = '#' + color.getHexString(); ctx.fill();
      });
      ctx.restore();
    }, 150);
    return () => { window.clearTimeout(timer); };
  }, [world, viewMode]);

  if (!world) return null;

  return (
    <canvas
      ref={canvasRef}
      width={240}
      height={120}
      className="cursor-crosshair opacity-90 hover:opacity-100 transition-opacity max-w-full"
    />
  );
};

const MiniMap: React.FC<MiniMapProps> = ({ world, viewMode }) => {
  const [isCollapsed, setIsCollapsed] = useState(false);

  if (!world) return null;

  return (
    <div className="absolute bottom-4 right-4 bg-black/80 border border-edge shadow-2xl overflow-hidden z-overlay transition-all duration-300">
      <button
        onClick={() => { setIsCollapsed(!isCollapsed); }}
        className="w-full flex items-center justify-between px-3 py-1.5 text-[10px] text-ink-muted font-bold uppercase hover:text-ink-strong transition-colors"
      >
        <span>2D Projection</span>
        {isCollapsed ? <ChevronUp size={12}/> : <ChevronDown size={12}/>}
      </button>
      {!isCollapsed && (
        <div className="p-1">
          <MiniMapCanvas world={world} viewMode={viewMode} />
        </div>
      )}
    </div>
  );
};

export default MiniMap;
