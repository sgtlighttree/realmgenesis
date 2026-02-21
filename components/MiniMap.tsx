import React, { useEffect, useRef, useState } from 'react';
import * as d3 from 'd3';
import { WorldData, ViewMode } from '../types';
import { getCellColor } from '../utils/colors';
import { ChevronDown, ChevronUp } from 'lucide-react';

interface MiniMapProps {
  world: WorldData | null;
  viewMode: ViewMode;
}

const MiniMap: React.FC<MiniMapProps> = ({ world, viewMode }) => {
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const [isCollapsed, setIsCollapsed] = useState(false);

  useEffect(() => {
    if (!world || !canvasRef.current || isCollapsed) return;
    const canvas = canvasRef.current;
    const ctx = canvas.getContext('2d');
    if (!ctx) return;
    const width = canvas.width; const height = canvas.height;
    ctx.fillStyle = '#111'; ctx.fillRect(0, 0, width, height);
    ctx.save();
    ctx.translate(width, 0); ctx.scale(-1, 1);
    const projection = d3.geoEquirectangular().fitSize([width, height], { type: "Sphere" } as d3.GeoPermissibleObjects);
    const pathGenerator = d3.geoPath(projection, ctx);
    world.cells.forEach((cell, i) => {
        // eslint-disable-next-line security/detect-object-injection, nosemgrep
        // codacy-disable-next-line
        if (!world.geoJson || !world.geoJson.features[i]) { return; }
        // eslint-disable-next-line security/detect-object-injection, nosemgrep
        // codacy-disable-next-line
        const feature = world.geoJson.features[i];
        if (!feature.geometry || feature.geometry.coordinates.length === 0) return;
        const color = getCellColor(cell, viewMode, world.params.seaLevel);
        ctx.beginPath(); pathGenerator(feature);
        ctx.fillStyle = '#' + color.getHexString(); ctx.fill();
    });
    ctx.restore();
  }, [world, viewMode, isCollapsed]);

  if (!world) return null;

  return (
    <div className="absolute bottom-4 right-4 bg-black/80 border border-gray-700 shadow-2xl overflow-hidden z-10 transition-all duration-300">
      <button 
        onClick={() => { setIsCollapsed(!isCollapsed); }}
        className="w-full flex items-center justify-between px-3 py-1.5 text-[10px] text-gray-400 font-bold uppercase hover:text-white transition-colors"
      >
        <span>2D Projection</span>
        {isCollapsed ? <ChevronUp size={12}/> : <ChevronDown size={12}/>}
      </button>
      {!isCollapsed && (
        <div className="p-1">
          <canvas 
            ref={canvasRef} 
            width={240} 
            height={120} 
            className="cursor-crosshair opacity-90 hover:opacity-100 transition-opacity"
          />
        </div>
      )}
    </div>
  );
};

export default MiniMap;