import React, { useState } from 'react';
import { BiomeType } from '../types';
import { BIOME_COLORS } from '../utils/colors';
import { ChevronDown, ChevronUp } from 'lucide-react';

export const BiomeLegend: React.FC = () => {
  const [isCollapsed, setIsCollapsed] = useState(false);
  return (
    <div className="absolute bottom-4 left-4 bg-gray-900/80 backdrop-blur border border-gray-700 rounded-lg shadow-xl z-10 transition-all duration-300">
      <button 
        onClick={() => setIsCollapsed(!isCollapsed)}
        className="w-full flex items-center justify-between px-3 py-2 text-xs font-bold text-gray-400 uppercase hover:text-white transition-colors"
      >
        <span>Biomes</span>
        {isCollapsed ? <ChevronUp size={12}/> : <ChevronDown size={12}/>}
      </button>
      {!isCollapsed && (
        <div className="px-3 pb-3 space-y-1 max-h-[250px] overflow-y-auto w-48">
          {(Object.keys(BIOME_COLORS) as BiomeType[]).map((biome) => (
            <div key={biome} className="flex items-center gap-2">
              <div 
                className="w-3 h-3 rounded-full border border-white/10" 
                style={{ backgroundColor: BIOME_COLORS[biome] }}
              />
              <span className="text-[10px] text-gray-300 whitespace-nowrap">{biome}</span>
            </div>
          ))}
        </div>
      )}
    </div>
  );
};