import React, { useState } from 'react';
import { BiomeType } from '../types';
import { BIOME_COLORS } from '../utils/colors';
import { ChevronDown, ChevronUp } from 'lucide-react';

/**
 * The swatch list on its own, with no chrome and no collapse affordance — the
 * shell's Read rail supplies both via `Panel`. `BiomeLegend` below composes it
 * with the floating-over-the-globe chrome the classic layout needs.
 */
export const BiomeLegendList: React.FC<{ columns?: 1 | 2 }> = ({ columns = 1 }) => (
  <div className={columns === 2 ? 'grid grid-cols-2 gap-x-3 gap-y-1' : 'space-y-1'}>
    {(Object.keys(BIOME_COLORS) as BiomeType[]).map((biome) => (
      // `min-w-0` is load-bearing: a flex item defaults to `min-width: auto`,
      // so the label could not shrink below its own content and pushed the grid
      // track wider than the 256px Read rail — which is what produced the
      // horizontal scrollbar on "Temperate Rainforest".
      <div key={biome} className="flex items-center gap-2 min-w-0">
        <div
          className="w-3 h-3 shrink-0 border border-white/10"
          style={{ backgroundColor: BIOME_COLORS[biome] }}
        />
        {/* Wraps rather than truncates. A legend exists to NAME things, so
            "Temperate Rainfo…" defeats the point; two lines cost a few pixels
            of height in a rail that already scrolls. */}
        <span className="text-[10px] text-ink-soft leading-tight">{biome}</span>
      </div>
    ))}
  </div>
);

export const BiomeLegend: React.FC = () => {
  const [isCollapsed, setIsCollapsed] = useState(false);
  return (
    <div className="absolute bottom-4 left-4 bg-surface/80 backdrop-blur border border-edge shadow-xl z-overlay transition-all duration-300">
      <button
        onClick={() => setIsCollapsed(!isCollapsed)}
        className="w-full flex items-center justify-between px-3 py-2 text-xs font-bold text-ink-muted uppercase hover:text-ink-strong transition-colors"
      >
        <span>Biomes</span>
        {isCollapsed ? <ChevronUp size={12}/> : <ChevronDown size={12}/>}
      </button>
      {!isCollapsed && (
        <div className="px-3 pb-3 w-48 max-h-[250px] overflow-y-auto">
          <BiomeLegendList />
        </div>
      )}
    </div>
  );
};
