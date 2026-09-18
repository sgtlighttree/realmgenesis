import React, { useState } from 'react';
import { BiomeType } from '../types';
import { BIOME_COLORS, Legend } from '../utils/colors';
import { ChevronDown, ChevronUp } from 'lucide-react';

/** One colour swatch + name, shared by the swatch legends and ramp OOB rows. */
const SwatchRow: React.FC<{ color: string; label: string }> = ({ color, label }) => (
  // `min-w-0` lets the label shrink inside the 2-col grid instead of forcing a
  // horizontal scrollbar (see BiomeLegendList note below).
  <div className="flex items-center gap-2 min-w-0">
    <div className="w-3 h-3 shrink-0 border border-white/10" style={{ backgroundColor: color }} />
    <span className="text-[10px] text-ink-soft leading-tight">{label}</span>
  </div>
);

/**
 * F5: the layer-aware legend. Renders whatever `legendFor(viewMode, world)`
 * returned — a swatch list for categorical layers, or a sampled gradient bar
 * with end labels, a sea-level tick, and out-of-band swatches for ramp layers.
 * The panel is only mounted when `legendFor` returns non-null, so satellite and
 * plates carry no legend at all.
 */
export const LayerLegend: React.FC<{ legend: Legend; columns?: 1 | 2 }> = ({ legend, columns = 1 }) => {
  if (legend.kind === 'swatches') {
    if (legend.items.length === 0) {
      return <p className="text-[10px] text-ink-muted italic leading-tight">No {legend.title.toLowerCase()} yet.</p>;
    }
    return (
      <div className={columns === 2 ? 'grid grid-cols-2 gap-x-3 gap-y-1' : 'space-y-1'}>
        {legend.items.map((it) => (
          <SwatchRow key={it.label} color={it.color} label={it.label} />
        ))}
      </div>
    );
  }

  const gradient = `linear-gradient(to right, ${legend.stops
    .map((s) => `${s.color} ${(s.pos * 100).toFixed(1)}%`)
    .join(', ')})`;
  const tickPct = legend.tick ? `${(legend.tick.pos * 100).toFixed(1)}%` : undefined;
  return (
    <div className="space-y-1.5">
      <div className="relative">
        <div className="h-3 w-full border border-white/10" style={{ background: gradient }} />
        {legend.tick && (
          <div
            className="absolute top-0 bottom-0 w-px bg-white/80"
            style={{ left: tickPct }}
            title={legend.tick.label}
          />
        )}
      </div>
      <div className="flex justify-between text-[10px] text-ink-soft leading-tight">
        <span>{legend.minLabel}</span>
        {legend.tick && <span className="text-ink-muted">{legend.tick.label}</span>}
        <span>{legend.maxLabel}</span>
      </div>
      {legend.oob && legend.oob.length > 0 && (
        <div className="space-y-1 pt-1 border-t border-edge/60">
          {legend.oob.map((s) => (
            <SwatchRow key={s.label} color={s.color} label={s.label} />
          ))}
        </div>
      )}
    </div>
  );
};

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
