import React from 'react';
import {
  Globe, Map as MapIcon, Hexagon, Waves, Route as RouteIcon, Sun,
  Flag, Pencil, ArrowUp, ArrowDown, Minus, Leaf, RefreshCw, Grid3x3,
} from 'lucide-react';

/**
 * shellKit — placement-agnostic content stubs + shared chrome for the F1
 * "design debug" prototype (reachable via ?shell=1).
 *
 * Nothing in here positions itself. Panels fill whatever container the shell
 * mounts them in. This mirrors the real F1 contract: App composes panels, the
 * shell only positions them. These are throwaway STUBS — no real state.
 */

export type Bucket = 'make' | 'view' | 'do' | 'read';

/** What every shell receives. App composes these; the shell only positions them. */
export interface ShellProps {
  make: React.ReactNode;
  view: React.ReactNode;
  read: ReadCard[];
  doTools: React.ReactNode;
  canvas: React.ReactNode;
  editing: boolean;
  onSetEditing: (v: boolean) => void;
}

export const BUCKET: Record<Bucket, { label: string }> = {
  make: { label: 'Make' },
  view: { label: 'View' },
  do: { label: 'Do' },
  read: { label: 'Read' },
};

// One chrome definition — the single source the "consistency cleanup" buys us.
// Change radius/border/fill here and every panel moves together. Solid, not
// glass: docked panels don't need to see through, and default blur is a tell.
export const PANEL = 'border border-gray-800 bg-gray-900';
const HEADER = 'flex items-center gap-2 px-3 py-2 border-b border-gray-800';
export const FOCUS = 'focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-blue-500/70 focus-visible:ring-offset-1 focus-visible:ring-offset-gray-900';
const CHIP = `text-[10px] px-2 py-1 border border-gray-700 bg-gray-800 text-gray-300 whitespace-nowrap transition-colors hover:border-gray-600 hover:text-white ${FOCUS}`;
const CHIP_ON = `text-[10px] px-2 py-1 border border-blue-500 bg-blue-600 text-white whitespace-nowrap ${FOCUS}`;

/** A panel wrapper with a consistent, quiet header — the title carries it; no
 *  eyebrow tag (bucket identity is conveyed by position, not a repeated label). */
export const Panel: React.FC<{
  title?: string;
  right?: React.ReactNode;
  className?: string;
  bodyClassName?: string;
  children?: React.ReactNode;
}> = ({ title, right, className = '', bodyClassName = '', children }) => (
  <section className={`flex flex-col min-h-0 ${PANEL} ${className}`}>
    {(title || right) && (
      <header className={HEADER}>
        {title && <h2 className="text-xs font-semibold text-gray-100 tracking-tight">{title}</h2>}
        {right && <div className="ml-auto flex items-center gap-1">{right}</div>}
      </header>
    )}
    <div className={`min-h-0 overflow-auto p-3 ${bodyClassName}`}>{children}</div>
  </section>
);

/** Consistent slider — single accent color (kills the 5-thumb-color tell). */
const FauxSlider: React.FC<{ label: string; value: string }> = ({ label, value }) => (
  <label className="block">
    <div className="flex justify-between text-[11px] text-gray-400 mb-1">
      <span>{label}</span>
      <span className="text-gray-300 tabular-nums">{value}</span>
    </div>
    <input type="range" defaultValue={50} readOnly
      className="w-full h-1 bg-gray-700 appearance-none cursor-pointer accent-blue-400" />
  </label>
);

const TabRow: React.FC<{ tabs: string[]; active?: string }> = ({ tabs, active }) => (
  <div className="flex gap-1 mb-3">
    {tabs.map(t => (
      <span key={t} className={t === active ? CHIP_ON : CHIP}>{t}</span>
    ))}
  </div>
);

/* ------------------------------------------------------------------ *
 *  STUB PANELS — one per bucket. Dumb content, no logic.
 * ------------------------------------------------------------------ */

export const MakePanel: React.FC = () => (
  <div className="flex flex-col h-full">
    <TabRow tabs={['SYS', 'GEO', 'CLIM', 'CIV', 'EXP']} active="GEO" />
    <div className="flex flex-col gap-3 flex-1">
      <FauxSlider label="Sea Level" value="55%" />
      <FauxSlider label="Planet Radius" value="6371 km" />
      <FauxSlider label="Tectonic Plates" value="12" />
      <FauxSlider label="Terrain Roughness" value="50%" />
      <FauxSlider label="Detail Octaves" value="3" />
    </div>
    <button className={`mt-3 flex items-center justify-center gap-2 w-full py-2.5 bg-blue-600 hover:bg-blue-500 text-white text-sm font-semibold transition-colors ${FOCUS}`}>
      <RefreshCw size={14} /> Generate World
    </button>
  </div>
);

const RENDER_MODES = [
  { icon: <Globe size={13} />, label: '3D' },
  { icon: <MapIcon size={13} />, label: 'Mercator' },
  { icon: <Hexagon size={13} />, label: 'Dymaxion' },
];
const LAYER_TOGGLES = [
  { icon: <Grid3x3 size={12} />, label: 'Grid' },
  { icon: <Waves size={12} />, label: 'Rivers' },
  { icon: <RouteIcon size={12} />, label: 'Roads' },
  { icon: <Sun size={12} />, label: 'Relief' },
  { icon: <Flag size={12} />, label: 'Labels' },
];

/** View panel — render mode + layer toggles. Lays out inline; the shell decides
 *  whether it's a top strip (wide) or a sheet (narrow). */
export const ViewPanel: React.FC = () => (
  <div className="flex flex-wrap items-center gap-2">
    <div className="inline-flex overflow-hidden border border-gray-700">
      {RENDER_MODES.map((m, i) => (
        <span key={m.label}
          className={`inline-flex items-center gap-1 px-2.5 py-1.5 text-[11px] ${i === 0 ? 'bg-blue-600 text-white' : 'bg-gray-800 text-gray-400'}`}>
          {m.icon}{m.label}
        </span>
      ))}
    </div>
    <div className="flex flex-wrap gap-1">
      {LAYER_TOGGLES.map(t => (
        <span key={t.label} className={`${CHIP} inline-flex items-center gap-1`}>{t.icon}{t.label}</span>
      ))}
    </div>
  </div>
);

const DO_TOOLS = [
  { icon: <ArrowUp size={11} />, label: 'Raise' },
  { icon: <ArrowDown size={11} />, label: 'Lower' },
  { icon: <Minus size={11} />, label: 'Flatten' },
  { icon: <Waves size={11} />, label: 'Smooth' },
  { icon: <Leaf size={11} />, label: 'Biome' },
  { icon: <Flag size={11} />, label: 'Political' },
  { icon: <Pencil size={11} />, label: 'Edit' },
];

/** Do panel — the edit toolset. Contextual; the shell only mounts it when
 *  editing is active. */
export const DoPanel: React.FC = () => (
  <div className="flex flex-wrap items-center gap-1">
    {DO_TOOLS.map((t, i) => (
      <span key={t.label}
        className={`${i === 0 ? CHIP_ON : CHIP} inline-flex items-center gap-1`}>{t.icon}{t.label}</span>
    ))}
  </div>
);

/** Read cards — each self-contained, shell decides stacking vs tabbing. */
export interface ReadCard { key: string; title: string; node: React.ReactNode; }

const BIOMES = [
  ['#3a4a9e', 'Deep Ocean'], ['#2f7fc4', 'Ocean'], ['#e8e8ee', 'Ice Cap'],
  ['#9fb0b6', 'Tundra'], ['#e0b552', 'Hot Desert'], ['#7fa05a', 'Steppe'],
  ['#3f7d3a', 'Temperate Forest'], ['#2f6d4a', 'Boreal Forest'],
];

export const READ_CARDS: ReadCard[] = [
  {
    key: 'inspector', title: 'Inspector',
    node: <p className="text-[11px] text-gray-400 italic py-6 text-center">Click a cell to inspect it.</p>,
  },
  {
    key: 'biomes', title: 'Biomes',
    node: (
      <ul className="space-y-1.5">
        {BIOMES.map(([c, name]) => (
          <li key={name} className="flex items-center gap-2 text-[11px] text-gray-300">
            <span className="w-3 h-3 " style={{ background: c }} />{name}
          </li>
        ))}
      </ul>
    ),
  },
  {
    key: 'minimap', title: '2D Projection',
    node: <div className="aspect-[2/1] bg-gradient-to-b from-blue-950 to-slate-800 border border-gray-800" />,
  },
];

/** Placeholder globe — a dumb sphere. Nothing rendered too hard. */
export const PlaceholderGlobe: React.FC = () => (
  <div className="w-full h-full grid place-items-center overflow-hidden">
    <div className="relative rounded-full"
      style={{
        width: 'min(72vmin, 640px)', aspectRatio: '1',
        background: 'radial-gradient(circle at 38% 32%, #26405e, #0c1626 72%)',
        boxShadow: 'inset -18px -22px 60px rgba(0,0,0,.55), 0 0 80px rgba(40,90,150,.18)',
      }}>
      <div className="absolute inset-[8%] rounded-full border border-white/10" />
      <div className="absolute inset-[8%_30%] rounded-full border border-white/10" />
      <span className="absolute inset-0 grid place-items-center text-[11px] font-mono text-white/40">
        placeholder globe
      </span>
    </div>
  </div>
);
