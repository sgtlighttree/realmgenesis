import React from 'react';
import {
  Globe, Map as MapIcon, Hexagon, Waves, Route as RouteIcon, Sun,
  Flag, Pencil, ArrowUp, ArrowDown, Minus, Leaf, RefreshCw, Grid3x3,
  ChevronDown,
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
export const PANEL = 'border border-edge-subtle bg-surface';
const HEADER = 'flex items-center gap-2 px-3 py-2 border-b border-edge-subtle';
export const FOCUS = 'focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-brand/70 focus-visible:ring-offset-1 focus-visible:ring-offset-surface';
const CHIP = `text-[10px] px-2 py-1 border border-edge bg-surface-raised text-ink-soft whitespace-nowrap transition-colors hover:border-edge-strong hover:text-ink-strong ${FOCUS}`;
const CHIP_ON = `text-[10px] px-2 py-1 border border-brand bg-brand-strong text-ink-strong whitespace-nowrap ${FOCUS}`;

/** A panel wrapper with a consistent, quiet header — the title carries it; no
 *  eyebrow tag (bucket identity is conveyed by position, not a repeated label).
 *
 *  `collapsible` restores the affordance the classic Legend and MiniMap had.
 *  An earlier pass argued docking made it unnecessary (the rail scrolls), which
 *  was wrong: the Read rail holds three tall cards, and users want to fold the
 *  reference ones away to see the globe. Implemented ONCE here rather than
 *  per-component, which is the whole point of the shared chrome.
 *
 *  Collapsing UNMOUNTS the body rather than hiding it — `MiniMapCanvas` runs a
 *  d3 redraw on every world change, and CSS-hiding would keep that work alive. */
export const Panel: React.FC<{
  title?: string;
  right?: React.ReactNode;
  className?: string;
  bodyClassName?: string;
  collapsible?: boolean;
  defaultCollapsed?: boolean;
  children?: React.ReactNode;
}> = ({
  title, right, className = '', bodyClassName = '',
  collapsible = false, defaultCollapsed = false, children,
}) => {
  const [collapsed, setCollapsed] = React.useState(defaultCollapsed);
  const canCollapse = collapsible && !!title;
  return (
    <section className={`flex flex-col min-h-0 ${PANEL} ${className}`}>
      {(title || right) && (
        <header className={HEADER}>
          {title && (
            canCollapse ? (
              <button
                onClick={() => { setCollapsed(v => !v); }}
                aria-expanded={!collapsed}
                className={`flex flex-1 items-center gap-1.5 text-left text-xs font-semibold text-ink-bright tracking-tight hover:text-ink-strong ${FOCUS}`}
              >
                <ChevronDown
                  size={12}
                  className={`shrink-0 text-ink-muted transition-transform duration-150 ${collapsed ? '-rotate-90' : ''}`}
                />
                {title}
              </button>
            ) : (
              <h2 className="text-xs font-semibold text-ink-bright tracking-tight">{title}</h2>
            )
          )}
          {right && <div className="ml-auto flex items-center gap-1">{right}</div>}
        </header>
      )}
      {!collapsed && (
        <div className={`min-h-0 overflow-auto p-3 ${bodyClassName}`}>{children}</div>
      )}
    </section>
  );
};

/** Consistent slider — single accent color (kills the 5-thumb-color tell). */
const FauxSlider: React.FC<{ label: string; value: string }> = ({ label, value }) => (
  <label className="block">
    <div className="flex justify-between text-[11px] text-ink-muted mb-1">
      <span>{label}</span>
      <span className="text-ink-soft tabular-nums">{value}</span>
    </div>
    <input type="range" defaultValue={50} readOnly
      className="w-full h-1 bg-surface-hover appearance-none cursor-pointer accent-brand-soft" />
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
    <button className={`mt-3 flex items-center justify-center gap-2 w-full py-2.5 bg-brand-strong hover:bg-brand text-ink-strong text-sm font-semibold transition-colors ${FOCUS}`}>
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
    <div className="inline-flex overflow-hidden border border-edge">
      {RENDER_MODES.map((m, i) => (
        <span key={m.label}
          className={`inline-flex items-center gap-1 px-2.5 py-1.5 text-[11px] ${i === 0 ? 'bg-brand-strong text-ink-strong' : 'bg-surface-raised text-ink-muted'}`}>
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
export interface ReadCard {
  key: string;
  title: string;
  node: React.ReactNode;
  /** Reference cards (Biomes, 2D Projection) fold away; the Inspector does not. */
  collapsible?: boolean;
}

const BIOMES = [
  ['#3a4a9e', 'Deep Ocean'], ['#2f7fc4', 'Ocean'], ['#e8e8ee', 'Ice Cap'],
  ['#9fb0b6', 'Tundra'], ['#e0b552', 'Hot Desert'], ['#7fa05a', 'Steppe'],
  ['#3f7d3a', 'Temperate Forest'], ['#2f6d4a', 'Boreal Forest'],
];

export const READ_CARDS: ReadCard[] = [
  {
    key: 'inspector', title: 'Inspector',
    node: <p className="text-[11px] text-ink-muted italic py-6 text-center">Click a cell to inspect it.</p>,
  },
  {
    key: 'biomes', title: 'Biomes',
    node: (
      <ul className="space-y-1.5">
        {BIOMES.map(([c, name]) => (
          <li key={name} className="flex items-center gap-2 text-[11px] text-ink-soft">
            <span className="w-3 h-3 " style={{ background: c }} />{name}
          </li>
        ))}
      </ul>
    ),
  },
  {
    key: 'minimap', title: '2D Projection',
    node: <div className="aspect-[2/1] bg-gradient-to-b from-blue-950 to-slate-800 border border-edge-subtle" />,
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
      <span className="absolute inset-0 grid place-items-center text-[11px] font-mono text-ink-strong/40">
        placeholder globe
      </span>
    </div>
  </div>
);
