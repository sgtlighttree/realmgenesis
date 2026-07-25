import React from 'react';
import {
  Globe, Satellite, Mountain, Eye, Thermometer, Droplets, Layers, Flag,
  Landmark, Palette, Church, Users, Grid, Waves, Route, Sun, LineChart,
} from 'lucide-react';

import { ViewMode, DisplayMode, LabelVisibility } from '../types';

/** Lucide icons take size/className; React.ElementType is too loose to accept them. */
type IconType = React.ComponentType<{ size?: number; className?: string }>;

/**
 * ViewControls — the "View" bucket's shared primitives.
 *
 * Deliberately exports PRIMITIVES plus one composed strip rather than a single
 * component with an orientation flag: the Sys tab stacks these vertically and
 * interleaves them with Make content (Render Mode, then Seed/points, then the
 * layer toggles, then View Layer), while the shell's top strip lays them out
 * inline. Layout is the part that genuinely differs per host, so each host
 * composes it; only the buttons, toggle definitions, and layer list are shared.
 */

export interface ViewControlsProps {
  viewMode: ViewMode;
  setViewMode: (m: ViewMode) => void;
  displayMode: DisplayMode;
  setDisplayMode: (m: DisplayMode) => void;
  showGrid: boolean;
  setShowGrid: (b: boolean) => void;
  showRivers: boolean;
  setShowRivers: (b: boolean) => void;
  showRoutes: boolean;
  setShowRoutes: (b: boolean) => void;
  showHillshade: boolean;
  setShowHillshade: (b: boolean) => void;
  showContours: boolean;
  setShowContours: (b: boolean) => void;
  labelVisibility: LabelVisibility;
  setLabelVisibility: React.Dispatch<React.SetStateAction<LabelVisibility>>;
}

/* ------------------------------------------------------------------ *
 *  Shared data
 * ------------------------------------------------------------------ */

export const DISPLAY_MODES: { mode: DisplayMode; label: string; short: string }[] = [
  { mode: 'globe', label: '3D Globe', short: '3D' },
  { mode: 'mercator', label: '2D Mercator', short: 'Mercator' },
  { mode: 'dymaxion', label: '2D Dymaxion', short: 'Dymaxion' },
];

export const VIEW_LAYERS: { mode: ViewMode; icon: IconType; label: string }[] = [
  { mode: 'biome', icon: Globe, label: 'Biomes' },
  { mode: 'satellite', icon: Satellite, label: 'Satellite' },
  { mode: 'height', icon: Mountain, label: 'Height' },
  { mode: 'height_bw', icon: Eye, label: 'Height BW' },
  { mode: 'temperature', icon: Thermometer, label: 'Temp' },
  { mode: 'moisture', icon: Droplets, label: 'Rain' },
  { mode: 'plates', icon: Layers, label: 'Plates' },
  { mode: 'political', icon: Flag, label: 'Borders' },
  { mode: 'province', icon: Landmark, label: 'Provinces' },
  { mode: 'culture', icon: Palette, label: 'Cultures' },
  { mode: 'religion', icon: Church, label: 'Religions' },
  { mode: 'population', icon: Users, label: 'Population' },
];

export const OVERLAY_KEYS: [keyof LabelVisibility, string][] = [
  ['borders', 'Faction Borders'],
  ['factions', 'Faction Names'],
  ['capitals', 'Capital Names'],
  ['provinces', 'Province Names'],
  ['towns', 'Town Names'],
  ['geography', 'Geographic Names'],
  ['markers', 'Markers'],
];

export interface LayerToggle {
  key: string;
  label: string;
  icon: IconType;
  checked: boolean;
  onChange: (b: boolean) => void;
  /** Tailwind text color when active — Roads is amber, the rest blue. */
  accent: string;
}

/** The five overlay layers, bound to live state. Order matches the Sys tab. */
export const buildLayerToggles = (p: ViewControlsProps): LayerToggle[] => [
  { key: 'grid', label: 'Lat/Long Grid', icon: Grid, checked: p.showGrid, onChange: p.setShowGrid, accent: 'text-blue-400' },
  { key: 'rivers', label: 'River Network', icon: Waves, checked: p.showRivers, onChange: p.setShowRivers, accent: 'text-blue-400' },
  { key: 'routes', label: 'Roads & Routes', icon: Route, checked: p.showRoutes, onChange: p.setShowRoutes, accent: 'text-amber-400' },
  { key: 'hillshade', label: 'Hillshading', icon: Sun, checked: p.showHillshade, onChange: p.setShowHillshade, accent: 'text-blue-400' },
  { key: 'contours', label: 'Contour Lines', icon: LineChart, checked: p.showContours, onChange: p.setShowContours, accent: 'text-blue-400' },
];

/* ------------------------------------------------------------------ *
 *  Primitives — markup preserved verbatim from the Sys tab so the
 *  classic panel renders byte-identically after the extraction.
 * ------------------------------------------------------------------ */

export const DisplayButton: React.FC<{
  mode: DisplayMode;
  label: string;
  displayMode: DisplayMode;
  setDisplayMode: (m: DisplayMode) => void;
}> = ({ mode, label, displayMode, setDisplayMode }) => (
  <button
    onClick={() => { setDisplayMode(mode); }}
    className={`px-2 py-1.5 text-xs transition-all flex-1 border ${
 displayMode === mode
 ? 'bg-blue-600 text-white border-blue-500 border-b-2'
 : 'bg-gray-800 text-gray-400 border-gray-700 hover:bg-gray-700 hover:text-white'
 }`}
  >
    {label}
  </button>
);

export const ViewButton: React.FC<{
  mode: ViewMode;
  icon: IconType;
  label: string;
  viewMode: ViewMode;
  setViewMode: (m: ViewMode) => void;
}> = ({ mode, icon: Icon, label, viewMode, setViewMode }) => (
  <button
    onClick={() => { setViewMode(mode); }}
    className={`flex items-center gap-2 px-2 py-1.5 text-xs transition-all flex-1 justify-center border ${
 viewMode === mode
 ? 'bg-blue-600 text-white border-blue-500 border-b-2'
 : 'bg-gray-800 text-gray-400 border-gray-700 hover:bg-gray-700 hover:text-white'
 }`}
  >
    <Icon size={14} />
    {label}
  </button>
);

/** Stacked label-left / checkbox-right row — the Sys tab shape. */
export const LayerToggleRow: React.FC<{ toggle: LayerToggle; className?: string }> = ({
  toggle, className = 'pt-2',
}) => {
  const Icon = toggle.icon;
  return (
    <div className={`flex items-center justify-between text-xs text-gray-400 ${className}`}>
      <div className="flex items-center gap-2">
        <Icon size={12} className={toggle.checked ? toggle.accent : 'text-gray-600'} />
        <label>{toggle.label}</label>
      </div>
      <input
        type="checkbox"
        checked={toggle.checked}
        onChange={(e) => { toggle.onChange(e.target.checked); }}
        className="bg-gray-700"
      />
    </div>
  );
};

/** The Map Overlays group (icon header + indented per-kind checkboxes). */
export const OverlayToggles: React.FC<Pick<ViewControlsProps, 'labelVisibility' | 'setLabelVisibility'>> = ({
  labelVisibility, setLabelVisibility,
}) => (
  <div className="pt-2">
    <div className="flex items-center gap-2 text-xs text-gray-400 mb-1">
      <Flag size={12} className={(labelVisibility.borders || labelVisibility.factions) ? 'text-blue-400' : 'text-gray-600'} />
      <label className="font-medium">Map Overlays</label>
    </div>
    <div className="ml-5 space-y-1">
      {OVERLAY_KEYS.map(([key, label]) => (
        <div key={key} className="flex items-center justify-between text-xs text-gray-400">
          <label>{label}</label>
          <input
            type="checkbox"
            checked={labelVisibility[key]}
            onChange={(e) => { setLabelVisibility(prev => ({ ...prev, [key]: e.target.checked })); }}
            className="bg-gray-700"
          />
        </div>
      ))}
    </div>
  </div>
);

/** The 12-layer grid — Sys tab shape. */
export const ViewLayerGrid: React.FC<Pick<ViewControlsProps, 'viewMode' | 'setViewMode'>> = ({
  viewMode, setViewMode,
}) => (
  <div className="grid grid-cols-2 gap-2">
    {VIEW_LAYERS.map(l => (
      <ViewButton key={l.mode} mode={l.mode} icon={l.icon} label={l.label}
        viewMode={viewMode} setViewMode={setViewMode} />
    ))}
  </div>
);

/* ------------------------------------------------------------------ *
 *  Compositions
 * ------------------------------------------------------------------ */

const CHIP_BASE = 'inline-flex items-center gap-1 text-[10px] px-2 py-1 rounded border whitespace-nowrap transition-colors focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-blue-500/70';

/** Compact chip form of a layer toggle — the strip shape. */
const LayerChip: React.FC<{ toggle: LayerToggle }> = ({ toggle }) => {
  const Icon = toggle.icon;
  return (
    <button
      onClick={() => { toggle.onChange(!toggle.checked); }}
      aria-pressed={toggle.checked}
      title={toggle.label}
      className={`${CHIP_BASE} ${toggle.checked
        ? 'border-blue-500 bg-blue-600 text-white'
        : 'border-gray-700 bg-gray-800 text-gray-400 hover:border-gray-600 hover:text-white'}`}
    >
      <Icon size={11} /> {toggle.label}
    </button>
  );
};

/**
 * ViewStrip — the horizontal composition for the wide shell's top strip.
 * Render mode as a segmented control, the 12 view layers as a select (a
 * 12-button grid does not fit a strip), and the layer toggles as chips.
 */
export const ViewStrip: React.FC<ViewControlsProps> = (p) => (
  <div className="flex flex-wrap items-center gap-2 min-w-0">
    <div className="inline-flex rounded-md overflow-hidden border border-gray-700 shrink-0">
      {DISPLAY_MODES.map(m => (
        <button
          key={m.mode}
          onClick={() => { p.setDisplayMode(m.mode); }}
          aria-pressed={p.displayMode === m.mode}
          className={`px-2.5 py-1.5 text-[11px] transition-colors ${p.displayMode === m.mode
            ? 'bg-blue-600 text-white'
            : 'bg-gray-800 text-gray-400 hover:bg-gray-700 hover:text-white'}`}
        >
          {m.short}
        </button>
      ))}
    </div>

    <select
      value={p.viewMode}
      onChange={(e) => { p.setViewMode(e.target.value as ViewMode); }}
      title="View layer"
      className="shrink-0 bg-gray-800 border border-gray-700 rounded text-[11px] text-gray-200 px-2 py-1.5 focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-blue-500/70"
    >
      {VIEW_LAYERS.map(l => (
        <option key={l.mode} value={l.mode}>{l.label}</option>
      ))}
    </select>

    <div className="flex flex-wrap items-center gap-1 min-w-0">
      {buildLayerToggles(p).map(t => <LayerChip key={t.key} toggle={t} />)}
    </div>
  </div>
);
