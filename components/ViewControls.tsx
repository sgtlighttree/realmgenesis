import React from 'react';
import {
  Globe, Satellite, Mountain, Eye, Thermometer, Droplets, Layers, Flag,
  Landmark, Palette, Church, Users, Grid, Waves, Route, Sun, LineChart,
  Pause, Play, Wind,
} from 'lucide-react';

import { ViewMode, DisplayMode, LabelVisibility } from '../types';
import Select, { SelectOption } from './Select';

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
  showCurrents: boolean;
  setShowCurrents: (b: boolean) => void;
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

/** Same list as VIEW_LAYERS, shaped for the themed Select. */
export const VIEW_LAYER_OPTIONS: SelectOption<ViewMode>[] =
  VIEW_LAYERS.map(l => ({ value: l.mode, label: l.label }));

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
  { key: 'grid', label: 'Lat/Long Grid', icon: Grid, checked: p.showGrid, onChange: p.setShowGrid, accent: 'text-brand-soft' },
  { key: 'rivers', label: 'River Network', icon: Waves, checked: p.showRivers, onChange: p.setShowRivers, accent: 'text-brand-soft' },
  { key: 'routes', label: 'Roads & Routes', icon: Route, checked: p.showRoutes, onChange: p.setShowRoutes, accent: 'text-warn-soft' },
  { key: 'hillshade', label: 'Hillshading', icon: Sun, checked: p.showHillshade, onChange: p.setShowHillshade, accent: 'text-brand-soft' },
  { key: 'contours', label: 'Contour Lines', icon: LineChart, checked: p.showContours, onChange: p.setShowContours, accent: 'text-brand-soft' },
  { key: 'currents', label: 'Ocean Currents', icon: Wind, checked: p.showCurrents, onChange: p.setShowCurrents, accent: 'text-brand-soft' },
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
 ? 'bg-brand-strong text-ink-strong border-brand border-b-2'
 : 'bg-surface-raised text-ink-muted border-edge hover:bg-surface-hover hover:text-ink-strong'
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
 ? 'bg-brand-strong text-ink-strong border-brand border-b-2'
 : 'bg-surface-raised text-ink-muted border-edge hover:bg-surface-hover hover:text-ink-strong'
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
    <div className={`flex items-center justify-between text-xs text-ink-muted ${className}`}>
      <div className="flex items-center gap-2">
        <Icon size={12} className={toggle.checked ? toggle.accent : 'text-ink-faint'} />
        <label>{toggle.label}</label>
      </div>
      <input
        type="checkbox"
        checked={toggle.checked}
        onChange={(e) => { toggle.onChange(e.target.checked); }}
        className="bg-surface-hover"
      />
    </div>
  );
};

/** The Map Overlays group (icon header + indented per-kind checkboxes). */
export const OverlayToggles: React.FC<Pick<ViewControlsProps, 'labelVisibility' | 'setLabelVisibility'>> = ({
  labelVisibility, setLabelVisibility,
}) => (
  <div className="pt-2">
    <div className="flex items-center gap-2 text-xs text-ink-muted mb-1">
      <Flag size={12} className={(labelVisibility.borders || labelVisibility.factions) ? 'text-brand-soft' : 'text-ink-faint'} />
      <label className="font-medium">Map Overlays</label>
    </div>
    <div className="ml-5 space-y-1">
      {OVERLAY_KEYS.map(([key, label]) => (
        <div key={key} className="flex items-center justify-between text-xs text-ink-muted">
          <label>{label}</label>
          <input
            type="checkbox"
            checked={labelVisibility[key]}
            onChange={(e) => { setLabelVisibility(prev => ({ ...prev, [key]: e.target.checked })); }}
            className="bg-surface-hover"
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

const CHIP_BASE = 'inline-flex items-center gap-1 text-[10px] px-2 py-1 border whitespace-nowrap transition-colors focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-brand/70';

/** Compact chip form of a layer toggle — the strip shape. */
const LayerChip: React.FC<{ toggle: LayerToggle }> = ({ toggle }) => {
  const Icon = toggle.icon;
  return (
    <button
      onClick={() => { toggle.onChange(!toggle.checked); }}
      aria-pressed={toggle.checked}
      title={toggle.label}
      className={`${CHIP_BASE} ${toggle.checked
        ? 'border-brand bg-brand-strong text-ink-strong'
        : 'border-edge bg-surface-raised text-ink-muted hover:border-edge-strong hover:text-ink-strong'}`}
    >
      <Icon size={11} /> {toggle.label}
    </button>
  );
};

/**
 * Rotation control for the strip. Scoped to `ViewStrip` rather than added to
 * `ViewControlsProps` because the Sys tab has no use for it — the globe's own
 * canvas overlay serves the folds that still carry one.
 */
export interface RotationControl {
  paused: boolean;
  onToggle: () => void;
  /** The Dymaxion overlay editor force-pauses rotation and owns it while open. */
  disabled?: boolean;
}

/**
 * ViewStrip — the horizontal composition for the wide shell's top strip.
 * Render mode as a segmented control, the 12 view layers as a select (a
 * 12-button grid does not fit a strip), and the layer toggles as chips.
 */
export const ViewStrip: React.FC<ViewControlsProps & { rotation?: RotationControl }> = ({ rotation, ...p }) => (
  <div className="flex flex-wrap items-center gap-2 min-w-0">
    <div className="inline-flex overflow-hidden border border-edge shrink-0">
      {DISPLAY_MODES.map(m => (
        <button
          key={m.mode}
          onClick={() => { p.setDisplayMode(m.mode); }}
          aria-pressed={p.displayMode === m.mode}
          className={`px-2.5 py-1.5 text-[11px] transition-colors ${p.displayMode === m.mode
            ? 'bg-brand-strong text-ink-strong'
            : 'bg-surface-raised text-ink-muted hover:bg-surface-hover hover:text-ink-strong'}`}
        >
          {m.short}
        </button>
      ))}
    </div>

    {/* Rotation sits immediately after render mode because it modifies that
        choice — it is a property of how the globe is PRESENTED, not something
        you do to the world, so it does not belong with Edit. Leading the strip
        with it would put a secondary control ahead of the primary one.
        It previously lived as a canvas overlay anchored to the viewer's own
        corner, which the wide shell's left-shifted canvas clipped out of view
        entirely — the reason it appeared to vanish. */}
    {rotation && (
      <button
        onClick={rotation.onToggle}
        disabled={rotation.disabled}
        aria-pressed={rotation.paused}
        aria-label={rotation.paused ? 'Resume globe rotation' : 'Pause globe rotation'}
        title={rotation.paused ? 'Resume rotation' : 'Pause rotation'}
        className={`shrink-0 inline-flex items-center justify-center px-2.5 py-1.5 border transition-colors ${
          rotation.disabled
            ? 'bg-surface-raised text-ink-faint border-edge opacity-40 cursor-not-allowed'
            : rotation.paused
              ? 'bg-surface-raised text-ink-strong border-edge-strong hover:bg-surface-hover'
              : 'bg-surface-raised text-ink-muted border-edge hover:bg-surface-hover hover:text-ink-strong'
        }`}
      >
        {rotation.paused ? <Play size={12} /> : <Pause size={12} />}
      </button>
    )}

    <Select
      value={p.viewMode}
      options={VIEW_LAYER_OPTIONS}
      onChange={p.setViewMode}
      label="View layer"
      className="shrink-0"
      triggerClassName="min-w-[7.5rem] justify-between"
    />

    <div className="flex flex-wrap items-center gap-1 min-w-0">
      {buildLayerToggles(p).map(t => <LayerChip key={t.key} toggle={t} />)}
    </div>
  </div>
);
