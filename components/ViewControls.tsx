import React from 'react';
import {
  Globe, Satellite, Mountain, Eye, Thermometer, Droplets, Layers, Flag,
  Landmark, Palette, Church, Users, Grid, Waves, Route, Sun, LineChart,
  Pause, Play, Wind, Circle, Hexagon,
} from 'lucide-react';

import { ViewMode, DisplayMode, LabelVisibility } from '../types';
import Select, { SelectOption } from './Select';
import { MAP_STYLES, MapStyleId } from '../utils/mapStyle';

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
  /**
   * A3 style axis. Sits beside viewMode because it IS its sibling: viewMode
   * decides what the map shows, mapStyleId decides how it is drawn.
   */
  mapStyleId: MapStyleId;
  setMapStyleId: (id: MapStyleId) => void;
  displayMode: DisplayMode;
  setDisplayMode: (m: DisplayMode) => void;
  showGrid: boolean;
  setShowGrid: (b: boolean) => void;
  smoothGlobe: boolean;
  setSmoothGlobe: (b: boolean) => void;
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
  showCellEdges: boolean;
  setShowCellEdges: (b: boolean) => void;
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
export const MAP_STYLE_OPTIONS: SelectOption<MapStyleId>[] =
  Object.values(MAP_STYLES).map(st => ({ value: st.id, label: st.name }));

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

/**
 * Faction borders as a LayerToggle, so the strip can offer it as a chip.
 *
 * Deliberately NOT inside `buildLayerToggles`: that list feeds the Sys tab's
 * stacked rows as well as the strip, and the Sys tab already renders borders in
 * the Map Overlays group (`OVERLAY_KEYS`). Adding it there would put the same
 * checkbox on screen twice. The strip has no such group — its chips are the only
 * quick access to any overlay — so it composes this one in on its own.
 *
 * The state also differs in kind: every `buildLayerToggles` entry is its own
 * boolean prop, while borders is one key of the `labelVisibility` record, hence
 * the functional updater rather than a bare setter.
 */
export const buildBordersToggle = (p: ViewControlsProps): LayerToggle => ({
  key: 'borders',
  label: 'Faction Borders',
  icon: Flag, // matches the Map Overlays group header in the Sys tab
  checked: p.labelVisibility.borders,
  onChange: (b) => { p.setLabelVisibility(prev => ({ ...prev, borders: b })); },
  accent: 'text-brand-soft',
});

/** The overlay layers, bound to live state. Order matches the Sys tab. */
export const buildLayerToggles = (p: ViewControlsProps): LayerToggle[] => [
  { key: 'grid', label: 'Lat/Long Grid', icon: Grid, checked: p.showGrid, onChange: p.setShowGrid, accent: 'text-brand-soft' },
  { key: 'smooth', label: 'Smooth Globe', icon: Circle, checked: p.smoothGlobe, onChange: p.setSmoothGlobe, accent: 'text-brand-soft' },
  { key: 'rivers', label: 'River Network', icon: Waves, checked: p.showRivers, onChange: p.setShowRivers, accent: 'text-brand-soft' },
  { key: 'routes', label: 'Roads & Routes', icon: Route, checked: p.showRoutes, onChange: p.setShowRoutes, accent: 'text-warn-soft' },
  { key: 'hillshade', label: 'Hillshading', icon: Sun, checked: p.showHillshade, onChange: p.setShowHillshade, accent: 'text-brand-soft' },
  { key: 'contours', label: 'Contour Lines', icon: LineChart, checked: p.showContours, onChange: p.setShowContours, accent: 'text-brand-soft' },
  { key: 'currents', label: 'Ocean Currents', icon: Wind, checked: p.showCurrents, onChange: p.setShowCurrents, accent: 'text-brand-soft' },
  { key: 'celledges', label: 'Cell Edges', icon: Hexagon, checked: p.showCellEdges, onChange: p.setShowCellEdges, accent: 'text-brand-soft' },
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

/**
 * Compact chip form of a layer toggle — the strip shape.
 *
 * `compact` drops the label to icon-only. The name still reaches the user two
 * ways: `title` (hover tooltip) and `aria-label` (screen readers, which cannot
 * read a tooltip). Both are set in BOTH modes, so the accessible name never
 * depends on available width.
 */
const LayerChip: React.FC<{ toggle: LayerToggle; compact?: boolean }> = ({ toggle, compact = false }) => {
  const Icon = toggle.icon;
  return (
    <button
      onClick={() => { toggle.onChange(!toggle.checked); }}
      aria-pressed={toggle.checked}
      aria-label={toggle.label}
      title={toggle.label}
      className={`${CHIP_BASE} ${compact ? 'px-1.5' : ''} ${toggle.checked
        ? 'border-brand bg-brand-strong text-ink-strong'
        : 'border-edge bg-surface-raised text-ink-muted hover:border-edge-strong hover:text-ink-strong'}`}
    >
      <Icon size={11} />{compact ? null : <> {toggle.label}</>}
    </button>
  );
};

/**
 * Collapses the strip to icon-only chips when the full-label form would not fit
 * on ONE row, which is the constraint: the strip must never wrap.
 *
 * The decision is content-driven, not a hard breakpoint, because the toggle list
 * grows (it gained Cell Edges this session) and any pixel threshold would rot
 * the next time one is added. A hidden mirror of the full-label row is measured
 * against the live container: the mirror never changes with `compact`, so the
 * measurement cannot oscillate the way `scrollWidth` on the live row would.
 */
function useCompactStrip(): {
  containerRef: React.RefObject<HTMLDivElement | null>;
  chipsRef: React.RefObject<HTMLDivElement | null>;
  mirrorRef: React.RefObject<HTMLDivElement | null>;
  compact: boolean;
} {
  const containerRef = React.useRef<HTMLDivElement | null>(null);
  const chipsRef = React.useRef<HTMLDivElement | null>(null);
  const mirrorRef = React.useRef<HTMLDivElement | null>(null);
  const [compact, setCompact] = React.useState(false);

  React.useLayoutEffect(() => {
    const container = containerRef.current;
    if (!container) return;

    const measure = (): void => {
      const chips = chipsRef.current;
      const mirror = mirrorRef.current;
      if (!chips || !mirror) return;

      // Sum the strip's other children explicitly — render-mode buttons,
      // rotation, the view-layer select. Deriving this by subtracting the chip
      // row from container.scrollWidth does NOT work: the mirror is absolutely
      // positioned but still counts toward scrollWidth, which inflated the
      // result enough to force compact mode at any width.
      const gap = parseFloat(getComputedStyle(container).columnGap) || 0;
      let others = 0;
      let laidOut = 0;
      for (const el of Array.from(container.children) as HTMLElement[]) {
        if (el === mirror) continue; // out of flow — never occupies a row slot
        laidOut++;
        if (el === chips) continue;
        others += el.offsetWidth;
      }

      // What one row would cost with every label shown. The mirror never
      // reflows with `compact`, so this total is invariant and cannot oscillate.
      const needed = others + mirror.scrollWidth + gap * Math.max(0, laidOut - 1);
      setCompact(needed > container.clientWidth + 1); // 1px absorbs rounding
    };

    measure();
    const ro = new ResizeObserver(measure);
    ro.observe(container);
    if (mirrorRef.current) ro.observe(mirrorRef.current);
    return () => { ro.disconnect(); };
  }, []);

  return { containerRef, chipsRef, mirrorRef, compact };
}

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
export const ViewStrip: React.FC<ViewControlsProps & { rotation?: RotationControl }> = ({ rotation, ...p }) => {
  // Borders sits next to Roads & Routes rather than at the end: both are
  // cell-bound line overlays, and a trailing chip is the first thing lost when
  // the row starts scrolling. The strip order therefore diverges from the Sys
  // tab's on purpose — see buildBordersToggle for why it is composed in here.
  const layers = buildLayerToggles(p);
  const routeAt = layers.findIndex(t => t.key === 'routes');
  const at = routeAt < 0 ? layers.length : routeAt + 1;
  const toggles = [...layers.slice(0, at), buildBordersToggle(p), ...layers.slice(at)];
  const { containerRef, chipsRef, mirrorRef, compact } = useCompactStrip();

  // `flex-1` on the root is load-bearing for the compact measurement, not
  // cosmetic: without it the root sizes to its own CONTENT, so clientWidth
  // reported the content width rather than the space available and the
  // full-label branch could never win. `min-w-0` keeps it able to shrink.
  return (
  <div ref={containerRef} className="relative flex flex-1 flex-nowrap items-center gap-2 min-w-0 overflow-hidden">
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

    {/* Style sits next to the view layer because it is the SIBLING axis, not a
        setting: viewMode picks what the map shows, style picks how it is drawn.
        It lived in the Sys tab under Auto-Update and was reported as missing —
        a display control filed among system options is a hidden one. Only
        rendered for the 2D modes: the 3D globe is deliberately out of scope for
        styles, so offering the control there would promise something it cannot
        do. */}
    {p.displayMode !== 'globe' && (
      <Select
        value={p.mapStyleId}
        options={MAP_STYLE_OPTIONS}
        onChange={p.setMapStyleId}
        label="Map style"
        className="shrink-0"
        triggerClassName="min-w-[6.5rem] justify-between"
      />
    )}

    {/* Hidden full-label mirror: the yardstick useCompactStrip measures against.
        Never reflows with `compact`, so the measurement cannot oscillate.
        aria-hidden + inert keeps it out of the a11y tree and tab order. */}
    <div
      ref={mirrorRef}
      aria-hidden="true"
      inert
      className="pointer-events-none invisible absolute left-0 top-0 -z-10 flex flex-nowrap items-center gap-1"
    >
      {toggles.map(t => <LayerChip key={t.key} toggle={t} />)}
    </div>

    {/* Below ~500px of strip width even icon-only chips exceed the row (the wide
        shell runs down to 768px viewport). The constraint is ONE row, so the
        chips scroll horizontally rather than wrapping or being clipped
        unreachable. overflow-x-auto shows no scrollbar until it is needed. */}
    <div ref={chipsRef} className="flex flex-nowrap items-center gap-1 min-w-0 overflow-x-auto">
      {toggles.map(t => <LayerChip key={t.key} toggle={t} compact={compact} />)}
    </div>
  </div>
  );
};
