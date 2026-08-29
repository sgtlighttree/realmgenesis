import React from 'react';
import {
  Globe, Satellite, Mountain, Eye, Thermometer, Droplets, Layers, Flag,
  Landmark, Palette, Church, Users, Grid, Waves, Route, Sun, LineChart,
  Pause, Play, Wind, Circle, Hexagon,
} from 'lucide-react';

import { ViewMode, DisplayMode, LabelVisibility } from '../types';
import Select, { SelectOption } from './Select';
import CheckboxMenu from './CheckboxMenu';
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
  // Equirectangular first among the flat ones: it is the projection the whole
  // pipeline is expressed in (the globe texture is baked in it, and every
  // lon/lat mapping assumes it), so it is the least surprising default 2D view.
  { mode: 'equirectangular', label: '2D Equirectangular', short: 'Equirect' },
  { mode: 'mercator', label: '2D Mercator', short: 'Mercator' },
  { mode: 'winkeltripel', label: '2D Winkel Tripel', short: 'Winkel' },
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

/**
 * DISPLAY_MODES shaped for the themed Select, using the FULL label rather than
 * the segmented control's `short`. The 2D/3D prefix is the load-bearing part —
 * it is the difference between a globe and a flat map, not decoration — and a
 * menu row has the width for it where a five-way segmented control did not.
 */
export const PROJECTION_OPTIONS: SelectOption<DisplayMode>[] =
  DISPLAY_MODES.map(m => ({ value: m.mode, label: m.label }));

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

/**
 * A control with a visible caption above it.
 *
 * The caption is the control's accessible name — it is passed down as
 * `labelledBy`, and `Select`/`CheckboxMenu` drop their `aria-label` when they
 * receive one. Two accessible names on one control is worse than none: a screen
 * reader announces whichever wins, and it is not necessarily the visible text
 * the user is looking at.
 */
/**
 * Sizing shared by every dropdown field.
 *
 * `flex-1` between a floor and a cap, rather than `shrink-0` at a fixed width.
 * The wide shell starts at 768px, where a left rail, a right panel and an Edit
 * button leave the strip a couple of hundred pixels, and the three sizings that
 * do NOT work were each tried in the browser:
 *
 * - Fixed `min-w` per trigger: the row overflowed and pushed Overlays under the
 *   right panel at 1024px.
 * - `min-w-0` with no floor: at 768px the four controls squeezed to empty stubs
 *   with their captions overlapping. Legible truncation needs a floor.
 * - No cap: four controls stretched across a 2560px display.
 *
 * So they share what there is, down to a floor, and the ROW scrolls below that
 * — the same last-resort the chip row used to have. A scroll there closes any
 * open menu, which is correct: the trigger has moved out from under it.
 */
const FIELD = 'flex-1 min-w-[5rem] max-w-[11rem]';

const Field: React.FC<{
  caption: string;
  children: (captionId: string) => React.ReactNode;
  className?: string;
}> = ({ caption, children, className = '' }) => {
  const id = React.useId();
  return (
    <div className={`flex flex-col gap-0.5 ${className}`}>
      <span id={id} className="text-[9px] uppercase tracking-wider text-ink-faint">{caption}</span>
      {children(id)}
    </div>
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
 *
 * ONE row of captioned dropdowns: projection, layer, style, overlays, rotation.
 *
 * It was two rows — a five-way segmented control for the projection, three
 * selects, and a row of nine overlay chips. The forcing function was the
 * projection: it grew from three entries to five (Equirectangular and Winkel
 * Tripel landed in S27f) and a segmented control does not have five slots in a
 * strip that also has to hold everything else. Making that a dropdown and
 * leaving the chips would have left one row of dropdowns above one row of
 * chips, which is two grammars for the same job, so the chips became a
 * multi-select menu at the same time.
 *
 * What that costs, and what pays it back: the chips showed which overlays were
 * ON without being opened. `CheckboxMenu` puts the active count on its trigger
 * for exactly that reason — see its doc.
 *
 * `useCompactStrip` and its hidden measuring mirror went with the chips. The
 * whole apparatus existed to decide whether nine labelled chips fitted on one
 * row; five dropdowns always do.
 */
export const ViewStrip: React.FC<ViewControlsProps & { rotation?: RotationControl }> = ({ rotation, ...p }) => {
  const rotationId = React.useId();
  // Borders sits next to Roads & Routes rather than at the end: both are
  // cell-bound line overlays. The strip order therefore diverges from the Sys
  // tab's on purpose — see buildBordersToggle for why it is composed in here.
  const layers = buildLayerToggles(p);
  const routeAt = layers.findIndex(t => t.key === 'routes');
  const at = routeAt < 0 ? layers.length : routeAt + 1;
  const toggles = [...layers.slice(0, at), buildBordersToggle(p), ...layers.slice(at)];

  return (
    <div className="flex flex-1 flex-nowrap items-end gap-2 min-w-0 overflow-x-auto">
      <Field caption="Projection" className={FIELD}>
        {id => (
          <Select
            value={p.displayMode}
            options={PROJECTION_OPTIONS}
            onChange={p.setDisplayMode}
            label="Map projection"
            labelledBy={id}
            triggerClassName="w-full justify-between"
          />
        )}
      </Field>

      <Field caption="Layer" className={FIELD}>
        {id => (
          <Select
            value={p.viewMode}
            options={VIEW_LAYER_OPTIONS}
            onChange={p.setViewMode}
            label="View layer"
            labelledBy={id}
            triggerClassName="w-full justify-between"
          />
        )}
      </Field>

      {/* Style sits next to the view layer because it is the SIBLING axis, not a
          setting: viewMode picks what the map shows, style picks how it is drawn.
          Shown in every projection: the globe samples a baked texture of the
          same style, so the control applies there too. */}
      <Field caption="Style" className={FIELD}>
        {id => (
          <Select
            value={p.mapStyleId}
            options={MAP_STYLE_OPTIONS}
            onChange={p.setMapStyleId}
            label="Map style"
            labelledBy={id}
            triggerClassName="w-full justify-between"
          />
        )}
      </Field>

      <Field caption="Overlays" className={FIELD}>
        {id => (
          <CheckboxMenu
            items={toggles.map(t => ({
              key: t.key, label: t.label, icon: t.icon,
              checked: t.checked, onChange: t.onChange,
            }))}
            label="Map overlays"
            labelledBy={id}
            triggerClassName="w-full justify-between"
          />
        )}
      </Field>

      {/* Rotation last rather than beside the projection it modifies: it is the
          only control here that is not a choice, and putting a lone icon button
          between two dropdowns broke the row's rhythm. It disappears entirely in
          every projection but the globe, so a trailing slot also means the other
          four never shift position when it goes. */}
      {rotation && (
        <Field caption="Rotation" className="shrink-0">
          {id => (
            <button
              id={rotationId}
              onClick={rotation.onToggle}
              disabled={rotation.disabled}
              aria-pressed={rotation.paused}
              // Caption AND the button's own text, so the accessible name is
              // "Rotation Playing" rather than a caption with no state in it.
              aria-labelledby={`${id} ${rotationId}`}
              title={rotation.paused ? 'Resume rotation' : 'Pause rotation'}
              className={`inline-flex items-center justify-center gap-1.5 px-2.5 py-1.5 border text-[11px] transition-colors ${
                rotation.disabled
                  ? 'bg-surface-raised text-ink-faint border-edge opacity-40 cursor-not-allowed'
                  : rotation.paused
                    ? 'bg-surface-raised text-ink-strong border-edge-strong hover:bg-surface-hover'
                    : 'bg-surface-raised text-ink-muted border-edge hover:bg-surface-hover hover:text-ink-strong'
              }`}
            >
              {rotation.paused ? <Play size={12} /> : <Pause size={12} />}
              {rotation.paused ? 'Paused' : 'Playing'}
            </button>
          )}
        </Field>
      )}
    </div>
  );
};
