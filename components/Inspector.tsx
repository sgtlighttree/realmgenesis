import React, { useEffect, useMemo, useState } from 'react';
import { Eye, EyeOff, ChevronDown, ChevronUp, Ruler, MapPin, Trash2, Star } from 'lucide-react';
import { WorldData, InspectMode, EditMode, TownData, MarkerKind } from '../types';

interface WorldDataUpdates {
  townName?: string;
  townPopulation?: number;
  factionName?: string;
  factionColor?: string;
  factionDescription?: string;
  provinceName?: string;
}

interface MarkerUpdates {
  kind?: MarkerKind;
  name?: string;
  note?: string;
}

const MARKER_KINDS: MarkerKind[] = ['dungeon', 'ruin', 'battlefield', 'portal', 'poi'];

interface InspectorProps {
  world: WorldData | null;
  cellId: number | null;
  inspectMode: InspectMode;
  collapsed: boolean;
  onToggleEnabled: () => void;
  onToggleCollapsed: () => void;
  editMode?: EditMode;
  onEditWorldData?: (cellId: number, updates: WorldDataUpdates) => void;
  onRenameProvince?: (factionId: number, provinceId: number, name: string) => void;
  onRenameTown?: (factionId: number, provinceId: number, cellId: number, name: string) => void;
  onRelocateCapital?: (factionId: number, townCellId: number) => void;
  rulerActive?: boolean;
  onToggleRuler?: () => void;
  markerMode?: boolean;
  onToggleMarkerMode?: () => void;
  selectedMarkerId?: number | null;
  onSelectMarker?: (id: number) => void;
  onUpdateMarker?: (id: number, updates: MarkerUpdates) => void;
  onDeleteMarker?: (id: number) => void;
}

const Inspector: React.FC<InspectorProps> = ({
  world,
  cellId,
  inspectMode,
  collapsed,
  onToggleEnabled,
  onToggleCollapsed,
  editMode = 'off',
  onEditWorldData,
  onRenameProvince,
  onRenameTown,
  onRelocateCapital,
  rulerActive = false,
  onToggleRuler,
  markerMode = false,
  onToggleMarkerMode,
  selectedMarkerId = null,
  onSelectMarker,
  onUpdateMarker,
  onDeleteMarker,
}) => {
  const [markersListOpen, setMarkersListOpen] = useState(false);
  // Local drafts for the province/town name inputs below — committed (and
  // validated) on blur rather than on every keystroke, so a rejected empty
  // value reverts to the real name instead of the field looking frozen.
  const [provinceDraft, setProvinceDraft] = useState<string | null>(null);
  const [townDraft, setTownDraft] = useState<string | null>(null);
  useEffect(() => {
    setProvinceDraft(null);
    setTownDraft(null);
  }, [cellId]);
  const markers = world?.markers ?? [];
  const selectedMarker = selectedMarkerId !== null ? markers.find(m => m.id === selectedMarkerId) ?? null : null;
  const cell = world && cellId !== null ? world.cells[cellId] : null;
  const enabled = inspectMode === 'click' || editMode === 'world-edit';
  const isEditing = editMode === 'world-edit' && cell !== null && onEditWorldData !== undefined;

  const factionMap = useMemo(() => {
    if (!world?.civData) return new Map();
    return new Map(world.civData.factions.map(f => [f.id, f]));
  }, [world?.civData]);

  // Lazy cell -> geographic feature name lookup. First match wins, and detection
  // orders specific kinds (ranges/deserts/forests) ahead of broad ones (islands),
  // so a mountain on an island reads as its range.
  const featureByCell = useMemo(() => {
    const map = new Map<number, string>();
    for (const feature of world?.features ?? []) {
      for (const id of feature.cellIds) {
        if (!map.has(id)) map.set(id, feature.name);
      }
    }
    return map;
  }, [world?.features]);

  const featureName = cell ? featureByCell.get(cell.id) ?? null : null;

  const cultureMap = useMemo(() => {
    const map = new Map<number, string>();
    for (const culture of world?.cultures ?? []) map.set(culture.id, culture.name);
    return map;
  }, [world?.cultures]);
  const cultureName = cell?.cultureId !== undefined ? cultureMap.get(cell.cultureId) ?? null : null;

  const religionMap = useMemo(() => {
    const map = new Map<number, string>();
    for (const religion of world?.religions ?? []) map.set(religion.id, religion.name);
    return map;
  }, [world?.religions]);
  const religionName = cell?.religionId !== undefined ? religionMap.get(cell.religionId) ?? null : null;

  const faction = cell?.regionId !== undefined ? factionMap.get(cell.regionId) : null;
  const province = (faction && cell?.provinceId !== undefined) ? faction.provinces[cell.provinceId] : null;
  const town = province?.towns.find((t: TownData) => t.cellId === cell?.id) ?? null;

  const locationName = useMemo(() => {
    if (!cell || !world?.civData || !enabled) return null;
    const { regionId, provinceId } = cell;
    if (regionId === undefined) return "Unclaimed Land";
    const f = factionMap.get(regionId);
    if (!f) return `Region ${regionId}`;
    const p = provinceId !== undefined && f.provinces[provinceId];
    return (
      <div className="flex flex-col">
        <span className="font-bold text-blue-200">{f.name}</span>
        {p && <span className="text-xs text-blue-100/70">{p.name}</span>}
      </div>
    );
  }, [cell, world?.civData, enabled, factionMap]);

  const emit = (updates: WorldDataUpdates) => {
    if (cell && onEditWorldData) onEditWorldData(cell.id, updates);
  };

  const inputCls = "bg-gray-900 border border-gray-700 text-white text-xs px-1.5 py-1 w-full focus:outline-none focus:border-blue-500";

  return (
    <div className="absolute top-6 left-1/2 -translate-x-1/2 flex flex-col items-center gap-2 pointer-events-none z-10">
      <div className={`bg-black/80 backdrop-blur text-white shadow-xl border border-white/20 transition-all duration-300 pointer-events-auto ${collapsed ? 'w-28' : selectedMarker || isEditing ? 'w-64' : 'min-w-[220px]'}`}>
        <div className="flex items-center justify-between p-2 border-b border-white/10">
          {!collapsed && (
            <span className="text-[10px] font-bold text-gray-400 uppercase tracking-tighter">
              {editMode === 'world-edit' ? 'World Edit' : 'Inspector'}
            </span>
          )}
          <div className="flex items-center gap-2">
            {onToggleMarkerMode && (
              <button
                onClick={onToggleMarkerMode}
                className={`p-1 transition-colors ${markerMode ? 'text-amber-400 hover:bg-amber-900/40' : 'text-gray-600 hover:bg-gray-800'}`}
                title={markerMode ? 'Disable Marker Placement' : 'Place Marker'}
              >
                <MapPin size={12} />
              </button>
            )}
            {onToggleRuler && (
              <button
                onClick={onToggleRuler}
                className={`p-1 transition-colors ${rulerActive ? 'text-amber-400 hover:bg-amber-900/40' : 'text-gray-600 hover:bg-gray-800'}`}
                title={rulerActive ? 'Disable Ruler' : 'Measure Distance'}
              >
                <Ruler size={12} />
              </button>
            )}
            {editMode !== 'world-edit' && (
              <button
                onClick={onToggleEnabled}
                className={`p-1 transition-colors ${enabled ? 'text-blue-400 hover:bg-blue-900/40' : 'text-gray-600 hover:bg-gray-800'}`}
                title={enabled ? "Disable Inspector" : "Enable Inspector"}
              >
                {enabled ? <Eye size={12} /> : <EyeOff size={12} />}
              </button>
            )}
            <button onClick={onToggleCollapsed} className="p-1 text-gray-400 hover:bg-gray-800">
              {collapsed ? <ChevronDown size={12} /> : <ChevronUp size={12} />}
            </button>
          </div>
        </div>

        {/* Marker editing panel — takes priority over cell views while a marker is selected */}
        {!collapsed && selectedMarker && (
          <div className="p-2 text-xs space-y-1.5">
            <div className="flex items-center justify-between border-b border-white/10 pb-1 mb-1">
              <span className="font-bold text-amber-300">Marker</span>
              {onDeleteMarker && (
                <button
                  onClick={() => onDeleteMarker(selectedMarker.id)}
                  className="p-0.5 text-red-400 hover:text-red-300"
                  title="Delete marker"
                >
                  <Trash2 size={12} />
                </button>
              )}
            </div>
            <select
              value={selectedMarker.kind}
              onChange={e => onUpdateMarker?.(selectedMarker.id, { kind: e.target.value as MarkerKind })}
              className={inputCls}
            >
              {MARKER_KINDS.map(k => <option key={k} value={k}>{k}</option>)}
            </select>
            <input
              type="text"
              value={selectedMarker.name}
              onChange={e => onUpdateMarker?.(selectedMarker.id, { name: e.target.value })}
              className={inputCls}
              placeholder="Marker name"
            />
            <textarea
              value={selectedMarker.note}
              onChange={e => onUpdateMarker?.(selectedMarker.id, { note: e.target.value })}
              rows={3}
              className={inputCls + ' resize-none'}
              placeholder="Notes..."
            />
          </div>
        )}

        {/* Read-only inspect view */}
        {!collapsed && !selectedMarker && enabled && cell && !isEditing && (
          <div className="p-2 text-xs">
            <div className="font-bold flex justify-between gap-4 mb-2 border-b border-white/10 pb-1 items-start">
              <div className="flex flex-col">
                <span>Cell {cell.id}</span>
                {locationName && <div className="mt-1">{locationName}</div>}
                {featureName && <div className="mt-1 text-[10px] text-gray-400 italic">Part of: <span className="text-gray-200">{featureName}</span></div>}
              </div>
              <span style={{ color: '#aaa' }}>{cell.biome}</span>
            </div>
            <div className="grid grid-cols-2 gap-x-4 gap-y-1">
              <div className="text-gray-400">Temp: <span className="text-white">{cell.temperature.toFixed(1)}°C</span></div>
              <div className="text-gray-400">Rain: <span className="text-white">{(cell.moisture * 100).toFixed(0)}%</span></div>
              <div className="text-gray-400">Elev: <span className="text-white">{(cell.height * 100).toFixed(0)}%</span></div>
              <div className="text-gray-400">Pop: <span className="text-white">{cell.population?.toLocaleString()}</span></div>
              {cultureName && (
                <div className="text-gray-400 col-span-2">Culture: <span className="text-white">{cultureName}</span></div>
              )}
              {religionName && (
                <div className="text-gray-400 col-span-2">Faith: <span className="text-white">{religionName}</span></div>
              )}
            </div>
          </div>
        )}

        {/* World-edit view with editable fields */}
        {!collapsed && !selectedMarker && isEditing && cell && (
          <div className="p-2 text-xs space-y-2">
            {/* Read-only physical stats */}
            <div className="grid grid-cols-2 gap-x-4 gap-y-1 pb-2 border-b border-white/10">
              <div className="text-gray-400">Temp: <span className="text-white">{cell.temperature.toFixed(1)}°C</span></div>
              <div className="text-gray-400">Rain: <span className="text-white">{(cell.moisture * 100).toFixed(0)}%</span></div>
              <div className="text-gray-400">Elev: <span className="text-white">{(cell.height * 100).toFixed(0)}%</span></div>
              <div className="text-gray-400">Biome: <span className="text-white">{cell.biome}</span></div>
            </div>

            {/* Faction fields */}
            {faction && (
              <div className="space-y-1.5">
                <p className="text-[10px] text-gray-500 uppercase tracking-wide">Faction</p>
                <div className="flex gap-1.5 items-center">
                  <input
                    type="color"
                    value={faction.color}
                    onChange={e => emit({ factionColor: e.target.value })}
                    className="w-7 h-6 border border-gray-700 bg-transparent cursor-pointer"
                    title="Faction color"
                  />
                  <input
                    type="text"
                    value={faction.name}
                    onChange={e => emit({ factionName: e.target.value })}
                    className={inputCls + ' flex-1'}
                    placeholder="Faction name"
                  />
                </div>
                {faction.description !== undefined && (
                  <textarea
                    value={faction.description}
                    onChange={e => emit({ factionDescription: e.target.value })}
                    rows={2}
                    className={inputCls + ' resize-none'}
                    placeholder="Faction description"
                  />
                )}
              </div>
            )}

            {/* Province fields */}
            {province && faction && (
              <div className="space-y-1.5">
                <p className="text-[10px] text-gray-500 uppercase tracking-wide">Province</p>
                <input
                  type="text"
                  value={provinceDraft ?? province.name}
                  onChange={e => setProvinceDraft(e.target.value)}
                  onBlur={() => {
                    if (provinceDraft !== null) onRenameProvince?.(faction.id, province.id, provinceDraft);
                    setProvinceDraft(null);
                  }}
                  className={inputCls}
                  placeholder="Province name"
                />
              </div>
            )}

            {/* Town fields */}
            {town && faction && province && (
              <div className="space-y-1.5">
                <div className="flex items-center justify-between">
                  <p className="text-[10px] text-gray-500 uppercase tracking-wide">{town.isCapital ? 'Capital' : 'Town'}</p>
                  {!town.isCapital && onRelocateCapital && (
                    <button
                      onClick={() => onRelocateCapital(faction.id, town.cellId)}
                      className="flex items-center gap-1 px-1.5 py-0.5 text-[9px] font-semibold uppercase tracking-wide bg-amber-900/40 text-amber-300 border border-amber-700 hover:bg-amber-900/70"
                      title="Promote this town to faction capital"
                    >
                      <Star size={9} /> Make Capital
                    </button>
                  )}
                </div>
                <input
                  type="text"
                  value={townDraft ?? town.name}
                  onChange={e => setTownDraft(e.target.value)}
                  onBlur={() => {
                    if (townDraft !== null) onRenameTown?.(faction.id, province.id, town.cellId, townDraft);
                    setTownDraft(null);
                  }}
                  className={inputCls}
                  placeholder="Town name"
                />
                <input
                  type="number"
                  value={town.population}
                  onChange={e => emit({ townPopulation: parseInt(e.target.value) || 0 })}
                  className={inputCls}
                  placeholder="Population"
                  min={0}
                />
              </div>
            )}

            {!faction && (
              <p className="text-gray-500 italic text-center py-1">Unclaimed — no faction data</p>
            )}
          </div>
        )}

        {!collapsed && !selectedMarker && !enabled && (
          <div className="p-4 text-[10px] text-gray-500 text-center italic">
            Inspector Disabled
          </div>
        )}
        {!collapsed && !selectedMarker && enabled && !cell && (
          <div className="p-4 text-[10px] text-gray-500 text-center italic">
            {markerMode ? 'Click the map to place a marker...' : editMode === 'world-edit' ? 'Click a cell to edit...' : 'Click a cell...'}
          </div>
        )}

        {/* Markers list — always available at the bottom when any markers exist */}
        {!collapsed && markers.length > 0 && (
          <div className="border-t border-white/10">
            <button
              onClick={() => { setMarkersListOpen(v => !v); }}
              className="w-full flex items-center justify-between p-2 text-[10px] text-gray-400 hover:bg-gray-800"
            >
              <span>Markers ({markers.length})</span>
              {markersListOpen ? <ChevronUp size={10} /> : <ChevronDown size={10} />}
            </button>
            {markersListOpen && (
              <div className="max-h-32 overflow-y-auto">
                {markers.map(m => (
                  <button
                    key={m.id}
                    onClick={() => onSelectMarker?.(m.id)}
                    className={`w-full flex items-center justify-between gap-2 px-2 py-1 text-[10px] text-left ${m.id === selectedMarkerId ? 'bg-blue-900/40 text-blue-200' : 'text-gray-300 hover:bg-gray-800'}`}
                  >
                    <span className="truncate">{m.name}</span>
                    <span className="text-gray-500 uppercase shrink-0">{m.kind}</span>
                  </button>
                ))}
              </div>
            )}
          </div>
        )}
      </div>
    </div>
  );
};

export default Inspector;
