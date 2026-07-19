import React, { useMemo } from 'react';
import { Eye, EyeOff, ChevronDown, ChevronUp, Ruler } from 'lucide-react';
import { WorldData, InspectMode, EditMode, TownData } from '../types';

interface WorldDataUpdates {
  townName?: string;
  townPopulation?: number;
  factionName?: string;
  factionColor?: string;
  factionDescription?: string;
  provinceName?: string;
}

interface InspectorProps {
  world: WorldData | null;
  cellId: number | null;
  inspectMode: InspectMode;
  collapsed: boolean;
  onToggleEnabled: () => void;
  onToggleCollapsed: () => void;
  editMode?: EditMode;
  onEditWorldData?: (cellId: number, updates: WorldDataUpdates) => void;
  rulerActive?: boolean;
  onToggleRuler?: () => void;
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
  rulerActive = false,
  onToggleRuler,
}) => {
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
      <div className={`bg-black/80 backdrop-blur text-white shadow-xl border border-white/20 transition-all duration-300 pointer-events-auto ${collapsed ? 'w-28' : isEditing ? 'w-64' : 'min-w-[220px]'}`}>
        <div className="flex items-center justify-between p-2 border-b border-white/10">
          {!collapsed && (
            <span className="text-[10px] font-bold text-gray-400 uppercase tracking-tighter">
              {editMode === 'world-edit' ? 'World Edit' : 'Inspector'}
            </span>
          )}
          <div className="flex items-center gap-2">
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

        {/* Read-only inspect view */}
        {!collapsed && enabled && cell && !isEditing && (
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
            </div>
          </div>
        )}

        {/* World-edit view with editable fields */}
        {!collapsed && isEditing && cell && (
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
            {province && (
              <div className="space-y-1.5">
                <p className="text-[10px] text-gray-500 uppercase tracking-wide">Province</p>
                <input
                  type="text"
                  value={province.name}
                  onChange={e => emit({ provinceName: e.target.value })}
                  className={inputCls}
                  placeholder="Province name"
                />
              </div>
            )}

            {/* Town fields */}
            {town && (
              <div className="space-y-1.5">
                <p className="text-[10px] text-gray-500 uppercase tracking-wide">{town.isCapital ? 'Capital' : 'Town'}</p>
                <input
                  type="text"
                  value={town.name}
                  onChange={e => emit({ townName: e.target.value })}
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

        {!collapsed && !enabled && (
          <div className="p-4 text-[10px] text-gray-500 text-center italic">
            Inspector Disabled
          </div>
        )}
        {!collapsed && enabled && !cell && (
          <div className="p-4 text-[10px] text-gray-500 text-center italic">
            {editMode === 'world-edit' ? 'Click a cell to edit...' : 'Click a cell...'}
          </div>
        )}
      </div>
    </div>
  );
};

export default Inspector;
