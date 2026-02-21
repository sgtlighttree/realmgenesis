import React, { useMemo } from 'react';
import { Eye, EyeOff, ChevronDown, ChevronUp } from 'lucide-react';
import { WorldData, InspectMode } from '../types';

interface InspectorProps {
  world: WorldData | null;
  cellId: number | null;
  inspectMode: InspectMode;
  collapsed: boolean;
  onToggleEnabled: () => void;
  onToggleCollapsed: () => void;
}

const Inspector: React.FC<InspectorProps> = ({
  world,
  cellId,
  inspectMode,
  collapsed,
  onToggleMode,
  onToggleEnabled,
  onToggleCollapsed
}) => {
  // eslint-disable-next-line security/detect-object-injection, nosemgrep
        // codacy-disable-next-line
  const cell = world && cellId !== null ? world.cells[cellId] : null;
  const enabled = inspectMode === 'click';

  const factionMap = useMemo(() => {
    if (!world?.civData) return new Map();
    return new Map(world.civData.factions.map(f => [f.id, f]));
  }, [world?.civData]);

  const locationName = useMemo(() => {
    if (!cell || !world?.civData || !enabled) return null;
    const { regionId, provinceId } = cell;
    if (regionId === undefined) return "Unclaimed Land";
    const faction = factionMap.get(regionId);
    if (!faction) return `Region ${regionId}`;
    // eslint-disable-next-line security/detect-object-injection, nosemgrep
        // codacy-disable-next-line
    const province = provinceId !== undefined && faction.provinces[provinceId];
    return (
      <div className="flex flex-col">
        <span className="font-bold text-blue-200">{faction.name}</span>
        {province && <span className="text-xs text-blue-100/70">{province.name}</span>}
      </div>
    );
  }, [cell, world?.civData, enabled, factionMap]);

  return (
    <div className="absolute top-6 left-1/2 -translate-x-1/2 flex flex-col items-center gap-2 pointer-events-none z-10">
      <div className={`bg-black/80 backdrop-blur text-white shadow-xl border border-white/20 transition-all duration-300 pointer-events-auto ${collapsed ? 'w-28' : 'min-w-[220px]'}`}>
        <div className="flex items-center justify-between p-2 border-b border-white/10">
          {!collapsed && <span className="text-[10px] font-bold text-gray-400 uppercase tracking-tighter">Inspector</span>}
          <div className="flex items-center gap-2">
            <button
              onClick={onToggleEnabled}
              className={`p-1 transition-colors ${enabled ? 'text-blue-400 hover:bg-blue-900/40' : 'text-gray-600 hover:bg-gray-800'}`}
              title={enabled ? "Disable Inspector (Performance Boost)" : "Enable Inspector"}
            >
              {enabled ? <Eye size={12} /> : <EyeOff size={12} />}
            </button>
            <button
              onClick={onToggleCollapsed}
              className="p-1 text-gray-400 hover:bg-gray-800"
            >
              {collapsed ? <ChevronDown size={12} /> : <ChevronUp size={12} />}
            </button>
          </div>
        </div>

        {!collapsed && enabled && cell && (
          <div className="p-2 text-xs">
            <div className="font-bold flex justify-between gap-4 mb-2 border-b border-white/10 pb-1 items-start">
              <div className="flex flex-col">
                <span>Cell {cell.id}</span>
                {locationName && <div className="mt-1">{locationName}</div>}
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
        {!collapsed && !enabled && (
          <div className="p-4 text-[10px] text-gray-500 text-center italic">
            Inspector Disabled
          </div>
        )}
        {!collapsed && enabled && !cell && (
          <div className="p-4 text-[10px] text-gray-500 text-center italic">
            Click a cell...
          </div>
        )}
      </div>
    </div>
  );
};

export default Inspector;
