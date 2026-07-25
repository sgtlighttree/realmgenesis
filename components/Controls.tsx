import React, { useState, useEffect, useRef } from 'react';
import { WorldParams, ViewMode, LoreData, LandStyle, CivData, DisplayMode, DymaxionSettings, DymaxionControlMode, LabelVisibility, MarkerData } from '../types';
import { RefreshCw, Globe, Mountain, Lock, Unlock, Shuffle, Layers, Zap, Save, Trash2, Image, Terminal, XCircle, ChevronDown, ChevronUp, FolderOpen, Box, Copy, Check, FileCode, FileJson } from 'lucide-react';
import { exportMap, saveMapConfig, loadMapConfig, saveMapToBrowser, getSavedMaps, deleteSavedMap, ExportResolution, ProjectionType } from '../utils/export';
import { downloadSVG, downloadGeoJSON } from '../utils/exportVector';
import { NAME_STYLES, NameStyle } from '../utils/namegen';
import { exportGLB } from '../utils/exportGLB';
import { WorldData } from '../types';
import DymaxionPreview2D from './DymaxionPreview2D';
import {
  ViewControlsProps, DISPLAY_MODES, DisplayButton, ViewLayerGrid,
  LayerToggleRow, OverlayToggles, buildLayerToggles,
} from './ViewControls';

interface ControlsProps {
  params: WorldParams;
  setParams: React.Dispatch<React.SetStateAction<WorldParams>>;
  onGenerate: (p?: WorldParams) => void;
  onLoadWorld: (p: WorldParams, civData?: CivData, markers?: MarkerData[]) => void;
  onCancel?: () => void;
  onUpdateCivs: (p?: WorldParams) => void;
  onUpdateProvinces: (p?: WorldParams) => void;
  viewMode: ViewMode;
  setViewMode: (m: ViewMode) => void;
  displayMode: DisplayMode;
  setDisplayMode: (m: DisplayMode) => void;
  lore: LoreData | null;
  loading: boolean;
  generatingLore: boolean;
  onGenerateLore: () => void;
  worldData: WorldData | null;
  logs: string[];
  genProgress: number;
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
  dymaxionSettings: DymaxionSettings;
  onDymaxionChange: React.Dispatch<React.SetStateAction<DymaxionSettings>>;
  apiKey: string;
  onApiKeyChange: (key: string) => void;
  onInspect: (id: number | null) => void;
  onEditFaction?: (factionId: number, updates: { name?: string; color?: string; description?: string }) => void;
  onMergeFactions?: (srcId: number, dstId: number) => void;
}

type Tab = 'geo' | 'climate' | 'political' | 'system' | 'export';

const ConsoleOutput: React.FC<{ logs: string[]; isOpen: boolean }> = ({ logs, isOpen }) => {
    const endRef = useRef<HTMLDivElement>(null);
    useEffect(() => {
        if (isOpen) {
            endRef.current?.scrollIntoView({ behavior: "smooth" });
        }
    }, [logs, isOpen]);

    if (!isOpen) return null;

    return (
        <div className="bg-black border border-gray-800 p-2 h-32 overflow-y-auto font-mono text-[10px] space-y-1 shadow-inner relative transition-all">
            {logs.length === 0 && <div className="text-gray-600 italic text-center mt-10">System Ready</div>}
            {logs.map((log, i) => (
                <div key={i} className="text-green-400 break-words border-b border-gray-900/50 pb-0.5 last:border-0">
                    <span className="text-gray-600 mr-2">[{i+1}]</span>
                    {log}
                </div>
            ))}
            <div ref={endRef} />
        </div>
    );
};

const Controls: React.FC<ControlsProps> = ({
  params,
  setParams,
  onGenerate,
  onLoadWorld,
  onCancel,
  onUpdateCivs,
  onUpdateProvinces,
  viewMode,
  setViewMode,
  displayMode,
  setDisplayMode,
  lore,
  loading,
  generatingLore,
  onGenerateLore,
  worldData,
  logs,
  genProgress,
  showGrid,
  setShowGrid,
  showRivers,
  setShowRivers,
  showRoutes,
  setShowRoutes,
  showHillshade,
  setShowHillshade,
  showContours,
  setShowContours,
  labelVisibility,
  setLabelVisibility,
  dymaxionSettings,
  onDymaxionChange,
  apiKey,
  onApiKeyChange,
  onInspect,
  onEditFaction,
  onMergeFactions,
}) => {
  const [activeTab, setActiveTab] = useState<Tab>('system');
  const [seedLocked, setSeedLocked] = useState(false);
  const [civSeedLocked, setCivSeedLocked] = useState(false);
  const [autoUpdate, setAutoUpdate] = useState(true);
  const [consoleOpen, setConsoleOpen] = useState(true);
  const [showCivParams, setShowCivParams] = useState(false);
  // Per-faction "merge into" target, keyed by source faction id. Two-step:
  // pick a target in the select, then click Merge to confirm (no window.confirm).
  const [mergeTargets, setMergeTargets] = useState<Record<number, number | ''>>({});
  
  // Export State
  const [expRes, setExpRes] = useState<ExportResolution>(4096);
  const [expProj, setExpProj] = useState<ProjectionType>('equirectangular');
  const [saveName, setSaveName] = useState('');
  const [savedMaps, setSavedMaps] = useState(getSavedMaps());
  const [showDymaxion2D, setShowDymaxion2D] = useState(false);
  const [seedCopied, setSeedCopied] = useState(false);
  const [showStats, setShowStats] = useState(false);

  const updateDymaxion = (patch: Partial<DymaxionSettings>) => {
      onDymaxionChange((prev) => ({ ...prev, ...patch }));
  };

  // Global keyboard shortcuts
  useEffect(() => {
    const TABS: Tab[] = ['system', 'geo', 'climate', 'political', 'export'];
    const onKey = (e: KeyboardEvent) => {
      const tag = (e.target as HTMLElement).tagName;
      if (tag === 'INPUT' || tag === 'TEXTAREA' || tag === 'SELECT') return;
      switch (e.key) {
        case 'g': case 'G': if (!loading) { handleGenerateClick(); } break;
        case 'r': case 'R': handleRandomizeSeed(); break;
        case 'Escape': onInspect(null); break;
        case '1': case '2': case '3': case '4': case '5':
          setActiveTab(TABS[parseInt(e.key) - 1]);
          break;
      }
    };
    window.addEventListener('keydown', onKey);
    return () => window.removeEventListener('keydown', onKey);
    // Re-register whenever anything the handlers close over changes, so the
    // G shortcut never generates from (or writes back) stale params
  }, [loading, params, seedLocked, civSeedLocked, onInspect]); // eslint-disable-line react-hooks/exhaustive-deps

  useEffect(() => {
     if (autoUpdate && !loading && params.points <= 20000) {
         const timer = setTimeout(() => {
             onGenerate();
         }, 400); 
         return () => { clearTimeout(timer); };
     }
  }, [
      // Dependency list for auto-update
      params.landStyle, 
      params.noiseScale, 
      params.roughness, 
      params.seaLevel, 
      params.plates, 
      params.warpStrength, 
      params.plateInfluence,
      params.ridgeBlend,
      params.mountainHeight,
      params.oceanDepth,
      params.cellJitter,
      params.detailLevel,
      params.erosionIterations,
      params.baseTemperature,
      params.poleTemperature,
      params.rainfallMultiplier,
      params.moistureTransport,
      params.temperatureVariance,
      params.axialTilt,
      autoUpdate
  ]);

  // Update default save name when active tab changes to export or system
  useEffect(() => {
      if (activeTab === 'export' || activeTab === 'system') {
          const now = new Date();
          const yymmdd = now.toISOString().slice(2,10).replace(/-/g, '');
          const hhmmss = now.toTimeString().slice(0,8).replace(/:/g, '');
          const defaultName = `${params.mapName || 'map'}_${yymmdd}_${hhmmss}`;
          if (!saveName || saveName.includes(params.mapName || 'map')) {
              setSaveName(defaultName);
          }
      }
  }, [activeTab, params.mapName]);

  const clamp = (val: number, min: number, max: number) => Math.max(min, Math.min(max, val));

  const handleChange = <K extends keyof WorldParams>(key: K, value: WorldParams[K]) => {  
    setParams({ ...params, [key]: value });
  };

  const handleNumberChange = <K extends keyof WorldParams>(key: K, rawValue: string, min: number, max: number, step?: number) => {
    let val = step ? parseFloat(rawValue) : parseInt(rawValue, 10);
    if (isNaN(val)) return;
    val = clamp(val, min, max);
    setParams({ ...params, [key]: val as WorldParams[K] });
  };

  const handleAdvancedChange = <K extends keyof WorldParams>(key: K, value: WorldParams[K]) => {  
      setParams({ ...params, [key]: value, landStyle: 'Custom' });
  };

  const handlePresetChange = (style: LandStyle) => {
      let updates: Partial<WorldParams> = { landStyle: style };
      switch(style) {
          case 'Continents':
              updates.noiseScale = 1.0;
              updates.ridgeBlend = 0.5;
              updates.maskType = 'None';
              updates.warpStrength = 0.2;
              updates.plateInfluence = 0.8;
              updates.erosionIterations = 2;
              updates.cellJitter = 0.5;
              updates.mountainHeight = 1.0;
              updates.oceanDepth = 1.0;
              break;
          case 'Pangea':
              updates.noiseScale = 0.8;
              updates.ridgeBlend = 0.4;
              updates.maskType = 'Pangea';
              updates.warpStrength = 0.5;
              updates.plateInfluence = 0.6;
              updates.erosionIterations = 2;
              updates.cellJitter = 0.4;
              updates.mountainHeight = 1.0;
              updates.oceanDepth = 1.0;
              break;
          case 'Archipelago':
              updates.noiseScale = 3.0;
              updates.ridgeBlend = 0.55;
              updates.maskType = 'None';
              updates.warpStrength = 0.7;
              updates.plateInfluence = 0.25;
              updates.erosionIterations = 4;
              updates.cellJitter = 0.8;
              updates.seaLevel = 0.65;
              updates.plates = 6;
              updates.roughness = 0.55;
              updates.mountainHeight = 0.8;
              updates.oceanDepth = 1.3;
              break;
          case 'Islands':
              updates.noiseScale = 1.8;
              updates.ridgeBlend = 0.25;
              updates.maskType = 'None';
              updates.warpStrength = 1.0;
              updates.plateInfluence = 0.4;
              updates.erosionIterations = 4;
              updates.cellJitter = 0.65;
              updates.seaLevel = 0.60;
              updates.plates = 8;
              updates.roughness = 0.55;
              updates.mountainHeight = 1.0;
              updates.oceanDepth = 1.1;
              break;
          case 'Custom':
              break;
      }
      setParams({ ...params, ...updates });
  };

  const handleRandomizeSeed = () => {
    if (seedLocked) return;
    setParams(prev => ({
      ...prev,
      seed: crypto.getRandomValues(new Uint32Array(1))[0].toString(36),
      ...(civSeedLocked ? {} : { civSeed: crypto.getRandomValues(new Uint32Array(1))[0].toString(36) }),
    }));
  };
  
  const handleRandomizeCivSeed = () => {
    if (!civSeedLocked) {
        handleChange('civSeed', crypto.getRandomValues(new Uint32Array(1))[0].toString(36));
    }
  };

  const handleGenerateClick = () => {
    let p = { ...params };
    if (!seedLocked) {
       p.seed = crypto.getRandomValues(new Uint32Array(1))[0].toString(36);
    }
    if (!civSeedLocked) {
       p.civSeed = crypto.getRandomValues(new Uint32Array(1))[0].toString(36);
    }
    setParams(p);
    setTimeout(() => {
        onGenerate(p);
    }, 0);
  };

  const handleRerollBorders = () => {
      let newCivSeed = params.civSeed;
      if (!civSeedLocked) {
          newCivSeed = crypto.getRandomValues(new Uint32Array(1))[0].toString(36);
          setParams({ ...params, civSeed: newCivSeed });
      }
      // Pass the updated params explicitly so the callback uses the new seed immediately
      setTimeout(() => onUpdateCivs({ ...params, civSeed: newCivSeed }), 50);
  };

  const handleRerollProvinces = () => {
      // Just triggers the province recalculation logic
      setTimeout(() => onUpdateProvinces(), 50);
  };
  
  const handleExport = async () => {
    if (!worldData) return;
    try {
        await exportMap(
          worldData,
          viewMode,
          expRes,
          expProj,
          expProj === 'dymaxion' ? { layout: dymaxionSettings.layout, lon: dymaxionSettings.lon, lat: dymaxionSettings.lat, roll: dymaxionSettings.roll } : undefined,
          labelVisibility,
          showHillshade,
          showContours,
          showRoutes
        );
    } catch(e) {
        console.error(e);
        alert("Export failed. Try a lower resolution.");
    }
  };

  const handleSaveBrowser = () => {
      if (!saveName) return;
      // Pass civData/markers if available to save lore and POIs
      if (saveMapToBrowser(saveName, params, worldData?.civData, worldData?.markers)) {
          setSavedMaps(getSavedMaps());
          // Generate next default name
          const now = new Date();
          const yymmdd = now.toISOString().slice(2,10).replace(/-/g, '');
          const hhmmss = now.toTimeString().slice(0,8).replace(/:/g, '');
          setSaveName(`${params.mapName || 'map'}_${yymmdd}_${hhmmss}`);
      } else {
          alert("Failed to save (storage full?)");
      }
  };

  const handleLoadBrowser = (entryParams: WorldParams, civData?: CivData, markers?: MarkerData[]) => {
      if (confirm("Load this map configuration? Unsaved changes will be lost.")) {
          setParams(entryParams);
          // Use the dedicated load function that handles restoration
          setTimeout(() => onLoadWorld(entryParams, civData, markers), 50);
      }
  };
  
  const handleDeleteBrowser = (name: string) => {
      if (confirm(`Delete ${name}?`)) {
          deleteSavedMap(name);
          setSavedMaps(getSavedMaps());
      }
  };

  const handleFileUpload = async (e: React.ChangeEvent<HTMLInputElement>) => {
      if (e.target.files && e.target.files[0]) {
          const loaded = await loadMapConfig(e.target.files[0]);
          if (loaded) {
              setParams(loaded.params);
              // Use the dedicated load function that handles restoration
              setTimeout(() => onLoadWorld(loaded.params, loaded.civData, loaded.markers), 50);
          } else {
              alert("Invalid config file");
          }
      }
  };

  // Render-mode / view-layer / overlay-toggle primitives are shared with the
  // shell's View strip (components/ViewControls.tsx) — one definition, each
  // host composing its own layout.
  const viewProps: ViewControlsProps = {
    viewMode, setViewMode, displayMode, setDisplayMode,
    showGrid, setShowGrid, showRivers, setShowRivers, showRoutes, setShowRoutes,
    showHillshade, setShowHillshade, showContours, setShowContours,
    labelVisibility, setLabelVisibility,
  };
  const layerToggles = buildLayerToggles(viewProps);

  return (
    <div className="w-full md:w-80 bg-gray-950 border-r border-gray-800 flex flex-col h-full overflow-hidden text-sm relative z-20">
      <div className="p-4 border-b border-gray-800">
        <h1 className="text-xl font-bold text-white flex items-center gap-2">
          <Globe className="text-blue-500" />
          RealmGenesis 3D
        </h1>
      </div>

      <div className="flex border-b border-gray-800">
         <button onClick={() => { setActiveTab('system'); }} className={`flex-1 py-3 text-[10px] font-semibold uppercase tracking-wide ${activeTab === 'system' ? 'text-blue-400 border-b-2 border-blue-500' : 'text-gray-500 hover:text-gray-300'}`}>Sys</button>
         <button onClick={() => { setActiveTab('geo'); }} className={`flex-1 py-3 text-[10px] font-semibold uppercase tracking-wide ${activeTab === 'geo' ? 'text-blue-400 border-b-2 border-blue-500' : 'text-gray-500 hover:text-gray-300'}`}>Geo</button>
         <button onClick={() => { setActiveTab('climate'); }} className={`flex-1 py-3 text-[10px] font-semibold uppercase tracking-wide ${activeTab === 'climate' ? 'text-blue-400 border-b-2 border-blue-500' : 'text-gray-500 hover:text-gray-300'}`}>Clim</button>
         <button onClick={() => { setActiveTab('political'); }} className={`flex-1 py-3 text-[10px] font-semibold uppercase tracking-wide ${activeTab === 'political' ? 'text-blue-400 border-b-2 border-blue-500' : 'text-gray-500 hover:text-gray-300'}`}>Civ</button>
         <button onClick={() => { setActiveTab('export'); }} className={`flex-1 py-3 text-[10px] font-semibold uppercase tracking-wide ${activeTab === 'export' ? 'text-blue-400 border-b-2 border-blue-500' : 'text-gray-500 hover:text-gray-300'}`}>Exp</button>
      </div>

      <div className="flex-1 overflow-y-auto p-4 space-y-6">
        
        {activeTab === 'system' && (
          <div className="space-y-4">
             <div className="space-y-1">
               <label className="text-xs text-gray-400 block">Render Mode</label>
               <div className="flex gap-2">
                 {DISPLAY_MODES.map(m => (
                   <DisplayButton key={m.mode} mode={m.mode} label={m.label}
                     displayMode={displayMode} setDisplayMode={setDisplayMode} />
                 ))}
               </div>
             </div>

             {/* Map Name Input */}
             <div className="space-y-1">
                 <label className="text-xs text-gray-400 block">Map Name</label>
                 <input 
                    type="text" 
                    value={params.mapName} 
                    onChange={(e) => { handleChange('mapName', e.target.value); }}
                    className="w-full bg-gray-900 border border-gray-700 px-2 py-1.5 text-white text-xs"
                    placeholder="My World"
                 />
             </div>

             {/* Seed Input */}
             <div className="bg-gray-900 p-3 border border-gray-800">
                <label className="text-xs text-gray-400 mb-1 block">Seed</label>
                <div className="flex gap-2">
                   <input 
                      type="text" 
                      value={params.seed} 
                      onChange={(e) => { handleChange('seed', e.target.value); }}
                      disabled={seedLocked}
                      className="bg-black border border-gray-700 px-2 py-1 text-white text-xs flex-1 disabled:opacity-50"
                   />
                   <button
                      onClick={() => {
                        void navigator.clipboard.writeText(params.seed);
                        setSeedCopied(true);
                        setTimeout(() => setSeedCopied(false), 1500);
                      }}
                      title="Copy seed to clipboard"
                      className="text-gray-400 hover:text-white transition-colors"
                   >
                      {seedCopied ? <Check size={14} className="text-green-400" /> : <Copy size={14} />}
                   </button>
                   <button
                      onClick={() => { setSeedLocked(!seedLocked); }}
                      className={`${seedLocked ? 'text-blue-500' : 'text-gray-400'} hover:text-white transition-colors`}
                   >
                      {seedLocked ? <Lock size={14}/> : <Unlock size={14}/>}
                   </button>
                   <button onClick={handleRandomizeSeed} disabled={seedLocked} className="text-gray-400 hover:text-white disabled:opacity-50">
                      <Shuffle size={14} />
                   </button>
                </div>
             </div>
             
             <div className="space-y-1">
              <div className="flex justify-between items-center text-xs text-gray-400">
                <label>Resolution</label>
                <input
                    type="number"
                    min="2000"
                    max="200000"
                    step="1000"
                    value={params.points}
                    onChange={(e) => { handleNumberChange('points', e.target.value, 2000, 200000); }}
                    className="w-24 bg-gray-900 border border-gray-700 px-1 py-0.5 text-right text-white text-xs"
                />
              </div>
              <input
                type="range"
                min="2000"
                max="200000"
                step="1000"
                value={Math.min(200000, params.points)}
                onChange={(e) => { handleChange('points', parseInt(e.target.value) as 1 | 2 | 3); }}
                className="w-full h-1 bg-gray-700 appearance-none cursor-pointer accent-blue-500"
              />
            </div>
             
             {layerToggles.map((t, i) => (
               <LayerToggleRow key={t.key} toggle={t}
                 className={i === 0 ? 'pt-2 border-t border-gray-800' : 'pt-2'} />
             ))}

             <OverlayToggles
               labelVisibility={labelVisibility}
               setLabelVisibility={setLabelVisibility}
             />

            <div className="flex items-center justify-between text-xs text-gray-400 pt-2">
                 <div className="flex items-center gap-2">
                    <Zap size={12} className={autoUpdate ? "text-yellow-400" : "text-gray-600"}/>
                    <label>Auto-Update (Low Res)</label>
                 </div>
                 <input 
                    type="checkbox"
                    checked={autoUpdate}
                    onChange={(e) => { setAutoUpdate(e.target.checked); }}
                    disabled={params.points > 20000}
                    className="bg-gray-700"
                 />
            </div>
            
            <div className="pt-4 border-t border-gray-800">
              <h3 className="text-xs font-semibold text-gray-500 mb-2">View Layer</h3>
              <ViewLayerGrid viewMode={viewMode} setViewMode={setViewMode} />
            </div>

            <div className="pt-4 border-t border-gray-800 space-y-3">
              <h3 className="text-xs font-semibold text-gray-500 mb-2">AI Settings (BYOK)</h3>
              <div className="bg-gray-900 p-3 border border-gray-800 space-y-2">
                <div className="flex items-center justify-between">
                  <label className="text-xs text-gray-400">Gemini API Key</label>
                  <a 
                    href="https://aistudio.google.com/app/apikey" 
                    target="_blank" 
                    rel="noopener noreferrer"
                    className="text-[10px] text-blue-400 hover:underline flex items-center gap-1"
                  >
                    Get Key <Layers size={8} />
                  </a>
                </div>
                <input 
                  type="password"
                  value={apiKey}
                  onChange={(e) => { onApiKeyChange(e.target.value); }}
                  placeholder="Paste your API key here..."
                  className="w-full bg-black border border-gray-700 px-2 py-1.5 text-white text-xs"
                />
                <p className="text-[9px] text-gray-500 italic">
                  Key is stored ephemerally in memory and will be lost on refresh.
                </p>
              </div>
            </div>

            {/* World Stats */}
            {worldData && (
              <div className="pt-4 border-t border-gray-800">
                <button
                  className="flex items-center justify-between w-full text-xs text-gray-400 hover:text-gray-200"
                  onClick={() => setShowStats(v => !v)}
                >
                  <span className="font-semibold text-gray-500">World Stats</span>
                  {showStats ? <ChevronUp size={12}/> : <ChevronDown size={12}/>}
                </button>
                {showStats && (() => {
                  const cells = worldData.cells;
                  const seaLevel = worldData.params.seaLevel;
                  const landCells = cells.filter(c => c.height >= seaLevel);
                  const landPct = ((landCells.length / cells.length) * 100).toFixed(1);
                  const biomeCounts: Record<string, number> = {};
                  cells.forEach(c => { biomeCounts[c.biome] = (biomeCounts[c.biome] || 0) + 1; });
                  const topBiomes = Object.entries(biomeCounts).sort((a, b) => b[1] - a[1]).slice(0, 5);
                  const totalPop = cells.reduce((s, c) => s + (c.population ?? 0), 0);
                  const popStr = totalPop >= 1e6 ? `${(totalPop/1e6).toFixed(1)}M` : totalPop >= 1e3 ? `${(totalPop/1e3).toFixed(0)}K` : String(totalPop);
                  return (
                    <div className="mt-2 space-y-1 text-[10px] text-gray-400">
                      <div className="flex justify-between"><span>Land coverage</span><span className="text-gray-200">{landPct}%</span></div>
                      <div className="flex justify-between"><span>Total cells</span><span className="text-gray-200">{cells.length.toLocaleString()}</span></div>
                      {totalPop > 0 && <div className="flex justify-between"><span>Population</span><span className="text-gray-200">{popStr}</span></div>}
                      {worldData.civData && <div className="flex justify-between"><span>Factions</span><span className="text-gray-200">{worldData.civData.factions.length}</span></div>}
                      {worldData.rivers && <div className="flex justify-between"><span>Rivers</span><span className="text-gray-200">{worldData.rivers.length}</span></div>}
                      <div className="mt-1 pt-1 border-t border-gray-800">
                        {topBiomes.map(([biome, count]) => (
                          <div key={biome} className="flex justify-between">
                            <span className="truncate">{biome}</span>
                            <span className="text-gray-200 ml-2">{((count/cells.length)*100).toFixed(1)}%</span>
                          </div>
                        ))}
                      </div>
                    </div>
                  );
                })()}
              </div>
            )}
          </div>
        )}

        {/* ... (Other Tabs omitted for brevity, logic remains identical to previous) ... */}
        {/* Keeping existing Geo, Climate, Political, Export tabs rendering logic as is */}
        {activeTab === 'geo' && (
           <div className="space-y-5">
              <div className="space-y-1">
                 <label className="text-xs text-gray-400 block mb-1">Terrain Preset</label>
                 <select 
                    value={params.landStyle}
                    onChange={(e) => { handlePresetChange(e.target.value as LandStyle); }}
                    className="w-full bg-gray-800 text-white text-xs border border-gray-700 p-2"
                 >
                    <option value="Continents">Continents</option>
                    <option value="Pangea">Pangea</option>
                    <option value="Archipelago">Archipelago</option>
                    <option value="Islands">Islands</option>
                    <option value="Custom">Custom</option>
                 </select>
              </div>
              <div className="space-y-1">
                <div className="flex justify-between text-xs text-gray-400">
                  <label>Sea Level</label>
                  <span>{(params.seaLevel * 100).toFixed(0)}%</span>
                </div>
                <input
                  type="range" min="0.1" max="0.9" step="0.05"
                  value={params.seaLevel}
                  onChange={(e) => { handleChange('seaLevel', parseFloat(e.target.value)); }}
                  className="w-full h-1 bg-gray-700 appearance-none cursor-pointer accent-blue-500"
                />
              </div>
               <div className="space-y-1">
                  <div className="flex justify-between text-xs text-gray-400">
                    <label>Planet Radius</label>
                    <span>{params.planetRadius} km</span>
                  </div>
                  <input
                    type="range" min="1000" max="20000" step="100"
                    value={params.planetRadius || 6371}
                    onChange={(e) => { handleChange('planetRadius', parseInt(e.target.value) as 1 | 2 | 3); }}
                    className="w-full h-1 bg-gray-700 appearance-none cursor-pointer accent-indigo-500"
                  />
              </div>
              <div className="space-y-1">
                <div className="flex justify-between text-xs text-gray-400">
                  <label>Tectonic Plates</label>
                  <span>{params.plates}</span>
                </div>
                <input
                  type="range" min="2" max="50" step="1"
                  value={params.plates}
                  onChange={(e) => { handleAdvancedChange('plates', parseInt(e.target.value) as 1 | 2 | 3); }}
                  className="w-full h-1 bg-gray-700 appearance-none cursor-pointer accent-rose-500"
                />
              </div>
              <div className="space-y-1">
                <div className="flex justify-between text-xs text-gray-400">
                  <label>Terrain Roughness</label>
                  <span>{(params.roughness * 100).toFixed(0)}%</span>
                </div>
                <input
                  type="range" min="0" max="1" step="0.1"
                  value={params.roughness}
                  onChange={(e) => { handleAdvancedChange('roughness', parseFloat(e.target.value)); }}
                  className="w-full h-1 bg-gray-700 appearance-none cursor-pointer accent-slate-400"
                />
              </div>
              <div className="space-y-1" title="FBM octave count for structural terrain noise. More octaves = finer nested detail; fewer = smoother, broader forms.">
                <div className="flex justify-between text-xs text-gray-400">
                  <label>Detail Octaves</label>
                  <span>{params.detailLevel}</span>
                </div>
                <input
                  type="range" min="1" max="6" step="1"
                  value={params.detailLevel}
                  onChange={(e) => { handleAdvancedChange('detailLevel', parseInt(e.target.value, 10)); }}
                  className="w-full h-1 bg-gray-700 appearance-none cursor-pointer accent-emerald-500"
                />
              </div>
              <div className="space-y-1" title="Randomizes the cell grid. 0 = regular Fibonacci lattice; 1 = fully jittered organic cells.">
                <div className="flex justify-between text-xs text-gray-400">
                  <label>Cell Jitter</label>
                  <span>{(params.cellJitter * 100).toFixed(0)}%</span>
                </div>
                <input
                  type="range" min="0" max="1" step="0.05"
                  value={params.cellJitter}
                  onChange={(e) => { handleAdvancedChange('cellJitter', parseFloat(e.target.value)); }}
                  className="w-full h-1 bg-gray-700 appearance-none cursor-pointer accent-lime-500"
                />
              </div>
              <div className="space-y-1" title="Controls terrain feature size. Lower = broader continents; higher = more fragmented detail.">
                <div className="flex justify-between text-xs text-gray-400">
                  <label>Feature Frequency</label>
                  <span>{params.noiseScale.toFixed(1)}</span>
                </div>
                <input
                  type="range" min="0.1" max="5.0" step="0.1"
                  value={params.noiseScale}
                  onChange={(e) => { handleAdvancedChange('noiseScale', parseFloat(e.target.value)); }}
                  className="w-full h-1 bg-gray-700 appearance-none cursor-pointer accent-green-500"
                />
              </div>
              <div className="space-y-1" title="0 = smooth rounded hills (FBM). 1 = sharp jagged mountain ridges (ridged noise).">
                <div className="flex justify-between text-xs text-gray-400">
                  <label>Ridge Intensity</label>
                  <span>{(params.ridgeBlend * 100).toFixed(0)}%</span>
                </div>
                <input
                  type="range" min="0" max="1" step="0.1"
                  value={params.ridgeBlend}
                  onChange={(e) => { handleAdvancedChange('ridgeBlend', parseFloat(e.target.value)); }}
                  className="w-full h-1 bg-gray-700 appearance-none cursor-pointer accent-yellow-500"
                />
              </div>
              <div className="space-y-1" title="Amplifies terrain above sea level using a power curve. >1.0 = taller peaks; <1.0 = flatter land.">
                <div className="flex justify-between text-xs text-gray-400">
                  <label>Mountain Heights</label>
                  <span>{params.mountainHeight.toFixed(2)}</span>
                </div>
                <input
                  type="range" min="0.5" max="2.0" step="0.05"
                  value={params.mountainHeight}
                  onChange={(e) => { handleAdvancedChange('mountainHeight', parseFloat(e.target.value)); }}
                  className="w-full h-1 bg-gray-700 appearance-none cursor-pointer accent-yellow-300"
                />
              </div>
              <div className="space-y-1" title="Amplifies ocean depth below sea level using a power curve. >1.0 = deeper trenches; <1.0 = shallower ocean floor.">
                <div className="flex justify-between text-xs text-gray-400">
                  <label>Sea / Trench Depth</label>
                  <span>{params.oceanDepth.toFixed(2)}</span>
                </div>
                <input
                  type="range" min="0.5" max="2.0" step="0.05"
                  value={params.oceanDepth}
                  onChange={(e) => { handleAdvancedChange('oceanDepth', parseFloat(e.target.value)); }}
                  className="w-full h-1 bg-gray-700 appearance-none cursor-pointer accent-amber-600"
                />
              </div>
              <div className="space-y-1" title="Domain warping — twists terrain shapes for more organic, swirled coastlines and mountain ranges.">
                <div className="flex justify-between text-xs text-gray-400">
                  <label>Swirl / Warp</label>
                  <span>{params.warpStrength.toFixed(1)}</span>
                </div>
                <input
                  type="range" min="0" max="2.0" step="0.1"
                  value={params.warpStrength}
                  onChange={(e) => { handleAdvancedChange('warpStrength', parseFloat(e.target.value)); }}
                  className="w-full h-1 bg-gray-700 appearance-none cursor-pointer accent-purple-500"
                />
              </div>
               <div className="space-y-1" title="How strongly tectonic plate boundaries shape mountain ranges and rifts. Capped at 1.0 internally.">
                <div className="flex justify-between text-xs text-gray-400">
                  <label>Tectonic Strength</label>
                  <span>{params.plateInfluence.toFixed(2)}</span>
                </div>
                <input
                  type="range" min="0.1" max="1.0" step="0.05"
                  value={params.plateInfluence}
                  onChange={(e) => { handleAdvancedChange('plateInfluence', parseFloat(e.target.value)); }}
                  className="w-full h-1 bg-gray-700 appearance-none cursor-pointer accent-red-500"
                />
              </div>
              <div className="space-y-1">
                <div className="flex justify-between text-xs text-gray-400">
                  <label>Hydraulic Erosion</label>
                  <span>{params.erosionIterations} Steps</span>
                </div>
                <input
                  type="range" min="0" max="50" step="1"
                  value={params.erosionIterations}
                  onChange={(e) => { handleAdvancedChange('erosionIterations', parseInt(e.target.value) as 1 | 2 | 3); }}
                  className="w-full h-1 bg-gray-700 appearance-none cursor-pointer accent-stone-400"
                />
              </div>
           </div>
        )}

        {/* Climate Tab content */}
        {activeTab === 'climate' && (
           <div className="space-y-5">
              <div className="space-y-1">
                <div className="flex justify-between text-xs text-gray-400">
                  <label>Axial Tilt (Visual & Climatic)</label>
                  <span>{params.axialTilt || 0}°</span>
                </div>
                <input
                  type="range" min="-90" max="90" step="1"
                  value={params.axialTilt || 0}
                  onChange={(e) => { handleChange('axialTilt', parseFloat(e.target.value)); }}
                  className="w-full h-1 bg-gray-700 appearance-none cursor-pointer accent-indigo-400"
                />
             </div>
             {/* ... (rest of climate sliders) ... */}
              <div className="space-y-1">
                <div className="flex justify-between text-xs text-gray-400">
                  <label>Equator Temp (°C)</label>
                  <span>{params.baseTemperature}°C</span>
                </div>
                <input
                  type="range" min="-10" max="50" step="1"
                  value={params.baseTemperature}
                  onChange={(e) => { handleChange('baseTemperature', parseFloat(e.target.value)); }}
                  className="w-full h-1 bg-gray-700 appearance-none cursor-pointer accent-orange-500"
                />
              </div>
              <div className="space-y-1">
                <div className="flex justify-between text-xs text-gray-400">
                  <label>Pole Temp (°C)</label>
                  <span>{params.poleTemperature}°C</span>
                </div>
                <input
                  type="range" min="-50" max="20" step="1"
                  value={params.poleTemperature}
                  onChange={(e) => { handleChange('poleTemperature', parseFloat(e.target.value)); }}
                  className="w-full h-1 bg-gray-700 appearance-none cursor-pointer accent-sky-500"
                />
              </div>
              <div className="space-y-1">
                <div className="flex justify-between text-xs text-gray-400">
                  <label>Global Rainfall</label>
                  <span>{params.rainfallMultiplier.toFixed(1)}x</span>
                </div>
                <input
                  type="range" min="0" max="3" step="0.1"
                  value={params.rainfallMultiplier}
                  onChange={(e) => { handleChange('rainfallMultiplier', parseFloat(e.target.value)); }}
                  className="w-full h-1 bg-gray-700 appearance-none cursor-pointer accent-cyan-500"
                />
              </div>
              <div className="space-y-1" title="How far wind carries moisture inland before it dissipates. Higher = wetter interiors; lower = stronger rain shadows.">
                <div className="flex justify-between text-xs text-gray-400">
                  <label>Wind Strength / Moisture Transport</label>
                  <span>{(params.moistureTransport * 100).toFixed(0)}%</span>
                </div>
                <input
                  type="range" min="0" max="1" step="0.1"
                  value={params.moistureTransport}
                  onChange={(e) => { handleChange('moistureTransport', parseFloat(e.target.value)); }}
                  className="w-full h-1 bg-gray-700 appearance-none cursor-pointer accent-blue-300"
                />
                <p className="text-[9px] text-gray-500">Affects rain shadows & moisture spread</p>
              </div>
              <div className="space-y-1" title="Adds simplex noise to temperature — creates local hot/cold anomalies beyond the baseline latitude gradient.">
                <div className="flex justify-between text-xs text-gray-400">
                  <label>Random Temp</label>
                  <span>{params.temperatureVariance}</span>
                </div>
                <input
                  type="range" min="0" max="20" step="1"
                  value={params.temperatureVariance}
                  onChange={(e) => { handleChange('temperatureVariance', parseFloat(e.target.value)); }}
                  className="w-full h-1 bg-gray-700 appearance-none cursor-pointer accent-red-400"
                />
              </div>
           </div>
        )}

        {/* Civ Tab content */}
        {activeTab === 'political' && (
           <div className="space-y-5">
              {/* Civ Seed Input */}
             <div className="bg-gray-900 p-3 border border-gray-800">
                <label className="text-xs text-gray-400 mb-1 block">Civ Seed</label>
                <div className="flex gap-2">
                   <input 
                      type="text" 
                      value={params.civSeed} 
                      onChange={(e) => { handleChange('civSeed', e.target.value); }}
                      disabled={civSeedLocked}
                      className="bg-black border border-gray-700 px-2 py-1 text-white text-xs flex-1 disabled:opacity-50"
                   />
                   <button 
                      onClick={() => { setCivSeedLocked(!civSeedLocked); }} 
                      className={`${civSeedLocked ? 'text-blue-500' : 'text-gray-400'} hover:text-white transition-colors`}
                   >
                      {civSeedLocked ? <Lock size={14}/> : <Unlock size={14}/>}
                   </button>
                   <button onClick={handleRandomizeCivSeed} disabled={civSeedLocked} className="text-gray-400 hover:text-white disabled:opacity-50">
                      <Shuffle size={14} />
                   </button>
                </div>
             </div>

              {/* Generation Parameters — collapsible */}
              <div>
                <button
                  onClick={() => setShowCivParams(v => !v)}
                  className="flex items-center justify-between w-full text-[10px] font-semibold text-gray-500 uppercase tracking-wide hover:text-gray-300 transition-colors"
                >
                  <span>Generation Parameters</span>
                  {showCivParams ? <ChevronUp size={12}/> : <ChevronDown size={12}/>}
                </button>
                {showCivParams && (
                  <div className="mt-3 space-y-4">
                    <div className="grid grid-cols-2 gap-2">
                      <button
                        onClick={handleRerollBorders}
                        disabled={loading}
                        className="flex items-center justify-center gap-1 text-[10px] bg-blue-900/40 text-blue-300 px-2 py-2 border border-blue-900/50 hover:bg-blue-800/40"
                      >
                        <Shuffle size={10} /> Reroll Borders
                      </button>
                      <button
                        onClick={handleRerollProvinces}
                        disabled={loading}
                        className="flex items-center justify-center gap-1 text-[10px] bg-teal-900/40 text-teal-300 px-2 py-2 border border-teal-900/50 hover:bg-teal-800/40"
                      >
                        <Layers size={10} /> Reroll Provs
                      </button>
                    </div>

              <div className="space-y-1" title="Number of culture regions (C1). Factions inherit their naming style from the culture their capital falls in.">
                <div className="flex justify-between text-xs text-gray-400">
                  <label>Cultures</label>
                  <span>{params.numCultures}</span>
                </div>
                <input
                  type="range" min="2" max="8"
                  value={params.numCultures}
                  onChange={(e) => { handleChange('numCultures', parseInt(e.target.value) as 1 | 2 | 3); }}
                  className="w-full h-1 bg-gray-700 appearance-none cursor-pointer accent-pink-400"
                />
              </div>
              <div className="space-y-1">
                <div className="flex justify-between text-xs text-gray-400">
                  <label>Factions</label>
                  <span>{params.numFactions}</span>
                </div>
                <input
                  type="range" min="2" max="20"
                  value={params.numFactions}
                  onChange={(e) => { handleChange('numFactions', parseInt(e.target.value) as 1 | 2 | 3); }}
                  className="w-full h-1 bg-gray-700 appearance-none cursor-pointer accent-purple-500"
                />
              </div>
              <div className="space-y-1" title="Minimum angular separation between faction capitals. Higher = capitals spawn further apart, producing more evenly distributed territories.">
                <div className="flex justify-between text-xs text-gray-400">
                  <label>Capital Spacing</label>
                  <span>{(params.capitalSpacing * 100).toFixed(0)}%</span>
                </div>
                <input
                  type="range" min="0" max="1" step="0.1"
                  value={params.capitalSpacing}
                  onChange={(e) => { handleChange('capitalSpacing', parseFloat(e.target.value)); }}
                  className="w-full h-1 bg-gray-700 appearance-none cursor-pointer accent-purple-400"
                />
              </div>
              <div className="space-y-1">
                <div className="flex justify-between text-xs text-gray-400">
                  <label>Province Size (Admin Density)</label>
                  <span>{(params.provinceSize || 0.5).toFixed(2)}</span>
                </div>
                <input
                  type="range" min="0.1" max="1.0" step="0.1"
                  value={params.provinceSize || 0.5}
                  onChange={(e) => { handleChange('provinceSize', parseFloat(e.target.value)); }}
                  className="w-full h-1 bg-gray-700 appearance-none cursor-pointer accent-teal-400"
                />
              </div>
              <div className="space-y-1" title="How unequal faction sizes can be. 0 = all factions roughly equal; 1 = some factions much larger than others.">
                <div className="flex justify-between text-xs text-gray-400">
                  <label>Country Size Variance</label>
                  <span>{(params.civSizeVariance * 100).toFixed(0)}%</span>
                </div>
                <input
                  type="range" min="0" max="1" step="0.1"
                  value={params.civSizeVariance}
                  onChange={(e) => { handleChange('civSizeVariance', parseFloat(e.target.value)); }}
                  className="w-full h-1 bg-gray-700 appearance-none cursor-pointer accent-orange-400"
                />
              </div>
              <div className="space-y-1" title="How easily factions cross water. Higher = more seafaring civilisations that readily claim distant islands.">
                <div className="flex justify-between text-xs text-gray-400">
                  <label>Seafaring (Water Crossing Cost)</label>
                  <span>{(1.0 - params.waterCrossingCost).toFixed(1)}</span>
                </div>
                <input
                  type="range" min="0.1" max="1.0" step="0.1"
                  value={params.waterCrossingCost}
                  onChange={(e) => { handleChange('waterCrossingCost', parseFloat(e.target.value)); }}
                  className="w-full h-1 bg-gray-700 appearance-none cursor-pointer accent-blue-400"
                />
              </div>
              <div className="space-y-1" title="How far from coastline a faction claims ocean cells as territorial waters.">
                <div className="flex justify-between text-xs text-gray-400">
                  <label>Territorial Waters (Range)</label>
                  <span>{params.territorialWaters?.toFixed(2)}</span>
                </div>
                <input
                  type="range" min="0.01" max="1.0" step="0.01"
                  value={params.territorialWaters ?? 0.2}
                  onChange={(e) => { handleChange('territorialWaters', parseFloat(e.target.value)); }}
                  className="w-full h-1 bg-gray-700 appearance-none cursor-pointer accent-cyan-600"
                />
              </div>

                  </div>
                )}
              </div>

              {/* Factions editor */}
              {worldData?.civData && (
                <div className="space-y-2 pt-2 border-t border-gray-800">
                  <h3 className="text-[10px] font-semibold text-gray-500 uppercase tracking-wide pt-1">Factions</h3>
                  {worldData.civData.factions.map(f => {
                    const otherFactions = worldData.civData!.factions.filter(o => o.id !== f.id);
                    const mergeTarget = mergeTargets[f.id] ?? '';
                    return (
                    <div key={f.id} className="flex flex-col gap-1.5 bg-gray-900 p-2 border border-gray-800">
                      <div className="flex items-center gap-2">
                        <input
                          type="color"
                          value={f.color}
                          onChange={e => onEditFaction?.(f.id, { color: e.target.value })}
                          className="w-7 h-6 border border-gray-700 bg-transparent cursor-pointer flex-shrink-0"
                          title="Faction color"
                        />
                        <input
                          type="text"
                          value={f.name}
                          onChange={e => onEditFaction?.(f.id, { name: e.target.value })}
                          className="flex-1 bg-black border border-gray-700 px-2 py-1 text-white text-xs focus:outline-none focus:border-blue-500"
                          placeholder="Faction name"
                        />
                      </div>
                      {onMergeFactions && otherFactions.length > 0 && (
                        <div className="flex items-center gap-1.5 pl-1">
                          <select
                            value={mergeTarget}
                            onChange={e => {
                              const val = e.target.value === '' ? '' : parseInt(e.target.value, 10);
                              setMergeTargets(prev => ({ ...prev, [f.id]: val }));
                            }}
                            className="flex-1 bg-black border border-gray-700 px-1.5 py-1 text-gray-300 text-[10px] focus:outline-none focus:border-blue-500"
                          >
                            <option value="">Merge into…</option>
                            {otherFactions.map(o => (
                              <option key={o.id} value={o.id}>{o.name}</option>
                            ))}
                          </select>
                          <button
                            onClick={() => {
                              if (mergeTarget === '') return;
                              onMergeFactions(f.id, mergeTarget);
                              setMergeTargets(prev => ({ ...prev, [f.id]: '' }));
                            }}
                            disabled={mergeTarget === ''}
                            className="px-2 py-1 text-[10px] font-semibold uppercase tracking-wide bg-red-900/50 text-red-300 border border-red-800 hover:bg-red-900 disabled:opacity-40 disabled:cursor-not-allowed disabled:hover:bg-red-900/50"
                            title={`Merge ${f.name} into the selected faction`}
                          >
                            Merge
                          </button>
                        </div>
                      )}
                    </div>
                    );
                  })}
                </div>
              )}

              {/* Name Style */}
              <div className="space-y-1 border-t border-gray-800 pt-3">
                  <label className="text-xs text-gray-400 block mb-1">Name Style</label>
                  <select
                     value={params.nameStyle || 'fantasy'}
                     onChange={(e) => { handleChange('nameStyle', e.target.value as NameStyle); }}
                     className="w-full bg-gray-800 text-white text-xs border border-gray-700 p-2"
                  >
                     {NAME_STYLES.map(style => (
                       <option key={style} value={style}>{style.charAt(0).toUpperCase() + style.slice(1)}</option>
                     ))}
                  </select>
                  <p className="text-[9px] text-gray-500">Applies on Reroll Borders or regeneration.</p>
              </div>

              {/* Lore Level */}
              <div className="space-y-1 border-t border-gray-800 pt-3">
                  <label className="text-xs text-gray-400 block mb-1">Lore Generation Detail</label>
                  <select 
                     value={params.loreLevel || 1}
                     onChange={(e) => { handleChange('loreLevel', parseInt(e.target.value) as 1 | 2 | 3); }}
                     className="w-full bg-gray-800 text-white text-xs border border-gray-700 p-2"
                  >
                     <option value={1}>Level 1: Factions & Capitals</option>
                     <option value={2}>Level 2: Provinces & Towns</option>
                     <option value={3}>Level 3: Backstories (Slow)</option>
                  </select>
              </div>

              <div className="bg-gray-800/50 p-3 border border-gray-700 mt-4">
                <div className="flex justify-between items-center mb-2">
                  <h2 className="text-xs font-semibold text-gray-300">AI Lore</h2>
                   <button 
                    onClick={onGenerateLore}
                    disabled={generatingLore || !apiKey}
                    className={`text-[10px] px-2 py-1 transition-colors ${
 apiKey 
 ? 'bg-blue-900/50 text-blue-300 hover:bg-blue-900 border border-blue-800' 
 : 'bg-gray-800 text-gray-500 cursor-not-allowed border border-gray-700'
 }`}
                  >
                    {generatingLore ? 'Thinking...' : 'Generate'}
                  </button>
                </div>
                {!apiKey && (
                  <p className="text-[9px] text-yellow-500/80 bg-yellow-500/10 p-1.5 border border-yellow-500/20 mb-2">
                    Enter a Gemini API Key in the "Sys" tab to enable AI lore.
                  </p>
                )}
                {lore ? (
                  <div className="space-y-2">
                    <h3 className="font-bold text-white text-xs">{lore.name}</h3>
                    <p className="text-[10px] text-gray-400 max-h-32 overflow-y-auto">
                      {lore.description}
                    </p>
                    {worldData?.civData && (
                        <div className="space-y-1 mt-2">
                            {worldData.civData.factions.map(f => (
                                <div key={f.id} className="text-[10px] bg-gray-900 p-1 border border-gray-700">
                                    <div style={{color: f.color}} className="font-bold">{f.name}</div>
                                    <div className="text-gray-500 pl-1">Cap: {f.provinces.flatMap(p => p.towns).find(t => t.cellId === f.capitalId)?.name || 'Unknown'}</div>
                                    {f.description && <div className="text-gray-400 italic pl-1 mt-1 border-t border-gray-800 pt-1">{f.description}</div>}
                                </div>
                            ))}
                        </div>
                    )}
                  </div>
                ) : (
                  <p className="text-[10px] text-gray-600 italic">Generate a world first.</p>
                )}
              </div>
           </div>
        )}

        {/* Export Tab Content */}
        {activeTab === 'export' && (
            <div className="space-y-6">
                 {/* ... (export tab same as before) ... */}
                 <div className="space-y-2">
                    <h3 className="text-xs font-bold text-gray-400 uppercase tracking-wider mb-2">Image Export</h3>
                    
                    <div className="space-y-1">
                        <label className="text-xs text-gray-400">Resolution</label>
                        <select 
                            value={expRes} 
                            onChange={(e) => { setExpRes(parseInt(e.target.value) as 1 | 2 | 3 as ExportResolution); }}
                            className="w-full bg-gray-800 text-white text-xs border border-gray-700 p-2"
                        >
                            <option value="2048">2K (2048px)</option>
                            <option value="4096">4K (4096px)</option>
                            <option value="8192">8K (8192px)</option>
                        </select>
                    </div>

                    <div className="space-y-1">
                        <label className="text-xs text-gray-400">Projection</label>
                        <select 
                            value={expProj} 
                            onChange={(e) => { setExpProj(e.target.value as ProjectionType); }}
                            className="w-full bg-gray-800 text-white text-xs border border-gray-700 p-2"
                        >
                            <option value="equirectangular">Equirectangular</option>
                            <option value="mercator">Mercator</option>
                            <option value="winkeltripel">Winkel Tripel</option>
                            <option value="robinson">Robinson</option>
                            <option value="mollweide">Mollweide</option>
                            <option value="orthographic">Orthographic</option>
                            <option value="dymaxion">Dymaxion (Icosahedron) (Experimental)</option>
                        </select>
                    </div>

                    {expProj === 'dymaxion' && (
                        <div className="border border-gray-800 p-3 space-y-3 bg-gray-900/40">
                            <div className="flex items-center justify-between">
                                <div className="text-xs font-semibold text-gray-300">Dymaxion Controls</div>
                                <label className="flex items-center gap-2 text-[10px] text-gray-400">
                                    <input
                                        type="checkbox"
                                        checked={dymaxionSettings.showOverlay}
                                        onChange={(e) => { updateDymaxion({ showOverlay: e.target.checked }); }}
                                        className="accent-blue-500"
                                    />
                                    Show Overlay
                                </label>
                            </div>

                            <div className="space-y-1">
                                <label className="text-xs text-gray-400">Manipulation Mode</label>
                                <div className="grid grid-cols-2 gap-2">
                                    <button
                                        onClick={() => { updateDymaxion({ mode: 'planet' as DymaxionControlMode }); }}
                                        className={`text-[10px] py-2 border ${dymaxionSettings.mode === 'planet' ? 'bg-blue-700/70 border-blue-500 text-white' : 'bg-gray-800 border-gray-700 text-gray-300'}`}
                                    >
                                        Rotate Planet
                                    </button>
                                    <button
                                        onClick={() => { updateDymaxion({ mode: 'overlay' as DymaxionControlMode }); }}
                                        className={`text-[10px] py-2 border ${dymaxionSettings.mode === 'overlay' ? 'bg-blue-700/70 border-blue-500 text-white' : 'bg-gray-800 border-gray-700 text-gray-300'}`}
                                    >
                                        Rotate Overlay
                                    </button>
                                </div>
                                <div className="text-[10px] text-gray-500">
                                    Drag the globe to rotate. Hold Shift while dragging to roll the overlay.
                                </div>
                            </div>

                            <div className="space-y-2">
                                <div className="flex justify-between text-xs text-gray-400">
                                    <label>Longitude</label>
                                    <span>{dymaxionSettings.lon}°</span>
                                </div>
                                <input
                                    type="range" min="-180" max="180" step="1"
                                    value={dymaxionSettings.lon}
                                    onChange={function(e) { updateDymaxion({ lon: parseInt(e.target.value) }); }}
                                    className="w-full h-1 bg-gray-700 appearance-none cursor-pointer accent-blue-400"
                                />
                            </div>

                            <div className="space-y-2">
                                <div className="flex justify-between text-xs text-gray-400">
                                    <label>Latitude</label>
                                    <span>{dymaxionSettings.lat}°</span>
                                </div>
                                <input
                                    type="range" min="-90" max="90" step="1"
                                    value={dymaxionSettings.lat}
                                    onChange={function(e) { updateDymaxion({ lat: parseInt(e.target.value) }); }}
                                    className="w-full h-1 bg-gray-700 appearance-none cursor-pointer accent-blue-400"
                                />
                            </div>

                            <div className="space-y-2">
                                <div className="flex justify-between text-xs text-gray-400">
                                    <label>Roll</label>
                                    <span>{dymaxionSettings.roll}°</span>
                                </div>
                                <input
                                    type="range" min="-180" max="180" step="1"
                                    value={dymaxionSettings.roll}
                                    onChange={function(e) { updateDymaxion({ roll: parseInt(e.target.value) }); }}
                                    className="w-full h-1 bg-gray-700 appearance-none cursor-pointer accent-blue-400"
                                />
                            </div>

                            <label className="flex items-center gap-2 text-[10px] text-gray-400">
                                <input
                                    type="checkbox"
                                    checked={dymaxionSettings.layout === 'blender'}
                                    onChange={(e) => { updateDymaxion({ layout: e.target.checked ? 'blender' : 'classic' }); }}
                                    className="accent-blue-500"
                                />
                                Blender UV Net (export only)
                            </label>

                            <div className="space-y-1">
                                <label className="text-[10px] text-gray-500">Orientation Presets</label>
                                <div className="grid grid-cols-4 gap-1">
                                    {([['N.Pole', 0, -90, 0], ['Pacific', -150, 0, 0], ['Atlantic', 0, 0, 0], ['Asia', 90, 0, 0]] as const).map(([label, lon, lat, roll]) => (
                                        <button
                                            key={label}
                                            onClick={() => { updateDymaxion({ lon, lat, roll }); }}
                                            className="text-[9px] bg-gray-800 hover:bg-gray-700 text-gray-300 py-1.5 border border-gray-700"
                                        >
                                            {label}
                                        </button>
                                    ))}
                                </div>
                            </div>

                            <button
                                onClick={() => { updateDymaxion({ lon: 0, lat: 0, roll: 0 }); }}
                                className="w-full text-[10px] bg-gray-800 hover:bg-gray-700 text-gray-200 py-2 border border-gray-700"
                            >
                                Reset Orientation
                            </button>

                            <label className="flex items-center gap-2 text-[10px] text-gray-400">
                                <input
                                    type="checkbox"
                                    checked={showDymaxion2D}
                                    onChange={(e) => { setShowDymaxion2D(e.target.checked); }}
                                    className="accent-blue-500"
                                />
                                Show 2D Preview
                            </label>

                            {showDymaxion2D && (
                                <DymaxionPreview2D
                                    world={worldData}
                                    viewMode={viewMode}
                                    settings={dymaxionSettings}
                                    onChange={onDymaxionChange}
                                />
                            )}
                        </div>
                    )}

                    <button
                        onClick={() => { void handleExport(); }}
                        disabled={!worldData}
                        className="w-full flex items-center justify-center gap-2 bg-green-700 hover:bg-green-600 text-white py-2 text-xs mt-2 disabled:opacity-50 border border-green-600"
                    >
                        <Image size={14}/> Download PNG
                    </button>
                    <button
                        onClick={() => { if (worldData) void exportMap(worldData, 'height_bw', expRes, 'equirectangular'); }}
                        disabled={!worldData}
                        className="w-full flex items-center justify-center gap-2 bg-gray-700 hover:bg-gray-600 text-white py-2 text-xs disabled:opacity-50 border border-gray-600"
                    >
                        <Mountain size={14}/> Export Heightmap (BW)
                    </button>
                </div>

                <div className="border-t border-gray-800 pt-4 space-y-2">
                    <h3 className="text-xs font-bold text-gray-400 uppercase tracking-wider mb-2">Vector Export</h3>
                    <p className="text-[10px] text-gray-500">Editable coastlines, borders, rivers, and labels for Inkscape/Illustrator, or geodesic GeoJSON for QGIS/web-GIS.</p>
                    <button
                        onClick={() => { if (worldData && expProj !== 'dymaxion') downloadSVG(worldData, viewMode, expProj); }}
                        disabled={!worldData || expProj === 'dymaxion'}
                        title={expProj === 'dymaxion' ? 'SVG export is raster-only for Dymaxion — choose another projection' : undefined}
                        className="w-full flex items-center justify-center gap-2 bg-teal-700 hover:bg-teal-600 text-white py-2 text-xs disabled:opacity-50 border border-teal-600"
                    >
                        <FileCode size={14}/> Download SVG
                    </button>
                    <button
                        onClick={() => { if (worldData) downloadGeoJSON(worldData); }}
                        disabled={!worldData}
                        className="w-full flex items-center justify-center gap-2 bg-teal-700 hover:bg-teal-600 text-white py-2 text-xs disabled:opacity-50 border border-teal-600"
                    >
                        <FileJson size={14}/> Download GeoJSON
                    </button>
                </div>

                <div className="border-t border-gray-800 pt-4 space-y-2">
                    <h3 className="text-xs font-bold text-gray-400 uppercase tracking-wider mb-2">3D Export</h3>
                    <p className="text-[10px] text-gray-500">Exports the current view as a GLB file. World mesh uses per-vertex colors. Rivers exported as line geometry. City markers included when civilization data is present.</p>
                    <button
                        onClick={() => { if (worldData) exportGLB(worldData, viewMode); }}
                        disabled={!worldData}
                        className="w-full flex items-center justify-center gap-2 bg-indigo-700 hover:bg-indigo-600 text-white py-2 text-xs disabled:opacity-50 border border-indigo-600"
                    >
                        <Box size={14}/> Export GLB
                    </button>
                </div>

                <div className="border-t border-gray-800 pt-4 space-y-3">
                    <h3 className="text-xs font-bold text-gray-400 uppercase tracking-wider">File Management</h3>
                    
                    <button
                        onClick={() => { if (params) { void saveMapConfig(params, worldData || undefined); } }}
                        className="w-full flex items-center justify-center gap-2 bg-gray-800 hover:bg-gray-700 text-white py-2 text-xs border border-gray-700"
                    >
                        <Save size={14} /> Save Config (JSON)
                    </button>
                    
                    <div className="relative">
                        <input 
                            type="file" 
                            accept=".json" 
                            onChange={(e) => { void handleFileUpload(e); }}
                            className="absolute inset-0 w-full h-full opacity-0 cursor-pointer"
                        />
                        <button className="w-full flex items-center justify-center gap-2 bg-gray-800 hover:bg-gray-700 text-white py-2 text-xs pointer-events-none border border-gray-700">
                            <FolderOpen size={14} /> Load Config (JSON)
                        </button>
                    </div>
                </div>

                <div className="border-t border-gray-800 pt-4 space-y-3">
                    <h3 className="text-xs font-bold text-gray-400 uppercase tracking-wider">Browser Storage</h3>
                    
                    <div className="flex gap-2">
                        <input 
                            type="text" 
                            placeholder="Save Name..." 
                            value={saveName}
                            onChange={(e) => { setSaveName(e.target.value); }}
                            className="flex-1 bg-gray-900 border border-gray-700 px-2 text-xs text-white"
                        />
                        <button 
                            onClick={handleSaveBrowser}
                            disabled={!saveName}
                            className="bg-blue-600 hover:bg-blue-500 text-white px-3"
                        >
                            <Save size={14}/>
                        </button>
                    </div>

                    <div className="space-y-1 max-h-40 overflow-y-auto">
                        {savedMaps.length === 0 && <p className="text-xs text-gray-600 italic">No saved maps.</p>}
                        {savedMaps.map(entry => (
                            <div key={entry.name} className="flex items-center justify-between bg-gray-900 p-2 border border-gray-800 group">
                                <div className="flex flex-col">
                                    <span className="text-xs font-bold text-gray-300">{entry.name}</span>
                                    <span className="text-[10px] text-gray-500">{new Date(entry.date).toLocaleDateString()}</span>
                                </div>
                                <div className="flex gap-1 opacity-0 group-hover:opacity-100 transition-opacity">
                                    <button onClick={() => { handleLoadBrowser(entry.params, entry.civData, entry.markers); }} className="text-blue-400 hover:text-white p-1"><FolderOpen size={12}/></button>
                                    <button onClick={() => { handleDeleteBrowser(entry.name); }} className="text-red-400 hover:text-white p-1"><Trash2 size={12}/></button>
                                </div>
                            </div>
                        ))}
                    </div>
                </div>
            </div>
        )}
      </div>

      <div className="p-4 border-t border-gray-800 space-y-2">
         {/* Console Output area */}
         <div className="mb-2">
             <div 
               className="flex items-center justify-between text-xs text-gray-500 mb-1 cursor-pointer hover:text-gray-300"
               onClick={() => { setConsoleOpen(!consoleOpen); }}
             >
                 <div className="flex items-center gap-1">
                    <Terminal size={10} />
                    <span>System Console</span>
                 </div>
                 {consoleOpen ? <ChevronDown size={10}/> : <ChevronUp size={10}/>}
             </div>
             <ConsoleOutput logs={logs} isOpen={consoleOpen} />
         </div>

         {/* Generation progress bar */}
         {loading && (
           <div className="w-full h-1 bg-gray-800 rounded overflow-hidden">
             <div
               className="h-full bg-blue-500 transition-all duration-300"
               style={{ width: `${genProgress * 100}%` }}
             />
           </div>
         )}

         {!loading ? (
             <button
              onClick={handleGenerateClick}
              className={`w-full py-3 font-semibold flex items-center justify-center gap-2 transition-all relative overflow-hidden bg-blue-600 hover:bg-blue-500 text-white border border-blue-500`}
            >
              <div className="relative flex items-center gap-2 z-10">
                  <RefreshCw size={16} />
                  Generate World
              </div>
            </button>
         ) : (
            <button
              onClick={onCancel}
              className={`w-full py-3 font-semibold flex items-center justify-center gap-2 transition-all relative overflow-hidden bg-red-600 hover:bg-red-500 text-white border border-red-500`}
            >
              <div className="relative flex items-center gap-2 z-10">
                  <XCircle size={16} />
                  Cancel Generation
              </div>
            </button>
         )}
      </div>
    </div>
  );
};

export default Controls;
