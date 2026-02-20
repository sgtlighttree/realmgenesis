import React, { useState, useEffect, useRef } from 'react';
import { WorldParams, ViewMode, LoreData, LandStyle, CivData, DisplayMode, DymaxionSettings, DymaxionLayout, DymaxionControlMode } from '../types';
import { RefreshCw, Globe, Thermometer, Droplets, Flag, Mountain, Lock, Unlock, Shuffle, Eye, Layers, Zap, Grid, Download, Save, FileJson, FolderOpen, Trash2, Image, Satellite, Waves, Terminal, XCircle, ChevronDown, ChevronUp } from 'lucide-react';
import { exportMap, saveMapConfig, loadMapConfig, saveMapToBrowser, getSavedMaps, deleteSavedMap, ExportResolution, ProjectionType } from '../utils/export';
import { WorldData } from '../types';
import DymaxionPreview2D from './DymaxionPreview2D';

interface ControlsProps {
  params: WorldParams;
  setParams: (p: WorldParams) => void;
  onGenerate: (p?: WorldParams) => void;
  onLoadWorld: (p: WorldParams, civData?: CivData) => void;
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
  showGrid: boolean;
  setShowGrid: (b: boolean) => void;
  showRivers: boolean;
  setShowRivers: (b: boolean) => void;
  dymaxionSettings: DymaxionSettings;
  onDymaxionChange: React.Dispatch<React.SetStateAction<DymaxionSettings>>;
  apiKey: string;
  onApiKeyChange: (key: string) => void;
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
        <div className="bg-black border border-gray-800 rounded-b-lg p-2 h-32 overflow-y-auto font-mono text-[10px] space-y-1 shadow-inner relative transition-all">
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
  showGrid,
  setShowGrid,
  showRivers,
  setShowRivers,
  dymaxionSettings,
  onDymaxionChange,
  apiKey,
  onApiKeyChange
}) => {
  const [activeTab, setActiveTab] = useState<Tab>('system');
  const [seedLocked, setSeedLocked] = useState(false);
  const [civSeedLocked, setCivSeedLocked] = useState(false);
  const [autoUpdate, setAutoUpdate] = useState(true);
  const [consoleOpen, setConsoleOpen] = useState(true);
  
  // Export State
  const [expRes, setExpRes] = useState<ExportResolution>(4096);
  const [expProj, setExpProj] = useState<ProjectionType>('equirectangular');
  const [saveName, setSaveName] = useState('');
  const [savedMaps, setSavedMaps] = useState(getSavedMaps());
  const [showDymaxion2D, setShowDymaxion2D] = useState(false);

  const updateDymaxion = (patch: Partial<DymaxionSettings>) => {
      onDymaxionChange((prev) => ({ ...prev, ...patch }));
  };

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

  const handleChange = <K extends keyof WorldParams>(key: K, value: WorldParams[K]) => { // eslint-disable-next-line @typescript-eslint/no-explicit-any
    setParams({ ...params, [key]: value });
  };

  const handleAdvancedChange = <K extends keyof WorldParams>(key: K, value: WorldParams[K]) => { // eslint-disable-next-line @typescript-eslint/no-explicit-any
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
              break;
          case 'Pangea':
              updates.noiseScale = 0.8;
              updates.ridgeBlend = 0.4;
              updates.maskType = 'Pangea';
              updates.warpStrength = 0.5;
              updates.plateInfluence = 0.6;
              updates.erosionIterations = 2;
              updates.cellJitter = 0.4;
              break;
          case 'Archipelago':
              updates.noiseScale = 2.5;
              updates.ridgeBlend = 0.8;
              updates.maskType = 'None';
              updates.warpStrength = 0.1;
              updates.plateInfluence = 1.2;
              updates.erosionIterations = 5;
              updates.cellJitter = 0.6;
              break;
          case 'Islands':
              updates.noiseScale = 2.0;
              updates.ridgeBlend = 0.3;
              updates.maskType = 'None';
              updates.warpStrength = 0.1;
              updates.plateInfluence = 0.8;
              updates.erosionIterations = 3;
              updates.cellJitter = 0.7;
              break;
          case 'Custom':
              break;
      }
      setParams({ ...params, ...updates });
  };

  const handleRandomizeSeed = () => {
    if (!seedLocked) {
      handleChange('seed', crypto.getRandomValues(new Uint32Array(1))[0].toString(36));
      if (!civSeedLocked) handleChange('civSeed', crypto.getRandomValues(new Uint32Array(1))[0].toString(36));
    }
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
          expProj === 'dymaxion' ? { layout: dymaxionSettings.layout, lon: dymaxionSettings.lon, lat: dymaxionSettings.lat, roll: dymaxionSettings.roll } : undefined
        );
    } catch(e) {
        console.error(e);
        alert("Export failed. Try a lower resolution.");
    }
  };

  const handleSaveBrowser = () => {
      if (!saveName) return;
      // Pass civData if available to save lore
      if (saveMapToBrowser(saveName, params, worldData?.civData)) {
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

  const handleLoadBrowser = (entryParams: WorldParams, civData?: CivData) => {
      if (confirm("Load this map configuration? Unsaved changes will be lost.")) {
          setParams(entryParams);
          // Use the dedicated load function that handles restoration
          setTimeout(() => onLoadWorld(entryParams, civData), 50);
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
              setTimeout(() => onLoadWorld(loaded.params, loaded.civData), 50);
          } else {
              alert("Invalid config file");
          }
      }
  };

  const ViewButton = ({ mode, icon: Icon, label }: { mode: ViewMode, icon: any, label: string }) => (
    <button
      onClick={() => { setViewMode(mode); }}
      className={`flex items-center gap-2 px-2 py-1.5 rounded-lg text-xs transition-all flex-1 justify-center ${
        viewMode === mode 
          ? 'bg-blue-600 text-white shadow-lg shadow-blue-900/50' 
          : 'bg-gray-800 text-gray-400 hover:bg-gray-700 hover:text-white'
      }`}
    >
      <Icon size={14} />
      {label}
    </button>
  );

  const DisplayButton = ({ mode, label }: { mode: DisplayMode; label: string }) => (
    <button
      onClick={() => { setDisplayMode(mode); }}
      className={`px-2 py-1.5 rounded-lg text-xs transition-all flex-1 ${
        displayMode === mode
          ? 'bg-blue-600 text-white shadow-lg shadow-blue-900/50'
          : 'bg-gray-800 text-gray-400 hover:bg-gray-700 hover:text-white'
      }`}
    >
      {label}
    </button>
  );

  return (
    <div className="w-80 bg-gray-950 border-r border-gray-800 flex flex-col h-full overflow-hidden text-sm relative z-20">
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
                 <DisplayButton mode="globe" label="3D Globe" />
                 <DisplayButton mode="mercator" label="2D Mercator" />
                 <DisplayButton mode="dymaxion" label="2D Dymaxion" />
               </div>
             </div>

             {/* Map Name Input */}
             <div className="space-y-1">
                 <label className="text-xs text-gray-400 block">Map Name</label>
                 <input 
                    type="text" 
                    value={params.mapName} 
                    onChange={(e) => { handleChange('mapName', e.target.value); }}
                    className="w-full bg-gray-900 border border-gray-700 rounded px-2 py-1.5 text-white text-xs"
                    placeholder="My World"
                 />
             </div>

             {/* Seed Input */}
             <div className="bg-gray-900 p-3 rounded-lg border border-gray-800">
                <label className="text-xs text-gray-400 mb-1 block">Seed</label>
                <div className="flex gap-2">
                   <input 
                      type="text" 
                      value={params.seed} 
                      onChange={(e) => { handleChange('seed', e.target.value); }}
                      disabled={seedLocked}
                      className="bg-black border border-gray-700 rounded px-2 py-1 text-white text-xs flex-1 disabled:opacity-50"
                   />
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
                    max="1000000"
                    step="1000"
                    value={params.points}
                    onChange={(e) => { handleChange('points', parseInt(e.target.value) as 1 | 2 | 3); }}
                    className="w-24 bg-gray-900 border border-gray-700 rounded px-1 py-0.5 text-right text-white text-xs"
                />
              </div>
              <input
                type="range"
                min="2000"
                max="200000"
                step="1000"
                value={Math.min(200000, params.points)}
                onChange={(e) => { handleChange('points', parseInt(e.target.value) as 1 | 2 | 3); }}
                className="w-full h-1 bg-gray-700 rounded-lg appearance-none cursor-pointer accent-blue-500"
              />
            </div>
             
             <div className="flex items-center justify-between text-xs text-gray-400 pt-2 border-t border-gray-800">
                 <div className="flex items-center gap-2">
                    <Grid size={12} className={showGrid ? "text-blue-400" : "text-gray-600"}/>
                    <label>Lat/Long Grid</label>
                 </div>
                 <input 
                    type="checkbox"
                    checked={showGrid}
                    onChange={(e) => { setShowGrid(e.target.checked); }}
                    className="rounded bg-gray-700"
                 />
            </div>

            <div className="flex items-center justify-between text-xs text-gray-400 pt-2">
                 <div className="flex items-center gap-2">
                    <Waves size={12} className={showRivers ? "text-blue-400" : "text-gray-600"}/>
                    <label>River Network</label>
                 </div>
                 <input 
                    type="checkbox"
                    checked={showRivers}
                    onChange={(e) => { setShowRivers(e.target.checked); }}
                    className="rounded bg-gray-700"
                 />
            </div>

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
                    className="rounded bg-gray-700"
                 />
            </div>
            
            <div className="pt-4 border-t border-gray-800">
              <h3 className="text-xs font-semibold text-gray-500 mb-2">View Layer</h3>
              <div className="grid grid-cols-2 gap-2">
                <ViewButton mode="biome" icon={Globe} label="Biomes" />
                <ViewButton mode="satellite" icon={Satellite} label="Satellite" />
                <ViewButton mode="height" icon={Mountain} label="Height" />
                <ViewButton mode="height_bw" icon={Eye} label="Height BW" />
                <ViewButton mode="temperature" icon={Thermometer} label="Temp" />
                <ViewButton mode="moisture" icon={Droplets} label="Rain" />
                <ViewButton mode="plates" icon={Layers} label="Plates" />
                <ViewButton mode="political" icon={Flag} label="Borders" />
              </div>
            </div>

            <div className="pt-4 border-t border-gray-800 space-y-3">
              <h3 className="text-xs font-semibold text-gray-500 mb-2">AI Settings (BYOK)</h3>
              <div className="bg-gray-900 p-3 rounded-lg border border-gray-800 space-y-2">
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
                  className="w-full bg-black border border-gray-700 rounded px-2 py-1.5 text-white text-xs"
                />
                <p className="text-[9px] text-gray-500 italic">
                  Key is stored ephemerally in memory and will be lost on refresh.
                </p>
              </div>
            </div>
          </div>
        )}

        {/* ... (Other Tabs omitted for brevity, logic remains identical to previous) ... */}
        {/* Keeping existing Geo, Climate, Political, Export tabs rendering logic as is */}
        {activeTab === 'geo' && (
           <div className="space-y-4">
              <div className="space-y-1">
                 <label className="text-xs text-gray-400 block mb-1">Terrain Preset</label>
                 <select 
                    value={params.landStyle}
                    onChange={(e) => { handlePresetChange(e.target.value as LandStyle); }}
                    className="w-full bg-gray-800 text-white text-xs border border-gray-700 rounded p-2"
                 >
                    <option value="Continents">Continents</option>
                    <option value="Pangea">Pangea</option>
                    <option value="Archipelago">Archipelago</option>
                    <option value="Islands">Islands</option>
                    <option value="Custom">Custom</option>
                 </select>
              </div>
              <div className="p-3 bg-gray-900 rounded border border-gray-800 space-y-3">
                  <h3 className="text-[10px] font-bold text-gray-400 uppercase tracking-wider">Advanced Terrain</h3>
                  
                  <div className="space-y-1">
                    <div className="flex justify-between text-xs text-gray-400">
                      <label>Sea Level</label>
                      <span>{(params.seaLevel * 100).toFixed(0)}%</span>
                    </div>
                    <input
                      type="range" min="0.1" max="0.9" step="0.05"
                      value={params.seaLevel}
                      onChange={(e) => { handleChange('seaLevel', parseFloat(e.target.value)); }}
                      className="w-full h-1 bg-gray-700 rounded-lg appearance-none cursor-pointer accent-blue-500"
                    />
                  </div>
                  {/* ... other geo sliders ... */}
                   <div className="space-y-1">
                      <div className="flex justify-between text-xs text-gray-400">
                        <label>Planet Radius</label>
                        <span>{params.planetRadius} km</span>
                      </div>
                      <input
                        type="range" min="1000" max="20000" step="100"
                        value={params.planetRadius || 6371}
                        onChange={(e) => { handleChange('planetRadius', parseInt(e.target.value) as 1 | 2 | 3); }}
                        className="w-full h-1 bg-gray-700 rounded-lg appearance-none cursor-pointer accent-indigo-500"
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
                      className="w-full h-1 bg-gray-700 rounded-lg appearance-none cursor-pointer accent-rose-500"
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
                      className="w-full h-1 bg-gray-700 rounded-lg appearance-none cursor-pointer accent-slate-400"
                    />
                  </div>
                  <div className="space-y-1">
                    <div className="flex justify-between text-xs text-gray-400">
                      <label>Feature Frequency</label>
                      <span>{params.noiseScale.toFixed(1)}</span>
                    </div>
                    <input
                      type="range" min="0.1" max="5.0" step="0.1"
                      value={params.noiseScale}
                      onChange={(e) => { handleAdvancedChange('noiseScale', parseFloat(e.target.value)); }}
                      className="w-full h-1 bg-gray-700 rounded-lg appearance-none cursor-pointer accent-green-500"
                    />
                  </div>
                  <div className="space-y-1">
                    <div className="flex justify-between text-xs text-gray-400">
                      <label>Ridge Intensity</label>
                      <span>{(params.ridgeBlend * 100).toFixed(0)}%</span>
                    </div>
                    <input
                      type="range" min="0" max="1" step="0.1"
                      value={params.ridgeBlend}
                      onChange={(e) => { handleAdvancedChange('ridgeBlend', parseFloat(e.target.value)); }}
                      className="w-full h-1 bg-gray-700 rounded-lg appearance-none cursor-pointer accent-yellow-500"
                    />
                  </div>
                  <div className="space-y-1">
                    <div className="flex justify-between text-xs text-gray-400">
                      <label>Swirl / Warp</label>
                      <span>{params.warpStrength.toFixed(1)}</span>
                    </div>
                    <input
                      type="range" min="0" max="2.0" step="0.1"
                      value={params.warpStrength}
                      onChange={(e) => { handleAdvancedChange('warpStrength', parseFloat(e.target.value)); }}
                      className="w-full h-1 bg-gray-700 rounded-lg appearance-none cursor-pointer accent-purple-500"
                    />
                  </div>
                   <div className="space-y-1">
                    <div className="flex justify-between text-xs text-gray-400">
                      <label>Tectonic Strength</label>
                      <span>{params.plateInfluence.toFixed(1)}x</span>
                    </div>
                    <input
                      type="range" min="0" max="2.0" step="0.1"
                      value={params.plateInfluence}
                      onChange={(e) => { handleAdvancedChange('plateInfluence', parseFloat(e.target.value)); }}
                      className="w-full h-1 bg-gray-700 rounded-lg appearance-none cursor-pointer accent-red-500"
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
                      className="w-full h-1 bg-gray-700 rounded-lg appearance-none cursor-pointer accent-stone-400"
                    />
                  </div>
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
                  className="w-full h-1 bg-gray-700 rounded-lg appearance-none cursor-pointer accent-indigo-400"
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
                  className="w-full h-1 bg-gray-700 rounded-lg appearance-none cursor-pointer accent-orange-500"
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
                  className="w-full h-1 bg-gray-700 rounded-lg appearance-none cursor-pointer accent-sky-500"
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
                  className="w-full h-1 bg-gray-700 rounded-lg appearance-none cursor-pointer accent-cyan-500"
                />
              </div>
              <div className="space-y-1">
                <div className="flex justify-between text-xs text-gray-400">
                  <label>Wind Strength / Moisture Transport</label>
                  <span>{(params.moistureTransport * 100).toFixed(0)}%</span>
                </div>
                <input
                  type="range" min="0" max="1" step="0.1"
                  value={params.moistureTransport}
                  onChange={(e) => { handleChange('moistureTransport', parseFloat(e.target.value)); }}
                  className="w-full h-1 bg-gray-700 rounded-lg appearance-none cursor-pointer accent-blue-300"
                />
                <p className="text-[9px] text-gray-500">Affects rain shadows & moisture spread</p>
              </div>
              <div className="space-y-1">
                <div className="flex justify-between text-xs text-gray-400">
                  <label>Random Temp</label>
                  <span>{params.temperatureVariance}</span>
                </div>
                <input
                  type="range" min="0" max="20" step="1"
                  value={params.temperatureVariance}
                  onChange={(e) => { handleChange('temperatureVariance', parseFloat(e.target.value)); }}
                  className="w-full h-1 bg-gray-700 rounded-lg appearance-none cursor-pointer accent-red-400"
                />
              </div>
           </div>
        )}

        {/* Civ Tab content */}
        {activeTab === 'political' && (
           <div className="space-y-5">
              {/* Civ Seed Input */}
             <div className="bg-gray-900 p-3 rounded-lg border border-gray-800">
                <label className="text-xs text-gray-400 mb-1 block">Civ Seed</label>
                <div className="flex gap-2">
                   <input 
                      type="text" 
                      value={params.civSeed} 
                      onChange={(e) => { handleChange('civSeed', e.target.value); }}
                      disabled={civSeedLocked}
                      className="bg-black border border-gray-700 rounded px-2 py-1 text-white text-xs flex-1 disabled:opacity-50"
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

              <div>
                <h3 className="text-[10px] font-bold text-gray-400 uppercase tracking-wider mb-2">Parameters</h3>
                <div className="grid grid-cols-2 gap-2">
                    <button 
                      onClick={handleRerollBorders}
                      disabled={loading}
                      className="flex items-center justify-center gap-1 text-[10px] bg-blue-900/40 text-blue-300 px-2 py-2 rounded border border-blue-900/50 hover:bg-blue-800/40"
                    >
                      <Shuffle size={10} /> Reroll Borders
                    </button>
                    <button 
                      onClick={handleRerollProvinces}
                      disabled={loading}
                      className="flex items-center justify-center gap-1 text-[10px] bg-teal-900/40 text-teal-300 px-2 py-2 rounded border border-teal-900/50 hover:bg-teal-800/40"
                    >
                      <Layers size={10} /> Reroll Provs
                    </button>
                </div>
              </div>

              {/* ... (rest of civ sliders) ... */}
              <div className="space-y-1">
                <div className="flex justify-between text-xs text-gray-400">
                  <label>Factions</label>
                  <span>{params.numFactions}</span>
                </div>
                <input
                  type="range" min="2" max="20"
                  value={params.numFactions}
                  onChange={(e) => { handleChange('numFactions', parseInt(e.target.value) as 1 | 2 | 3); }}
                  className="w-full h-1 bg-gray-700 rounded-lg appearance-none cursor-pointer accent-purple-500"
                />
              </div>
              <div className="space-y-1">
                <div className="flex justify-between text-xs text-gray-400">
                  <label>Capital Spacing</label>
                  <span>{(params.capitalSpacing * 100).toFixed(0)}%</span>
                </div>
                <input
                  type="range" min="0" max="1" step="0.1"
                  value={params.capitalSpacing}
                  onChange={(e) => { handleChange('capitalSpacing', parseFloat(e.target.value)); }}
                  className="w-full h-1 bg-gray-700 rounded-lg appearance-none cursor-pointer accent-purple-400"
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
                  className="w-full h-1 bg-gray-700 rounded-lg appearance-none cursor-pointer accent-teal-400"
                />
              </div>
              <div className="space-y-1">
                <div className="flex justify-between text-xs text-gray-400">
                  <label>Country Size Variance</label>
                  <span>{(params.civSizeVariance * 100).toFixed(0)}%</span>
                </div>
                <input
                  type="range" min="0" max="1" step="0.1"
                  value={params.civSizeVariance}
                  onChange={(e) => { handleChange('civSizeVariance', parseFloat(e.target.value)); }}
                  className="w-full h-1 bg-gray-700 rounded-lg appearance-none cursor-pointer accent-orange-400"
                />
              </div>
              <div className="space-y-1">
                <div className="flex justify-between text-xs text-gray-400">
                  <label>Seafaring (Water Crossing Cost)</label>
                  <span>{(1.0 - params.waterCrossingCost).toFixed(1)}</span>
                </div>
                <input
                  type="range" min="0.1" max="1.0" step="0.1"
                  value={params.waterCrossingCost}
                  onChange={(e) => { handleChange('waterCrossingCost', parseFloat(e.target.value)); }}
                  className="w-full h-1 bg-gray-700 rounded-lg appearance-none cursor-pointer accent-blue-400"
                />
              </div>
              <div className="space-y-1">
                <div className="flex justify-between text-xs text-gray-400">
                  <label>Territorial Waters (Range)</label>
                  <span>{params.territorialWaters?.toFixed(2)}</span>
                </div>
                <input
                  type="range" min="0.01" max="1.0" step="0.01"
                  value={params.territorialWaters ?? 0.2}
                  onChange={(e) => { handleChange('territorialWaters', parseFloat(e.target.value)); }}
                  className="w-full h-1 bg-gray-700 rounded-lg appearance-none cursor-pointer accent-cyan-600"
                />
              </div>

              {/* Lore Level */}
              <div className="space-y-1 border-t border-gray-800 pt-3">
                  <label className="text-xs text-gray-400 block mb-1">Lore Generation Detail</label>
                  <select 
                     value={params.loreLevel || 1}
                     onChange={(e) => { handleChange('loreLevel', parseInt(e.target.value) as 1 | 2 | 3); }}
                     className="w-full bg-gray-800 text-white text-xs border border-gray-700 rounded p-2"
                  >
                     <option value={1}>Level 1: Factions & Capitals</option>
                     <option value={2}>Level 2: Provinces & Towns</option>
                     <option value={3}>Level 3: Backstories (Slow)</option>
                  </select>
              </div>

              <div className="bg-gray-800/50 rounded-xl p-3 border border-gray-700 mt-4">
                <div className="flex justify-between items-center mb-2">
                  <h2 className="text-xs font-semibold text-gray-300">AI Lore</h2>
                   <button 
                    onClick={onGenerateLore}
                    disabled={generatingLore || !apiKey}
                    className={`text-[10px] px-2 py-1 rounded transition-colors ${
                      apiKey 
                        ? 'bg-blue-900/50 text-blue-300 hover:bg-blue-900 border border-blue-800' 
                        : 'bg-gray-800 text-gray-500 cursor-not-allowed border border-gray-700'
                    }`}
                  >
                    {generatingLore ? 'Thinking...' : 'Generate'}
                  </button>
                </div>
                {!apiKey && (
                  <p className="text-[9px] text-yellow-500/80 bg-yellow-500/10 p-1.5 rounded border border-yellow-500/20 mb-2">
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
                                <div key={f.id} className="text-[10px] bg-gray-900 p-1 rounded border border-gray-700">
                                    <div style={{color: f.color}} className="font-bold">{f.name}</div>
                                    <div className="text-gray-500 pl-1">Cap: {f.provinces[0]?.towns[0]?.name || 'Unknown'}</div>
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
                            className="w-full bg-gray-800 text-white text-xs border border-gray-700 rounded p-2"
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
                            className="w-full bg-gray-800 text-white text-xs border border-gray-700 rounded p-2"
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
                        <div className="border border-gray-800 rounded-lg p-3 space-y-3 bg-gray-900/40">
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
                                        className={`text-[10px] py-2 rounded border ${dymaxionSettings.mode === 'planet' ? 'bg-blue-700/70 border-blue-500 text-white' : 'bg-gray-800 border-gray-700 text-gray-300'}`}
                                    >
                                        Rotate Planet
                                    </button>
                                    <button
                                        onClick={() => { updateDymaxion({ mode: 'overlay' as DymaxionControlMode }); }}
                                        className={`text-[10px] py-2 rounded border ${dymaxionSettings.mode === 'overlay' ? 'bg-blue-700/70 border-blue-500 text-white' : 'bg-gray-800 border-gray-700 text-gray-300'}`}
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
                                    onChange={(e) => updateDymaxion({ lon: parseInt(e.target.value) })}
                                    className="w-full h-1 bg-gray-700 rounded-lg appearance-none cursor-pointer accent-blue-400"
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
                                    onChange={(e) => updateDymaxion({ lat: parseInt(e.target.value) })}
                                    className="w-full h-1 bg-gray-700 rounded-lg appearance-none cursor-pointer accent-blue-400"
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
                                    onChange={(e) => updateDymaxion({ roll: parseInt(e.target.value) })}
                                    className="w-full h-1 bg-gray-700 rounded-lg appearance-none cursor-pointer accent-blue-400"
                                />
                            </div>

                            <button
                                onClick={() => { updateDymaxion({ lon: 0, lat: 0, roll: 0 }); }}
                                className="w-full text-[10px] bg-gray-800 hover:bg-gray-700 text-gray-200 py-2 rounded border border-gray-700"
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
                        className="w-full flex items-center justify-center gap-2 bg-green-700 hover:bg-green-600 text-white py-2 rounded text-xs mt-2 disabled:opacity-50"
                    >
                        <Image size={14}/> Download PNG
                    </button>
                </div>

                <div className="border-t border-gray-800 pt-4 space-y-3">
                    <h3 className="text-xs font-bold text-gray-400 uppercase tracking-wider">File Management</h3>
                    
                    <button
                        onClick={() => { if (params) { void saveMapConfig(params, worldData || undefined); } }}
                        className="w-full flex items-center justify-center gap-2 bg-gray-800 hover:bg-gray-700 text-white py-2 rounded text-xs"
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
                        <button className="w-full flex items-center justify-center gap-2 bg-gray-800 hover:bg-gray-700 text-white py-2 rounded text-xs pointer-events-none">
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
                            className="flex-1 bg-gray-900 border border-gray-700 rounded px-2 text-xs text-white"
                        />
                        <button 
                            onClick={handleSaveBrowser}
                            disabled={!saveName}
                            className="bg-blue-600 hover:bg-blue-500 text-white px-3 rounded"
                        >
                            <Save size={14}/>
                        </button>
                    </div>

                    <div className="space-y-1 max-h-40 overflow-y-auto">
                        {savedMaps.length === 0 && <p className="text-xs text-gray-600 italic">No saved maps.</p>}
                        {savedMaps.map(entry => (
                            <div key={entry.name} className="flex items-center justify-between bg-gray-900 p-2 rounded border border-gray-800 group">
                                <div className="flex flex-col">
                                    <span className="text-xs font-bold text-gray-300">{entry.name}</span>
                                    <span className="text-[10px] text-gray-500">{new Date(entry.date).toLocaleDateString()}</span>
                                </div>
                                <div className="flex gap-1 opacity-0 group-hover:opacity-100 transition-opacity">
                                    <button onClick={() => { handleLoadBrowser(entry.params, entry.civData); }} className="text-blue-400 hover:text-white p-1"><FolderOpen size={12}/></button>
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

         {!loading ? (
             <button
              onClick={handleGenerateClick}
              className={`w-full py-3 rounded-lg font-semibold flex items-center justify-center gap-2 transition-all relative overflow-hidden bg-blue-600 hover:bg-blue-500 text-white shadow-lg shadow-blue-900/30`}
            >
              <div className="relative flex items-center gap-2 z-10">
                  <RefreshCw size={16} />
                  Generate World
              </div>
            </button>
         ) : (
            <button
              onClick={onCancel}
              className={`w-full py-3 rounded-lg font-semibold flex items-center justify-center gap-2 transition-all relative overflow-hidden bg-red-600 hover:bg-red-500 text-white shadow-lg shadow-red-900/30`}
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
