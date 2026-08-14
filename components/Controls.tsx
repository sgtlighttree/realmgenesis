import React, { useState, useEffect, useRef } from 'react';
import { WorldParams, ViewMode, LoreData, LandStyle, CivData, DisplayMode, DymaxionSettings, DymaxionControlMode, LabelVisibility, MarkerData } from '../types';
import { RefreshCw, Globe, Mountain, Lock, Unlock, Shuffle, Layers, Zap, Save, Trash2, Image, Terminal, XCircle, ChevronDown, ChevronUp, FolderOpen, Box, Copy, Check, FileCode, FileJson } from 'lucide-react';
import { exportMap, saveMapConfig, loadMapConfig, saveMapToBrowser, getSavedMaps, deleteSavedMap, ExportResolution, ProjectionType } from '../utils/export';
import { downloadSVG, downloadGeoJSON } from '../utils/exportVector';
import { NAME_STYLES } from '../utils/namegen';
import { STAR_CLASSES, STAR_CLASS_LABELS } from '../utils/planetary';
import { exportGLB } from '../utils/exportGLB';
import { WorldData } from '../types';
import DymaxionPreview2D from './DymaxionPreview2D';
import Select, { SelectOption } from './Select';
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
  /**
   * Root classes. Defaults to the classic left-rail chrome; the shell passes
   * layout-neutral classes because its Panel already draws the box.
   */
  className?: string;
  /** The shell draws its own brand header; classic keeps this one. */
  showHeader?: boolean;
  /**
   * Render mode, layer toggles, and the view-layer grid. The shell renders
   * these in its View strip, so it turns them off here rather than showing
   * every one of them twice. Classic has no View bucket and keeps them.
   * Map Overlays is NOT included: it has no equivalent in the strip.
   */
  showViewControls?: boolean;
}

type Tab = 'geo' | 'climate' | 'political' | 'system' | 'export';

const ConsoleOutput: React.FC<{ logs: string[]; isOpen: boolean }> = ({ logs, isOpen }) => {
    const boxRef = useRef<HTMLDivElement>(null);
    // Scroll THIS box, not the end marker. `scrollIntoView` walks up and scrolls
    // every scrollable ancestor, so inside the mobile sheet it dragged the whole
    // Make panel down to the console — the panel would open on the log output
    // instead of Render Mode. Setting scrollTop keeps the effect local.
    useEffect(() => {
        if (!isOpen) return;
        const el = boxRef.current;
        if (el) el.scrollTop = el.scrollHeight;
    }, [logs, isOpen]);

    if (!isOpen) return null;

    return (
        <div ref={boxRef} className="bg-black border border-edge-subtle p-2 h-32 overflow-y-auto font-mono text-[10px] space-y-1 shadow-inner relative transition-all">
            {logs.length === 0 && <div className="text-ink-faint italic text-center mt-10">System Ready</div>}
            {logs.map((log, i) => (
                <div key={i} className="text-green-400 break-words border-b border-surface/50 pb-0.5 last:border-0">
                    <span className="text-ink-faint mr-2">[{i+1}]</span>
                    {log}
                </div>
            ))}
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
  showHeader = true,
  showViewControls = true,
  className = 'w-full md:w-80 bg-surface-sunken border-r border-edge-subtle flex flex-col h-full overflow-hidden text-sm relative z-chrome',
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
      params.tectonicStrength,
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
      params.starClass,
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
              updates.tectonicStrength = 0.8;
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
              updates.tectonicStrength = 0.6;
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
              updates.tectonicStrength = 0.25;
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
              updates.tectonicStrength = 0.4;
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
    <div className={className}>
      {showHeader && (
        <div className="p-4 border-b border-edge-subtle">
          <h1 className="text-xl font-bold text-ink-strong flex items-center gap-2">
            <Globe className="text-brand" />
            RealmGenesis 3D
          </h1>
        </div>
      )}

      <div className="flex border-b border-edge-subtle">
         <button onClick={() => { setActiveTab('system'); }} className={`flex-1 py-3 text-[10px] font-semibold uppercase tracking-wide ${activeTab === 'system' ? 'text-brand-soft border-b-2 border-brand' : 'text-ink-muted hover:text-ink-soft'}`}>Sys</button>
         <button onClick={() => { setActiveTab('geo'); }} className={`flex-1 py-3 text-[10px] font-semibold uppercase tracking-wide ${activeTab === 'geo' ? 'text-brand-soft border-b-2 border-brand' : 'text-ink-muted hover:text-ink-soft'}`}>Geo</button>
         <button onClick={() => { setActiveTab('climate'); }} className={`flex-1 py-3 text-[10px] font-semibold uppercase tracking-wide ${activeTab === 'climate' ? 'text-brand-soft border-b-2 border-brand' : 'text-ink-muted hover:text-ink-soft'}`}>Clim</button>
         <button onClick={() => { setActiveTab('political'); }} className={`flex-1 py-3 text-[10px] font-semibold uppercase tracking-wide ${activeTab === 'political' ? 'text-brand-soft border-b-2 border-brand' : 'text-ink-muted hover:text-ink-soft'}`}>Civ</button>
         <button onClick={() => { setActiveTab('export'); }} className={`flex-1 py-3 text-[10px] font-semibold uppercase tracking-wide ${activeTab === 'export' ? 'text-brand-soft border-b-2 border-brand' : 'text-ink-muted hover:text-ink-soft'}`}>Exp</button>
      </div>

      <div className="flex-1 overflow-y-auto p-3 space-y-4">
        
        {activeTab === 'system' && (
          <div className="space-y-4">
             {showViewControls && (
               <div className="space-y-1">
                 <label className="text-xs text-ink-muted block">Render Mode</label>
                 <div className="flex gap-2">
                   {DISPLAY_MODES.map(m => (
                     <DisplayButton key={m.mode} mode={m.mode} label={m.label}
                       displayMode={displayMode} setDisplayMode={setDisplayMode} />
                   ))}
                 </div>
               </div>
             )}

             {/* Map Name Input */}
             <div className="space-y-1">
                 <label className="text-xs text-ink-muted block">Map Name</label>
                 <input 
                    type="text" 
                    value={params.mapName} 
                    onChange={(e) => { handleChange('mapName', e.target.value); }}
                    className="w-full bg-surface border border-edge px-2 py-1.5 text-ink-strong text-xs"
                    placeholder="My World"
                 />
             </div>

             {/* Seed Input */}
             <div className="bg-surface p-3 border border-edge-subtle">
                <label className="text-xs text-ink-muted mb-1 block">Seed</label>
                <div className="flex gap-2">
                   <input 
                      type="text" 
                      value={params.seed} 
                      onChange={(e) => { handleChange('seed', e.target.value); }}
                      disabled={seedLocked}
                      className="bg-black border border-edge px-2 py-1 text-ink-strong text-xs flex-1 min-w-0 disabled:opacity-50"
                   />
                   <button
                      onClick={() => {
                        void navigator.clipboard.writeText(params.seed);
                        setSeedCopied(true);
                        setTimeout(() => setSeedCopied(false), 1500);
                      }}
                      title="Copy seed to clipboard"
                      aria-label={seedCopied ? 'Seed copied to clipboard' : 'Copy seed to clipboard'}
                      className="text-ink-muted hover:text-ink-strong transition-colors"
                   >
                      {seedCopied ? <Check size={14} className="text-green-400" /> : <Copy size={14} />}
                   </button>
                   <button
                      onClick={() => { setSeedLocked(!seedLocked); }}
                      aria-label={seedLocked ? 'Unlock seed' : 'Lock seed'}
                      aria-pressed={seedLocked}
                      className={`${seedLocked ? 'text-brand' : 'text-ink-muted'} hover:text-ink-strong transition-colors`}
                   >
                      {seedLocked ? <Lock size={14}/> : <Unlock size={14}/>}
                   </button>
                   <button onClick={handleRandomizeSeed} disabled={seedLocked} aria-label="Randomize seed" className="text-ink-muted hover:text-ink-strong disabled:opacity-50">
                      <Shuffle size={14} />
                   </button>
                </div>
             </div>
             
             <div className="space-y-1">
              <div className="flex justify-between items-center text-xs text-ink-muted">
                <label>Resolution</label>
                <input
                    type="number"
                    min="2000"
                    max="200000"
                    step="1000"
                    value={params.points}
                    onChange={(e) => { handleNumberChange('points', e.target.value, 2000, 200000); }}
                    className="w-24 bg-surface border border-edge px-1 py-0.5 text-right text-ink-strong text-xs"
                />
              </div>
              <input
                type="range"
                min="2000"
                max="200000"
                step="1000"
                value={Math.min(200000, params.points)}
                onChange={(e) => { handleChange('points', parseInt(e.target.value) as 1 | 2 | 3); }}
className="w-full h-1 bg-surface-hover appearance-none cursor-pointer accent-brand-soft"
                />
              </div>

              {showViewControls && layerToggles.map((t, i) => (
               <LayerToggleRow key={t.key} toggle={t}
                 className={i === 0 ? 'pt-2 border-t border-edge-subtle' : 'pt-2'} />
             ))}

             <OverlayToggles
               labelVisibility={labelVisibility}
               setLabelVisibility={setLabelVisibility}
             />

            <div className="flex items-center justify-between text-xs text-ink-muted pt-2">
                 <div className="flex items-center gap-2">
                    <Zap size={12} className={autoUpdate ? "text-yellow-400" : "text-ink-faint"}/>
                    <label>Auto-Update (Low Res)</label>
                 </div>
                 <input 
                    type="checkbox"
                    checked={autoUpdate}
                    onChange={(e) => { setAutoUpdate(e.target.checked); }}
                    disabled={params.points > 20000}
                    className="bg-surface-hover"
                 />
            </div>
            
            {showViewControls && (
              <div className="pt-3 border-t border-edge-subtle">
                <h3 className="text-xs font-semibold text-ink-muted mb-2">View Layer</h3>
                <ViewLayerGrid viewMode={viewMode} setViewMode={setViewMode} />
              </div>
            )}

            <div className="pt-4 border-t border-edge-subtle space-y-3">
              <h3 className="text-xs font-semibold text-ink-muted mb-2">AI Settings (BYOK)</h3>
              <div className="bg-surface p-3 border border-edge-subtle space-y-2">
                <div className="flex items-center justify-between">
                  <label className="text-xs text-ink-muted">Gemini API Key</label>
                  <a 
                    href="https://aistudio.google.com/app/apikey" 
                    target="_blank" 
                    rel="noopener noreferrer"
                    className="text-[10px] text-brand-soft hover:underline flex items-center gap-1"
                  >
                    Get Key <Layers size={8} />
                  </a>
                </div>
                <input 
                  type="password"
                  value={apiKey}
                  onChange={(e) => { onApiKeyChange(e.target.value); }}
                  placeholder="Paste your API key here..."
                  className="w-full bg-black border border-edge px-2 py-1.5 text-ink-strong text-xs"
                />
                <p className="text-[10px] text-ink-muted italic">
                  Key is stored ephemerally in memory and will be lost on refresh.
                </p>
              </div>
            </div>

            {/* World Stats */}
            {worldData && (
              <div className="pt-4 border-t border-edge-subtle">
                <button
                  className="flex items-center justify-between w-full text-xs text-ink-muted hover:text-ink"
                  onClick={() => setShowStats(v => !v)}
                >
                  <span className="font-semibold text-ink-muted">World Stats</span>
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
                    <div className="mt-2 space-y-1 text-[10px] text-ink-muted">
                      <div className="flex justify-between"><span>Land coverage</span><span className="text-ink">{landPct}%</span></div>
                      <div className="flex justify-between"><span>Total cells</span><span className="text-ink">{cells.length.toLocaleString()}</span></div>
                      {totalPop > 0 && <div className="flex justify-between"><span>Population</span><span className="text-ink">{popStr}</span></div>}
                      {worldData.civData && <div className="flex justify-between"><span>Factions</span><span className="text-ink">{worldData.civData.factions.length}</span></div>}
                      {worldData.rivers && <div className="flex justify-between"><span>Rivers</span><span className="text-ink">{worldData.rivers.length}</span></div>}
                      <div className="mt-1 pt-1 border-t border-edge-subtle">
                        {topBiomes.map(([biome, count]) => (
                          <div key={biome} className="flex justify-between">
                            <span className="truncate">{biome}</span>
                            <span className="text-ink ml-2">{((count/cells.length)*100).toFixed(1)}%</span>
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
                 <label className="text-xs text-ink-muted block mb-1">Terrain Preset</label>
                 <Select
                    value={params.landStyle}
                    options={[{ value: 'Continents', label: 'Continents' }, { value: 'Pangea', label: 'Pangea' }, { value: 'Archipelago', label: 'Archipelago' }, { value: 'Islands', label: 'Islands' }, { value: 'Custom', label: 'Custom' }] as SelectOption<LandStyle>[]}
                    onChange={(v) => { handlePresetChange(v); }}
                    label="Terrain preset"
                    className="w-full"
                    triggerClassName="w-full justify-between bg-surface-raised border-edge px-2 py-2 text-xs"
                 />
              </div>
              <div className="space-y-1">
                <div className="flex justify-between text-xs text-ink-muted">
                  <label>Sea Level</label>
                  <span>{(params.seaLevel * 100).toFixed(0)}%</span>
                </div>
                <input
                  type="range" min="0.1" max="0.9" step="0.05"
                  value={params.seaLevel}
                  onChange={(e) => { handleChange('seaLevel', parseFloat(e.target.value)); }}
                  className="w-full h-1 bg-surface-hover appearance-none cursor-pointer accent-brand-soft"
                />
              </div>
               <div className="space-y-1">
                  <div className="flex justify-between text-xs text-ink-muted">
                    <label>Planet Radius</label>
                    <span>{params.planetRadius} km</span>
                  </div>
                  <input
                    type="range" min="1000" max="20000" step="100"
                    value={params.planetRadius || 6371}
                    onChange={(e) => { handleChange('planetRadius', parseInt(e.target.value) as 1 | 2 | 3); }}
                    className="w-full h-1 bg-surface-hover appearance-none cursor-pointer accent-brand-soft"
                  />
              </div>
              <div className="space-y-1">
                <div className="flex justify-between text-xs text-ink-muted">
                  <label>Tectonic Plates</label>
                  <span>{params.plates}</span>
                </div>
                <input
                  type="range" min="2" max="50" step="1"
                  value={params.plates}
                  onChange={(e) => { handleAdvancedChange('plates', parseInt(e.target.value) as 1 | 2 | 3); }}
                  className="w-full h-1 bg-surface-hover appearance-none cursor-pointer accent-brand-soft"
                />
              </div>
              <div className="space-y-1">
                <div className="flex justify-between text-xs text-ink-muted">
                  <label>Terrain Roughness</label>
                  <span>{(params.roughness * 100).toFixed(0)}%</span>
                </div>
                <input
                  type="range" min="0" max="1" step="0.1"
                  value={params.roughness}
                  onChange={(e) => { handleAdvancedChange('roughness', parseFloat(e.target.value)); }}
                  className="w-full h-1 bg-surface-hover appearance-none cursor-pointer accent-brand-soft"
                />
              </div>
              <div className="space-y-1" title="FBM octave count for structural terrain noise. More octaves = finer nested detail; fewer = smoother, broader forms.">
                <div className="flex justify-between text-xs text-ink-muted">
                  <label>Detail Octaves</label>
                  <span>{params.detailLevel}</span>
                </div>
                <input
                  type="range" min="1" max="6" step="1"
                  value={params.detailLevel}
                  onChange={(e) => { handleAdvancedChange('detailLevel', parseInt(e.target.value, 10)); }}
                  className="w-full h-1 bg-surface-hover appearance-none cursor-pointer accent-brand-soft"
                />
              </div>
              <div className="space-y-1" title="Randomizes the cell grid. 0 = regular Fibonacci lattice; 1 = fully jittered organic cells.">
                <div className="flex justify-between text-xs text-ink-muted">
                  <label>Cell Jitter</label>
                  <span>{(params.cellJitter * 100).toFixed(0)}%</span>
                </div>
                <input
                  type="range" min="0" max="1" step="0.05"
                  value={params.cellJitter}
                  onChange={(e) => { handleAdvancedChange('cellJitter', parseFloat(e.target.value)); }}
                  className="w-full h-1 bg-surface-hover appearance-none cursor-pointer accent-brand-soft"
                />
              </div>
              <div className="space-y-1" title="Controls terrain feature size. Lower = broader continents; higher = more fragmented detail.">
                <div className="flex justify-between text-xs text-ink-muted">
                  <label>Feature Frequency</label>
                  <span>{params.noiseScale.toFixed(1)}</span>
                </div>
                <input
                  type="range" min="0.1" max="5.0" step="0.1"
                  value={params.noiseScale}
                  onChange={(e) => { handleAdvancedChange('noiseScale', parseFloat(e.target.value)); }}
                  className="w-full h-1 bg-surface-hover appearance-none cursor-pointer accent-brand-soft"
                />
              </div>
              <div className="space-y-1" title="0 = smooth hills (FBM). 1 = sharp jagged mountain ridges (ridged noise).">
                <div className="flex justify-between text-xs text-ink-muted">
                  <label>Ridge Intensity</label>
                  <span>{(params.ridgeBlend * 100).toFixed(0)}%</span>
                </div>
                <input
                  type="range" min="0" max="1" step="0.1"
                  value={params.ridgeBlend}
                  onChange={(e) => { handleAdvancedChange('ridgeBlend', parseFloat(e.target.value)); }}
                  className="w-full h-1 bg-surface-hover appearance-none cursor-pointer accent-brand-soft"
                />
              </div>
              <div className="space-y-1" title="Amplifies terrain above sea level using a power curve. >1.0 = taller peaks; <1.0 = flatter land.">
                <div className="flex justify-between text-xs text-ink-muted">
                  <label>Mountain Heights</label>
                  <span>{params.mountainHeight.toFixed(2)}</span>
                </div>
                <input
                  type="range" min="0.5" max="2.0" step="0.05"
                  value={params.mountainHeight}
                  onChange={(e) => { handleAdvancedChange('mountainHeight', parseFloat(e.target.value)); }}
                  className="w-full h-1 bg-surface-hover appearance-none cursor-pointer accent-brand-soft"
                />
              </div>
              <div className="space-y-1" title="Amplifies ocean depth below sea level using a power curve. >1.0 = deeper trenches; <1.0 = shallower ocean floor.">
                <div className="flex justify-between text-xs text-ink-muted">
                  <label>Sea / Trench Depth</label>
                  <span>{params.oceanDepth.toFixed(2)}</span>
                </div>
                <input
                  type="range" min="0.5" max="2.0" step="0.05"
                  value={params.oceanDepth}
                  onChange={(e) => { handleAdvancedChange('oceanDepth', parseFloat(e.target.value)); }}
                  className="w-full h-1 bg-surface-hover appearance-none cursor-pointer accent-brand-soft"
                />
              </div>
              <div className="space-y-1" title="Domain warping — twists terrain shapes for more organic, swirled coastlines and mountain ranges.">
                <div className="flex justify-between text-xs text-ink-muted">
                  <label>Swirl / Warp</label>
                  <span>{params.warpStrength.toFixed(1)}</span>
                </div>
                <input
                  type="range" min="0" max="2.0" step="0.1"
                  value={params.warpStrength}
                  onChange={(e) => { handleAdvancedChange('warpStrength', parseFloat(e.target.value)); }}
                  className="w-full h-1 bg-surface-hover appearance-none cursor-pointer accent-brand-soft"
                />
              </div>
               <div className="space-y-1" title="How strongly tectonic plate boundaries shape mountain ranges and rifts. Capped at 1.0 internally.">
                <div className="flex justify-between text-xs text-ink-muted">
                  <label>Tectonic Strength</label>
                  <span>{params.tectonicStrength.toFixed(2)}</span>
                </div>
                <input
                  type="range" min="0.1" max="1.0" step="0.05"
                  value={params.tectonicStrength}
                  onChange={(e) => { handleAdvancedChange('tectonicStrength', parseFloat(e.target.value)); }}
                  className="w-full h-1 bg-surface-hover appearance-none cursor-pointer accent-brand-soft"
                />
              </div>
              <div className="space-y-1">
                <div className="flex justify-between text-xs text-ink-muted">
                  <label>Hydraulic Erosion</label>
                  <span>{params.erosionIterations} Steps</span>
                </div>
                <input
                  type="range" min="0" max="50" step="1"
                  value={params.erosionIterations}
                  onChange={(e) => { handleAdvancedChange('erosionIterations', parseInt(e.target.value) as 1 | 2 | 3); }}
                  className="w-full h-1 bg-surface-hover appearance-none cursor-pointer accent-brand-soft"
                />
              </div>
              <div className="space-y-1" title="How strongly mountain belts tend to align with continental margins (V3).">
                <div className="flex justify-between text-xs text-ink-muted">
                  <label>Margin Coupling</label>
                  <span>{params.marginCoupling?.toFixed(2) ?? '0.30'}</span>
                </div>
                <input
                  type="range" min="0" max="1" step="0.05"
                  value={params.marginCoupling ?? 0.3}
                  onChange={(e) => { handleAdvancedChange('marginCoupling', parseFloat(e.target.value)); }}
                  className="w-full h-1 bg-surface-hover appearance-none cursor-pointer accent-brand-soft"
                />
              </div>
              <div className="space-y-1" title="Number of simulation timesteps for the tectonic model (V3).">
                <div className="flex justify-between text-xs text-ink-muted">
                  <label>Timesteps</label>
                  <span>{params.numTimesteps ?? 20}</span>
                </div>
                <input
                  type="range" min="5" max="60" step="1"
                  value={params.numTimesteps ?? 20}
                  onChange={(e) => { handleAdvancedChange('numTimesteps', parseInt(e.target.value)); }}
                  className="w-full h-1 bg-surface-hover appearance-none cursor-pointer accent-brand-soft"
                />
              </div>
              <div className="space-y-1" title="Macro-cell resolution for the tectonic simulation (V3).">
                <div className="flex justify-between text-xs text-ink-muted">
                  <label>Macro-Cells</label>
                  <span>{params.simulationResolution ?? 10000}</span>
                </div>
                <input
                  type="range" min="5000" max="20000" step="1000"
                  value={params.simulationResolution ?? 10000}
                  onChange={(e) => { handleAdvancedChange('simulationResolution', parseInt(e.target.value)); }}
                  className="w-full h-1 bg-surface-hover appearance-none cursor-pointer accent-brand-soft"
                />
              </div>
              <div className="space-y-1" title="How irregularly plate seeds are distributed. 0 = uniform Fibonacci, 1 = chaotic, >1 = strongly varied plate sizes.">
                <div className="flex justify-between text-xs text-ink-muted">
                  <label>Plate Jitter</label>
                  <span>{params.plateJitter?.toFixed(2) ?? '0.30'}</span>
                </div>
                <input
                  type="range" min="0" max="3" step="0.05"
                  value={params.plateJitter ?? 0.3}
                  onChange={(e) => { handleAdvancedChange('plateJitter', parseFloat(e.target.value)); }}
                  className="w-full h-1 bg-surface-hover appearance-none cursor-pointer accent-brand-soft"
                />
              </div>
              <div className="space-y-1" title="How jagged/fractal plate boundaries are. 0 = clean great-circle arcs, 1 = highly irregular, >1 = extreme interlocking fracture.">
                <div className="flex justify-between text-xs text-ink-muted">
                  <label>Boundary Roughness</label>
                  <span>{params.boundaryRoughness?.toFixed(2) ?? '0.30'}</span>
                </div>
                <input
                  type="range" min="0" max="3" step="0.05"
                  value={params.boundaryRoughness ?? 0.3}
                  onChange={(e) => { handleAdvancedChange('boundaryRoughness', parseFloat(e.target.value)); }}
                  className="w-full h-1 bg-surface-hover appearance-none cursor-pointer accent-brand-soft"
                />
              </div>
              <div className="space-y-1" title="Seafloor spreading rate. Lower = older, deeper ocean floor away from mid-ocean ridges (GDH1 bathymetry).">
                <div className="flex justify-between text-xs text-ink-muted">
                  <label>Spreading Rate</label>
                  <span>{params.spreadRate?.toFixed(3) ?? '0.008'}</span>
                </div>
                <input
                  type="range" min="0.004" max="0.02" step="0.001"
                  value={params.spreadRate ?? 0.008}
                  onChange={(e) => { handleAdvancedChange('spreadRate', parseFloat(e.target.value)); }}
                  className="w-full h-1 bg-surface-hover appearance-none cursor-pointer accent-brand-soft"
                />
              </div>
              <div className="space-y-1" title="Overall ocean-floor depth. Lower = shallower seas (floor rises toward the coast); higher = deeper abyss. Coastline stays put — this scales mean water depth, unlike Ocean Depth which reshapes the trench-vs-shelf contrast.">
                <div className="flex justify-between text-xs text-ink-muted">
                  <label>Seafloor Depth</label>
                  <span>{params.seafloorDepth?.toFixed(2) ?? '1.00'}</span>
                </div>
                <input
                  type="range" min="0.3" max="2" step="0.05"
                  value={params.seafloorDepth ?? 1.0}
                  onChange={(e) => { handleAdvancedChange('seafloorDepth', parseFloat(e.target.value)); }}
                  className="w-full h-1 bg-surface-hover appearance-none cursor-pointer accent-brand-soft"
                />
              </div>
              <div className="space-y-1" title="How many shear-driven microplates to inject along high-strain boundaries. 0 = none, higher = more elongated slivers breaking up round plates.">
                <div className="flex justify-between text-xs text-ink-muted">
                  <label>Microplates</label>
                  <span>{params.microplateIntensity?.toFixed(2) ?? '0.35'}</span>
                </div>
                <input
                  type="range" min="0" max="1" step="0.05"
                  value={params.microplateIntensity ?? 0.35}
                  onChange={(e) => { handleAdvancedChange('microplateIntensity', parseFloat(e.target.value)); }}
                  className="w-full h-1 bg-surface-hover appearance-none cursor-pointer accent-brand-soft"
                />
              </div>
              <div className="space-y-1" title="Seed-chain length for plate growth. 0 = round Voronoi blobs, higher = elongated, band-shaped plates along their motion direction.">
                <div className="flex justify-between text-xs text-ink-muted">
                  <label>Plate Elongation</label>
                  <span>{params.plateElongation?.toFixed(2) ?? '0.40'}</span>
                </div>
                <input
                  type="range" min="0" max="1" step="0.05"
                  value={params.plateElongation ?? 0.4}
                  onChange={(e) => { handleAdvancedChange('plateElongation', parseFloat(e.target.value)); }}
                  className="w-full h-1 bg-surface-hover appearance-none cursor-pointer accent-brand-soft"
                />
              </div>
           </div>
        )}

        {/* Climate Tab content */}
        {activeTab === 'climate' && (
           <div className="space-y-5">
              <div className="space-y-1">
                <div className="flex justify-between text-xs text-ink-muted">
                  <label>Axial Tilt (Visual & Climatic)</label>
                  <span>{params.axialTilt || 0}°</span>
                </div>
                <input
                  type="range" min="-90" max="90" step="1"
                  value={params.axialTilt || 0}
                  onChange={(e) => { handleChange('axialTilt', parseFloat(e.target.value)); }}
                  className="w-full h-1 bg-surface-hover appearance-none cursor-pointer accent-brand-soft"
                />
             </div>
              {/* D1: Season is render-only — it recolors (temperature, snow line,
                  biome edges) without regenerating. Neutral (0.5) shows the
                  canonical annual-mean world. Read out as the subsolar latitude. */}
              <div className="space-y-1">
                <div className="flex justify-between text-xs text-ink-muted">
                  <label>Season (subsolar latitude)</label>
                  <span>{(() => {
                    const decl = (params.axialTilt || 0) * Math.sin(2 * Math.PI * ((params.season ?? 0.5)));
                    return Math.abs(decl) < 0.05 ? 'Equinox (annual mean)' : `Sun ${decl >= 0 ? 'N' : 'S'} ${Math.abs(decl).toFixed(1)}°`;
                  })()}</span>
                </div>
                <input
                  type="range" min="0" max="1" step="0.01"
                  value={params.season ?? 0.5}
                  onChange={(e) => { handleChange('season', parseFloat(e.target.value)); }}
                  className="w-full h-1 bg-surface-hover appearance-none cursor-pointer accent-brand-soft"
                  disabled={!params.axialTilt}
                />
             </div>
              {/* D5: host star spectral class — scales global insolation
                  (temperature → biomes → sea-ice all follow). G is Sun-like. */}
              <div className="space-y-1">
                <label className="text-xs text-ink-muted block mb-1">Host Star Class</label>
                <Select
                   value={params.starClass || 'G'}
                   options={STAR_CLASSES.map(sc => ({ value: sc, label: STAR_CLASS_LABELS[sc] }))}
                   onChange={(v) => { handleChange('starClass', v); }}
                   label="Host star class"
                   className="w-full"
                   triggerClassName="w-full justify-between bg-surface-raised border-edge px-2 py-2 text-xs"
                />
             </div>
             {/* ... (rest of climate sliders) ... */}
              <div className="space-y-1">
                <div className="flex justify-between text-xs text-ink-muted">
                  <label>Equator Temp (°C)</label>
                  <span>{params.baseTemperature}°C</span>
                </div>
                <input
                  type="range" min="-10" max="50" step="1"
                  value={params.baseTemperature}
                  onChange={(e) => { handleChange('baseTemperature', parseFloat(e.target.value)); }}
                  className="w-full h-1 bg-surface-hover appearance-none cursor-pointer accent-brand-soft"
                />
              </div>
              <div className="space-y-1">
                <div className="flex justify-between text-xs text-ink-muted">
                  <label>Pole Temp (°C)</label>
                  <span>{params.poleTemperature}°C</span>
                </div>
                <input
                  type="range" min="-50" max="20" step="1"
                  value={params.poleTemperature}
                  onChange={(e) => { handleChange('poleTemperature', parseFloat(e.target.value)); }}
                  className="w-full h-1 bg-surface-hover appearance-none cursor-pointer accent-brand-soft"
                />
              </div>
              <div className="space-y-1">
                <div className="flex justify-between text-xs text-ink-muted">
                  <label>Global Rainfall</label>
                  <span>{params.rainfallMultiplier.toFixed(1)}x</span>
                </div>
                <input
                  type="range" min="0" max="3" step="0.1"
                  value={params.rainfallMultiplier}
                  onChange={(e) => { handleChange('rainfallMultiplier', parseFloat(e.target.value)); }}
                  className="w-full h-1 bg-surface-hover appearance-none cursor-pointer accent-brand-soft"
                />
              </div>
              <div className="space-y-1" title="How far wind carries moisture inland before it dissipates. Higher = wetter interiors; lower = stronger rain shadows.">
                <div className="flex justify-between text-xs text-ink-muted">
                  <label>Wind Strength / Moisture Transport</label>
                  <span>{(params.moistureTransport * 100).toFixed(0)}%</span>
                </div>
                <input
                  type="range" min="0" max="1" step="0.1"
                  value={params.moistureTransport}
                  onChange={(e) => { handleChange('moistureTransport', parseFloat(e.target.value)); }}
                  className="w-full h-1 bg-surface-hover appearance-none cursor-pointer accent-brand-soft"
                />
                <p className="text-[10px] text-ink-muted">Affects rain shadows & moisture spread</p>
              </div>
              <div className="space-y-1" title="Adds simplex noise to temperature — creates local hot/cold anomalies beyond the baseline latitude gradient.">
                <div className="flex justify-between text-xs text-ink-muted">
                  <label>Random Temp</label>
                  <span>{params.temperatureVariance}</span>
                </div>
                <input
                  type="range" min="0" max="20" step="1"
                  value={params.temperatureVariance}
                  onChange={(e) => { handleChange('temperatureVariance', parseFloat(e.target.value)); }}
                  className="w-full h-1 bg-surface-hover appearance-none cursor-pointer accent-brand-soft"
                />
              </div>
           </div>
        )}

        {/* Civ Tab content */}
        {activeTab === 'political' && (
           <div className="space-y-5">
              {/* Civ Seed Input */}
             <div className="bg-surface p-3 border border-edge-subtle">
                <label className="text-xs text-ink-muted mb-1 block">Civ Seed</label>
                <div className="flex gap-2">
                   <input 
                      type="text" 
                      value={params.civSeed} 
                      onChange={(e) => { handleChange('civSeed', e.target.value); }}
                      disabled={civSeedLocked}
                      className="bg-black border border-edge px-2 py-1 text-ink-strong text-xs flex-1 min-w-0 disabled:opacity-50"
                   />
                   <button 
                      onClick={() => { setCivSeedLocked(!civSeedLocked); }}
                      aria-label={civSeedLocked ? 'Unlock civilization seed' : 'Lock civilization seed'}
                      aria-pressed={civSeedLocked}
                      className={`${civSeedLocked ? 'text-brand' : 'text-ink-muted'} hover:text-ink-strong transition-colors`}
                   >
                      {civSeedLocked ? <Lock size={14}/> : <Unlock size={14}/>}
                   </button>
                   <button onClick={handleRandomizeCivSeed} disabled={civSeedLocked} aria-label="Randomize civilization seed" className="text-ink-muted hover:text-ink-strong disabled:opacity-50">
                      <Shuffle size={14} />
                   </button>
                </div>
             </div>

              {/* Generation Parameters — collapsible */}
              <div>
                <button
                  onClick={() => setShowCivParams(v => !v)}
                  className="flex items-center justify-between w-full text-[10px] font-semibold text-ink-muted uppercase tracking-wide hover:text-ink-soft transition-colors"
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
                <div className="flex justify-between text-xs text-ink-muted">
                  <label>Cultures</label>
                  <span>{params.numCultures}</span>
                </div>
                <input
                  type="range" min="2" max="8"
                  value={params.numCultures}
                  onChange={(e) => { handleChange('numCultures', parseInt(e.target.value) as 1 | 2 | 3); }}
                  className="w-full h-1 bg-surface-hover appearance-none cursor-pointer accent-brand-soft"
                />
              </div>
              <div className="space-y-1">
                <div className="flex justify-between text-xs text-ink-muted">
                  <label>Factions</label>
                  <span>{params.numFactions}</span>
                </div>
                <input
                  type="range" min="2" max="20"
                  value={params.numFactions}
                  onChange={(e) => { handleChange('numFactions', parseInt(e.target.value) as 1 | 2 | 3); }}
                  className="w-full h-1 bg-surface-hover appearance-none cursor-pointer accent-brand-soft"
                />
              </div>
              <div className="space-y-1" title="Minimum angular separation between faction capitals. Higher = capitals spawn further apart, producing more evenly distributed territories.">
                <div className="flex justify-between text-xs text-ink-muted">
                  <label>Capital Spacing</label>
                  <span>{(params.capitalSpacing * 100).toFixed(0)}%</span>
                </div>
                <input
                  type="range" min="0" max="1" step="0.1"
                  value={params.capitalSpacing}
                  onChange={(e) => { handleChange('capitalSpacing', parseFloat(e.target.value)); }}
                  className="w-full h-1 bg-surface-hover appearance-none cursor-pointer accent-brand-soft"
                />
              </div>
              <div className="space-y-1">
                <div className="flex justify-between text-xs text-ink-muted">
                  <label>Province Size (Admin Density)</label>
                  <span>{(params.provinceSize || 0.5).toFixed(2)}</span>
                </div>
                <input
                  type="range" min="0.1" max="1.0" step="0.1"
                  value={params.provinceSize || 0.5}
                  onChange={(e) => { handleChange('provinceSize', parseFloat(e.target.value)); }}
                  className="w-full h-1 bg-surface-hover appearance-none cursor-pointer accent-brand-soft"
                />
              </div>
              <div className="space-y-1" title="How unequal faction sizes can be. 0 = all factions roughly equal; 1 = some factions much larger than others.">
                <div className="flex justify-between text-xs text-ink-muted">
                  <label>Country Size Variance</label>
                  <span>{(params.civSizeVariance * 100).toFixed(0)}%</span>
                </div>
                <input
                  type="range" min="0" max="1" step="0.1"
                  value={params.civSizeVariance}
                  onChange={(e) => { handleChange('civSizeVariance', parseFloat(e.target.value)); }}
                  className="w-full h-1 bg-surface-hover appearance-none cursor-pointer accent-brand-soft"
                />
              </div>
              <div className="space-y-1" title="How easily factions cross water. Higher = more seafaring civilisations that readily claim distant islands.">
                <div className="flex justify-between text-xs text-ink-muted">
                  <label>Seafaring (Water Crossing Cost)</label>
                  <span>{(1.0 - params.waterCrossingCost).toFixed(1)}</span>
                </div>
                <input
                  type="range" min="0.1" max="1.0" step="0.1"
                  value={params.waterCrossingCost}
                  onChange={(e) => { handleChange('waterCrossingCost', parseFloat(e.target.value)); }}
                  className="w-full h-1 bg-surface-hover appearance-none cursor-pointer accent-brand-soft"
                />
              </div>
              <div className="space-y-1" title="How far from coastline a faction claims ocean cells as territorial waters.">
                <div className="flex justify-between text-xs text-ink-muted">
                  <label>Territorial Waters (Range)</label>
                  <span>{params.territorialWaters?.toFixed(2)}</span>
                </div>
                <input
                  type="range" min="0.01" max="1.0" step="0.01"
                  value={params.territorialWaters ?? 0.2}
                  onChange={(e) => { handleChange('territorialWaters', parseFloat(e.target.value)); }}
                  className="w-full h-1 bg-surface-hover appearance-none cursor-pointer accent-brand-soft"
                />
              </div>

                  </div>
                )}
              </div>

              {/* Factions editor */}
              {worldData?.civData && (
                <div className="space-y-2 pt-2 border-t border-edge-subtle">
                  <h3 className="text-[10px] font-semibold text-ink-muted uppercase tracking-wide pt-1">Factions</h3>
                  {worldData.civData.factions.map(f => {
                    const otherFactions = worldData.civData!.factions.filter(o => o.id !== f.id);
                    const mergeTarget = mergeTargets[f.id] ?? '';
                    return (
                    <div key={f.id} className="flex flex-col gap-1.5 bg-surface p-2 border border-edge-subtle">
                      <div className="flex items-center gap-2">
                        <input
                          type="color"
                          value={f.color}
                          onChange={e => onEditFaction?.(f.id, { color: e.target.value })}
                          className="w-7 h-6 border border-edge bg-transparent cursor-pointer flex-shrink-0"
                          title="Faction color"
                        />
                        <input
                          type="text"
                          value={f.name}
                          onChange={e => onEditFaction?.(f.id, { name: e.target.value })}
                          className="flex-1 bg-black border border-edge px-2 py-1 text-ink-strong text-xs focus:outline-none focus:border-brand"
                          placeholder="Faction name"
                        />
                      </div>
                      {onMergeFactions && otherFactions.length > 0 && (
                        <div className="flex items-center gap-1.5 pl-1">
                          <Select
                            value={mergeTarget === '' ? '' : String(mergeTarget)}
                            options={[
                              { value: '', label: 'Merge into\u2026' },
                              ...otherFactions.map(o => ({ value: String(o.id), label: o.name })),
                            ]}
                            onChange={(v) => {
                              const val = v === '' ? '' : parseInt(v, 10);
                              setMergeTargets(prev => ({ ...prev, [f.id]: val }));
                            }}
                            label={`Merge ${f.name} into another faction`}
                            className="flex-1 min-w-0"
                            triggerClassName="w-full justify-between bg-black border-edge px-1.5 py-1 text-[10px] text-ink-soft"
                          />
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
              <div className="space-y-1 border-t border-edge-subtle pt-3">
                  <label className="text-xs text-ink-muted block mb-1">Name Style</label>
                  <Select
                     value={params.nameStyle || 'fantasy'}
                     options={NAME_STYLES.map(style => ({ value: style, label: style.charAt(0).toUpperCase() + style.slice(1) }))}
                     onChange={(v) => { handleChange('nameStyle', v); }}
                     label="Name style"
                     className="w-full"
                     triggerClassName="w-full justify-between bg-surface-raised border-edge px-2 py-2 text-xs"
                  />
                  <p className="text-[10px] text-ink-muted">Applies on Reroll Borders or regeneration.</p>
              </div>

              {/* Lore Level */}
              <div className="space-y-1 border-t border-edge-subtle pt-3">
                  <label className="text-xs text-ink-muted block mb-1">Lore Generation Detail</label>
                  <Select
                     value={String(params.loreLevel || 1)}
                     options={[
                       { value: '1', label: 'Level 1: Factions & Capitals' },
                       { value: '2', label: 'Level 2: Provinces & Towns' },
                       { value: '3', label: 'Level 3: Backstories (Slow)' },
                     ]}
                     onChange={(v) => { handleChange('loreLevel', parseInt(v, 10) as 1 | 2 | 3); }}
                     label="Lore generation detail"
                     className="w-full"
                     triggerClassName="w-full justify-between bg-surface-raised border-edge px-2 py-2 text-xs"
                  />
              </div>

              <div className="bg-surface-raised/50 p-3 border border-edge mt-4">
                <div className="flex justify-between items-center mb-2">
                  <h2 className="text-xs font-semibold text-ink-soft">AI Lore</h2>
                   <button 
                    onClick={onGenerateLore}
                    disabled={generatingLore || !apiKey}
                    className={`text-[10px] px-2 py-1 transition-colors ${
 apiKey 
 ? 'bg-blue-900/50 text-blue-300 hover:bg-blue-900 border border-blue-800' 
 : 'bg-surface-raised text-ink-muted cursor-not-allowed border border-edge'
 }`}
                  >
                    {generatingLore ? 'Thinking...' : 'Generate'}
                  </button>
                </div>
                {!apiKey && (
                  <p className="text-[10px] text-yellow-500/80 bg-yellow-500/10 p-1.5 border border-yellow-500/20 mb-2">
                    Enter a Gemini API Key in the "Sys" tab to enable AI lore.
                  </p>
                )}
                {lore ? (
                  <div className="space-y-2">
                    <h3 className="font-bold text-ink-strong text-xs">{lore.name}</h3>
                    <p className="text-[10px] text-ink-muted max-h-32 overflow-y-auto">
                      {lore.description}
                    </p>
                    {worldData?.civData && (
                        <div className="space-y-1 mt-2">
                            {worldData.civData.factions.map(f => (
                                <div key={f.id} className="text-[10px] bg-surface p-1 border border-edge">
                                    <div style={{color: f.color}} className="font-bold">{f.name}</div>
                                    <div className="text-ink-muted pl-1">Cap: {f.provinces.flatMap(p => p.towns).find(t => t.cellId === f.capitalId)?.name || 'Unknown'}</div>
                                    {f.description && <div className="text-ink-muted italic pl-1 mt-1 border-t border-edge-subtle pt-1">{f.description}</div>}
                                </div>
                            ))}
                        </div>
                    )}
                  </div>
                ) : (
                  <p className="text-[10px] text-ink-faint italic">Generate a world first.</p>
                )}
              </div>
           </div>
        )}

        {/* Export Tab Content */}
        {activeTab === 'export' && (
            <div className="space-y-6">
                 {/* ... (export tab same as before) ... */}
                 <div className="space-y-2">
                    <h3 className="text-xs font-bold text-ink-muted uppercase tracking-wider mb-2">Image Export</h3>
                    
                    <div className="space-y-1">
                        <label className="text-xs text-ink-muted">Resolution</label>
                        <Select
                            value={String(expRes)}
                            options={[
                              { value: '2048', label: '2K (2048px)' },
                              { value: '4096', label: '4K (4096px)' },
                              { value: '8192', label: '8K (8192px)' },
                            ]}
                            onChange={(v) => { setExpRes(parseInt(v, 10) as ExportResolution); }}
                            label="Export resolution"
                            className="w-full"
                            triggerClassName="w-full justify-between bg-surface-raised border-edge px-2 py-2 text-xs"
                        />
                    </div>

                    <div className="space-y-1">
                        <label className="text-xs text-ink-muted">Projection</label>
                        <Select
                            value={expProj}
                            options={[
                              { value: 'equirectangular', label: 'Equirectangular' },
                              { value: 'mercator', label: 'Mercator' },
                              { value: 'winkeltripel', label: 'Winkel Tripel' },
                              { value: 'robinson', label: 'Robinson' },
                              { value: 'mollweide', label: 'Mollweide' },
                              { value: 'orthographic', label: 'Orthographic' },
                              { value: 'dymaxion', label: 'Dymaxion (Icosahedron) (Experimental)' },
                            ] as SelectOption<ProjectionType>[]}
                            onChange={setExpProj}
                            label="Export projection"
                            className="w-full"
                            triggerClassName="w-full justify-between bg-surface-raised border-edge px-2 py-2 text-xs"
                        />
                    </div>

                    {expProj === 'dymaxion' && (
                        <div className="border border-edge-subtle p-3 space-y-3 bg-surface/40">
                            <div className="flex items-center justify-between">
                                <div className="text-xs font-semibold text-ink-soft">Dymaxion Controls</div>
                                <label className="flex items-center gap-2 text-[10px] text-ink-muted">
                                    <input
                                        type="checkbox"
                                        checked={dymaxionSettings.showOverlay}
                                        onChange={(e) => { updateDymaxion({ showOverlay: e.target.checked }); }}
                                        className="accent-brand-soft"
                                    />
                                    Show Overlay
                                </label>
                            </div>

                            <div className="space-y-1">
                                <label className="text-xs text-ink-muted">Manipulation Mode</label>
                                <div className="grid grid-cols-2 gap-2">
                                    <button
                                        onClick={() => { updateDymaxion({ mode: 'planet' as DymaxionControlMode }); }}
                                        className={`text-[10px] py-2 border ${dymaxionSettings.mode === 'planet' ? 'bg-brand-deep/70 border-brand text-ink-strong' : 'bg-surface-raised border-edge text-ink-soft'}`}
                                    >
                                        Rotate Planet
                                    </button>
                                    <button
                                        onClick={() => { updateDymaxion({ mode: 'overlay' as DymaxionControlMode }); }}
                                        className={`text-[10px] py-2 border ${dymaxionSettings.mode === 'overlay' ? 'bg-brand-deep/70 border-brand text-ink-strong' : 'bg-surface-raised border-edge text-ink-soft'}`}
                                    >
                                        Rotate Overlay
                                    </button>
                                </div>
                                <div className="text-[10px] text-ink-muted">
                                    Drag the globe to rotate. Hold Shift while dragging to roll the overlay.
                                </div>
                            </div>

                            <div className="space-y-2">
                                <div className="flex justify-between text-xs text-ink-muted">
                                    <label>Longitude</label>
                                    <span>{dymaxionSettings.lon}°</span>
                                </div>
                                <input
                                    type="range" min="-180" max="180" step="1"
                                    value={dymaxionSettings.lon}
                                    onChange={function(e) { updateDymaxion({ lon: parseInt(e.target.value) }); }}
                                    className="w-full h-1 bg-surface-hover appearance-none cursor-pointer accent-brand-soft"
                                />
                            </div>

                            <div className="space-y-2">
                                <div className="flex justify-between text-xs text-ink-muted">
                                    <label>Latitude</label>
                                    <span>{dymaxionSettings.lat}°</span>
                                </div>
                                <input
                                    type="range" min="-90" max="90" step="1"
                                    value={dymaxionSettings.lat}
                                    onChange={function(e) { updateDymaxion({ lat: parseInt(e.target.value) }); }}
                                    className="w-full h-1 bg-surface-hover appearance-none cursor-pointer accent-brand-soft"
                                />
                            </div>

                            <div className="space-y-2">
                                <div className="flex justify-between text-xs text-ink-muted">
                                    <label>Roll</label>
                                    <span>{dymaxionSettings.roll}°</span>
                                </div>
                                <input
                                    type="range" min="-180" max="180" step="1"
                                    value={dymaxionSettings.roll}
                                    onChange={function(e) { updateDymaxion({ roll: parseInt(e.target.value) }); }}
                                    className="w-full h-1 bg-surface-hover appearance-none cursor-pointer accent-brand-soft"
                                />
                            </div>

                            <label className="flex items-center gap-2 text-[10px] text-ink-muted">
                                <input
                                    type="checkbox"
                                    checked={dymaxionSettings.layout === 'blender'}
                                    onChange={(e) => { updateDymaxion({ layout: e.target.checked ? 'blender' : 'classic' }); }}
                                    className="accent-brand-soft"
                                />
                                Blender UV Net (export only)
                            </label>

                            <div className="space-y-1">
                                <label className="text-[10px] text-ink-muted">Orientation Presets</label>
                                <div className="grid grid-cols-4 gap-1">
                                    {([['N.Pole', 0, -90, 0], ['Pacific', -150, 0, 0], ['Atlantic', 0, 0, 0], ['Asia', 90, 0, 0]] as const).map(([label, lon, lat, roll]) => (
                                        <button
                                            key={label}
                                            onClick={() => { updateDymaxion({ lon, lat, roll }); }}
                                            className="text-[10px] bg-surface-raised hover:bg-surface-hover text-ink-soft py-1.5 border border-edge"
                                        >
                                            {label}
                                        </button>
                                    ))}
                                </div>
                            </div>

                            <button
                                onClick={() => { updateDymaxion({ lon: 0, lat: 0, roll: 0 }); }}
                                className="w-full text-[10px] bg-surface-raised hover:bg-surface-hover text-ink py-2 border border-edge"
                            >
                                Reset Orientation
                            </button>

                            <label className="flex items-center gap-2 text-[10px] text-ink-muted">
                                <input
                                    type="checkbox"
                                    checked={showDymaxion2D}
                                    onChange={(e) => { setShowDymaxion2D(e.target.checked); }}
                                    className="accent-brand-soft"
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
                        className="w-full flex items-center justify-center gap-2 bg-green-700 hover:bg-green-600 text-ink-strong py-2 text-xs mt-2 disabled:opacity-50 border border-green-600"
                    >
                        <Image size={14}/> Download PNG
                    </button>
                    <button
                        onClick={() => { if (worldData) void exportMap(worldData, 'height_bw', expRes, 'equirectangular'); }}
                        disabled={!worldData}
                        className="w-full flex items-center justify-center gap-2 bg-surface-hover hover:bg-gray-600 text-ink-strong py-2 text-xs disabled:opacity-50 border border-edge-strong"
                    >
                        <Mountain size={14}/> Export Heightmap (BW)
                    </button>
                </div>

                <div className="border-t border-edge-subtle pt-4 space-y-2">
                    <h3 className="text-xs font-bold text-ink-muted uppercase tracking-wider mb-2">Vector Export</h3>
                    <p className="text-[10px] text-ink-muted">Editable coastlines, borders, rivers, and labels for Inkscape/Illustrator, or geodesic GeoJSON for QGIS/web-GIS.</p>
                    <button
                        onClick={() => { if (worldData && expProj !== 'dymaxion') downloadSVG(worldData, viewMode, expProj); }}
                        disabled={!worldData || expProj === 'dymaxion'}
                        title={expProj === 'dymaxion' ? 'SVG export is raster-only for Dymaxion — choose another projection' : undefined}
                        className="w-full flex items-center justify-center gap-2 bg-teal-700 hover:bg-teal-600 text-ink-strong py-2 text-xs disabled:opacity-50 border border-teal-600"
                    >
                        <FileCode size={14}/> Download SVG
                    </button>
                    <button
                        onClick={() => { if (worldData) downloadGeoJSON(worldData); }}
                        disabled={!worldData}
                        className="w-full flex items-center justify-center gap-2 bg-teal-700 hover:bg-teal-600 text-ink-strong py-2 text-xs disabled:opacity-50 border border-teal-600"
                    >
                        <FileJson size={14}/> Download GeoJSON
                    </button>
                </div>

                <div className="border-t border-edge-subtle pt-4 space-y-2">
                    <h3 className="text-xs font-bold text-ink-muted uppercase tracking-wider mb-2">3D Export</h3>
                    <p className="text-[10px] text-ink-muted">Exports the current view as a GLB file. World mesh uses per-vertex colors. Rivers exported as line geometry. City markers included when civilization data is present.</p>
                    <button
                        onClick={() => { if (worldData) exportGLB(worldData, viewMode); }}
                        disabled={!worldData}
                        className="w-full flex items-center justify-center gap-2 bg-indigo-700 hover:bg-indigo-600 text-ink-strong py-2 text-xs disabled:opacity-50 border border-indigo-600"
                    >
                        <Box size={14}/> Export GLB
                    </button>
                </div>

                <div className="border-t border-edge-subtle pt-4 space-y-3">
                    <h3 className="text-xs font-bold text-ink-muted uppercase tracking-wider">File Management</h3>
                    
                    <button
                        onClick={() => { if (params) { void saveMapConfig(params, worldData || undefined); } }}
                        className="w-full flex items-center justify-center gap-2 bg-surface-raised hover:bg-surface-hover text-ink-strong py-2 text-xs border border-edge"
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
                        <button className="w-full flex items-center justify-center gap-2 bg-surface-raised hover:bg-surface-hover text-ink-strong py-2 text-xs pointer-events-none border border-edge">
                            <FolderOpen size={14} /> Load Config (JSON)
                        </button>
                    </div>
                </div>

                <div className="border-t border-edge-subtle pt-4 space-y-3">
                    <h3 className="text-xs font-bold text-ink-muted uppercase tracking-wider">Browser Storage</h3>
                    
                    <div className="flex gap-2">
                        <input 
                            type="text" 
                            placeholder="Save Name..." 
                            value={saveName}
                            onChange={(e) => { setSaveName(e.target.value); }}
                            className="flex-1 bg-surface border border-edge px-2 text-xs text-ink-strong"
                        />
                        <button 
                            onClick={handleSaveBrowser}
                            disabled={!saveName}
                            aria-label="Save map to browser storage"
                            className="bg-brand-strong hover:bg-brand text-ink-strong px-3"
                        >
                            <Save size={14}/>
                        </button>
                    </div>

                    <div className="space-y-1 max-h-40 overflow-y-auto">
                        {savedMaps.length === 0 && <p className="text-xs text-ink-faint italic">No saved maps.</p>}
                        {savedMaps.map(entry => (
                            <div key={entry.name} className="flex items-center justify-between bg-surface p-2 border border-edge-subtle group">
                                <div className="flex flex-col">
                                    <span className="text-xs font-bold text-ink-soft">{entry.name}</span>
                                    <span className="text-[10px] text-ink-muted">{new Date(entry.date).toLocaleDateString()}</span>
                                </div>
                                <div className="flex gap-1 opacity-0 group-hover:opacity-100 transition-opacity">
                                    {/* Named per ENTRY, not just "Load": the list repeats these
                                        two icons per row, so a bare "load button" leaves a reader
                                        no way to tell which saved map it is on. */}
                                    <button onClick={() => { handleLoadBrowser(entry.params, entry.civData, entry.markers); }} aria-label={`Load saved map ${entry.name}`} className="text-brand-soft hover:text-ink-strong p-1"><FolderOpen size={12}/></button>
                                    <button onClick={() => { handleDeleteBrowser(entry.name); }} aria-label={`Delete saved map ${entry.name}`} className="text-danger-soft hover:text-ink-strong p-1"><Trash2 size={12}/></button>
                                </div>
                            </div>
                        ))}
                    </div>
                </div>
            </div>
        )}
      </div>

      <div className="p-4 border-t border-edge-subtle space-y-2">
         {/* Console Output area */}
         <div className="mb-2">
             {/* A real <button>, not a clickable <div>: as a div this was not
                 focusable, had no role, and could not be operated from the
                 keyboard at all. The text content supplies the accessible name. */}
             <button
               type="button"
               className="w-full flex items-center justify-between text-xs text-ink-muted mb-1 cursor-pointer hover:text-ink-soft"
               onClick={() => { setConsoleOpen(!consoleOpen); }}
               aria-expanded={consoleOpen}
             >
                 <span className="flex items-center gap-1">
                    <Terminal size={10} />
                    <span>System Console</span>
                 </span>
                 {consoleOpen ? <ChevronDown size={10}/> : <ChevronUp size={10}/>}
             </button>
             <ConsoleOutput logs={logs} isOpen={consoleOpen} />
         </div>

         {/* Generation progress bar */}
         {loading && (
           <div className="w-full h-1 bg-surface-raised overflow-hidden">
             <div
               className="h-full bg-brand transition-all duration-300"
               style={{ width: `${genProgress * 100}%` }}
             />
           </div>
         )}

         {!loading ? (
             <button
              onClick={handleGenerateClick}
              className={`w-full py-3 font-semibold flex items-center justify-center gap-2 transition-all relative overflow-hidden bg-brand-strong hover:bg-brand text-ink-strong border border-brand`}
            >
              <div className="relative flex items-center gap-2 z-overlay">
                  <RefreshCw size={16} />
                  Generate World
              </div>
            </button>
         ) : (
            <button
              onClick={onCancel}
              className={`w-full py-3 font-semibold flex items-center justify-center gap-2 transition-all relative overflow-hidden bg-red-600 hover:bg-danger text-ink-strong border border-danger`}
            >
              <div className="relative flex items-center gap-2 z-overlay">
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
