import React, { useState, useEffect, useCallback, useRef } from 'react';
import Controls from './components/Controls';
import WorldViewer from './components/WorldViewer';
import Map2D from './components/Map2D';
import MiniMap from './components/MiniMap';
import Inspector from './components/Inspector';
import { BiomeLegend } from './components/Legend';
import { WorldData, WorldParams, ViewMode, LoreData, CivData, DisplayMode, InspectMode, DymaxionSettings } from './types';
import { generateWorld, recalculateCivs, recalculateProvinces } from './utils/worldGen';
import { generateWorldLore, setRuntimeApiKey } from './services/gemini';
import { Menu, X } from 'lucide-react';

const DEFAULT_PARAMS: WorldParams = {
  mapName: 'map',
  points: 5000,
  planetRadius: 6371, 
  axialTilt: 23.5,
  plates: 12,
  seaLevel: 0.55,
  roughness: 0.5,
  detailLevel: 2, 
  landStyle: 'Continents',
  cellJitter: 0.5,
  noiseScale: 0.4,
  ridgeBlend: 0.1,
  maskType: 'None',
  warpStrength: 0.5,
  plateInfluence: 0.5,
  erosionIterations: 2,
  baseTemperature: 30, 
  poleTemperature: -30, 
  rainfallMultiplier: 1.0,
  moistureTransport: 0.5,
  temperatureVariance: 5,
  numFactions: 6,
  civSeed: 'realmgenesis_civs',
  borderRoughness: 0.2, 
  civSizeVariance: 0.5,
  waterCrossingCost: 0.8,
  territorialWaters: 0.15,
  capitalSpacing: 0.5,
  provinceSize: 0.5,
  loreLevel: 1,
  seed: 'realmgenesis',
};

const App: React.FC = () => {
  const [params, setParams] = useState<WorldParams>(DEFAULT_PARAMS);
  const [world, setWorld] = useState<WorldData | null>(null);
  const [viewMode, setViewMode] = useState<ViewMode>('biome');
  const [displayMode, setDisplayMode] = useState<DisplayMode>('globe');
  const [inspectMode, setInspectMode] = useState<InspectMode>('click');
  const [inspectorCollapsed, setInspectorCollapsed] = useState(false);
  const [inspectedCellId, setInspectedCellId] = useState<number | null>(null);
  const [isGenerating, setIsGenerating] = useState(false);
  const [logs, setLogs] = useState<string[]>([]);
  const [lore, setLore] = useState<LoreData | null>(null);
  const [isLoreLoading, setIsLoreLoading] = useState(false);
  const [showGrid, setShowGrid] = useState(false);
  const [showRivers, setShowRivers] = useState(true);
  const [sidebarOpen, setSidebarOpen] = useState(true);
  const [dymaxionSettings, setDymaxionSettings] = useState<DymaxionSettings>({
    layout: 'classic',
    lon: 0,
    lat: 0,
    roll: 0,
    showOverlay: false,
    mode: 'planet',
  });

  const [apiKey, setApiKey] = useState('');

  useEffect(() => {
    setRuntimeApiKey(apiKey);
  }, [apiKey]);

  // Controller reference to persist across renders
  const abortControllerRef = useRef<AbortController | null>(null);

  const addLog = useCallback((msg: string) => {
      setLogs(prev => {
          const lastLog = prev[prev.length - 1];
          // Check for repetitive progress messages to update in-place
          if (lastLog && lastLog.startsWith("Rivers: Drainage processed") && msg.startsWith("Rivers: Drainage processed")) {
              const newLogs = [...prev];
              newLogs[newLogs.length - 1] = msg;
              return newLogs;
          }
          // Standard append, keeping history limit
          return [...prev.slice(-49), msg];
      });
  }, []);

  const handleGenerate = useCallback(async (overrideParams?: WorldParams) => {
    // Abort previous if running
    if (abortControllerRef.current) {
        abortControllerRef.current.abort();
    }
    
    // Create new controller
    const controller = new AbortController();
    abortControllerRef.current = controller;

    setIsGenerating(true);
    setLore(null);
    setLogs(['--- Starting Generation ---']);
    const p = overrideParams || params;
    
    // Defer execution to let UI update
    await new Promise(r => setTimeout(r, 100));
    
    try {
        const newWorld = await generateWorld(p, (msg) => { addLog(msg); }, controller.signal);
        setWorld(newWorld);
        addLog('World Generation Complete.');
    } catch (e: any) { // eslint-disable-next-line @typescript-eslint/no-explicit-any 
        if (e.message === "Generation Cancelled") {
            addLog('Cancelled by user.');
        } else {
            console.error("Generation failed", e);
            addLog(`Error: ${(e as Error).message}`);
        }
    }
    finally { 
        // Only clear generating state if this was the active controller
        if (abortControllerRef.current === controller) {
            setIsGenerating(false); 
            abortControllerRef.current = null;
        }
    }
  }, [params, addLog]);

  const handleLoadWorld = useCallback(async (newParams: WorldParams, savedCivData?: CivData) => {
     if (abortControllerRef.current) abortControllerRef.current.abort();
     const controller = new AbortController();
     abortControllerRef.current = controller;

     setIsGenerating(true);
     setLore(null);
     setLogs(['--- Loading Map from File ---']);
     setParams(newParams);

     await new Promise(r => setTimeout(r, 100));

     try {
         // 1. Regenerate World Geometry & Civs based on Seed
         const newWorld = await generateWorld(newParams, (msg) => { addLog(msg); }, controller.signal);
         
         // 2. Restore Saved Metadata (Names, Descriptions, Colors)
         if (savedCivData && newWorld.civData) {
              addLog("Restoring historical records...");
              savedCivData.factions.forEach(savedF => {
                  // Match by ID since seed is identical
                  const genF = newWorld.civData?.factions.find(f => f.id === savedF.id);
                  if (genF) {
                      genF.name = savedF.name;
                      genF.color = savedF.color;
                      genF.description = savedF.description;
                      
                      // Restore province names
                      savedF.provinces.forEach((savedP, idx) => {
                          // eslint-disable-next-line security/detect-object-injection
                          if (genF.provinces[idx]) {
                              // eslint-disable-next-line security/detect-object-injection
                              genF.provinces[idx].name = savedP.name;
                              // Restore town names
                              savedP.towns.forEach((savedT, tIdx) => {
                                  // eslint-disable-next-line security/detect-object-injection
                                  if (genF.provinces[idx].towns[tIdx]) {
                                      // eslint-disable-next-line security/detect-object-injection
                                      genF.provinces[idx].towns[tIdx].name = savedT.name;
                                  }
                              });
                          }
                      });
                  }
              });
         }

         setWorld(newWorld);
         addLog('Map Loaded Successfully.');
     } catch (e: any) { // eslint-disable-next-line @typescript-eslint/no-explicit-any
        if (e.message === "Generation Cancelled") addLog('Cancelled.');
        else {
            console.error("Load failed", e);
            addLog(`Load Error: ${(e as Error).message}`);
        }
     } finally {
         if (abortControllerRef.current === controller) {
             setIsGenerating(false);
             abortControllerRef.current = null;
         }
     }
  }, [addLog]);

  const handleCancel = useCallback(() => {
      if (abortControllerRef.current) {
          abortControllerRef.current.abort();
          addLog("Cancelling...");
      }
  }, [addLog]);

  const handleUpdateCivs = useCallback(async (overrideParams?: WorldParams) => {
      if (!world) return;
      setIsGenerating(true); 
      addLog('--- Updating Civilizations ---');
      const p = overrideParams || params;
      await new Promise(r => setTimeout(r, 50));
      try {
          const updatedWorld = recalculateCivs(world, p, (msg) => { addLog(msg); });
          setWorld({ ...updatedWorld });
          if (viewMode !== 'political') setViewMode('political');
          addLog('Civilizations Updated.');
      } catch(e) { 
          console.error("Civ update failed", e); 
          addLog(`Error: ${(e as Error).message}`);
      }
      finally { setIsGenerating(false); }
  }, [world, params, viewMode, addLog]);

  const handleUpdateProvinces = useCallback(async (overrideParams?: WorldParams) => {
    if (!world || !world.civData) return;
    setIsGenerating(true); 
    addLog('--- Updating Provinces ---');
    const p = overrideParams || params;
    await new Promise(r => setTimeout(r, 50));
    try {
        const updatedWorld = recalculateProvinces(world, p);
        setWorld({ ...updatedWorld });
        if (viewMode !== 'political') setViewMode('political');
        addLog('Provinces Updated.');
    } catch(e) { 
        console.error("Province update failed", e);
        addLog(`Error: ${(e as Error).message}`); 
    }
    finally { setIsGenerating(false); }
  }, [world, params, viewMode, addLog]);

  useEffect(() => {
    setInspectedCellId(null);
  }, [world?.params.seed, displayMode]);

  const toggleInspectEnabled = useCallback(() => {
    setInspectMode(prev => (prev === 'off' ? 'click' : 'off'));
    if (inspectMode === 'click') {
      setInspectedCellId(null);
    }
  }, [inspectMode]);

  useEffect(() => { handleGenerate(); }, []);

  const handleGenerateLore = async () => {
    if (!world) return;
    setIsLoreLoading(true);
    addLog('Contacting Gemini API for lore...');
    try {
      const newLore = await generateWorldLore(world);
      setLore(newLore);
      setWorld({ ...world });
      addLog('Lore Received.');
    } catch (e: any) { // eslint-disable-next-line @typescript-eslint/no-explicit-any 
        console.error("Lore gen failed", e);
        addLog(`Lore Error: ${e.message}`);
    }
    finally { setIsLoreLoading(false); }
  };

  return (
    <div className="flex flex-col md:flex-row w-full h-full bg-black overflow-hidden font-sans text-gray-200">
      {/* Sidebar / Bottom Drawer */}
      <div className={`fixed inset-x-0 bottom-0 z-30 md:relative md:inset-auto md:w-80 md:h-full
 bg-gray-950 border-t md:border-t-0 md:border-r border-gray-800 transition-transform duration-300
 ${sidebarOpen ? 'translate-y-0 md:translate-x-0' : 'translate-y-full md:-translate-x-full'}
 max-h-[85vh] md:max-h-full flex flex-col shadow-2xl`}>
        <Controls 
          params={params} setParams={setParams}
          onGenerate={handleGenerate} 
          onLoadWorld={handleLoadWorld}
          onCancel={handleCancel}
          onUpdateCivs={handleUpdateCivs} onUpdateProvinces={handleUpdateProvinces}
          viewMode={viewMode} setViewMode={setViewMode}
          displayMode={displayMode} setDisplayMode={setDisplayMode}
          loading={isGenerating} logs={logs}
          lore={lore} generatingLore={isLoreLoading} onGenerateLore={handleGenerateLore}
          worldData={world} 
          showGrid={showGrid} setShowGrid={setShowGrid}
          showRivers={showRivers} setShowRivers={setShowRivers}
          dymaxionSettings={dymaxionSettings}
          onDymaxionChange={setDymaxionSettings}
          apiKey={apiKey}
          onApiKeyChange={setApiKey}
        />
        <button 
          onClick={() => { setSidebarOpen(false); }}
          className="md:hidden absolute -top-12 right-4 bg-gray-900 text-white p-2 border border-gray-700 shadow-lg"
        >
          <X size={20} />
        </button>
      </div>

      {/* Floating menu button - Top Left to avoid overlapping with bottom overlays */}
      {!sidebarOpen && (
        <button 
          onClick={() => { setSidebarOpen(true); }}
          className="fixed top-4 left-4 z-40 bg-blue-600 text-white p-3 shadow-2xl hover:bg-blue-500 md:hidden border border-white/10"
        >
          <Menu size={24} />
        </button>
      )}

      {/* Main Content Area */}
      <main className="flex-1 relative h-full overflow-hidden">
        {displayMode === 'globe' ? (
          <WorldViewer
            world={world}
            viewMode={viewMode}
            showGrid={showGrid}
            showRivers={showRivers}
            inspectMode={inspectMode}
            onInspect={setInspectedCellId}
            dymaxionSettings={dymaxionSettings}
            onDymaxionChange={setDymaxionSettings}
          />
        ) : (
          <Map2D
            world={world}
            viewMode={viewMode}
            inspectMode={inspectMode}
            onInspect={setInspectedCellId}
            projectionType={displayMode === 'dymaxion' ? 'dymaxion' : 'mercator'}
            dymaxionSettings={dymaxionSettings}
            showGrid={showGrid}
            showRivers={showRivers}
          />
        )}

        {/* Overlay HUD elements */}
        {displayMode === 'globe' && (
          <div className="absolute top-4 left-24 bg-black/50 backdrop-blur-md p-3 border border-white/10 text-left pointer-events-none z-10 hidden md:block">
           <h3 className="text-white text-xs font-bold">{world ? `Seed: ${params.seed}` : 'No World'}</h3>
           <p className="text-gray-400 text-[10px]">{world ? `${world.cells.length.toLocaleString()} Cells` : ''}</p>
          </div>
        )}

        <BiomeLegend />
        {displayMode === 'globe' && <MiniMap world={world} viewMode={viewMode} />}
        <Inspector
          world={world}
          cellId={inspectedCellId}
          inspectMode={inspectMode}
          collapsed={inspectorCollapsed}
          onToggleEnabled={toggleInspectEnabled}
          onToggleCollapsed={() => { setInspectorCollapsed(v => !v); }}
        />
      </main>
    </div>
  );
};

export default App;
