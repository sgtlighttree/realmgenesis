import React, { useState, useEffect, useCallback, useRef, useMemo } from 'react';
import Controls from './components/Controls';
import WorldViewer from './components/WorldViewer';
import Map2D from './components/Map2D';
import MiniMap from './components/MiniMap';
import Inspector from './components/Inspector';
import { BiomeLegend } from './components/Legend';
import { WorldData, WorldParams, ViewMode, LoreData, CivData, DisplayMode, InspectMode, DymaxionSettings, EditMode, PaintStyle, UndoSnapshot, BiomeType, POLITICAL_ERASER_ID, LabelVisibility, DEFAULT_LABEL_VISIBILITY } from './types';
import { generateWorld, recalculateCivs, recalculateProvinces } from './utils/worldGen';
import { getCellsInRadius, applyTerrainStroke, applyFlattenStroke, applySmoothStroke, applyPoliticalStroke, applyBiomeStroke, refreshBiomes } from './utils/paintUtils';
import { generateWorldLore, setRuntimeApiKey } from './services/gemini';
import EditToolbar from './components/EditToolbar';
import { Menu, X } from 'lucide-react';

const DEFAULT_PARAMS: WorldParams = {
  mapName: 'map',
  points: 5000,
  planetRadius: 6371, 
  axialTilt: 23.5,
  plates: 12,
  seaLevel: 0.55,
  roughness: 0.5,
  detailLevel: 3, // FBM octaves; 3 = historical default terrain character
  landStyle: 'Continents',
  cellJitter: 0.5,
  noiseScale: 0.4,
  ridgeBlend: 0.1,
  mountainHeight: 1.0,
  oceanDepth: 1.0,
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

const isValidProvinceId = (world: WorldData, factionId: number, provinceId: number | undefined): provinceId is number => {
  if (provinceId === undefined) return false;
  const faction = world.civData?.factions.find(f => f.id === factionId);
  return faction?.provinces.some(p => p.id === provinceId) ?? false;
};

const resolvePoliticalProvinceId = (world: WorldData, startCellId: number, factionId: number): number | undefined => {
  const startCell = world.cells[startCellId];
  if (!startCell) return undefined;
  if (startCell.regionId === factionId && isValidProvinceId(world, factionId, startCell.provinceId)) {
    return startCell.provinceId;
  }

  let bestProvinceId: number | undefined;
  let bestDot = -Infinity;
  for (const cell of world.cells) {
    if (cell.regionId !== factionId || !isValidProvinceId(world, factionId, cell.provinceId)) continue;
    const d = startCell.center.x * cell.center.x + startCell.center.y * cell.center.y + startCell.center.z * cell.center.z;
    if (d > bestDot) {
      bestDot = d;
      bestProvinceId = cell.provinceId;
    }
  }

  if (bestProvinceId !== undefined) return bestProvinceId;
  return world.civData?.factions.find(f => f.id === factionId)?.provinces[0]?.id;
};

const recalculatePoliticalTotals = (world: WorldData): void => {
  if (!world.civData) return;
  const factionById = new Map(world.civData.factions.map(f => [f.id, f]));
  world.civData.factions.forEach(faction => {
    faction.totalPopulation = 0;
    faction.provinces.forEach(province => { province.totalPopulation = 0; });
  });

  world.cells.forEach(cell => {
    if (cell.regionId === undefined || cell.provinceId === undefined) return;
    const faction = factionById.get(cell.regionId);
    const province = faction?.provinces.find(p => p.id === cell.provinceId);
    if (!faction || !province) return;
    province.totalPopulation += cell.population || 0;
  });

  world.civData.factions.forEach(faction => {
    faction.totalPopulation = faction.provinces.reduce((sum, province) => sum + province.totalPopulation, 0);
  });
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
  const [genProgress, setGenProgress] = useState(0);
  const [logs, setLogs] = useState<string[]>([]);
  const [lore, setLore] = useState<LoreData | null>(null);
  const [isLoreLoading, setIsLoreLoading] = useState(false);
  const [showGrid, setShowGrid] = useState(false);
  const [showRivers, setShowRivers] = useState(true);
  const [labelVisibility, setLabelVisibility] = useState<LabelVisibility>(DEFAULT_LABEL_VISIBILITY);
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

  // Edit mode state
  const [editMode, setEditMode] = useState<EditMode>('off');
  const [paintStyle, setPaintStyle] = useState<PaintStyle>('adaptive');
  const [brushSize, setBrushSize] = useState<number>(1);
  const [paintStrength, setPaintStrength] = useState<number>(0.5);
  const [paintFaction, setPaintFaction] = useState<number>(0);
  const [paintBiome, setPaintBiome] = useState<BiomeType>(BiomeType.OCEAN);
  const [sampleHeight, setSampleHeight] = useState<number | null>(null);
  const [adaptiveBiomes, setAdaptiveBiomes] = useState<boolean>(true);
  const [undoStack, setUndoStack] = useState<UndoSnapshot[]>([]);
  const isPainting = useRef<boolean>(false);
  const currentStrokeSnapshot = useRef<UndoSnapshot['cells'] | null>(null);
  const currentPoliticalProvince = useRef<number | undefined>(undefined);

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
    setGenProgress(0);
    setLore(null);
    setLogs(['--- Starting Generation ---']);
    const p = overrideParams || params;

    // Defer execution to let UI update
    await new Promise(r => setTimeout(r, 100));

    try {
        const newWorld = await generateWorld(p, (msg) => { addLog(msg); }, controller.signal, (s, t) => setGenProgress(s / t));
        setWorld(newWorld);
        setGenProgress(1);
        addLog('World Generation Complete.');
    } catch (e: any) {  
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
         // Wrapped separately: corrupt metadata degrades to "terrain only"
         // instead of failing the whole load
         try {
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
                          if (genF.provinces[idx]) {
                              genF.provinces[idx].name = savedP.name;
                              savedP.towns.forEach((savedT, tIdx) => {
                                  if (genF.provinces[idx].towns[tIdx]) {
                                      genF.provinces[idx].towns[tIdx].name = savedT.name;
                                  }
                              });
                          }
                      });
                  }
              });
         }
         } catch (metaErr) {
             console.error("civData restore failed", metaErr);
             addLog('Warning: saved names/colors could not be restored; loaded terrain only.');
         }

         setWorld(newWorld);
         addLog('Map Loaded Successfully.');
     } catch (e: any) {  
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
    } catch (e: any) {  
        console.error("Lore gen failed", e);
        addLog(`Lore Error: ${e.message}`);
    }
    finally { setIsLoreLoading(false); }
  };

  const factionColors = useMemo<Map<number, string>>(() => {
    const m = new Map<number, string>();
    world?.civData?.factions.forEach(f => m.set(f.id, f.color));
    return m;
  }, [world]);

  const handlePaint = useCallback((cellId: number, phase: 'start' | 'stroke' | 'end', isRightClick = false) => {
    if (!world || editMode === 'off' || editMode === 'world-edit') return;
    const center = world.cells[cellId];
    if (!center) return;
    const brush = getCellsInRadius(center, brushSize, world.cells);

    if (phase === 'start') {
      if (editMode === 'terrain-flatten' && isRightClick) {
        setSampleHeight(center.height);
        return;
      }
      currentPoliticalProvince.current = editMode === 'political' && paintFaction !== POLITICAL_ERASER_ID
        ? resolvePoliticalProvinceId(world, center.id, paintFaction)
        : undefined;
      // Create the snapshot Map and push its reference to undoStack immediately.
      // During 'stroke' we mutate the same Map in place — the stack entry grows automatically.
      // This means 'end' never needs to push, making undo reliable even if 'end' is missed.
      const snapMap: UndoSnapshot['cells'] = new Map(
        brush.map(({ cell }) => [
          cell.id,
          { height: cell.height, biome: cell.biome, regionId: cell.regionId, provinceId: cell.provinceId },
        ])
      );
      currentStrokeSnapshot.current = snapMap;
      setUndoStack(prev => [...prev.slice(-19), { cells: snapMap }]);
      isPainting.current = true;
    }

    if (phase === 'stroke' && isPainting.current) {
      // Add newly-touched cells to the snapshot BEFORE modifying them
      if (currentStrokeSnapshot.current) {
        brush.forEach(({ cell }) => {
          if (!currentStrokeSnapshot.current!.has(cell.id)) {
            currentStrokeSnapshot.current!.set(cell.id, {
              height: cell.height,
              biome: cell.biome,
              regionId: cell.regionId,
              provinceId: cell.provinceId,
            });
          }
        });
      }
      if (editMode === 'terrain-raise' || editMode === 'terrain-lower') {
        applyTerrainStroke(brush, brushSize, editMode === 'terrain-raise' ? 'raise' : 'lower', paintStyle, world.cells, paintStrength);
        if (adaptiveBiomes) refreshBiomes(brush.map(b => b.cell), world.params.seaLevel);
      } else if (editMode === 'terrain-smooth') {
        applySmoothStroke(brush, brushSize, paintStrength, world.cells);
        if (adaptiveBiomes) refreshBiomes(brush.map(b => b.cell), world.params.seaLevel);
      } else if (editMode === 'terrain-flatten' && sampleHeight !== null) {
        applyFlattenStroke(brush, brushSize, sampleHeight, paintStrength);
        if (adaptiveBiomes) refreshBiomes(brush.map(b => b.cell), world.params.seaLevel);
      } else if (editMode === 'political') {
        applyPoliticalStroke(
          brush,
          paintFaction === POLITICAL_ERASER_ID ? undefined : paintFaction,
          currentPoliticalProvince.current,
        );
      } else if (editMode === 'biome') {
        applyBiomeStroke(brush, paintBiome);
      }
      setWorld({ ...world });
    }

    if (phase === 'end') {
      if (editMode === 'political' && isPainting.current) {
        recalculatePoliticalTotals(world);
        setWorld({ ...world });
      }
      currentStrokeSnapshot.current = null; // snapshot already lives in undoStack
      currentPoliticalProvince.current = undefined;
      isPainting.current = false;
    }
  }, [world, editMode, paintStyle, brushSize, paintStrength, paintFaction, paintBiome, sampleHeight, adaptiveBiomes]);

  const handleUndo = useCallback(() => {
    if (!world || undoStack.length === 0) return;
    const snap = undoStack[undoStack.length - 1];
    snap.cells.forEach(({ height, biome, regionId, provinceId }, id) => {
      world.cells[id].height = height;
      world.cells[id].biome = biome;
      world.cells[id].regionId = regionId;
      world.cells[id].provinceId = provinceId;
    });
    recalculatePoliticalTotals(world);
    setWorld({ ...world });
    setUndoStack(prev => prev.slice(0, -1));
  }, [world, undoStack]);

  useEffect(() => {
    const onKey = (e: KeyboardEvent) => {
      if ((e.ctrlKey || e.metaKey) && e.key === 'z' && !e.shiftKey) {
        e.preventDefault();
        handleUndo();
      }
    };
    window.addEventListener('keydown', onKey);
    return () => window.removeEventListener('keydown', onKey);
  }, [handleUndo]);

  const handleEditWorldData = useCallback((cellId: number, updates: {
    townName?: string; townPopulation?: number;
    factionName?: string; factionColor?: string; factionDescription?: string;
    provinceName?: string;
  }) => {
    if (!world?.civData) return;
    const cell = world.cells[cellId];
    const faction = world.civData.factions.find(f => f.id === cell.regionId);
    if (faction) {
      if (updates.factionName !== undefined) faction.name = updates.factionName;
      if (updates.factionColor !== undefined) faction.color = updates.factionColor;
      if (updates.factionDescription !== undefined) faction.description = updates.factionDescription;
      if (cell.provinceId !== undefined) {
        const prov = faction.provinces.find(p => p.id === cell.provinceId);
        if (prov) {
          if (updates.provinceName !== undefined) prov.name = updates.provinceName;
          const town = prov.towns.find(t => t.cellId === cellId);
          if (town) {
            if (updates.townName !== undefined) town.name = updates.townName;
            if (updates.townPopulation !== undefined) town.population = updates.townPopulation;
          }
        }
      }
    }
    setWorld({ ...world });
  }, [world]);

  const handleEditFaction = useCallback((factionId: number, updates: {
    name?: string; color?: string; description?: string;
  }) => {
    if (!world?.civData) return;
    const faction = world.civData.factions.find(f => f.id === factionId);
    if (!faction) return;
    if (updates.name !== undefined) faction.name = updates.name;
    if (updates.color !== undefined) faction.color = updates.color;
    if (updates.description !== undefined) faction.description = updates.description;
    setWorld({ ...world });
  }, [world]);

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
          loading={isGenerating} logs={logs} genProgress={genProgress}
          lore={lore} generatingLore={isLoreLoading} onGenerateLore={handleGenerateLore}
          worldData={world} 
          showGrid={showGrid} setShowGrid={setShowGrid}
          showRivers={showRivers} setShowRivers={setShowRivers}
          labelVisibility={labelVisibility} setLabelVisibility={setLabelVisibility}
          dymaxionSettings={dymaxionSettings}
          onDymaxionChange={setDymaxionSettings}
          apiKey={apiKey}
          onApiKeyChange={setApiKey}
          onInspect={setInspectedCellId}
          onEditFaction={handleEditFaction}
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
            labelVisibility={labelVisibility}
            inspectMode={inspectMode}
            onInspect={setInspectedCellId}
            selectedCellId={inspectedCellId}
            dymaxionSettings={dymaxionSettings}
            onDymaxionChange={setDymaxionSettings}
            editMode={editMode}
            onPaint={handlePaint}
            factionColors={factionColors}
            brushSize={brushSize}
          />
        ) : (
          <Map2D
            world={world}
            viewMode={viewMode}
            inspectMode={inspectMode}
            onInspect={setInspectedCellId}
            highlightCellId={inspectedCellId}
            projectionType={displayMode === 'dymaxion' ? 'dymaxion' : 'mercator'}
            dymaxionSettings={dymaxionSettings}
            showGrid={showGrid}
            showRivers={showRivers}
            labelVisibility={labelVisibility}
            editMode={editMode}
            onPaint={handlePaint}
            factionColors={factionColors}
            brushSize={brushSize}
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
          editMode={editMode}
          onEditWorldData={handleEditWorldData}
        />
        {world && (
          <EditToolbar
            editMode={editMode} setEditMode={setEditMode}
            paintStyle={paintStyle} setPaintStyle={setPaintStyle}
            brushSize={brushSize} setBrushSize={setBrushSize}
            paintStrength={paintStrength} setPaintStrength={setPaintStrength}
            adaptiveBiomes={adaptiveBiomes} setAdaptiveBiomes={setAdaptiveBiomes}
            paintFaction={paintFaction} setPaintFaction={setPaintFaction}
            paintBiome={paintBiome} setPaintBiome={setPaintBiome}
            sampleHeight={sampleHeight}
            undoCount={undoStack.length} onUndo={handleUndo}
            world={world}
          />
        )}
      </main>
    </div>
  );
};

export default App;
