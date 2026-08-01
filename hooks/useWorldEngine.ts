import { useState, useEffect, useCallback, useRef, useMemo } from 'react';
import { WorldData, WorldParams, ViewMode, LoreData, CivData, DisplayMode, InspectMode, DymaxionSettings, EditMode, PaintStyle, UndoSnapshot, BiomeType, POLITICAL_ERASER_ID, LabelVisibility, DEFAULT_LABEL_VISIBILITY, Point, MarkerData } from '../types';
import { recalculateCivs, recalculateProvinces } from '../utils/worldGen';
import { generateWorldInWorker } from '../utils/worldGenClient';
import { computeRoutes } from '../utils/routes';
import { mergeFactions, renameProvince, renameTown, relocateCapital } from '../utils/civEdit';
import { getCellsInRadius, applyTerrainStroke, applyFlattenStroke, applySmoothStroke, applyPoliticalStroke, applyBiomeStroke, refreshBiomes } from '../utils/paintUtils';
import { generateWorldLore, setRuntimeApiKey } from '../services/gemini';
import { greatCircleDistanceKm, sampleGreatCircleArc } from '../utils/measure';

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
  tectonicStrength: 0.5,
  erosionIterations: 2,
  marginCoupling: 0.3,
  numTimesteps: 20,
  simulationResolution: 10000,
  plateJitter: 0.3,
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
  numCultures: 4,
  nameStyle: 'fantasy',
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

export function useWorldEngine() {
  const [params, setParams] = useState<WorldParams>(DEFAULT_PARAMS);
  const [world, setWorld] = useState<WorldData | null>(null);
  const [viewMode, setViewMode] = useState<ViewMode>('biome');
  const [displayMode, setDisplayMode] = useState<DisplayMode>('globe');
  const [inspectMode, setInspectMode] = useState<InspectMode>('click');
  const [inspectorCollapsed, setInspectorCollapsed] = useState(false);
  const [inspectedCellId, setInspectedCellId] = useState<number | null>(null);
  const [rulerActive, setRulerActive] = useState(false);
  const [rulerCells, setRulerCells] = useState<[number, number | null] | null>(null);
  const [markerMode, setMarkerMode] = useState(false);
  const [selectedMarkerId, setSelectedMarkerId] = useState<number | null>(null);
  const [isGenerating, setIsGenerating] = useState(false);
  const [genProgress, setGenProgress] = useState(0);
  const [logs, setLogs] = useState<string[]>([]);
  const [lore, setLore] = useState<LoreData | null>(null);
  const [isLoreLoading, setIsLoreLoading] = useState(false);
  const [showGrid, setShowGrid] = useState(false);
  const [showRivers, setShowRivers] = useState(true);
  const [showRoutes, setShowRoutes] = useState(false);
  const [showHillshade, setShowHillshade] = useState(false);
  const [showContours, setShowContours] = useState(false);
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
    setUndoStack([]);
    setLogs(['--- Starting Generation ---']);
    const p = overrideParams || params;

    // The main thread no longer blocks during generation (it runs in a worker),
    // so this delay is no longer needed to let the UI update. It is kept as an
    // intentional debounce: a rapid second handleGenerate call inside this 100ms
    // window aborts the first before any worker is constructed, which suppresses
    // worker churn during rapid slider drags. Removing it trades that off for
    // lower latency on a single generate call.
    await new Promise(r => setTimeout(r, 100));

    try {
        const newWorld = await generateWorldInWorker(p, (msg) => { addLog(msg); }, controller.signal, (s, t) => setGenProgress(s / t));
        // Markers are sphere-anchored (not cellId-based), so they stay
        // meaningful across regeneration — carry them over from the world
        // being replaced. Functional updater avoids a stale `world` closure.
        setWorld(prev => {
          newWorld.markers = prev?.markers;
          return newWorld;
        });
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

  const handleLoadWorld = useCallback(async (newParams: WorldParams, savedCivData?: CivData, savedMarkers?: MarkerData[]) => {
     if (abortControllerRef.current) abortControllerRef.current.abort();
     const controller = new AbortController();
     abortControllerRef.current = controller;

     setIsGenerating(true);
     setLore(null);
     setUndoStack([]);
     setLogs(['--- Loading Map from File ---']);
     setParams(newParams);

     await new Promise(r => setTimeout(r, 100));

     try {
         // 1. Regenerate World Geometry & Civs based on Seed
         const newWorld = await generateWorldInWorker(newParams, (msg) => { addLog(msg); }, controller.signal);
         
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

         newWorld.markers = savedMarkers;
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
    setRulerCells(null);
  }, [world?.params.seed, displayMode]);

  // C3: routes are computed LAZILY here — only when the Roads & Routes toggle
  // is on and the current world has none cached. computeRoutes is O(towns · A*)
  // and runs several seconds near the 200k-cell cap, so gating it keeps a
  // routes-off generation (the default) free. A fresh generate / civ-update
  // clears world.routes (see recalculateProvinces), so this recomputes against
  // the new town graph. Mutating + shallow-copying matches the paint pattern.
  useEffect(() => {
    if (!showRoutes || !world || !world.civData || world.routes) return;
    let cancelled = false;
    setIsGenerating(true);
    addLog('Computing roads & trade routes...');
    const id = setTimeout(() => {
      if (cancelled) return;
      const routes = computeRoutes(world, world.params);
      if (cancelled) return;
      world.routes = routes;
      setWorld({ ...world });
      addLog(`Routes: ${routes.filter(r => r.kind === 'road').length} roads, ${routes.filter(r => r.kind === 'searoute').length} sea.`);
      setIsGenerating(false);
    }, 30);
    return () => { cancelled = true; clearTimeout(id); setIsGenerating(false); };
  }, [showRoutes, world, addLog]);

  const toggleInspectEnabled = useCallback(() => {
    setInspectMode(prev => (prev === 'off' ? 'click' : 'off'));
    if (inspectMode === 'click') {
      setInspectedCellId(null);
    }
  }, [inspectMode]);

  // Ruler and marker placement are mutually exclusive click-capture modes;
  // both force inspectMode to 'click' for the same reason: the 3D/2D pickers
  // only ever fire onInspect from click handling when inspectMode === 'click'
  // (hover mode fires it continuously on pointer move, which would spam
  // ruler points / drop a marker on every hover).
  const toggleRuler = useCallback(() => {
    setRulerActive(prev => {
      const next = !prev;
      if (next) { setInspectMode('click'); setMarkerMode(false); }
      return next;
    });
    setRulerCells(null);
  }, []);

  const toggleMarkerMode = useCallback(() => {
    setMarkerMode(prev => {
      const next = !prev;
      if (next) { setInspectMode('click'); setRulerActive(false); }
      return next;
    });
  }, []);

  // Intercepts the shared onInspect plumbing: while the ruler or marker mode
  // is active, cell clicks build up a measurement / drop a pin instead of
  // updating the Inspector's selection. Children (WorldViewer/Map2D/Controls)
  // are unaware of either mode — they just keep calling onInspect(cellId).
  const handleInspect = useCallback((cellId: number | null) => {
    if (markerMode) {
      if (cellId === null || !world) return;
      const cell = world.cells[cellId];
      if (!cell) return;
      const markers = world.markers ?? [];
      const nextId = markers.reduce((max, m) => Math.max(max, m.id), -1) + 1;
      const len = Math.hypot(cell.center.x, cell.center.y, cell.center.z) || 1;
      const newMarker: MarkerData = {
        id: nextId,
        kind: 'poi',
        name: `Marker ${nextId}`,
        note: '',
        position: { x: cell.center.x / len, y: cell.center.y / len, z: cell.center.z / len },
      };
      world.markers = [...markers, newMarker];
      setWorld({ ...world });
      setSelectedMarkerId(nextId);
      return; // stay in markerMode for multiple placements
    }
    if (rulerActive) {
      if (cellId === null) return; // misses/hover-outs never touch ruler state
      setRulerCells(prev => (!prev || prev[1] !== null) ? [cellId, null] : [prev[0], cellId]);
      return;
    }
    setInspectedCellId(cellId);
  }, [rulerActive, markerMode, world]);

  const updateMarker = useCallback((id: number, patch: Partial<Pick<MarkerData, 'kind' | 'name' | 'note'>>) => {
    if (!world?.markers) return;
    const marker = world.markers.find(m => m.id === id);
    if (!marker) return;
    Object.assign(marker, patch);
    setWorld({ ...world });
  }, [world]);

  const deleteMarker = useCallback((id: number) => {
    if (!world?.markers) return;
    world.markers = world.markers.filter(m => m.id !== id);
    setWorld({ ...world });
    setSelectedMarkerId(prev => (prev === id ? null : prev));
  }, [world]);

  const rulerArc = useMemo<Point[] | null>(() => {
    if (!world || !rulerCells || rulerCells[1] === null) return null;
    const a = world.cells[rulerCells[0]]?.center;
    const b = world.cells[rulerCells[1]]?.center;
    if (!a || !b) return null;
    return sampleGreatCircleArc(a, b);
  }, [world, rulerCells]);

  const rulerDistanceKm = useMemo<number | null>(() => {
    if (!world || !rulerCells || rulerCells[1] === null) return null;
    const a = world.cells[rulerCells[0]]?.center;
    const b = world.cells[rulerCells[1]]?.center;
    if (!a || !b) return null;
    return greatCircleDistanceKm(a, b, world.params.planetRadius);
  }, [world, rulerCells]);

  useEffect(() => { handleGenerate(); }, []);

  /**
   * Generation is DESTRUCTIVE: it replaces the world, discarding any terrain,
   * biome, or border painting. The undo stack is the only evidence that hand
   * work exists, so it is the gate.
   *
   * The gate lives here rather than in the Generate button so that every entry
   * point inherits it — today that is the button in both shells; a future
   * keyboard shortcut or command-palette entry gets it for free. Auto-update
   * deliberately keeps calling `handleGenerate` directly and stays ungated: it
   * fires on every slider change, so prompting would make it unusable.
   */
  const [pendingGenerate, setPendingGenerate] = useState<{ params?: WorldParams } | null>(null);

  const requestGenerate = useCallback((overrideParams?: WorldParams) => {
    if (undoStack.length > 0) setPendingGenerate({ params: overrideParams });
    else void handleGenerate(overrideParams);
  }, [undoStack.length, handleGenerate]);

  const confirmGenerate = useCallback(() => {
    const p = pendingGenerate;
    setPendingGenerate(null);
    void handleGenerate(p?.params);
  }, [pendingGenerate, handleGenerate]);

  const cancelGenerate = useCallback(() => { setPendingGenerate(null); }, []);

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

  const cultureColors = useMemo<Map<number, string>>(() => {
    const m = new Map<number, string>();
    world?.cultures?.forEach(c => m.set(c.id, c.color));
    return m;
  }, [world]);

  const religionColors = useMemo<Map<number, string>>(() => {
    const m = new Map<number, string>();
    world?.religions?.forEach(r => m.set(r.id, r.color));
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

  const handleMergeFactions = useCallback((srcId: number, dstId: number) => {
    if (!world) return;
    if (mergeFactions(world, srcId, dstId)) {
      recalculatePoliticalTotals(world);
      setWorld({ ...world });
    }
  }, [world]);

  const handleRenameProvince = useCallback((factionId: number, provinceId: number, name: string) => {
    if (!world) return;
    if (renameProvince(world, factionId, provinceId, name)) {
      setWorld({ ...world });
    }
  }, [world]);

  const handleRenameTown = useCallback((factionId: number, provinceId: number, cellId: number, name: string) => {
    if (!world) return;
    if (renameTown(world, factionId, provinceId, cellId, name)) {
      setWorld({ ...world });
    }
  }, [world]);

  const handleRelocateCapital = useCallback((factionId: number, townCellId: number) => {
    if (!world) return;
    if (relocateCapital(world, factionId, townCellId)) {
      setWorld({ ...world });
    }
  }, [world]);
  return {
    params, setParams, world, setWorld, viewMode, setViewMode, displayMode, setDisplayMode, inspectMode, setInspectMode, inspectorCollapsed, setInspectorCollapsed, inspectedCellId, setInspectedCellId, rulerActive, setRulerActive, rulerCells, setRulerCells, markerMode, setMarkerMode, selectedMarkerId, setSelectedMarkerId, isGenerating, setIsGenerating, genProgress, setGenProgress, logs, setLogs, lore, setLore, isLoreLoading, setIsLoreLoading, showGrid, setShowGrid, showRivers, setShowRivers, showRoutes, setShowRoutes, showHillshade, setShowHillshade, showContours, setShowContours, labelVisibility, setLabelVisibility, sidebarOpen, setSidebarOpen, dymaxionSettings, setDymaxionSettings, apiKey, setApiKey, editMode, setEditMode, paintStyle, setPaintStyle, brushSize, setBrushSize, paintStrength, setPaintStrength, paintFaction, setPaintFaction, paintBiome, setPaintBiome, sampleHeight, setSampleHeight, adaptiveBiomes, setAdaptiveBiomes, undoStack, setUndoStack, addLog, handleGenerate, requestGenerate, pendingGenerate, confirmGenerate, cancelGenerate, handleLoadWorld, handleCancel, handleUpdateCivs, handleUpdateProvinces, toggleInspectEnabled, toggleRuler, toggleMarkerMode, handleInspect, updateMarker, deleteMarker, rulerArc, rulerDistanceKm, handleGenerateLore, factionColors, cultureColors, religionColors, handlePaint, handleUndo, handleEditWorldData, handleEditFaction, handleMergeFactions, handleRenameProvince, handleRenameTown, handleRelocateCapital,
  };
}

export type WorldEngine = ReturnType<typeof useWorldEngine>;
