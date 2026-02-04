import React, { useState, useEffect, useCallback } from 'react';
import Controls from './components/Controls';
import WorldViewer from './components/WorldViewer';
import MiniMap from './components/MiniMap';
import { BiomeLegend } from './components/Legend';
import { WorldData, WorldParams, ViewMode, LoreData } from './types';
import { generateWorld, recalculateCivs, recalculateProvinces } from './utils/worldGen';
import { generateWorldLore } from './services/gemini';
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
  civSeed: 'terraforge_civs',
  borderRoughness: 0.2, 
  civSizeVariance: 0.5,
  waterCrossingCost: 0.8, 
  capitalSpacing: 0.5,
  provinceSize: 0.5,
  loreLevel: 1,
  seed: 'terraforge',
};

const App: React.FC = () => {
  const [params, setParams] = useState<WorldParams>(DEFAULT_PARAMS);
  const [world, setWorld] = useState<WorldData | null>(null);
  const [viewMode, setViewMode] = useState<ViewMode>('biome');
  const [isGenerating, setIsGenerating] = useState(false);
  const [progress, setProgress] = useState<{msg: string, percent: number}>({msg: '', percent: 0});
  const [lore, setLore] = useState<LoreData | null>(null);
  const [isLoreLoading, setIsLoreLoading] = useState(false);
  const [showGrid, setShowGrid] = useState(false);
  const [sidebarOpen, setSidebarOpen] = useState(true);

  const handleGenerate = useCallback(async (overrideParams?: WorldParams) => {
    setIsGenerating(true);
    setLore(null);
    setProgress({msg: 'Starting...', percent: 0});
    const p = overrideParams || params;
    await new Promise(r => setTimeout(r, 50));
    try {
        const newWorld = await generateWorld(p, (msg, pct) => setProgress({msg, percent: pct}));
        setWorld(newWorld);
    } catch (e) { console.error("Generation failed", e); }
    finally { setIsGenerating(false); setProgress({msg: '', percent: 100}); }
  }, [params]);

  const handleUpdateCivs = useCallback(async () => {
      if (!world) return;
      setIsGenerating(true); setProgress({msg: 'Updating Factions...', percent: 0});
      await new Promise(r => setTimeout(r, 50));
      try {
          const updatedWorld = recalculateCivs(world, params, (msg, pct) => setProgress({msg, percent: pct}));
          setWorld({ ...updatedWorld });
          if (viewMode !== 'political') setViewMode('political');
      } catch(e) { console.error("Civ update failed", e); }
      finally { setIsGenerating(false); setProgress({msg: '', percent: 100}); }
  }, [world, params, viewMode]);

  const handleUpdateProvinces = useCallback(async () => {
    if (!world || !world.civData) return;
    setIsGenerating(true); setProgress({msg: 'Reshuffling Provinces...', percent: 0});
    await new Promise(r => setTimeout(r, 50));
    try {
        const updatedWorld = recalculateProvinces(world, params);
        setWorld({ ...updatedWorld });
        if (viewMode !== 'political') setViewMode('political');
    } catch(e) { console.error("Province update failed", e); }
    finally { setIsGenerating(false); setProgress({msg: '', percent: 100}); }
  }, [world, params, viewMode]);

  useEffect(() => { handleGenerate(); }, []);

  const handleGenerateLore = async () => {
    if (!world) return;
    setIsLoreLoading(true);
    try {
      const newLore = await generateWorldLore(world);
      setLore(newLore);
      setWorld({ ...world });
    } catch (e) { console.error("Lore gen failed", e); }
    finally { setIsLoreLoading(false); }
  };

  return (
    <div className="flex flex-col md:flex-row w-full h-full bg-black overflow-hidden font-sans text-gray-200">
      {/* Sidebar / Bottom Drawer */}
      <div className={`
        fixed inset-x-0 bottom-0 z-30 md:relative md:inset-auto md:w-80 md:h-full
        bg-gray-950 border-t md:border-t-0 md:border-r border-gray-800 transition-transform duration-300
        ${sidebarOpen ? 'translate-y-0 md:translate-x-0' : 'translate-y-full md:-translate-x-full'}
        max-h-[85vh] md:max-h-full flex flex-col shadow-2xl
      `}>
        <Controls 
          params={params} setParams={setParams}
          onGenerate={handleGenerate} onUpdateCivs={handleUpdateCivs} onUpdateProvinces={handleUpdateProvinces}
          viewMode={viewMode} setViewMode={setViewMode}
          loading={isGenerating} progress={progress}
          lore={lore} generatingLore={isLoreLoading} onGenerateLore={handleGenerateLore}
          worldData={world} showGrid={showGrid} setShowGrid={setShowGrid}
        />
        <button 
          onClick={() => setSidebarOpen(false)}
          className="md:hidden absolute -top-12 right-4 bg-gray-900 text-white p-2 rounded-full border border-gray-700 shadow-lg"
        >
          <X size={20} />
        </button>
      </div>

      {/* Floating menu button - Top Left to avoid overlapping with bottom overlays */}
      {!sidebarOpen && (
        <button 
          onClick={() => setSidebarOpen(true)}
          className="fixed top-4 left-4 z-40 bg-blue-600 text-white p-3 rounded-full shadow-2xl hover:bg-blue-500 md:hidden border border-white/10"
        >
          <Menu size={24} />
        </button>
      )}

      {/* Main Content Area */}
      <main className="flex-1 relative h-full overflow-hidden">
        <WorldViewer world={world} viewMode={viewMode} showGrid={showGrid} />
        
        {/* Overlay HUD elements */}
        <div className="absolute top-4 left-24 bg-black/50 backdrop-blur-md p-3 rounded-lg border border-white/10 text-left pointer-events-none z-10 hidden md:block">
           <h3 className="text-white text-xs font-bold">{world ? `Seed: ${params.seed}` : 'No World'}</h3>
           <p className="text-gray-400 text-[10px]">{world ? `${world.cells.length.toLocaleString()} Cells` : ''}</p>
        </div>

        <BiomeLegend />
        <MiniMap world={world} viewMode={viewMode} />
      </main>
    </div>
  );
};

export default App;