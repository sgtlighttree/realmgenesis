import React from 'react';
import { Menu, X } from 'lucide-react';
import Controls from '../Controls';
import WorldViewer from '../WorldViewer';
import Map2D from '../Map2D';
import MiniMap from '../MiniMap';
import Inspector from '../Inspector';
import { BiomeLegend } from '../Legend';
import EditToolbar from '../EditToolbar';
import { PlaceholderGlobe } from './shellKit';
import { useWorldEngine } from '../../hooks/useWorldEngine';

/**
 * ShellApp — the F1 redesign entry (?shell=1). Consumes the SAME useWorldEngine
 * hook as the classic App, so it shares one brain. v1 is a clean reframe of the
 * real render with the ?globe=0 escape hatch (swap the Three.js globe for a
 * placeholder during UI iteration). The docked bucket-model layout
 * (WideShell/NarrowShell) is wired in a later pass; ?shell=stub still serves the
 * DesignShell prototype for that reference.
 */
const globeDisabled = (): boolean =>
  new URLSearchParams(window.location.search).get('globe') === '0';

const ShellApp: React.FC = () => {
  const {
    params, setParams, world, viewMode, setViewMode, displayMode, setDisplayMode,
    inspectMode, inspectorCollapsed, setInspectorCollapsed, inspectedCellId,
    rulerActive, markerMode, selectedMarkerId, setSelectedMarkerId, isGenerating,
    genProgress, logs, lore, isLoreLoading, showGrid, setShowGrid, showRivers,
    setShowRivers, showRoutes, setShowRoutes, showHillshade, setShowHillshade,
    showContours, setShowContours, labelVisibility, setLabelVisibility, sidebarOpen,
    setSidebarOpen, dymaxionSettings, setDymaxionSettings, apiKey, setApiKey,
    editMode, setEditMode, paintStyle, setPaintStyle, brushSize, setBrushSize,
    paintStrength, setPaintStrength, paintFaction, setPaintFaction, paintBiome,
    setPaintBiome, sampleHeight, adaptiveBiomes, setAdaptiveBiomes, undoStack,
    handleGenerate, handleLoadWorld, handleCancel, handleUpdateCivs,
    handleUpdateProvinces, toggleInspectEnabled, toggleRuler, toggleMarkerMode,
    handleInspect, updateMarker, deleteMarker, rulerArc, rulerDistanceKm,
    handleGenerateLore, factionColors, cultureColors, religionColors, handlePaint,
    handleUndo, handleEditWorldData, handleEditFaction, handleMergeFactions,
    handleRenameProvince, handleRenameTown, handleRelocateCapital,
  } = useWorldEngine();

  const noGlobe = globeDisabled();

  return (
    <div className="flex flex-col md:flex-row w-full h-full bg-black overflow-hidden font-sans text-gray-200">
      {/* Make rail (shell-styled) / bottom drawer on mobile */}
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
          showRoutes={showRoutes} setShowRoutes={setShowRoutes}
          showHillshade={showHillshade} setShowHillshade={setShowHillshade}
          showContours={showContours} setShowContours={setShowContours}
          labelVisibility={labelVisibility} setLabelVisibility={setLabelVisibility}
          dymaxionSettings={dymaxionSettings}
          onDymaxionChange={setDymaxionSettings}
          apiKey={apiKey}
          onApiKeyChange={setApiKey}
          onInspect={handleInspect}
          onEditFaction={handleEditFaction}
          onMergeFactions={handleMergeFactions}
        />
        <button
          onClick={() => { setSidebarOpen(false); }}
          className="md:hidden absolute -top-12 right-4 bg-gray-900 text-white p-2 border border-gray-700 shadow-lg"
        >
          <X size={20} />
        </button>
      </div>

      {!sidebarOpen && (
        <button
          onClick={() => { setSidebarOpen(true); }}
          className="fixed top-4 left-4 z-40 bg-blue-600 text-white p-3 shadow-2xl hover:bg-blue-500 md:hidden border border-white/10"
        >
          <Menu size={24} />
        </button>
      )}

      {/* Canvas */}
      <main className="flex-1 relative h-full overflow-hidden">
        {noGlobe ? (
          <PlaceholderGlobe />
        ) : displayMode === 'globe' ? (
          <WorldViewer
            world={world}
            viewMode={viewMode}
            showGrid={showGrid}
            showRivers={showRivers}
            showRoutes={showRoutes}
            showHillshade={showHillshade}
            showContours={showContours}
            labelVisibility={labelVisibility}
            inspectMode={inspectMode}
            onInspect={handleInspect}
            selectedCellId={inspectedCellId}
            dymaxionSettings={dymaxionSettings}
            onDymaxionChange={setDymaxionSettings}
            editMode={editMode}
            onPaint={handlePaint}
            factionColors={factionColors}
            cultureColors={cultureColors}
            religionColors={religionColors}
            brushSize={brushSize}
            rulerArc={rulerArc}
          />
        ) : (
          <Map2D
            world={world}
            viewMode={viewMode}
            inspectMode={inspectMode}
            onInspect={handleInspect}
            highlightCellId={inspectedCellId}
            projectionType={displayMode === 'dymaxion' ? 'dymaxion' : 'mercator'}
            dymaxionSettings={dymaxionSettings}
            showGrid={showGrid}
            showRivers={showRivers}
            showRoutes={showRoutes}
            showHillshade={showHillshade}
            showContours={showContours}
            labelVisibility={labelVisibility}
            editMode={editMode}
            onPaint={handlePaint}
            factionColors={factionColors}
            cultureColors={cultureColors}
            religionColors={religionColors}
            brushSize={brushSize}
            rulerArc={rulerArc}
          />
        )}

        {!noGlobe && displayMode === 'globe' && (
          <div className="absolute top-4 left-24 bg-black/50 backdrop-blur-md p-3 border border-white/10 text-left pointer-events-none z-10 hidden md:block">
            <h3 className="text-white text-xs font-bold">{world ? `Seed: ${params.seed}` : 'No World'}</h3>
            <p className="text-gray-400 text-[10px]">{world ? `${world.cells.length.toLocaleString()} Cells` : ''}</p>
          </div>
        )}

        {!noGlobe && <BiomeLegend />}
        {!noGlobe && displayMode === 'globe' && <MiniMap world={world} viewMode={viewMode} />}
        {!noGlobe && (
          <Inspector
            world={world}
            cellId={inspectedCellId}
            inspectMode={inspectMode}
            collapsed={inspectorCollapsed}
            onToggleEnabled={toggleInspectEnabled}
            onToggleCollapsed={() => { setInspectorCollapsed(v => !v); }}
            editMode={editMode}
            onEditWorldData={handleEditWorldData}
            onRenameProvince={handleRenameProvince}
            onRenameTown={handleRenameTown}
            onRelocateCapital={handleRelocateCapital}
            rulerActive={rulerActive}
            onToggleRuler={toggleRuler}
            markerMode={markerMode}
            onToggleMarkerMode={toggleMarkerMode}
            selectedMarkerId={selectedMarkerId}
            onSelectMarker={setSelectedMarkerId}
            onUpdateMarker={updateMarker}
            onDeleteMarker={deleteMarker}
          />
        )}
        {!noGlobe && rulerActive && (
          <div className="absolute bottom-6 left-1/2 -translate-x-1/2 pointer-events-none z-10">
            <div className="bg-black/80 backdrop-blur text-white text-xs font-bold px-3 py-1.5 border border-white/20 shadow-xl">
              {rulerDistanceKm !== null
                ? `${Math.round(rulerDistanceKm).toLocaleString()} km`
                : 'Click two points'}
            </div>
          </div>
        )}
        {!noGlobe && world && (
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

        {/* Placeholder-globe iteration banner */}
        {noGlobe && (
          <div className="absolute top-4 left-1/2 -translate-x-1/2 z-10 text-[11px] font-mono text-gray-500 bg-gray-950/80 border border-gray-800 rounded px-3 py-1.5 pointer-events-none">
            globe=0 · placeholder mode (real globe disabled for UI iteration)
          </div>
        )}
      </main>
    </div>
  );
};

export default ShellApp;
