import React from 'react';
import Controls from './components/Controls';
import WorldViewer from './components/WorldViewer';
import Map2D from './components/Map2D';
import MiniMap from './components/MiniMap';
import Inspector from './components/Inspector';
import { BiomeLegend } from './components/Legend';
import EditToolbar from './components/EditToolbar';
import { Menu, X } from 'lucide-react';
import ConfirmDialog from './components/ConfirmDialog';
import { useWorldEngine } from './hooks/useWorldEngine';

const App: React.FC = () => {
  const {
    params, setParams, world, viewMode, setViewMode, mapStyleId, setMapStyleId, displayMode, setDisplayMode, inspectMode, inspectorCollapsed, setInspectorCollapsed, inspectedCellId, rulerActive, markerMode, selectedMarkerId, setSelectedMarkerId, isGenerating, genProgress, logs, lore, isLoreLoading, showGrid, setShowGrid, smoothGlobe, setSmoothGlobe, showRivers, setShowRivers, showRoutes, setShowRoutes, showHillshade, setShowHillshade, showContours, setShowContours, showCurrents, setShowCurrents, labelVisibility, setLabelVisibility, sidebarOpen, setSidebarOpen, dymaxionSettings, setDymaxionSettings, apiKey, setApiKey, editMode, setEditMode, paintStyle, setPaintStyle, brushSize, setBrushSize, paintStrength, setPaintStrength, paintFaction, setPaintFaction, paintBiome, setPaintBiome, sampleHeight, adaptiveBiomes, setAdaptiveBiomes, undoStack, requestGenerate, pendingGenerate, confirmGenerate, cancelGenerate, handleLoadWorld, handleCancel, handleUpdateCivs, handleUpdateProvinces, toggleInspectEnabled, toggleRuler, toggleMarkerMode, handleInspect, updateMarker, deleteMarker, rulerArc, rulerDistanceKm, handleGenerateLore, factionColors, cultureColors, religionColors, handlePaint, handleUndo, handleEditWorldData, handleEditFaction, handleMergeFactions, handleRenameProvince, handleRenameTown, handleRelocateCapital,
  } = useWorldEngine();

  return (
    <div className="flex flex-col md:flex-row w-full h-full bg-black overflow-hidden font-sans text-ink">
      {/* Sidebar / Bottom Drawer */}
      <div className={`fixed inset-x-0 bottom-0 z-sheet md:relative md:inset-auto md:w-80 md:h-full
 bg-surface-sunken border-t md:border-t-0 md:border-r border-edge-subtle transition-transform duration-300
 ${sidebarOpen ? 'translate-y-0 md:translate-x-0' : 'translate-y-full md:-translate-x-full'}
 max-h-[85vh] md:max-h-full flex flex-col shadow-2xl`}>
        <Controls 
          params={params} setParams={setParams}
          onGenerate={requestGenerate} 
          onLoadWorld={handleLoadWorld}
          onCancel={handleCancel}
          onUpdateCivs={handleUpdateCivs} onUpdateProvinces={handleUpdateProvinces}
          viewMode={viewMode} setViewMode={setViewMode}
          mapStyleId={mapStyleId} setMapStyleId={setMapStyleId}
          displayMode={displayMode} setDisplayMode={setDisplayMode}
          loading={isGenerating} logs={logs} genProgress={genProgress}
          lore={lore} generatingLore={isLoreLoading} onGenerateLore={handleGenerateLore}
          worldData={world} 
          showGrid={showGrid} setShowGrid={setShowGrid}
          smoothGlobe={smoothGlobe} setSmoothGlobe={setSmoothGlobe}
          showRivers={showRivers} setShowRivers={setShowRivers}
          showRoutes={showRoutes} setShowRoutes={setShowRoutes}
          showHillshade={showHillshade} setShowHillshade={setShowHillshade}
          showContours={showContours} setShowContours={setShowContours}
          showCurrents={showCurrents} setShowCurrents={setShowCurrents}
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
          aria-label="Close settings panel"
          className="md:hidden absolute -top-12 right-4 bg-surface text-ink-strong p-2 border border-edge shadow-lg"
        >
          <X size={20} />
        </button>
      </div>

      {/* Floating menu button - Top Left to avoid overlapping with bottom overlays */}
      {!sidebarOpen && (
        <button 
          onClick={() => { setSidebarOpen(true); }}
          aria-label="Open settings panel"
          aria-expanded={sidebarOpen}
          className="fixed top-4 left-4 z-modal bg-brand-strong text-ink-strong p-3 shadow-2xl hover:bg-brand md:hidden border border-white/10"
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
            smoothGlobe={smoothGlobe}
            showRivers={showRivers}
            showRoutes={showRoutes}
            showHillshade={showHillshade}
            showContours={showContours}
            showCurrents={showCurrents}
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
            showCurrents={showCurrents}
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

        {/* Overlay HUD elements */}
        {displayMode === 'globe' && (
          <div className="absolute top-4 left-24 bg-black/50 backdrop-blur-md p-3 border border-white/10 text-left pointer-events-none z-overlay hidden md:block">
           <h3 className="text-ink-strong text-xs font-bold">{world ? `Seed: ${params.seed}` : 'No World'}</h3>
           <p className="text-ink-muted text-[10px]">{world ? `${world.cells.length.toLocaleString()} Cells` : ''}</p>
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
        {rulerActive && (
          <div className="absolute bottom-6 left-1/2 -translate-x-1/2 pointer-events-none z-overlay">
            <div className="bg-black/80 backdrop-blur text-ink-strong text-xs font-bold px-3 py-1.5 border border-white/20 shadow-xl">
              {rulerDistanceKm !== null
                ? `${Math.round(rulerDistanceKm).toLocaleString()} km`
                : 'Click two points'}
            </div>
          </div>
        )}
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
      <ConfirmDialog
        open={pendingGenerate !== null}
        title="Discard your edits?"
        body={`Generating replaces the world. ${undoStack.length} painted ${undoStack.length === 1 ? 'change' : 'changes'} will be lost, and this cannot be undone.`}
        confirmLabel="Discard & Generate"
        destructive
        onConfirm={confirmGenerate}
        onCancel={cancelGenerate}
      />
    </div>
  );
};

export default App;
