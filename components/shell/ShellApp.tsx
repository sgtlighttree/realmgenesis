import React, { useState, useEffect, useCallback } from 'react';

import Controls from '../Controls';
import WorldViewer from '../WorldViewer';
import Map2D from '../Map2D';
import Inspector from '../Inspector';
import EditToolbar from '../EditToolbar';
import { BiomeLegendList } from '../Legend';
import { MiniMapCanvas } from '../MiniMap';
import DymaxionNetPreview from '../DymaxionNetPreview';
import { ViewStrip } from '../ViewControls';
import WideShell from './WideShell';
import NarrowShell from './NarrowShell';
import { PlaceholderGlobe, ReadCard } from './shellKit';
import ConfirmDialog from '../ConfirmDialog';
import { useWorldEngine } from '../../hooks/useWorldEngine';

/**
 * ShellApp — the F1 redesign entry (?shell=1). WIRING ONLY: it consumes the
 * same `useWorldEngine` hook as the classic App, composes each panel WITH its
 * props into a finished element, and hands those to whichever shell is active
 * as named slots. The shell only positions pre-built elements and never sees
 * their props, which is what keeps the no-Context / prop-drilling invariant
 * intact across two layout trees.
 *
 * `?globe=0` swaps the Three.js globe for a placeholder (fast UI iteration
 * without the WebGL cost). It skips RENDERING, not generation — a world still
 * generates on mount.
 */
const globeDisabled = (): boolean =>
  new URLSearchParams(window.location.search).get('globe') === '0';

/** Which fold to render. Matches Tailwind's `md` breakpoint. */
const useIsNarrow = (): boolean => {
  const [narrow, setNarrow] = useState(() => window.matchMedia('(max-width: 767px)').matches);
  useEffect(() => {
    const mq = window.matchMedia('(max-width: 767px)');
    const onChange = (e: MediaQueryListEvent) => { setNarrow(e.matches); };
    mq.addEventListener('change', onChange);
    return () => { mq.removeEventListener('change', onChange); };
  }, []);
  return narrow;
};

const ShellApp: React.FC = () => {
  const engine = useWorldEngine();
  const {
    params, setParams, world, viewMode, setViewMode, mapStyleId, setMapStyleId, displayMode, setDisplayMode,
    inspectMode, inspectedCellId, rulerActive, markerMode, selectedMarkerId,
    setSelectedMarkerId, isGenerating, genProgress, logs, lore, isLoreLoading,
    showGrid, setShowGrid, smoothGlobe, setSmoothGlobe, showRivers, setShowRivers, showRoutes, setShowRoutes,
    showHillshade, setShowHillshade, showContours, setShowContours,
    showCurrents, setShowCurrents,
    showCellEdges, setShowCellEdges,
    labelVisibility, setLabelVisibility, dymaxionSettings, setDymaxionSettings,
    apiKey, setApiKey, editMode, setEditMode, paintStyle, setPaintStyle,
    brushSize, setBrushSize, paintStrength, setPaintStrength, paintFaction,
    setPaintFaction, paintBiome, setPaintBiome, sampleHeight, adaptiveBiomes,
    setAdaptiveBiomes, undoStack, requestGenerate, pendingGenerate, confirmGenerate,
    cancelGenerate, handleLoadWorld, handleCancel,
    handleUpdateCivs, handleUpdateProvinces, toggleInspectEnabled, toggleRuler,
    toggleMarkerMode, handleInspect, updateMarker, deleteMarker, rulerArc,
    rulerDistanceKm, handleGenerateLore, factionColors, cultureColors,
    religionColors, handlePaint, handleUndo, handleEditWorldData,
    handleEditFaction, handleMergeFactions, handleRenameProvince,
    handleRenameTown, handleRelocateCapital,
  } = engine;

  const noGlobe = globeDisabled();
  const isNarrow = useIsNarrow();

  // Edit mode is the shell's contextual "Do" bucket. `editing` and the Esc
  // listener stay HERE, not in a shell: NarrowShell derives editing from its
  // open tab and WideShell from its Edit button, but both just call
  // onSetEditing, so the single source of truth (and the setEditMode('off')
  // side effect) lives in one place for both folds.
  const [editing, setEditing] = useState(false);
  const onSetEditing = useCallback((v: boolean) => {
    setEditing(v);
    if (!v) setEditMode('off');
  }, [setEditMode]);

  useEffect(() => {
    const onKey = (e: KeyboardEvent) => { if (e.key === 'Escape') onSetEditing(false); };
    window.addEventListener('keydown', onKey);
    return () => { window.removeEventListener('keydown', onKey); };
  }, [onSetEditing]);

  // Globe rotation is lifted out of WorldViewer so the wide fold can render the
  // control in the top strip. It is ephemeral presentation state, so it lives
  // here in the shell rather than in the engine hook — same rule as `editing`.
  const [paused, setPaused] = useState(false);
  const togglePause = useCallback(() => { setPaused(v => !v); }, []);

  /* ---------------- Slot composition ---------------- */

  const make = (
    <Controls
      className="w-full flex flex-col h-full overflow-hidden text-sm"
      showHeader={false}
      showViewControls={false}
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
      showCellEdges={showCellEdges} setShowCellEdges={setShowCellEdges}
      labelVisibility={labelVisibility} setLabelVisibility={setLabelVisibility}
      dymaxionSettings={dymaxionSettings}
      onDymaxionChange={setDymaxionSettings}
      apiKey={apiKey}
      onApiKeyChange={setApiKey}
      onInspect={handleInspect}
      onEditFaction={handleEditFaction}
      onMergeFactions={handleMergeFactions}
    />
  );

  const view = (
    <ViewStrip
      viewMode={viewMode} setViewMode={setViewMode}
      mapStyleId={mapStyleId} setMapStyleId={setMapStyleId}
      displayMode={displayMode} setDisplayMode={setDisplayMode}
      showGrid={showGrid} setShowGrid={setShowGrid}
      smoothGlobe={smoothGlobe} setSmoothGlobe={setSmoothGlobe}
      showRivers={showRivers} setShowRivers={setShowRivers}
      showRoutes={showRoutes} setShowRoutes={setShowRoutes}
      showHillshade={showHillshade} setShowHillshade={setShowHillshade}
      showContours={showContours} setShowContours={setShowContours}
      showCurrents={showCurrents} setShowCurrents={setShowCurrents}
      showCellEdges={showCellEdges} setShowCellEdges={setShowCellEdges}
      labelVisibility={labelVisibility} setLabelVisibility={setLabelVisibility}
      // Only the wide strip adopts the rotation control; the narrow fold's View
      // sheet is behind a tab, and its canvas is unshifted, so the viewer's own
      // overlay is both reachable and correctly placed there.
      rotation={!isNarrow && !noGlobe && displayMode === 'globe' ? {
        paused,
        onToggle: togglePause,
        disabled: dymaxionSettings.mode === 'overlay',
      } : undefined}
    />
  );

  // Inspector carries its own header of working tools (ruler, marker, eye), so
  // its card gets no Panel title — that header IS the header. Docked, it needs
  // neither the float anchoring nor its own box, and `onToggleCollapsed` is
  // omitted deliberately: the rail scrolls, so there is nothing to collapse
  // into and the chevron would render as a dead control.
  const read: ReadCard[] = [
    {
      key: 'inspector',
      title: '',
      node: (
        <Inspector
          className="w-full"
          cardClassName="text-ink-strong w-full"
          world={world}
          cellId={inspectedCellId}
          inspectMode={inspectMode}
          collapsed={false}
          onToggleEnabled={toggleInspectEnabled}
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
      ),
    },
    { key: 'biomes', title: 'Biomes', collapsible: true, node: <BiomeLegendList columns={2} /> },
  ];
  // MiniMapCanvas renders null without a world, which would leave an empty
  // titled box; and it must be UNMOUNTED (not hidden) to stop its redraw.
  if (world && !noGlobe) {
    read.push({
      key: 'minimap',
      title: '2D Projection',
      collapsible: true,
      // With the Dymaxion cage on, show the live ACTUAL Dymaxion net projection
      // here (updates as the cage moves) instead of the plain minimap.
      node: dymaxionSettings.showOverlay
        ? <DymaxionNetPreview world={world} viewMode={viewMode} settings={dymaxionSettings} />
        : <MiniMapCanvas world={world} viewMode={viewMode} />,
    });
  }

  const doTools = world ? (
    <EditToolbar
      className="flex flex-col items-center gap-1 select-none w-full"
      rowChrome=""
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
  ) : (
    <p className="text-[11px] text-ink-muted italic">Generate a world to start editing.</p>
  );

  const canvas = (
    <div className="absolute inset-0">
      {noGlobe ? (
        <PlaceholderGlobe />
      ) : displayMode === 'globe' ? (
        <WorldViewer
          world={world}
          viewMode={viewMode}
          showGrid={showGrid}
          smoothGlobe={smoothGlobe}
          showRivers={showRivers}
          showRoutes={showRoutes}
          showHillshade={showHillshade}
          mapStyleId={mapStyleId}
          showContours={showContours}
          showCurrents={showCurrents}
          showCellEdges={showCellEdges}
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
          paused={paused}
          onPausedChange={setPaused}
          // Wide renders this in the top strip instead. The canvas overlay was
          // anchored bottom-LEFT of the viewer, and WideShell shifts the canvas
          // left by the Read rail's width inside an overflow-hidden column — so
          // the control was painting at x≈36 in a column that clips below 288
          // and had been invisible since that shift landed.
          showPauseControl={isNarrow}
          overlayClassName="absolute bottom-3 left-3 z-overlay flex gap-2"
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
          mapStyleId={mapStyleId}
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

      {/* Seed caption. Second casualty of the same clipping as the pause
          control: everything here is anchored to the CANVAS, which WideShell
          shifts left by 16.5rem, so a left-anchored overlay lands outside the
          visible column. Narrow's canvas is unshifted and keeps plain `left-3`
          (offset to `left-16` only when the viewer's own pause overlay is
          beside it). Wide adds the shift back so the caption sits 12px inside
          the visible edge; the 16.5rem must track WideShell's canvas offset. */}
      {world && (
        <div className={`absolute bottom-3 pointer-events-none text-[10px] font-mono text-ink-muted ${
          isNarrow
            ? (!noGlobe && displayMode === 'globe' ? 'left-16' : 'left-3')
            : 'left-[17.25rem]'
        }`}>
          {params.seed} · {world.cells.length.toLocaleString()} cells
        </div>
      )}

      {rulerActive && (
        <div className="absolute bottom-6 left-1/2 -translate-x-1/2 pointer-events-none z-overlay">
          <div className="bg-black/80 backdrop-blur text-ink-strong text-xs font-bold px-3 py-1.5 border border-white/20 shadow-xl">
            {rulerDistanceKm !== null
              ? `${Math.round(rulerDistanceKm).toLocaleString()} km`
              : 'Click two points'}
          </div>
        </div>
      )}

      {noGlobe && (
        <div className="absolute top-3 left-1/2 -translate-x-1/2 z-overlay text-[11px] font-mono text-ink-muted bg-surface-sunken/80 border border-edge-subtle px-3 py-1.5 pointer-events-none">
          globe=0 · placeholder mode (real globe disabled for UI iteration)
        </div>
      )}
    </div>
  );

  const Shell = isNarrow ? NarrowShell : WideShell;

  return (
    <>
    <Shell
      make={make}
      view={view}
      read={read}
      doTools={doTools}
      canvas={canvas}
      editing={editing}
      onSetEditing={onSetEditing}
    />
    <ConfirmDialog
      open={pendingGenerate !== null}
      title="Discard your edits?"
      body={`Generating replaces the world. ${undoStack.length} painted ${undoStack.length === 1 ? 'change' : 'changes'} will be lost, and this cannot be undone.`}
      confirmLabel="Discard & Generate"
      destructive
      onConfirm={confirmGenerate}
      onCancel={cancelGenerate}
    />
    </>
  );
};

export default ShellApp;
