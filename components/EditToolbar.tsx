import React, { useEffect } from 'react';
import { Undo2, Pencil, ArrowUp, ArrowDown, Minus, Waves, Leaf, Flag, Globe, Eraser } from 'lucide-react';
import { EditMode, PaintStyle, BiomeType, WorldData, POLITICAL_ERASER_ID } from '../types';
import { BIOME_COLORS } from '../utils/colors';

interface EditToolbarProps {
  editMode: EditMode;
  setEditMode: (m: EditMode) => void;
  paintStyle: PaintStyle;
  setPaintStyle: (s: PaintStyle) => void;
  brushSize: number;
  setBrushSize: (n: number) => void;
  paintStrength: number;
  setPaintStrength: (n: number) => void;
  adaptiveBiomes: boolean;
  setAdaptiveBiomes: (v: boolean) => void;
  paintFaction: number;
  setPaintFaction: (id: number) => void;
  paintBiome: BiomeType;
  setPaintBiome: (b: BiomeType) => void;
  sampleHeight: number | null;
  undoCount: number;
  onUndo: () => void;
  world: WorldData;
}

// Lakes are hydrology-derived (surfaced from depressions), not hand-painted.
const BIOME_LIST = Object.values(BiomeType).filter(
  b => b !== BiomeType.LAKE && b !== BiomeType.SALT_LAKE,
);

const ModeBtn: React.FC<{
  active: boolean;
  onClick: () => void;
  icon: React.ReactNode;
  label: string;
  title: string;
}> = ({ active, onClick, icon, label, title }) => (
  <button
    onClick={onClick}
    title={title}
    className={`flex items-center gap-1 px-2 py-1.5 text-[10px] font-medium border transition-colors whitespace-nowrap
      ${active
        ? 'bg-blue-600 text-white border-blue-500'
        : 'bg-gray-800 text-gray-400 border-gray-700 hover:bg-gray-700 hover:text-white'
      }`}
  >
    {icon}
    {label}
  </button>
);

const BrushSlider: React.FC<{ brushSize: number; setBrushSize: (n: number) => void }> = ({ brushSize, setBrushSize }) => (
  <div className="flex items-center gap-1.5">
    <span className="text-[10px] text-gray-400">Brush</span>
    <input type="range" min="0" max="5" step="1"
      value={brushSize} onChange={e => setBrushSize(parseInt(e.target.value))}
      className="w-16 h-1 bg-gray-700 appearance-none cursor-pointer accent-blue-400" />
    <span className="text-[10px] text-gray-300 w-3 text-center">{brushSize}</span>
  </div>
);

const StrengthSlider: React.FC<{ paintStrength: number; setPaintStrength: (n: number) => void }> = ({ paintStrength, setPaintStrength }) => (
  <div className="flex items-center gap-1.5">
    <span className="text-[10px] text-gray-400">Strength</span>
    <input type="range" min="0.1" max="1.0" step="0.05"
      value={paintStrength} onChange={e => setPaintStrength(parseFloat(e.target.value))}
      className="w-16 h-1 bg-gray-700 appearance-none cursor-pointer accent-blue-400" />
  </div>
);

const Sep = () => <div className="w-px h-4 bg-gray-700 flex-shrink-0" />;

const getContrastText = (hex: string) => {
  const value = hex.replace('#', '');
  if (value.length !== 6) return '#ffffff';
  const r = parseInt(value.slice(0, 2), 16);
  const g = parseInt(value.slice(2, 4), 16);
  const b = parseInt(value.slice(4, 6), 16);
  const luminance = (0.299 * r + 0.587 * g + 0.114 * b) / 255;
  return luminance > 0.58 ? '#111827' : '#ffffff';
};

const EditToolbar: React.FC<EditToolbarProps> = ({
  editMode, setEditMode,
  paintStyle, setPaintStyle,
  brushSize, setBrushSize,
  paintStrength, setPaintStrength,
  adaptiveBiomes, setAdaptiveBiomes,
  paintFaction, setPaintFaction,
  paintBiome, setPaintBiome,
  sampleHeight,
  undoCount, onUndo,
  world,
}) => {
  const isTerrainMode = editMode === 'terrain-raise' || editMode === 'terrain-lower'
    || editMode === 'terrain-flatten' || editMode === 'terrain-smooth';
  const factions = world.civData?.factions ?? [];

  useEffect(() => {
    if (factions.length > 0 && paintFaction !== POLITICAL_ERASER_ID && !factions.find(f => f.id === paintFaction)) {
      setPaintFaction(factions[0].id);
    }
  }, [factions, paintFaction, setPaintFaction]);

  return (
    <div className="absolute bottom-20 left-1/2 -translate-x-1/2 z-20 flex flex-col items-center gap-1 pointer-events-auto select-none">

      {/* Sub-controls */}
      {editMode !== 'off' && (
        <div className="flex items-center gap-2 bg-black/85 backdrop-blur border border-white/10 px-3 py-1.5 shadow-xl flex-wrap justify-center max-w-[560px]">

          {/* Terrain sub-controls */}
          {isTerrainMode && (
            <>
              {(editMode === 'terrain-raise' || editMode === 'terrain-lower') && (
                <div className="flex gap-1">
                  <button onClick={() => setPaintStyle('adaptive')} title="Gentle sculpt — blends edges with neighbors"
                    className={`text-[10px] px-2 py-1 border transition-colors ${paintStyle === 'adaptive' ? 'bg-teal-700 text-white border-teal-500' : 'bg-gray-800 text-gray-400 border-gray-700 hover:bg-gray-700 hover:text-white'}`}>
                    Adaptive
                  </button>
                  <button onClick={() => setPaintStyle('freeform')} title="Strong direct paint — no edge blending"
                    className={`text-[10px] px-2 py-1 border transition-colors ${paintStyle === 'freeform' ? 'bg-orange-700 text-white border-orange-500' : 'bg-gray-800 text-gray-400 border-gray-700 hover:bg-gray-700 hover:text-white'}`}>
                    Freeform
                  </button>
                </div>
              )}
              {editMode === 'terrain-flatten' && (
                sampleHeight !== null
                  ? <span className="text-[10px] text-yellow-400">Target: {(sampleHeight * 100).toFixed(0)}%</span>
                  : <span className="text-[10px] text-gray-400">Right-click cell to set target height</span>
              )}
              <Sep />
              <StrengthSlider paintStrength={paintStrength} setPaintStrength={setPaintStrength} />
              <Sep />
              <BrushSlider brushSize={brushSize} setBrushSize={setBrushSize} />
              <Sep />
              <label className="flex items-center gap-1 text-[10px] text-gray-400 cursor-pointer">
                <input type="checkbox" checked={adaptiveBiomes} onChange={e => setAdaptiveBiomes(e.target.checked)}
                  className="w-3 h-3 accent-teal-500" />
                Adaptive Biomes
              </label>
            </>
          )}

          {/* Biome sub-controls */}
          {editMode === 'biome' && (
            <>
              <span className="text-[10px] text-gray-300 font-medium">{paintBiome}</span>
              <Sep />
              <div className="flex gap-1 flex-wrap max-w-[220px]">
                {BIOME_LIST.map(b => (
                  <button key={b} onClick={() => setPaintBiome(b)} title={b}
                    className={`w-5 h-5 border-2 transition-all ${paintBiome === b ? 'scale-125 border-white' : 'border-transparent hover:border-gray-400'}`}
                    style={{ backgroundColor: BIOME_COLORS[b] }} />
                ))}
              </div>
              <Sep />
              <BrushSlider brushSize={brushSize} setBrushSize={setBrushSize} />
            </>
          )}

          {/* Political sub-controls */}
          {editMode === 'political' && (
            <>
              <BrushSlider brushSize={brushSize} setBrushSize={setBrushSize} />
              {factions.length > 0 && (
                <>
                  <Sep />
                  <div className="flex gap-1.5 flex-wrap">
                    <button
                      onClick={() => setPaintFaction(POLITICAL_ERASER_ID)}
                      title="Eraser: mark cells as unclaimed"
                      className={`w-6 h-6 border-2 transition-all flex items-center justify-center bg-gray-950 text-gray-200 ${paintFaction === POLITICAL_ERASER_ID ? 'scale-125 border-white' : 'border-gray-700 hover:border-gray-400'}`}
                    >
                      <Eraser size={13} />
                    </button>
                    {factions.map((f, idx) => {
                      const textColor = getContrastText(f.color);
                      return (
                        <button key={f.id} onClick={() => setPaintFaction(f.id)} title={`${idx + 1}. ${f.name}`}
                          className={`w-6 h-6 border-2 transition-all flex items-center justify-center text-[10px] font-bold leading-none ${paintFaction === f.id ? 'scale-125 border-white' : 'border-transparent hover:border-gray-400'}`}
                          style={{ backgroundColor: f.color, color: textColor }}>
                          <span style={{ textShadow: textColor === '#ffffff' ? '0 1px 2px #000' : '0 1px 2px #fff' }}>
                            {idx + 1}
                          </span>
                        </button>
                      );
                    })}
                  </div>
                </>
              )}
            </>
          )}

          {/* World-edit hint */}
          {editMode === 'world-edit' && (
            <span className="text-[10px] text-gray-400">Click a cell to edit its data in the Inspector</span>
          )}

          <Sep />
          <button onClick={onUndo} disabled={undoCount === 0}
            title={`Undo last stroke (${undoCount}) — Ctrl+Z`}
            className="flex items-center gap-1 text-[10px] text-gray-400 hover:text-white disabled:opacity-30 disabled:cursor-not-allowed transition-colors">
            <Undo2 size={12} />
            {undoCount > 0 && <span>{undoCount}</span>}
          </button>
        </div>
      )}

      {/* Space hint */}
      {editMode !== 'off' && (
        <div className="text-[9px] text-gray-600 text-center">
          Hold <kbd className="bg-gray-800 border border-gray-700 px-1 py-0.5 text-gray-400">Space</kbd> + drag to orbit / pan
        </div>
      )}

      {/* Mode button row — grouped by category */}
      <div className="flex items-center gap-0.5 bg-black/85 backdrop-blur border border-white/10 px-2 py-1.5 shadow-xl">
        {/* Off + Undo */}
        <ModeBtn active={editMode === 'off'} onClick={() => setEditMode('off')} icon={<Globe size={11} />} label="Off" title="Disable editing" />
        <Sep />
        {/* Terrain group */}
        <span className="text-[9px] text-gray-600 px-1 uppercase tracking-wide hidden sm:block">Terrain</span>
        <ModeBtn active={editMode === 'terrain-raise'} onClick={() => setEditMode('terrain-raise')} icon={<ArrowUp size={11} />} label="Raise" title="Raise terrain" />
        <ModeBtn active={editMode === 'terrain-lower'} onClick={() => setEditMode('terrain-lower')} icon={<ArrowDown size={11} />} label="Lower" title="Lower terrain" />
        <ModeBtn active={editMode === 'terrain-flatten'} onClick={() => setEditMode('terrain-flatten')} icon={<Minus size={11} />} label="Flatten" title="Right-click to sample height, drag to flatten" />
        <ModeBtn active={editMode === 'terrain-smooth'} onClick={() => setEditMode('terrain-smooth')} icon={<Waves size={11} />} label="Smooth" title="Smooth terrain — average heights with neighbors" />
        <Sep />
        {/* Other tools */}
        <ModeBtn active={editMode === 'biome'} onClick={() => setEditMode('biome')} icon={<Leaf size={11} />} label="Biome" title="Force-paint biome type" />
        <ModeBtn active={editMode === 'political'} onClick={() => setEditMode('political')} icon={<Flag size={11} />} label="Political" title="Repaint faction borders" />
        <ModeBtn active={editMode === 'world-edit'} onClick={() => setEditMode('world-edit')} icon={<Pencil size={11} />} label="Edit" title="Edit names & demographics" />
      </div>
    </div>
  );
};

export default EditToolbar;
