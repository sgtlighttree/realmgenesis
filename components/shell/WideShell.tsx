import React from 'react';
import { Globe, Pencil } from 'lucide-react';
import { Panel, PANEL, FOCUS, ShellProps } from './shellKit';

/**
 * WideShell — the "A · Tidy" desktop layout.
 * Make = solid left rail. View = top strip. Read = docked stack, right.
 * Do = contextual bar, summoned only while editing. Canvas fills the middle.
 * Positions pre-built panels; never touches their props.
 *
 * SPACING CONTRACT (4pt scale, and the reason the rail is flush):
 * the Make rail was previously `p-3` around a `Panel` around Controls' own
 * `p-4` — three nested paddings costing ~58px of the 288px rail versus
 * classic's single padding, which is what produced the cramped column and its
 * horizontal scrollbar. The rail is now flush to the viewport edge like the
 * classic sidebar: one border, one padding, owned by Controls. Floating
 * surfaces over the canvas keep a uniform 8px inset (`top-2`/`right-2`) and 8px
 * between siblings; panel interiors use 12px.
 */
const WideShell: React.FC<ShellProps> = ({
  make, view, read, doTools, canvas, editing, onSetEditing,
}) => (
  <div className="absolute inset-0 flex bg-black text-gray-200 font-sans">

    {/* Make — flush left rail. No outer padding and no Panel wrapper: Controls
        supplies its own chrome and scroll region, so anything here is nesting. */}
    <aside className="w-72 shrink-0 h-full flex flex-col border-r border-gray-800 bg-gray-950">
      <div className="flex items-center gap-2 px-3 py-2.5 border-b border-gray-800">
        <Globe size={16} className="text-blue-400" />
        <span className="font-semibold text-[13px] tracking-tight">RealmGenesis 3D</span>
      </div>
      <div className="flex-1 min-h-0">{make}</div>
    </aside>

    {/* Canvas + docked panels */}
    <div className="relative flex-1 min-w-0">
      <div className="absolute inset-0">{canvas}</div>

      {/* View — top strip. Right edge clears the Read rail (16rem + 2×8px). */}
      <div className="absolute top-2 left-2 right-[16.5rem] z-20">
        <div className={`${PANEL} flex items-center gap-2 px-2 py-1.5`}>
          {view}
          <button
            onClick={() => onSetEditing(!editing)}
            aria-pressed={editing}
            className={`ml-auto shrink-0 inline-flex items-center gap-1.5 px-2.5 py-1.5 text-[11px] font-medium border transition-colors ${FOCUS}
              ${editing
                ? 'bg-amber-500 text-black border-amber-400'
                : 'bg-gray-800 text-gray-300 border-gray-700 hover:bg-gray-700 hover:text-white'}`}
          >
            <Pencil size={12} /> Edit
          </button>
        </div>
      </div>

      {/* Read — docked stack, right edge */}
      <div className="absolute top-2 right-2 bottom-2 w-64 flex flex-col gap-2 overflow-y-auto z-20">
        {read.map(card => (
          <Panel key={card.key} title={card.title} className="shrink-0">
            {card.node}
          </Panel>
        ))}
      </div>

      {/* Do — contextual bottom bar. Max width clears the Read rail. */}
      {editing && (
        <div className="absolute bottom-3 left-1/2 -translate-x-1/2 max-w-[calc(100%_-_17.5rem)] z-20 rg-rise">
          <div className={`${PANEL} px-2 py-1.5 shadow-xl`}>{doTools}</div>
        </div>
      )}
    </div>
  </div>
);

export default WideShell;
