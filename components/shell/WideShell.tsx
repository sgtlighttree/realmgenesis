import React from 'react';
import { Globe, Pencil } from 'lucide-react';
import { Panel, PANEL, FOCUS, ShellProps } from './shellKit';

/**
 * WideShell — the "A · Tidy" desktop layout.
 * Make = solid left rail. View = top strip. Read = floating stack, right.
 * Do = contextual bar, summoned only while editing. Canvas fills the middle.
 * Positions pre-built panels; never touches their props.
 */
const WideShell: React.FC<ShellProps> = ({
  make, view, read, doTools, canvas, editing, onSetEditing,
}) => (
  <div className="absolute inset-0 flex bg-black text-gray-200 font-sans">

    {/* Make — solid left rail */}
    <aside className="w-72 shrink-0 h-full flex flex-col gap-3 p-3 border-r border-gray-800 bg-gray-950">
      <div className="flex items-center gap-2 px-1 py-1">
        <Globe size={18} className="text-blue-400" />
        <span className="font-semibold text-sm tracking-tight">RealmGenesis 3D</span>
      </div>
      <Panel className="flex-1 min-h-0">{make}</Panel>
    </aside>

    {/* Canvas + floating panels */}
    <div className="relative flex-1 min-w-0">
      <div className="absolute inset-0">{canvas}</div>

      {/* View — top strip */}
      <div className="absolute top-3 left-3 right-[17rem]">
        <div className={`${PANEL} flex items-center gap-3 px-3 py-2`}>
          {view}
          <button
            onClick={() => onSetEditing(!editing)}
            aria-pressed={editing}
            className={`ml-auto shrink-0 inline-flex items-center gap-1.5 px-3 py-1.5 rounded-md text-[11px] font-medium border transition-colors ${FOCUS}
              ${editing
                ? 'bg-amber-500 text-black border-amber-400'
                : 'bg-gray-800 text-gray-300 border-gray-700 hover:bg-gray-700 hover:text-white'}`}
          >
            <Pencil size={12} /> Edit
          </button>
        </div>
      </div>

      {/* Read — floating stack, right edge */}
      <div className="absolute top-3 right-3 bottom-3 w-64 flex flex-col gap-3 overflow-y-auto">
        {read.map(card => (
          <Panel key={card.key} title={card.title} className="shrink-0">
            {card.node}
          </Panel>
        ))}
      </div>

      {/* Do — contextual bottom bar */}
      {editing && (
        <div className="absolute bottom-4 left-1/2 -translate-x-1/2 max-w-[calc(100%_-_18rem)] rg-rise">
          <div className={`${PANEL} px-3 py-2 shadow-xl`}>{doTools}</div>
        </div>
      )}
    </div>
  </div>
);

export default WideShell;
