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
  <div className="absolute inset-0 flex bg-black text-ink font-sans">

    {/* Make — flush left rail. No outer padding and no Panel wrapper: Controls
        supplies its own chrome and scroll region, so anything here is nesting. */}
    <aside className="w-72 shrink-0 h-full flex flex-col border-r border-edge-subtle bg-surface-sunken">
      <div className="flex items-center gap-2 px-3 py-2.5 border-b border-edge-subtle">
        <Globe size={16} className="text-brand-soft" />
        <span className="font-semibold text-[13px] tracking-tight">RealmGenesis 3D</span>
      </div>
      <div className="flex-1 min-h-0">{make}</div>
    </aside>

    {/* Canvas + docked panels */}
    <div className="relative flex-1 min-w-0 overflow-hidden">
      {/* The globe must centre on the VISIBLE gap between the rails, not on the
          element box. Insetting the canvas from the right did that but left a
          dead black gutter under the Read rail. Instead the canvas keeps full
          coverage and is shifted LEFT by the rail's footprint: it still paints
          to the right edge, the overspill on the left is clipped by the parent,
          and the centre lands on the visible centre. */}
      <div className="absolute inset-y-0 left-[-16.5rem] right-0">{canvas}</div>

      {/* View — top strip. Right edge clears the Read rail (16rem + 2×8px). */}
      <div className="absolute top-2 left-2 right-[17rem] z-chrome">
        <div className={`${PANEL} flex items-center gap-2 px-2 py-1.5`}>
          {view}
          <button
            onClick={() => onSetEditing(!editing)}
            aria-pressed={editing}
            className={`ml-auto shrink-0 inline-flex items-center gap-1.5 px-2.5 py-1.5 text-[11px] font-medium border transition-colors ${FOCUS}
              ${editing
                ? 'bg-warn text-black border-warn-soft'
                : 'bg-surface-raised text-ink-soft border-edge hover:bg-surface-hover hover:text-ink-strong'}`}
          >
            <Pencil size={12} /> Edit
          </button>
        </div>
      </div>

      {/* Read — docked stack, right edge */}
      <div className="absolute top-2 right-2 bottom-2 w-64 flex flex-col gap-2 overflow-y-auto z-chrome">
        {read.map(card => (
          <Panel key={card.key} title={card.title} collapsible={card.collapsible} className="shrink-0">
            {card.node}
          </Panel>
        ))}
      </div>

      {/* Do — contextual bottom bar. It lives inside the inset canvas column,
          so it centres with the globe and needs no rail allowance of its own. */}
      {editing && (
        <div className="absolute bottom-3 left-0 right-[16.5rem] z-chrome flex justify-center px-4">
          <div className={`${PANEL} px-2 py-1.5 shadow-xl rg-rise`}>{doTools}</div>
        </div>
      )}
    </div>
  </div>
);

export default WideShell;
