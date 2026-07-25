import React, { useEffect, useLayoutEffect, useRef, useState } from 'react';
import { Sliders, Eye, Pencil, PanelRight } from 'lucide-react';
import { Panel, PANEL, FOCUS, BUCKET, Bucket, ShellProps } from './shellKit';

/**
 * NarrowShell — the "C · Studio" mobile fold. The four docks collapse into a
 * bottom tab bar; each tab opens a sheet over the globe (one panel at a time).
 * The Do tab IS edit mode. `openTab` is ephemeral presentation state and lives
 * here, not in App (per the state-ownership contract).
 */
type Tab = Bucket | null;

const TABS: { key: Bucket; icon: React.ReactNode }[] = [
  { key: 'make', icon: <Sliders size={16} /> },
  { key: 'view', icon: <Eye size={16} /> },
  { key: 'do', icon: <Pencil size={16} /> },
  { key: 'read', icon: <PanelRight size={16} /> },
];

const NarrowShell: React.FC<ShellProps> = ({
  make, view, read, doTools, canvas, onSetEditing,
}) => {
  const [openTab, setOpenTab] = useState<Tab>(null);
  const sheetBodyRef = useRef<HTMLDivElement>(null);

  // The sheet body is the SAME DOM node for every tab — only its children swap
  // — so scrollTop persists across tab switches and Make would open scrolled to
  // wherever the last tab was left, hiding Render Mode / Seed / Resolution.
  //
  // Reset twice: once before paint, and once on the next frame. The single
  // pre-paint reset alone left ~50px of drift, because the panel's content is
  // still settling (scroll anchoring re-applies an offset after layout).
  useLayoutEffect(() => {
    const el = sheetBodyRef.current;
    if (!el) return;
    el.scrollTop = 0;
    const raf = requestAnimationFrame(() => { el.scrollTop = 0; });
    return () => { cancelAnimationFrame(raf); };
  }, [openTab]);

  // The Do sheet is the edit context; opening it enters edit mode.
  useEffect(() => { onSetEditing(openTab === 'do'); }, [openTab, onSetEditing]);

  const sheetBody: Record<Bucket, React.ReactNode> = {
    make, view, do: doTools,
    read: (
      <div className="flex flex-col gap-3">
        {read.map(card => (
          <Panel key={card.key} title={card.title} className="shrink-0">{card.node}</Panel>
        ))}
      </div>
    ),
  };

  return (
    <div className="absolute inset-0 flex flex-col bg-black text-gray-200 font-sans">
      {/* Canvas hero */}
      <div className="relative flex-1 min-h-0">
        <div className="absolute inset-0">{canvas}</div>

        <div className="absolute top-3 left-3 text-xs font-semibold tracking-tight text-white/80">
          RealmGenesis 3D
        </div>

        {/* Sheet — one panel at a time, over the globe */}
        {openTab && (
          <div className={`absolute inset-x-2 bottom-2 max-h-[52%] flex flex-col ${PANEL} shadow-2xl rg-rise`}>
            <div className="flex items-center gap-2 px-3 py-2.5 border-b border-gray-800/80">
              <span className="text-xs font-semibold text-gray-100 tracking-tight">{BUCKET[openTab].label}</span>
              <button onClick={() => setOpenTab(null)} aria-label="Close panel"
                className={`ml-auto -mr-1 grid place-items-center w-6 h-6 rounded text-gray-500 hover:text-white hover:bg-gray-800 text-lg leading-none ${FOCUS}`}>×</button>
            </div>
            <div ref={sheetBodyRef} className="min-h-0 overflow-auto p-3">{sheetBody[openTab]}</div>
          </div>
        )}
      </div>

      {/* Bottom tab bar — the four docks folded */}
      <nav className="shrink-0 flex border-t border-gray-800 bg-gray-950">
        {TABS.map(t => {
          const active = openTab === t.key;
          return (
            <button key={t.key}
              onClick={() => setOpenTab(prev => (prev === t.key ? null : t.key))}
              aria-pressed={active}
              className={`flex-1 flex flex-col items-center gap-1 py-2.5 border-t-2 transition-colors ${FOCUS}
                ${active
                  ? 'text-blue-400 border-blue-500 bg-gray-900/60'
                  : 'text-gray-500 border-transparent hover:text-gray-300'}`}
            >
              {t.icon}
              <span className="text-[9px] font-mono uppercase tracking-wide">{BUCKET[t.key].label}</span>
            </button>
          );
        })}
      </nav>
    </div>
  );
};

export default NarrowShell;
