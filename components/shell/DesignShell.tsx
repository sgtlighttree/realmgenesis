import React, { useEffect, useState } from 'react';
import WideShell from './WideShell';
import NarrowShell from './NarrowShell';
import {
  MakePanel, ViewPanel, DoPanel, READ_CARDS, PlaceholderGlobe, ShellProps,
} from './shellKit';

/**
 * DesignShell — the F1 "design debug" harness (reachable via ?shell=1).
 *
 * Owns two things a real App would own: which shell is active (breakpoint) and
 * whether editing is on (Esc exits). Composes the stub panels ONCE and hands the
 * finished elements to whichever shell positions them — the composition pattern
 * that keeps panel props one hop from the owner and out of the shells entirely.
 */
type Override = 'auto' | 'wide' | 'narrow';
const BREAKPOINT = 768;

const DesignShell: React.FC = () => {
  const [override, setOverride] = useState<Override>('auto');
  const [width, setWidth] = useState(() => window.innerWidth);
  const [editing, setEditing] = useState(false);

  useEffect(() => {
    const onResize = () => setWidth(window.innerWidth);
    window.addEventListener('resize', onResize);
    return () => window.removeEventListener('resize', onResize);
  }, []);

  useEffect(() => {
    const onKey = (e: KeyboardEvent) => { if (e.key === 'Escape') setEditing(false); };
    window.addEventListener('keydown', onKey);
    return () => window.removeEventListener('keydown', onKey);
  }, []);

  const effective: 'wide' | 'narrow' =
    override === 'auto' ? (width < BREAKPOINT ? 'narrow' : 'wide') : override;

  const shellProps: ShellProps = {
    make: <MakePanel />,
    view: <ViewPanel />,
    read: READ_CARDS,
    doTools: <DoPanel />,
    canvas: <PlaceholderGlobe />,
    editing,
    onSetEditing: setEditing,
  };

  return (
    <div className="w-full h-full flex flex-col bg-neutral-900 text-gray-200 font-sans">
      {/* Dev bar — not part of the design, just the harness controls. */}
      <div className="shrink-0 h-9 flex items-center gap-3 px-3 border-b border-gray-800 bg-gray-950 text-[11px]">
        <span className="font-mono text-gray-400 uppercase tracking-wider">F1 design debug</span>
        <div className="flex rounded-md overflow-hidden border border-gray-700">
          {(['auto', 'wide', 'narrow'] as Override[]).map(o => (
            <button key={o} onClick={() => setOverride(o)}
              className={`px-2.5 py-1 capitalize transition-colors ${override === o ? 'bg-blue-600 text-white' : 'bg-gray-800 text-gray-400 hover:text-white'}`}>
              {o}
            </button>
          ))}
        </div>
        <span className="text-gray-500">→ {effective === 'wide' ? 'A · Tidy' : 'C · Studio'}</span>
        {editing && <span className="ml-auto text-amber-400">editing — Esc to exit</span>}
      </div>

      {/* Stage */}
      <div className="relative flex-1 min-h-0 flex justify-center">
        <div className={effective === 'narrow'
          ? 'relative w-full max-w-[420px] border-x border-gray-800'
          : 'relative flex-1 min-w-0'}>
          {effective === 'wide'
            ? <WideShell {...shellProps} />
            : <NarrowShell {...shellProps} />}
        </div>
      </div>
    </div>
  );
};

export default DesignShell;
