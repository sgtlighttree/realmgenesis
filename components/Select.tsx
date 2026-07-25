import React, { useCallback, useEffect, useLayoutEffect, useRef, useState } from 'react';
import { createPortal } from 'react-dom';
import { Check, ChevronDown } from 'lucide-react';

/**
 * Select — a themed listbox replacing the native `<select>`.
 *
 * Why replace it at all: a native `<select>` renders its menu with OS chrome,
 * so on this near-black UI it drops a light macOS popup in the middle of the
 * app. That is a theming defect, not a missing feature.
 *
 * So the APPEARANCE is ours and the BEHAVIOUR stays exactly what users expect
 * from a select: type-ahead, Home/End, arrow navigation, Enter/Space to commit,
 * Esc to cancel, focus returned to the trigger on close, and full ARIA
 * listbox semantics. Reinventing the affordance would be worse than the
 * native control; only the paint is reinvented here.
 *
 * The menu is PORTALED to `document.body` and positioned `fixed`. An absolutely
 * positioned menu would be clipped by the View strip's own overflow and by the
 * mobile sheet's `overflow-auto` — the classic dropdown-in-a-scroll-container
 * bug. Portal + fixed escapes both stacking contexts.
 */

export interface SelectOption<T extends string> {
  value: T;
  label: string;
}

interface SelectProps<T extends string> {
  value: T;
  options: SelectOption<T>[];
  onChange: (value: T) => void;
  /** Accessible name — there is no visible <label> in the compact strip. */
  label: string;
  className?: string;
  /** Extra classes for the trigger button. */
  triggerClassName?: string;
}

const MENU_MAX_H = 320;
const GAP = 4;

function Select<T extends string>({
  value, options, onChange, label, className = '', triggerClassName = '',
}: SelectProps<T>) {
  const [open, setOpen] = useState(false);
  const [activeIndex, setActiveIndex] = useState(0);
  const [menuStyle, setMenuStyle] = useState<React.CSSProperties>({});

  const triggerRef = useRef<HTMLButtonElement>(null);
  const menuRef = useRef<HTMLDivElement>(null);
  const optionRefs = useRef<(HTMLDivElement | null)[]>([]);
  const typeahead = useRef({ query: '', at: 0 });

  const selectedIndex = Math.max(0, options.findIndex(o => o.value === value));
  const selected = options[selectedIndex];
  const listId = React.useId();

  const position = useCallback(() => {
    const trigger = triggerRef.current;
    if (!trigger) return;
    const r = trigger.getBoundingClientRect();
    const below = window.innerHeight - r.bottom - GAP;
    const above = r.top - GAP;
    // Flip up when the space below cannot hold a usefully tall menu.
    const flip = below < Math.min(MENU_MAX_H, 160) && above > below;
    const maxHeight = Math.min(MENU_MAX_H, (flip ? above : below) - 4);
    setMenuStyle({
      position: 'fixed',
      left: Math.round(Math.min(r.left, window.innerWidth - r.width - 8)),
      minWidth: Math.round(r.width),
      maxHeight: Math.max(96, Math.round(maxHeight)),
      ...(flip
        ? { bottom: Math.round(window.innerHeight - r.top + GAP) }
        : { top: Math.round(r.bottom + GAP) }),
    });
  }, []);

  useLayoutEffect(() => {
    if (open) position();
  }, [open, position]);

  // Close on outside interaction; reposition-by-closing on scroll/resize so the
  // menu can never detach from a trigger that has moved out from under it.
  useEffect(() => {
    if (!open) return;
    const onPointerDown = (e: PointerEvent) => {
      const t = e.target as Node;
      if (menuRef.current?.contains(t) || triggerRef.current?.contains(t)) return;
      setOpen(false);
    };
    const onScrollOrResize = () => { setOpen(false); };
    document.addEventListener('pointerdown', onPointerDown, true);
    window.addEventListener('resize', onScrollOrResize);
    window.addEventListener('scroll', onScrollOrResize, true);
    return () => {
      document.removeEventListener('pointerdown', onPointerDown, true);
      window.removeEventListener('resize', onScrollOrResize);
      window.removeEventListener('scroll', onScrollOrResize, true);
    };
  }, [open]);

  // Keep the active option in view during keyboard navigation.
  useEffect(() => {
    if (!open) return;
    optionRefs.current[activeIndex]?.scrollIntoView({ block: 'nearest' });
  }, [open, activeIndex]);

  const openMenu = (index = selectedIndex) => {
    setActiveIndex(index);
    setOpen(true);
  };

  const close = (refocus = true) => {
    setOpen(false);
    if (refocus) triggerRef.current?.focus();
  };

  const commit = (index: number) => {
    const opt = options[index];
    if (opt) onChange(opt.value);
    close();
  };

  const onKeyDown = (e: React.KeyboardEvent) => {
    // Printable single characters drive type-ahead, matching a native select.
    if (e.key.length === 1 && !e.metaKey && !e.ctrlKey && !e.altKey) {
      const now = Date.now();
      const t = typeahead.current;
      t.query = now - t.at > 700 ? e.key : t.query + e.key;
      t.at = now;
      const q = t.query.toLowerCase();
      const from = open ? activeIndex : selectedIndex;
      // Search after the cursor first so repeated keys cycle matches.
      const order = [...options.slice(from + 1), ...options.slice(0, from + 1)];
      const hit = order.find(o => o.label.toLowerCase().startsWith(q));
      if (hit) {
        const i = options.indexOf(hit);
        if (open) setActiveIndex(i); else onChange(hit.value);
      }
      e.preventDefault();
      return;
    }

    switch (e.key) {
      case 'ArrowDown':
      case 'ArrowUp': {
        e.preventDefault();
        const delta = e.key === 'ArrowDown' ? 1 : -1;
        if (!open) { openMenu(Math.min(options.length - 1, Math.max(0, selectedIndex + delta))); return; }
        setActiveIndex(i => Math.min(options.length - 1, Math.max(0, i + delta)));
        return;
      }
      case 'Home':
        if (open) { e.preventDefault(); setActiveIndex(0); }
        return;
      case 'End':
        if (open) { e.preventDefault(); setActiveIndex(options.length - 1); }
        return;
      case 'Enter':
      case ' ':
        e.preventDefault();
        if (open) commit(activeIndex); else openMenu();
        return;
      case 'Escape':
        if (open) { e.preventDefault(); e.stopPropagation(); close(); }
        return;
      case 'Tab':
        if (open) setOpen(false);
        return;
      default:
    }
  };

  return (
    <div className={`relative ${className}`}>
      <button
        ref={triggerRef}
        type="button"
        role="combobox"
        aria-haspopup="listbox"
        aria-expanded={open}
        aria-controls={open ? listId : undefined}
        aria-label={label}
        onClick={() => { open ? close(false) : openMenu(); }}
        onKeyDown={onKeyDown}
        className={`inline-flex items-center gap-1.5 rounded border border-gray-700 bg-gray-800 px-2 py-1.5 text-[11px] text-gray-200 transition-colors hover:border-gray-600 hover:text-white focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-blue-500/70 ${triggerClassName}`}
      >
        <span className="truncate">{selected?.label ?? ''}</span>
        <ChevronDown
          size={12}
          className={`shrink-0 text-gray-400 transition-transform duration-150 ${open ? 'rotate-180' : ''}`}
        />
      </button>

      {open && createPortal(
        <div
          ref={menuRef}
          id={listId}
          role="listbox"
          aria-label={label}
          tabIndex={-1}
          style={{ ...menuStyle, zIndex: 60 }}
          className="overflow-y-auto rounded-md border border-gray-700 bg-gray-900 py-1 shadow-2xl rg-rise"
        >
          {options.map((opt, i) => {
            const isSelected = opt.value === value;
            return (
              <div
                key={opt.value}
                ref={el => { optionRefs.current[i] = el; }}
                role="option"
                aria-selected={isSelected}
                onPointerEnter={() => { setActiveIndex(i); }}
                onClick={() => { commit(i); }}
                className={`flex cursor-pointer items-center gap-2 px-2.5 py-1.5 text-[11px] whitespace-nowrap transition-colors ${
                  i === activeIndex ? 'bg-blue-600 text-white' : 'text-gray-300'
                }`}
              >
                <Check
                  size={11}
                  className={`shrink-0 ${isSelected ? 'opacity-100' : 'opacity-0'}`}
                />
                {opt.label}
              </div>
            );
          })}
        </div>,
        document.body,
      )}
    </div>
  );
}

export default Select;
