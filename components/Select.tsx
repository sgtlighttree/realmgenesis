import React, { useEffect, useRef, useState } from 'react';
import { createPortal } from 'react-dom';
import { Check, ChevronDown } from 'lucide-react';

import { usePopoverMenu } from './usePopoverMenu';

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

/**
 * Re-exported from `usePopoverMenu`, where it now lives alongside the scroll
 * listener that uses it. Kept here because this is the module that documented
 * the fix and the module its test imports from.
 */
export { shouldCloseOnScroll } from './usePopoverMenu';

export interface SelectOption<T extends string> {
  value: T;
  label: string;
}

interface SelectProps<T extends string> {
  value: T;
  options: SelectOption<T>[];
  onChange: (value: T) => void;
  /**
   * Accessible name. Used as `aria-label` UNLESS `labelledBy` is given: a
   * control with a visible caption must take its name from that caption, or it
   * has two accessible names and a screen reader announces the wrong one.
   */
  label: string;
  /** Id of a visible caption element. Replaces `aria-label` when present. */
  labelledBy?: string;
  className?: string;
  /** Extra classes for the trigger button. */
  triggerClassName?: string;
}

function Select<T extends string>({
  value, options, onChange, label, labelledBy, className = '', triggerClassName = '',
}: SelectProps<T>) {
  const { open, close, setOpen, triggerRef, menuRef, menuStyle } = usePopoverMenu();
  const [activeIndex, setActiveIndex] = useState(0);

  const optionRefs = useRef<(HTMLDivElement | null)[]>([]);
  const typeahead = useRef({ query: '', at: 0 });

  const selectedIndex = Math.max(0, options.findIndex(o => o.value === value));
  const selected = options[selectedIndex];
  const listId = React.useId();

  // Keep the active option in view during keyboard navigation.
  useEffect(() => {
    if (!open) return;
    optionRefs.current[activeIndex]?.scrollIntoView({ block: 'nearest' });
  }, [open, activeIndex]);

  // Seeds the active index, which is listbox behaviour and so stays here rather
  // than in the shared hook.
  const openMenu = (index = selectedIndex) => {
    setActiveIndex(index);
    setOpen(true);
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
        aria-label={labelledBy ? undefined : label}
        aria-labelledby={labelledBy}
        onClick={() => { open ? close(false) : openMenu(); }}
        onKeyDown={onKeyDown}
        className={`inline-flex items-center gap-1.5 border border-edge bg-surface-raised px-2 py-1.5 text-[11px] text-ink transition-colors hover:border-edge-strong hover:text-ink-strong focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-brand/70 ${triggerClassName}`}
      >
        <span className="truncate">{selected?.label ?? ''}</span>
        <ChevronDown
          size={12}
          className={`shrink-0 text-ink-muted transition-transform duration-150 ${open ? 'rotate-180' : ''}`}
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
          // `overscroll-contain` stops scroll CHAINING: at the top or bottom of the
          // list, further wheel/touch scrolling would otherwise propagate to the
          // page and reach the R3F canvas, zooming the globe while the menu is open.
          className="overflow-y-auto overscroll-contain border border-edge bg-surface py-1 shadow-2xl rg-rise"
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
                  i === activeIndex ? 'bg-brand-strong text-ink-strong' : 'text-ink-soft'
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
