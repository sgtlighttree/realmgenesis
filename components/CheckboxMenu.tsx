import React from 'react';
import { createPortal } from 'react-dom';
import { ChevronDown } from 'lucide-react';

import { usePopoverMenu } from './usePopoverMenu';

/** Lucide icons take size/className; React.ElementType is too loose to accept them. */
type IconType = React.ComponentType<{ size?: number; className?: string }>;

export interface CheckboxMenuItem {
  key: string;
  label: string;
  icon?: IconType;
  checked: boolean;
  onChange: (checked: boolean) => void;
}

interface CheckboxMenuProps {
  items: CheckboxMenuItem[];
  /** Accessible name. Replaced by `labelledBy` when a visible caption exists. */
  label: string;
  labelledBy?: string;
  className?: string;
  triggerClassName?: string;
}

/**
 * A dropdown of independent checkboxes — the multi-select sibling of `Select`.
 *
 * It exists because `Select` cannot be it. `Select` is a listbox: one of N, and
 * choosing commits and closes. Nine overlay toggles are none of those things, so
 * bolting a mode onto it would have meant a flag threaded through its active
 * index, its type-ahead and its commit path. They share the part that is
 * genuinely the same — `usePopoverMenu`, which carries the portal positioning
 * and the hard-won scroll-containment fix — and nothing else.
 *
 * Two consequences of "independent", both deliberate:
 *
 * - **Clicking an item does NOT close the menu.** You come here to turn on three
 *   things, not one.
 * - **The trigger's text IS the active count.** This control replaced a row of
 *   chips where you could see at a glance that Rivers and Borders were on.
 *   A dropdown reading only "Overlays" would tell the user strictly less than
 *   what it replaced; "2 shown" is what keeps that from being a regression.
 *
 * ARIA is `menu`/`menuitemcheckbox`, not `listbox`/`option` — an option is a
 * choice among alternatives, and these are not alternatives.
 */
const CheckboxMenu: React.FC<CheckboxMenuProps> = ({
  items, label, labelledBy, className = '', triggerClassName = '',
}) => {
  const { open, close, setOpen, triggerRef, menuRef, menuStyle } = usePopoverMenu();
  const [activeIndex, setActiveIndex] = React.useState(0);
  const itemRefs = React.useRef<(HTMLDivElement | null)[]>([]);
  const listId = React.useId();

  const activeCount = items.filter(i => i.checked).length;

  React.useEffect(() => {
    if (!open) return;
    itemRefs.current[activeIndex]?.scrollIntoView({ block: 'nearest' });
  }, [open, activeIndex]);

  const toggle = (index: number): void => {
    const item = items[index];
    if (item) item.onChange(!item.checked);
  };

  const onKeyDown = (e: React.KeyboardEvent): void => {
    switch (e.key) {
      case 'ArrowDown':
      case 'ArrowUp': {
        e.preventDefault();
        const delta = e.key === 'ArrowDown' ? 1 : -1;
        if (!open) { setActiveIndex(0); setOpen(true); return; }
        setActiveIndex(i => Math.min(items.length - 1, Math.max(0, i + delta)));
        return;
      }
      case 'Home':
        if (open) { e.preventDefault(); setActiveIndex(0); }
        return;
      case 'End':
        if (open) { e.preventDefault(); setActiveIndex(items.length - 1); }
        return;
      case 'Enter':
      case ' ':
        e.preventDefault();
        // Toggles and STAYS OPEN. This is the one place the behaviour has to
        // diverge from Select, whose Enter commits and closes.
        if (open) toggle(activeIndex); else { setActiveIndex(0); setOpen(true); }
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
        aria-haspopup="menu"
        aria-expanded={open}
        aria-controls={open ? listId : undefined}
        aria-label={labelledBy ? undefined : label}
        aria-labelledby={labelledBy}
        onClick={() => { open ? close(false) : setOpen(true); }}
        onKeyDown={onKeyDown}
        className={`inline-flex items-center gap-1.5 border border-edge bg-surface-raised px-2 py-1.5 text-[11px] text-ink transition-colors hover:border-edge-strong hover:text-ink-strong focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-brand/70 ${triggerClassName}`}
      >
        {/* The summary IS the trigger's text. A noun here would only repeat the
            caption above it; what the caption cannot say is how many are on. */}
        <span className={`truncate ${activeCount > 0 ? 'text-ink' : 'text-ink-faint'}`}>
          {activeCount === 0 ? 'None' : `${activeCount} shown`}
        </span>
        <ChevronDown
          size={12}
          className={`shrink-0 text-ink-muted transition-transform duration-150 ${open ? 'rotate-180' : ''}`}
        />
      </button>

      {open && createPortal(
        <div
          ref={menuRef}
          id={listId}
          role="menu"
          aria-label={label}
          tabIndex={-1}
          style={{ ...menuStyle, zIndex: 60 }}
          // `overscroll-contain` stops scroll CHAINING: at the top or bottom of
          // the list, further wheel/touch scrolling would otherwise propagate to
          // the page and reach the R3F canvas, zooming the globe while open.
          className="overflow-y-auto overscroll-contain border border-edge bg-surface py-1 shadow-2xl rg-rise"
        >
          {items.map((item, i) => {
            const Icon = item.icon;
            return (
              <div
                key={item.key}
                ref={el => { itemRefs.current[i] = el; }}
                role="menuitemcheckbox"
                aria-checked={item.checked}
                onPointerEnter={() => { setActiveIndex(i); }}
                onClick={() => { toggle(i); }}
                className={`flex cursor-pointer items-center gap-2 px-2.5 py-1.5 text-[11px] whitespace-nowrap transition-colors ${
                  i === activeIndex ? 'bg-brand-strong text-ink-strong' : 'text-ink-soft'
                }`}
              >
                <input
                  type="checkbox"
                  checked={item.checked}
                  readOnly
                  tabIndex={-1}
                  // The row owns the click; the box is the visual state only, so
                  // it must not consume the event or fire a second toggle.
                  className="pointer-events-none bg-surface-hover"
                />
                {Icon ? <Icon size={11} className={item.checked ? 'text-brand-soft' : 'text-ink-faint'} /> : null}
                {item.label}
              </div>
            );
          })}
        </div>,
        document.body,
      )}
    </div>
  );
};

export default CheckboxMenu;
