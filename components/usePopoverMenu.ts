import { useCallback, useEffect, useLayoutEffect, useRef, useState } from 'react';

/**
 * The positioning and dismissal machinery shared by every dropdown in the app.
 *
 * Extracted from `Select` when a second dropdown was needed (the overlays menu,
 * which is multi-select and so cannot be a listbox). It deliberately carries
 * ONLY the parts that are the same for any popover: where the menu goes, and
 * what closes it. Listbox semantics — an active index, type-ahead, commit — stay
 * in `Select`, because a checkbox menu needs different ones. Pulling those in
 * here is how a shared hook turns into a second implementation with a flag.
 */

const MENU_MAX_H = 320;
const GAP = 4;

/**
 * Should a scroll event close the menu?
 *
 * The scroll listener is registered in the CAPTURE phase, deliberately: the menu
 * is portaled and positioned `fixed`, so a scroll in ANY ancestor scroll
 * container can move the trigger out from under it and leave the menu stranded
 * in mid-air. Capture is the only way to see those scrolls.
 *
 * The cost of capture is that it also sees the menu scrolling ITSELF — the
 * option list is `overflow-y-auto`. Without this guard, trackpad-scrolling the
 * options or dragging their scrollbar closed the dropdown instantly, and the
 * remainder of the gesture fell through to the 2D/3D viewport underneath,
 * rotating the globe or panning the map.
 *
 * A scroll inside the menu never moves the trigger, so it must never close.
 * The sibling `pointerdown` handler already applies the same containment test;
 * this brings the scroll handler in line with it.
 */
export const shouldCloseOnScroll = (
  target: unknown,
  menu: { contains: (node: never) => boolean } | null,
): boolean => {
  if (!menu || target == null) return true;
  return !menu.contains(target as never);
};

export interface PopoverMenu {
  open: boolean;
  /** Open the menu. */
  openMenu: () => void;
  /** Close it; returns focus to the trigger unless told not to. */
  close: (refocus?: boolean) => void;
  setOpen: (open: boolean) => void;
  triggerRef: React.RefObject<HTMLButtonElement | null>;
  menuRef: React.RefObject<HTMLDivElement | null>;
  /** `fixed` placement for the portaled menu, recomputed on every open. */
  menuStyle: React.CSSProperties;
}

export function usePopoverMenu(): PopoverMenu {
  const [open, setOpen] = useState(false);
  const [menuStyle, setMenuStyle] = useState<React.CSSProperties>({});
  const triggerRef = useRef<HTMLButtonElement>(null);
  const menuRef = useRef<HTMLDivElement>(null);

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
    const onScroll = (e: Event) => {
      if (shouldCloseOnScroll(e.target, menuRef.current)) setOpen(false);
    };
    const onResize = () => { setOpen(false); };
    document.addEventListener('pointerdown', onPointerDown, true);
    window.addEventListener('resize', onResize);
    window.addEventListener('scroll', onScroll, true);
    return () => {
      document.removeEventListener('pointerdown', onPointerDown, true);
      window.removeEventListener('resize', onResize);
      window.removeEventListener('scroll', onScroll, true);
    };
  }, [open]);

  const openMenu = useCallback(() => { setOpen(true); }, []);
  const close = useCallback((refocus = true) => {
    setOpen(false);
    if (refocus) triggerRef.current?.focus();
  }, []);

  return { open, openMenu, close, setOpen, triggerRef, menuRef, menuStyle };
}
