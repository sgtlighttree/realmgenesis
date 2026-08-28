import { describe, it, expect } from 'vitest';

import { shouldCloseOnScroll } from '../components/Select';

// A stand-in for a DOM node with `contains`. The real predicate only ever calls
// `contains`, so this needs no DOM — the repo has no jsdom and does not need it
// for this.
const nodeWith = (children: string[]) => ({
  contains: (n: unknown) => children.includes(n as string),
});

describe('shouldCloseOnScroll', () => {
  // The bug: Select closes on ANY scroll event because the listener is
  // registered in the capture phase, so it also catches the menu scrolling
  // ITSELF. Trackpad-scrolling the option list, or dragging its scrollbar,
  // closed the dropdown — and the rest of the gesture then fell through to the
  // 2D/3D viewport underneath.
  it('does NOT close when the scroll came from inside the menu', () => {
    const menu = nodeWith(['menu', 'option-3']);
    expect(shouldCloseOnScroll('menu', menu)).toBe(false);
    expect(shouldCloseOnScroll('option-3', menu)).toBe(false);
  });

  // The listener exists for a real reason: the menu is portaled and positioned
  // `fixed`, so a scroll that moves the TRIGGER would leave the menu stranded
  // in mid-air. That case must still close.
  it('closes when an outside scroll container moves', () => {
    const menu = nodeWith(['menu']);
    expect(shouldCloseOnScroll('sidebar', menu)).toBe(true);
    expect(shouldCloseOnScroll('document', menu)).toBe(true);
  });

  it('closes when there is no menu or no target', () => {
    expect(shouldCloseOnScroll('anything', null)).toBe(true);
    expect(shouldCloseOnScroll(null, nodeWith(['menu']))).toBe(true);
    expect(shouldCloseOnScroll(undefined, nodeWith(['menu']))).toBe(true);
  });
});
