import { describe, it, expect } from 'vitest';

import { SCREEN_STROKE_PX, screenStroke } from '../components/map2d/strokeWidths';

// THE BUG: widths were written `Math.max(FLOOR, W / s)` under a transform that
// already carries `s` (VectorOverlay's `setTransform(dpr * s, ...)`). The floor
// is therefore in TRANSFORMED units while it was meant as a screen minimum, so
// once `W / s` fell below it the width pinned at the floor and the on-screen
// weight became `floor * s` — growing linearly with zoom. The coastline crossed
// over at s = 1.5 and rendered ~4x too heavy by full zoom, swamping the labels,
// which are constant size. Same defect the label `Math.max(8, ...)` floor had.
describe('screenStroke', () => {
  it('holds a constant on-screen width at every zoom', () => {
    for (const px of Object.values(SCREEN_STROKE_PX)) {
      for (const s of [0.6, 1, 1.5, 2, 4, 6]) {
        expect(screenStroke(px, s) * s, `${px}px at scale ${s}`).toBeCloseTo(px, 9);
      }
    }
  });

  it('never returns a non-positive width', () => {
    for (const s of [0, -1, 0.6, 6]) expect(screenStroke(1.5, s)).toBeGreaterThan(0);
  });

  // The regression in one assertion: the old expression at full zoom.
  it('is thinner than the old floored expression once the floor bound', () => {
    const s = 6;
    const old = Math.max(1, 1.5 / s);            // = 1 transformed unit
    const now = screenStroke(SCREEN_STROKE_PX.coast, s);
    expect(old * s).toBeCloseTo(6, 9);           // 6 screen px — the thick outline
    expect(now * s).toBeCloseTo(SCREEN_STROKE_PX.coast, 9);
    expect(now).toBeLessThan(old);
  });

  // A coastline must not out-weigh the text it sits beside. Labels render at
  // 12-14 screen px (LABEL_CONFIG base sizes), so a stroke in the same order of
  // magnitude reads as a border rather than a coast.
  it('keeps every stroke well under label weight', () => {
    for (const px of Object.values(SCREEN_STROKE_PX)) expect(px).toBeLessThanOrEqual(3);
  });
});
