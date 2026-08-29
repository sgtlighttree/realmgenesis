import { useEffect, useState } from 'react';

import { ensureLabelFonts, LabelTheme } from './labelTheme';

/**
 * Load a theme's faces and report when a repaint is due.
 *
 * **Why this exists at all.** Canvas2D does not wait for webfonts. `ctx.font`
 * accepts a family that has not finished loading, falls back to the next entry
 * in the stack, and draws without error — so a map painted before the fonts
 * land looks EXACTLY like an unstyled one, and nothing anywhere reports a
 * problem. Nobody would find that by reading the code; it presents as "the
 * style did nothing", which is a diagnosis this project has already spent two
 * sessions on.
 *
 * The returned number changes once per theme when its faces resolve. Put it in
 * the deps of whatever paints, so the fallback frame is replaced.
 */
export const useLabelFonts = (theme: LabelTheme): number => {
  const [tick, setTick] = useState(0);

  useEffect(() => {
    let live = true;
    ensureLabelFonts(theme).then(() => {
      // A theme with no webfonts resolves immediately; bumping anyway would
      // force one wasted repaint on every style change.
      if (live && theme.faces.length > 0) setTick(t => t + 1);
    });
    return () => { live = false; };
  }, [theme]);

  return tick;
};
