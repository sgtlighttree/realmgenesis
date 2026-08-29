import { useEffect, useState } from 'react';

import { ensureLabelFonts, LabelTheme } from './labelTheme';

/**
 * Load a theme's faces, and hand back a theme whose IDENTITY changes once they
 * are ready.
 *
 * **Why this exists.** Canvas2D does not wait for webfonts. `ctx.font` accepts a
 * family that has not finished loading, falls back to the next entry in the
 * stack, and draws without error — so a map painted before the fonts land looks
 * EXACTLY like an unstyled one, and nothing anywhere reports a problem. It
 * presents as "the style did nothing", which is a diagnosis this project has
 * already spent sessions on.
 *
 * **Why it returns a theme rather than a ready flag.** The first version
 * returned a counter that callers listed in their dependency arrays without
 * reading. That works, but `react-hooks/exhaustive-deps` correctly calls an
 * unread dependency unnecessary, and silencing it needs a disable comment at
 * every call site — which is a standing invitation for someone to delete the
 * "pointless" dep and quietly reintroduce the fallback-face bug. Returning the
 * theme makes the dependency real: callers already pass it to the draw call, so
 * a new identity repaints them by the ordinary rules.
 *
 * The returned object is `===` to the input until the faces resolve.
 */
export const useLabelFonts = (theme: LabelTheme): LabelTheme => {
  const [loaded, setLoaded] = useState<LabelTheme | null>(null);

  useEffect(() => {
    let live = true;
    setLoaded(null);
    ensureLabelFonts(theme).then(() => {
      // A theme with no webfonts is ready on arrival; a fresh identity there
      // would force one wasted repaint on every style change.
      if (live && theme.faces.length > 0) setLoaded({ ...theme });
    });
    return () => { live = false; };
  }, [theme]);

  return loaded ?? theme;
};
