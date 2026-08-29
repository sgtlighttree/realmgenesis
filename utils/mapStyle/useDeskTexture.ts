import { useEffect, useState } from 'react';

/** Decoded textures, kept for the page's life: a style switch must not re-fetch. */
const cache = new Map<string, HTMLImageElement>();

/**
 * Load a style's desk texture, and hand back the image only once it is
 * decoded — never before.
 *
 * **Why this exists.** Same failure as the webfonts (see `useLabelFonts`): an
 * `HTMLImageElement` that has not finished decoding draws NOTHING through
 * `drawImage`/`createPattern` and throws nothing either, so an early paint
 * looks exactly like a style with no texture and nothing reports a problem.
 * Returning the image rather than a ready flag makes the dependency real —
 * callers already pass it to the draw call, so a new identity repaints them by
 * the ordinary rules.
 *
 * Returns `null` for a style with no texture, and for one whose texture fails
 * to load: a missing desk graphic must fall back to the palette colour, not
 * blank the viewport.
 */
export const useDeskTexture = (src?: string): HTMLImageElement | null => {
  const [image, setImage] = useState<HTMLImageElement | null>(() => (src ? cache.get(src) ?? null : null));

  useEffect(() => {
    if (!src) { setImage(null); return; }
    const cached = cache.get(src);
    if (cached) { setImage(cached); return; }

    let live = true;
    setImage(null);
    const img = new Image();
    img.decoding = 'async';
    const ready = () => {
      if (!img.naturalWidth) return;
      cache.set(src, img);
      if (live) setImage(img);
    };
    // `decode()` resolves only when the bitmap is paintable, which is the whole
    // point; it rejects on some browsers for a cross-origin or already-complete
    // image, so `onload` stays as the fallback.
    img.onload = ready;
    img.onerror = () => { if (live) setImage(null); };
    img.src = src;
    img.decode?.().then(ready).catch(() => { /* onload covers it */ });

    return () => { live = false; };
  }, [src]);

  return image;
};
