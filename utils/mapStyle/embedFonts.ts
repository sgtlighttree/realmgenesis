import cinzel400 from '@fontsource/cinzel/files/cinzel-latin-400-normal.woff2?url';
import cinzel700 from '@fontsource/cinzel/files/cinzel-latin-700-normal.woff2?url';
import fell400 from '@fontsource/im-fell-english/files/im-fell-english-latin-400-normal.woff2?url';
import fellItalic from '@fontsource/im-fell-english/files/im-fell-english-latin-400-italic.woff2?url';

import { LabelTheme } from './labelTheme';

interface EmbeddableFace {
  family: string;
  weight: number;
  style: 'normal' | 'italic';
  url: string;
}

/**
 * The webfonts an exported SVG can carry with it.
 *
 * Imported as `?url` rather than `?inline` on purpose: `?inline` would base64
 * roughly 150KB of woff2 into the main bundle for every visitor, whether or not
 * they ever export an SVG. As URLs the files are emitted as assets and fetched
 * only when someone actually exports.
 */
const EMBEDDABLE_FACES: EmbeddableFace[] = [
  { family: 'Cinzel', weight: 400, style: 'normal', url: cinzel400 },
  { family: 'Cinzel', weight: 700, style: 'normal', url: cinzel700 },
  { family: 'IM Fell English', weight: 400, style: 'normal', url: fell400 },
  { family: 'IM Fell English', weight: 400, style: 'italic', url: fellItalic },
];

const toBase64 = (buf: ArrayBuffer): string => {
  const bytes = new Uint8Array(buf);
  let s = '';
  // Chunked: String.fromCharCode(...bytes) blows the argument limit on a
  // 60KB font.
  for (let i = 0; i < bytes.length; i += 0x8000) {
    s += String.fromCharCode(...bytes.subarray(i, i + 0x8000));
  }
  return btoa(s);
};

/**
 * `@font-face` rules, with the woff2 inlined as data URIs, for the faces a
 * theme actually names.
 *
 * **Why an SVG needs this.** `exportSVG` writes `font-family` by name, and an
 * SVG is not a web page — it carries no stylesheet and cannot reach the app's
 * fonts. Opened anywhere that does not happen to have Cinzel and IM Fell
 * installed, every label silently falls back to a system serif and the export
 * does not look like the map that produced it. Embedding is the only way the
 * file is self-contained.
 *
 * Returns '' when nothing matches (the default theme has no webfonts) or when a
 * fetch fails: a label in a fallback face beats a failed export.
 */
export const embeddedFontCss = async (theme: LabelTheme): Promise<string> => {
  if (theme.faces.length === 0) return '';
  const wanted = EMBEDDABLE_FACES.filter(f =>
    theme.displayFamily.includes(f.family) || theme.geoFamily.includes(f.family));
  if (wanted.length === 0) return '';

  const rules = await Promise.all(wanted.map(async (f) => {
    try {
      const res = await fetch(f.url);
      if (!res.ok) return '';
      const b64 = toBase64(await res.arrayBuffer());
      return `@font-face{font-family:'${f.family}';font-style:${f.style};`
        + `font-weight:${f.weight};src:url(data:font/woff2;base64,${b64}) format('woff2');}`;
    } catch {
      return '';
    }
  }));
  return rules.filter(Boolean).join('');
};
