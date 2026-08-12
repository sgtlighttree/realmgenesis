// Leaf palette module. Imports NOTHING but types, on purpose: utils/worldGen.ts
// pulls this in, and worldGen runs inside a Web Worker where three.js is dead
// weight re-parsed on every spawn. Keep this file dependency-free.

// Perceptually distinct palette for faction coloring — maximum contrast across hue, lightness
export const FACTION_COLORS = [
  '#e53935', // vivid red
  '#43a047', // vivid green
  '#1e88e5', // vivid blue
  '#fb8c00', // vivid orange
  '#8e24aa', // vivid purple
  '#00acc1', // cyan
  '#f06292', // pink
  '#6d4c41', // brown
  '#c0ca33', // lime
  '#546e7a', // steel blue
  '#00897b', // teal
  '#fdd835', // yellow
  '#d81b60', // deep pink
  '#039be5', // light blue
  '#558b2f', // dark green
  '#6200ea', // deep purple
  '#ff6f00', // amber
  '#78909c', // blue-grey
];

// Muted, desaturated palette for the culture layer (C1) — deliberately
// softer than FACTION_COLORS so faction borders (drawn on top, in political
// mode) stay the visually dominant political layer.
export const CULTURE_COLORS = [
  '#a1887f', // muted taupe
  '#90a4ae', // muted blue-grey
  '#a5d6a7', // muted sage green
  '#ce93d8', // muted lavender
  '#ffcc80', // muted amber
  '#80cbc4', // muted teal
  '#ef9a9a', // muted rose
  '#c5cae9', // muted periwinkle
];

// Saturated, distinct palette for organized religions (C2) — deliberately
// more vivid than CULTURE_COLORS so a spreading faith reads clearly against
// the muted culture-layer shading of the folk faith it displaces.
export const RELIGION_COLORS = [
  '#ffd700', // gold
  '#7b1fa2', // deep violet
  '#c62828', // crimson
  '#00838f', // deep teal
  '#3949ab', // indigo
  '#ef6c00', // burnt orange
];

// Precomputed darkenForFolk(CULTURE_COLORS[i]), index-aligned.
//
// Generated once from the THREE implementation rather than ported, because
// THREE.Color applies an sRGB -> working-colorspace conversion in setHex that a
// hand-written HSL port would not reproduce exactly. tests/palette.test.ts
// recomputes these via THREE and fails if a three.js upgrade shifts them.
export const FOLK_COLORS = [
  '#857671',
  '#7e8a91',
  '#92b293',
  '#ab88b1',
  '#ddae64',
  '#79a5a1',
  '#cd8282',
  '#9ea5ce',
];
