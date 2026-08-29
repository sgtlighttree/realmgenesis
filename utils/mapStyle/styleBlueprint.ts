import { ViewMode } from '../../types';
import {
  coastlinePass, glyphPass, hillshadePass, landPass, oceanFillPass, oceanHatchPass,
  paperPass,
} from './passes';
import { BLUEPRINT_LABEL_THEME } from './labelTheme';
import { BLUEPRINT_OVERLAY_INK } from './overlayInk';
import { FillPolicy, MapStyle, StylePalette } from './types';

/**
 * Cyanotype. White line work on process blue, the way a world would look if it
 * had been drawn on a board and run through a blueprint machine.
 *
 * The one thing a cyanotype inverts relative to every other style here: the
 * PAPER is the dark tone and the ink is the light one. `shadow` and `highlight`
 * follow — the shade is a deeper blue and the highlight a pale one, so relief
 * still reads as a light source rather than as two unrelated tints.
 */
const BLUEPRINT_PALETTE: StylePalette = {
  paper: '#1a5183',
  ink: '#e6f0f8',
  inkLight: '#8bb0cf',
  // A wide gap, not a shade or two. Land is BARE paper under this style's fill
  // policy, so the only things separating a continent from the sea are this
  // step and the coastline drawn over it. The first pair were four percent
  // apart in luminance and the map read as one navy field with white outlines
  // on it — the same mistake parchment's `sea` comment records.
  sea: '#0a2338',
  seaHatch: '#2e618e',
  coast: '#ffffff',
  shadow: '#08243c',
  highlight: '#bcd9f2',
  // Pale, but NOT paper-white: the coastline is pure white under this style, and
  // an ice cap at #dceaf6 swallowed it whole — the polar coasts simply had no
  // visible line. Backed off until the coast reads over the ice.
  ice: '#bcd4e8',
  // A drafting table rather than a desk: neutral dark grey, no warmth. Wood
  // under a blueprint would read as a period map that lost its colour.
  desk: '#1a1f26',
  deskEdge: '#080a0d',
  deskShadow: '#03060a',
};

const BARE: ViewMode[] = ['satellite', 'biome', 'height_bw'];
const CATEGORICAL: ViewMode[] = ['political', 'province', 'culture', 'religion', 'plates'];

/** Spec §4, same split as parchment: a ramp mode's fill IS its content. */
const blueprintFillPolicy = (mode: ViewMode): FillPolicy => {
  if (BARE.includes(mode)) return 'bare';
  if (CATEGORICAL.includes(mode)) return 'categorical';
  return 'ramp';
};

export const styleBlueprint: MapStyle = {
  id: 'blueprint',
  name: 'Blueprint',
  palette: BLUEPRINT_PALETTE,
  labelTheme: BLUEPRINT_LABEL_THEME,
  overlayInk: BLUEPRINT_OVERLAY_INK,
  fillPolicy: blueprintFillPolicy,
  passes: [
    paperPass(BLUEPRINT_PALETTE, 'blueprint'),
    oceanFillPass(BLUEPRINT_PALETTE),
    landPass(blueprintFillPolicy, BLUEPRINT_PALETTE),
    // Lighter than parchment's 0.55. The palette's shadow is already a strong
    // dark blue against a mid-blue paper, so the same opacity reads as a bruise.
    hillshadePass(BLUEPRINT_PALETTE, 0.4, true),
    oceanHatchPass(BLUEPRINT_PALETTE),
    coastlinePass(BLUEPRINT_PALETTE),
    glyphPass(BLUEPRINT_PALETTE),
  ],
};
