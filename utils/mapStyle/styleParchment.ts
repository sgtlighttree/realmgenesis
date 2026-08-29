import { ViewMode } from '../../types';
import {
  coastlinePass, glyphPass, hillshadePass, landPass, oceanFillPass, oceanHatchPass, paperPass,
} from './passes';
import { FillPolicy, MapStyle, StylePalette } from './types';

/**
 * Warm, low-contrast paper. `shadow` and `highlight` are sepia rather than
 * black/white so hillshading reads as drawn relief instead of dirt on the page.
 */
const PARCHMENT_PALETTE: StylePalette = {
  paper: '#e8d9b5',
  ink: '#3b2f1c',
  inkLight: '#8a7550',
  sea: '#c9bf9a',
  seaHatch: '#a89b74',
  coast: '#2e2414',
  shadow: '#5c4a2e',
  highlight: '#fbf1d8',
};

const BARE: ViewMode[] = ['satellite', 'biome', 'height_bw'];
const CATEGORICAL: ViewMode[] = ['political', 'province', 'culture', 'religion', 'plates'];

/**
 * Spec §4. Continuous-ramp modes (height, temperature, moisture, population)
 * keep their own fill: their entire information content IS the fill, so bare
 * paper would render them blank. Glyphs are suppressed there too — they would
 * fight the ramp.
 */
const parchmentFillPolicy = (mode: ViewMode): FillPolicy => {
  if (BARE.includes(mode)) return 'bare';
  if (CATEGORICAL.includes(mode)) return 'categorical';
  return 'ramp';
};

export const styleParchment: MapStyle = {
  id: 'parchment',
  name: 'Parchment',
  palette: PARCHMENT_PALETTE,
  fillPolicy: parchmentFillPolicy,
  passes: [
    paperPass(PARCHMENT_PALETTE, 'parchment'),
    oceanFillPass(PARCHMENT_PALETTE),
    landPass(parchmentFillPolicy, PARCHMENT_PALETTE),
    hillshadePass(PARCHMENT_PALETTE, 0.5),
    // Hatch AFTER hillshade (spec §7.4: shading sits under the hatching), and
    // over an ocean-only composite region so bare-paper land is never covered.
    oceanHatchPass(PARCHMENT_PALETTE),
    coastlinePass(PARCHMENT_PALETTE),
    glyphPass(PARCHMENT_PALETTE),
  ],
};
