import { ViewMode } from '../../types';
import { FillPolicy, MapStyle } from './types';

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
  palette: {
    paper: '#e8d9b5',
    ink: '#3b2f1c',
    inkLight: '#8a7550',
    sea: '#c9bf9a',
    seaHatch: '#a89b74',
    coast: '#2e2414',
  },
  fillPolicy: parchmentFillPolicy,
};
