import { ViewMode } from '../../types';
import {
  coastlinePass, glyphPass, hillshadePass, landPass, oceanFillPass, paperPass,
} from './passes';
import { INKWASH_LABEL_THEME } from './labelTheme';
import { INKWASH_OVERLAY_INK } from './overlayInk';
import { FillPolicy, MapStyle, StylePalette } from './types';

/**
 * Pen and a single grey wash — a drawn map with no engraving in it.
 *
 * The distinction from Parchment is a real one, not a colour swap. Parchment
 * imitates an ENGRAVED plate: hatched sea, warm paper, shading kept off the
 * water so it does not fight the hatch. This imitates a WASH drawing, where a
 * brush carries the elevation and the pen only draws edges. So it has no ocean
 * hatch at all, its paper is neutral rather than warm, and its shading runs
 * harder and over the water as well as the land.
 */
const INKWASH_PALETTE: StylePalette = {
  paper: '#f4f1ea',
  ink: '#22221f',
  inkLight: '#8a8880',
  // With no hatch to carry it, the sea is separated from the paper by TONE
  // ALONE, so the step has to be a real one. #dbe2e1 was not: the ocean and the
  // bare-paper land read as the same field.
  sea: '#c6d2d1',
  seaHatch: '#a9b4b2',
  coast: '#1a1a17',
  // Neutral grey-green, the colour a diluted ink wash actually dries — and
  // DILUTE. #454a47 on near-white paper put so much contrast between adjacent
  // cells that the Voronoi polygons themselves became the dominant texture:
  // the land read as a honeycomb, not as relief. Parchment gets away with a
  // darker shadow because its paper is beige, so the same step is a smaller
  // one. A near-white sheet needs a weaker wash, not a stronger one.
  shadow: '#646b67',
  highlight: '#ffffff',
  ice: '#fbfdfe',
  // Cool dark grey. The desk under a wash drawing is a studio surface, not the
  // wood-and-leather table parchment sits on.
  desk: '#23262b',
  deskEdge: '#0a0b0d',
  deskShadow: '#0a0b0d',
};

const BARE: ViewMode[] = ['satellite', 'biome', 'height_bw'];
const CATEGORICAL: ViewMode[] = ['political', 'province', 'culture', 'religion', 'plates'];

/** Spec §4, same split as parchment: a ramp mode's fill IS its content. */
const inkWashFillPolicy = (mode: ViewMode): FillPolicy => {
  if (BARE.includes(mode)) return 'bare';
  if (CATEGORICAL.includes(mode)) return 'categorical';
  return 'ramp';
};

export const styleInkWash: MapStyle = {
  id: 'inkwash',
  name: 'Ink & Wash',
  palette: INKWASH_PALETTE,
  labelTheme: INKWASH_LABEL_THEME,
  overlayInk: INKWASH_OVERLAY_INK,
  fillPolicy: inkWashFillPolicy,
  passes: [
    paperPass(INKWASH_PALETTE, 'inkwash'),
    oceanFillPass(INKWASH_PALETTE),
    landPass(inkWashFillPolicy, INKWASH_PALETTE),
    // Land only, and no harder than parchment despite the wash being this
    // style's entire elevation story: the strength comes from the shadow tone,
    // and pushing the opacity as well brought the cell polygons straight back.
    //
    // Shading the sea bed here was tried and reverted. `hillshadePass` already
    // records why parchment turns it off, and the same thing happens on neutral
    // paper: the sea came out as mottled grey blotches that read as a stain on
    // the sheet, and they wiped out what little tonal separation the coastline
    // had to work with. A wash drawing washes the land and leaves the water
    // flat; that turns out to be the practical rule as well as the period one.
    hillshadePass(INKWASH_PALETTE, 0.6, true),
    coastlinePass(INKWASH_PALETTE),
    glyphPass(INKWASH_PALETTE),
  ],
};
