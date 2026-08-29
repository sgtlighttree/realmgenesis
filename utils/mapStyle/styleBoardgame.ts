import { ViewMode } from '../../types';
import {
  coastlinePass, glyphPass, landPass, oceanFillPass, paperPass,
} from './passes';
import { BOARDGAME_LABEL_THEME } from './labelTheme';
import { BOARDGAME_OVERLAY_INK } from './overlayInk';
import { FillPolicy, MapStyle, StylePalette } from './types';

/**
 * A printed board. Flat territory colour, heavy black rules, no texture of any
 * kind — the look of something screen-printed on card and read from a metre
 * away across a table.
 *
 * It is the only style here that is not imitating a drawing, and that shows in
 * what it LEAVES OUT rather than what it adds: no paper grain, no ocean hatch,
 * no relief shading. Every one of those is a mark made by a hand, and a printed
 * board has none.
 */
const BOARDGAME_PALETTE: StylePalette = {
  paper: '#f2e9d5',
  ink: '#1a1a18',
  inkLight: '#8a8272',
  // A flat printed blue. Deliberately duller than the territory colours: on a
  // board the sea is the space between the pieces, not a region you play in.
  sea: '#8fb0c4',
  seaHatch: '#8fb0c4',
  coast: '#141414',
  // Present but unused — nothing in this style's pass list shades. Kept honest
  // rather than left as a copy of another style's tints.
  shadow: '#141414',
  highlight: '#ffffff',
  ice: '#f4fbfd',
  // A table, literally: this is the one style whose subject is a thing that
  // sits on one. Warm mid-dark, lighter than the other desks so the board reads
  // as an object placed on wood rather than as a page in a dark room.
  desk: '#3a2c1e',
  deskEdge: '#150f09',
  deskShadow: '#0f0a06',
};

/**
 * Everything categorical or terrain-like takes the muted flat fill; only the
 * continuous ramps keep their own.
 *
 * `bare` is deliberately absent. Bare paper is a drawing technique — it asks
 * the reader to infer terrain from glyphs and outlines — and a board prints
 * every region in a colour so you can see at a glance who holds it. The 45%
 * mute toward the cream paper is what turns the app's screen-saturated palette
 * into something that looks printed rather than lit.
 */
const BARE_ON_OTHER_STYLES: ViewMode[] = ['satellite', 'biome', 'height_bw'];
const CATEGORICAL: ViewMode[] = ['political', 'province', 'culture', 'religion', 'plates'];

const boardgameFillPolicy = (mode: ViewMode): FillPolicy => {
  if (BARE_ON_OTHER_STYLES.includes(mode) || CATEGORICAL.includes(mode)) return 'categorical';
  return 'ramp';
};

export const styleBoardgame: MapStyle = {
  id: 'boardgame',
  name: 'Boardgame',
  palette: BOARDGAME_PALETTE,
  labelTheme: BOARDGAME_LABEL_THEME,
  overlayInk: BOARDGAME_OVERLAY_INK,
  fillPolicy: boardgameFillPolicy,
  passes: [
    // Grain 0: a printed board has no tooth. This is the reason paperPass takes
    // the opacity at all.
    paperPass(BOARDGAME_PALETTE, 'boardgame', 0),
    oceanFillPass(BOARDGAME_PALETTE),
    landPass(boardgameFillPolicy, BOARDGAME_PALETTE),
    // Double weight. The outline is the whole look, and on a board it is also
    // the thing you physically move pieces across.
    coastlinePass(BOARDGAME_PALETTE, 2),
    glyphPass(BOARDGAME_PALETTE),
  ],
};
