import { DEFAULT_LABEL_THEME } from './labelTheme';
import { DEFAULT_OVERLAY_INK } from './overlayInk';
import { MapStyle } from './types';

/**
 * The pre-A3 look, expressed as a style so the registry has a neutral member.
 * Every mode is 'ramp': each view mode paints its own fill exactly as before,
 * and no style pass draws anything. Selecting this style is a visual no-op.
 */
export const styleDefault: MapStyle = {
  id: 'default',
  name: 'Default',
  palette: {
    paper: '#000000',
    ink: '#ffffff',
    inkLight: '#999999',
    sea: '#050505',
    seaHatch: '#050505',
    coast: '#ffffff',
    shadow: '#000000',
    highlight: '#ffffff',
    ice: '#ffffff',
    // The pre-A3 look had no backdrop and this style draws nothing anyway (its
    // pass list is empty by design); these exist to satisfy the palette.
    desk: '#000000',
    deskShadow: '#000000',
  },
  labelTheme: DEFAULT_LABEL_THEME,
  overlayInk: DEFAULT_OVERLAY_INK,
  fillPolicy: () => 'ramp',
  // Empty by design: this style IS the pre-A3 look, and every view mode already
  // paints its own fill. An empty pass list is also the signal a render path
  // uses to run its legacy loop instead — see MapStyle.passes.
  passes: [],
};
