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
  },
  fillPolicy: () => 'ramp',
  // Empty by design: this style IS the pre-A3 look, and every view mode already
  // paints its own fill. An empty pass list is also the signal a render path
  // uses to run its legacy loop instead — see MapStyle.passes.
  passes: [],
};
