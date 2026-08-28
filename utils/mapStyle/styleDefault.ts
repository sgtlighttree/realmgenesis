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
  },
  fillPolicy: () => 'ramp',
};
