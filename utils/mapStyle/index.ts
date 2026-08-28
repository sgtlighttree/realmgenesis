import { MapStyle, MapStyleId } from './types';
import { styleDefault } from './styleDefault';
import { styleParchment } from './styleParchment';

export * from './types';

export const MAP_STYLES: Record<MapStyleId, MapStyle> = {
  default: styleDefault,
  parchment: styleParchment,
};

/** Unknown ids fall back to `default` — a saved/stale id must never blank the map. */
export const getMapStyle = (id: MapStyleId): MapStyle => MAP_STYLES[id] ?? styleDefault;
