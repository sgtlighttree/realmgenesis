import { describe, it, expect } from 'vitest';

import { passesGlobeLod, GLOBE_MIN_CAM_DIST, GLOBE_MAX_CAM_DIST } from '../components/overlays/tenants';
import { VILLAGE_MARKER_MAX_DIST, HAMLET_MARKER_MAX_DIST } from '../components/markerLod';
import { MapLabel, MAP2D_MAX_SCALE, MAP2D_MIN_SCALE, ZOOM_THRESHOLDS } from '../utils/labels';
import { DEFAULT_LABEL_VISIBILITY } from '../types';

const at = (kind: MapLabel['kind']): MapLabel => ({
  kind, name: 'X', position: { x: 1, y: 0, z: 0 }, priority: 1,
});

// THE BUG THIS FILE EXISTS FOR: the hamlet cutoff was 1.15 while OrbitControls
// clamps minDistance to 1.2, so no camera position could ever satisfy it —
// hamlet labels AND hamlet markers were dead on the globe by construction. A
// cutoff is only meaningful inside the camera's reachable range.
describe('globe LOD cutoffs are reachable', () => {
  const tiers: MapLabel['kind'][] = ['capital', 'city', 'province', 'town', 'village', 'hamlet'];

  it('shows every settlement tier at the closest the camera can get', () => {
    for (const kind of tiers) {
      expect(passesGlobeLod(at(kind), GLOBE_MIN_CAM_DIST), `${kind} unreachable`).toBe(true);
    }
  });

  it('culls the small tiers at the furthest the camera can get', () => {
    for (const kind of ['town', 'village', 'hamlet'] as MapLabel['kind'][]) {
      expect(passesGlobeLod(at(kind), GLOBE_MAX_CAM_DIST), `${kind} never culls`).toBe(false);
    }
    // Capitals are the one tier that must survive a fully zoomed-out globe.
    expect(passesGlobeLod(at('capital'), GLOBE_MAX_CAM_DIST)).toBe(true);
  });

  it('reveals the tiers in order as the camera closes in', () => {
    const firstVisible = (kind: MapLabel['kind']): number => {
      for (let d = GLOBE_MAX_CAM_DIST; d >= GLOBE_MIN_CAM_DIST; d -= 0.01) {
        if (passesGlobeLod(at(kind), d)) return d;
      }
      return GLOBE_MIN_CAM_DIST;
    };
    const order = tiers.map(firstVisible);
    // capital appears furthest out, hamlet nearest: strictly descending.
    for (let i = 1; i < order.length; i++) {
      expect(order[i], `${tiers[i]} does not appear after ${tiers[i - 1]}`)
        .toBeLessThan(order[i - 1] + 1e-9);
    }
  });

  // The marker and its name must appear at the same zoom, or a dot floats with
  // no label (or a label floats over nothing).
  it('keeps marker cutoffs identical to the label cutoffs', () => {
    expect(passesGlobeLod(at('village'), VILLAGE_MARKER_MAX_DIST)).toBe(true);
    expect(passesGlobeLod(at('village'), VILLAGE_MARKER_MAX_DIST + 0.01)).toBe(false);
    expect(passesGlobeLod(at('hamlet'), HAMLET_MARKER_MAX_DIST)).toBe(true);
    expect(passesGlobeLod(at('hamlet'), HAMLET_MARKER_MAX_DIST + 0.01)).toBe(false);
    expect(VILLAGE_MARKER_MAX_DIST).toBeGreaterThanOrEqual(GLOBE_MIN_CAM_DIST);
    expect(HAMLET_MARKER_MAX_DIST).toBeGreaterThanOrEqual(GLOBE_MIN_CAM_DIST);
  });
});

// The 2D map has the same failure mode: ZOOM_THRESHOLDS entries outside
// Map2D's clamp would hide a tier at every zoom the user can reach.
describe('2D zoom thresholds are reachable', () => {
  it('puts every label kind inside Map2D\'s zoom range', () => {
    for (const [kind, threshold] of Object.entries(ZOOM_THRESHOLDS)) {
      expect(threshold, `${kind} is above max zoom`).toBeLessThanOrEqual(MAP2D_MAX_SCALE);
      expect(threshold, `${kind} is below min zoom`).toBeGreaterThanOrEqual(0);
    }
  });

  it('shows factions and capitals at full zoom-out, and the small tiers only on zoom-in', () => {
    expect(ZOOM_THRESHOLDS.faction).toBeLessThanOrEqual(MAP2D_MIN_SCALE);
    expect(ZOOM_THRESHOLDS.capital).toBeLessThanOrEqual(MAP2D_MIN_SCALE);
    for (const kind of ['town', 'village', 'hamlet'] as const) {
      expect(ZOOM_THRESHOLDS[kind], `${kind} shows when zoomed out`).toBeGreaterThan(MAP2D_MIN_SCALE);
    }
  });

  it('reveals the settlement tiers in order', () => {
    const ladder = ['capital', 'city', 'town', 'village', 'hamlet'] as const;
    for (let i = 1; i < ladder.length; i++) {
      expect(ZOOM_THRESHOLDS[ladder[i]], `${ladder[i]} appears before ${ladder[i - 1]}`)
        .toBeGreaterThan(ZOOM_THRESHOLDS[ladder[i - 1]]);
    }
  });

  // The tiers the density slider actually adds must default to visible, or the
  // slider looks inert — the S39 bug this file was extended for.
  it('defaults the tiers the density slider adds to visible', () => {
    for (const key of ['towns', 'villages', 'hamlets', 'cities', 'capitals'] as const) {
      expect(DEFAULT_LABEL_VISIBILITY[key], `${key} defaults hidden`).toBe(true);
    }
  });
});
