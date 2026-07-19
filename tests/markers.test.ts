import { describe, it, expect } from 'vitest';
import { collectLabels, LABEL_CONFIG } from '../utils/labels';
import { validateMarkers } from '../utils/export';
import { DEFAULT_LABEL_VISIBILITY, MarkerData, WorldData } from '../types';
import { makeParams } from './helpers';

// Minimal WorldData fixture — collectLabels only touches cells/civData for
// faction/province/town labels, so an empty cell set is enough to isolate
// marker emission (which is independent of civData, like geo features).
const makeWorld = (markers: MarkerData[] | undefined): WorldData => ({
  cells: [],
  params: makeParams(),
  geoJson: { type: 'FeatureCollection', features: [] },
  markers,
});

const marker = (overrides: Partial<MarkerData> = {}): MarkerData => ({
  id: 0,
  kind: 'poi',
  name: 'Old Watchtower',
  note: 'Abandoned since the war.',
  position: { x: 0.6, y: 0.2, z: 0.7746 },
  ...overrides,
});

describe('collectLabels — markers', () => {
  it('emits one label per marker, independent of civData', () => {
    const world = makeWorld([marker()]);
    expect(world.civData).toBeUndefined();
    const labels = collectLabels(world);
    expect(labels).toHaveLength(1);
    expect(labels[0].kind).toBe('marker');
    expect(labels[0].name).toBe('Old Watchtower');
    expect(labels[0].priority).toBe(1.5);
    expect(labels[0].position).toEqual({ x: 0.6, y: 0.2, z: 0.7746 });
  });

  it('emits nothing when markers is absent or empty', () => {
    expect(collectLabels(makeWorld(undefined))).toHaveLength(0);
    expect(collectLabels(makeWorld([]))).toHaveLength(0);
  });

  it('emits multiple markers alongside each other', () => {
    const world = makeWorld([marker({ id: 0, name: 'A' }), marker({ id: 1, name: 'B' })]);
    const labels = collectLabels(world).filter(l => l.kind === 'marker');
    expect(labels.map(l => l.name).sort()).toEqual(['A', 'B']);
  });

  it('is gated by the markers visibility key in LABEL_CONFIG', () => {
    expect(LABEL_CONFIG.marker.visibilityKey).toBe('markers');
  });
});

describe('LabelVisibility defaults', () => {
  it('markers default to visible', () => {
    expect(DEFAULT_LABEL_VISIBILITY.markers).toBe(true);
  });
});

describe('validateMarkers', () => {
  const valid: MarkerData[] = [
    marker({ id: 0, kind: 'dungeon', name: 'The Deep', note: '' }),
    marker({ id: 1, kind: 'ruin', name: 'Old Fort', note: 'Half-collapsed.', position: { x: -0.5, y: 0.5, z: 0.7071 } }),
  ];

  it('round-trips a well-formed marker list through JSON serialization', () => {
    const serialized = JSON.parse(JSON.stringify(valid));
    expect(validateMarkers(serialized)).toEqual(valid);
  });

  it('accepts every marker kind in the union', () => {
    for (const kind of ['dungeon', 'ruin', 'battlefield', 'portal', 'poi'] as const) {
      const result = validateMarkers([marker({ kind })]);
      expect(result).toHaveLength(1);
      expect(result![0].kind).toBe(kind);
    }
  });

  it('drops entries with an unrecognized kind, keeping the rest', () => {
    const result = validateMarkers([marker({ id: 0, kind: 'castle' as never }), marker({ id: 1 })]);
    expect(result).toHaveLength(1);
    expect(result![0].id).toBe(1);
  });

  it('drops entries with a non-finite or malformed position', () => {
    const nan = validateMarkers([marker({ position: { x: NaN, y: 0, z: 0 } })]);
    expect(nan).toHaveLength(0);

    const infinite = validateMarkers([marker({ position: { x: 0, y: Infinity, z: 0 } })]);
    expect(infinite).toHaveLength(0);

    const missingAxis = validateMarkers([{ ...marker(), position: { x: 0, y: 0 } as never }]);
    expect(missingAxis).toHaveLength(0);
  });

  it('drops entries missing a string name', () => {
    const result = validateMarkers([{ ...marker(), name: undefined }]);
    expect(result).toHaveLength(0);
  });

  it('drops entries with a non-finite id', () => {
    const result = validateMarkers([marker({ id: NaN })]);
    expect(result).toHaveLength(0);
  });

  it('returns undefined (not an empty array) when the top level is not an array', () => {
    expect(validateMarkers(undefined)).toBeUndefined();
    expect(validateMarkers(null)).toBeUndefined();
    expect(validateMarkers({ 0: marker() })).toBeUndefined();
    expect(validateMarkers('markers')).toBeUndefined();
  });

  it('drops non-object entries within an otherwise-valid array', () => {
    const result = validateMarkers([marker({ id: 0 }), null, 'nope', 42]);
    expect(result).toHaveLength(1);
  });
});

// Regeneration carry-over (markers persist across handleGenerate/App.tsx
// re-generation because they're sphere-anchored, not cellId-based) is wired
// in App.tsx's handleGenerate — App-level state-management behavior with no
// pure-engine surface to unit test here.
