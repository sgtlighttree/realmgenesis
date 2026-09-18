import { describe, it, expect } from 'vitest';
import { legendFor, getCellColor, TEMP_MIN_C, TEMP_MAX_C, BIOME_COLORS } from '../utils/colors';
import { formatElevation } from '../utils/datum';
import { BiomeType, Cell, ViewMode, WorldData } from '../types';

// F5: the legend panel must follow the LAYER dropdown, keyed to the same colour
// source as the map fill. `legendFor(viewMode, world)` is that single source —
// the legend and `getCellColor` cannot drift because ramp stops are SAMPLED
// from `getCellColor`. Synthetic cells only; never generateWorld (keeps this in
// the seconds tier, off the workstation CPU — see CLAUDE.md testing notes).

const SEA = 0.55;

const world = (over: Partial<WorldData> = {}): WorldData =>
  ({
    cells: [],
    params: { seaLevel: SEA, maxElevationM: 9000 } as WorldData['params'],
    geoJson: { type: 'FeatureCollection', features: [] } as unknown as WorldData['geoJson'],
    ...over,
  }) as WorldData;

const cell = (over: Partial<Cell>): Cell => ({
  id: 0,
  center: { x: 0, y: 0, z: 0 },
  vertices: [],
  neighbors: [],
  height: 0.7,
  plateId: 0,
  temperature: 15,
  moisture: 0.5,
  biome: BiomeType.GRASSLAND,
  ...over,
});

const hex = (mode: ViewMode, c: Cell) => `#${getCellColor(c, mode, { seaLevel: SEA }).getHexString()}`;

// The full ViewMode union. A new mode added to types.ts without a legendFor
// branch trips the guard below — the paramLiveness analogue for legends.
const ALL_MODES: ViewMode[] = [
  'biome', 'height', 'height_bw', 'temperature', 'moisture', 'plates',
  'political', 'population', 'province', 'satellite', 'culture', 'religion',
];

describe('legendFor — coverage guard', () => {
  it('returns a legend or an explicit null for every ViewMode', () => {
    for (const mode of ALL_MODES) {
      const legend = legendFor(mode, world());
      if (legend === null) continue;
      expect(legend.title, `${mode} legend needs a title`).toBeTruthy();
      expect(legend.kind === 'ramp' || legend.kind === 'swatches').toBe(true);
    }
  });

  it('hides the panel (null) on satellite and plates', () => {
    expect(legendFor('satellite', world())).toBeNull();
    expect(legendFor('plates', world())).toBeNull();
  });
});

describe('legendFor — ramp modes sample getCellColor (no drift)', () => {
  it('temperature ramp endpoints equal the fill at the domain edges', () => {
    const legend = legendFor('temperature', world());
    expect(legend?.kind).toBe('ramp');
    if (legend?.kind !== 'ramp') return;
    const first = legend.stops[0];
    const last = legend.stops[legend.stops.length - 1];
    expect(first.pos).toBe(0);
    expect(last.pos).toBe(1);
    expect(first.color).toBe(hex('temperature', cell({ temperature: TEMP_MIN_C })));
    expect(last.color).toBe(hex('temperature', cell({ temperature: TEMP_MAX_C })));
  });

  it('temperature labels carry the real domain in °C', () => {
    const legend = legendFor('temperature', world());
    if (legend?.kind !== 'ramp') throw new Error('expected ramp');
    expect(legend.minLabel).toContain(String(TEMP_MIN_C));
    expect(legend.maxLabel).toContain(String(TEMP_MAX_C));
  });

  it('moisture ramp lists ocean as an out-of-band swatch, not on the bar', () => {
    const legend = legendFor('moisture', world());
    if (legend?.kind !== 'ramp') throw new Error('expected ramp');
    // Ocean short-circuits before the ramp in getCellColor; it must appear as an
    // OOB swatch (it covers most of the map) rather than being invisible.
    const ocean = legend.oob?.find((s) => /ocean|water/i.test(s.label));
    expect(ocean).toBeDefined();
    expect(ocean?.color).toBe(hex('moisture', cell({ height: 0.2 })));
  });

  it('population ramp lists ocean and uninhabited as out-of-band swatches', () => {
    const legend = legendFor('population', world());
    if (legend?.kind !== 'ramp') throw new Error('expected ramp');
    const labels = (legend.oob ?? []).map((s) => s.label.toLowerCase()).join('|');
    expect(labels).toMatch(/ocean|water/);
    expect(labels).toMatch(/uninhabited|none|empty/);
  });

  it('height ramp spans the full field with a sea-level tick in real metres', () => {
    const legend = legendFor('height', world({ params: { seaLevel: SEA, maxElevationM: 9000 } as WorldData['params'] }));
    if (legend?.kind !== 'ramp') throw new Error('expected ramp');
    expect(legend.stops[0].pos).toBe(0);
    expect(legend.stops[legend.stops.length - 1].pos).toBe(1);
    expect(legend.tick).toBeDefined();
    expect(legend.tick?.pos).toBeCloseTo(SEA, 5);
    // Elevation labels go through datum.ts (frac^2 curve), matching the Inspector.
    expect(legend.maxLabel).toBe(formatElevation(1, SEA, 9000));
  });
});

describe('legendFor — categorical modes', () => {
  const civWorld = world({
    civData: {
      factions: [
        { id: 0, name: 'Aurelia', color: '#ff0000', capitalId: 0, provinces: [], totalPopulation: 0 },
        { id: 1, name: 'Borealis', color: '#00ff00', capitalId: 1, provinces: [], totalPopulation: 0 },
      ],
      settlements: [],
    },
    cultures: [{ id: 0, name: 'Highlanders', color: '#0000ff', nameStyle: 'fantasy', homeCellId: 0 }],
    religions: [{ id: 0, name: 'The Old Faith', kind: 'folk', color: '#ffff00', cultureId: 0, holyCellId: null }],
  });

  it('biome legend lists all biome colours plus a sea-ice swatch', () => {
    const legend = legendFor('biome', world());
    if (legend?.kind !== 'swatches') throw new Error('expected swatches');
    // Every real biome, plus sea-ice which paints ocean white at cold seasons
    // yet has no BIOME_COLORS entry.
    expect(legend.items.length).toBe(Object.keys(BIOME_COLORS).length + 1);
    expect(legend.items.some((i) => /sea.?ice/i.test(i.label))).toBe(true);
  });

  it('political legend lists faction names, colours, and an Unclaimed swatch', () => {
    const legend = legendFor('political', civWorld);
    if (legend?.kind !== 'swatches') throw new Error('expected swatches');
    const byLabel = new Map(legend.items.map((i) => [i.label, i.color]));
    expect(byLabel.get('Aurelia')).toBe('#ff0000');
    expect(byLabel.get('Borealis')).toBe('#00ff00');
    expect(legend.items.some((i) => /unclaimed/i.test(i.label))).toBe(true);
  });

  it('province legend is faction-keyed (variants are derived shades)', () => {
    const legend = legendFor('province', civWorld);
    if (legend?.kind !== 'swatches') throw new Error('expected swatches');
    expect(legend.items.some((i) => i.label === 'Aurelia')).toBe(true);
  });

  it('culture and religion legends list their own names and colours', () => {
    const culture = legendFor('culture', civWorld);
    const religion = legendFor('religion', civWorld);
    if (culture?.kind !== 'swatches' || religion?.kind !== 'swatches') throw new Error('expected swatches');
    expect(culture.items.some((i) => i.label === 'Highlanders' && i.color === '#0000ff')).toBe(true);
    expect(religion.items.some((i) => i.label === 'The Old Faith' && i.color === '#ffff00')).toBe(true);
  });
});
