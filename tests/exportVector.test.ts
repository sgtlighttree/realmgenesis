import { describe, it, expect, beforeAll } from 'vitest';
import { generateWorld } from '../utils/worldGen';
import { exportSVG, exportGeoJSON } from '../utils/exportVector';
import { WorldData } from '../types';
import { makeParams } from './helpers';

let world: WorldData;

beforeAll(async () => {
  world = await generateWorld(makeParams());
}, 30000);

const cellsWithGeometry = (w: WorldData): number =>
  w.cells.filter((_, i) => w.geoJson.features[i]?.geometry).length;

describe('exportSVG', () => {
  it('starts with an <svg> tag and contains all five group ids', () => {
    const svg = exportSVG(world, 'biome', 'equirectangular');
    expect(svg.startsWith('<svg')).toBe(true);
    expect(svg).toContain('id="cells"');
    expect(svg).toContain('id="coastlines"');
    expect(svg).toContain('id="rivers"');
    expect(svg).toContain('id="borders"');
    expect(svg).toContain('id="labels"');
  });

  it('emits one cell <path> per cell that has geometry', () => {
    const svg = exportSVG(world, 'biome', 'equirectangular');
    const cellsGroup = svg.match(/<g id="cells"[^>]*>(.*?)<\/g>/s)?.[1] ?? '';
    const pathCount = (cellsGroup.match(/<path /g) ?? []).length;
    expect(pathCount).toBe(cellsWithGeometry(world));
    expect(pathCount).toBeGreaterThan(0);
  });

  it('includes label <text> elements', () => {
    const svg = exportSVG(world, 'political', 'equirectangular');
    expect(svg).toMatch(/<text[^>]*>[^<]+<\/text>/);
  });

  it('escapes XML-unsafe characters in label names', () => {
    const mutated: WorldData = {
      ...world,
      civData: world.civData
        ? { factions: world.civData.factions.map((f, i) => (i === 0 ? { ...f, name: 'Fire & <Ice>' } : f)) }
        : world.civData,
    };
    const svg = exportSVG(mutated, 'political', 'equirectangular');
    // Faction labels render uppercase (LABEL_CONFIG.faction.uppercase).
    expect(svg).toContain('FIRE &amp; &lt;ICE&gt;');
    expect(svg).not.toContain('Fire & <Ice>');
  });

  it('produces no NaN anywhere in the output', () => {
    const svg = exportSVG(world, 'satellite', 'equirectangular');
    expect(svg).not.toContain('NaN');
  });

  it('works across supported projections without NaN', () => {
    (['mercator', 'winkeltripel', 'robinson', 'mollweide', 'orthographic'] as const).forEach(proj => {
      const svg = exportSVG(world, 'biome', proj);
      expect(svg.startsWith('<svg')).toBe(true);
      expect(svg).not.toContain('NaN');
    });
  });
});

describe('exportGeoJSON', () => {
  it('parses as a valid FeatureCollection', () => {
    const json = exportGeoJSON(world);
    const parsed = JSON.parse(json);
    expect(parsed.type).toBe('FeatureCollection');
    expect(Array.isArray(parsed.features)).toBe(true);
  });

  it('emits one polygon feature per cell with geometry', () => {
    const parsed = JSON.parse(exportGeoJSON(world));
    const polygons = parsed.features.filter((f: { geometry: { type: string } }) => f.geometry.type === 'Polygon');
    expect(polygons.length).toBe(cellsWithGeometry(world));
  });

  it('every coordinate pair is finite and within lon/lat bounds', () => {
    const parsed = JSON.parse(exportGeoJSON(world));
    const checkCoord = (coord: unknown): void => {
      expect(Array.isArray(coord)).toBe(true);
      const c = coord as number[];
      if (typeof c[0] === 'number') {
        expect(Number.isFinite(c[0])).toBe(true);
        expect(Number.isFinite(c[1])).toBe(true);
        expect(c[0]).toBeGreaterThanOrEqual(-180);
        expect(c[0]).toBeLessThanOrEqual(180);
        expect(c[1]).toBeGreaterThanOrEqual(-90);
        expect(c[1]).toBeLessThanOrEqual(90);
      } else {
        (coord as unknown[]).forEach(checkCoord);
      }
    };
    parsed.features.forEach((f: { geometry: { coordinates: unknown } }) => checkCoord(f.geometry.coordinates));
  });

  it('includes river, border, coastline, and label features with kind properties', () => {
    const parsed = JSON.parse(exportGeoJSON(world));
    const kinds = new Set(parsed.features.map((f: { properties: { kind?: string } }) => f.properties.kind));
    expect(kinds.has('river')).toBe(true);
    expect(kinds.has('border')).toBe(true);
    expect(kinds.has('coastline')).toBe(true);
    // Label features carry their own MapLabel kind (faction/capital/town/etc),
    // not a generic 'label' kind — check by presence of a name property + Point geometry instead.
    const labelFeatures = parsed.features.filter(
      (f: { geometry: { type: string }; properties: { name?: string } }) =>
        f.geometry.type === 'Point' && typeof f.properties.name === 'string',
    );
    expect(labelFeatures.length).toBeGreaterThan(0);
  });

  it('resolves factionName for a cell that has a regionId', () => {
    const parsed = JSON.parse(exportGeoJSON(world));
    const withFaction = parsed.features.find(
      (f: { geometry: { type: string }; properties: { regionId?: number; factionName?: string } }) =>
        f.geometry.type === 'Polygon' && f.properties.regionId !== undefined,
    );
    expect(withFaction).toBeDefined();
    expect(typeof withFaction.properties.factionName).toBe('string');
    expect(withFaction.properties.factionName.length).toBeGreaterThan(0);
  });

  it('does not include a top-level crs member (RFC 7946-clean)', () => {
    const parsed = JSON.parse(exportGeoJSON(world));
    expect(parsed.crs).toBeUndefined();
  });
});

describe('route export (C3)', () => {
  // 'route-test' is the seed routes.test.ts confirms yields road routes.
  let routeWorld: Awaited<ReturnType<typeof generateWorld>>;
  beforeAll(async () => {
    routeWorld = await generateWorld(makeParams({ seed: 'route-test', civSeed: 'route-test' }));
  });

  it('SVG emits an id="routes" group with the road stroke color', () => {
    expect((routeWorld.routes ?? []).some(r => r.kind === 'road')).toBe(true);
    const svg = exportSVG(routeWorld, 'biome', 'equirectangular');
    expect(svg).toContain('id="routes"');
    expect(svg).toContain('stroke="#c8a25a"'); // road color
  });

  it('GeoJSON includes road route features with a "road" kind', () => {
    const parsed = JSON.parse(exportGeoJSON(routeWorld));
    const kinds = new Set(parsed.features.map((f: { properties: { kind?: string } }) => f.properties.kind));
    expect(kinds.has('road')).toBe(true);
  });
});
