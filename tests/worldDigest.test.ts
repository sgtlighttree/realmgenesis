import { describe, it, expect } from 'vitest';
import { generateWorld } from '../utils/worldGen';
import { digestWorld, diffDigests } from './helpers/worldDigest';
import { makeParams } from './helpers';
import { WorldData, Cell } from '../types';

const clone = (w: WorldData): WorldData =>
  ({ ...w, cells: w.cells.map(c => ({ ...c, center: { ...c.center }, vertices: c.vertices.map(v => ({ ...v })), neighbors: [...c.neighbors] })) });

describe('worldDigest catches what the existing determinism suite cannot', () => {
  it('reports no differences for an unmodified copy', async () => {
    const w = await generateWorld(makeParams());
    expect(diffDigests(digestWorld(w), digestWorld(clone(w)))).toEqual([]);
  }, 30000);

  it('catches a Float64 -> Float32 downcast of height', async () => {
    const w = await generateWorld(makeParams());
    const m = clone(w);
    m.cells.forEach(c => { c.height = Math.fround(c.height); });
    expect(diffDigests(digestWorld(w), digestWorld(m))).toContain('cell.height');
  }, 30000);

  it('catches a 1e-12 perturbation that toFixed(6) hides', async () => {
    const w = await generateWorld(makeParams());
    const m = clone(w);
    m.cells[0].temperature += 1e-12;
    expect(diffDigests(digestWorld(w), digestWorld(m))).toContain('cell.temperature');
  }, 30000);

  it('catches undefined collapsing to 0', async () => {
    const w = await generateWorld(makeParams());
    const m = clone(w);
    m.cells.forEach(c => { if (c.regionId === undefined) (c as Cell).regionId = 0; });
    expect(diffDigests(digestWorld(w), digestWorld(m))).toContain('cell.regionId');
  }, 30000);

  it('catches a dropped neighbor link', async () => {
    const w = await generateWorld(makeParams());
    const m = clone(w);
    m.cells[0].neighbors.pop();
    expect(diffDigests(digestWorld(w), digestWorld(m))).toContain('cell.neighbors');
  }, 30000);

  it('catches a changed geoJson coordinate', async () => {
    const w = await generateWorld(makeParams());
    const m = clone(w);
    m.geoJson = JSON.parse(JSON.stringify(w.geoJson));
    const g = m.geoJson.features.find(f => f.geometry)!;
    g.geometry!.coordinates[0][0][0] += 1e-9;
    expect(diffDigests(digestWorld(w), digestWorld(m))).toContain('geoJson');
  }, 30000);

  it('catches a lost vertices ring', async () => {
    const w = await generateWorld(makeParams());
    const m = clone(w);
    m.cells[5].vertices = [];
    expect(diffDigests(digestWorld(w), digestWorld(m))).toContain('cell.vertices');
  }, 30000);
});
