import { describe, it, expect } from 'vitest';
import { generateWorld } from '../utils/worldGen';
import { WorldParams } from '../types';
import { makeParams, terrainSignature, civSignature, nameSignature, cultureSignature } from './helpers';

// Every tunable parameter must provably influence the generated world.
// This guards against "dead slider" regressions (a UI control bound to a
// param the engine never reads — see AUDIT.md finding A3).
//
// Display-only params are explicitly allowlisted:
//   mapName        – filename/label only
//   planetRadius   – documented display-only
//   loreLevel      – only affects the Gemini prompt
//   civSeed/seed   – covered by their own cases below

type Perturbation = Partial<WorldParams>;

// Terrain-layer params: expect height/biome/temp/moisture changes
const TERRAIN_PERTURBATIONS: Record<string, Perturbation> = {
  seed: { seed: 'other_seed' },
  points: { points: 320 },
  plates: { plates: 12 },
  seaLevel: { seaLevel: 0.45 },
  roughness: { roughness: 0.9 },
  detailLevel: { detailLevel: 1 },
  landStyle: { landStyle: 'Islands' },
  cellJitter: { cellJitter: 0.1 },
  noiseScale: { noiseScale: 1.2 },
  ridgeBlend: { ridgeBlend: 0.7 },
  mountainHeight: { mountainHeight: 1.8 },
  oceanDepth: { oceanDepth: 1.8 },
  maskType: { maskType: 'Pangea' },
  warpStrength: { warpStrength: 1.5 },
  tectonicStrength: { tectonicStrength: 0.9 },
  erosionIterations: { erosionIterations: 5 },
  baseTemperature: { baseTemperature: 10 },
  poleTemperature: { poleTemperature: -10 },
  rainfallMultiplier: { rainfallMultiplier: 2.5 },
  moistureTransport: { moistureTransport: 0.9 },
  temperatureVariance: { temperatureVariance: 15 },
  axialTilt: { axialTilt: 60 },
};

// Civ-layer params: expect regionId/provinceId/population changes
const CIV_PERTURBATIONS: Record<string, Perturbation> = {
  civSeed: { civSeed: 'other_civs' },
  numFactions: { numFactions: 2 },
  borderRoughness: { borderRoughness: 0.9 },
  civSizeVariance: { civSizeVariance: 1.0 },
  waterCrossingCost: { waterCrossingCost: 0.1 },
  territorialWaters: { territorialWaters: 0.9 },
  provinceSize: { provinceSize: 0.1 },
};

describe('every tunable param influences the world', () => {
  it('terrain params change the terrain signature', async () => {
    const baseline = await generateWorld(makeParams());
    const baseSig = terrainSignature(baseline);

    for (const [name, perturbation] of Object.entries(TERRAIN_PERTURBATIONS)) {
      const world = await generateWorld(makeParams(perturbation));
      expect(terrainSignature(world), `param "${name}" appears to be dead — output unchanged`).not.toBe(baseSig);
    }
  }, 120000);

  it('civ params change the civilization signature', async () => {
    const baseline = await generateWorld(makeParams());
    const baseSig = civSignature(baseline);

    for (const [name, perturbation] of Object.entries(CIV_PERTURBATIONS)) {
      const world = await generateWorld(makeParams(perturbation));
      expect(civSignature(world), `param "${name}" appears to be dead — output unchanged`).not.toBe(baseSig);
    }
  }, 120000);

  // capitalSpacing only *binds* when capitals are dense enough for the minimum
  // separation to reject a candidate. At the default faction count under V3
  // terrain the capitals already spread past the threshold, so the constraint
  // is inert — the param is live (verified at 8/12 factions), just non-binding
  // there. Isolate it at a binding density instead of the default seed.
  it('capitalSpacing changes capital placement at binding density', async () => {
    const tight = await generateWorld(makeParams({ numFactions: 12, capitalSpacing: 0.1 }));
    const loose = await generateWorld(makeParams({ numFactions: 12, capitalSpacing: 1.0 }));
    expect(civSignature(loose), 'param "capitalSpacing" appears to be dead — output unchanged')
      .not.toBe(civSignature(tight));
  }, 120000);

  // nameStyle is a labels-only param: it must change the generated names but
  // must NOT perturb terrain or civ geometry, so it gets its own case rather
  // than joining CIV_PERTURBATIONS (which asserts civSignature *changes*).
  it('nameStyle changes names while leaving geometry untouched', async () => {
    const baseline = await generateWorld(makeParams({ nameStyle: 'fantasy' }));
    const restyled = await generateWorld(makeParams({ nameStyle: 'norse' }));

    expect(nameSignature(restyled), 'param "nameStyle" appears to be dead — names unchanged')
      .not.toBe(nameSignature(baseline));
    expect(civSignature(restyled), 'nameStyle must not alter civ geometry')
      .toBe(civSignature(baseline));
    expect(terrainSignature(restyled), 'nameStyle must not alter terrain')
      .toBe(terrainSignature(baseline));
  }, 120000);

  // numCultures is a culture-layer-only param (C1): it must reshape culture
  // assignment (and therefore, indirectly, faction/province/town names — see
  // recalculateCivs/recalculateProvinces), but recalculateCultures runs on
  // its own RNG side-stream before civRng is ever created, so it must NOT
  // perturb civ geometry (regionId/provinceId/population) or terrain at all.
  // The civSignature/terrainSignature assertions here double as the
  // enforcement of that determinism constraint.
  it('numCultures changes culture assignment while leaving civ/terrain geometry untouched', async () => {
    const baseline = await generateWorld(makeParams({ numCultures: 4 }));
    const changed = await generateWorld(makeParams({ numCultures: 8 }));

    expect(cultureSignature(changed), 'param "numCultures" appears to be dead — culture assignment unchanged')
      .not.toBe(cultureSignature(baseline));
    expect(civSignature(changed), 'numCultures must not alter civ geometry')
      .toBe(civSignature(baseline));
    expect(terrainSignature(changed), 'numCultures must not alter terrain')
      .toBe(terrainSignature(baseline));
  }, 120000);

  // V3 terrain model params are inert when V3_ENABLED is false (default).
  // Remove .skip when V3 goes live.
  it.skip('V3 params change the terrain signature when V3 is active', async () => {
    const baseline = await generateWorld(makeParams());
    const baseSig = terrainSignature(baseline);

    const v3Perturbations: Record<string, Perturbation> = {
      marginCoupling: { marginCoupling: 0.8 },
      numTimesteps: { numTimesteps: 10 },
      simulationResolution: { simulationResolution: 5000 },
      plateJitter: { plateJitter: 0.8 },
      boundaryRoughness: { boundaryRoughness: 0.8 },
    };

    for (const [name, perturbation] of Object.entries(v3Perturbations)) {
      const world = await generateWorld(makeParams(perturbation));
      expect(terrainSignature(world), `V3 param "${name}" appears to be dead — output unchanged`).not.toBe(baseSig);
    }
  }, 120000);
});
