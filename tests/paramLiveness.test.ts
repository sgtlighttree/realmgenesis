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
  civSizeVariance: { civSizeVariance: 1.0 },
  waterCrossingCost: { waterCrossingCost: 0.1 },
  territorialWaters: { territorialWaters: 0.9 },
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

  // borderRoughness is inert at the default 4-faction/300-cell test world
  // under D7 plateElongation terrain (verified deterministic across reruns —
  // not a flake); it is live at higher faction density (verified at
  // numFactions 8/12/16). Isolate it there instead of the default seed.
  it('borderRoughness changes civ borders at binding density', async () => {
    const smooth = await generateWorld(makeParams({ numFactions: 8, borderRoughness: 0.0 }));
    const rough = await generateWorld(makeParams({ numFactions: 8, borderRoughness: 0.9 }));
    expect(civSignature(rough), 'param "borderRoughness" appears to be dead — output unchanged')
      .not.toBe(civSignature(smooth));
  }, 120000);

  // provinceSize, like capitalSpacing, only binds when factions are large enough
  // to subdivide into multiple provinces. The default 300-cell world (4 tiny
  // factions) leaves it inert, so isolate it at a higher land/faction density.
  it('provinceSize changes province subdivision at binding density', async () => {
    const small = await generateWorld(makeParams({ points: 1000, numFactions: 5, provinceSize: 0.1 }));
    const large = await generateWorld(makeParams({ points: 1000, numFactions: 5, provinceSize: 0.9 }));
    expect(civSignature(large), 'param "provinceSize" appears to be dead — output unchanged')
      .not.toBe(civSignature(small));
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

  // V3 terrain model params. V3 is the only terrain path, so these must all
  // influence the terrain signature.
  it('V3 params change the terrain signature', async () => {
    const baseline = await generateWorld(makeParams());
    const baseSig = terrainSignature(baseline);

    const v3Perturbations: Record<string, Perturbation> = {
      marginCoupling: { marginCoupling: 0.8 },
      numTimesteps: { numTimesteps: 10 },
      simulationResolution: { simulationResolution: 5000 },
      plateJitter: { plateJitter: 0.8 },
      boundaryRoughness: { boundaryRoughness: 0.8 },
      spreadRate: { spreadRate: 0.02 },
      seafloorDetail: { seafloorDetail: 1.0 },
      microplateIntensity: { microplateIntensity: 0.9 },
    };

    for (const [name, perturbation] of Object.entries(v3Perturbations)) {
      const world = await generateWorld(makeParams(perturbation));
      expect(terrainSignature(world), `V3 param "${name}" appears to be dead — output unchanged`).not.toBe(baseSig);
    }

    expect(terrainSignature(await generateWorld(makeParams({ plateElongation: 0.0 }))),
      'param "plateElongation" appears to be dead — output unchanged')
      .not.toBe(terrainSignature(await generateWorld(makeParams({ plateElongation: 1.0 }))));
  }, 120000);
});
