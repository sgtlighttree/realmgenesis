# WorldParams Reference

Every generation parameter, its default, and what it does. Type +
inline-comment source: `types.ts`. Defaults: `DEFAULT_PARAMS` in
`hooks/useWorldEngine.ts` (mirrored in `tests/helpers.ts`). Sliders:
`components/Controls.tsx`.

> **Every tunable key here must influence output** — `tests/paramLiveness.test.ts`
> fails otherwise. When you add a param, add it to that test. When a param reads
> "dead," check whether the default test world simply can't exercise it before
> assuming a bug (see [testing.md](testing.md)).

## System
| Param | Default | Range | Effect |
|-------|---------|-------|--------|
| `mapName` | `'map'` | string | Display name + export filename. |
| `seed` | `'realmgenesis'` | string | Terrain RNG seed (hashed to uint32). Drives all terrain-derived names. |
| `points` | `5000` | 2,000–200,000 | Voronoi cell count. Higher = more detail, slower. UI + import share the cap. |
| `planetRadius` | `6371` | km | **Display only** — feeds the geodesic ruler/scale; no simulation effect. |
| `axialTilt` | `23.5` | −90–90° | **D1:** obliquity — the amplitude of the seasonal declination `δ(s)=tilt·sin(2πs)`. Enters generation via the orbit-averaged annual-mean temperature (Jensen on the quadratic latitude curve), no longer a static climate offset. Also drives the per-season excursion in the render layer. |

## Geography
| Param | Default | Range | Effect |
|-------|---------|-------|--------|
| `landStyle` | `'Continents'` | Continents / Archipelago / Islands / Pangea / Custom | Terrain preset. Moving an advanced slider sets it to `Custom`. |
| `cellJitter` | `0.5` | 0–1 | Randomizes Fibonacci-sphere points; 0 = regular grid. |
| `seaLevel` | `0.55` | 0–1 | Height threshold separating ocean from land. Hard threshold for biomes. |
| `plates` | `12` | 2–30 | Number of tectonic plates. |

## Terrain (advanced)
| Param | Default | Range | Effect |
|-------|---------|-------|--------|
| `noiseScale` | `0.4` | 0.1–5.0 | FBM feature frequency; lower = broader continents. |
| `roughness` | `0.5` | 0–1 | Display-resolution structural relief amplitude (centered so 0.5 = ×1.0). |
| `detailLevel` | `3` | 1–6 | FBM octave count for structural noise. |
| `ridgeBlend` | `0.1` | 0–1 | 0 = smooth FBM, 1 = sharp ridged mountains. |
| `warpStrength` | `0.5` | 0–2 | Domain-warp intensity for organic shapes; also warps plate-boundary sampling. |
| `mountainHeight` | `1.0` | 0.5–2.0 | Stage-9b power-curve on heights **above** `seaLevel`; >1 = taller peaks. |
| `oceanDepth` | `1.0` | 0.5–2.0 | Stage-9b power-curve on depths **below** `seaLevel` — reshapes trench-vs-shelf **contrast**. |
| `seafloorDepth` | `1.0` | 0.3–2.0 | Stage-9b **linear** multiplier on depth below `seaLevel` — shifts **mean** water depth, relative shape preserved. Coastline fixed. Complements `oceanDepth`. |
| `maskType` | `'None'` | None / Pangea | Optional supercontinent height mask. |
| `erosionIterations` | `2` | 0–50 | Hydraulic + thermal erosion pass count. |

> `mountainHeight`/`oceanDepth`/`seafloorDepth` all live in the same Stage-9b
> remap, applied after final normalization, before climate. `seafloorDepth`
> replaced the former `seafloorDetail` texture knob. See
> [ENGINEERING-NOTES.md](ENGINEERING-NOTES.md).

## V3 tectonics
| Param | Default | Range | Effect |
|-------|---------|-------|--------|
| `tectonicStrength` | `0.5` | 0–2 | How strongly tectonics deform crust (uplift/rift magnitude). Renamed from `plateInfluence`; **no clamp**. |
| `simulationResolution` | `10000` | 5,000–20,000 | Macro-cell count for the tectonic sim (projected onto display cells). |
| `numTimesteps` | `20` | 10–60 | Kinematic simulation timesteps. |
| `marginCoupling` | `0.3` | 0–1 | Geometric correlation between mountain belts and continental margins. |
| `plateJitter` | `1.5` | 0–3 | Irregularity of plate seed placement → size/position variety. **Primary de-blob lever.** |
| `boundaryRoughness` | `1.5` | 0–3 | Jaggedness of plate boundaries (0 = straight arcs). **Primary de-blob lever.** |
| `spreadRate` | `0.008` | 0.004–0.02 | Chord-per-Ma seafloor spreading rate; smaller = older/deeper floor (GDH1). |
| `microplateIntensity` | `0.35` | 0–1 | Shear-driven microplates injected along high-strain boundaries. 0 = none (plate layout byte-identical). |
| `plateElongation` | `0.4` | 0–1 | Seed-chain length → mild plate elongation. **Near-inert at the macro silhouette** — kept cheap; see [ENGINEERING-NOTES.md](ENGINEERING-NOTES.md). |

## Climate
| Param | Default | Range | Effect |
|-------|---------|-------|--------|
| `baseTemperature` | `30` | −30–60°C | Equatorial temperature before elevation adjustment. |
| `poleTemperature` | `−30` | −60–10°C | Polar temperature. |
| `rainfallMultiplier` | `1.0` | 0.1–3.0 | Scales moisture globally. |
| `moistureTransport` | `0.5` | 0–1 | How far wind carries moisture inland. |
| `temperatureVariance` | `5` | 0–20 | Simplex noise added to temperature. |
| `season` | `0.5` | 0–1 | **D1: render-only** — orbital position; 0.5 = neutral (shows the canonical annual-mean world). Off-neutral shifts shown temperature, snow line, and biome edges (and D3 sea-ice extent) via a closed-form excursion in the color path. **Does not regenerate**, is **excluded from paramLiveness**, and is synced into `world.params` by an effect in `useWorldEngine`. |
| `starClass` | `'G'` | O / B / A / F / G / K / M | **D5: generation param** — host star spectral class. Scales global insolation → temperature in Kelvin (Stefan-Boltzmann), cascading into biomes + D3 sea-ice. G = 1.0 exact no-op. Live in paramLiveness; **regenerates** on change (in the auto-update dep list). |

> **season vs starClass — the render-only / generation split.** `season` never
> re-runs generation (it recolors live) and is kept out of paramLiveness on
> purpose. `starClass` *is* a generation param (changing it re-runs the world) and
> is asserted live in paramLiveness. Both sit in `WorldParams`; do not move
> `season` into the regen path or `starClass` out of it.

## Political
| Param | Default | Range | Effect |
|-------|---------|-------|--------|
| `numFactions` | `6` | 1–20 | Number of factions. |
| `civSeed` | `'realmgenesis_civs'` | string | Separate RNG seed for civ placement — re-roll civs without renaming geography. |
| `borderRoughness` | `0.2` | 0–1 | Noise on Dijkstra costs for irregular borders. |
| `civSizeVariance` | `0.5` | 0–1 | Per-faction size-factor spread; cheaper movement wins Dijkstra races → larger factions. |
| `waterCrossingCost` | `0.8` | 0–1 | Dijkstra cost multiplier for crossing water. |
| `territorialWaters` | `0.15` | 0–1 | Max graph distance from land to claim water cells. |
| `capitalSpacing` | `0.5` | 0–1 | Min angular separation between capitals (`spacing²·4/numFactions`). |
| `provinceSize` | `0.5` | 0.1–1.0 | Province target size (small/many → large/few). |

## Culture & naming
| Param | Default | Range | Effect |
|-------|---------|-------|--------|
| `numCultures` | `4` | 2–8 | Home count for the culture layer (C1). |
| `nameStyle` | `'fantasy'` | `NameStyle` | Default offline namebase style (A2); per-culture styles override. |

## Meta
| Param | Default | Options | Effect |
|-------|---------|---------|--------|
| `loreLevel` | `1` | 1 / 2 / 3 | Gemini prompt depth: 1 = names, 2 = + provinces/towns, 3 = + backstories. |

## Validation note
`utils/export.ts` bounds-checks `tectonicStrength: [0, 2.0]` on import. Before
Session 13's follow-up this was a dead `plateInfluence` key (the V2 name); it was
renamed to the live param, so the import validator now actually guards the value.
