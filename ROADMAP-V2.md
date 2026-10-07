# RealmGenesis V2: Roadmap

Started 2026-10-07. V1 (spherical Voronoi cells, see `ROADMAP.md`) is retired as the main paradigm. Its Voronoi/cell machinery may still be useful later (e.g. political or province layers), so it's kept, not deleted.

**V2 is incubating in `~/tcotg`** (the Chronicle of the Gale map run) and moves here once it has been refined there. Until then, tcotg is the source of truth; this file records the direction so it isn't lost.

## Paradigm

- **Raster-first, on the sphere.** Elevation and other layers are equirectangular rasters; all noise is sampled in 3D on the unit sphere (no seams, no pole pinch).
- **Drawn intent → generated detail.** The user draws rough coastlines (intent shapes); the generator produces realistic coast and relief that stays within a per-region **wander** distance of the drawn line (signed distance to the intent coast + clamped noise). Topology guard: no split landmasses, no sealed straits.
- **Generation zones.** Drawn polygons carrying params (wander, roughness, detail wavelength, aridity, biome override, arcs on/off), feathered and prioritised, instead of one global settings panel. Volcanic arcs are drawable line features.
- **World-coherent noise.** One world seed and one wavelength stack shared by every region; zones only scale it, never reseed it.
- **Relief separate from the coastline.** The land/sea mask comes from the coast field; mountains/hillshade/biome come from a relief field applied on land only.
- **Claims → borders.** Political borders = land within editable claim polygons; overlaps split by nearest settlement (hard border); everything else fuzzy frontier.

## Pieces (as prototyped in tcotg)

| Piece | tcotg location | V2 home |
|---|---|---|
| Coast/relief generator (Python, numpy + PIL) | `tools/map/terrain/coastgen.py` | generation backend |
| Map editor (vanilla JS, single file) | `tools/map/editor.html` | reference UX for the TS/three.js frontend |
| Local server (stdlib; guarded `/regen`) | `tools/map/serve.py` | dev server / backend API |
| Transforms (borders, merges, spreads) | `tools/map/transforms/` | library functions |
| Data contract | `tools/map/README.md` ("Generation contract v3") | API/schema docs |

## Outputs (tile-ready, standard formats)

- **World master:** 8192×4096 equirectangular (~3.45 km/px at the equator for R=4500 km); elevation as 16-bit PNG in metres plus float EXR; JSON sidecar with bounds, resolution and radius.
- **Regional detail DEMs:** e.g. 100–200 m/px over ~500 km windows, same format.
- Cut on a simple tile grid so a viewer can stream them.

## V2 milestones (after the move from tcotg)

1. Port the generator as a backend (Python service, or WASM if it must stay browser-only).
2. TS + three.js frontend: 2D editor (port the tcotg UX) and a 3D globe.
3. Tile pyramid + streaming: world master at low zoom, regional DEMs at high zoom, seam-free.
4. Whole-planet pass: base relief + bathymetry everywhere; intents raise land.
5. Climate, rivers and biomes on the raster; SVG/cartographic export.
6. Re-introduce V1 Voronoi as an optional political/province layer if useful.

## Non-goals for now

- No viewer/tile work inside tcotg; it stays a canon tool.
- No changes to V1 code from this plan.
