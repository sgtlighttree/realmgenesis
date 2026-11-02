# Climate solver research (2026-11-02)

Web survey by a subagent (Haiku); digest below. **Confidence:** sources were found, but the claims were not checked against them by the orchestrator. **Licences are unverified** except where marked. Reimplement from papers and design docs; don't copy code until the licence is checked.

## Goal

A one-shot climate run on a procedural planet: seconds to about a minute, at 1° (later 0.5°). It produces fixed monthly fields (temperature, precipitation, wind, surface currents, sea ice), then a Köppen classification. Planet parameters are tunable: radius, tilt, rotation rate, stellar flux, sea level.

## Closest prior art

- **raguilar011095/planet_heightmap_generation** (JS, licence unverified). A fixed 1° grid with 12 months in a few seconds. Includes seasonal energy balance, sea ice, pressure belts (ITCZ, Hadley/Ferrel, monsoon lows), wind-driven gyres with western boundary currents, rain budget with orographic/rain-shadow terms, and Köppen. It claims about 50% Köppen match on Earth. Design doc: `plans/CLIMATE_REBUILD_REPORT.md`. **The strongest tier B template.** https://github.com/raguilar011095/planet_heightmap_generation
- **HexWorld** `climate.js` (MIT per the subagent). Wind-lane moisture advection with orographic uplift and subsidence drain. Documented pitfalls: descending air must drain the column, and small bumps must not trigger uplift. https://github.com/Artekka/hexworld
- **Tectonics.js** (Carl Davidson). Analytic temperature, orbital insolation. Issue #9 records precipitation approaches that failed. https://github.com/davidson16807/tectonics.js/issues/9
- **Island-Proc-Gen** (Rust, wgpu, pre-alpha). Same architecture class. https://github.com/AmanThuL/Island-Proc-Gen

## Physics references

- **Insolation and EBMs:** climlab (Brian Rose). The `DailyInsolation` and `EBM_seasonal` math. https://climlab.readthedocs.io/en/latest/models.html
- **Wind-driven gyres:** Sverdrup interior balance `β·v = curl(τ)/(ρH)`. Stommel adds linear bottom friction, which gives the western boundary current. https://uhslc.soest.hawaii.edu/ocn620/notes/western-boundary-currents.html and https://epic.awi.de/id/eprint/2330/1/Olb8888a.pdf
- **Orographic precipitation:** Smith & Barstad (2004), linear, FFT-based. Use as a later refinement over the wind-lane sweep. https://journals.ametsoc.org/view/journals/atsc/61/12/1520-0469_2004_061_1377_altoop_2.0.co_2.xml
- **Köppen-Geiger rules:** Peel, Finlayson & McMahon (2007). https://hess.copernicus.org/articles/11/1633/2007/
- **Tier A rules:** Geoff Eddy's Climate Cookbook. https://astrographer.wordpress.com/bookshelf/geoff-eddys-climate-cookbook/

## Rejected for one-shot use (useful for validation)

ExoPlaSim / PlaSim (GPL-2.0, runs of simulated years), SPEEDY (Earth-specific boundary data), Isca (heavy FMS build), PLASIM-GENIE (about 2 weeks for 1000 years), CLIMBER-X (GPL-3.0, Earth-specific), SpeedyWeather.jl (Julia; fast, but not portable). Koppenpasta (Worldbuilding Pasta) runs ExoPlaSim to produce Köppen maps, which makes it a good visual reference for what "right" looks like.

## Tier B pipeline (adopted, see ARCHITECTURE.md §Climate)

1. Insolation per latitude and month from obliquity, eccentricity and flux. Coriolis from rotation rate.
2. Seasonal surface energy balance: per-cell heat capacity (land/ocean), snow and ice albedo, diffusion, lapse rate.
3. Pressure belts, monsoon lows and highs, then winds from pressure gradient + Coriolis + friction.
4. Sverdrup/Stommel gyres per basin. Ekman upwelling.
5. Evaporation (Clausius-Clapeyron), vapour advection along the wind, rain-out on ascent, drain on descent.
6. Sea ice with seasonal growth and melt, feeding albedo.
7. Iterate steps 2–6 two or three times.
8. Köppen classification.

Steps 2, 4 and 5 map to compute shaders.

## Tier A (reroll preview)

Latitude bands + lapse rate + ocean moderation + a band-wind rain-shadow test, then Köppen. Runtime in milliseconds.
