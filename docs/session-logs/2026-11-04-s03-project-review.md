# s03: project-wide review (2026-11-04, local M1)

Context: s02 wrap-up (`2026-11-03-s02-local-m1.md`, "Wrap-up"). Matt unhappy with terrain *detail*; asked for a broad review before more building.

## Matt's answers to the 12 grilling questions (near-verbatim)

1. **Reference look:** "Basically a Photoreal Earth. I don't expect VFX/Hollywood grade levels of detail for terrain but I expect very naturalistic terrain but only at a lower resolution."
2. **Zoom range:** infinite zoom is "honestly just a habit. We can cap it to a reasonable limit. Your call, then I'll let you know where to push the floor."
3. **Product:** "Worldbuilding tool, for sure. Basically an all in one app for generating a fantasy world", hence the roadmap's civilizations and fantasy map.
4. **Science vs look:** "Science first, but also malleable." Wants to **import vector/raster coastlines** (e.g. from other fantasy map generators), align/transform/scale them on the sphere, and have the generator adapt.
5. **Interaction:** drawing isn't necessarily core. "Editing coastlines and regenerating terrain would be good."
6. **Finite resolution:** raster world (~5 km everywhere, sharp only in baked regions) is acceptable: "Also makes optimization easier."
7. **Stack:** offered the s01 cloud transcript (he asked Rust vs Java for cross-platform, forgot to ask about Python). Would like to save disk by dropping the Rust tooling, but sees it's a major rewrite.
8. **Web/mobile:** only mattered when the project was architected in the cloud session. **Native first for now.**
9. **External tools:** exports are for **Blender** (planetary rendering), **Illustrator** (vector maps) and **game engines** (interactive). (Didn't directly answer "Gaea/WM does regions"; my reading below.)
10. **Keep:** "the bones are there for the terrain generator, it just doesn't look right when we add more detail." Will flag other things as they come up.
11. **Waits:** overnight is overkill; an Ultra reroll on the M1 is fine as is (~40 s measured in s02).
12. **Usable:** "when I feel like it is." Terrain first: it's the first thing he sees and it feeds everything downstream (civs, borders, cultures, cities, à la Azgaar's).

## Findings this session

- **Disk (n=1, 2026-11-04):** 38 GB free on `/` (short of the 40 GB ideal). Rust footprint: `~/.rustup` 1.4 GB, `~/.cargo` 188 MB, `target/` 1.9 GB (≈3.5 GB total; `target/` is reclaimable any time). Untracked in the repo: `node_modules/` **372 MB** (contents look like V1's TS deps: `@dimforge`, `@google`, eslint; probably left over from a `main` checkout, unconfirmed), `spike/` 56 MB (s02 spike PNGs), `heightmap-1.png/json` 8.4 MB (an app heightmap export?). Asked Matt before touching them.
- **Macro grid confirmed:** `rg-tectonics/src/grid.rs`: cube-sphere nodes, 8-connected. So the macro planform's straight runs/45° turns are D8-on-cube-sphere, as s02 concluded.
- Committed the lever-round spike code (`02213a0`, was uncommitted at the s02 wrap).

## Decisions / readings (pending Matt's confirm where marked)

- **Q7 → keep Rust (recommendation).** Rewriting tectonics, erosion, lakes, climate and the app to save ~2 GB net is not worth it: a Python + numba + node/three.js stack is itself ~1–1.5 GB. The s01 transcript isn't needed: Rust was a deliberate choice over Java. Disk lever instead: `cargo clean` / prune `target/`, quarantine the stray `node_modules/`.
- **Q8 → native-first relaxes two settled rules:** "GPU-first because the web has no threads" and "CPU fallbacks single-threaded on wasm32; rayon native-only". The web-artifact preview loop becomes optional. CLAUDE.md edits **proposed to Matt, not made**.
- **Q9 reading (confirm):** regional detail stays in-house (all-in-one tool); Blender/Illustrator/engines are export *consumers*, not the regional generator.
- **Q2 → zoom floor proposal:** fixed **20 km camera altitude** to start (finest stored data will be ~150–300 m baked regions, ~5 km elsewhere; at 20 km a 150 m texel spans ~10 screen px, still readable with per-pixel normals). Matt pushes from there.
- **Q4/Q5 is the architecture-changing answer:** today coastlines *emerge* from tectonics; Matt wants terrain *conditioned on a land mask he supplies or edits*. Direction doc: `docs/design/v2-direction-s03.md`.

## Fresh-context review of the direction: go with changes

All points adopted; detail in `docs/design/v2-direction-s03.md` → "Fresh-context review".
- **Correction:** the macro is already Rho8 (`solver.rs:121`), not plain D8, so "D8 griddiness" was imprecise. The planform is still bound to the 8 lattice neighbours.
- **Refuted (my proposal): a flat 20 km zoom floor.** At 20 km a 5 km texel is about 300 px. New rule: one texel spans at most ~8 px, which gives ~800 km over unbaked land and ~25 km over 150 m baked regions.
- **Order changed** to R2 (raster renderer on today's macro + river overlay) → R3 design → R1 jitter-only experiment → R4 → R5. Rejected: R1 first. R2 doesn't depend on jitter, and a 3D view is the better judge.

## Matt's go (s03)

- **Go on R2** (raster renderer on today's macro).
- **Q9 confirmed:** "RG becomes a fantasy oriented, user friendly, but less feature rich FOSS alternative to Gaea." Feature parity with Gaea/Terragen is long-term.
- **Native-first:** CLAUDE.md "Code rules" and "Preview loop" updated; HANDOFF settled decisions updated (native-first, Rust stays, product line, zoom floor).
- **Disk:** `node_modules/` (V1 leftovers) moved to `~/.Trash/node_modules-realmgenesis-v1-20261104`; `cargo clean` freed 1.9 GB (`target/probe-cache` worlds included; they rebuild on the next probe run). Free space now **40 GB**.
- Side note (side agent, agrees with the s02/s03 finding): raster rendering alone won't fix the griddy planform. That's R1's job. Saved `.rgworld` files store committed fields, so an old seed-1 save should still load after R1 (unverified; schema permitting), but regenerating seed 1 won't reproduce it.

## R2 design (raster renderer)

- Draft: `docs/design/raster-renderer.md`. The key finding from reading the code: the quadtree, camera-relative instances, skirts, geomorph and exaggeration already exist (`rg-render/src/terrain.rs`, `rg-terrain/src/tile.rs`), so R2 only changes **where heights come from**: a baked 6×2049² `DisplayRaster` (R32F height + Rg16F gradient) on the GPU instead of per-vertex CPU sampling.
- **Expectation flag for Matt:** at the same data, R2 looks about like today's macro-only Ultra (which he already called "pretty decent" in s02). Its payoff is no stalls, per-pixel normals at every zoom, and the plumbing for R4/R5. Visible *detail* needs R1/R4.
- Fresh-context review: pending.
- **R2 review: go with changes.** Adopted a simpler design: heights stay per-vertex on the CPU but come from the `DisplayRaster` (cheap Catmull-Rom), plus a fragment-only baked gradient texture for per-pixel normals, plus one surface for tiles, camera and hover, plus a level cap at L6 and the zoom floor on `distance_m` using the raster texel.
  - **Rejected:** vertex-shader displacement for now. It needs the bind-group changes, a skirt-parity workaround, a corrected (odd, odd) diagonal geomorph, and a min/max pyramid for the culling/split errors, all for no visible gain at a max level of 6.
- Matt: "clear them" → `spike/` and `heightmap-1.*` moved to `~/.Trash/realmgenesis-s02-leftovers-20261104/`.
