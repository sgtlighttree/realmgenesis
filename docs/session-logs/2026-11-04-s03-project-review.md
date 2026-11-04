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
