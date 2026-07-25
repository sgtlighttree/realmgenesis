# Adversarial research pass — crust advection, cost, and design red-team

> **Provenance:** commissioned 2026-07-25 during the D6 brainstorm, as a
> deliberate second opinion on [the first research pass](./2026-07-25-realistic-terrain-generation.md)
> (Antigravity/Gemini). Run by a Claude `sonnet-medium` subagent with web search,
> briefed to fill gaps and to *disagree* rather than survey. Body is the agent's,
> unedited.
>
> **Why a second pass was commissioned:** the first report's treatment of crust
> advection was explicitly self-marked as *inference*, not sourced — and that one
> line ("advect crust by each plate's motion") is what the whole simulation loop
> stands on. If it is implemented as naive resampling, the crust field diffuses
> and continents dissolve into mush after ~30 timesteps.
>
> **It found three things that changed the design**, the second of which
> invalidated a decision we had already made:
>
> 1. Nobody advects crust by interpolated resampling. Implementations either
>    carry Lagrangian markers rotated exactly and reassign cells by
>    nearest-seed, or never move crust at all and classify boundaries
>    algebraically. Confirmed across five independent implementations.
> 2. **"Accumulate uplift over 20–40 timesteps" does NOT fix the faceted ridge,
>    and can sharpen it.** With small per-step rotation, the same cell-graph
>    edge is re-selected as the boundary every step, so uplift piles onto one
>    edge — a taller, thinner wall exactly on the Voronoi cut. This was our
>    headline seam fix and it was wrong.
> 3. Nobody runs a multi-step CPU simulation at 200k cells. The scaling pattern
>    is: simulate coarse (5k–20k), then project once onto the display mesh.
>
> **Caveat on confidence:** the agent could NOT access Cortial et al.'s full
> text (paywalled; HAL and LIRIS mirrors failed three times). Its Cortial-specific
> claims are inference from secondary sources. The cross-implementation pattern
> in Gap 1 is well-sourced and does not depend on Cortial.
>
> It also independently confirmed that the first pass's reference [4] is
> fabricated, and supplied correct canonical URLs.

---

## TOP 3 THINGS THAT CHANGE THE DESIGN

1. **Do not advect crust by resampling a field. Advect *markers* (rotated plate seed points / vertices), then reassign cells to their nearest rotated marker each step.** Every real implementation I could find that avoids diffusion does this — none of them interpolate a scalar field forward. This is a structural change to "advect crust by each plate's motion," not a tuning knob.

2. **"Accumulate uplift over 20–40 timesteps" does not by itself fix the faceted-ridge problem, and can make it worse.** If plate membership is a nearest-neighbor test against the same rotated generators every step, and per-step rotation is small, the same cell-graph edge gets picked as "the boundary" repeatedly — you get a taller, sharper wall on the identical edge, not a softened one. You need explicit boundary jitter/domain-warping or coarse-to-fine resampling to break the cell-graph alignment; timestep count alone doesn't do it.

3. **Simulate at a coarse resolution and project once, don't run the boundary/uplift loop at 200k cells.** No implementation I found runs a many-step CPU simulation at 200k-cell resolution. The consistent real-world pattern (World Orogen, Cortial's own follow-up paper, tectonics.js) is: simulate plates on a much smaller point/mesh set (10k–30k), then nearest-neighbor/barycentric-project the resulting crust+uplift field onto the full display mesh in a single final pass. This is cheaper, avoids the worst of the reassignment cost, and is consistent with your own "sub-cell detail is rendering only" principle — it should be "sub-cell *and* simulation-resolution detail are rendering-time concerns," not just sub-cell noise.

---

## GAP 1 — Crust advection

**Direct access to Cortial et al. 2019 full text failed** — Wiley is paywalled, the HAL and LIRIS PDF mirrors either exceeded fetch size limits or returned binary/blocked content on three attempts. I could not confirm their exact data structure from primary text. What I could confirm from the abstract/secondary sources: it operates on continental/oceanic crust as fields deformed by "plate subduction or collisions," with a companion paper describing a separate "amplification" pass for detail — meaning even Cortial's own team splits simulation resolution from final-detail resolution. ([Wiley abstract](https://onlinelibrary.wiley.com/doi/abs/10.1111/cgf.13614), [Real-time hyper-amplification of planets](https://link.springer.com/article/10.1007/s00371-020-01923-4)) — treat any Cortial-specific claim below as unconfirmed unless flagged otherwise.

**What I could confirm, robustly, from five independent working implementations**, is the general answer to your question — and it's a consistent pattern across all of them:

- **A Blender addon reimplementing Cortial's technique** moves actual mesh vertices via quaternion rotation per plate (Lagrangian, no interpolation). To avoid degeneration it does a **full remesh** periodically and **recomputes plate-membership as nearest-neighbor**, explicitly not interpolated: "I re-meshed the entire sphere and then recalculated all the custom vertex data that I need." ([Blender Artists thread](https://blenderartists.org/t/procedural-remeshing-for-tectonic-planets-addon-that-im-currently-developing/1242087))
- **World Orogen** doesn't advect a field iteratively at all: it does farthest-point seeding + flood fill once on a fixed ~20k-region mesh, classifies boundaries, propagates stress via BFS once, then **projects the result onto up to 2.56M display cells via nearest-neighbor walks**. No multi-step physics loop. ([orogen.studio](https://www.orogen.studio/), [repo](https://github.com/raguilar011095/planet_heightmap_generation))
- **Nick McDonald's "Clustered Convection"** represents each plate as a cloud of point "segments" (masses, thickness, density) moved as rigid bodies. Each timestep it **regenerates a GPU Voronoi label map from the moved points** — hard nearest-neighbor assignment, explicitly not interpolation — then applies a smoothing/cascade shader *after* the fact, so boundary identity itself never blurs. ([nickmcd.me](https://nickmcd.me/2020/12/03/clustered-convection-for-simulating-plate-tectonics/))
- **davidar.io's GPU planet sim** moves plates by **whole integer pixels per step**, specifically to sidestep sub-pixel/interpolated movement. ([davidar.io/post/sim-glsl](https://davidar.io/post/sim-glsl))
- **Red Blob Games' sphere generator** doesn't advect crust at all — plates carry only a direction vector, and boundary type/elevation is inferred algebraically from comparing neighboring plates' vectors. Crust identity never physically moves. ([redblobgames.com/x/1843-planet-generation](https://www.redblobgames.com/x/1843-planet-generation/))
- The real-geophysics analog (GPlates-adjacent literature) uses **massless Lagrangian tracer particles** advected by plate rotation, with the age/crust grid reconstructed from tracers rather than a resampled field — this is the formal name for the pattern above (particle-in-cell / tracer methods). ([arXiv:1910.03351](https://arxiv.org/pdf/1910.03351), text extraction was partial — PDF parsed as mostly binary, treat with reduced confidence, but the title and abstract-level framing is clear from the search snippet).

**Conclusion (well-supported, cross-source):** the diffusion problem you're worried about is real and is specifically caused by *interpolated resampling* (bilinear/trilinear lookups at back-rotated positions). Nobody who avoids it does so via a fancier interpolation scheme (BFECC/MacCormack) — that's atmospheric-science machinery, still interpolates each step, and would still slowly diffuse over 30 steps unless paired with real flux-conservation, which is a lot of complexity for no clear payoff here. Instead, every practical implementation either (a) carries **Lagrangian markers/vertices/points** rotated exactly (no interpolation) and reconstructs the grid via **nearest-neighbor/Voronoi reassignment**, or (b) **doesn't move crust at all** and instead classifies boundaries algebraically from plate motion vectors.

**Recommendation for your design:** keep a small array of per-plate Euler poles/quaternions and, conceptually, one seed point per plate (or a modest set of markers per plate if you want sub-plate deformation later). Each timestep, for boundary classification, inverse-rotate each cell's position into each plate's rest frame and test membership by nearest-rotated-seed — this **is** "crust fixed in a per-plate local reference frame, transform the sampling," and it's exactly what the working implementations converge on. What breaks if you do this naively: triple junctions (three plates meeting) are ambiguous under pairwise nearest-neighbor tests and will show a visible notch/singularity — worth a specific tie-break rule. Also: a plate's own crust identity is trivially preserved this way (nothing diffuses), but any *scalar* you accumulate on top (uplift, thickness) is still susceptible to diffusion if you ever blend it across a reassignment boundary — keep uplift attached to the crust identity you're tracking, not resampled independently.

## GAP 2 — Cost and timestep budget

Concrete numbers found:

- **tectonics.js** (pure JS, typed arrays, browser): "tens of thousands of grid cells with reasonable framerates," and separately, ~10,000 cells at ~20fps in Chrome for continuous real-time simulation. ([tectonics.js blog](https://davidson16807.github.io/tectonics.js/blog/))
- **World Orogen**: plate simulation itself runs on a **~20k-region coarse mesh**; only the *rendering* projection goes up to 2.56M cells, via nearest-neighbor, done once (not per-timestep). No iterative multi-step loop at high resolution. ([orogen.studio repo](https://github.com/raguilar011095/planet_heightmap_generation))
- **Clustered Convection**: 256×256 grid (~65k cells), real-time but GPU-shader-driven; author explicitly notes it's "significantly slower" without GPU. ([nickmcd.me](https://nickmcd.me/2020/12/03/clustered-convection-for-simulating-plate-tectonics/))
- **RNGalaxy** (a from-scratch reimplementation attempt, Red-Blob-Games-style): 30,000 Fibonacci points, and notably **no multi-step timestep simulation at all** — single-pass classification. ([github.com/robertoost/RNGalaxy](https://github.com/robertoost/RNGalaxy))

**No implementation I found runs a 20–40-step CPU-only simulation at 200k cells.** The nearest analog to what you're proposing (tectonics.js, continuous JS simulation) caps out around 10k cells for smooth interactivity; the ones that use higher counts either go GPU (Clustered Convection) or drop iterative simulation for a single classification pass at coarse resolution then project to fine (World Orogen — this is your best precedent).

**Order-of-magnitude reasoning (my own estimate, not cited):** if plate reassignment is done as "rotate ~10–30 plate seeds, nearest-neighbor each of 200k cells to them," that's on the order of 200k × 30 ≈ 6M distance ops per step, ×30 steps ≈ 180M ops — plausibly a few hundred ms to low seconds in a single-threaded JS worker with typed arrays (no GC pressure, simple float math). That's tolerable off the main thread. But if instead you do point-in-spherical-polygon tests against actual plate boundary polygons (not nearest-seed), cost balloons with polygon complexity and this estimate no longer holds — use nearest-rotated-seed, not polygon containment.

Given the precedent above, I'd recommend: **run the 20–40-step plate-motion/boundary/uplift loop at a coarse resolution (5k–20k macro-cells, independent of your final cell count), then do one nearest-neighbor projection onto the full 200k-cell mesh at the end**, adding structural/isostasy noise at full resolution in that single final pass. This is cheaper, matches every implementation that scaled past 10k cells, and doesn't require you to solve the O(200k × N_plates × N_steps) cost problem at all.

**Does uplift need every timestep?** No precedent runs it every step at full resolution. World Orogen computes stress once via BFS, not iteratively. If you keep the multi-step design, boundary *classification* (convergent/divergent/transform) benefits from being recomputed every step since it depends on evolving relative geometry, but *uplift accumulation* could reasonably be computed as a closed-form integral (relative velocity × dt × Nsteps, evaluated once at the coarse-resolution boundary) rather than literally looping — cutting most of the cost with no obvious visual loss, since you're already discarding sub-cell detail as "rendering only."

## GAP 3 — Red-team

- **Multi-step accumulation does not remove lattice alignment; it can sharpen it.** If plate membership is decided by nearest-rotated-seed every step and per-step rotation is small relative to cell spacing, the same cell-graph edge is re-selected as "the boundary" on most steps early on, so uplift piles onto the identical edge — you can end up with a *taller, thinner* wall precisely aligned to the Voronoi cut, which is arguably a worse version of your current faceted-ridge bug, not a fix. **Fix:** perturb the boundary test with per-cell noise/domain-warp so the effective boundary wanders across a band of cells rather than snapping to one graph edge, and/or ensure the coarse-to-fine projection (Gap 2) is what actually determines the *visual* boundary width, not the coarse simulation's cell graph.

- **Triple junctions.** Pairwise boundary classification (convergent/divergent/transform between two plates) is undefined or ambiguous where three plates meet. These points are geometrically prominent — expect visible notches, spikes, or discontinuous uplift right at plate corners unless you add an explicit triple-junction rule.

- **Two independent fields, structural incoherence:** the actual failure mode isn't "ocean crust sits under a mid-ocean ridge" (your own rule set already handles that correctly — divergent boundaries are supposed to overwrite with new oceanic crust, which is geologically right, that's literally how rifting works). The real gap is **intra-plate crust-type boundaries that have no associated plate boundary at all** — since crust type is seeded independently of plate id, you will get continental/oceanic transitions sitting in the *interior* of a single plate, far from any classified boundary, which currently has zero uplift/rift treatment in your pipeline. In reality this is a passive margin (e.g. the Atlantic coast) and should get a gentle shelf slope from isostasy alone — but if your compose step assumes all interesting relief comes from boundary uplift, these will render as unmotivated flat cliffs or invisible edges. Add an explicit passive-margin case.

- **Structural mismatch, not a bug you can fix with tuning:** because crust and plates are seeded independently, no number of timesteps makes coastlines and mountain belts *causally* related the way they are on Earth (where plates ARE defined by where crust broke). Expect mountain ranges/trenches that cut across continents in ways that look arbitrary relative to the coastline, especially at low seed counts or few timesteps, because the belts are tracking a different random field than the coastline is. This is presented as a deliberate simplification in your design — worth stating explicitly in the design doc as an accepted tradeoff, not something 20–40 steps will "cook out."

- **Isostasy specifics:** a simple linear thickness→elevation relation will make any long-lived collision belt grow **unboundedly** over 30 steps unless clamped — real crustal thickening saturates (flexural rigidity, not naive Airy isostasy) well before Himalaya-scale thickness. Clamp or use a saturating (e.g. log/asymptotic) response. Also watch for **double-counting**: if "accumulated uplift" and "isostasy from thickness" are both driven by the same convergent-collision event, you may be adding two upward terms for one physical cause — make sure uplift represents the *kinematic* forcing and isostasy represents the *static* buoyancy response to the resulting thickness, not two independent bonuses stacked on the same cells.

- **Anything else likely to look wrong:** Euler-pole rotation can produce inconsistent apparent motion across a single plate if the pole sits inside the plate's own territory (near-zero velocity near the pole, fast motion far from it) — expect uneven uplift intensity along one plate's boundary depending on rotation-pole placement, which may read as arbitrary rather than tectonically motivated.

## GAP 4 — Citation check

- **"Erleben, K. et al., Lattice-aligned artifact mitigation via Lloyd's Relaxation and Voronoi smoothing" — could not find this anywhere, and I could not find any combination of Erleben's real published work (rigid-body/physics simulation, University of Copenhagen) with terrain/Voronoi artifact mitigation.** Treat as fabricated; strike it. If you want a real citation for the underlying idea (Lloyd relaxation reducing irregular-cell artifacts), cite Lloyd's actual 1982 paper ("Least Squares Quantization in PCM," IEEE Trans. Information Theory) or, more practically, Red Blob Games' own polygon map generation writeup, which explicitly discusses Lloyd relaxation for this purpose.
- **Red Blob Games — polygonal map generation:** current canonical URL is **https://www.redblobgames.com/maps/mapgen2/** (the maintained JS rewrite); the original 2010 article lives at the older Stanford mirror (http://www-cs-students.stanford.edu/~amitp/game-programming/polygon-map-generation/, still referenced by forks on GitHub). The sphere-specific piece relevant to this project is **https://www.redblobgames.com/x/1843-planet-generation/**.
- **mapgen4:** **https://www.redblobgames.com/maps/mapgen4/** — confirmed.
- **Cortial et al. 2019, correct citation:** Y. Cortial, A. Peytavie, E. Galin, E. Guérin, "Procedural Tectonic Planets," *Computer Graphics Forum*, vol. 38, no. 2 (Eurographics 2019 proceedings), DOI [10.1111/cgf.13614](https://onlinelibrary.wiley.com/doi/abs/10.1111/cgf.13614). Free PDF mirrors: [HAL](https://hal.science/hal-02136820v1/file/2019-Procedural-Tectonic-Planets.pdf), [LIRIS](https://perso.liris.cnrs.fr/eric.galin/Articles/2019-planets.pdf) (both fetch-blocked for me during this pass — access failed 3 times via different methods; I could not independently verify their crust-advection method from primary text, only via the cross-implementation pattern in Gap 1).
