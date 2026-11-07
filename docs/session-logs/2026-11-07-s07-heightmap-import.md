# s07 (2026-11-07, local M1): heightmap import, review and build

Continues HANDOFF → "a fresh-context review of `docs/design/heightmap-import.md`, then build".

## Fresh review (sonnet-medium, read-only): go with changes, all adopted
- Full list in the design doc, "Revised design (s07)". Key points:
  - **Blocker:** calibration (p99 rescale, mean fit or quantile remap), `finish`'s 3,800 m soft limit, the elevation taper and refine's global detail scale all rewrite imported heights. Fix: a per-node pin weight.
  - **Guided = inverted uplift** `U = K·A^m·S` on the import's filled network: the import becomes a fixed point of the implicit step. My own read before the review agreed: 30 Myr ≫ channel response time (~0.1–0.3 Myr at 30 km cells), so "start from the import" alone is forgotten. **Rejected:** "U ∝ imported height + pin ±x%" (no pin mechanism in the solver; steady state imposes its own χ-profile).
  - Exact and Guided are one mechanism (λ, steps). Loose is deferred.
  - Heights sample a mip pyramid (aliasing breaks Normal/Ultra coherence); coverage keeps level 0 (byte-identical old worlds).
  - Streaming row decode for V2 (~500 MB → tens of MB).
- **Checked:** both Klerga PNGs are **8-bit** (V1 7616² RGBA, V2 11200×5600 RGBA), so no 16-bit decode is needed. V1 also has a 32-bit TIF (deferred).
