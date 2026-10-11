# V2 incubator snapshot

A **read-only snapshot** of the V2 raster tooling as prototyped in `~/tcotg` (The Chronicle of the Gale map run). tcotg stays the source of truth until the move described in `../ROADMAP-V2.md`. Re-copy from tcotg rather than editing here.

- **Snapshot:** 2026-10-11, tcotg commit `7c02c70` (map locked as tag `map-lock-2026-10-11c` there).
- **What's here:** `tools/` (generator, close-up renderer, editor, server, transforms, venv launcher) and `FINDINGS.md` (what the run taught us, for the V2 design).
- **What's not here:** tcotg's canon data (`library/maps/map-state.json`, generated layers, snapshots). The scripts assume tcotg's layout (`library/maps/…` relative to the repo root, `parents[3]`), so they won't run from this folder as-is.

## Tools

| File | Role |
|---|---|
| `tools/terrain/coastgen.py` | Generator. v3 = coast follows the drawn intent within `wander_km` (world-coherent noise, topology guard, lakes, climate/biomes). v2 `--wild-group` = noise-led pass whose output gets baked back into the intent. |
| `tools/terrain/closeup.py` | Regional close-up **of the locked map**: rasterise the locked coast at ~0.2 km, add a separate fine-octave stack, gate against the locked coast. Never regenerates. |
| `tools/editor.html`, `tools/serve.py` | Vanilla-JS editor + stdlib server (rev-guarded saves, guarded `/regen`). Clean view = reading map. |
| `tools/transforms/` | `bake_generated.py` (generated → intent), `2026-10-07c_borders.py` (claims → political regions), `2026-10-11_barsdorf_north.py` (sphere-correct latitude move, as a template). |
| `tools/py`, `setup-venv.sh`, `requirements.txt` | Pinned venv launcher (numpy 2.5.4, Pillow 12.3.0 on python@3.14). |
| `tools/terrain/configs-examples/` | Wild-pass configs showing `protect_small_km2`, `protect_only`, `protect_detail`, basins, conjugate margin. |
| `tools/README.md` | Editor docs + "Generation contract v3". |
