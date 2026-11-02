# V2 findings from the tcotg map run (through 2026-10-11)

Verified in tcotg unless marked otherwise. Full logs: tcotg `meta/map-run.md` (rounds 1–25) and `tech-docs/session-logs/`.

## Workflow that worked

**Draw intent → wild pass → bake → tight v3 → lock.**
1. The user sketches rough outlines.
2. A noise-led "wild" pass invents realistic coast.
3. The user picks a variant, and it is **baked** back into editable intent.
4. A tight v3 pass (wander ~15 km) follows it from then on.
5. Locking = a git tag + snapshots. Later changes need explicit sign-off, and every regen is checked byte-for-byte against the locked coast.

## Generator lessons

- **Wild passes are blob-shaped.** A noise-dominant field can't hold a tight archipelago (thin islands fuse or vanish). Two fixes:
  - **`protect_small_km2`:** intents under a size keep their drawn coast. They need a **sea feather** (else they merge with wild neighbours), and they must be **removed from the wild prior** (else their warped copy survives as a "ghost" islet beside them).
  - **`protect_only`:** for archipelagos, every intent is protected and all land comes from protection. The wild field is mirrored to sea, which kills the open-sea islands and the slivers in the feather band.
- **Coast detail on protected shapes:** displace the signed distance field with fBm. Scale amplitude by island size (`sqrt(area/ref)`, clamped; the nearest zone's scale wins). Add a **channel guard**: no land growth near the midline between two islands.
- **fBm is normalised:** the effective displacement is ~0.25 × the nominal amplitude. Settings that look large on paper are invisible on the map. Tune by measured coast deviation, not by the parameter.
- **Area solve:** a multiplicative sea-level step oscillates forever on flat lowlands. Bisect once the target is bracketed, and keep the closest iteration if it never converges.
- **Never extend the global octave stacks** for one region: it silently changes every group's output. Close-ups get their own fine stack, seeded from the world seed.
- **Point features need local resolution.** A 0.1° climate grid (~9 km) can't tell the Bastion (4 km from the peak) from the Threshold (7 km). Derive them locally (close-up) or give named features an analytic profile. The Spire uses `height · (1 − d/base)^p` (ease-in-quart) plus an ice radius.
- **Window edges:** land cut by the window edge traces open contours. Close the mask with a one-cell sea frame for contouring, and skip frame segments when drawing.
- **Sphere-correct moves:**
  - A longitude shift is a rotation about the polar axis.
  - A latitude move is a rotation about the equatorial axis perpendicular to the group's central meridian. Watch the sign: the first attempt went south.
  - A plain +lat shift stretches shapes near the pole.
- **Determinism:** outputs were byte-identical across Python 3.14/3.15 and numpy 2.5.1/2.5.4 (Asahi v3, NE wild). Pin the interpreter path (`/opt/homebrew/opt/python@3.14`), never bare `python3`.

## Gates (and how they lied)

- **A protect gate's max deviation is unreliable on thin rings** when intent points are sampled by raster cell (sparse). Use per-component deviation instead (generated components that overlap each intent).
- **The lost-intent check passes falsely for intents nested inside other intents.** Two stale islets drawn inside a reshaped main island looked like a "merge". V2 needs an *intent-inside-intent* check.
- **Lakes drawn over land** make "pin on land" checks pass for pins in lakes. Add an explicit pins-in-lakes gate.
- **Cross-group contact:** after any regen, intersect each group's generated rings with every other group's. Distance helpers that cap at their search pad report the cap, not the true gap.
- **Bearing checks need the canon frame.** One region was rotated 50° to fit a coast, so canon directions ("west into the hills") map to rotated map bearings. Store a per-region canon frame, or the checks fail on correct geography.

## Data model notes

- **Lakes** = `region` features of kind `lake`, drawn as water over land. This keeps land single-ring; holes in land rings would have broken the schema everywhere.
- **Line features don't follow pin moves.** A sketch import moved pins and stranded the river and road on their old anchors. V2: lines should reference pin ids at their endpoints, or move with their region.
- **Political borders are outputs** (claims → nearest settlement). Re-run them after any land change.

## Open needs for V2

- **An engine-facing semantic layer.** A narrative engine can't use raw coordinates, GeoJSON or images. It needs a gazetteer: each place's terrain, elevation and climate, its neighbours, distances and travel times by mode (walk, horse, cart, drake), routes (rivers, roads) and canon-frame bearings. It should be generated from the map, not hand-kept.
- **Tile outputs** (world master + regional DEMs) as in `ROADMAP-V2.md`. The close-up renderer is the prototype for regional detail.
