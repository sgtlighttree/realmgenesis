# Map editor

Local worldbuilding map editor (equirectangular, lon -180..180, lat -90..90). Python stdlib + vanilla JS, no deps.

**Python:** the editor server is stdlib-only, but the generators (coastgen, bake, borders, transforms) and the editor's `/regen` need numpy + Pillow. Run everything through `tools/map/py`, which uses the project venv `tools/map/.venv` (Python 3.14 via `/opt/homebrew/opt/python@3.14`, pinned in `requirements.txt`; Matt stays on 3.14 for now, bare `python3` is unreliable) whatever bare `python3` points at. First time, or after a brew Python change: `tools/map/setup-venv.sh` (moves any old venv to `tools/map/.quarantine/`). Start the server through the launcher too: `/regen` runs coastgen with the server's own interpreter (`sys.executable`). 3.14 and 3.15 give byte-identical outputs (checked 2026-10-11 on Asahi v3 and NE wild).

    tools/map/py tools/map/serve.py --host 0.0.0.0   # always (Matt, 2026-10-06); prints the phone/LAN URL
    # then open http://127.0.0.1:8765 here, or http://<mac-ip>:8765 / http://<host>.local:8765 on the phone

Options: `--port N`, `--state <path>` (history/ goes next to it). `GET /state` returns an `X-Map-Rev` header; `POST /state` must send it back, and gets 409 (nothing written) if the file changed on disk since. The editor polls every 3 s: silently reloads external changes, or shows a Reload / Keep mine banner if you also have unsaved edits or a save hits 409.

`GET /generated/` lists `library/maps/generated/*.json` (`{files:[{name,size,rev}]}`); `GET /generated/<name>.json` serves one (read-only; strict `[A-Za-z0-9_-]` names, no symlinks, same Host guard). The editor loads all of them as a read-only layer and re-polls every 5 s.

Port busy = already running (it exits, kills nothing). Stop with Ctrl+C.
State: `library/maps/map-state.json` (autosaved ~1 s after each change, atomic write). Each overwrite first copies the old file to `library/maps/history/`; only the last 50 stay, older ones move to `history/quarantine/`.

Keys: V select, P pin, N note, G polygon, L line, E vertices, D delete, M measure, Esc cancel, Cmd/Ctrl+Z undo, Shift+Cmd+Z redo.
S marquee, A lasso (Shift adds, Alt subtracts). Left-panel "Marquee / lasso picks" checkboxes (Land, Regions, Lines, Pins, Notes; saved in localStorage) choose what they catch; Whole features = any vertex inside (Alt/Ctrl = fully inside only), Vertices = vertices across all checked features, then V to drag/scale/rotate them. Select: Shift-click multi-select; box handles scale (Shift = uniform) / rotate (equal-area frame, so scaling is area-honest at any latitude); arrows nudge 0.1 deg (Shift 1 deg). Side panel shows live km2.
Transform box: 4 corner handles + 4 side handles (N, S, E, W = one-axis scale) + rotate handle. Scaling is Photoshop-style: the anchor is the OPPOSITE corner / side midpoint (it stays put on screen); Alt/Option = scale about the box centre; Shift = uniform (corners and sides); Alt+Shift = both. Pressing or releasing Alt/Shift mid-drag switches live. The frame is the equal-area (sinusoidal) frame centred on the anchor, so 1.3x uniform = 1.69x area at any latitude. Rotation pivots about the centre. Handles win over Alt-click-delete-vertex. Side/corner hit areas are 40/32 px in touch edit mode; the edit toolbar chips "Anchor: Opposite/Center" and "Uniform" stand in for Alt/Shift.
Proportional editing (O toggles; left-panel "Proportional edit", chip "Prop" on phones): when on, move/scale/rotate of the selection also moves unselected points within radius R (km): p' = p + w(d) (T(p) - p) in the local frame, d = great-circle km to the nearest selected vertex/edge. Falloffs smooth (default), sphere, root, linear, sharp, constant. R: `[` / `]`, mouse wheel while dragging (no zoom then), number/slider in the panel, slider chip on phones; drawn as a sampled geodesic circle. Scope = the Land/Regions/Lines/Pins/Notes checkboxes (plus visible layers), never hidden/locked tags or the generated layer; "Same feature only" restricts to the selected feature(s). One undo reverts the whole drag. R/falloff/same-only are remembered in localStorage `mapPropEdit`.
Selection tools (left panel): "Densify to N km" (default 20; inserts points so no edge exceeds N km great-circle length, placed linearly in lon/lat so the drawn outline does not change) and "Simplify" (Douglas-Peucker, km tolerance, default 5, measured against the drawn straight edge; keeps >=3 / 2 vertices). They act on the selected features, or on the selected vertices (edges between two selected vertices / runs of selected vertices). Both undoable.
Vertices: click selects (Shift toggles), drag moves all selected, click edge midpoint to insert, Delete removes selected (min 3 / 2), Alt-click deletes one. Drag the diamond to move the Breach pole.

## Phone / touch

- **View / Edit toggle** (top bar, remembered per device in localStorage `mapMode`). Default is View when the viewport is <=820 px, the UA is a phone/tablet, or `pointer: coarse`; desktop defaults to Edit. View is read-only: no autosave POST can fire (guarded in `save()` itself), edit tools and handles are gone, detail fields are read-only. The 3 s poll still silently reloads external changes. Switching to View flushes pending edits first and refuses if a conflict/error is outstanding.
- **Gestures** (Pointer Events, `touch-action: none`): one finger pans, two fingers pinch-zoom about the midpoint and pan, double-tap zooms in. Tap selects a pin/note/feature; details open in a bottom sheet (peek = name + status; tap or swipe the handle to expand to note + cites).
- **Narrow layout (<=820 px)**: compact top bar (mode, Layers, List, Home, search). Layers and Features are slide-over drawers; safe-area insets respected; buttons >=44 px.
- **Search**: the filter box finds by name/kind/status/type; Enter (or tapping a result) selects and flies to it. Works on desktop too.
- **Labels**: pin/note labels keep a constant on-screen size (12 px touch, 11 px desktop); overlapping ones are hidden, selected first, then canon/ruled, then conflict, then proposed.
- **Edit on touch** (narrow only): bottom tool bar (Select, Pin, Note, Polygon, Line, Vertices, Marquee, Lasso, Measure, Undo, Redo). Chips above it stand in for modifiers: Add (Shift), Subtract (Alt), Fully inside. Long-press = Alt: on a vertex deletes it, on a selected feature removes it from the selection, before a marquee/lasso drag switches that drag to subtract. Polygon/Line show Cancel / Done buttons. Vertex/handle hit areas are >=28 px. In Select, dragging an unselected feature pans (tap selects it); dragging an already-selected one moves it.

## Groups and canon checks

- **Groups** (top-bar button; a slide-over drawer on phones): lists every tag found in `groups` with a feature count. Per tag: *show* (hides every feature carrying the tag), *lock* (not selectable/editable, clicks pass through, marquee/lasso skip it), *Select all* (selects every visible, unlocked feature with the tag, switches to Select, so move/scale/rotate/nudge act on all together; edit mode only). *Add selection to tag...* / *Remove tag from selection* use the tag box (tap a tag name to fill it; undoable, autosaved). A feature with several tags is hidden/locked if ANY of its tags is. Hidden/locked sets are per device in localStorage `mapGroupState` (never in the state file).
- **Canon checks** (top-bar button, red badge with the count of failing checks, panel updates live while dragging). Warn only: never blocks or locks anything. Click a row to zoom to the pair.

## Generated coasts (coastgen v2)

    tools/map/py tools/map/terrain/coastgen.py --preset wi --seed 23 --config tools/map/terrain/configs/wi-v2.json
    # -> library/maps/generated/wi.json  +  /tmp/mapcheck/coastgen-wi-v2.png   (reads map-state.json read-only; numpy + PIL only, ~12 s, ~0.4 GB)

Generated coasts are OUTPUT: edit the intent shapes (land kind `intent`) and pins in the editor, then regenerate; hand edits to the output are not possible (read-only layer) and a rerun replaces the file.
- **File** `generated/<group>.json` (compact, 4 dp, budget 1.5 MB, the generator shrinks the full-res tolerance if it would not fit): `{version, group, generator{seed, params_hash, target_km2, achieved_km2, lod_vertices, stats...}, features[{id,type:"land",kind:"generated",group,generator,seed,params,area_km2,bbox,lod:[full, ~1/4, ~1/16 rings]}], biome{...}}`. LODs are nested Douglas-Peucker subsets, crossing-checked. The editor picks the LOD by zoom (<4 px/deg lowest, <14 middle, else full), culls by bbox, draws it above intent shapes, never selects or edits it.
- **Quality**: full-res marching-squares contour + one Chaikin pass + tolerance DP (`simplify_eps_cells`, default 0.06 cell), no global vertex cap. The sea level is iterated until the FINAL simplified, rounded area is within `area_tol` (0.5%) of the target. `min_island_km2` plus `min_island_rules: [{bbox:[lon0,lat0,lon1,lat1], min_km2}]` (first rule whose bbox holds the island centroid wins); islands holding a land pin are always kept. Pins are checked against the generated land (land pins on land, `sea` pins in water); counts are printed.
- **Biome tint** (`biome` block in the config; preview + optional editor layer, not a climate engine): temperature = 27 - 28 sin^2(lat) - 6.5 C/km x elevation; moisture = latitude bands (wet equator, trade-wind desert ~25 deg, westerly belt ~50 deg) + `aridity` shift, x(1 - shadow_k x rain shadow behind ridges upwind) + lift_k x windward lift (wind from the east below 30 deg and above 60 deg, from the west between; `wind: "W"|"E"` forces it); `zones: [{name, ring, biome}]` override the class in a polygon (e.g. warm current: `"subtropical"`). Classes: desert, steppe, mediterranean scrub, temperate forest, subtropical/tropical forest, alpine, ice. Stored as a coarse 0.1 deg RLE raster (`rle` = [class, run] pairs, rows north to south, class 0 = sea; the editor clips it to the coastline). Layer toggle "Biome tint" (off by default).
- The editor shows land of kind `intent` as a dashed outline with a light fill under its own "Intent shapes" toggle.

## Schema (all geometry is [lon, lat] degrees)

    { "version": 1,
      "planet":  {"radius_km": 5600},
      "speeds":  {"walk_mph": 3, "horse_mph": 9, "cart_mph": 2.5, "riding_drake_mph": 180, "primeval_cruise_mph": 500},
      "breach":  {"pole": {"lat": 40, "lon": 10}, "note": "..."},
      "features": [ {"id", "type", "kind", "name", "status", "note", "cites", "geometry"} ] }

- `type`: `land` | `region` (geometry = one ring `[[lon,lat],...]`, implicitly closed) | `line` (open `[[lon,lat],...]`) | `pin` | `note` (geometry = `{"lon","lat"}`).
- `kind`: free text (island, kingdom, river, road, canal, city...).
- `status`: `canon` | `proposed` | `conflict` | `ruled`.
- `cites`: free text, e.g. `ep37 L1017`.
- `groups` (optional, any feature): `["dolma", "asahi", ...]` tags. Country/jurisdiction tags (`dolma`, `asahi`, `wi`, `barsdorf`, `bandit-lands`) go on pins, regions, lines and notes; landmass tags (`haldei`, `ne-peninsula`, `asahi-arch`, `wi`, `barsdorf`) go on land (intent) polygons. A feature can carry several.
- `constraints` (optional, top level): `[{"id", "a", "b", "target_mi", "tol_mi", "kind": "eq"|"max"|"bearing", "bearing_deg"?}]`, `a`/`b` = feature ids (pins; a region/land id uses its centroid). Great-circle miles use `planet.radius_km`. `eq`: red if `|miles - target_mi| > tol_mi`. `max`: red if `miles > target_mi + tol_mi`. `bearing`: initial bearing a->b in degrees clockwise from north, red if the wrapped difference from `bearing_deg` exceeds `tol_deg` (if present, else `tol_mi` is read as degrees). A missing endpoint id shows MISSING and counts red.
- The editor preserves unknown fields on every feature and at the top level on save (the whole document is round-tripped; only the fields you edit change).
- Breach = great circle of all points 90 deg from `breach.pole` (drawn with a cosmetic wobble).
- `speeds` feed the Measure tool's travel times; `planet.radius_km` feeds distances and areas.

## Units (2026-10-07)
Mapping tools are metric: distances in km (measure also shows mi and nm), constraints use `target_km`/`tol_km` (bearing: `tol_deg`; legacy `*_mi` still read). Speeds: ground `*_kmh`, air `*_kt` (knots); legacy `*_mph` still read. Prose stays imperial-first.

## Borders
`tools/map/transforms/2026-10-07c_borders.py` rebuilds political regions from the land: islands by landmass tag, the main body by nearest Dolma/Asahi pin within a claim radius (`CLAIM_KM`), the bandit-lands region kept. Unclaimed land stays frontier. Regions carry `edge: "fuzzy"|"hard"`: fuzzy regions draw with a soft blurred edge; where two claims meet, a hard `line` with `kind: "border"` is traced. Re-run after reshaping land. The Groups panel has per-tag colour swatches and a "border" toggle.

## Generation contract v3 (2026-10-07) — shared by coastgen, serve.py and the editor
- CLI: `tools/map/py tools/map/terrain/coastgen.py --group <tag> --state <path> [--png <path>|none]`. `<tag>` is a landmass tag (`wi`, `haldei`, `asahi-arch`, `barsdorf`, `ne-peninsula`). Reads `state.gen` + genzones + arcs + pins from the state file (read-only). Writes ONLY `library/maps/generated/<tag>.json` (+ preview). Prints exactly one JSON line on stdout: `{group, area_km2, intent_area_km2, coast_dev_p95_km, coast_dev_max_km, topology_ok, topology_issues[], pins_land_ok, pins_water_ok, pins_coastal_ok, runtime_s, peak_mb}`. Exit code 0 even when topology fails (the summary says so); it writes nothing if topology fails.
- `state.gen` = `{"defaults": {wander_km, rough, detail_wl_km, aridity, feather_km}, "groups": {"<tag>": {...overrides, "seed": int}}}`.
- Genzone = region feature `{"kind": "genzone", "groups": ["<landmass tag>"], "params": {wander_km?, rough?, detail_wl_km?, aridity?, biome?, arcs_on?}, "feather_km": number, "priority": int}`. Higher priority wins; params blend across `feather_km`.
- Arc = line feature `{"kind": "arc", "groups": ["<landmass tag>"], "params": {width_km, amp, bead_km, bead_sigma, bead_max}}` (opt-in volcanic chains; they move with the land).
- Pin flags: `"coastal": true` → must end within 3 km of the generated coast; land pins on land; pins tagged `sea` in water.
- Coast = signed distance (km) to the intent coast + noise clamped to ±wander_km (so deviation ≤ wander by construction). Relief/biome is a separate field, applied on land only.

## Move is a sphere rotation (2026-10-07)
Dragging (any group move, incl. vertices) rotates every point about the axis `start x cursor` (unit vectors) so the drag-start point lands on the cursor: true km and km² are kept, and the shape visibly stretches on the equirectangular screen when dragged poleward (width ~1/cos(lat)). Proportional editing rotates each affected point by `w * theta` about the same axis. Arrow nudge = small rotation east (about the polar axis; step is the arc step at the selection centroid, so `Delta lon = step / cos(lat)`, capped at 10 deg) or north (about the local east axis through the centroid). Scale/rotate handles are unchanged (equal-area local frame).

## Generation panel (Generation button, desktop and mobile drawer)
- `state.gen` is created in memory with defaults (`wander_km 40, rough 1, detail_wl_km 170, aridity 0, feather_km 30`; groups seeded `wi 15`, `haldei 70`, `asahi-arch 70` wander) if missing, and saves through the normal autosave.
- One block per landmass tag (tags on `land` features): sliders write per-group overrides (`*` = overridden, `reset` = inherit), a seed field, and **Regenerate <group>**. Zones of that group are listed under it.
- Genzones (region `kind: "genzone"`): selecting one shows wander_km, rough, detail_wl_km, aridity, feather_km (top-level), priority, biome (none or desert/steppe/mediterranean/temperate/subtropical/tropical/forest/alpine/ice), arcs_on in the right panel. Drawn with polygon type "genzone" (or the panel's New genzone button), auto-tagged with the landmass under it. Layer "Generation zones": dotted amber outline + `~wander km` label. Arc lines (`kind: "arc"`) draw dashed amber.
- Regenerate flushes the autosave first, then `POST /regen`; the panel polls `/regen/status`, shows area vs intent, coast deviation p95/max, topology, pin results, and reloads the generated layer when the run ends. Edit mode only.

## /regen endpoint and its safety rules (serve.py)
- `POST /regen` body `{"group": "<tag>", "rev": "<X-Map-Rev>"}` (Content-Type application/json, 4 KB max) -> 202. Runs `[sys.executable, coastgen, "--group", group, "--state", STATE]` with no shell, stdin closed, nothing else request-supplied.
- 400 bad body or group not in the whitelist (landmass tags of `land` features in the state file, `^[a-z0-9][a-z0-9-]*$`); 415 wrong content type; 409 a run is active, or `rev` != the current state-file revision (no regen on stale data); 501 coastgen missing. One run at a time, 180 s timeout (child killed).
- `GET /regen/status` -> `{running, group, started, last_summary, stderr_tail}`; `last_summary` is the last JSON line on coastgen stdout (or `{error}`).
- The endpoint never writes map-state.json; coastgen writes only `library/maps/generated/<tag>.json`. Same host check as the rest of the server.
- `COASTGEN_PATH=/path/to/stub.py tools/map/py tools/map/serve.py ...` overrides the script (for tests while coastgen is being rewritten). Note coastgen derives its output dir itself, so only use the real script against the real state file.
