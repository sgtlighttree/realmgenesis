#!/usr/bin/env python3
"""Bake a generated coast back into editable intent outlines (2026-10-07).

Workflow: generate a group loose/wild -> bake -> the realistic coast becomes the new intent -> generate it tight.
  tools/map/py tools/map/transforms/bake_generated.py --group wi [--lod 1] [--eps-km 3] [--min-km2 60] [--src /tmp/mapcheck/<tag>-v2wild.json] [--dry-run]

- Reads library/maps/generated/<group>.json (LOD 0 = full, 1 = ~1/4, 2 = ~1/16 vertices), simplifies each ring with
  Douglas-Peucker at --eps-km, and drops islands under --min-km2.
- Replaces every land feature of kind "intent" carrying the group tag. Each new island keeps the id/name/tags of the
  old intent whose centroid it contains (else gets a fresh id). Pins, regions, claims and genzones are untouched.
- Backs up map-state.json to library/maps/history/ first. It can't tell whether the generated file is stale relative
  to the intents, so regenerate the group right before baking. Use --dry-run to see counts first.
"""
import argparse, json, math, os, sys, tempfile, time, shutil
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tools/map/terrain"))
from coastgen import dp_tol  # noqa: E402

STATE = ROOT / "library/maps/map-state.json"


def pip(x, y, ring):
    c = False
    for i in range(len(ring)):
        x1, y1 = ring[i]; x2, y2 = ring[i - 1]
        if (y1 > y) != (y2 > y) and x < (x2 - x1) * (y - y1) / (y2 - y1) + x1: c = not c
    return c


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--group", required=True); ap.add_argument("--lod", type=int, default=1)
    ap.add_argument("--eps-km", type=float, default=3.0); ap.add_argument("--min-km2", type=float, default=60.0)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--src", default=None, help="generated file to bake (default library/maps/generated/<group>.json), e.g. a /tmp wild preview")
    a = ap.parse_args()
    s = json.loads(STATE.read_text()); R = s["planet"]["radius_km"]; DEG = R * math.pi / 180
    gen = json.loads((Path(a.src) if a.src else ROOT / "library/maps/generated" / f"{a.group}.json").read_text())
    old = [f for f in s["features"] if f["type"] == "land" and f.get("kind") == "intent" and a.group in (f.get("groups") or [])]
    cents = [(f, np.array(f["geometry"]).mean(0)) for f in old]
    new, used = [], set()
    for gf in gen["features"]:
        if gf.get("kind") == "lake": continue   # lakes bake separately (region kind "lake")
        if gf.get("area_km2", 0) < a.min_km2: continue
        ring = np.array(gf["lod"][min(a.lod, len(gf["lod"]) - 1)], float)
        lat0 = ring[:, 1].mean(); k = math.cos(math.radians(lat0))
        P = np.stack([ring[:, 0] * k * DEG, ring[:, 1] * DEG], 1)       # local km (fine for DP)
        ring = ring[dp_tol(P, a.eps_km)]
        if len(ring) < 4: continue
        geo = [[round(float(x), 4), round(float(y), 4)] for x, y in ring]
        match = next((f for f, c in cents if f["id"] not in used and pip(c[0], c[1], geo)), None)
        if match:
            used.add(match["id"]); nf = {k2: v for k2, v in match.items() if k2 != "geometry"}
        else:
            nf = {"id": f"intent-{a.group}-baked-{len(new) + 1:03d}", "type": "land", "kind": "intent", "status": "proposed",
                  "name": f"{gen.get('label', a.group)} isle {len(new) + 1}", "groups": [a.group], "cites": "", "note": ""}
        tag = f"(baked from generated coast {gen['features'][0].get('params', '')}, {time.strftime('%Y-%m-%d')})"
        if "baked from generated" not in nf.get("note", ""): nf["note"] = (nf.get("note", "") + " " + tag).strip()
        nf["geometry"] = geo; new.append(nf)
    # lakes: replace this group's lake regions with the generated ones (drawn as water over land; land rings stay single-ring)
    old_lakes = [f for f in s["features"] if f.get("type") == "region" and f.get("kind") == "lake" and a.group in (f.get("groups") or [])]
    lakes = []
    for gf in gen["features"]:
        if gf.get("kind") != "lake": continue
        ring = np.array(gf["lod"][0], float); k = math.cos(math.radians(ring[:, 1].mean()))
        ring = ring[dp_tol(np.stack([ring[:, 0] * k * DEG, ring[:, 1] * DEG], 1), min(a.eps_km, 1.5))]
        if len(ring) < 4: continue
        lakes.append({"id": f"lake-{a.group}-{len(lakes) + 1:03d}", "type": "region", "kind": "lake", "status": "proposed",
                      "name": f"{gen.get('label', a.group)} lake {len(lakes) + 1}", "groups": [a.group], "cites": "",
                      "note": f"baked from generated basins, {time.strftime('%Y-%m-%d')}", "geometry": [[round(float(x), 4), round(float(y), 4)] for x, y in ring]})
    new_lakes_msg = f", lakes {len(old_lakes)} -> {len(lakes)}"
    verts = sum(len(f["geometry"]) for f in new)
    print(f"{a.group}: {len(old)} intents -> {len(new)} baked ({len(used)} kept their id/name), {verts} vertices, eps {a.eps_km} km, LOD {a.lod}{new_lakes_msg}")
    if a.dry_run: return
    hist = STATE.parent / "history"; hist.mkdir(exist_ok=True)
    shutil.copy2(STATE, hist / f"map-state-{time.strftime('%Y%m%d-%H%M%S')}-pre-bake-{a.group}.json")
    oldids = {f["id"] for f in old} | {f["id"] for f in old_lakes}
    s["features"] = [f for f in s["features"] if f["id"] not in oldids] + new + lakes
    fd, tmp = tempfile.mkstemp(dir=STATE.parent)
    with os.fdopen(fd, "w") as fh: fh.write(json.dumps(s, indent=2, ensure_ascii=False) + "\n")
    os.replace(tmp, STATE)


if __name__ == "__main__":
    main()
