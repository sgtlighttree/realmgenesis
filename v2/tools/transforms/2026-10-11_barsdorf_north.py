#!/usr/bin/env python3
"""2026-10-11: move Barsdorf north by DLAT degrees (Matt: "push barsdorf a little north", picked 3 deg; it ran warm poleward of 55).
Sphere-correct: rotation about the equatorial axis perpendicular to the group's central meridian, so shapes keep true size
(a plain +lat shift would stretch them east-west near the pole). Moves barsdorf intents/regions/lines + pins; political regions
are outputs (re-run 2026-10-07c_borders.py); generated coast is stale until regen. Backup + atomic write.
  tools/map/py tools/map/transforms/2026-10-11_barsdorf_north.py"""
import json, math, os, shutil, tempfile, time
from pathlib import Path
import numpy as np
DLAT = 3.0
STATE = Path(__file__).resolve().parents[3] / "library/maps/map-state.json"
s = json.loads(STATE.read_text())
mine = [f for f in s["features"] if "barsdorf" in (f.get("groups") or []) and f.get("kind") != "political"]
lons = np.concatenate([np.asarray(f["geometry"], float)[:, 0] for f in mine if f["type"] == "land"])
c = math.radians(0.5 * (lons.min() + lons.max()))
k = np.array([math.sin(c), -math.cos(c), 0.0])   # west-pointing equatorial axis: +angle tips the central meridian north
a = math.radians(DLAT); K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
Rm = np.eye(3) + math.sin(a) * K + (1 - math.cos(a)) * K @ K
def rot(lon, lat):
    lo, la = np.radians(np.asarray(lon, float)), np.radians(np.asarray(lat, float))
    v = Rm @ np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])
    return np.degrees(np.arctan2(v[1], v[0])), np.degrees(np.arcsin(np.clip(v[2], -1, 1)))
n = 0
for f in mine:
    if f["type"] == "pin":
        lo, la = rot([f["geometry"]["lon"]], [f["geometry"]["lat"]]); f["geometry"]["lon"], f["geometry"]["lat"] = round(float(lo[0]), 4), round(float(la[0]), 4); n += 1
    elif f["type"] in ("land", "region", "line"):
        g = np.asarray(f["geometry"], float); lo, la = rot(g[:, 0], g[:, 1])
        f["geometry"] = [[round(float(x), 4), round(float(y), 4)] for x, y in zip(lo, la)]; n += 1
hist = STATE.parent / "history"; hist.mkdir(exist_ok=True)
shutil.copy2(STATE, hist / f"map-state-{time.strftime('%Y%m%d-%H%M%S')}-pre-barsdorf-north.json")
fd, tmp = tempfile.mkstemp(dir=STATE.parent)
with os.fdopen(fd, "w") as fh: fh.write(json.dumps(s, indent=2, ensure_ascii=False) + "\n")
os.replace(tmp, STATE); print(f"moved {n} barsdorf features {DLAT:+} deg north (axis lon {math.degrees(c) - 90:.1f})")
