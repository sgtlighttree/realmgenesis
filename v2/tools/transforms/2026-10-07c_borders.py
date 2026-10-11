#!/usr/bin/env python3
"""Redraw country regions from current land shapes (2026-10-07).

Land is rasterised on a lon/lat grid. Each land cell gets a nation:
  - islands by their landmass tag (asahi-arch -> asahi, wi -> wi, barsdorf -> barsdorf,
    ne-peninsula -> dolma);
  - the main Haldei body by nearest country-tagged pin (km), with the existing
    NW bandit-lands region kept as an override.
Each nation mask is traced (coastgen.marching_squares), simplified, and written back
as region features (one per ring), replacing the old country regions. Sub-regions
(Old Reach) and the WI footprint are left alone. Re-runnable.
"""
import json, math, os, sys, tempfile, time, shutil
from pathlib import Path
import numpy as np
from PIL import Image, ImageDraw

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tools/map/terrain"))
from coastgen import marching_squares, dp_tol  # noqa: E402

STATE = ROOT / "library/maps/map-state.json"
RES = 0.05                      # degrees per cell
EPS = 0.04                      # DP tolerance, degrees
MIN_CELLS = 6                   # drop specks
# Medieval claims: a nation holds land within CLAIM_KM of its own pins (main body only).
# Where two claims meet, the border is hard (traced as a line); other edges are fuzzy frontier.
CLAIM_KM = {"dolma": 650, "asahi": 300}
FUZZY = {"dolma", "asahi", "bandit-lands"}
# Default owner of land NOT covered by any claim (islands only; the main body is frontier unless claimed).
# NE archipelago: loosely held, bandit country (ruling 2026-10-07).
LAND_TO_NATION = {"asahi-arch": "asahi", "wi": "wi", "barsdorf": "barsdorf", "ne-peninsula": "bandit-lands"}
NAMES = {"dolma": "Dolma", "asahi": "Asahi", "wi": "Emirate of Majid", "barsdorf": "Confederation of Barsdorf (fallen)",
         "bandit-lands": "NW bandit-lands"}
REPLACE = {"Dolma", "Asahi — Archipelago", "Asahi — Mainland", "Emirate of Majid",
           "Confederation of Barsdorf (fallen)", "NW bandit-lands", "NE bandit-isles"}


def main():
    s = json.loads(STATE.read_text())
    feats = s["features"]
    lands = [f for f in feats if f["type"] == "land"]
    allpts = np.array([p for f in lands for p in f["geometry"]])
    lon0, lat0 = allpts[:, 0].min() - 1, allpts[:, 1].min() - 1
    lon1, lat1 = allpts[:, 0].max() + 1, allpts[:, 1].max() + 1
    nx, ny = int((lon1 - lon0) / RES) + 1, int((lat1 - lat0) / RES) + 1
    px = lambda ring: [((lo - lon0) / RES, (lat1 - la) / RES) for lo, la in ring]

    def raster(rings):
        im = Image.new("L", (nx, ny), 0); d = ImageDraw.Draw(im)
        for r in rings: d.polygon(px(r), fill=1)
        return np.array(im, bool)

    nations = ["dolma", "asahi", "wi", "barsdorf", "bandit-lands"]
    lab = np.full((ny, nx), -1, np.int16)
    claims = [f for f in feats if f["type"] == "region" and f.get("kind") == "claim"]
    oldreach = next((f for f in feats if f["type"] == "region" and f["name"] == "Old Reach"), None)
    if oldreach: oldreach["kind"] = "subregion"   # always Dolman (ruling 2026-10-07)
    # main body (haldei-tagged land): nearest country pin
    body = raster([f["geometry"] for f in lands if "haldei" in (f.get("groups") or [])])
    seeds = [(f["geometry"]["lon"], f["geometry"]["lat"], nations.index(t))
             for f in feats if f["type"] == "pin" for t in (f.get("groups") or []) if t in ("dolma", "asahi")]
    seeds = [sd for sd in seeds if body[min(ny - 1, int((lat1 - sd[1]) / RES)), min(nx - 1, int((sd[0] - lon0) / RES))]]
    yy, xx = np.nonzero(body)
    lon = lon0 + xx * RES; lat = lat1 - yy * RES
    DEG = s["planet"]["radius_km"] * math.pi / 180
    best = np.full(len(xx), 1e18); who = np.zeros(len(xx), np.int16)
    for slo, sla, k in seeds:
        dx = (lon - slo) * np.cos(np.radians((lat + sla) / 2)); dy = lat - sla
        dd = dx * dx + dy * dy
        m = dd < best; best[m] = dd[m]; who[m] = k
    reach = np.array([CLAIM_KM.get(nations[k], 1e9) for k in who]) / DEG
    who = np.where(best <= reach ** 2, who, -1)                # beyond every claim: unclaimed frontier
    if not claims:   # ---- bootstrap claims once from radii + rulings, then they are Matt's to sculpt
        lab[yy, xx] = who
        bandit = next((f for f in feats if f["type"] == "region" and f["name"] == "NW bandit-lands"), None)
        if bandit: lab[body & raster([bandit["geometry"]])] = nations.index("bandit-lands")
        if oldreach: lab[body & raster([oldreach["geometry"]])] = nations.index("dolma")
        cap = next(f for f in feats if f["id"] == "pin-capital")["geometry"]
        cols = body.sum(0); bx = np.nonzero(cols)[0]
        c0 = int((cap["lon"] + 8 - lon0) / RES); c1 = int(bx.max() - 8 / RES)   # stay off the body's far tip
        nz = [c for c in range(c0, max(c0 + 1, c1)) if cols[c] > 0]
        if nz:
            neck = min(nz, key=lambda c: cols[c])          # the isthmus: narrowest column east of the core
            east = np.zeros_like(body); east[:, neck:] = True
            lab[body & east & (lab < 0)] = nations.index("dolma")
            print(f"isthmus at lon {lon0 + neck * RES:.1f}; NE peninsula claimed by Dolma")
        for k, nat in enumerate(nations[:1] + nations[1:2] + nations[4:5]):
            kk = nations.index(nat); M = lab == kk
            if M.sum() < MIN_CELLS: continue
            Mp = np.pad(M, 1)
            for j, r in enumerate(sorted(marching_squares(Mp.astype(np.float32), 0.5, Mp), key=len, reverse=True)):
                r = r - 1.0
                if len(r) < 12: continue
                ll = np.stack([lon0 + r[:, 0] * RES, lat1 - r[:, 1] * RES], 1); ll = ll[dp_tol(ll, EPS * 2)]
                claims.append({"id": f"claim-{nat}-{j + 1:02d}", "type": "region", "kind": "claim", "name": f"{NAMES[nat]} claim",
                               "status": "proposed", "groups": [nat],
                               "note": "Editable claim. Political borders = land within the claim (re-run 2026-10-07c_borders.py). Overlapping claims split by nearest pin (hard border).",
                               "cites": "", "geometry": [[round(float(a), 4), round(float(b), 4)] for a, b in ll]})
        feats.extend(claims)
        print("bootstrapped claims:", len(claims))
    # ---- claims decide: every land cell inside a claim belongs to it (overlaps -> nearest own pin)
    lab[:] = -1
    allland = raster([f["geometry"] for f in lands])
    cmask = {}
    for c in claims:
        nat = (c.get("groups") or [None])[0]
        if nat in nations: cmask[nat] = cmask.get(nat, np.zeros_like(body)) | raster([c["geometry"]])
    cnt = sum(m.astype(np.int8) for m in cmask.values()) if cmask else np.zeros(body.shape, np.int8)
    for nat, m in cmask.items(): lab[m & allland & (cnt == 1)] = nations.index(nat)
    both = allland & (cnt > 1)
    if both.any():
        yy2, xx2 = np.nonzero(both); lon2 = lon0 + xx2 * RES; lat2 = lat1 - yy2 * RES
        bd = np.full(len(xx2), 1e18); wh = np.full(len(xx2), -1, np.int16)
        pins = [(f["geometry"]["lon"], f["geometry"]["lat"], nations.index(t)) for f in feats if f["type"] == "pin"
                for t in (f.get("groups") or []) if t in cmask]
        for slo, sla, k in pins:
            inside = cmask[nations[k]][yy2, xx2]
            dd = ((lon2 - slo) * np.cos(np.radians((lat2 + sla) / 2))) ** 2 + (lat2 - sla) ** 2
            m = inside & (dd < bd); bd[m] = dd[m]; wh[m] = k
        lab[yy2, xx2] = wh
    if oldreach:   # canon: the Old Reach is always Dolman, whatever the claims say
        lab[allland & raster([oldreach["geometry"]])] = nations.index("dolma")
    for f in lands:
        tags = f.get("groups") or []
        nat = next((LAND_TO_NATION[t] for t in tags if t in LAND_TO_NATION), None)
        if nat and "haldei" not in tags:
            m = raster([f["geometry"]]); lab[m & (lab < 0) & (cnt == 0)] = nations.index(nat)

    new = []
    for k, nat in enumerate(nations):
        M = lab == k
        if M.sum() < MIN_CELLS: continue
        Mp = np.pad(M, 1)
        rings = marching_squares(Mp.astype(np.float32), 0.5, Mp)
        i = 0
        for r in rings:
            r = r - 1.0
            if len(r) < 4: continue
            area = 0.5 * abs(np.dot(r[:, 0], np.roll(r[:, 1], 1)) - np.dot(r[:, 1], np.roll(r[:, 0], 1)))
            if area < MIN_CELLS: continue
            ll = np.stack([lon0 + r[:, 0] * RES, lat1 - r[:, 1] * RES], 1)
            keep = dp_tol(ll, EPS); ll = ll[keep]
            if len(ll) < 4: continue
            i += 1
            cap_lon = next(f for f in feats if f["id"] == "pin-capital")["geometry"]["lon"]
            nm = ("NE bandit-isles" if float(ll[:, 0].mean()) > cap_lon + 15 else "NW bandit-lands") if nat == "bandit-lands" else NAMES[nat]
            new.append({"id": f"reg-{nat}-{i:03d}", "type": "region", "kind": "political", "name": nm,
                        "status": "proposed", "groups": [nat], "edge": "fuzzy" if nat in FUZZY else "hard",
                        "note": "Auto border from land + nearest country pin (2026-10-07). Re-run tools/map/transforms/2026-10-07c_borders.py after reshaping land.",
                        "cites": "", "geometry": [[round(float(a), 4), round(float(b), 4)] for a, b in ll]})
        print(f"{nat:<13} {i} ring(s), {int(M.sum())} cells")

    # hard borders: where two nations' cells touch on land, trace the boundary as lines
    hard = []
    for a in range(len(nations)):
        for b in range(a + 1, len(nations)):
            A, B = lab == a, lab == b
            if not A.any() or not B.any(): continue
            touch = np.zeros_like(A)
            for sy, sx in ((0, 1), (1, 0), (1, 1), (1, -1)):
                Bs = np.roll(np.roll(B, sy, 0), sx, 1); As = np.roll(np.roll(A, sy, 0), sx, 1)
                touch |= (A & Bs) | (B & As)
            if touch.sum() < 3: continue
            Ap = np.pad(A, 1)
            for r in marching_squares(Ap.astype(np.float32), 0.5, Ap):
                r = r - 1.0; iy = np.clip(np.round(r[:, 1]).astype(int), 0, ny - 1); ix = np.clip(np.round(r[:, 0]).astype(int), 0, nx - 1)
                near = np.zeros(len(r), bool)
                for oy in (-1, 0, 1):
                    for ox in (-1, 0, 1):
                        near |= touch[np.clip(iy + oy, 0, ny - 1), np.clip(ix + ox, 0, nx - 1)]
                run = []
                for j in list(range(len(r))) + [0]:
                    if near[j]: run.append(r[j])
                    elif len(run) >= 3: hard.append(np.array(run)); run = []
                    else: run = []
                if len(run) >= 3: hard.append(np.array(run))
            nm = f"Border: {NAMES[nations[a]]} / {NAMES[nations[b]]}"
            for r in hard:
                ll = np.stack([lon0 + r[:, 0] * RES, lat1 - r[:, 1] * RES], 1)
                if len(ll) > 4: ll = ll[np.unique(np.r_[0, dp_tol(ll, EPS / 2), len(ll) - 1])]
                new.append({"id": f"brd-{nations[a]}-{nations[b]}-{len(new):03d}", "type": "line", "kind": "border",
                            "name": nm, "status": "proposed", "groups": [nations[a], nations[b]],
                            "note": "Auto border (contested edge, hard). Re-run 2026-10-07c_borders.py after reshaping land.",
                            "cites": "", "geometry": [[round(float(x), 4), round(float(y), 4)] for x, y in ll]})
            hard = []
    print("hard border lines:", sum(1 for f in new if f["type"] == "line"))
    hist = STATE.parent / "history"; hist.mkdir(exist_ok=True)
    shutil.copy2(STATE, hist / f"map-state-{time.strftime('%Y%m%d-%H%M%S')}-pre-borders.json")
    kept = [f for f in feats if not ((f["type"] == "region" and (f["name"] in REPLACE or f["id"].startswith("reg-") and f.get("kind") == "political" and "Auto border" in f.get("note", "")))
                                     or (f["type"] == "line" and f.get("kind") == "border" and "Auto border" in f.get("note", "")))]
    s["features"] = kept + new
    fd, tmp = tempfile.mkstemp(dir=STATE.parent)
    with os.fdopen(fd, "w") as fh: fh.write(json.dumps(s, indent=2, ensure_ascii=False) + "\n")
    os.replace(tmp, STATE)
    print(f"features written: {len(new)}; removed {len(feats) - len(kept)}")


if __name__ == "__main__":
    main()
