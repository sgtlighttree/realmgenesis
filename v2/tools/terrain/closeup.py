#!/usr/bin/env python3
"""Regional close-up of the LOCKED map (2026-10-11).

A zoom, not a regen: the locked generated coast (library/maps/generated/<group>.json) is rasterised at a fine grid and
displaced only by a separate fine-octave stack (wavelengths <= 16 km, amplitude ~1.5 km), so the close-up can't drift from
the overview. Relief = the v3 formula on the same world-coherent stacks + a fine ridged stack + the Spire spec; biomes =
the locked biome raster with smoothed class boundaries. The global COAST_WLS / REL_*_WLS stacks are untouched.

Writes ONLY library/maps/closeups/<name>.{png,json}; never generated/ or map-state.json.
  tools/map/py tools/map/terrain/closeup.py --name dolma-heartland --group haldei --bbox 52.1 31.9 55.3 35.0 --res 0.2
Gates (printed + stored): close-up coast vs locked coast (p95/max km), land/water/coastal pins, a named strait stays open
(--strait-pin, min width km), and pin land checks (land at given km E/W of a pin).
"""
import argparse, json, math, os, sys, time
from pathlib import Path
import numpy as np
from PIL import Image, ImageDraw, ImageFont, ImageFilter
sys.path.insert(0, str(Path(__file__).resolve().parent))
import coastgen as cg   # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
FINE_COAST_WLS = [16.0, 8.0, 4.0, 2.0, 1.0]     # close-up only; finer than the global stack's 30 km floor
FINE_RID_WLS = [24.0, 12.0, 6.0, 3.0]
FINE_GAIN = 0.6

def fine_fields(pj, x0, y0, res, nx, ny, seed, R, rows=64):
    """World-coherent fine coast noise (~[-0.5,0.5]) and fine ridged relief (~[0,1]) on the close-up grid."""
    cn = np.empty((ny, nx), np.float32); rn = np.empty((ny, nx), np.float32); xs = x0 + np.arange(nx) * res
    nc = sum(FINE_GAIN ** o for o in range(len(FINE_COAST_WLS))); nr = sum(0.5 ** o for o in range(len(FINE_RID_WLS)))
    for a in range(0, ny, rows):
        b = min(ny, a + rows); X, Y = np.meshgrid(xs, y0 + np.arange(a, b) * res)
        lon, lat = pj.inv(X.ravel(), Y.ravel()); P = cg.sph_xyz(lon, lat)
        vc = cg.octave_values(P, FINE_COAST_WLS, seed + 701, R); vr = cg.octave_values(P, FINE_RID_WLS, seed + 801, R, ridged=True)
        cn[a:b] = (0.5 * sum(FINE_GAIN ** o * v for o, v in enumerate(vc)) / nc).reshape(b - a, nx)
        rn[a:b] = (sum(0.5 ** o * v for o, v in enumerate(vr)) / nr).reshape(b - a, nx)
    return cn, rn

def decode_biome(b):
    flat = np.empty(b["w"] * b["h"], np.uint8); i = 0; r = b["rle"]
    for k in range(0, len(r), 2): flat[i:i + r[k + 1]] = r[k]; i += r[k + 1]
    return flat.reshape(b["h"], b["w"])[::-1]   # rows south -> north

def soft_biome(cls, bb, res_deg, lon, lat, sigma_px=1.2):
    """Smooth class boundaries: per-class masks blurred on the 0.1 deg grid, bilinear-sampled, argmax."""
    best = np.zeros(lon.shape, np.float32); out = np.zeros(lon.shape, np.uint8)
    px = (lon - bb[0]) / res_deg - 0.5; py = (lat - bb[1]) / res_deg - 0.5
    for c in np.unique(cls):
        if c == 0: continue
        m = Image.fromarray(((cls == c) * 255).astype(np.uint8)).filter(ImageFilter.GaussianBlur(sigma_px))
        v = cg.bilin(np.asarray(m, np.float32) / 255.0, px, py)
        sel = v > best; best = np.where(sel, v, best); out = np.where(sel, c, out)
    return out

def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--name", required=True); ap.add_argument("--group", default="haldei")
    ap.add_argument("--bbox", type=float, nargs=4, required=True, metavar=("LON0", "LAT0", "LON1", "LAT1"))
    ap.add_argument("--res", type=float, default=0.2, help="km per cell"); ap.add_argument("--exag", type=float, default=20.0, help="hillshade vertical exaggeration"); ap.add_argument("--amp", type=float, default=6.0, help="fine coast displacement, km")
    ap.add_argument("--strait-pin", default="Emerald Straits"); ap.add_argument("--strait-min-km", type=float, default=4.0)
    ap.add_argument("--land-check", default="pin-weir:320@8,15;140@8,15", help="pin:bearing@km,km;... land expected (map bearings; the Dolma pocket is rotated 50 deg cw from canon, so canon W->320, E->140)"); ap.add_argument("--sub-inset", default="Capital of Dolma (name needed):5", help="pin name:half-km for a magnified corner inset ('' = none)")
    ap.add_argument("--state", default=str(ROOT / "library/maps/map-state.json")); ap.add_argument("--out-dir", default=str(ROOT / "library/maps/closeups"))
    a = ap.parse_args(); t0 = time.time()
    out_dir = Path(a.out_dir).resolve()
    if out_dir == (ROOT / "library/maps/generated").resolve(): raise SystemExit("close-ups never write generated/")
    state = json.load(open(a.state)); R = float(state["planet"]["radius_km"]); seed = int(state["gen"]["world_seed"])
    gen = json.load(open(ROOT / f"library/maps/generated/{a.group}.json"))
    lo0, la0, lo1, la1 = a.bbox; pj = cg.Proj(0.5 * (lo0 + lo1), 0.5 * (la0 + la1), R); res = a.res
    edge = np.vstack([np.stack([np.linspace(lo0, lo1, 50), np.full(50, v)], 1) for v in (la0, la1)] + [np.stack([np.full(50, v), np.linspace(la0, la1, 50)], 1) for v in (lo0, lo1)])
    ex, ey = pj.fwd(edge[:, 0], edge[:, 1]); x0, x1, y0, y1 = ex.min(), ex.max(), ey.min(), ey.max()
    nx = int((x1 - x0) / res) + 1; ny = int((y1 - y0) / res) + 1
    print(f"[{a.name}] grid {nx}x{ny} @ {res} km ({nx*ny/1e6:.2f}M), centre {pj.lon0:.2f}E {pj.lat0:.2f}N")
    pad = 1.0
    def near(f):
        b = f.get("bbox"); return b and b[0] < lo1 + pad and b[2] > lo0 - pad and b[1] < la1 + pad and b[3] > la0 - pad
    def km_px(r):
        x, y = pj.fwd(*cg.densify_ll(np.asarray(r, float), 0.01).T); return np.stack([(x - x0) / res, (y - y0) / res], 1)
    land_r = [f["lod"][0] for f in gen["features"] if f.get("type") == "land" and near(f)]
    lake_r = [f["lod"][0] for f in gen["features"] if f.get("kind") == "lake" and near(f)]
    inside = cg.rasterize_cov([km_px(r) for r in land_r], ny, nx) >= 0.5
    if lake_r: inside &= ~(cg.rasterize_cov([km_px(r) for r in lake_r], ny, nx) >= 0.5)
    sdf0 = cg.sdf_from_mask(inside, res)   # locked coast, km (+ inland)
    cn, rn = fine_fields(pj, x0, y0, res, nx, ny, seed, R)
    C = sdf0 + np.float32(a.amp) * 2.0 * cn   # cn ~[-0.5,0.5] -> displacement ~ +-amp
    pins, notes = cg.assign_pins(state, a.group, pj, cg.group_rings_all(state), 0.0)
    inwin = lambda p: lo0 <= p["lon"] <= lo1 and la0 <= p["lat"] <= la1
    pins = {k: [p for p in v if inwin(p)] for k, v in pins.items()}
    cg.force_pins(C, pins, x0, y0, res, sig_km=1.5)
    L = C > 0
    for m, keep_if in ((L, True), (~L, False)):   # drop specks < 1 km2 (land specks at sea, puddles on land)
        lab, n = cg.label_mask(m, 8); sz = np.bincount(lab.ravel())
        small = (sz < 1.0 / res ** 2); small[0] = False
        L = np.where(small[lab], not keep_if, L)
    sdfL = cg.sdf_from_mask(L, res)
    # ---- gates
    inner = lambda px, py: (px > 10) & (px < nx - 11) & (py > 10) & (py < ny - 11)   # ignore the window edge
    lp = np.vstack([cg.densify_km(np.stack(pj.fwd(*np.asarray(r, float).T), 1), res) for r in land_r]); lpx, lpy = (lp[:, 0] - x0) / res, (lp[:, 1] - y0) / res
    k = inner(lpx, lpy); d1 = np.abs(cg.bilin(sdfL, lpx[k], lpy[k]))
    Lc = L.copy(); Lc[0, :] = Lc[-1, :] = Lc[:, 0] = Lc[:, -1] = False   # one-cell sea frame so edge-cut land traces closed rings
    Cc = np.where(Lc, np.maximum(C, 1e-3), np.minimum(C, -1e-3)).astype(np.float32)
    raw = cg.marching_squares32(Cc, 0.0, Lc); del Cc, Lc
    rings_km = [np.stack([x0 + r[:, 0] * res, y0 + r[:, 1] * res], 1) for r in raw]
    fp = np.vstack(rings_km); k2 = inner((fp[:, 0] - x0) / res, (fp[:, 1] - y0) / res)
    d2 = np.abs(cg.bilin(sdf0, (fp[k2, 0] - x0) / res, (fp[k2, 1] - y0) / res))
    dev = np.concatenate([d1, d2]); gates = dict(coast_dev_p95_km=round(float(np.percentile(dev, 95)), 2), coast_dev_max_km=round(float(dev.max()), 2))
    val = lambda p: float(cg.bilin(sdfL, np.array([(p["x"] - x0) / res]), np.array([(p["y"] - y0) / res]))[0])
    gates["land_pins_failed"] = [p["name"] for p in pins["land"] if val(p) <= 0]
    gates["water_pins_failed"] = [p["name"] for p in pins["water"] if val(p) >= 0]
    gates["coastal_pins_km"] = {p["name"]: round(abs(val(p)), 2) for p in pins["coastal"]}
    sp = [p for p in pins["water"] if p["name"] == a.strait_pin]
    if sp: w = 2 * -val(sp[0]); gates["strait"] = dict(pin=a.strait_pin, width_at_pin_km=round(w, 2), ok=w >= a.strait_min_km)
    pid = {f["id"]: f for f in state["features"]}
    if a.land_check:
        pk, spec_ = a.land_check.split(":", 1); g = pid[pk]["geometry"]; lc = []; kd = R * math.pi / 180
        for part in spec_.split(";"):
            brg, kms = part.split("@"); t = math.radians(float(brg))
            for dk in [float(v) for v in kms.split(",")]:
                lon = g["lon"] + dk * math.sin(t) / (kd * math.cos(math.radians(g["lat"]))); lat = g["lat"] + dk * math.cos(t) / kd
                x, y = pj.fwd(lon, lat); lc.append(dict(at=f"{dk:.0f}km @{float(brg):.0f}", land=bool(cg.bilin(sdfL, np.array([(x - x0) / res]), np.array([(y - y0) / res]))[0] > 0)))
        gates["land_checks"] = {pid[pk]["name"]: lc}
    gates["ok"] = (gates["coast_dev_p95_km"] <= 2.0 and not gates["land_pins_failed"] and not gates["water_pins_failed"]
                   and all(v <= 3.0 for v in gates["coastal_pins_km"].values()) and gates.get("strait", {}).get("ok", True)
                   and all(c["land"] for v in gates.get("land_checks", {}).values() for c in v))
    print(json.dumps(gates))
    # ---- relief (v3 formula on the world stacks + fine ridged + Spire), metres
    base, rid = cg.relief_noise(pj, x0, y0, res, nx, ny, seed, R); base = base * 2.5
    up = 1 - np.exp(-np.clip(sdfL, 0, None) / 140.0); up_s = 1 - np.exp(-np.clip(sdfL, 0, None) / 4.0)
    rel = 0.12 + 0.55 * up + 0.6 * up * rid + 0.2 * up * (0.5 + 0.5 * base) + 0.22 * up_s * (rn - 0.3)
    ev = np.where(L, np.clip(rel, 0, None) * 1000.0, 0).astype(np.float32)
    spec = cg.climate_spec(state); sps = spec.get("spire", {})
    X, Y = np.meshgrid(x0 + np.arange(nx) * res, y0 + np.arange(ny) * res); LON, LAT = pj.inv(X, Y); del X, Y
    d_sp = None
    if "lon" in sps:
        d_sp = cg._gc_km(LON, LAT, sps["lon"], sps["lat"], R)
        ev = np.where(L, np.maximum(ev, cg.spire_height(d_sp, sps)), ev)
    # ---- biome
    b = gen["biome"]; cls = soft_biome(decode_biome(b), b["bbox"], b["res_deg"], LON, LAT)
    cls = np.where(cls == 0, 4, cls)
    # the overview climate grid (~9 km cells) can't resolve the Spire: drop its alpine/ice there, re-derive from the fine relief
    cls = np.where(np.isin(cls, (6, 7)), 4, cls)
    cls = np.where(L & (ev >= 1200.0), 6, cls)   # coastgen's alpine_min_m
    if d_sp is not None: cls = np.where(d_sp + 3.0 * cn <= sps.get("ice_radius_km", 15.0), 7, cls)   # edge ~+-0.75 km
    # ---- render
    pal = np.array([[96, 96, 96]] + [[int(c[i:i + 2], 16) for i in (1, 3, 5)] for _, c in cg.BIOMES], np.float32) / 255
    gy, gx = np.gradient(ev * a.exag / (res * 1000.0))
    shade = np.clip(0.7 + 0.45 * (gx - gy) / math.sqrt(2), 0.25, 1.25).astype(np.float32)
    land_c = pal[cls] * (0.5 + 0.55 * shade[..., None])
    depth = np.clip(-sdfL / 40.0, 0, 1)[..., None]
    sea_c = np.array([0.62, 0.82, 0.90], np.float32) * (1 - depth) + np.array([0.12, 0.28, 0.48], np.float32) * depth
    img = Image.fromarray((np.clip(np.flipud(np.where(L[..., None], land_c, sea_c)), 0, 1) * 255).astype(np.uint8), "RGB")
    del land_c, sea_c, shade
    try:
        font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 15); small = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 12)
    except Exception: font = small = ImageFont.load_default()
    d = ImageDraw.Draw(img, "RGBA")
    to_px = lambda x, y: ((x - x0) / res, (ny - 1) - (y - y0) / res)
    for K in rings_km:   # coastline only: break the ring where it runs along the window frame
        seg = []
        for p in list(K) + [K[0]]:
            px_, py_ = to_px(*p)
            if 1.5 < px_ < nx - 2.5 and 1.5 < py_ < ny - 2.5: seg.append((px_, py_))
            else:
                if len(seg) > 1: d.line(seg, fill=(20, 20, 20, 255), width=1)
                seg = []
        if len(seg) > 1: d.line(seg, fill=(20, 20, 20, 255), width=1)
    for f in state["features"]:   # rivers / roads / canals drawn in the state, clipped to the window
        if f.get("type") != "line" or f.get("kind") not in ("river", "road", "canal"): continue
        g = np.asarray(f["geometry"], float)
        if not ((g[:, 0] >= lo0) & (g[:, 0] <= lo1) & (g[:, 1] >= la0) & (g[:, 1] <= la1)).any(): continue
        dg = np.vstack([a_ + (b_ - a_) * t_ for a_, b_ in zip(g[:-1], g[1:]) for t_ in np.linspace(0, 1, 40, endpoint=False)[:, None]] + [g[-1:]])   # open polyline
        x, y = pj.fwd(dg[:, 0], dg[:, 1]); P_ = [to_px(a_, b_) for a_, b_ in zip(x, y)]
        if f["kind"] == "river":
            d.line(P_, fill=(40, 90, 170, 255), width=4, joint="curve"); d.line(P_, fill=(110, 170, 230, 255), width=2, joint="curve")
        else:
            for n_ in range(0, len(P_) - 1, 2): d.line(P_[n_:n_ + 2], fill=(250, 235, 200, 235), width=2)
            lx, ly = P_[len(P_) // 2]; d.text((lx + 4, ly + 4), f.get("name", "").split(" (")[0], fill=(250, 235, 200, 255), font=small, stroke_width=2, stroke_fill=(0, 0, 0, 255))
    for gl in np.arange(math.ceil(lo0 * 2) / 2, lo1, 0.5):   # graticule, 0.5 deg
        x, y = pj.fwd(np.full(60, gl), np.linspace(la0, la1, 60)); d.line([to_px(*p) for p in zip(x, y)], fill=(255, 255, 255, 70), width=1)
        px_, py_ = to_px(x[0], y[0]); d.text((px_ + 3, min(ny - 16, py_ - 16)), f"{gl:g}°E", fill=(255, 255, 255, 200), font=small)
    for gl in np.arange(math.ceil(la0 * 2) / 2, la1, 0.5):
        x, y = pj.fwd(np.linspace(lo0, lo1, 60), np.full(60, gl)); d.line([to_px(*p) for p in zip(x, y)], fill=(255, 255, 255, 70), width=1)
        px_, py_ = to_px(x[0], y[0]); d.text((max(3, px_ + 3), py_ - 15), f"{gl:g}°N", fill=(255, 255, 255, 200), font=small)
    # pins: greedy label placement over 8 candidate offsets, avoiding earlier labels and dots
    allp = [(p, "land") for p in pins["land"]] + [(p, "water") for p in pins["water"]] + [(p, "coastal") for p in pins["coastal"] if p not in pins["land"]]
    for f in state["features"]:   # in-window pins assign_pins skipped (owned by no group / unplaced kinds are excluded upstream)
        if f.get("type") == "pin" and f.get("kind") != "unplaced" and inwin(f["geometry"]) and not any(q["name"] == f.get("name") for q, _ in allp):
            x, y = pj.fwd(f["geometry"]["lon"], f["geometry"]["lat"]); allp.append((dict(name=f.get("name"), lon=f["geometry"]["lon"], lat=f["geometry"]["lat"], x=float(x), y=float(y)), "other"))
    sub = None
    if a.sub_inset:
        sn, sh = a.sub_inset.rsplit(":", 1); sp_ = [p for p, _ in allp if p["name"] == sn]
        if sp_: sub = (sp_[0], float(sh))
    def in_sub(p): return sub and abs(p["x"] - sub[0]["x"]) <= sub[1] and abs(p["y"] - sub[0]["y"]) <= sub[1]
    boxes = []
    def place(dd, p, kind, fnt, scale_to_px, label_ok=True):
        x, y = scale_to_px(p["x"], p["y"]); col = {"land": (220, 40, 40), "water": (40, 110, 240), "coastal": (250, 205, 0), "other": (230, 230, 230)}[kind]
        rr = 4; dd.ellipse([x - rr, y - rr, x + rr, y + rr], fill=col + (255,), outline=(255, 255, 255, 255)); boxes.append((x - rr, y - rr, x + rr, y + rr))
        if not label_ok: return
        t = p["name"]; tw, th = dd.textbbox((0, 0), t, font=fnt)[2:]; W_, H_ = dd.im.size
        def ov(bx): return sum(max(0, min(bx[2], q[2]) - max(bx[0], q[0])) * max(0, min(bx[3], q[3]) - max(bx[1], q[1])) for q in boxes)
        cands = [(x + ox, y + oy, x + ox + tw, y + oy + th) for ox, oy in ((7, -th // 2), (-tw - 7, -th // 2), (-tw // 2, -th - 7), (-tw // 2, 7),
                 (7, -th - 4), (7, 4), (-tw - 7, -th - 4), (-tw - 7, 4), (7, -2 * th - 6), (7, th + 6), (-tw - 7, -2 * th - 6), (-tw - 7, th + 6))]
        cands = [b for b in cands if b[0] >= 2 and b[1] >= 2 and b[2] <= W_ - 2 and b[3] <= H_ - 2] or cands
        bx = min(cands, key=ov)
        boxes.append(bx); dd.text(bx[:2], t, fill=(255, 255, 255, 255), font=fnt, stroke_width=2, stroke_fill=(0, 0, 0, 255))
    crop = None
    if sub:   # magnified corner inset for the capital cluster: crop the bare map first, reserve its corner for labels
        cpx, cpy = to_px(sub[0]["x"], sub[0]["y"]); h = sub[1] / res; S = max(2, int(380 / (2 * h)))
        crop = img.crop((int(cpx - h), int(cpy - h), int(cpx + h), int(cpy + h))).resize((int(2 * h) * S, int(2 * h) * S), Image.LANCZOS)
        ix0 = img.width - crop.width - 12; iy0 = img.height - crop.height - 12; boxes.append((ix0 - 4, iy0 - 16, img.width, img.height))
    for p, kind in allp: place(d, p, kind, font, to_px, label_ok=not in_sub(p))
    if sub:
        d.rectangle([cpx - h, cpy - h, cpx + h, cpy + h], outline=(255, 255, 255, 230), width=2)
        dd = ImageDraw.Draw(crop, "RGBA"); boxes_main = boxes; boxes = [(-1e9, -1e9, 1e9, 0), (-1e9, crop.height, 1e9, 1e9), (-1e9, -1e9, 0, 1e9), (crop.width, -1e9, 1e9, 1e9), (0, 0, 340, 30)]
        sp2 = lambda x, y: (((x - x0) / res - (cpx - h)) * S, (((ny - 1) - (y - y0) / res) - (cpy - h)) * S)
        for p, kind in allp:
            if in_sub(p): place(dd, p, kind, font, sp2)
        dd.rectangle([0, 0, crop.width - 1, crop.height - 1], outline=(255, 255, 255, 255), width=3)
        dd.text((8, 6), f"{sub[0]['name'].split(' /')[0]} area ({2*sub[1]:.0f} km box)", fill=(255, 255, 255, 255), font=font, stroke_width=2, stroke_fill=(0, 0, 0, 255))
        img.paste(crop, (ix0, iy0)); boxes = boxes_main; d = ImageDraw.Draw(img, "RGBA")
    sb = 25.0 / res; bx0, by0 = 20, ny - 40   # scale bar
    d.rectangle([bx0 - 6, by0 - 22, bx0 + sb + 50, by0 + 14], fill=(0, 0, 0, 140))
    d.line([(bx0, by0), (bx0 + sb, by0)], fill=(255, 255, 255, 255), width=4); d.text((bx0 + sb + 6, by0 - 8), "25 km", fill=(255, 255, 255, 255), font=small)
    foot = Image.new("RGB", (img.width, 52), (22, 22, 22)); fd = ImageDraw.Draw(foot)
    fd.text((8, 6), f"{a.name} · zoom of the locked map ({a.group}) · {res} km/px · coast p95 {gates['coast_dev_p95_km']} km from locked · gates {'ok' if gates['ok'] else 'FAIL'}", fill=(255, 255, 255), font=small)
    if sps.get("profile") == "ease_in":
        def el(pin_id):
            g_ = pid[pin_id]["geometry"]; dk = float(cg._gc_km(np.array(g_["lon"]), np.array(g_["lat"]), sps["lon"], sps["lat"], R))
            return f"{pid[pin_id]['name'].split(' (')[0]} ({dk:.1f} km) ~{float(cg.spire_height(np.array(dk), sps)):,.0f} m"
        sn = f"Spire: ease-in^{sps.get('power', 4):g}, {sps.get('base_km', 18):g} km base, ice {sps.get('ice_radius_km', 15):g} km · " + " · ".join(el(k) for k in ("pin-bastion", "pin-threshold") if k in pid)
    else: sn = "Spire relief provisional (legacy exp profile)"
    fd.text((8, 26), sn + " · river, coast road from the map state", fill=(255, 200, 120), font=small)
    out = Image.new("RGB", (img.width, img.height + foot.height)); out.paste(img, (0, 0)); out.paste(foot, (0, img.height))
    out_dir.mkdir(parents=True, exist_ok=True)
    png = out_dir / f"{a.name}.png"; out.save(png, optimize=True)
    rings_ll = [[[round(float(lo), 5), round(float(la), 5)] for lo, la in zip(*pj.inv(K[cg.dp_tol(K, 0.15 * res)][:, 0], K[cg.dp_tol(K, 0.15 * res)][:, 1]))] for K in rings_km]
    meta = dict(name=a.name, group=a.group, bbox=a.bbox, res_km=res, fine_amp_km=a.amp, fine_coast_wls=FINE_COAST_WLS, fine_rid_wls=FINE_RID_WLS, seed=seed,
                source=f"library/maps/generated/{a.group}.json", gates=gates, notes=notes, rings=rings_ll, runtime_s=round(time.time() - t0, 1))
    (out_dir / f"{a.name}.json").write_text(json.dumps(meta, separators=(",", ":")) + "\n")
    print(f"wrote {png} ({out.width}x{out.height}), {len(rings_km)} rings, {time.time()-t0:.1f}s")

if __name__ == "__main__":
    main()
