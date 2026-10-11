#!/usr/bin/env python3
"""coastgen: v3 = signed-distance coast generator (see tools/map/README.md, Generation contract v3); v2 = raster-heightfield coastline generator for the TCoTG map (numpy + PIL only).

Pipeline: intent polygons -> blurred base field -> 3D (unit-sphere) fBm + domain warp + ridged
(optionally anisotropic) + arc terms -> pins (Gaussian bumps) -> sea level by area quantile ->
lake fill -> hand-rolled marching squares -> per-island min-area filter ->
(v2: tolerance DP at full res, area-corrected sea level, nested LODs, biome raster) -> library/maps/generated/<group>.json + PNG preview.

Schema in/out: [lon, lat] degrees, single-ring polygons. Never writes map-state.json.

    tools/map/py tools/map/terrain/coastgen.py --preset wi --seed 23 --config tools/map/terrain/configs/wi-v2.json
"""
import argparse, copy, hashlib, heapq, json, math, os, resource, sys, time
try:
    import numpy as np
    from PIL import Image, ImageDraw, ImageFont
except ImportError as e:   # bare python3 without the deps (e.g. after a brew Python switch)
    raise SystemExit(f"coastgen: {e}. Run map tools through tools/map/py (venv; setup: tools/map/setup-venv.sh)")

R_KM_DEFAULT = 4500.0
VERSION = "0.2"

# ----------------------------------------------------------------------------- geometry helpers
def sph_xyz(lon, lat):
    lo, la = np.radians(lon), np.radians(lat)
    c = np.cos(la)
    return np.stack([c * np.cos(lo), c * np.sin(lo), np.sin(la)], axis=0)

def ring_area_km2(ring, R):
    """Spherical-ish area (cylindrical equal-area trapezoid sum), km^2, positive."""
    a = np.asarray(ring, float)
    lo, la = np.radians(a[:, 0]), np.radians(a[:, 1])
    lo2, la2 = np.roll(lo, -1), np.roll(la, -1)
    d = (lo2 - lo + np.pi) % (2 * np.pi) - np.pi
    return abs(float(np.sum(d * (np.sin(la) + np.sin(la2))))) / 2 * R * R

def point_in_ring(pts, ring):
    """Even-odd, vectorised over pts (N,2). Planar lon/lat."""
    pts = np.asarray(pts, float); r = np.asarray(ring, float)
    x, y = pts[:, 0][:, None], pts[:, 1][:, None]
    x1, y1 = r[:, 0][None], r[:, 1][None]
    x2, y2 = np.roll(r[:, 0], -1)[None], np.roll(r[:, 1], -1)[None]
    cond = (y1 > y) != (y2 > y)
    with np.errstate(divide="ignore", invalid="ignore"):
        xi = x1 + (y - y1) * (x2 - x1) / (y2 - y1)
    return (np.sum(cond & (x < xi), axis=1) % 2) == 1

def _segs_cross(A, B):
    """A,B: (n,2,2) edge arrays. Proper-crossing matrix (n_a, n_b)."""
    def cr(o, p, q):
        return (p[..., 0] - o[..., 0]) * (q[..., 1] - o[..., 1]) - (p[..., 1] - o[..., 1]) * (q[..., 0] - o[..., 0])
    p1, p2 = A[:, None, 0], A[:, None, 1]
    p3, p4 = B[None, :, 0], B[None, :, 1]
    d1, d2 = cr(p3, p4, p1), cr(p3, p4, p2)
    d3, d4 = cr(p1, p2, p3), cr(p1, p2, p4)
    return (d1 * d2 < 0) & (d3 * d4 < 0)

def ring_self_intersects(ring):
    r = np.asarray(ring, float); n = len(r)
    E = np.stack([r, np.roll(r, -1, axis=0)], axis=1)
    M = _segs_cross(E, E)
    idx = np.arange(n)
    adj = (np.abs(idx[:, None] - idx[None]) <= 1) | (np.abs(idx[:, None] - idx[None]) == n - 1)
    return bool(np.any(M & ~adj))

def rings_cross(r1, r2):
    a = np.asarray(r1, float); b = np.asarray(r2, float)
    Ea = np.stack([a, np.roll(a, -1, axis=0)], axis=1)
    Eb = np.stack([b, np.roll(b, -1, axis=0)], axis=1)
    return bool(np.any(_segs_cross(Ea, Eb)))

# ----------------------------------------------------------------------------- 3D gradient noise
_U = np.uint32
def _hash(ix, iy, iz, seed):
    with np.errstate(over="ignore"):
        h = (ix.astype(_U) * _U(73856093)) ^ (iy.astype(_U) * _U(19349663)) ^ (iz.astype(_U) * _U(83492791)) ^ _U(seed & 0xFFFFFFFF)
        h ^= h >> _U(13); h = h * _U(1274126177); h ^= h >> _U(16)
        h = h * _U(2246822519); h ^= h >> _U(15)
    return h

def perlin3(p, seed):
    """p: (3,N) float64 -> (N,) approx [-1,1]."""
    i = np.floor(p); f = p - i; ii = i.astype(np.int64)
    u = f * f * f * (f * (f * 6 - 15) + 10)
    out = 0.0
    for dx in (0, 1):
        wx = u[0] if dx else 1 - u[0]
        for dy in (0, 1):
            wy = u[1] if dy else 1 - u[1]
            for dz in (0, 1):
                wz = u[2] if dz else 1 - u[2]
                h = _hash(ii[0] + dx, ii[1] + dy, ii[2] + dz, seed)
                gx = (h & _U(255)).astype(np.float64) / 127.5 - 1
                gy = ((h >> _U(8)) & _U(255)).astype(np.float64) / 127.5 - 1
                gz = ((h >> _U(16)) & _U(255)).astype(np.float64) / 127.5 - 1
                out = out + wx * wy * wz * (gx * (f[0] - dx) + gy * (f[1] - dy) + gz * (f[2] - dz))
    return out * 1.25

_ROT = np.array([[0.36, 0.48, -0.80], [-0.80, 0.60, 0.00], [0.48, 0.64, 0.60]])  # orthonormal-ish
def fbm3(p, freq, octaves, gain, seed, ridged=False, lac=2.0):
    """p unit vectors (3,N); freq in cycles per unit-sphere-radius. Normalised ~[-1,1] (ridged: ~[0,1])."""
    q = p * freq; amp = 1.0; tot = 0.0; acc = 0.0
    for o in range(octaves):
        n = perlin3(q + 17.31 * o, seed + 101 * o)
        if ridged:
            n = 1.0 - np.abs(n); n = n * n
        acc = acc + amp * n; tot += amp
        q = (_ROT @ q) * lac; amp *= gain
    return acc / tot

# ----------------------------------------------------------------------------- grid + blur
class Grid:
    def __init__(self, bbox, d):
        self.lon0, self.lat0, lon1, lat1 = bbox
        self.d = d
        self.nx = int(round((lon1 - self.lon0) / d)) + 1
        self.ny = int(round((lat1 - self.lat0) / d)) + 1
    def lon(self): return self.lon0 + np.arange(self.nx) * self.d
    def lat(self): return self.lat0 + np.arange(self.ny) * self.d

def gblur(F, sx, sy):
    """Separable gaussian blur via FFT, zero padded. sx, sy in cells."""
    ny, nx = F.shape
    px, py = int(3 * sx) + 2, int(3 * sy) + 2
    P = np.zeros((ny + 2 * py, nx + 2 * px)); P[py:py + ny, px:px + nx] = F
    fy = np.fft.fftfreq(P.shape[0])[:, None]; fx = np.fft.rfftfreq(P.shape[1])[None]
    K = np.exp(-2 * np.pi ** 2 * ((sx * fx) ** 2 + (sy * fy) ** 2))
    out = np.fft.irfft2(np.fft.rfft2(P) * K, s=P.shape)
    return out[py:py + ny, px:px + nx].astype(np.float32)

def rasterize_intent(g, polys):
    """polys: [{'ring':[[lon,lat]..], 'weight':w}] -> float32 (ny,nx), lat increasing with row."""
    im = Image.new("F", (g.nx, g.ny), 0.0); dr = ImageDraw.Draw(im)
    for p in sorted(polys, key=lambda p: p.get("weight", 1.0)):
        xy = [((lo - g.lon0) / g.d, (g.ny - 1) - (la - g.lat0) / g.d) for lo, la in p["ring"]]
        dr.polygon(xy, fill=float(p.get("weight", 1.0)))
    return np.flipud(np.asarray(im, dtype=np.float32)).copy()

# ----------------------------------------------------------------------------- arcs
def _catmull(pts, per=40):
    P = np.vstack([2 * pts[0] - pts[1], pts, 2 * pts[-1] - pts[-2]])
    out = []
    for i in range(1, len(P) - 2):
        p0, p1, p2, p3 = P[i - 1], P[i], P[i + 1], P[i + 2]
        for t in np.linspace(0, 1, per, endpoint=False):
            out.append(0.5 * ((2 * p1) + (-p0 + p2) * t + (2 * p0 - 5 * p1 + 4 * p2 - p3) * t * t + (-p0 + 3 * p1 - 3 * p2 + p3) * t ** 3))
    out.append(P[-2])
    return np.array(out)

def densify_arc(pts, kind, R):
    """pts [[lon,lat]..] -> dense unit xyz (M,3), cumulative km (M,)."""
    xyz = sph_xyz(np.array([p[0] for p in pts]), np.array([p[1] for p in pts])).T
    if kind == "gc" or len(xyz) == 2:
        out = []
        for a, b in zip(xyz[:-1], xyz[1:]):
            om = math.acos(max(-1, min(1, float(a @ b)))) or 1e-9
            n = max(2, int(om * R / 10))
            for t in np.linspace(0, 1, n, endpoint=False):
                out.append((math.sin((1 - t) * om) * a + math.sin(t * om) * b) / math.sin(om))
        out.append(xyz[-1]); D = np.array(out)
    else:
        D = _catmull(xyz, 40)
    D = D / np.linalg.norm(D, axis=1, keepdims=True)
    seg = np.linalg.norm(np.diff(D, axis=0), axis=1) * R
    return D, np.concatenate([[0], np.cumsum(seg)])

class ArcSet:
    def __init__(self, arcs, R, seed):
        self.R = R; self.arcs = []
        rng = np.random.default_rng(seed + 991)
        for a in arcs:
            D, s = densify_arc(a["pts"], a.get("kind", "spline"), R)
            L = s[-1]
            spacing = a.get("bead_km", 110.0)
            nb = max(3, int(L / spacing) + 2)
            vals = rng.standard_normal(nb) * a.get("bead_sigma", 1.0)
            amp_s = np.exp(vals - 0.3) if a.get("bead_sigma", 1.0) > 0 else np.ones(nb)
            amp_s = np.minimum(amp_s, a.get("bead_max", 3.5))
            # smooth interpolation over arclength
            xs = np.linspace(0, L, nb)
            fine = np.interp(s, xs, amp_s)
            fine = np.convolve(np.pad(fine, 3, mode="edge"), np.ones(7) / 7, mode="valid")
            taper = np.minimum(1, np.minimum(s, L - s) / max(1.0, a.get("taper_km", 150.0)))
            self.arcs.append(dict(D=D.astype(np.float32), amp=(fine * taper).astype(np.float32),
                                  w=a.get("amp", 1.0), sig=a.get("width_km", 45.0), spec=a))
    def eval(self, P):
        """P (3,N) unit vectors -> (N,) arc term."""
        out = np.zeros(P.shape[1], np.float32); Pf = P.T.astype(np.float32)
        for a in self.arcs:
            dot = Pf @ a["D"].T
            j = np.argmax(dot, axis=1)
            d2 = np.maximum(2 - 2 * dot[np.arange(len(j)), j], 0)
            d = np.sqrt(d2) * self.R
            out += a["w"] * a["amp"][j] * np.exp(-0.5 * (d / a["sig"]) ** 2)
        return out

# ----------------------------------------------------------------------------- field
DEFAULT_STYLE = dict(
    intent_w=1.0, intent_sig_km=(300.0, 70.0), intent_mix=0.6,
    warp1_amp_km=380.0, warp1_wl_km=1500.0, warp2_amp_km=75.0, warp2_wl_km=300.0,
    base_w=0.8, base_wl_km=1300.0, base_oct=4,
    ridged_w=0.8, ridged_wl_km=600.0, ridged_oct=5, ridged_gate=0.4,
    aniso=1.0, trend_deg=0.0,
    arc_w=1.0,
    detail_w=0.5, detail_wl_km=170.0, detail_oct=5, gain=0.55,
    edge_taper_deg=2.0,
)

def compute_field(g, cfg, R, log):
    st = {**DEFAULT_STYLE, **cfg.get("style", {})}
    seed = cfg["seed"]
    polys = cfg.get("intent", [])
    mask = rasterize_intent(g, polys) if polys else np.zeros((g.ny, g.nx), np.float32)
    kmx = R * math.radians(g.d)
    latc = math.radians(0.5 * (g.lat0 + g.lat0 + (g.ny - 1) * g.d))
    s1, s2 = st["intent_sig_km"]
    sy1, sy2 = s1 / kmx, s2 / kmx
    sx1, sx2 = sy1 / max(0.3, math.cos(latc)), sy2 / max(0.3, math.cos(latc))
    I = st["intent_mix"] * gblur(mask, sx1, sy1) + (1 - st["intent_mix"]) * gblur(mask, sx2, sy2)
    I = (2 * I - 1).astype(np.float32)  # roughly -1..1, 0 at 50% blurred coverage
    log(f"intent raster {g.nx}x{g.ny}, polys={len(polys)}")

    arcs = ArcSet(cfg.get("arcs", []), R, seed) if cfg.get("arcs") else None
    lonv, latv = g.lon(), g.lat()
    LON, LAT = np.meshgrid(lonv, latv)
    N = g.nx * g.ny
    flat_lon, flat_lat = LON.ravel(), LAT.ravel()
    E = np.empty(N, np.float32); RID = np.empty(N, np.float32)
    # trend direction (global tangent at region centre)
    clon = 0.5 * (lonv[0] + lonv[-1]); clat = 0.5 * (latv[0] + latv[-1])
    c = sph_xyz(np.array([clon]), np.array([clat]))[:, 0]
    east = np.array([-math.sin(math.radians(clon)), math.cos(math.radians(clon)), 0.0])
    north = np.cross(c, east)
    b = math.radians(st["trend_deg"])  # bearing clockwise from north
    u = math.cos(b) * north + math.sin(b) * east
    aniso = st["aniso"]
    A1, A2 = st["warp1_amp_km"] / R, st["warp2_amp_km"] / R
    f = lambda wl: R / wl
    CH = 16384
    Ifl = I
    for a0 in range(0, N, CH):
        sl_ = slice(a0, min(N, a0 + CH))
        P = sph_xyz(flat_lon[sl_], flat_lat[sl_])
        fw1 = f(st["warp1_wl_km"]); fw2 = f(st["warp2_wl_km"])
        w1 = np.stack([fbm3(P, fw1, 3, 0.5, seed + 11 * k + 1) for k in range(3)])
        w2 = np.stack([fbm3(P, fw2, 3, 0.5, seed + 23 * k + 7) for k in range(3)])
        P1 = P + A1 * w1; P1 /= np.linalg.norm(P1, axis=0)
        P2 = P + A2 * w2; P2 /= np.linalg.norm(P2, axis=0)
        # intent sampled at big-warped position
        lo1 = np.degrees(np.arctan2(P1[1], P1[0])); la1 = np.degrees(np.arcsin(np.clip(P1[2], -1, 1)))
        fx = (lo1 - g.lon0) / g.d; fy = (la1 - g.lat0) / g.d
        ix = np.clip(np.floor(fx).astype(int), 0, g.nx - 2); iy = np.clip(np.floor(fy).astype(int), 0, g.ny - 2)
        tx = np.clip(fx - ix, 0, 1); ty = np.clip(fy - iy, 0, 1)
        Iv = (Ifl[iy, ix] * (1 - tx) * (1 - ty) + Ifl[iy, ix + 1] * tx * (1 - ty) + Ifl[iy + 1, ix] * (1 - tx) * ty + Ifl[iy + 1, ix + 1] * tx * ty)
        outside = (fx < 0) | (fx > g.nx - 1) | (fy < 0) | (fy > g.ny - 1)
        Iv = np.where(outside, -1.0, Iv)
        base = fbm3(P1, f(st["base_wl_km"]), st["base_oct"], st["gain"], seed + 301)
        Pr = P1 + (1.0 / aniso - 1.0) * (u[:, None] * (u @ P1)[None]) if aniso != 1.0 else P1
        rid = fbm3(Pr, f(st["ridged_wl_km"]), st["ridged_oct"], st["gain"], seed + 401, ridged=True)
        det = fbm3(P2, f(st["detail_wl_km"]), st["detail_oct"], st["gain"], seed + 501)
        gate = st["ridged_gate"] + np.clip(Iv * 0.5 + 0.5, 0, 1)
        e = st["intent_w"] * Iv + st["base_w"] * base + st["ridged_w"] * (rid - 0.45) * gate + st["detail_w"] * det
        if arcs is not None:
            arc = arcs.eval(P2)
            e = e + st["arc_w"] * arc * (0.75 + 0.5 * (det * 0.5 + 0.5))
        E[sl_] = e; RID[sl_] = rid
    E = E.reshape(g.ny, g.nx); RID = RID.reshape(g.ny, g.nx)
    # edge taper: push to ocean near bbox borders so no land is clipped
    tp = st["edge_taper_deg"] / g.d
    xi = np.minimum(np.arange(g.nx), np.arange(g.nx)[::-1])[None]; yi = np.minimum(np.arange(g.ny), np.arange(g.ny)[::-1])[:, None]
    dist = np.minimum(xi, yi).astype(np.float32)
    ramp = np.clip(dist / max(tp, 1), 0, 1)
    E = E - (1 - ramp) ** 2 * 3.0
    E = (E - E.mean()) / E.std()
    return E.astype(np.float32), RID

# ----------------------------------------------------------------------------- sea level + pins
def cell_area_weights(g, R):
    lat = np.radians(g.lat())
    return (R * math.radians(g.d)) ** 2 * np.cos(lat)  # km2 per cell per row

def solve_sea_level(E, wrow, target):
    """Sea level so that land (E>sl) area == target (cell-based)."""
    W = np.broadcast_to(wrow[:, None], E.shape).ravel()
    e = E.ravel(); o = np.argsort(-e)
    cs = np.cumsum(W[o])
    k = int(np.searchsorted(cs, target))
    k = min(k, len(o) - 1)
    return float(0.5 * (e[o[k]] + e[o[min(k + 1, len(o) - 1)]]))

def _bump(g, lon, lat, sig_km, R):
    kmy = R * math.radians(g.d); lonv, latv = g.lon(), g.lat()
    rad = 4.5 * sig_km
    x0 = max(0, int((lon - g.lon0) / g.d - rad / (kmy * max(0.2, math.cos(math.radians(lat)))) ))
    x1 = min(g.nx, int((lon - g.lon0) / g.d + rad / (kmy * max(0.2, math.cos(math.radians(lat)))) ) + 2)
    y0 = max(0, int((lat - g.lat0) / g.d - rad / kmy)); y1 = min(g.ny, int((lat - g.lat0) / g.d + rad / kmy) + 2)
    LO, LA = np.meshgrid(lonv[x0:x1], latv[y0:y1])
    dx = (LO - lon) * math.cos(math.radians(lat)) * R * math.pi / 180
    dy = (LA - lat) * R * math.pi / 180
    return (y0, y1, x0, x1), np.exp(-0.5 * (dx * dx + dy * dy) / sig_km ** 2).astype(np.float32)

def _sample(E, g, lon, lat):
    fx = (lon - g.lon0) / g.d; fy = (lat - g.lat0) / g.d
    ix = int(max(0, min(g.nx - 2, math.floor(fx)))); iy = int(max(0, min(g.ny - 2, math.floor(fy))))
    tx, ty = fx - ix, fy - iy
    return float(E[iy, ix] * (1 - tx) * (1 - ty) + E[iy, ix + 1] * tx * (1 - ty) + E[iy + 1, ix] * (1 - tx) * ty + E[iy + 1, ix + 1] * tx * ty)

def apply_pins(E, g, land_pins, water_pins, wrow, target, R, st, log):
    sig_l, sig_w = st.get("pin_land_sig_km", 28.0), st.get("pin_water_sig_km", 38.0)
    margin = st.get("pin_margin", 0.45)
    sl = solve_sea_level(E, wrow, target)
    for it in range(10):
        changed = 0
        for p in land_pins:
            ev = min(_sample(E, g, p["lon"] + dx, p["lat"] + dy) for dx in (-0.12, 0, 0.12) for dy in (-0.12, 0, 0.12))
            if ev < sl + margin:
                (y0, y1, x0, x1), B = _bump(g, p["lon"], p["lat"], sig_l, R)
                E[y0:y1, x0:x1] += (sl + margin - ev + 0.15) * B; changed += 1
        for p in water_pins:
            ev = max(_sample(E, g, p["lon"] + dx, p["lat"] + dy) for dx in (-0.12, 0, 0.12) for dy in (-0.12, 0, 0.12))
            if ev > sl - margin:
                (y0, y1, x0, x1), B = _bump(g, p["lon"], p["lat"], sig_w, R)
                E[y0:y1, x0:x1] -= (ev - (sl - margin) + 0.15) * B; changed += 1
        sl = solve_sea_level(E, wrow, target)
        log(f"  pins iter {it}: adjusted {changed}, sl={sl:.3f}")
        if not changed: break
    return sl

# ----------------------------------------------------------------------------- marching squares
def fill_lakes(L):
    """Water not 4-connected to the border becomes land. Returns new bool mask."""
    im = Image.fromarray(np.where(L, 255, 0).astype(np.uint8)).copy()
    ImageDraw.floodfill(im, (0, 0), 128)
    a = np.asarray(im)
    return a != 128

def marching_squares(E, sl, L):
    """Closed rings (x,y float grid coords), land on the left (CCW outer). L is the land mask (lakes filled)."""
    ny, nx = E.shape
    E = np.where(L & (E <= sl), sl + 1e-4, E)           # filled lake cells: no crossings inside
    E = np.where(~L & (E > sl), sl - 1e-4, E)           # forced-ocean border
    C = [L[:-1, :-1], L[:-1, 1:], L[1:, 1:], L[1:, :-1]]  # corners v0..v3
    yy, xx = np.mgrid[0:ny - 1, 0:nx - 1]
    NH = ny * nx
    eid = [yy * nx + xx, NH + yy * nx + (xx + 1), (yy + 1) * nx + xx, NH + yy * nx + xx]
    ex = [C[k] & ~C[(k + 1) % 4] for k in range(4)]
    en = [~C[k] & C[(k + 1) % 4] for k in range(4)]
    saddle = (C[0] == C[2]) & (C[1] == C[3]) & (C[0] != C[1])
    center = (E[:-1, :-1] + E[:-1, 1:] + E[1:, 1:] + E[1:, :-1]) * 0.25 > sl
    use_next = saddle & center
    nxt = np.full(2 * NH, -1, np.int64)
    for k in range(4):
        def partner(direction):
            pj = np.full(C[0].shape, -1, np.int64)
            for step in (3, 2, 1):
                j = (k + direction * step) % 4
                pj = np.where(en[j], j, pj)
            return pj
        pj = np.where(use_next, partner(+1), partner(-1))
        m = ex[k] & (pj >= 0)
        if not m.any(): continue
        to = np.choose(np.clip(pj, 0, 3), [eid[0], eid[1], eid[2], eid[3]])
        nxt[eid[k][m]] = to[m]
    used = np.nonzero(nxt >= 0)[0]
    posx = np.zeros(2 * NH, np.float32); posy = np.zeros(2 * NH, np.float32)
    ids = np.concatenate([used, nxt[used]]); ids = np.unique(ids)
    hm = ids < NH
    h = ids[hm]; y, x = h // nx, h % nx
    t = (sl - E[y, x]) / (E[y, x + 1] - E[y, x]); posx[h] = x + t; posy[h] = y
    v = ids[~hm] - NH; y, x = v // nx, v % nx
    t = (sl - E[y, x]) / (E[y + 1, x] - E[y, x]); posx[v + NH] = x; posy[v + NH] = y + t
    nl = nxt.tolist(); seen = np.zeros(2 * NH, bool); rings = []
    for s in used.tolist():
        if seen[s]: continue
        cur = s; seq = []
        while cur >= 0 and not seen[cur]:
            seen[cur] = True; seq.append(cur); cur = nl[cur]
        if cur == s and len(seq) >= 4:
            q = np.array(seq); rings.append(np.stack([posx[q], posy[q]], axis=1))
    return rings


def dp_by_count(P, n_target):
    """Closed ring P (N,2) planar -> kept index array (sorted), ~n_target vertices, by DP significance."""
    N = len(P)
    if n_target >= N: return np.arange(N)
    d0 = np.hypot(*(P - P[0]).T); i1 = int(np.argmax(d0))
    def chain(a, b):
        if b < a: b += N
        return a, b
    def maxdist(a, b):
        idx = np.arange(a, b + 1) % N
        Q = P[idx]; A, B = Q[0], Q[-1]; d = B - A; L = math.hypot(*d)
        if L < 1e-12: dd = np.hypot(*(Q - A).T)
        else: dd = np.abs(d[0] * (Q[:, 1] - A[1]) - d[1] * (Q[:, 0] - A[0])) / L
        j = int(np.argmax(dd)); return float(dd[j]), (a + j)
    kept = {0, i1}; heap = []
    for a, b in (chain(0, i1), chain(i1, 0)):
        if b - a >= 2:
            dmax, j = maxdist(a, b); heapq.heappush(heap, (-dmax, a, b, j))
    while heap and len(kept) < n_target:
        nd, a, b, j = heapq.heappop(heap)
        kept.add(j % N)
        for aa, bb in ((a, j), (j, b)):
            if bb - aa >= 2:
                dmax, jj = maxdist(aa, bb); heapq.heappush(heap, (-dmax, aa, bb, jj))
    return np.array(sorted(kept))

def relief_field(E, RID, sl, g, cfg, R):
    H = E - sl
    land = np.clip(H, 0, None)
    gate = np.clip(H / 0.8, 0, 1) ** 0.8
    rel = land ** 0.9 * 0.55 + 0.55 * gate * (RID - 0.35)
    rel = np.where(H > 0, np.clip(rel, 0, None), H * 0.25)
    for cn in cfg.get("canyons", []):  # carve (relief only; coast untouched)
        D, s = densify_arc(cn["pts"], "spline", R)
        LO, LA = np.meshgrid(g.lon(), g.lat()); P = sph_xyz(LO.ravel(), LA.ravel()).T.astype(np.float32)
        d = np.empty(len(P), np.float32)
        for a0 in range(0, len(P), 20000):
            dot = P[a0:a0 + 20000] @ D.astype(np.float32).T
            d[a0:a0 + 20000] = np.sqrt(np.maximum(2 - 2 * dot.max(axis=1), 0)) * R
        d = d.reshape(g.ny, g.nx)
        w = cn.get("width_km", 18.0)
        rel = rel - cn.get("depth", 0.5) * np.exp(-0.5 * (d / w) ** 2) * (H > 0.02)
    return rel.astype(np.float32)

# ----------------------------------------------------------------------------- v2: tolerance simplify, LODs, topology
def dp_tol(P, eps):
    """Closed ring P (N,2) planar -> sorted kept indices (Douglas-Peucker, tolerance eps)."""
    N = len(P)
    if N <= 4: return np.arange(N)
    d0 = np.hypot(*(P - P[0]).T); i1 = int(np.argmax(d0))
    keep = np.zeros(N, bool); keep[0] = keep[i1] = True
    Pc = np.vstack([P, P[:1]])
    stack = [(0, i1), (i1, N)]
    while stack:
        a, b = stack.pop()
        if b - a < 2: continue
        A, B = Pc[a], Pc[b]; d = B - A; L = math.hypot(d[0], d[1])
        Q = Pc[a + 1:b]
        if L < 1e-12: dd = np.hypot(Q[:, 0] - A[0], Q[:, 1] - A[1])
        else: dd = np.abs(d[0] * (Q[:, 1] - A[1]) - d[1] * (Q[:, 0] - A[0])) / L
        j = int(np.argmax(dd))
        if dd[j] > eps:
            k = a + 1 + j; keep[k] = True
            stack.append((a, k)); stack.append((k, b))
    return np.nonzero(keep)[0]

def planar(r):
    return np.stack([r[:, 0] * math.cos(math.radians(float(r[:, 1].mean()))), r[:, 1]], axis=1)

def rnd_ring(r, nd=4):
    q = np.round(np.asarray(r, float), nd)
    keep = np.ones(len(q), bool); keep[1:] = np.any(q[1:] != q[:-1], axis=1)
    q = q[keep]
    if len(q) > 1 and np.all(q[0] == q[-1]): q = q[:-1]
    return q

def _cross_pairs(a0, a1, b0, b1):
    def cr(o, p, q): return (p[:, 0] - o[:, 0]) * (q[:, 1] - o[:, 1]) - (p[:, 1] - o[:, 1]) * (q[:, 0] - o[:, 0])
    return (cr(b0, b1, a0) * cr(b0, b1, a1) < 0) & (cr(a0, a1, b0) * cr(a0, a1, b1) < 0)

def find_crossings(rings):
    """Grid-bucketed proper-crossing search over all edges of all (implicitly closed) rings.
    Returns a set of (i, j) ring-index pairs (i <= j; i == j is a self crossing)."""
    if not rings: return set()
    A = np.vstack(rings); B = np.vstack([np.roll(r, -1, axis=0) for r in rings])
    RID = np.concatenate([np.full(len(r), i) for i, r in enumerate(rings)])
    K = np.concatenate([np.arange(len(r)) for r in rings]); nr = np.array([len(r) for r in rings])
    NE = len(A)
    if NE < 4: return set()
    cs = max(float(np.median(np.hypot(*(B - A).T))) * 4, 1e-5)
    x0 = np.minimum(A[:, 0], B[:, 0]); x1 = np.maximum(A[:, 0], B[:, 0]); y0 = np.minimum(A[:, 1], B[:, 1]); y1 = np.maximum(A[:, 1], B[:, 1])
    ix0 = np.floor(x0 / cs).astype(np.int64); iy0 = np.floor(y0 / cs).astype(np.int64)
    ix1 = np.floor(x1 / cs).astype(np.int64); iy1 = np.floor(y1 / cs).astype(np.int64)
    mx, my = ix0.min(), iy0.min(); ix0 -= mx; ix1 -= mx; iy0 -= my; iy1 -= my
    nx, ny = ix1 - ix0, iy1 - iy0
    longe = np.nonzero((nx > 12) | (ny > 12))[0]
    isl = np.zeros(NE, bool); isl[longe] = True
    keys = []; eids = []
    for dx in range(int(min(nx.max(), 12)) + 1):
        for dy in range(int(min(ny.max(), 12)) + 1):
            idx = np.nonzero((nx >= dx) & (ny >= dy) & ~isl)[0]
            if len(idx): keys.append((ix0[idx] + dx) * 4000037 + (iy0[idx] + dy)); eids.append(idx)
    pa, pb = [], []
    if keys:
        keys = np.concatenate(keys); eids = np.concatenate(eids)
        o = np.argsort(keys, kind="stable"); keys = keys[o]; eids = eids[o]
        _, cnt = np.unique(keys, return_counts=True)
        for off in range(1, int(cnt.max())):
            m = keys[:-off] == keys[off:]
            if not m.any(): break
            pa.append(eids[:-off][m]); pb.append(eids[off:][m])
    for e in longe:   # very long edges: brute force against everything
        pa.append(np.full(NE, e)); pb.append(np.arange(NE))
    if not pa: return set()
    pa = np.concatenate(pa).astype(np.int64); pb = np.concatenate(pb).astype(np.int64)
    lo = np.minimum(pa, pb); hi = np.maximum(pa, pb); m = lo != hi
    u = np.unique(lo[m] * NE + hi[m]); lo = u // NE; hi = u % NE
    same = RID[lo] == RID[hi]; dk = np.abs(K[lo] - K[hi])
    adj = same & ((dk <= 1) | (dk == nr[RID[lo]] - 1))
    lo = lo[~adj]; hi = hi[~adj]
    out = set()
    for s in range(0, len(lo), 1_000_000):
        l, h = lo[s:s + 1_000_000], hi[s:s + 1_000_000]
        hit = _cross_pairs(A[l], B[l], A[h], B[h])
        for i, j in zip(RID[l][hit].tolist(), RID[h][hit].tolist()): out.add((min(i, j), max(i, j)))
    return out

def simplify_full(rings, eps, max_verts):
    """Full-res LOD: tolerance DP (eps in degrees, grown until <= max_verts), then 4-dp rounding. Returns (rings, eps_used)."""
    Ps = [planar(r) for r in rings]
    while True:
        out = []
        for r, P in zip(rings, Ps):
            q = rnd_ring(r[dp_tol(P, eps)])
            out.append(q if len(q) >= 3 else None)
        tot = sum(len(q) for q in out if q is not None)
        if tot <= max_verts or eps > 1.0: return out, eps
        eps *= 1.25

def build_lods(full, fracs, log, pl=None):
    """full: list of (N,2) rings. Returns lods[l] = list of rings for l = 0..len(fracs); nested DP subsets, crossings repaired."""
    levels = [list(full)]
    prev_idx = [np.arange(len(r)) for r in full]
    pl = pl or planar
    P = [pl(r) for r in full]
    for fr in fracs:
        tgt = [min(len(r), max(6, int(round(len(r) * fr)))) if len(r) > 8 else len(r) for r in full]
        idx = [None] * len(full)
        def make(i):
            base = prev_idx[i]
            if tgt[i] >= len(base): return base
            return base[dp_by_count(P[i][base], tgt[i])]
        for i in range(len(full)): idx[i] = make(i)
        for _ in range(8):
            cur = [full[i][idx[i]] for i in range(len(full))]
            bad = {i for pair in find_crossings(cur) for i in pair}
            if not bad: break
            for i in bad:
                tgt[i] = min(len(prev_idx[i]), int(tgt[i] * 1.35) + 2); idx[i] = make(i)
        else:
            log(f"  warn: LOD {fr} still has {len(bad)} crossing rings, using the finer level for them")
            for i in bad: idx[i] = prev_idx[i]
        levels.append([full[i][idx[i]] for i in range(len(full))]); prev_idx = idx
    return levels

# ----------------------------------------------------------------------------- v2: contour + area-corrected solve
def land_rings(E, sl, g, R, cfg, land_pins):
    L = E > sl
    L[0, :] = L[-1, :] = False; L[:, 0] = L[:, -1] = False
    L = fill_lakes(L)
    out = []
    dmin = cfg.get("min_island_km2", 90.0); rules = cfg.get("min_island_rules", [])
    for r in marching_squares(E, sl, L):
        sa = 0.5 * np.sum(r[:, 0] * np.roll(r[:, 1], -1) - np.roll(r[:, 0], -1) * r[:, 1])
        if sa <= 0: continue
        ll = np.stack([g.lon0 + r[:, 0] * g.d, g.lat0 + r[:, 1] * g.d], axis=1)
        a = ring_area_km2(ll, R)
        if a < 1.0: continue
        c = ll.mean(axis=0); need = dmin
        for ru in rules:
            b = ru["bbox"]
            if b[0] <= c[0] <= b[2] and b[1] <= c[1] <= b[3]: need = ru["min_km2"]; break
        keep = a >= need
        if not keep and cfg.get("keep_pin_islands", True) and land_pins:
            keep = any(point_in_ring(np.array([[p["lon"], p["lat"]]]), ll)[0] for p in land_pins)
        if keep:
            for _ in range(cfg.get("smooth_iters", 1)):   # Chaikin corner cutting: de-angularise the 0.05 deg contour
                nx_ = np.roll(ll, -1, axis=0); ll = np.stack([0.75 * ll + 0.25 * nx_, 0.25 * ll + 0.75 * nx_], axis=1).reshape(-1, 2)
            out.append(ll)
    return out

def lake_rings(E, sl, g, R, cfg):
    """Enclosed water (not 4-connected to the grid edge) >= basins.min_lake_km2 -> outer rings [lon,lat] (islands inside lakes dropped).
    Land rings keep lakes filled (single-ring schema); lakes are separate features drawn on top."""
    bc = cfg.get("basins") or {}; amin = bc.get("min_lake_km2", 2500.0)
    L = E > sl; L[0, :] = L[-1, :] = True; L[:, 0] = L[:, -1] = True
    wl, nw = label_mask(~L, 4)
    if nw == 0: return []
    L[0, :] = L[-1, :] = False; L[:, 0] = L[:, -1] = False
    edge = set(np.unique(np.concatenate([wl[1], wl[-2], wl[:, 1], wl[:, -2]])).tolist()) - {0}
    wrow = cell_area_weights(g, R)
    area = np.bincount(wl.ravel(), weights=np.broadcast_to(wrow[:, None], wl.shape).ravel(), minlength=nw + 1)
    keep = [k for k in range(1, nw + 1) if k not in edge and area[k] >= amin]
    if not keep: return []
    K = np.isin(wl, keep)
    D = K.copy()
    for _ in range(2):
        D = D | np.roll(D, 1, 0) | np.roll(D, -1, 0) | np.roll(D, 1, 1) | np.roll(D, -1, 1)
    F = np.where(D, sl - E, -1.0).astype(np.float32)
    out = []
    for r in marching_squares(F, 0.0, K):
        sa = 0.5 * np.sum(r[:, 0] * np.roll(r[:, 1], -1) - np.roll(r[:, 0], -1) * r[:, 1])
        if sa <= 0: continue
        ll = np.stack([g.lon0 + r[:, 0] * g.d, g.lat0 + r[:, 1] * g.d], axis=1)
        if ring_area_km2(ll, R) < amin * 0.5: continue
        nx_ = np.roll(ll, -1, axis=0); ll = np.stack([0.75 * ll + 0.25 * nx_, 0.25 * ll + 0.75 * nx_], axis=1).reshape(-1, 2)
        out.append(ll)
    return out

def params_hash(cfg):
    c = {k: v for k, v in cfg.items() if k != "output"}
    return hashlib.sha1(json.dumps(c, sort_keys=True, default=str).encode()).hexdigest()[:12]

# ----------------------------------------------------------------------------- v2: biome tint
BIOMES = [("desert", "#d8b46a"), ("steppe", "#c9c37c"), ("mediterranean scrub", "#92a35e"), ("temperate forest", "#4c8a4c"),
          ("subtropical / tropical forest", "#1f6e3d"), ("alpine", "#9c9288"), ("ice", "#eef5f8"), ("cold desert", "#b3aa94"), ("cold steppe", "#9ea58b"), ("boreal forest", "#3d6656"), ("tundra", "#a7b0a0")]
BIOME_KEYS = {"desert": 1, "steppe": 2, "mediterranean": 3, "temperate": 4, "subtropical": 5, "tropical": 5, "forest": 4, "alpine": 6, "ice": 7, "boreal": 10, "taiga": 10, "tundra": 11}

def biome_id(name):
    n = str(name).lower()
    for k, v in BIOME_KEYS.items():
        if n.startswith(k): return v
    raise ValueError(f"unknown biome '{name}'")

def lat_moisture(lat):
    """Latitude-band moisture 0..1: wet equator, trade-wind desert ~25, westerly belt ~50, dry above ~70."""
    return 0.5 + 0.35 * np.cos(2 * np.pi * np.abs(lat) / 50.0)

def prevailing_up(lat, mode):
    """+1: air arrives from the east (upwind = higher lon), -1: from the west."""
    a = np.abs(lat)
    if mode == "W": return -np.ones_like(lat)
    if mode == "E": return np.ones_like(lat)
    return np.where((a >= 30) & (a < 60), -1.0, 1.0)

def land_mask_poly(polys, g):
    m = Image.new("L", (g.nx, g.ny), 0); md = ImageDraw.Draw(m)
    for r in polys: md.polygon([((lo - g.lon0) / g.d, (g.ny - 1) - (la - g.lat0) / g.d) for lo, la in r], fill=255)
    return np.flipud(np.asarray(m) > 127)

def compute_biome(g, E, RID, sl, polys, cfg, R, log):
    bc = cfg.get("biome", {})
    res = bc.get("res_deg", 0.1); f = max(1, int(round(res / g.d))); res = f * g.d
    aridity = bc.get("aridity", 0.0)
    rel = relief_field(E, RID, sl, g, cfg, R)
    land = land_mask_poly(polys, g).astype(np.float32)
    elev = np.clip(rel, 0, None) * bc.get("elev_scale_m", 1500.0) * land
    ny2, nx2 = g.ny // f, g.nx // f
    def red(a): return a[:ny2 * f, :nx2 * f].reshape(ny2, f, nx2, f).mean(axis=(1, 3))
    lf = red(land); ev = red(elev) / np.maximum(lf, 1e-3)          # land-weighted mean elevation (m), rows south->north
    lat_c = g.lat0 - g.d / 2 + (np.arange(ny2) + 0.5) * res
    lat2 = lat_c[:, None] * np.ones((1, nx2))
    lon_c = g.lon0 - g.d / 2 + (np.arange(nx2) + 0.5) * res
    za = None
    zones = bc.get("zones", [])
    if zones:
        zi = Image.new("L", (nx2, ny2), 0); zd = ImageDraw.Draw(zi)
        for z in zones:
            zd.polygon([((lo - (g.lon0 - g.d / 2)) / res, ny2 - (la - (g.lat0 - g.d / 2)) / res) for lo, la in z["ring"]], fill=biome_id(z["biome"]))
        za = np.flipud(np.asarray(zi))
    spec = cfg.get("climate_spec")
    cls, info = climate_classify(lat_c, lon_c, ev, lf, -aridity, bc, R, res, spec, override=za)   # v2 config aridity is a moisture offset
    Mlat = lat_moisture(lat2)
    # mask: land + 1 cell dilation (editor clips to the coastline)
    lm = lf > 0
    d = lm.copy(); d[1:] |= lm[:-1]; d[:-1] |= lm[1:]; d[:, 1:] |= lm[:, :-1]; d[:, :-1] |= lm[:, 1:]
    out = np.where(d, cls, 0).astype(np.uint8)
    wt = lf * np.cos(np.radians(lat2))
    pct = {BIOMES[i - 1][0]: round(100 * float(wt[cls == i].sum() / wt.sum()), 1) for i in range(1, len(BIOMES) + 1)}
    cfg["_climate_report"] = climate_report(lat_c, lon_c, cls, info, lf, spec, res)
    sa = (lat2 <= bc.get("south_arc_max_lat", 8.0)) & (lf > 0); w2 = wt[sa]
    south = dict(land_km2_cells=int(sa.sum()), pct_humid_band_by_latitude=round(100 * float(w2[(Mlat[sa] >= 0.6)].sum() / max(w2.sum(), 1e-9)), 1),
                 pct_forest_or_scrub=round(100 * float(w2[np.isin(cls[sa], [3, 4, 5])].sum() / max(w2.sum(), 1e-9)), 1),
                 pct_desert_steppe=round(100 * float(w2[np.isin(cls[sa], [1, 2])].sum() / max(w2.sum(), 1e-9)), 1))
    flat = np.flipud(out).ravel()                                  # rows north -> south
    brk = np.nonzero(np.diff(flat))[0] + 1; st = np.concatenate([[0], brk]); ln = np.diff(np.concatenate([st, [len(flat)]]))
    rle = np.stack([flat[st], ln], axis=1).ravel().astype(int).tolist()
    bbox = [round(g.lon0 - g.d / 2, 4), round(g.lat0 - g.d / 2, 4), round(g.lon0 - g.d / 2 + nx2 * res, 4), round(g.lat0 - g.d / 2 + ny2 * res, 4)]
    blob = dict(bbox=bbox, res_deg=round(res, 4), w=nx2, h=ny2, order="rows north to south, [class, run] pairs",
                classes=[dict(id=i + 1, name=n, color=c) for i, (n, c) in enumerate(BIOMES)], rle=rle, aridity=aridity,
                zones=[z.get("name", z["biome"]) for z in zones], pct_land=pct, south_arc=south)
    return blob, out, pct, south

# ----------------------------------------------------------------------------- shared climate (v2 + v3, 2026-10-07 round 17)
PIN_TABLE = ["pin-capital", "pin-villa", "pin-bastion", "pin-threshold", "pin-weir", "pin-majidtown", "pin-emerald", "pin-gannet", "pin-ironpeaks", "pin-saltflats", "pin-spire"]

def climate_spec(state):
    """state.gen spire/first_cold/climate + resolved pin coordinates -> plain dict (goes into cfg; JSON-safe)."""
    gen = state.get("gen", {}) or {}
    pins = {f["id"]: f for f in state.get("features", []) if f.get("type") == "pin"}
    out = {k: dict(gen.get(k) or {}) for k in ("spire", "first_cold", "climate")}
    for k in ("spire", "first_cold"):
        pid = out[k].get("pin")
        if pid in pins: out[k]["lon"], out[k]["lat"] = pins[pid]["geometry"]["lon"], pins[pid]["geometry"]["lat"]
    out["pins"] = [dict(id=i, name=pins[i].get("name", i), lon=pins[i]["geometry"]["lon"], lat=pins[i]["geometry"]["lat"]) for i in PIN_TABLE if i in pins]
    return out

def _gc_km(LO, LA, lon, lat, R):
    P = sph_xyz(LO.ravel(), LA.ravel()); c = sph_xyz(np.array([lon]), np.array([lat]))[:, 0]
    return (np.arccos(np.clip(c @ P, -1, 1)) * R).reshape(LO.shape)

def spire_height(d_km, sp):
    """Spire relief (m) at great-circle distance d_km. profile "ease_in" (Matt 2026-10-11: wide base, rising slowly then sharply):
    height_m * (1 - d/base_km)^power inside base_km, 0 beyond. Legacy: height_m * exp(-d/profile_km)."""
    H = sp.get("height_m", 6000.0)
    if sp.get("profile") == "ease_in":
        return H * np.clip(1.0 - d_km / sp.get("base_km", 18.0), 0.0, 1.0) ** sp.get("power", 4.0)
    return H * np.exp(-d_km / sp.get("profile_km", 6.0))

def spire_extent_km(sp):
    """Distance beyond which the Spire adds no relief (for 'highest non-Spire point' style exclusions)."""
    return sp.get("base_km", 18.0) if sp.get("profile") == "ease_in" else 3 * (sp.get("profile_km", 6.0) + 10)

def climate_classify(lat_c, lon_c, ev, lf, dry, bc, R, bres, spec, override=None):
    """Biome grid rows south->north. ev: land elevation m; lf: land fraction; dry: moisture offset (>0 drier).
    Returns cls (uint8), dict(T, M, ev, cont) — T in C at the cell's (adjusted) elevation."""
    spec = spec or {}; sp = spec.get("spire", {}); fc = spec.get("first_cold", {}); cl = spec.get("climate", {})
    ny2, nx2 = ev.shape
    LO, LA = np.meshgrid(lon_c, lat_c); lat2 = LA
    # relief: soft cap for every peak, then the Spire on top (relief metres only; the coast never sees this)
    cap = sp.get("relief_cap_m")
    if cap:
        c0 = 0.8 * cap; ev = np.where(ev > c0, c0 + (cap - c0) * np.tanh((ev - c0) / (cap - c0)), ev)
    d_sp = None
    if "lon" in sp:
        d_sp = _gc_km(LO, LA, sp["lon"], sp["lat"], R)
        ev = np.where(lf > 0, np.maximum(ev, spire_height(d_sp, sp)), ev)
    sin2 = np.sin(np.radians(lat2)) ** 2
    T = bc.get("t_equator", 27.0) - bc.get("t_pole_drop", 28.0) * sin2 - bc.get("lapse_c_per_km", 6.5) * ev / 1000.0
    # prevailing wind: westerlies 30-60, easterlies elsewhere; blended over +-wind_blend_deg so there's no ruler line at 30/60
    a = np.abs(lat_c); b = cl.get("wind_blend_deg", 4.0)
    mode = bc.get("wind", "auto")
    wW = np.ones_like(a) if mode == "W" else np.zeros_like(a) if mode == "E" else smooth01((a - 30 + b) / (2 * b)) * (1 - smooth01((a - 60 + b) / (2 * b)))
    # continentality: km since the last sea cell, scanning from the upwind side
    kmc = (R * math.radians(bres) * np.cos(np.radians(lat_c)))[:, None]
    sea = lf < 0.5; ii = np.arange(nx2)[None, :] * np.ones((ny2, 1), int)
    lastW = np.maximum.accumulate(np.where(sea, ii, -1), axis=1)                     # wind from the west: sea at j <= i
    fW = np.where(lastW < 0, ii + 1, ii - lastW) * kmc
    kk = np.arange(nx2)[None, :]; lastR = np.maximum.accumulate(np.where(sea[:, ::-1], kk, -1), axis=1)   # wind from the east: sea at j >= i
    fE = np.where(lastR < 0, kk + 1, kk - lastR)[:, ::-1] * kmc
    fetch = wW[:, None] * fW + (1 - wW[:, None]) * fE
    cont = (1 - np.exp(-fetch / cl.get("cont_L_km", 1500.0))) * (lf > 0) if cl else np.zeros_like(ev)
    T = T - cl.get("cont_t_c", 0.0) * cont * np.clip(np.abs(lat2) / 50.0, 0.2, 1.2)
    if "lon" in fc:
        T = T - fc.get("amp_c", 0.0) * np.exp(-0.5 * (_gc_km(LO, LA, fc["lon"], fc["lat"], R) / fc.get("sigma_km", 400.0)) ** 2)
    Mlat = lat_moisture(lat2); M0 = np.clip(Mlat - dry, 0, 1)
    up = np.where(wW >= 0.5, -1.0, 1.0)[:, None]
    cols = np.arange(nx2)[None, :] * np.ones((ny2, 1), int); rows = np.arange(ny2)[:, None] * np.ones((1, nx2), int)
    nS = int(bc.get("shadow_range_km", 250.0) / (R * math.radians(bres))) + 1; Ls = bc.get("shadow_decay_km", 160.0)
    bar = np.zeros_like(ev); upm = np.zeros_like(ev); nl = 6
    for s in range(1, nS + 1):
        ci = (cols + (up * s).astype(int)); ok = (ci >= 0) & (ci < nx2)
        eu = np.where(ok, ev[rows, np.clip(ci, 0, nx2 - 1)], 0.0)
        bar = np.maximum(bar, eu * np.exp(-s * kmc / Ls))
        if s <= nl: upm += eu / nl
    shadow = np.clip((bar - ev) / bc.get("shadow_ref_m", 900.0), 0, 1)
    lift = np.clip((ev - upm) / bc.get("lift_ref_m", 600.0), 0, 1) * (lf > 0)
    M = M0 * (1 - bc.get("shadow_k", 0.55) * shadow) + bc.get("lift_k", 0.32) * lift + bc.get("elev_k", 0.12) * np.clip(ev / 2500.0, 0, 1)
    M = M * (1 - cl.get("cont_m", 0.0) * cont) * (1 + cl.get("pet_k", 0.0) * np.clip(cl.get("pet_ref_c", 18.0) - T, 0, None))
    M = np.clip(M, 0, 1)
    th = bc.get("thr", {}); t_des, t_ste, t_wet = th.get("desert", 0.2), th.get("steppe", 0.38), th.get("wet", 0.58)
    cold = T < bc.get("cold_dry_t", 8.0)
    cls = np.where(M < t_des, np.where(cold, 8, 1), np.where(M < t_ste, np.where(cold, 9, 2), np.where(M < t_wet,
          np.where((T >= bc.get("med_t", 13.0)) & (T < 26), 3, np.where(T >= 26, 2, 4)), np.where(T < 18, 4, 5)))).astype(np.uint8)
    # cold: alpine/ice only above the treeline band; cold lowlands go boreal -> tundra; polar ice anywhere
    high = ev >= bc.get("alpine_min_m", 1200.0)
    cls = np.where(~high & (T < bc.get("boreal_t", 3.0)) & np.isin(cls, [3, 4, 5]), 10, cls)
    cls = np.where(~high & (T < bc.get("tundra_t", -5.0)), 11, cls)
    cls = np.where(high & (T < bc.get("alpine_t", 6.0)), 6, cls)
    cls = np.where((high & (T < bc.get("ice_t", -6.0))) | (T < bc.get("polar_ice_t", -12.0)), 7, cls).astype(np.uint8)
    if override is not None: cls = np.where(override > 0, override, cls).astype(np.uint8)
    if d_sp is not None: cls = np.where(d_sp <= sp.get("ice_radius_km", 15.0), 7, cls).astype(np.uint8)   # the First Cold's ice, regardless of realism
    return cls, dict(T=T, M=M, ev=ev, cont=cont, d_sp=d_sp)

def climate_report(lat_c, lon_c, cls, info, lf, spec, bres):
    """Pin table (T, biome, elevation) for the canon pins inside this grid + the highest non-Spire point."""
    rows = []
    for p in (spec or {}).get("pins", []):
        j = int(round((p["lat"] - lat_c[0]) / bres)); i = int(round((p["lon"] - lon_c[0]) / bres))
        if 0 <= j < len(lat_c) and 0 <= i < len(lon_c):
            rows.append(dict(pin=p["name"], T_c=round(float(info["T"][j, i]), 1), biome=BIOMES[cls[j, i] - 1][0] if cls[j, i] else "sea",
                             elev_m=round(float(info["ev"][j, i])), M=round(float(info["M"][j, i]), 2), cont=round(float(info["cont"][j, i]), 2)))
    out = dict(pin_climate=rows)
    land = lf > 0
    if land.any():
        ev = np.where(land, info["ev"], -1e9)
        if info.get("d_sp") is not None: ev_ns = np.where(info["d_sp"] > spire_extent_km(spec.get("spire", {})), ev, -1e9)
        else: ev_ns = ev
        j, i = np.unravel_index(int(np.argmax(ev_ns)), ev_ns.shape)
        out["highest_non_spire"] = dict(elev_m=round(float(ev_ns[j, i])), lon=round(float(lon_c[i]), 2), lat=round(float(lat_c[j]), 2))
        if info.get("d_sp") is not None and (info["d_sp"] < 50).any():
            out["spire_cell_max_m"] = round(float(np.where(info["d_sp"] < 50, ev, -1e9).max()))
            out["spire_ice_cells"] = int(((info["d_sp"] <= spec["spire"].get("ice_radius_km", 15)) & land).sum())
    return out

# ----------------------------------------------------------------------------- preview
def render_preview(path, g, E, rel, sl, polys, intent, land_pins, water_pins, pin_ok, graticule, scale, title, R, biome=None):
    W, Hh = g.nx * scale, g.ny * scale
    m = Image.new("L", (W, Hh), 0); md = ImageDraw.Draw(m)
    def px(ring): return [((lo - g.lon0) / g.d * scale, (g.ny - 1 - (la - g.lat0) / g.d) * scale) for lo, la in ring]
    for r in polys: md.polygon(px(r), fill=255)
    mask = np.asarray(m) > 127
    def up(a): return np.asarray(Image.fromarray(np.flipud(a).astype(np.float32), "F").resize((W, Hh), Image.BICUBIC))
    rl = up(rel); El = up(E - sl)
    gy, gx = np.gradient(rl * 40.0)
    shade = np.clip(0.5 + 0.9 * (-gx - gy) / 1.4, 0.0, 1.0)
    h = np.clip(rl / 1.8, 0, 1)
    low = np.array([0.40, 0.55, 0.30]); mid = np.array([0.72, 0.64, 0.42]); high = np.array([0.55, 0.42, 0.34]); top = np.array([0.95, 0.94, 0.92])
    t = h[..., None]
    land_c = np.where(t < 0.35, low + (mid - low) * (t / 0.35), np.where(t < 0.7, mid + (high - mid) * ((t - 0.35) / 0.35), high + (top - high) * ((t - 0.7) / 0.3)))
    if biome is not None:
        f = int(round(biome["res_deg"] / g.d)); pal = np.array([[0, 0, 0]] + [[int(c["color"][i:i + 2], 16) for i in (1, 3, 5)] for c in biome["classes"]], np.uint8)
        arr = np.flipud(biome["cls"]); bimg = Image.fromarray(pal[arr], "RGB").resize((arr.shape[1] * f * scale, arr.shape[0] * f * scale), Image.BILINEAR)
        full = np.zeros((Hh, W, 3), np.float32) + np.array([0.5, 0.5, 0.5], np.float32)
        yo = (g.ny - arr.shape[0] * f) * scale; ba = np.asarray(bimg, np.float32) / 255
        full[yo:yo + ba.shape[0], :ba.shape[1]] = ba[:Hh - yo, :W]
        land_c = full * 0.9 + land_c * 0.1
        land_c = land_c * (0.62 + 0.55 * shade[..., None])
    else:
        land_c = land_c * (0.55 + 0.7 * shade[..., None])
    depth = np.clip(-El / 2.0, 0, 1)[..., None]
    sea_c = np.array([0.55, 0.78, 0.88]) * (1 - depth) + np.array([0.08, 0.20, 0.36]) * depth
    img = np.where(mask[..., None], land_c, sea_c)
    im = Image.fromarray((np.clip(img, 0, 1) * 255).astype(np.uint8), "RGB")
    d = ImageDraw.Draw(im, "RGBA")
    try: font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 13 * scale // 2 + 4)
    except Exception: font = ImageFont.load_default()
    step = graticule
    lo = math.ceil(g.lon0 / step) * step
    while lo <= g.lon0 + (g.nx - 1) * g.d:
        x = (lo - g.lon0) / g.d * scale; d.line([(x, 0), (x, Hh)], fill=(255, 255, 255, 60), width=1); d.text((x + 3, 3), f"{lo:g}", fill=(255, 255, 255, 200), font=font); lo += step
    la = math.ceil(g.lat0 / step) * step
    while la <= g.lat0 + (g.ny - 1) * g.d:
        y = (g.ny - 1 - (la - g.lat0) / g.d) * scale; d.line([(0, y), (W, y)], fill=(255, 255, 255, 60), width=1); d.text((3, y + 2), f"{la:g}", fill=(255, 255, 255, 200), font=font); la += step
    for r in intent:
        p = px(r); d.line(p + [p[0]], fill=(255, 255, 255, 110), width=1)
    for r in polys:
        p = px(r); d.line(p + [p[0]], fill=(20, 20, 20, 255), width=max(1, scale // 2))
    rad = 4 * scale // 2 + 2
    for p, ok in zip(land_pins, pin_ok[0]):
        x, y = px([[p["lon"], p["lat"]]])[0]
        d.ellipse([x - rad, y - rad, x + rad, y + rad], fill=(230, 40, 40, 255) if ok else (255, 230, 0, 255), outline=(255, 255, 255, 255))
    for p, ok in zip(water_pins, pin_ok[1]):
        x, y = px([[p["lon"], p["lat"]]])[0]
        d.rectangle([x - rad, y - rad, x + rad, y + rad], fill=(40, 120, 255, 255) if ok else (255, 230, 0, 255), outline=(255, 255, 255, 255))
    for p in land_pins + water_pins:
        if p.get("label"):
            x, y = px([[p["lon"], p["lat"]]])[0]; d.text((x + rad + 2, y - rad), p["label"], fill=(255, 255, 255, 255), font=font, stroke_width=1, stroke_fill=(0, 0, 0, 255))
    if biome is not None:
        lx, ly = 8, 40
        d.rectangle([lx - 4, ly - 4, lx + 205, ly + 20 * len(biome["classes"]) + 2], fill=(0, 0, 0, 140))
        for i, c in enumerate(biome["classes"]):
            col = tuple(int(c["color"][k:k + 2], 16) for k in (1, 3, 5))
            d.rectangle([lx, ly + i * 20, lx + 14, ly + i * 20 + 14], fill=col + (255,), outline=(255, 255, 255, 255))
            d.text((lx + 22, ly + i * 20 - 1), f"{c['name']} {biome['pct_land'].get(c['name'], 0):.0f}%", fill=(255, 255, 255, 255), font=font)
    d.rectangle([0, Hh - 22 * scale // 2 - 6, W, Hh], fill=(0, 0, 0, 150))
    d.text((6, Hh - 22 * scale // 2 - 2), title, fill=(255, 255, 255, 255), font=font)
    im.save(path)

# ----------------------------------------------------------------------------- main generate
def generate(cfg, log=print, prep=None):
    """prep (wild-group only): callable(g, R, log) -> dict(W, SDF, EXCL, s_km, bump_ok) — protection mask + midline exclusion."""
    t0 = time.time()
    R = cfg.get("radius_km", R_KM_DEFAULT)
    bb = cfg["bbox"]; res = cfg.get("res_deg", 0.05)
    g = Grid(bb, res); st = {**DEFAULT_STYLE, **cfg.get("style", {})}
    out_cfg = cfg.get("output", {})
    target = cfg["target_km2"]; tol = cfg.get("area_tol", 0.005)
    group = cfg.get("group", cfg.get("name", "region").lower())
    log(f"grid {g.nx}x{g.ny} ({g.nx*g.ny/1e6:.2f}M cells), target {target:,.0f} km2, seed {cfg['seed']}")
    E, RID = compute_field(g, cfg, R, log)
    log(f"field done {time.time()-t0:.1f}s")
    PR = prep(g, R, log) if prep else None
    if PR is not None and PR.get("EXCL") is not None:
        E = E - 4.0 * PR["EXCL"]   # nearer another group's land -> sea
    if PR is not None and PR.get("BASIN") is not None:
        E = E - PR["BASIN"]        # inland basins (lakes / inland seas)
    wrow = cell_area_weights(g, R)
    land_pins = cfg.get("land_pins", []); water_pins = cfg.get("water_pins", [])
    bump_l, bump_w = land_pins, water_pins
    if PR is not None:   # pins inside protection keep the intent coast; verified on the output
        bump_l = [p for p in land_pins if PR["bump_ok"](p)]; bump_w = [p for p in water_pins if PR["bump_ok"](p)]
    budget = cfg.get("file_budget_bytes", 1_500_000); bpv = 19.5; fracs = cfg.get("lod_fracs", [0.25, 0.0625])
    lodmul = 1 + sum(fracs)
    biome_reserve = cfg.get("biome_reserve_bytes", 160_000)
    max_verts = int((budget - biome_reserve - 20_000) / (bpv * lodmul))
    eps = cfg.get("simplify_eps_cells", 0.15) * g.d
    E0 = E.copy()
    for attempt in range(4):
        corr = 1.0; sl = None; full = None; hist = []; best = None; lo = hi = None   # lo/hi: corr brackets (area under / over target)
        for it in range(9):
            E = E0.copy()
            sl = apply_pins(E, g, bump_l, bump_w, wrow, target * corr, R, st, lambda s: None)
            if PR is not None and PR.get("CJ") is not None:
                cw = PR["CJ"]["w"]; E = (1 - cw) * E + cw * (sl + PR["CJ"]["target"])
            if PR is not None and PR.get("W") is not None:   # protected coast = intent coast (sdf = 0) whatever sl is
                W = PR["W"]; sd = PR["SDF"]; t = smooth01(np.abs(sd) / (3 * PR["s_km"]))
                mag = (1 - t) * np.abs(sd) / PR["s_km"] + t * np.maximum(np.abs(E - sl), 0.5)   # linear at the coast, wild relief inland
                if cfg.get("protect_only"):   # all land comes from protection zones: the wild field is sea everywhere (mirrored, keeps relief)
                    Ew = np.where(E > sl - 0.5, 2 * (sl - 0.5) - E, E)
                    E = (1 - W) * Ew + W * (sl + np.sign(sd) * mag)
                else:
                    E = (1 - W) * E + W * (sl + np.sign(sd) * mag)
            raw = land_rings(E, sl, g, R, cfg, land_pins)
            raw.sort(key=lambda r: -ring_area_km2(r, R))
            full, eps_used = simplify_full(raw, eps, max_verts)
            keep = [i for i, q in enumerate(full) if q is not None]
            full = [full[i] for i in keep]; raw = [raw[i] for i in keep]
            tot = sum(ring_area_km2(q, R) for q in full)
            raw_tot = sum(ring_area_km2(r, R) for r in raw)
            hist.append((corr, tot))
            log(f"iter {it}: corr={corr:.4f} sl={sl:.3f} polys={len(full)} verts={sum(len(q) for q in full)} eps={eps_used/g.d:.2f}c area={tot:,.0f} ({100*(tot/target-1):+.2f}%, raw {100*(raw_tot/target-1):+.2f}%) {time.time()-t0:.1f}s")
            err = abs(tot / target - 1)
            if best is None or err < best[0]: best = (err, E, sl, full, raw, eps_used, tot)
            if err <= tol: break
            if tot < target: lo = corr if lo is None else max(lo, corr)
            else: hi = corr if hi is None else min(hi, corr)
            # bracketed -> bisect (flat lowland flips the area across target under the multiplicative step)
            corr = math.sqrt(lo * hi) if lo is not None and hi is not None else corr * (target / tot) ** (1.0 if it < 3 else 0.9)
        if best[0] < abs(tot / target - 1):   # unconverged: keep the closest iteration, not the last
            _, E, sl, full, raw, eps_used, tot = best
            log(f"  area: kept closest iteration ({100*(tot/target-1):+.2f}%)")
        bad = {i for pair in find_crossings(full) for i in pair}
        for i in bad: full[i] = rnd_ring(raw[i])
        levels = build_lods(full, fracs, log)
        areas = [ring_area_km2(q, R) for q in full]
        area_tot = float(sum(areas))
        order = sorted(range(len(full)), key=lambda i: -areas[i])
        # pins on the final (full-res) generated land
        def inland(p):
            pt = np.array([[p["lon"], p["lat"]]])
            for q in full:
                if q[:, 0].min() <= p["lon"] <= q[:, 0].max() and q[:, 1].min() <= p["lat"] <= q[:, 1].max() and point_in_ring(pt, q)[0]: return True
            return False
        lp = [inland(p) for p in land_pins]; wp = [not inland(p) for p in water_pins]
        biome = None; bpct = {}; south = {}
        if cfg.get("biome"):
            blob, bcls, bpct, south = compute_biome(g, E, RID, sl, full, cfg, R, log)
            biome = blob; biome["cls"] = bcls
        gen = dict(name="coastgen", version=VERSION, seed=cfg["seed"], params_hash=params_hash(cfg), region=cfg.get("name", "region"), bbox=bb, res_deg=res,
                   radius_km=R, target_km2=round(target), achieved_km2=round(area_tot), simplify_eps_deg=round(eps_used, 5),
                   min_island_km2=cfg.get("min_island_km2", 90.0), min_island_rules=cfg.get("min_island_rules", []), style=st, arcs=cfg.get("arcs", []),
                   lod_vertices=[sum(len(r) for r in lv) for lv in levels])
        feats = []
        for n_, i in enumerate(order):
            q = [lv[i] for lv in levels]
            feats.append({"id": f"gen-{group}-{n_+1:03d}", "type": "land", "kind": "generated", "name": f"{cfg.get('label','Generated')} island {n_+1:03d}",
                          "status": "proposed", "group": group, "generator": "coastgen", "seed": cfg["seed"], "params": gen["params_hash"], "area_km2": round(areas[i]),
                          "bbox": [round(float(full[i][:, 0].min()), 4), round(float(full[i][:, 1].min()), 4), round(float(full[i][:, 0].max()), 4), round(float(full[i][:, 1].max()), 4)],
                          "lod": [[[float(x), float(y)] for x, y in r] for r in q]})
        result = dict(version=2, group=group, label=cfg.get("label", group), generator=gen, features=feats)
        if biome is not None:
            bj = {k: v for k, v in biome.items() if k != "cls"}; result["biome"] = bj
        text = json.dumps(result, separators=(",", ":"))
        size = len(text.encode())
        log(f"file {size/1e6:.3f} MB (budget {budget/1e6:.2f}) verts/LOD {gen['lod_vertices']} attempt {attempt}")
        if size <= budget: break
        max_verts = int(max_verts * budget / size * 0.93)
    lakes = []
    if cfg.get("basins"):
        for n_, lr in enumerate(sorted(lake_rings(E, sl, g, R, cfg), key=lambda r: -ring_area_km2(r, R))):
            q0 = rnd_ring(lr[dp_tol(lr, eps_used)]) if len(lr) > 8 else rnd_ring(lr)
            q1 = rnd_ring(lr[dp_tol(lr, eps_used * 4)]) if len(lr) > 16 else q0
            lakes.append({"id": f"gen-{group}-lake-{n_+1:03d}", "type": "lake", "kind": "lake", "name": f"{cfg.get('label','Generated')} lake {n_+1:03d}",
                          "group": group, "area_km2": round(ring_area_km2(lr, R)),
                          "bbox": [round(float(lr[:, 0].min()), 4), round(float(lr[:, 1].min()), 4), round(float(lr[:, 0].max()), 4), round(float(lr[:, 1].max()), 4)],
                          "lod": [[[float(x), float(y)] for x, y in q] for q in (q0, q1, q1)]})
        result["features"].extend(lakes)
        log(f"lakes: {len(lakes)} ({sum(l['area_km2'] for l in lakes):,} km2)")
    stats = dict(runtime_s=round(time.time() - t0, 1), peak_mb=round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6),
                 islands=len(full), file_bytes=size, lod_vertices=gen["lod_vertices"], area_km2=round(area_tot), target_km2=round(target),
                 area_err_pct=round(100 * (area_tot / target - 1), 2), land_pins_ok=f"{sum(lp)}/{len(lp)}", water_pins_ok=f"{sum(wp)}/{len(wp)}",
                 crossing_fallbacks=len(bad), self_intersections_full=len({i for p in find_crossings(full) for i in p}),
                 failed_pins=[p.get("label", "?") for p, ok in zip(land_pins, lp) if not ok] + [p.get("label", "?") for p, ok in zip(water_pins, wp) if not ok],
                 biome_pct=bpct, south_arc=south, lakes=[(l["name"], l["area_km2"]) for l in lakes])
    stats.update(cfg.pop("_climate_report", {}) or {})
    result["generator"]["stats"] = {k: v for k, v in stats.items() if k not in ("biome_pct", "south_arc")}
    text = json.dumps(result, separators=(",", ":")); stats["file_bytes"] = len(text.encode())
    if out_cfg.get("json"):
        tmp = out_cfg["json"] + ".tmp"; open(tmp, "w").write(text); os.replace(tmp, out_cfg["json"])
    if out_cfg.get("png"):
        rel = relief_field(E, RID, sl, g, cfg, R)
        title = (f"{cfg.get('label','')} seed {cfg['seed']} | {stats['islands']} islands, {stats['lod_vertices'][0]} verts, {area_tot/1e6:.2f}M km2 ({stats['area_err_pct']:+.2f}%) | "
                 f"pins L{stats['land_pins_ok']} W{stats['water_pins_ok']} | {size/1e6:.2f} MB")
        intent = [p["ring"] for p in cfg.get("intent", [])]
        render_preview(out_cfg["png"], g, E, rel, sl, full, intent, land_pins, water_pins, (lp, wp), cfg.get("graticule_deg", 5), out_cfg.get("scale", 1), title, R, biome)
        if lakes:
            im = Image.open(out_cfg["png"]).convert("RGB"); dr = ImageDraw.Draw(im); sc = out_cfg.get("scale", 1)
            for l in lakes:
                dr.polygon([((lo - g.lon0) / g.d * sc, (g.ny - 1 - (la - g.lat0) / g.d) * sc) for lo, la in l["lod"][0]], fill=(62, 120, 168), outline=(20, 30, 45))
            im.save(out_cfg["png"])
    stats["runtime_s"] = round(time.time() - t0, 1)
    stats["peak_mb"] = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6)
    if PR is not None and PR.get("report"):
        stats.update(PR["report"](full, lp, wp))
        lk = [np.asarray(l["lod"][0], float) for l in lakes]
        stats["land_pins_in_lakes"] = [p.get("label", "?") for p in land_pins if any(point_in_ring(np.array([[p["lon"], p["lat"]]]), r)[0] for r in lk)]
    return result, stats

# ----------------------------------------------------------------------------- Western Isles preset
WI_PRESET = dict(
    name="WI", label="Western Isles",
    bbox=[-48.0, -4.0, 8.0, 46.0], res_deg=0.05,
    style=dict(),
    land_pins=[], water_pins=[],
    arcs=[
        dict(pts=[[-47,12],[-39,6.5],[-29,4],[-18,3.5],[-9,6],[-3,10],[3,17]], width_km=48, amp=1.5, bead_km=100),
        dict(pts=[[-5,44],[1,37],[1.5,29],[-2,21],[-6.5,16.5]], width_km=45, amp=1.5, bead_km=100),
        dict(pts=[[-47,38],[-40,33],[-31.6,29.6],[-27,23],[-25,17],[-29,11]], width_km=45, amp=1.4, bead_km=100),
        dict(pts=[[-30,45],[-22,41],[-14,38.5],[-6,40]], width_km=40, amp=1.2, bead_km=90),
        dict(pts=[[-19.5,28.9],[-12,24],[-9,19]], width_km=30, amp=1.0, bead_sigma=0, taper_km=50),
        dict(pts=[[-19.5,28.9],[-24,33],[-26,38]], width_km=28, amp=0.9, bead_sigma=0, taper_km=50),
    ],
    canyons=[dict(pts=[[-19.5,29.0],[-20.4,31],[-19.6,33],[-21,35.5]], width_km=14, depth=0.6)],
    min_island_km2=90.0, group="wi",
)

def build_wi_cfg(state_path):
    s = json.load(open(state_path))
    R = s["planet"]["radius_km"]
    cfg = copy.deepcopy(WI_PRESET); cfg["radius_km"] = R
    lands = [f for f in s["features"] if f["type"] == "land" and f["name"].startswith("Western Isles")]
    areas = [(ring_area_km2(f["geometry"], R), f) for f in lands]
    areas.sort(key=lambda t: -t[0])
    cfg["target_km2"] = sum(a for a, _ in areas)
    cfg["intent"] = [dict(ring=f["geometry"], weight=1.0 if a > 150000 else 0.6, name=f["name"]) for a, f in areas]
    bb = cfg["bbox"]
    for f in s["features"]:
        if f["type"] != "pin": continue
        g = f["geometry"]
        if bb[0] < g["lon"] < bb[2] and bb[1] < g["lat"] < bb[3]:
            p = dict(lon=g["lon"], lat=g["lat"], label=f["name"])
            cfg["water_pins" if f["kind"] == "sea" else "land_pins"].append(p)
    return cfg

# ----------------------------------------------------------------------------- v2 wild on any group (2026-10-07): protection mask + midline exclusion
# Wild = v2 field (noise-dominant, intent as a prior) on one landmass group, for baking into intents. Inside protected genzones
# (`protect`: strength 0..1, or wander_km <= PROTECT_WANDER_KM) the field is blended toward the exact signed distance to the intent
# coast, so the drawn coast survives verbatim; cells nearer another group's land than this group's are forced to sea.
PROTECT_WANDER_KM = 10.0
PROTECT_S_KM = 25.0          # protected field slope: 1 E unit per 25 km
EXCL_FEATHER_KM = 40.0

def edt_km(mask, g, R):
    """Distance (km) to the nearest True cell on a lon/lat grid; x spacing per row = R*d*cos(lat) (local-flat metric)."""
    sy = R * math.radians(g.d); sx = sy * np.cos(np.radians(g.lat()))[:, None]
    f = _edt1d_sq(np.where(mask, 0.0, INF_SQ))
    f = np.where(f >= INF_SQ, INF_SQ, f * sx ** 2 / sy ** 2)
    f = _edt1d_sq(np.ascontiguousarray(f.T)).T * sy ** 2
    return np.sqrt(np.minimum(f, INF_SQ)).astype(np.float32)

def _ll_mask(g, rings):
    return land_mask_poly([np.asarray(r, float) for r in rings], g)

def _km_local(lon, lat, lonc, latc, R):
    k = R * math.pi / 180
    return np.stack([(np.asarray(lon) - lonc) * math.cos(math.radians(latc)) * k, (np.asarray(lat) - latc) * k], -1)

def _near_dist(Q, P, chunk=4096):
    """min distance from each Q (n,2) to the point cloud P (m,2), km."""
    out = np.empty(len(Q), np.float32)
    for a in range(0, len(Q), chunk):
        d = Q[a:a + chunk, None, :] - P[None]; out[a:a + chunk] = np.sqrt((d * d).sum(-1)).min(1)
    return out

def near_dist_ll(q, p, R, pad_deg=0.6, chunk=512):
    """min great-circle-ish (chord) distance, km, from lon/lat points q (n,2) to lon/lat points p (m,2); spatially chunked.
    Distances beyond ~pad_deg*R*pi/180 are capped at that (only near-coast values matter)."""
    out = np.full(len(q), pad_deg * R * math.pi / 180, np.float32)
    if not len(q) or not len(p): return out
    o = np.lexsort((q[:, 1], q[:, 0])); P3 = sph_xyz(p[:, 0], p[:, 1]).T
    for a in range(0, len(o), chunk):
        idx = o[a:a + chunk]; qq = q[idx]
        lo, hi = qq.min(0) - pad_deg, qq.max(0) + pad_deg
        sel = (p[:, 0] >= lo[0]) & (p[:, 0] <= hi[0]) & (p[:, 1] >= lo[1]) & (p[:, 1] <= hi[1])
        if not sel.any(): continue
        Q3 = sph_xyz(qq[:, 0], qq[:, 1]).T; pp = P3[sel]
        for b in range(0, len(pp), 8192):
            d = Q3[:, None, :] - pp[None, b:b + 8192]; out[idx] = np.minimum(out[idx], np.sqrt((d * d).sum(-1)).min(1) * R)
    return out

def group_rings_all(state):
    gr = {}
    for f in state["features"]:
        if f.get("type") == "land" and f.get("kind") == "intent":
            for t in f.get("groups") or []: gr.setdefault(t, []).append(np.asarray(f["geometry"], float))
    return gr

def build_group_wild_cfg(state_path, tag, pad=(5.0, 4.0)):
    s = json.load(open(state_path)); R = float(s["planet"]["radius_km"])
    feats = group_intent(s, tag)
    if not feats: raise SystemExit(f"no intent land carries group '{tag}'")
    areas = sorted(((ring_area_km2(f["geometry"], R), f) for f in feats), key=lambda t: -t[0])
    allp = np.vstack([np.asarray(f["geometry"], float) for f in feats])
    bb = [round(float(allp[:, 0].min()) - pad[0], 2), round(float(allp[:, 1].min()) - pad[1], 2),
          round(float(allp[:, 0].max()) + pad[0], 2), round(float(allp[:, 1].max()) + pad[1], 2)]
    cfg = dict(name=tag, label=V3_LABELS.get(tag, tag), group=tag, radius_km=R, bbox=bb, res_deg=0.05, style={}, min_island_km2=90.0,
               target_km2=sum(a for a, _ in areas),
               intent=[dict(ring=f["geometry"], weight=1.0 if a > 150000 else 0.6, name=f["name"], id=f["id"]) for a, f in areas],   # config intent_small_weight overrides 0.6
               arcs=[dict(pts=f["geometry"], **(f.get("params") or {})) for f in s["features"]
                     if f.get("type") == "line" and f.get("kind") == "arc" and tag in (f.get("groups") or [])])
    lonc, latc = 0.5 * (bb[0] + bb[2]), 0.5 * (bb[1] + bb[3])
    pj = Proj(lonc, latc, R)
    pins, notes = assign_pins(s, tag, pj, group_rings_all(s), 0.0)
    lp = {(p["lon"], p["lat"]): p for p in pins["land"]}
    cfg["land_pins"] = [dict(lon=p["lon"], lat=p["lat"], label=p["name"]) for p in lp.values()]
    cfg["water_pins"] = [dict(lon=p["lon"], lat=p["lat"], label=p["name"]) for p in pins["water"]]
    cfg["coastal_pins"] = [dict(lon=p["lon"], lat=p["lat"], label=p["name"]) for p in pins["coastal"]]
    zones = []
    for z in group_zones(s, tag, {}):
        st_ = z.get("protect"); wk = (z.get("params") or {}).get("wander_km")
        strength = float(st_) if st_ not in (None, False) else (1.0 if (wk is not None and wk <= PROTECT_WANDER_KM) else 0.0)
        if strength > 0: zones.append(dict(name=z["name"], ring=z["geometry"], feather_km=float(z.get("feather_km", 40.0)), strength=strength, wander_km=wk))
    cfg["protect_zones"] = zones
    return cfg, s, notes

def _detail_scale(area_km2, pd):
    """protect_detail amplitude scale by island size: sqrt(area / ref_km2) clipped to [min_scale, 1] (big coasts get big bays, islets survive)."""
    if not pd or not pd.get("ref_km2"): return 1.0
    return float(min(1.0, max(pd.get("min_scale", 0.3), math.sqrt(area_km2 / pd["ref_km2"]))))

def apply_protect_overrides(cfg):
    """config `protect_overrides`: {zone name: {feather_km, strength}} (wild runs only; genzone feathers stay as drawn)."""
    ov = cfg.get("protect_overrides") or {}
    for z in cfg.get("protect_zones", []):
        z.update(ov.get(z["name"], {}))

def make_wild_prep(cfg, state):
    """Returns prep(g, R, log) for generate(): W (protection), SDF (exact, km, + inland, only where W > 0), EXCL, report()."""
    tag = cfg["group"]; gr = group_rings_all(state)
    own = gr.get(tag, []); others = [r for t, rs in gr.items() if t != tag for r in rs]
    def prep(g, R, log):
        t0 = time.time(); out = {}
        LON, LAT = np.meshgrid(g.lon(), g.lat())
        # exclusion: cells nearer another group's land than this group's
        if others:
            d_own = edt_km(_ll_mask(g, own), g, R); d_oth = edt_km(_ll_mask(g, others), g, R)
            out["EXCL"] = smooth01((d_own - d_oth) / EXCL_FEATHER_KM + 0.5).astype(np.float32); del d_own, d_oth
        # protection weight
        W = np.zeros((g.ny, g.nx), np.float32); zinfo = []; AMP = np.zeros((g.ny, g.nx), np.float32); DNEAR = np.full((g.ny, g.nx), np.inf, np.float32); DSEC = np.full((g.ny, g.nx), np.inf, np.float32)
        own_m = _ll_mask(g, own)
        SDF = (edt_km(~own_m, g, R) - edt_km(own_m, g, R)).astype(np.float32)   # raster signed distance, + inland
        for z in cfg.get("protect_zones", []):
            m = _ll_mask(g, [z["ring"]]); d = edt_km(m, g, R)
            wz = z["strength"] * smooth01(1.0 - d / max(z["feather_km"], 1.0)); W = np.maximum(W, wz.astype(np.float32))
            sf = z.get("sea_feather_km", 0.0)
            if sf > z["feather_km"]:   # intent sea near the zone can't fill (keeps straits/bays open); land there stays wild
                ws = z["strength"] * smooth01(1.0 - (d - z["feather_km"]) / (sf - z["feather_km"])) * smooth01(-SDF / 20.0); W = np.maximum(W, ws.astype(np.float32))
            zinfo.append((z, m))
            if "detail_scale" in z:   # nearest zone's size scale wins (islets beside a big island keep theirs)
                nr = d < DNEAR; DSEC = np.where(nr, DNEAR, np.minimum(DSEC, d)); AMP[nr] = z["detail_scale"]; DNEAR[nr] = d[nr]
        coast = {}
        if W.max() > 0:
            ys, xs = np.nonzero((W > 0) & (np.abs(SDF) < 20.0))   # exact distance near the coast only
            qll = np.stack([LON[ys, xs], LAT[ys, xs]], 1)
            pts = np.vstack([densify_ll(r, 0.004) for r in own])   # ~0.3 km spacing
            dist = near_dist_ll(qll, pts, R)
            inside = np.zeros(len(qll), bool)
            for r in own:
                bx = (qll[:, 0] >= r[:, 0].min()) & (qll[:, 0] <= r[:, 0].max()) & (qll[:, 1] >= r[:, 1].min()) & (qll[:, 1] <= r[:, 1].max())
                if bx.any(): inside[bx] ^= point_in_ring(qll[bx], r)
            SDF[ys, xs] = np.where(inside, dist, -dist)
            coast = dict(pts=pts)
            pd = cfg.get("protect_detail")
            if pd:   # crenellate protected coasts: fractal offset of the drawn coast (shape kept, detail added)
                ys, xs = np.nonzero((W > 0) & (np.abs(SDF) < 3 * pd.get("amp_km", 4.0) + 10.0))
                P = sph_xyz(LON[ys, xs], LAT[ys, xs])
                n = fbm3(P, R / pd.get("wl_km", 40.0), pd.get("oct", 5), pd.get("gain", 0.55), int(cfg.get("seed", 7)) + 4243)
                sc = np.where(AMP[ys, xs] > 0, AMP[ys, xs], 1.0)   # per-zone size scale (zones without one: 1)
                disp = pd.get("amp_km", 4.0) * sc * n
                ck = pd.get("channel_km")
                if ck:   # channel guard: no land growth near the midline between two islands (straits stay open)
                    gap = DSEC[ys, xs] - DNEAR[ys, xs]
                    disp = np.where(disp > 0, disp * smooth01(np.nan_to_num(gap, posinf=1e6) / ck), disp)
                SDF[ys, xs] += disp.astype(np.float32)
                log(f"protect detail: amp {pd.get('amp_km', 4.0)} km, wl {pd.get('wl_km', 40.0)} km, cells {len(ys)}")
        out.update(W=W if W.max() > 0 else None, SDF=SDF, s_km=PROTECT_S_KM)
        bs = cfg.get("basins")
        if bs:   # low-frequency basin noise, only deep inland and never in protected zones -> enclosed lakes / inland seas
            pins_ll = {f["id"]: f["geometry"] for f in state["features"] if f.get("type") == "pin"}
            lo_f, la_f = LON.ravel(), LAT.ravel(); BAS = np.zeros((g.ny, g.nx), np.float32)
            for li, ly in enumerate([bs] + list(bs.get("layers", []))):   # layers stack (max); each may be limited to a lon range / kept clear of a pin
                seed = int(cfg.get("seed", 7)) + 911 + 37 * li; fq = R / ly.get("wl_km", 700.0)
                Bn = np.empty(g.ny * g.nx, np.float32)
                for a0 in range(0, len(Bn), 65536):
                    P = sph_xyz(lo_f[a0:a0 + 65536], la_f[a0:a0 + 65536]); Bn[a0:a0 + 65536] = fbm3(P, fq, ly.get("oct", 4), 0.5, seed)
                Bn = Bn.reshape(g.ny, g.nx); thr = ly.get("thr", 0.3)
                inland = smooth01((SDF - ly.get("inland_km", 250.0)) / ly.get("ramp_km", 150.0))
                m = ly.get("amp", 3.0) * np.clip((Bn - thr) / max(1e-3, 1 - thr), 0, 1) ** 0.7 * inland
                if "lon_range" in ly:
                    l0, l1 = ly["lon_range"]; fe = ly.get("lon_feather_deg", 3.0)
                    m = m * smooth01((LON - l0) / fe + 0.5) * smooth01((l1 - LON) / fe + 0.5)
                if "lat_range" in ly:
                    a0_, a1_ = ly["lat_range"]; fe = ly.get("lat_feather_deg", 2.0)
                    m = m * smooth01((LAT - a0_) / fe + 0.5) * smooth01((a1_ - LAT) / fe + 0.5)
                if ly.get("avoid_pin") in pins_ll:
                    pg = pins_ll[ly["avoid_pin"]]; m = m * smooth01((_gc_km(LON, LAT, pg["lon"], pg["lat"], R) - ly.get("avoid_km", 500.0)) / 150.0 + 0.5)
                BAS = np.maximum(BAS, m.astype(np.float32))
                log(f"basins[{li}]: wl {ly.get('wl_km', 700)} km, thr {thr}, amp {ly.get('amp', 3.0)}, active cells {int((m > 0.05).sum())}")
            clear = bs.get("pin_clear_km", 80.0)
            if clear and cfg.get("land_pins"):   # never drown a canon land pin
                pm = np.zeros((g.ny, g.nx), bool)
                for p_ in cfg["land_pins"]:
                    pm[int(round((p_["lat"] - g.lat0) / g.d)) % g.ny, int(round((p_["lon"] - g.lon0) / g.d)) % g.nx] = True
                BAS = BAS * smooth01((edt_km(pm, g, R) - clear) / 60.0 + 0.5)
            out["BASIN"] = (BAS * (W < 0.05)).astype(np.float32)
        cj = cfg.get("conjugate")
        if cj and gr.get(cj["group"]):   # conjugate margin: nudge our coast to run parallel to the other plate's facing coast
            dB = edt_km(_ll_mask(g, gr[cj["group"]]), g, R); reach = cj.get("reach_km", 700.0); band = cj.get("band_km", 200.0)
            near = (np.abs(SDF) < 6.0) & (dB < reach)
            gap = float(cj["gap_km"]) if cj.get("gap_km") else (float(np.median(dB[near])) if near.any() else None)
            if gap is not None:
                wc = cj.get("weight", 0.35) * smooth01(1.0 - np.abs(SDF) / band) * smooth01(1.0 - dB / reach)
                out["CJ"] = dict(w=wc.astype(np.float32), target=((dB - gap) / PROTECT_S_KM).astype(np.float32))
                log(f"conjugate margin vs {cj['group']}: gap {gap:.0f} km, band {band:.0f}, reach {reach:.0f}, weight {cj.get('weight', 0.35)}, cells {int((wc > 0.01).sum())}")
                out["conj_gap_km"] = gap
        def wat(lon, lat):
            ix = int(round((lon - g.lon0) / g.d)); iy = int(round((lat - g.lat0) / g.d))
            return float(W[min(max(iy, 0), g.ny - 1), min(max(ix, 0), g.nx - 1)])
        out["bump_ok"] = lambda p: wat(p["lon"], p["lat"]) <= 0.5
        def report(full, lp_ok, wp_ok):
            rep = {}
            # gate: generated coast vs intent coast inside each protected core (W >= 0.99)
            if coast:
                gz = []; gpts = np.vstack([densify_ll(q, 0.004) for q in full]); ipts = np.vstack([densify_ll(r, 0.004) for r in own])
                for z, m in zinfo:
                    ix = np.clip(np.round((gpts[:, 0] - g.lon0) / g.d).astype(int), 0, g.nx - 1); iy = np.clip(np.round((gpts[:, 1] - g.lat0) / g.d).astype(int), 0, g.ny - 1)
                    core = (W[iy, ix] >= 0.99 * z["strength"]) & m[iy, ix]
                    jx = np.clip(np.round((ipts[:, 0] - g.lon0) / g.d).astype(int), 0, g.nx - 1); jy = np.clip(np.round((ipts[:, 1] - g.lat0) / g.d).astype(int), 0, g.ny - 1)
                    gk = gpts[core]; ik = ipts[m[jy, jx]]
                    d1 = near_dist_ll(gk, ik, R, pad_deg=2.0) if len(gk) and len(ik) else np.zeros(1)
                    d2 = near_dist_ll(ik, gpts, R, pad_deg=2.0) if len(ik) else np.zeros(1)
                    dmax = float(max(d1.max(), d2.max())); gate = z.get("wander_km") or PROTECT_WANDER_KM
                    gz.append(dict(zone=z["name"], max_dev_km=round(dmax, 2), p95_km=round(float(np.percentile(np.concatenate([d1, d2]), 95)), 2), gate_km=gate, ok=dmax <= gate))
                rep["protect_gate"] = gz
            # land pins: on the generated component that overlaps the pin's own intent component
            bad = []
            for p in cfg.get("land_pins", []):
                pt = np.array([[p["lon"], p["lat"]]])
                src = [r for r in own if point_in_ring(pt, r)[0]]
                gen = [q for q in full if point_in_ring(pt, q)[0]]
                if not src or not gen: bad.append(f"{p['label']}: {'no source intent' if not src else 'at sea'}"); continue
                a_src = ring_area_km2(src[0], R); a_gen = ring_area_km2(gen[0], R)
                if a_gen < 0.3 * a_src: bad.append(f"{p['label']}: on a {a_gen:,.0f} km2 islet, source {a_src:,.0f} km2")
            rep["land_pin_component_issues"] = bad
            # coastal pins: distance to the generated coast
            cp = []
            if cfg.get("coastal_pins"):
                gpts = np.vstack([densify_ll(q, 0.01) for q in full])
                for p in cfg["coastal_pins"]:
                    k0 = _km_local(gpts[:, 0], gpts[:, 1], p["lon"], p["lat"], R)
                    cp.append(dict(pin=p["label"], coast_km=round(float(np.sqrt((k0 ** 2).sum(1)).min()), 1), W=round(wat(p["lon"], p["lat"]), 2)))
            rep["coastal_pins"] = cp
            # intents with no generated counterpart (centroid not on land and no overlap)
            lost = []
            for it in cfg["intent"]:
                r = np.asarray(it["ring"], float); c = r.mean(0)[None]
                if not any(point_in_ring(c, q)[0] for q in full) and not any(point_in_ring(r[::max(1, len(r) // 20)], q).any() for q in full):
                    lost.append(it["name"])
            rep["lost_intents"] = lost
            # canon land checks: generated land at given distances/bearings from a pin
            lc = []; k = R * math.pi / 180
            for chk in cfg.get("land_checks", []):
                for dname in chk.get("dirs", ["W", "E"]):
                    for dk in chk.get("km", [30, 50]):
                        bx, by = {"N": (0, 1), "S": (0, -1), "E": (1, 0), "W": (-1, 0)}[dname]
                        pt = np.array([[chk["lon"] + bx * dk / (k * math.cos(math.radians(chk["lat"]))), chk["lat"] + by * dk / k]])
                        lc.append(dict(pin=chk["pin"], at=f"{dk}km {dname}", land=bool(any(point_in_ring(pt, q)[0] for q in full))))
            rep["land_checks"] = lc
            return rep
        out["report"] = report
        log(f"wild prep: protect zones {len(zinfo)}, protected cells {int((W > 0).sum())}, excl {'on' if others else 'off'}, {time.time()-t0:.1f}s")
        return out
    return prep

# ============================================================================= v3: coast-following generator
# Signed-distance coast (km, azimuthal-equidistant local grid) + clamped fBm wander, relief on land only.
V3_VERSION = "0.3"
V3_DEFAULTS = dict(wander_km=15.0, rough=1.0, detail_wl_km=60.0, aridity=0.0, feather_km=30.0, arcs_on=False)
V3_LABELS = {"wi": "Western Isles", "haldei": "Haldei", "asahi-arch": "Asahi archipelago", "barsdorf": "Barsdorf", "ne-peninsula": "NE peninsula"}
NOISE_GAIN = 1.0   # set below by calibration constant
CAP_ALPHA = 0.65   # wander is also capped at CAP_ALPHA * local half-width (necks, straits)
INF_SQ = 1e12

def vlog(*a):
    print(*a, file=sys.stderr, flush=True)

# ----------------------------------------------------------------------------- exact EDT (Felzenszwalb-Huttenlocher, vectorised over rows)
def _edt1d_sq(f):
    """f (B,n) float64 (0 at seeds, INF_SQ elsewhere) -> squared distance (cells^2) along axis 1."""
    B, n = f.shape
    if n == 1: return f.copy()
    rows = np.arange(B)
    v = np.zeros((B, n), np.int32); z = np.empty((B, n + 1)); k = np.zeros(B, np.int64)
    z[:, 0] = -np.inf; z[:, 1] = np.inf
    for q in range(1, n):
        fq = f[:, q] + q * q
        act = np.ones(B, bool)
        while True:
            vk = v[rows, k].astype(np.float64)
            s = (fq - (f[rows, vk.astype(np.int64)] + vk * vk)) / (2.0 * (q - vk))
            pop = act & (s <= z[rows, k])
            if not pop.any(): break
            k = k - pop; act = pop
        k = k + 1
        v[rows, k] = q; z[rows, k] = s; z[rows, k + 1] = np.inf
    out = np.empty((B, n)); k = np.zeros(B, np.int64)
    for q in range(n):
        while True:
            m = z[rows, k + 1] < q
            if not m.any(): break
            k = k + m
        vk = v[rows, k].astype(np.float64)
        out[:, q] = (q - vk) ** 2 + f[rows, vk.astype(np.int64)]
    return out

def edt(mask):
    """Exact Euclidean distance (in cells) from every cell to the nearest True cell. float32. 0 on True."""
    f = np.where(mask, 0.0, INF_SQ)
    f = _edt1d_sq(f)
    f = _edt1d_sq(np.ascontiguousarray(f.T)).T
    return np.sqrt(np.minimum(f, INF_SQ)).astype(np.float32)

# ----------------------------------------------------------------------------- connected components (run-length union-find)
def label_mask(m, conn=8):
    """bool (ny,nx) -> (int32 labels 0=bg, n). Run-length graph + min-label hooking with pointer jumping."""
    ny, nx = m.shape
    pad = np.zeros((ny, nx + 2), np.int8); pad[:, 1:-1] = m
    d = np.diff(pad, axis=1)
    ry, rs = np.nonzero(d == 1); _, re = np.nonzero(d == -1)
    nr = len(rs)
    if nr == 0: return np.zeros((ny, nx), np.int32), 0
    K = nx + 2; e = 1 if conn == 8 else 0
    sk = ry * K + rs; ek = ry * K + re
    hi = np.searchsorted(sk, (ry + 1) * K + (re + e), "left")
    lo = np.searchsorted(ek, (ry + 1) * K + (rs - e), "right")
    cnt = np.maximum(hi - lo, 0); cnt[ry >= ny - 1] = 0
    tot = int(cnt.sum())
    if tot:
        ea = np.repeat(np.arange(nr), cnt)
        off = np.arange(tot) - np.repeat(np.cumsum(cnt) - cnt, cnt)
        eb = np.repeat(lo, cnt) + off
        lab = np.arange(nr)
        while True:
            a = lab[ea]; b = lab[eb]
            if not (a != b).any(): break
            mn = np.minimum(a, b)
            np.minimum.at(lab, a, mn); np.minimum.at(lab, b, mn)
            while True:
                l2 = lab[lab]
                if (l2 == lab).all(): break
                lab = l2
    else:
        lab = np.arange(nr)
    _, comp = np.unique(lab, return_inverse=True); comp = comp.astype(np.int32) + 1
    img = np.zeros((ny, nx + 1), np.int32)
    np.add.at(img, (ry, rs), comp); np.add.at(img, (ry, re), -comp)
    out = np.cumsum(img, axis=1)[:, :nx].astype(np.int32)
    return out, int(comp.max())

# ----------------------------------------------------------------------------- projection (azimuthal equidistant about a centre)
class Proj:
    def __init__(self, lon0, lat0, R):
        self.l0, self.p0, self.R = math.radians(lon0), math.radians(lat0), R
        self.lon0, self.lat0 = lon0, lat0
    def fwd(self, lon, lat):
        lo, la = np.radians(np.asarray(lon, float)), np.radians(np.asarray(lat, float))
        dl = lo - self.l0
        cc = np.clip(math.sin(self.p0) * np.sin(la) + math.cos(self.p0) * np.cos(la) * np.cos(dl), -1, 1)
        c = np.arccos(cc); sc = np.sin(c)
        k = np.where(sc < 1e-9, 1.0, c / np.where(sc < 1e-9, 1.0, sc))
        return self.R * k * np.cos(la) * np.sin(dl), self.R * k * (math.cos(self.p0) * np.sin(la) - math.sin(self.p0) * np.cos(la) * np.cos(dl))
    def inv(self, x, y):
        x = np.asarray(x, np.float64); y = np.asarray(y, np.float64)
        rho = np.hypot(x, y); c = rho / self.R; sc, cc = np.sin(c), np.cos(c)
        rs = np.where(rho < 1e-9, 1.0, rho)
        la = np.arcsin(np.clip(cc * math.sin(self.p0) + np.where(rho < 1e-9, 0.0, y * sc * math.cos(self.p0) / rs), -1, 1))
        lo = self.l0 + np.arctan2(x * sc, rho * math.cos(self.p0) * cc - y * math.sin(self.p0) * sc)
        return (np.degrees(lo) + 180) % 360 - 180, np.degrees(la)

def densify_ll(r, step_deg=0.1):
    """Closed lon/lat ring -> denser ring (straight in lon/lat, as the editor draws it)."""
    r = np.asarray(r, float); out = []
    for a, b in zip(r, np.roll(r, -1, axis=0)):
        n = max(1, int(max(abs(b[0] - a[0]), abs(b[1] - a[1])) / step_deg))
        t = np.arange(n)[:, None] / n
        out.append(a + (b - a) * t)
    return np.vstack(out)

def densify_km(P, step):
    """Closed km ring (N,2) -> points at <= step spacing."""
    out = []
    for a, b in zip(P, np.roll(P, -1, axis=0)):
        n = max(1, int(math.ceil(math.hypot(*(b - a)) / step)))
        out.append(a + (b - a) * (np.arange(n)[:, None] / n))
    return np.vstack(out)

def rasterize_cov(rings_px, ny, nx, S=3):
    """rings in grid px coords (x col, y row-up, pixel centres at integers) -> coverage float32 (ny,nx) via SxS supersampling."""
    im = Image.new("L", (nx * S, ny * S), 0); dr = ImageDraw.Draw(im)
    off = (S - 1) / 2.0
    for r in rings_px:
        dr.polygon([(float(x) * S + off, (ny - 1 - float(y)) * S + off) for x, y in r], fill=255)
    a = np.asarray(im, dtype=np.uint8).reshape(ny, S, nx, S)
    cov = a.mean(axis=(1, 3), dtype=np.float32) / 255.0
    return np.flipud(cov).copy()

def sdf_from_mask(inside, res):
    """Signed distance in km: + inland, - at sea; zero crossing midway between cell centres."""
    dl = edt(inside); ds = edt(~inside)
    return np.where(inside, ds - 0.5, -(dl - 0.5)).astype(np.float32) * np.float32(res)

def half_width_ladder(sdf, res, rmax):
    """Local half-width H(x) (km): largest r_k with x inside the morphological opening (r_k) of its own class
    (land cells vs land opening, sea cells vs sea opening). One EDT per rung."""
    rs = [max(res, 1.5)]
    while rs[-1] < rmax: rs.append(rs[-1] * 1.5)
    H = np.full(sdf.shape, rs[0] * 0.5, np.float32)
    for sign in (1, -1):
        s = sdf * sign
        cls = s > 0
        for r in rs:
            S = s > r
            if not S.any(): break
            f = max(1, min(4, int(r / (4 * res))))
            ok = cls[::f, ::f] & (edt(S[::f, ::f]) * res * f <= r + res * f)
            if f > 1: ok = np.repeat(np.repeat(ok, f, 0), f, 1)[:sdf.shape[0], :sdf.shape[1]]
            H = np.where(ok, np.maximum(H, r), H)
    return H, rs

def bilin(F, px, py):
    """Sample F (ny,nx) at fractional grid coords (arrays), edge-clamped."""
    ny, nx = F.shape
    px = np.clip(px, 0, nx - 1.000001); py = np.clip(py, 0, ny - 1.000001)
    ix = np.floor(px).astype(np.int64); iy = np.floor(py).astype(np.int64); tx = px - ix; ty = py - iy
    return (F[iy, ix] * (1 - tx) * (1 - ty) + F[iy, ix + 1] * tx * (1 - ty) + F[iy + 1, ix] * (1 - tx) * ty + F[iy + 1, ix + 1] * tx * ty)

def upsample_lin(F, dc, ny, nx):
    """Coarse grid F (sampled at full-res indices 0, dc, 2dc ...) -> full (ny,nx) float32, separable linear."""
    cy, cx = F.shape
    fx = np.clip(np.arange(nx) / dc, 0, cx - 1.000001); fy = np.clip(np.arange(ny) / dc, 0, cy - 1.000001)
    ix = np.floor(fx).astype(np.int64); iy = np.floor(fy).astype(np.int64); tx = (fx - ix).astype(np.float32); ty = (fy - iy).astype(np.float32)
    tmp = F[:, ix] * (1 - tx)[None] + F[:, ix + 1] * tx[None]
    return (tmp[iy] * (1 - ty)[:, None] + tmp[iy + 1] * ty[:, None]).astype(np.float32)

# ----------------------------------------------------------------------------- noise on the unit sphere, sampled on a decimated km grid
def coarse_axes(nx, ny, dc):
    return np.arange(0, nx, dc), np.arange(0, ny, dc)

def noise_coarse(pj, x0, y0, res, nx, ny, dc, wl_km, octaves, gain, seed, ridged=False, R=4500.0):
    """fBm (v2 helpers) on the coarse grid; finest cell ~ wl_km. Returns (cy,cx) float32 ~[-1,1] (ridged ~[0,1])."""
    ci, cj = coarse_axes(nx, ny, dc)
    X, Y = np.meshgrid(x0 + ci * res, y0 + cj * res)
    lon, lat = pj.inv(X.ravel(), Y.ravel())
    freq = R / (wl_km * 2 ** (octaves - 1))
    out = np.empty(lon.size, np.float32); CH = 200000
    for a in range(0, lon.size, CH):
        P = sph_xyz(lon[a:a + CH], lat[a:a + CH])
        out[a:a + CH] = fbm3(P, freq, octaves, gain, seed, ridged=ridged)
    return out.reshape(X.shape)

WORLD_SEED_DEFAULT = 20261006
COAST_WLS = [480.0, 240.0, 120.0, 60.0, 30.0]   # fixed world stack (finest cell, km); never rescaled per group/bbox
COAST_GAIN = 0.55
REL_BASE_WLS = [700.0, 350.0, 175.0]
REL_RID_WLS = [450.0, 225.0, 112.0, 56.0]

def _rotpow(n):
    out = [np.eye(3)]
    for _ in range(n): out.append(_ROT @ out[-1])
    return out

def octave_values(P, wls, seed, R, ridged=False):
    """P unit vectors (3,N) -> list of per-octave noise arrays (N,), a pure function of the sphere position."""
    rots = _rotpow(len(wls)); out = []
    for o, wl in enumerate(wls):
        q = (rots[o] @ P) * (R / wl)
        n = perlin3(q + 17.31 * o, (seed + 101 * o) & 0xFFFFFFFF)
        if ridged: n = (1.0 - np.abs(n)) ** 2
        out.append(n)
    return out

def detail_emphasis(wl_o, detail_wl):
    """Weight of an octave of cell size wl_o when the zone's detail_wl_km is detail_wl (octaves finer than it fade out)."""
    return smooth01(1.0 + np.log2(wl_o / detail_wl) / 1.5)

def coast_noise(pj, x0, y0, res, nx, ny, wl_map, seed, R, rows=48):
    """World-coherent coast noise (~[-0.5,0.5]), full grid res. wl_map (ny,nx) = detail_wl_km per cell (only re-weights octaves)."""
    out = np.empty((ny, nx), np.float32); norm = sum(COAST_GAIN ** o for o in range(len(COAST_WLS)))
    xs = x0 + np.arange(nx) * res
    for a in range(0, ny, rows):
        b = min(ny, a + rows)
        X, Y = np.meshgrid(xs, y0 + np.arange(a, b) * res)
        lon, lat = pj.inv(X.ravel(), Y.ravel()); P = sph_xyz(lon, lat)
        vals = octave_values(P, COAST_WLS, seed, R); w = wl_map[a:b].ravel().astype(np.float64)
        acc = 0.0
        for o, v in enumerate(vals): acc = acc + COAST_GAIN ** o * detail_emphasis(COAST_WLS[o], w) * v
        out[a:b] = (acc / norm).reshape(b - a, nx)
    return out

def relief_noise(pj, x0, y0, res, nx, ny, seed, R, rows=48):
    """World-coherent (base fBm, ridged fBm) for relief, full grid res."""
    base = np.empty((ny, nx), np.float32); rid = np.empty((ny, nx), np.float32)
    xs = x0 + np.arange(nx) * res
    nb = sum(0.5 ** o for o in range(len(REL_BASE_WLS))); nr = sum(0.5 ** o for o in range(len(REL_RID_WLS)))
    for a in range(0, ny, rows):
        b = min(ny, a + rows)
        X, Y = np.meshgrid(xs, y0 + np.arange(a, b) * res)
        lon, lat = pj.inv(X.ravel(), Y.ravel()); P = sph_xyz(lon, lat)
        vb = octave_values(P, REL_BASE_WLS, seed + 301, R); vr = octave_values(P, REL_RID_WLS, seed + 401, R, ridged=True)
        base[a:b] = (sum(0.5 ** o * v for o, v in enumerate(vb)) / nb).reshape(b - a, nx)
        rid[a:b] = (sum(0.5 ** o * v for o, v in enumerate(vr)) / nr).reshape(b - a, nx)
    return base, rid

def noise_coherence_check(R=4500.0, seed=WORLD_SEED_DEFAULT):
    """Two adjacent overlapping bboxes (shared projection, shifted by whole cells) and a res-2 vs res-4 pair on coincident nodes."""
    pj = Proj(60.0, 40.0, R); wl = lambda ny, nx: np.full((ny, nx), 60.0, np.float32)
    A = coast_noise(pj, -400.0, -300.0, 4.0, 200, 150, wl(150, 200), seed, R)
    B = coast_noise(pj, -400.0 + 120 * 4.0, -300.0 + 20 * 4.0, 4.0, 200, 150, wl(150, 200), seed, R)   # shifted +120 cols, +20 rows
    d1 = float(np.abs(A[20:, 120:] - B[:130, :80]).max())
    C2 = coast_noise(pj, 0.0, 0.0, 2.0, 100, 80, wl(80, 100), seed, R); C4 = coast_noise(pj, 0.0, 0.0, 4.0, 50, 40, wl(40, 50), seed, R)
    d2 = float(np.abs(C2[::2, ::2] - C4).max())
    rb, rr = relief_noise(pj, -400.0, -300.0, 4.0, 200, 150, seed, R); rb2, rr2 = relief_noise(pj, -400.0 + 480.0, -300.0 + 80.0, 4.0, 200, 150, seed, R)
    d3 = float(max(np.abs(rb[20:, 120:] - rb2[:130, :80]).max(), np.abs(rr[20:, 120:] - rr2[:130, :80]).max()))
    return dict(coast_shifted_bbox_max_abs_diff=d1, coast_res2_vs_res4_max_abs_diff=d2, relief_shifted_bbox_max_abs_diff=d3, ok=bool(max(d1, d2, d3) < 1e-5))

def smooth01(t):
    t = np.clip(t, 0, 1); return t * t * (3 - 2 * t)

def parse_params(state, tag):
    gen = state.get("gen", {}) or {}
    P = dict(V3_DEFAULTS); P.update(gen.get("defaults", {}) or {}); P.update((gen.get("groups", {}) or {}).get(tag, {}) or {})
    return P

def group_intent(state, tag):
    return [f for f in state["features"] if f.get("type") == "land" and f.get("kind") == "intent" and tag in (f.get("groups") or [])]

def group_zones(state, tag, defaults):
    zs = [f for f in state["features"] if f.get("type") == "region" and f.get("kind") == "genzone" and tag in (f.get("groups") or [])]
    zs.sort(key=lambda f: f.get("priority", 0))
    return zs

def zone_weight(pj, zone, x0, y0, res, nx, ny, feather):
    """Full inside the polygon, ramping to 0 over `feather` km outside. Returns (window slices, weight array) or None."""
    ring = np.asarray(zone["geometry"], float)
    zx, zy = pj.fwd(*densify_ll(ring, 0.1).T)
    pad = feather + 3 * res
    i0 = max(0, int((zx.min() - pad - x0) / res)); i1 = min(nx, int((zx.max() + pad - x0) / res) + 2)
    j0 = max(0, int((zy.min() - pad - y0) / res)); j1 = min(ny, int((zy.max() + pad - y0) / res) + 2)
    if i1 <= i0 or j1 <= j0: return None
    px = (zx - x0) / res - i0; py = (zy - y0) / res - j0
    cov = rasterize_cov([np.stack([px, py], 1)], j1 - j0, i1 - i0)
    ins = cov >= 0.5
    if not ins.any(): return None
    sd = sdf_from_mask(np.pad(ins, 1, constant_values=False), res)[1:-1, 1:-1]
    w = np.clip(1.0 + sd / max(feather, res), 0, 1) if feather > res else (sd > 0).astype(np.float32)
    return (slice(j0, j1), slice(i0, i1)), smooth01(w).astype(np.float32)

# ----------------------------------------------------------------------------- arcs
class ArcBumps:
    """Opt-in volcanic chains. Bump (km): A = (amp*bead*gauss(d/width) - 0.5) * width * 0.8; C = max(C, A * arcs_on)."""
    def __init__(self, arcs, pj, seed):
        self.arcs = []
        rng = np.random.default_rng(seed + 991)
        for a in arcs:
            pts = np.asarray(a["pts"], float)
            x, y = pj.fwd(pts[:, 0], pts[:, 1]); K = np.stack([x, y], 1)
            if len(K) >= 3 and a.get("kind", "spline") != "gc":
                P = np.vstack([2 * K[0] - K[1], K, 2 * K[-1] - K[-2]]); dense = []
                for i in range(1, len(P) - 2):
                    p0, p1, p2, p3 = P[i - 1], P[i], P[i + 1], P[i + 2]
                    for t in np.linspace(0, 1, 40, endpoint=False):
                        dense.append(0.5 * ((2 * p1) + (-p0 + p2) * t + (2 * p0 - 5 * p1 + 4 * p2 - p3) * t * t + (-p0 + 3 * p1 - 3 * p2 + p3) * t ** 3))
                dense.append(P[-2]); D = np.array(dense)
            else:
                D = K
            D = densify_open(D, 3.0)
            s = np.concatenate([[0], np.cumsum(np.hypot(*np.diff(D, axis=0).T))]); L = s[-1]
            nb = max(3, int(L / a.get("bead_km", 110.0)) + 2)
            sg = a.get("bead_sigma", 1.0)
            vals = rng.standard_normal(nb) * sg
            amp_s = np.minimum(np.exp(vals - 0.3) if sg > 0 else np.ones(nb), a.get("bead_max", 3.5))
            fine = np.interp(s, np.linspace(0, L, nb), amp_s)
            fine = np.convolve(np.pad(fine, 3, mode="edge"), np.ones(7) / 7, mode="valid")
            taper = np.minimum(1, np.minimum(s, L - s) / max(1.0, a.get("taper_km", 150.0)))
            self.arcs.append(dict(D=D, amp=fine * taper, w=a.get("amp", 1.0), sig=a.get("width_km", 45.0)))
    def eval(self, x0, y0, res, nx, ny, arcs_on):
        """Returns (A, u): A the bump in km (-inf where none), u the unit strength (for relief)."""
        A = np.full((ny, nx), -np.inf, np.float32); U = np.zeros((ny, nx), np.float32)
        for a in self.arcs:
            D = a["D"]; sig = a["sig"]; pad = 4.5 * sig
            i0 = max(0, int((D[:, 0].min() - pad - x0) / res)); i1 = min(nx, int((D[:, 0].max() + pad - x0) / res) + 2)
            j0 = max(0, int((D[:, 1].min() - pad - y0) / res)); j1 = min(ny, int((D[:, 1].max() + pad - y0) / res) + 2)
            if i1 <= i0 or j1 <= j0: continue
            X, Y = np.meshgrid(x0 + np.arange(i0, i1) * res, y0 + np.arange(j0, j1) * res)
            q = np.stack([X.ravel(), Y.ravel()], 1).astype(np.float32); Df = D.astype(np.float32)
            d = np.empty(len(q), np.float32); jn = np.empty(len(q), np.int64)
            step = max(1, 4_000_000 // len(Df))
            for b in range(0, len(q), step):
                dd = ((q[b:b + step, None, 0] - Df[None, :, 0]) ** 2 + (q[b:b + step, None, 1] - Df[None, :, 1]) ** 2)
                jj = dd.argmin(1); jn[b:b + step] = jj; d[b:b + step] = np.sqrt(dd[np.arange(len(jj)), jj])
            u = (a["w"] * a["amp"][jn] * np.exp(-0.5 * (d / sig) ** 2)).reshape(X.shape) * arcs_on[j0:j1, i0:i1]
            Aw = ((u - 0.5) * sig * 0.8).astype(np.float32)
            A[j0:j1, i0:i1] = np.maximum(A[j0:j1, i0:i1], Aw); U[j0:j1, i0:i1] = np.maximum(U[j0:j1, i0:i1], u.astype(np.float32))
        return A, U

def densify_open(D, step):
    out = []
    for a, b in zip(D[:-1], D[1:]):
        n = max(1, int(math.ceil(math.hypot(*(b - a)) / step)))
        out.append(a + (b - a) * (np.arange(n)[:, None] / n))
    out.append(D[-1:]); return np.vstack(out)

NOISE_GAIN = 3.0
POND_MIN_KM2 = 400.0        # enclosed generated water overlapping intent sea: smaller = filled, >= this = a cut-off gulf (issue)
SPLIT_MIN_KM2 = 2000.0      # intent land components at least this big must map to exactly one generated component

# ----------------------------------------------------------------------------- pins
def assign_pins(state, tag, pj, groups_rings, reach_km):
    """Each pin goes to the group whose intent land contains it, else the nearest. Returns dict(land, water, coastal) lists of
    dict(name, lon, lat, x, y, side)."""
    R = pj.R; mine = {"land": [], "water": [], "coastal": []}; notes = []
    own = groups_rings
    for f in state["features"]:
        if f.get("type") != "pin": continue
        grp = f.get("groups") or []
        if f.get("kind") == "unplaced" or "unplaced" in grp: continue
        g = f["geometry"]; p = np.array([[g["lon"], g["lat"]]])
        owner = None
        for t, rings in own.items():
            if any(point_in_ring(p, r)[0] for r in rings): owner = (t, True); break
        if owner is None:
            xyz = sph_xyz(p[:, 0], p[:, 1])[:, 0]; best = None
            for t, rings in own.items():
                for r in rings:
                    d = densify_ll(r, 0.25); dd = np.linalg.norm(sph_xyz(d[:, 0], d[:, 1]).T - xyz, axis=1).min() * R
                    if best is None or dd < best[0]: best = (dd, t)
            owner = (best[1], False, best[0])
        if owner[0] != tag: continue
        x, y = pj.fwd(g["lon"], g["lat"]); inside = owner[1]
        rec = dict(name=f.get("name", f.get("id", "?")), lon=g["lon"], lat=g["lat"], x=float(x), y=float(y), inside=inside)
        if f.get("kind") == "sea":
            if inside: notes.append(f"sea pin '{rec['name']}' lies inside intent land: forced to water locally")
            mine["water"].append(rec)
        elif f.get("coastal"):
            rec["side"] = 1 if inside else -1; mine["coastal"].append(rec)
            if inside: mine["land"].append(rec)
        elif inside:
            mine["land"].append(rec)
        else:
            if "sea" in grp: notes.append(f"pin '{rec['name']}' (kind {f.get('kind')}, tagged sea) is outside intent land, not enforced as land or water")
    return mine, notes

def force_pins(C, pins, x0, y0, res, sig_km=8.0, log=None):
    """Local Gaussian dents/bumps so land pins sit on land, sea pins in water, coastal pins at C ~ +-1 km. In place."""
    margin = max(1.5, 0.8 * res)
    ny, nx = C.shape; sig = sig_km; rad = 4 * sig
    items = [(p, "coast") if p in pins["coastal"] else (p, "land") for p in pins["land"] + [q for q in pins["coastal"] if q not in pins["land"]]] + [(p, "water") for p in pins["water"]]
    for it in range(14):
        changed = 0
        for p, kind in items:
            px, py = (p["x"] - x0) / res, (p["y"] - y0) / res
            cv = float(bilin(C, np.array([px]), np.array([py]))[0])
            if kind == "coast":
                t = 1.8 * p["side"]
                if abs(cv - t) <= 0.9: continue
            elif kind == "land":
                if cv >= margin: continue
                t = margin + 0.5
            else:
                if cv <= -margin: continue
                t = -margin - 0.5
            i0 = max(0, int(px - rad / res)); i1 = min(nx, int(px + rad / res) + 2); j0 = max(0, int(py - rad / res)); j1 = min(ny, int(py + rad / res) + 2)
            X, Y = np.meshgrid((np.arange(i0, i1) - px) * res, (np.arange(j0, j1) - py) * res)
            C[j0:j1, i0:i1] += (t - cv) * np.exp(-0.5 * (X * X + Y * Y) / sig ** 2).astype(np.float32); changed += 1
        if not changed: break
    return C

# ----------------------------------------------------------------------------- marching squares, int32 / low memory
def marching_squares32(E, sl, L):
    """As v2 marching_squares (land on the left, CCW outer) with int32 ids and sparse links."""
    ny, nx = E.shape
    E = np.where(L & (E <= sl), np.float32(sl + 1e-4), E); E = np.where(~L & (E > sl), np.float32(sl - 1e-4), E)
    C = [L[:-1, :-1], L[:-1, 1:], L[1:, 1:], L[1:, :-1]]
    yy, xx = np.mgrid[0:ny - 1, 0:nx - 1]; yy = yy.astype(np.int32); xx = xx.astype(np.int32)
    NH = ny * nx
    eid = [yy * nx + xx, NH + yy * nx + (xx + 1), (yy + 1) * nx + xx, NH + yy * nx + xx]
    ex = [C[k] & ~C[(k + 1) % 4] for k in range(4)]
    en = [~C[k] & C[(k + 1) % 4] for k in range(4)]
    saddle = (C[0] == C[2]) & (C[1] == C[3]) & (C[0] != C[1])
    center = (E[:-1, :-1] + E[:-1, 1:] + E[1:, 1:] + E[1:, :-1]) * 0.25 > sl
    use_next = saddle & center
    src = []; dst = []
    for k in range(4):
        def partner(direction):
            pj = np.full(C[0].shape, -1, np.int8)
            for step in (3, 2, 1):
                j = (k + direction * step) % 4
                pj = np.where(en[j], np.int8(j), pj)
            return pj
        pj = np.where(use_next, partner(+1), partner(-1))
        m = ex[k] & (pj >= 0)
        if not m.any(): continue
        to = np.choose(np.clip(pj, 0, 3).astype(np.int64), [eid[0], eid[1], eid[2], eid[3]])
        src.append(eid[k][m]); dst.append(to[m])
    if not src: return []
    src = np.concatenate(src); dst = np.concatenate(dst)
    ids = np.unique(np.concatenate([src, dst]))
    posx = {};
    hm = ids < NH
    h = ids[hm]; y, x = h // nx, h % nx
    t = (sl - E[y, x]) / (E[y, x + 1] - E[y, x]); px_h = x + t; py_h = y.astype(np.float32)
    v = ids[~hm] - NH; y, x = v // nx, v % nx
    t = (sl - E[y, x]) / (E[y + 1, x] - E[y, x]); px_v = x.astype(np.float32); py_v = y + t
    pos = np.empty((len(ids), 2), np.float32)
    pos[np.nonzero(hm)[0]] = np.stack([px_h, py_h], 1); pos[np.nonzero(~hm)[0]] = np.stack([px_v, py_v], 1)
    sidx = np.searchsorted(ids, src); didx = np.searchsorted(ids, dst)
    nxt = np.full(len(ids), -1, np.int64); nxt[sidx] = didx
    nl = nxt.tolist(); seen = np.zeros(len(ids), bool); rings = []
    for s in sidx.tolist():
        if seen[s]: continue
        cur = s; seq = []
        while cur >= 0 and not seen[cur]:
            seen[cur] = True; seq.append(cur); cur = nl[cur]
        if cur == s and len(seq) >= 4: rings.append(pos[np.array(seq)].astype(np.float64))
    return rings

# ----------------------------------------------------------------------------- topology + cleanup on the RAW mask
def topology_check(C, sdf, il, ni, res, pj, x0, y0, land_pin_cells, min_island_km2, lake=None):
    """C>0 raw mask. Returns (issues, protect, L): issues = [{kind,text,lon,lat}], protect = [(ys,xs)] cells to calm down, L = cleaned land mask."""
    ny, nx = C.shape; cell = res * res; issues = []; protect = []
    def where(ys, xs):
        lo, la = pj.inv(x0 + float(np.mean(xs)) * res, y0 + float(np.mean(ys)) * res); return round(float(lo), 2), round(float(la), 2)
    Lraw = C > 0; Lraw[:3, :] = False; Lraw[-3:, :] = False; Lraw[:, :3] = False; Lraw[:, -3:] = False
    wl, nw = label_mask(~Lraw, 4)
    sea = wl[0, 0]
    if nw > 1:
        cnt = np.bincount(wl.ravel(), minlength=nw + 1)
        isea = (sdf < 0)
        ov = np.bincount(wl[isea].ravel(), minlength=nw + 1)
        lov = np.bincount(wl[lake].ravel(), minlength=nw + 1) if lake is not None else np.zeros(nw + 1, np.int64)
        lakes_kept = [k for k in range(1, nw + 1) if k != sea and cnt[k] and lov[k] >= 0.5 * cnt[k]]
        for k in range(1, nw + 1):
            if k == sea or cnt[k] == 0 or k in lakes_kept: continue
            area = cnt[k] * cell
            if area >= POND_MIN_KM2 and ov[k] >= 0.5 * cnt[k]:
                ys, xs = np.nonzero(wl == k); lo, la = where(ys, xs)
                issues.append(dict(kind="cutoff", text=f"water body of intent sea cut off from the open sea ({area:,.0f} km2) near lon {lo}, lat {la}", lon=lo, lat=la)); protect.append((ys, xs))
        keepw = (wl == sea) | (wl == 0) | np.isin(wl, lakes_kept)
        L = Lraw | ~keepw
    else:
        L = Lraw.copy()
    ll, nl = label_mask(L, 8)
    cnt = np.bincount(ll.ravel(), minlength=nl + 1)
    keep = cnt * cell >= min_island_km2; keep[0] = False
    icnt0 = np.bincount(il.ravel(), minlength=ni + 1)
    mo = (il > 0) & (ll > 0)
    if mo.any():
        pr = np.unique(ll[mo].astype(np.int64) * (ni + 1) + il[mo]); gg = pr // (ni + 1); ii_ = pr % (ni + 1)
        sm = gg[(icnt0[ii_] * cell < SPLIT_MIN_KM2) & (cnt[gg] >= 2)]; keep[sm] = True   # small intent islands are kept whole
    if land_pin_cells is not None and len(land_pin_cells[0]):
        pl = np.unique(ll[land_pin_cells]); keep[pl[pl > 0]] = True
    L = keep[ll]; ll = np.where(L, ll, 0)
    cnt = np.where(keep, cnt, 0)
    m = (il > 0) & (ll > 0)
    pairs = np.unique(il[m].astype(np.int64) * (nl + 1) + ll[m]) if m.any() else np.array([], np.int64)
    pi = pairs // (nl + 1); pg = pairs % (nl + 1)
    icnt = np.bincount(il.ravel(), minlength=ni + 1)
    for i in range(1, ni + 1):
        if icnt[i] * cell < SPLIT_MIN_KM2: continue
        gs = pg[pi == i]
        ys, xs = np.nonzero(il == i); lo, la = where(ys, xs)
        if len(gs) == 0:
            issues.append(dict(kind="vanished", text=f"intent land component ({icnt[i]*cell:,.0f} km2) near lon {lo}, lat {la} has no generated land", lon=lo, lat=la)); protect.append((ys, xs))
        elif len(gs) > 1:
            big = gs[np.argmax(cnt[gs])]; fr = [g for g in gs if g != big]
            fy, fx = np.nonzero(np.isin(ll, fr)); lo, la = where(fy, fx)
            fd = []
            for g_ in fr[:5]:
                gy, gx = np.nonzero(ll == g_); flo, fla = where(gy, gx); fd.append(f"{cnt[g_]*cell:,.0f} km2 at lon {flo}, lat {fla}")
            issues.append(dict(kind="split", text=f"intent land component ({icnt[i]*cell:,.0f} km2) split into {len(gs)} generated components; fragments: " + "; ".join(fd), lon=lo, lat=la)); protect.append((fy, fx))
    for g in np.unique(pg):
        ins = pi[pg == g]
        ins = [i for i in ins if icnt[i] >= 3]
        if len(ins) > 1:
            by, bx = np.nonzero((ll == g) & (sdf < 0));
            if len(by) == 0: by, bx = np.nonzero(ll == g)
            lo, la = where(by, bx)
            issues.append(dict(kind="bridge", text=f"one generated component bridges {len(ins)} separate intent components near lon {lo}, lat {la}", lon=lo, lat=la)); protect.append((by, bx))
    return issues, protect, L

def calm_around(damp, cells, x0, y0, res, radius_km, factor=0.5):
    ys, xs = cells; ny, nx = damp.shape
    if len(ys) == 0: return
    pad = int(radius_km / res) + 2
    j0 = max(0, int(ys.min()) - pad); j1 = min(ny, int(ys.max()) + pad + 1); i0 = max(0, int(xs.min()) - pad); i1 = min(nx, int(xs.max()) + pad + 1)
    m = np.zeros((j1 - j0, i1 - i0), bool); m[ys - j0, xs - i0] = True
    d = edt(m) * res
    w = np.clip(1 - d / radius_km, 0, 1)
    damp[j0:j1, i0:i1] *= (1 - (1 - factor) * smooth01(w * 2)).astype(np.float32)

# ----------------------------------------------------------------------------- metrics
def ring_dist_to_polylines(P, Q):
    """Min distance from points P (n,2) to the closed polyline set Q = list of (m,2) rings. Brute force, chunked, via segment geometry."""
    out = np.full(len(P), np.inf)
    for r in Q:
        A = r; B = np.roll(r, -1, axis=0); AB = B - A; L2 = np.maximum((AB ** 2).sum(1), 1e-12)
        step = max(1, 3_000_000 // len(A))
        for s in range(0, len(P), step):
            p = P[s:s + step, None, :]
            t = np.clip(((p - A[None]) * AB[None]).sum(2) / L2[None], 0, 1)
            d = np.sqrt(((p - A[None] - t[..., None] * AB[None]) ** 2).sum(2)).min(1)
            out[s:s + step] = np.minimum(out[s:s + step], d)
    return out

# ----------------------------------------------------------------------------- biome (lon/lat raster, same file format as v2)
def biome_v3(pj, x0, y0, res, rel_m, Lmask, arid, bz, bbox_ll, bc, R, log, spec=None):
    """rel_m: land elevation (m) on the km grid (0 at sea). Lmask: bool land. arid: aridity map (km grid). bz: uint8 biome override."""
    bres = bc.get("res_deg", 0.1)
    bx0 = math.floor(bbox_ll[0] / bres) * bres - 2 * bres; by0 = math.floor(bbox_ll[1] / bres) * bres - 2 * bres
    nx2 = int(math.ceil((bbox_ll[2] - bx0) / bres)) + 2; ny2 = int(math.ceil((bbox_ll[3] - by0) / bres)) + 2
    lon_c = bx0 + (np.arange(nx2) + 0.5) * bres; lat_c = by0 + (np.arange(ny2) + 0.5) * bres
    LO, LA = np.meshgrid(lon_c, lat_c)
    ny, nx = Lmask.shape
    Lf = Lmask.astype(np.float32); lf = np.zeros(LO.shape, np.float32)
    offs = (-bres / 3, 0.0, bres / 3)
    for dx in offs:
        for dy in offs:
            X, Y = pj.fwd(LO + dx, LA + dy)
            lf += bilin(Lf, (X - x0) / res, (Y - y0) / res) * (1 / 9)
    X, Y = pj.fwd(LO, LA); PX, PY = (X - x0) / res, (Y - y0) / res
    inb = (PX >= 0) & (PX <= nx - 1) & (PY >= 0) & (PY <= ny - 1)
    lf = np.where(inb, lf, 0).astype(np.float32)
    ev = np.where(inb, bilin(rel_m, PX, PY), 0.0)
    ar = np.where(inb, bilin(arid, PX, PY), 0.0)
    zb = np.where(inb, bz[np.clip(np.round(PY).astype(int), 0, ny - 1), np.clip(np.round(PX).astype(int), 0, nx - 1)], 0)
    lat2 = lat_c[:, None] * np.ones((1, nx2))
    cls, info = climate_classify(lat_c, lon_c, ev, lf, ar, bc, R, bres, spec, override=zb.astype(np.uint8))   # v3: aridity > 0 = drier
    lm = lf > 0.02
    d = lm.copy(); d[1:] |= lm[:-1]; d[:-1] |= lm[1:]; d[:, 1:] |= lm[:, :-1]; d[:, :-1] |= lm[:, 1:]
    out = np.where(d, cls, 0).astype(np.uint8)
    wt = lf * np.cos(np.radians(lat2)); wsum = max(float(wt.sum()), 1e-9)
    pct = {BIOMES[i - 1][0]: round(100 * float(wt[(cls == i)].sum() / wsum), 1) for i in range(1, len(BIOMES) + 1)}
    biome_v3.last_report = climate_report(lat_c, lon_c, cls, info, lf, spec, bres)
    flat = np.flipud(out).ravel(); brk = np.nonzero(np.diff(flat))[0] + 1; st = np.concatenate([[0], brk]); ln = np.diff(np.concatenate([st, [len(flat)]]))
    rle = np.stack([flat[st], ln], axis=1).ravel().astype(int).tolist()
    bbox = [round(bx0, 4), round(by0, 4), round(bx0 + nx2 * bres, 4), round(by0 + ny2 * bres, 4)]
    blob = dict(bbox=bbox, res_deg=round(bres, 4), w=nx2, h=ny2, order="rows north to south, [class, run] pairs",
                classes=[dict(id=i + 1, name=n, color=c) for i, (n, c) in enumerate(BIOMES)], rle=rle, aridity=round(float(np.mean(ar[lm])) if lm.any() else 0.0, 3),
                zones=[], pct_land=pct)
    return blob, out, (bx0, by0, bres)

# ----------------------------------------------------------------------------- preview
def render_preview_v3(path, tag, pj, x0, y0, res, C, rel, Lmask, intent_km, gen_km, pins, zones_km, insets, biome_cls, bgeo, title, stats_txt):
    ny, nx = C.shape
    pal = np.array([[96, 96, 96]] + [[int(c[i:i + 2], 16) for i in (1, 3, 5)] for _, c in BIOMES], np.float32) / 255
    # biome class at each km cell (lon/lat raster lookup), chunked
    cls = np.zeros((ny, nx), np.uint8); bx0, by0, bres = bgeo; ch2, cw2 = biome_cls.shape
    for a in range(0, ny, 256):
        X, Y = np.meshgrid(x0 + np.arange(nx) * res, y0 + np.arange(a, min(ny, a + 256)) * res)
        lo, la = pj.inv(X, Y)
        ii = np.clip(((lo - bx0) / bres).astype(int), 0, cw2 - 1); jj = np.clip(((la - by0) / bres).astype(int), 0, ch2 - 1)
        cls[a:a + 256] = biome_cls[jj, ii]
    cls = np.where(cls == 0, 4, cls)
    gy, gx = np.gradient(rel * (1.5 / res) * 10.0)
    shade = np.clip(0.55 + 0.9 * (gx - gy) / 1.4, 0.15, 1.2).astype(np.float32)
    land_c = pal[cls] * (0.55 + 0.6 * shade[..., None]) * (0.9 + 0.1 * np.clip(rel, 0, 1.2)[..., None])
    depth = np.clip(-C / 250.0, 0, 1)[..., None]
    sea_c = np.array([0.55, 0.78, 0.88], np.float32) * (1 - depth) + np.array([0.08, 0.20, 0.36], np.float32) * depth
    img = np.where(Lmask[..., None], land_c, sea_c)
    full = Image.fromarray((np.clip(np.flipud(img), 0, 1) * 255).astype(np.uint8), "RGB"); del img, land_c, sea_c, cls
    try: font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 14)
    except Exception: font = ImageFont.load_default()
    def tx(P, ox, oy, s):  # km points -> image px (row 0 at top)
        return [((x - x0) / res * s - ox, (ny - 1 - (y - y0) / res) * s - oy) for x, y in P]
    def draw(im, ox, oy, s, lw, rad, label_pins):
        d = ImageDraw.Draw(im, "RGBA")
        for z in zones_km:
            p = tx(z, ox, oy, s); d.line(p + [p[0]], fill=(0, 220, 255, 220), width=lw)
        for r in intent_km:
            p = tx(r, ox, oy, s); d.line(p + [p[0]], fill=(255, 255, 255, 235), width=lw)
        for r in gen_km:
            p = tx(r, ox, oy, s); d.line(p + [p[0]], fill=(15, 15, 15, 255), width=lw)
        for kind, lst in pins.items():
            for p in lst:
                x, y = tx([(p["x"], p["y"])], ox, oy, s)[0]
                col = {"land": (230, 40, 40, 255), "water": (40, 120, 255, 255), "coastal": (255, 220, 0, 255)}[kind]
                if p.get("ok") is False: col = (255, 0, 255, 255)
                if kind == "water": d.rectangle([x - rad, y - rad, x + rad, y + rad], fill=col, outline=(255, 255, 255, 255))
                else: d.ellipse([x - rad, y - rad, x + rad, y + rad], fill=col, outline=(255, 255, 255, 255))
                if label_pins and kind != "land": d.text((x + rad + 2, y - rad), p["name"][:22], fill=(255, 255, 255, 255), font=font, stroke_width=1, stroke_fill=(0, 0, 0, 255))
        return d
    W = 1800; f = W / nx if nx > W else 1.0
    main = full.resize((int(nx * f), int(ny * f)), Image.BOX if f < 1 else Image.NEAREST)
    draw(main, 0, 0, f, 1, 3, True)
    tiles = []
    for name, cx, cy, half in insets:
        s = max(1, int(round(560 / (2 * half / res))))
        i0 = int((cx - half - x0) / res); j0 = int((ny - 1) - (cy + half - y0) / res); n = int(2 * half / res)
        i0 = max(0, min(nx - n, i0)); j0 = max(0, min(ny - n, j0))
        t = full.crop((i0, j0, i0 + n, j0 + n)).resize((n * s, n * s), Image.BILINEAR)
        draw(t, i0 * s, j0 * s, s, 2, 5, True)
        ImageDraw.Draw(t).text((6, 4), f"{name} ({2*half:.0f} km box)", fill=(255, 255, 255), font=font, stroke_width=1, stroke_fill=(0, 0, 0))
        tiles.append(t)
    strip_h = max([t.height for t in tiles], default=0)
    outim = Image.new("RGB", (main.width, main.height + strip_h + 28), (24, 24, 24))
    outim.paste(main, (0, 0))
    x = 0
    for t in tiles:
        t2 = t.resize((min(t.width, main.width // max(1, len(tiles))), min(t.width, main.width // max(1, len(tiles))))) if t.width * len(tiles) > main.width else t
        outim.paste(t2, (x, main.height)); x += t2.width + 4
    dd = ImageDraw.Draw(outim)
    dd.text((6, outim.height - 24), title + "   " + stats_txt, fill=(255, 255, 255), font=font)
    outim.save(path)

# ----------------------------------------------------------------------------- v3 driver
def run_group(tag, state_path, out_path, png_path, seed=None, res_override=None, extra=None):
    t0 = time.time(); extra = extra or {}
    state = json.load(open(state_path)); R = float(state["planet"]["radius_km"])
    P = parse_params(state, tag); seed_cli = seed; seed = int((state.get('gen') or {}).get('world_seed', WORLD_SEED_DEFAULT))
    feats = group_intent(state, tag)
    if not feats: raise SystemExit(f"no intent land carries group '{tag}'")
    rings_ll = [np.asarray(f["geometry"], float) for f in feats]
    xyz = sph_xyz(*np.vstack([densify_ll(r, 0.5) for r in rings_ll]).T).mean(axis=1); xyz /= np.linalg.norm(xyz)
    lon0 = math.degrees(math.atan2(xyz[1], xyz[0])); lat0 = math.degrees(math.asin(xyz[2]))
    pj = Proj(lon0, lat0, R)
    zones = group_zones(state, tag, P)
    zparams = [z.get("params", {}) or {} for z in zones]
    wmax = max([P["wander_km"]] + [zp.get("wander_km", 0.0) for zp in zparams])
    km_rings = [np.stack(pj.fwd(*densify_ll(r, 0.1).T), 1) for r in rings_ll]
    allp = np.vstack(km_rings); margin = 3 * wmax + 150
    x0, x1 = allp[:, 0].min() - margin, allp[:, 0].max() + margin; y0, y1 = allp[:, 1].min() - margin, allp[:, 1].max() + margin
    res = res_override or max(2.0, math.sqrt((x1 - x0) * (y1 - y0) / 6.0e6))
    nx = int((x1 - x0) / res) + 1; ny = int((y1 - y0) / res) + 1
    vlog(f"[{tag}] centre {lon0:.2f}E {lat0:.2f}N, grid {nx}x{ny} ({nx*ny/1e6:.2f}M) @ {res:.2f} km, seed {seed}, wander {P['wander_km']} (max {wmax}), zones {len(zones)}")
    # intent raster + SDF
    cov = rasterize_cov([np.stack([(r[:, 0] - x0) / res, (r[:, 1] - y0) / res], 1) for r in km_rings], ny, nx)
    inside = cov >= 0.5; del cov
    lake_feats = [f for f in state["features"] if f.get("type") == "region" and f.get("kind") == "lake" and tag in (f.get("groups") or [])]
    lakeM = None
    if lake_feats:   # drawn/baked lakes are intent water inside land (land rings stay single-ring)
        lk_km = [np.stack(pj.fwd(*densify_ll(np.asarray(f["geometry"], float), 0.1).T), 1) for f in lake_feats]
        lakeM = rasterize_cov([np.stack([(r[:, 0] - x0) / res, (r[:, 1] - y0) / res], 1) for r in lk_km], ny, nx) >= 0.5
        inside &= ~lakeM
        vlog(f"  lakes: {len(lake_feats)} carved from the intent")
    sdf = sdf_from_mask(inside, res)
    il, ni = label_mask(inside, 8)
    vlog(f"  sdf done {time.time()-t0:.1f}s, intent components {ni}")
    # local half-width cap
    H, rungs = half_width_ladder(sdf, res, wmax / CAP_ALPHA)
    cap = (CAP_ALPHA * H).astype(np.float32); del H
    vlog(f"  half-width ladder {len(rungs)} rungs, {time.time()-t0:.1f}s")
    # per-cell params (group defaults -> genzones in priority order, feathered)
    wander = np.full((ny, nx), P["wander_km"], np.float32); rough = np.full((ny, nx), P["rough"], np.float32)
    wl = np.full((ny, nx), P["detail_wl_km"], np.float32); arid = np.full((ny, nx), P["aridity"], np.float32)
    aon = np.full((ny, nx), 1.0 if P.get("arcs_on") else 0.0, np.float32); bz = np.zeros((ny, nx), np.uint8)
    zones_km = []
    for z, zp in zip(zones, zparams):
        zr = np.asarray(z["geometry"], float); zones_km.append(np.stack(pj.fwd(*densify_ll(zr, 0.1).T), 1))
        zw = zone_weight(pj, z, x0, y0, res, nx, ny, float(z.get("feather_km", P["feather_km"])))
        if zw is None: continue
        sl_, w = zw
        for key, arr in (("wander_km", wander), ("rough", rough), ("detail_wl_km", wl), ("aridity", arid), ("arcs_on", aon)):
            if key in zp:
                v = (1.0 if zp[key] else 0.0) if key == "arcs_on" else float(zp[key])
                arr[sl_] = arr[sl_] * (1 - w) + v * w
        if zp.get("biome"):
            bz[sl_][w > 0.5] = biome_id(zp["biome"])
    # noise: ONE world seed + fixed octave stack (pure function of lon/lat); params only scale / re-weight it
    wseed = int((state.get("gen") or {}).get("world_seed", WORLD_SEED_DEFAULT))
    if "seed" in P or any("seed" in (v or {}) for v in ((state.get("gen") or {}).get("groups") or {}).values()) or seed_cli is not None:
        vlog("  note: per-group/CLI seeds are ignored (world-coherent noise, world_seed %d)" % wseed)
    nz = coast_noise(pj, x0, y0, res, nx, ny, wl, wseed, R)
    nclamp = np.clip(rough * NOISE_GAIN * nz, -1, 1).astype(np.float32); del nz, rough
    vlog(f"  coast noise done {time.time()-t0:.1f}s")
    # arcs (opt-in)
    arcs_src = [dict(pts=f["geometry"], **(f.get("params", {}) or {})) for f in state["features"] if f.get("type") == "line" and f.get("kind") == "arc" and tag in (f.get("groups") or [])]
    if not arcs_src: arcs_src = (extra.get("arcs") or [])
    arcA = arcU = None
    if arcs_src and aon.max() > 0:
        arcA, arcU = ArcBumps(arcs_src, pj, seed).eval(x0, y0, res, nx, ny, aon)
        vlog(f"  arcs: {len(arcs_src)} on, bump unit km = (amp*bead*gauss - 0.5)*width*0.8")
    # pins
    groups_rings = {}
    for f in state["features"]:
        if f.get("type") == "land" and f.get("kind") == "intent":
            for g_ in f.get("groups") or []: groups_rings.setdefault(g_, []).append(np.asarray(f["geometry"], float))
    pins, pin_notes = assign_pins(state, tag, pj, groups_rings, wmax)
    for p in pins["land"] + pins["water"] + pins["coastal"]:
        pass
    for n_ in pin_notes: vlog("  note:", n_)
    vlog(f"  pins: land {len(pins['land'])}, water {len(pins['water'])}, coastal {len(pins['coastal'])}")
    min_isl = float(P.get("min_island_km2", extra.get("min_island_km2", 150.0)))
    damp = np.ones((ny, nx), np.float32)
    ok_topo = False; issues = []
    for rep in range(8):
        weff = np.minimum(wander * damp, cap)
        C = (sdf + weff * nclamp).astype(np.float32)
        if arcA is not None: C = np.where(np.isfinite(arcA), np.maximum(C, arcA), C)
        force_pins(C, pins, x0, y0, res)
        lpc = None
        if pins["land"]:
            ii = np.clip(np.round((np.array([p["x"] for p in pins["land"]]) - x0) / res).astype(int), 0, nx - 1); jj = np.clip(np.round((np.array([p["y"] for p in pins["land"]]) - y0) / res).astype(int), 0, ny - 1)
            lpc = (jj, ii)
        issues, protect, L = topology_check(C, sdf, il, ni, res, pj, x0, y0, lpc, min_isl, lake=lakeM)
        near = np.abs(sdf) < 2 * wander
        capped_share = float(((wander * damp > cap) & near).sum() / max(1, near.sum()))
        vlog(f"  round {rep}: issues {len(issues)}, capped share near coast {100*capped_share:.1f}%, {time.time()-t0:.1f}s")
        if not issues: ok_topo = True; break
        for cells in protect: calm_around(damp, cells, x0, y0, res, 90.0)
    stats = dict(group=tag)
    intent_area = float(sum(ring_area_km2(r, R) for r in rings_ll)); stats["intent_area_km2"] = round(intent_area)
    if not ok_topo:
        stats.update(area_km2=0, coast_dev_p95_km=None, coast_dev_max_km=None, topology_ok=False, topology_issues=[i["text"] for i in issues],
                     pins_land_ok=None, pins_water_ok=None, pins_coastal_ok=None)
        return stats, t0
    # rings from the cleaned mask
    raw = marching_squares32(C, 0.0, L)
    kept = []
    for r in raw:
        sa = 0.5 * np.sum(r[:, 0] * np.roll(r[:, 1], -1) - np.roll(r[:, 0], -1) * r[:, 1])
        if sa <= 0: continue
        K = np.stack([x0 + r[:, 0] * res, y0 + r[:, 1] * res], 1)
        for _ in range(int(P.get("smooth_iters", 1))):
            nx_ = np.roll(K, -1, axis=0); K = np.stack([0.75 * K + 0.25 * nx_, 0.25 * K + 0.75 * nx_], axis=1).reshape(-1, 2)
        kept.append(K)
    kept.sort(key=lambda K: -abs(0.5 * np.sum(K[:, 0] * np.roll(K[:, 1], -1) - np.roll(K[:, 0], -1) * K[:, 1])))
    raw_ll = [np.stack(pj.inv(K[:, 0], K[:, 1]), 1) for K in kept]
    # simplify (km space tolerance), 4-dp, LODs, size budget
    budget = int(extra.get("file_budget_bytes", 1_500_000)); fracs = [0.25, 0.0625]; bpv = 19.5
    biome_reserve = 160_000
    max_verts = int((budget - biome_reserve - 20_000) / (bpv * (1 + sum(fracs))))
    eps = 0.15 * res
    # biome / relief
    base, rid = relief_noise(pj, x0, y0, res, nx, ny, seed, R)
    base = base * 2.5; rid = rid * 1.0
    up = 1 - np.exp(-np.clip(C, 0, None) / 140.0)
    rel = (0.12 + 0.55 * up + 0.6 * up * rid + 0.2 * up * (0.5 + 0.5 * base)).astype(np.float32); del base, rid, up
    if arcU is not None: rel = rel + 0.9 * np.clip(arcU, 0, 1.5)
    rel = np.where(L, rel, 0).astype(np.float32)          # land only; never feeds back into L
    bc = dict(extra.get("biome", {})); bz_used = bz
    while True:
        full = []
        for K, rl in zip(kept, raw_ll):
            q = rnd_ring(np.stack(pj.inv(*K[dp_tol(K, eps)].T), 1)); full.append(q if len(q) >= 3 else None)
        tot = sum(len(q) for q in full if q is not None)
        if tot <= max_verts or eps > 40: break
        eps *= 1.25
    idx = [i for i, q in enumerate(full) if q is not None]
    full = [full[i] for i in idx]; raw_ll = [raw_ll[i] for i in idx]; kept = [kept[i] for i in idx]
    bad = {i for pair in find_crossings(full) for i in pair}
    for i in bad: full[i] = rnd_ring(raw_ll[i])
    plkm = lambda r: np.stack(pj.fwd(r[:, 0], r[:, 1]), 1)
    levels_ = build_lods(full, fracs, vlog, pl=plkm)
    areas = [ring_area_km2(q, R) for q in full]; area_tot = float(sum(areas))
    # metrics on the final output
    gen_km = [plkm(q) for q in full]
    cov_g = rasterize_cov([np.stack([(K[:, 0] - x0) / res, (K[:, 1] - y0) / res], 1) for K in gen_km], ny, nx)
    sdf_g = sdf_from_mask(cov_g >= 0.5, res); del cov_g
    gp = np.vstack([densify_km(K, 1.0) for K in gen_km])
    ip = np.vstack([densify_km(K, 1.0) for K in km_rings])
    d1 = np.abs(bilin(sdf, (gp[:, 0] - x0) / res, (gp[:, 1] - y0) / res)); d2 = np.abs(bilin(sdf_g, (ip[:, 0] - x0) / res, (ip[:, 1] - y0) / res))
    dev = np.concatenate([d1, d2]); p95 = float(np.percentile(dev, 95)); dmax = float(dev.max())
    del sdf_g
    # pins on the output
    def inland(p):
        pt = np.array([[p["lon"], p["lat"]]])
        for q in full:
            if q[:, 0].min() <= p["lon"] <= q[:, 0].max() and q[:, 1].min() <= p["lat"] <= q[:, 1].max() and point_in_ring(pt, q)[0]: return True
        return False
    for p in pins["land"]: p["ok"] = inland(p)
    for p in pins["water"]: p["ok"] = not inland(p)
    for p in pins["coastal"]:
        p["dist_km"] = float(ring_dist_to_polylines(np.array([[p["x"], p["y"]]]), gen_km)[0]); p["ok_c"] = p["dist_km"] <= 3.0
    lp_ok = all(p["ok"] for p in pins["land"]); wp_ok = all(p["ok"] for p in pins["water"]); cp_ok = all(p["ok_c"] for p in pins["coastal"])
    for p in pins["land"]:
        if not p["ok"]: vlog("  pin FAIL land:", p["name"])
    for p in pins["water"]:
        if not p["ok"]: vlog("  pin FAIL water:", p["name"])
    for p in pins["coastal"]:
        vlog(f"  coastal pin {p['name']}: {p['dist_km']:.2f} km from coast {'ok' if p['ok_c'] else 'FAIL'}")
        if not p["ok_c"]: p["ok"] = False
    # biome
    bbox_ll = [min(float(q[:, 0].min()) for q in full), min(float(q[:, 1].min()) for q in full), max(float(q[:, 0].max()) for q in full), max(float(q[:, 1].max()) for q in full)]
    blob, bcls, bgeo = biome_v3(pj, x0, y0, res, rel * float(bc.get("elev_scale_m", 1000.0)), L, arid, bz_used, bbox_ll, bc, R, vlog, spec=climate_spec(state))
    stats.update(getattr(biome_v3, "last_report", {}))
    bj = dict(blob)
    label = P.get("label", V3_LABELS.get(tag, tag))
    phash = hashlib.sha1(json.dumps([P, [z.get("params") for z in zones], seed], sort_keys=True, default=str).encode()).hexdigest()[:12]
    gen = dict(name="coastgen", version=V3_VERSION, seed=seed, params_hash=phash, region=label, bbox=[round(v, 4) for v in bbox_ll], res_km=round(res, 3),
               radius_km=R, target_km2=round(intent_area), achieved_km2=round(area_tot), simplify_eps_km=round(eps, 3), params=P, genzones=len(zones),
               lod_vertices=[sum(len(r) for r in lv) for lv in levels_], coast_dev_p95_km=round(p95, 2), coast_dev_max_km=round(dmax, 2), capped_share=round(capped_share, 3))
    order = sorted(range(len(full)), key=lambda i: -areas[i]); fl = []
    for n_, i in enumerate(order):
        q = [lv[i] for lv in levels_]
        fl.append({"id": f"gen-{tag}-{n_+1:03d}", "type": "land", "kind": "generated", "name": f"{label} island {n_+1:03d}", "status": "proposed", "group": tag, "generator": "coastgen",
                   "seed": seed, "params": phash, "area_km2": round(areas[i]),
                   "bbox": [round(float(full[i][:, 0].min()), 4), round(float(full[i][:, 1].min()), 4), round(float(full[i][:, 0].max()), 4), round(float(full[i][:, 1].max()), 4)],
                   "lod": [[[float(x), float(y)] for x, y in r] for r in q]})
    if lakeM is not None:   # generated lake shores: enclosed water of the cleaned mask that sits on a drawn lake
        wl_, nw_ = label_mask(~L, 4)
        lov_ = np.bincount(wl_[lakeM].ravel(), minlength=nw_ + 1); cnt_ = np.bincount(wl_.ravel(), minlength=nw_ + 1)
        Kl = np.isin(wl_, [k for k in range(1, nw_ + 1) if k != wl_[0, 0] and cnt_[k] and lov_[k] >= 0.5 * cnt_[k]])
        D = Kl.copy()
        for _ in range(2): D = D | np.roll(D, 1, 0) | np.roll(D, -1, 0) | np.roll(D, 1, 1) | np.roll(D, -1, 1)
        nl_ = 0
        for r in marching_squares32(np.where(D, -C, -1.0).astype(np.float32), 0.0, Kl):
            sa = 0.5 * np.sum(r[:, 0] * np.roll(r[:, 1], -1) - np.roll(r[:, 0], -1) * r[:, 1])
            if sa <= 0: continue
            K = np.stack([x0 + r[:, 0] * res, y0 + r[:, 1] * res], 1)
            q0 = rnd_ring(np.stack(pj.inv(*K[dp_tol(K, eps)].T), 1)); q1 = rnd_ring(np.stack(pj.inv(*K[dp_tol(K, eps * 4)].T), 1)) if len(K) > 16 else q0
            if len(q0) < 4: continue
            nl_ += 1
            fl.append({"id": f"gen-{tag}-lake-{nl_:03d}", "type": "lake", "kind": "lake", "name": f"{label} lake {nl_:03d}", "group": tag,
                       "area_km2": round(ring_area_km2(q0, R)), "bbox": [round(float(q0[:, 0].min()), 4), round(float(q0[:, 1].min()), 4), round(float(q0[:, 0].max()), 4), round(float(q0[:, 1].max()), 4)],
                       "lod": [[[float(x), float(y)] for x, y in q] for q in (q0, q1, q1)]})
        stats["lakes"] = nl_; stats["lakes_drawn"] = len(lake_feats)
    result = dict(version=2, group=tag, label=label, generator=gen, features=fl, biome=bj)
    text = json.dumps(result, separators=(",", ":")); size = len(text.encode())
    vlog(f"  file {size/1e6:.3f} MB, verts/LOD {gen['lod_vertices']}, eps {eps:.2f} km")
    if size > budget: vlog(f"  WARN over budget {size} > {budget}")
    stats.update(area_km2=round(area_tot), coast_dev_p95_km=round(p95, 2), coast_dev_max_km=round(dmax, 2), topology_ok=True, topology_issues=[],
                 pins_land_ok=bool(lp_ok), pins_water_ok=bool(wp_ok), pins_coastal_ok=bool(cp_ok))
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    tmp = out_path + ".tmp"; open(tmp, "w").write(text); os.replace(tmp, out_path)
    if png_path:
        os.makedirs(os.path.dirname(os.path.abspath(png_path)), exist_ok=True)
        ins = []
        for it in (P.get("insets") or []):
            ins.append((it.get("name", "inset"), *[float(v) for v in pj.fwd(it["lon"], it["lat"])], float(it.get("half_km", 120))))
        if not ins:
            for z, zk in list(zip(zones, zones_km))[:2]:
                cx, cy = zk.mean(0); ins.append((z.get("name", "zone"), float(cx), float(cy), max(80.0, float(np.ptp(zk[:, 0]) * 0.7))))
        title = f"{label} seed {seed} | wander {P['wander_km']} km | {len(full)} isles | {area_tot/1e6:.3f}M km2 vs intent {intent_area/1e6:.3f}M ({100*(area_tot/intent_area-1):+.2f}%)"
        stxt = f"dev p95 {p95:.1f} / max {dmax:.1f} km | pins L{int(lp_ok)} W{int(wp_ok)} C{int(cp_ok)} | {size/1e6:.2f} MB"
        render_preview_v3(png_path, tag, pj, x0, y0, res, C, rel, L, km_rings, gen_km, pins, zones_km, ins, bcls, bgeo, title, stxt)
    return stats, t0

def main_v3(a):
    tag = a.group
    out = a.out or os.path.normpath(os.path.join(GEN_DIR, tag + ".json"))
    png = None if a.png == "none" else (a.png or f"/tmp/mapcheck/coastgen-{tag}-v3.png")
    extra = json.load(open(a.config)) if a.config else {}
    stats, t0 = run_group(tag, a.state, out, png, seed=a.seed, res_override=a.res, extra=extra)
    stats["runtime_s"] = round(time.time() - t0, 1); stats["peak_mb"] = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6)
    keys = ["group", "area_km2", "intent_area_km2", "coast_dev_p95_km", "coast_dev_max_km", "topology_ok", "topology_issues", "pins_land_ok", "pins_water_ok", "pins_coastal_ok", "runtime_s", "peak_mb",
            "lakes", "lakes_drawn", "highest_non_spire", "spire_cell_max_m", "spire_ice_cells", "pin_climate"]
    print(json.dumps({k: stats.get(k) for k in keys if k in stats or k in keys[:12]}))

GEN_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "library", "maps", "generated")

def main():
    ap = argparse.ArgumentParser(description="coastgen: v3 (--group <tag>, coast follows the intent coast) or v2 (--preset wi) -> library/maps/generated/<group>.json")
    ap.add_argument("--wild-group", default=None, help="v2 wild on any group tag, with protection mask + midline exclusion; config defaults to configs/<tag>-v2wild.json, output to /tmp/mapcheck")
    ap.add_argument("--group", default=None, help="v3: landmass tag (wi, haldei, asahi-arch, barsdorf, ne-peninsula)")
    ap.add_argument("--preset", default="wi", choices=["wi"], help="v2 (legacy) preset")
    ap.add_argument("--state", default=os.path.join(os.path.dirname(__file__), "..", "..", "..", "library", "maps", "map-state.json"), help="read-only: intent shapes + pins")
    ap.add_argument("--seed", type=int, default=None, help="v3: overrides state.gen seed; v2 default 7")
    ap.add_argument("--out", default=None, help="output JSON (default library/maps/generated/<group>.json)")
    ap.add_argument("--png", default=None, help="preview PNG (default /tmp/mapcheck/coastgen-<group>-v2.png; 'none' to skip)")
    ap.add_argument("--res", type=float, default=None)
    ap.add_argument("--noise-check", action="store_true", help="v3: verify world-coherent noise on overlapping bboxes and exit")
    ap.add_argument("--scale", type=int, default=1)
    ap.add_argument("--style", default=None, help="JSON dict of style overrides")
    ap.add_argument("--config", default=None, help="JSON file merged over the preset (arcs, style, biome, ...)")
    a = ap.parse_args()
    if a.noise_check:
        print(json.dumps(noise_coherence_check())); return
    if a.group: return main_v3(a)
    if a.wild_group: return main_wild(a)
    if a.seed is None: a.seed = 7
    cfg = build_wi_cfg(a.state)
    if a.config:
        over = json.load(open(a.config))
        for k, v in over.items():
            if k == "style": cfg["style"].update(v)
            else: cfg[k] = v
    if a.style: cfg["style"].update(json.loads(a.style))
    cfg["seed"] = a.seed
    if a.res: cfg["res_deg"] = a.res
    cfg["climate_spec"] = climate_spec(json.load(open(a.state)))
    group = cfg.get("group", a.preset)
    out = a.out or os.path.normpath(os.path.join(GEN_DIR, group + ".json"))
    png = None if a.png == "none" else (a.png or f"/tmp/mapcheck/coastgen-{group}-v2.png")
    cfg["output"] = dict(json=out, png=png, scale=a.scale)
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    if png: os.makedirs(os.path.dirname(png) or ".", exist_ok=True)
    res, stats = generate(cfg)
    print(json.dumps(stats))

def main_wild(a):
    tag = a.wild_group
    cfg, state, notes = build_group_wild_cfg(a.state, tag)
    cpath = a.config or os.path.join(os.path.dirname(__file__), "configs", f"{tag}-v2wild.json")
    if os.path.exists(cpath):
        for k, v in json.load(open(cpath)).items():
            if k == "style": cfg["style"].update(v)
            elif k == "bbox" or k.startswith("_"): continue
            else: cfg[k] = v
    if a.style: cfg["style"].update(json.loads(a.style))
    apply_protect_overrides(cfg)
    ps = cfg.get("protect_small_km2")
    if ps:   # islets below this area keep their drawn coast (wild noise erases thin/small intents); relief inland stays wild
        for it in cfg["intent"]:
            if ring_area_km2(it["ring"], cfg["radius_km"]) < ps:
                cfg["protect_zones"].append(dict(name=f"small:{it['name']}", ring=it["ring"], feather_km=cfg.get("protect_small_feather_km", 25.0),
                                                 sea_feather_km=cfg.get("protect_small_sea_feather_km", 0.0), strength=1.0,
                                                 detail_scale=_detail_scale(ring_area_km2(it["ring"], cfg["radius_km"]), cfg.get("protect_detail")), wander_km=cfg.get("protect_small_gate_km", 15.0)))
    if "intent_small_weight" in cfg:
        for it in cfg["intent"]:
            if it["weight"] < 1.0: it["weight"] = cfg["intent_small_weight"]
    if ps:   # protected islets leave the wild prior: their warped copy would survive as a ghost islet beside the protected one
        for it in cfg["intent"]:
            if ring_area_km2(it["ring"], cfg["radius_km"]) < ps: it["weight"] = 0.0
    if cfg.get("conjugate"):   # the facing plate must be inside the grid
        ext = cfg["conjugate"].get("reach_km", 700.0) / (cfg["radius_km"] * math.pi / 180)
        bb = cfg["bbox"]; cfg["bbox"] = [bb[0], max(-85.0, bb[1] - ext), bb[2], min(85.0, bb[3] + ext)] if cfg["conjugate"].get("side") == "both" else [bb[0], bb[1], bb[2], min(85.0, bb[3] + ext)]
    pins_by_id = {f["id"]: f["geometry"] for f in state["features"] if f.get("type") == "pin"}
    for z in cfg.get("protect_extra", []):   # wild-only protect disks around canon pins (the bake makes them intent)
        pg = pins_by_id[z["pin"]]; rk = z.get("radius_km", 35.0); k = cfg["radius_km"] * math.pi / 180
        ring = [[pg["lon"] + rk * math.cos(t) / (k * math.cos(math.radians(pg["lat"]))), pg["lat"] + rk * math.sin(t) / k] for t in np.linspace(0, 2 * math.pi, 48, endpoint=False)]
        cfg["protect_zones"].append(dict(name=z.get("name", z["pin"]), ring=ring, feather_km=z.get("feather_km", 60.0), sea_feather_km=z.get("sea_feather_km", 0.0),
                                         strength=z.get("strength", 1.0), wander_km=z.get("gate_km", 10.0)))
    for chk in cfg.get("land_checks", []):
        pg = pins_by_id[chk["pin"]]; chk["lon"], chk["lat"] = pg["lon"], pg["lat"]
    cfg["seed"] = a.seed if a.seed is not None else int((state.get("gen") or {}).get("world_seed", 7))
    if a.res: cfg["res_deg"] = a.res
    out = a.out or f"/tmp/mapcheck/{tag}-v2wild.json"
    png = None if a.png == "none" else (a.png or f"/tmp/mapcheck/{tag}-v2wild.png")
    if os.path.abspath(out).startswith(os.path.abspath(GEN_DIR)): raise SystemExit("wild previews never write library/maps/generated (the bake reads it); use --out /tmp/...")
    cfg["output"] = dict(json=out, png=png, scale=a.scale)
    for d in (out, png):
        if d: os.makedirs(os.path.dirname(d) or ".", exist_ok=True)
    for n in notes: vlog("  note: " + n)
    vlog(f"[{tag} wild] bbox {cfg['bbox']}, intents {len(cfg['intent'])}, arcs {len(cfg['arcs'])}, pins L{len(cfg['land_pins'])} W{len(cfg['water_pins'])} C{len(cfg['coastal_pins'])}, protect {[z['name'] for z in cfg['protect_zones']]}")
    cfg["climate_spec"] = climate_spec(state)
    res, stats = generate(cfg, log=vlog, prep=make_wild_prep(cfg, state))
    print(json.dumps(stats))

if __name__ == "__main__":
    main()
