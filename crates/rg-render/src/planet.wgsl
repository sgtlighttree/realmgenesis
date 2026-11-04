// Planet renderer: a 3D globe and a 2D equirectangular map, sharing one layer texture.
//
// Layer texture: equirectangular RGBA8 sRGB albedo (lon -180..180 left→right, lat +90 top).
// Height texture: equirectangular R32Float metres above sea level, read with textureLoad
// (manual bilinear) since 32-bit float textures aren't filterable everywhere.

struct Globals {
    view_proj: mat4x4<f32>,
    sun_dir: vec4<f32>,  // xyz towards the sun (planet frame)
    eye: vec4<f32>,      // xyz camera position (globe)
    params: vec4<f32>,   // x: relief shading strength, y: ocean specular on/off, z: radius m, w: exaggeration (globe) / hillshade exaggeration (map)
    map: vec4<f32>,      // x: centre lon (rad), y: centre lat (rad), z: view width (rad of lon), w: aspect (w/h)
    lod: vec4<f32>,      // x: px per radian, y: cell px where tiles show their parent surface, z: morph range px, w: lat/lon grid spacing in degrees (0 = off)
    river: vec4<f32>,    // xyz: river tint (linear) over the coverage stored in the layer's alpha, w: viewport width px
};

@group(0) @binding(0) var<uniform> g: Globals;
@group(0) @binding(1) var albedo_tex: texture_2d<f32>;
@group(0) @binding(2) var albedo_smp: sampler;
@group(0) @binding(3) var height_tex: texture_2d<f32>;
// Displayed surface gradient per cube face, [dh/d east, dh/d north] (rg_terrain::gradient).
@group(0) @binding(4) var grad_tex: texture_2d_array<f32>;
@group(0) @binding(5) var grad_smp: sampler;

const PI: f32 = 3.14159265358979;

fn lonlat_to_uv(lon: f32, lat: f32) -> vec2<f32> {
    return vec2<f32>(fract(lon / (2.0 * PI) + 0.5), clamp(0.5 - lat / PI, 0.0, 1.0));
}

// Layer colour with rivers tinted in (coverage in the layer's alpha). Rivers no longer fade
// out when zoomed in: the drained detail network that used to carve them in 3D is gone
// (R2, s03), so this overlay is the only place they show.
fn albedo(lon: f32, lat: f32) -> vec3<f32> {
    let t = textureSampleLevel(albedo_tex, albedo_smp, lonlat_to_uv(lon, lat), 0.0);
    return mix(t.rgb, g.river.xyz, 1.0 - t.a);
}

// The local (east, north) frame, exactly as `rg_terrain::gradient::east_north` bakes it.
fn east_of(up: vec3<f32>) -> vec3<f32> {
    let e = cross(vec3<f32>(0.0, 0.0, 1.0), up);
    let l = length(e);
    if (l < 1e-6) {
        return vec3<f32>(1.0, 0.0, 0.0);
    }
    return e / l;
}

// Surface normal from the baked gradient at face position (s, t) on `face`, with the
// vertical exaggeration applied to the slope.
fn surface_normal(up: vec3<f32>, st: vec2<f32>, face: u32, exag: f32) -> vec3<f32> {
    let n = f32(textureDimensions(grad_tex).x);
    // Corner-aligned samples: s = 0 is the centre of texel 0, s = 1 the centre of the last.
    let tc = (st * (n - 1.0) + 0.5) / n;
    let grad = textureSampleLevel(grad_tex, grad_smp, tc, face, 0.0).xy;
    let east = east_of(up);
    let north = cross(up, east);
    return normalize(up - exag * (grad.x * east + grad.y * north));
}

fn height_at(uv: vec2<f32>) -> f32 {
    let dims = vec2<f32>(textureDimensions(height_tex));
    let p = uv * dims - 0.5;
    let i = floor(p);
    let f = p - i;
    let w = i32(dims.x);
    let h = i32(dims.y);
    let x0 = (i32(i.x) % w + w) % w;
    let x1 = (x0 + 1) % w;
    let y0 = clamp(i32(i.y), 0, h - 1);
    let y1 = clamp(i32(i.y) + 1, 0, h - 1);
    let a = textureLoad(height_tex, vec2<i32>(x0, y0), 0).r;
    let b = textureLoad(height_tex, vec2<i32>(x1, y0), 0).r;
    let c = textureLoad(height_tex, vec2<i32>(x0, y1), 0).r;
    let d = textureLoad(height_tex, vec2<i32>(x1, y1), 0).r;
    return mix(mix(a, b, f.x), mix(c, d, f.x), f.y);
}

// ---------------------------------------------------------------- lat/lon grid
//
// Anti-aliased ~1 px lines every `g.lod.w` degrees; `graticule` returns how far to mix
// towards white. It uses fwidth, so call it from uniform control flow (before any early
// return). The spacing always divides 360, so longitude wrap-around lands on a line.
fn grid_line(x: f32, n_wrap: f32) -> vec2<f32> {
    // x is in grid steps. Returns (line coverage 0..1, 1 if this is the zero line).
    // fwidth(x) spikes where atan2 wraps at the antimeridian; fract(x) stays continuous
    // there, so the smaller of the two is the true screen-space rate.
    let fw = max(min(fwidth(x), fwidth(fract(x))), 1e-6);
    let d = abs(fract(x + 0.5) - 0.5);
    var cover = clamp(1.5 - d / fw, 0.0, 1.0);
    // Lines closer than ~3 px just look like a tint: fade them out.
    cover = cover * (1.0 - smoothstep(0.2, 0.4, fw));
    let k = round(x);
    let zero = select(0.0, 1.0, k - n_wrap * round(k / n_wrap) == 0.0);
    return vec2<f32>(cover, zero);
}

fn graticule(lon: f32, lat: f32) -> f32 {
    let step = g.lod.w;
    if (step <= 0.0) {
        return 0.0;
    }
    let deg = 180.0 / PI;
    let lon_line = grid_line(lon * deg / step, 360.0 / step);
    let lat_line = grid_line(lat * deg / step, 1.0e9);
    // Meridians converge on the poles: fade them out rather than smear.
    let pole_fade = 1.0 - smoothstep(82.0, 87.0, abs(lat) * deg);
    let lon_a = lon_line.x * pole_fade * mix(0.35, 0.6, lon_line.y);
    let lat_a = lat_line.x * mix(0.35, 0.6, lat_line.y);
    return max(lon_a, lat_a);
}

// ---------------------------------------------------------------- globe (streamed tiles)
//
// Camera-relative: g.view_proj has the eye at the origin, and each tile instance carries
// (tile centre − eye) computed in f64 on the CPU. g.eye is the absolute eye (f32, ~0.5 m
// precision), used only for directions (local up, lat/lon), never for positions.
// g.params.w is the vertical exaggeration in globe mode.

struct TileIn {
    @location(0) pos: vec3<f32>,
    @location(1) height_m: f32,
    @location(2) face_uv: vec3<f32>,
    @location(3) parent_height_m: f32,
    @location(4) skirt: f32,
    @location(5) offset: vec3<f32>,
    @location(6) morph: f32,
    @location(7) skirt_m: f32,
    @location(8) cell_m: f32,
};

struct TileOut {
    @builtin(position) clip: vec4<f32>,
    @location(0) rel: vec3<f32>,
    @location(1) height_m: f32,
    @location(2) st: vec2<f32>,
    @location(3) up: vec3<f32>,
    @location(4) @interpolate(flat) face: u32,
};

@vertex
fn vs_tile(v: TileIn) -> TileOut {
    let local = v.offset + v.pos;              // sea-level point relative to the eye
    let up = normalize(g.eye.xyz + local);     // direction only: f32 is plenty
    let exag = g.params.w;
    // Geomorph (CDLOD-style): where this tile's cells project small enough that its parent
    // would be drawn instead, show the parent's surface exactly; blend to full detail as
    // cells grow. Neighbouring tiles of different levels then meet without steps. Combined
    // with the time-based morph-in for freshly streamed tiles.
    let dist = max(length(local), 1.0);
    let cell_px = v.cell_m / dist * g.lod.x;
    let dist_morph = clamp((cell_px - g.lod.y) / max(g.lod.z, 0.01), 0.0, 1.0);
    let h = mix(v.parent_height_m, v.height_m, min(v.morph, dist_morph));
    let lift = max(h, 0.0) * exag - v.skirt * v.skirt_m;
    var o: TileOut;
    o.rel = local + up * lift;
    o.clip = g.view_proj * vec4<f32>(o.rel, 1.0);
    o.height_m = h;
    o.st = v.face_uv.xy;
    o.face = u32(v.face_uv.z + 0.5);
    o.up = up;
    return o;
}

@fragment
fn fs_tile(i: TileOut) -> @location(0) vec4<f32> {
    let up = normalize(i.up);
    let lon = atan2(up.y, up.x);
    let lat = asin(clamp(up.z, -1.0, 1.0));
    let grid = graticule(lon, lat);
    var base = albedo(lon, lat);
    let normal = surface_normal(up, i.st, i.face, g.params.w);
    // Steep slopes show bare rock whatever the biome (the biome layer is ~km-scale).
    if (i.height_m > 0.0) {
        let slope = 1.0 - dot(normal, up);
        let rock = smoothstep(0.12, 0.35, slope) * g.params.x;
        base = mix(base, vec3<f32>(0.20, 0.18, 0.16), rock * 0.7);
    }
    let sun = normalize(g.sun_dir.xyz);
    let dist = length(i.rel);
    let view = -i.rel / max(dist, 1e-3);
    let ocean = i.height_m < 0.0;
    var n = up;
    if (!ocean) {
        n = normalize(mix(up, normal, g.params.x));
    }
    // Wrapped diffuse keeps the day side readable; exposure lifts the sRGB-authored
    // palettes (which are dark in linear space) back to how they look on a flat map.
    let diffuse = clamp((dot(n, sun) + 0.25) / 1.25, 0.0, 1.0);
    let day = smoothstep(-0.2, 0.2, dot(up, sun));
    var lit = base * (0.16 + 1.45 * diffuse * day);
    if (ocean && g.params.y > 0.5) {
        let hlf = normalize(sun + view);
        lit += vec3<f32>(1.0, 0.95, 0.85) * pow(max(dot(n, hlf), 0.0), 400.0) * 0.35 * day;
    }
    // Atmosphere: a rim at grazing angles from orbit, aerial haze with distance near the
    // ground (scaled to the planet radius so it works for any planet size).
    let rim = pow(1.0 - max(dot(up, view), 0.0), 3.0);
    let sky = vec3<f32>(0.42, 0.60, 0.92) * (0.25 + 0.75 * day);
    lit += vec3<f32>(0.30, 0.50, 0.90) * rim * 0.30 * (0.2 + 0.8 * day);
    // Haze only matters near the ground: from orbit a ray crosses little air (the rim
    // covers the limb). Fade it in below ~3% of the radius in altitude.
    let altitude = length(g.eye.xyz) - g.params.z;
    let near_ground = 1.0 - smoothstep(0.004 * g.params.z, 0.03 * g.params.z, altitude);
    let haze = (1.0 - exp(-dist / (g.params.z * 0.02))) * near_ground;
    lit = mix(lit, sky, clamp(haze, 0.0, 1.0) * 0.6);
    lit = mix(lit, vec3<f32>(1.0), grid);
    return vec4<f32>(lit, 1.0);
}

// ---------------------------------------------------------------- map

struct MapOut {
    @builtin(position) clip: vec4<f32>,
    @location(0) ndc: vec2<f32>,
};

@vertex
fn vs_map(@builtin(vertex_index) vi: u32) -> MapOut {
    // Fullscreen triangle.
    let xy = vec2<f32>(f32((vi << 1u) & 2u), f32(vi & 2u)) * 2.0 - 1.0;
    var o: MapOut;
    o.clip = vec4<f32>(xy, 0.0, 1.0);
    o.ndc = xy;
    return o;
}

@fragment
fn fs_map(i: MapOut) -> @location(0) vec4<f32> {
    let width = g.map.z;
    let height = width / g.map.w;
    let lon = g.map.x + i.ndc.x * 0.5 * width;
    let lat = g.map.y + i.ndc.y * 0.5 * height;
    let grid = graticule(lon, lat);
    if (abs(lat) > PI * 0.5) {
        return vec4<f32>(0.055, 0.067, 0.078, 1.0);
    }
    let uv = lonlat_to_uv(lon, lat);
    var c = albedo(lon, lat);
    // Hillshade from the height texture (land only), lit from the north-west as on
    // printed maps.
    let h = height_at(uv);
    if (h > 0.0 && g.params.x > 0.0) {
        let dims = vec2<f32>(textureDimensions(height_tex));
        let du = 1.0 / dims.x;
        let dv = 1.0 / dims.y;
        let px_m = g.params.z * 2.0 * PI / dims.x; // metres per texel along the equator
        let coslat = max(cos(lat), 0.05);
        let dx = (height_at(uv + vec2<f32>(du, 0.0)) - height_at(uv - vec2<f32>(du, 0.0))) / (2.0 * px_m * coslat);
        let dy = (height_at(uv - vec2<f32>(0.0, dv)) - height_at(uv + vec2<f32>(0.0, dv))) / (2.0 * px_m);
        let ex = g.params.w;
        let n = normalize(vec3<f32>(-dx * ex, -dy * ex, 1.0));
        let light = normalize(vec3<f32>(-0.6, 0.6, 0.7));
        let shade = dot(n, light) / light.z; // 1.0 on flat ground
        c = c * mix(1.0, clamp(shade, 0.35, 1.6), g.params.x);
    }
    c = mix(c, vec3<f32>(1.0), grid);
    return vec4<f32>(c, 1.0);
}
