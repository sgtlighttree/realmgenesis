// Planet renderer: a 3D globe and a 2D equirectangular map, sharing one layer texture.
//
// Layer texture: equirectangular RGBA8 sRGB albedo (lon -180..180 left→right, lat +90 top).
// Height texture: equirectangular R32Float metres above sea level, read with textureLoad
// (manual bilinear) since 32-bit float textures aren't filterable everywhere.

struct Globals {
    view_proj: mat4x4<f32>,
    sun_dir: vec4<f32>,  // xyz towards the sun (planet frame)
    eye: vec4<f32>,      // xyz camera position (globe)
    params: vec4<f32>,   // x: relief shading strength, y: ocean specular on/off, z: radius m, w: map exaggeration
    map: vec4<f32>,      // x: centre lon (rad), y: centre lat (rad), z: view width (rad of lon), w: aspect (w/h)
};

@group(0) @binding(0) var<uniform> g: Globals;
@group(0) @binding(1) var albedo_tex: texture_2d<f32>;
@group(0) @binding(2) var albedo_smp: sampler;
@group(0) @binding(3) var height_tex: texture_2d<f32>;

const PI: f32 = 3.14159265358979;

fn lonlat_to_uv(lon: f32, lat: f32) -> vec2<f32> {
    return vec2<f32>(fract(lon / (2.0 * PI) + 0.5), clamp(0.5 - lat / PI, 0.0, 1.0));
}

fn albedo(lon: f32, lat: f32) -> vec3<f32> {
    return textureSampleLevel(albedo_tex, albedo_smp, lonlat_to_uv(lon, lat), 0.0).rgb;
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

// ---------------------------------------------------------------- globe

struct VsIn {
    @location(0) pos: vec3<f32>,
    @location(1) height_m: f32,
    @location(2) normal: vec3<f32>,
};

struct GlobeOut {
    @builtin(position) clip: vec4<f32>,
    @location(0) world: vec3<f32>,
    @location(1) height_m: f32,
    @location(2) normal: vec3<f32>,
};

@vertex
fn vs_globe(v: VsIn) -> GlobeOut {
    var o: GlobeOut;
    o.clip = g.view_proj * vec4<f32>(v.pos, 1.0);
    o.world = v.pos;
    o.height_m = v.height_m;
    o.normal = v.normal;
    return o;
}

@fragment
fn fs_globe(i: GlobeOut) -> @location(0) vec4<f32> {
    let up = normalize(i.world);
    let lon = atan2(up.y, up.x);
    let lat = asin(clamp(up.z, -1.0, 1.0));
    let base = albedo(lon, lat);
    let sun = normalize(g.sun_dir.xyz);
    let view = normalize(g.eye.xyz - i.world);
    let ocean = i.height_m < 0.0;
    var n = up;
    if (!ocean) {
        n = normalize(mix(up, normalize(i.normal), g.params.x));
    }
    let diffuse = max(dot(n, sun), 0.0);
    let day = smoothstep(-0.15, 0.25, dot(up, sun));
    var lit = base * (0.10 + 0.90 * diffuse * day);
    if (ocean && g.params.y > 0.5) {
        let hlf = normalize(sun + view);
        lit += vec3<f32>(0.9, 0.85, 0.7) * pow(max(dot(n, hlf), 0.0), 120.0) * 0.5 * day;
    }
    let rim = pow(1.0 - max(dot(up, view), 0.0), 3.0);
    lit += vec3<f32>(0.30, 0.50, 0.90) * rim * 0.30 * (0.2 + 0.8 * day);
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
    return vec4<f32>(c, 1.0);
}
