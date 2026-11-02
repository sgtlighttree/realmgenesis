struct Globals {
    view_proj: mat4x4<f32>,
    sun_dir: vec4<f32>,      // xyz: direction towards the sun
    params: vec4<f32>,       // x: max land height (m), y: max ocean depth (m), z: show grid
    eye: vec4<f32>,
};

@group(0) @binding(0) var<uniform> g: Globals;

struct VsIn {
    @location(0) pos: vec3<f32>,
    @location(1) height_m: f32,
    @location(2) normal: vec3<f32>,
};

struct VsOut {
    @builtin(position) clip: vec4<f32>,
    @location(0) world: vec3<f32>,
    @location(1) height_m: f32,
    @location(2) normal: vec3<f32>,
};

@vertex
fn vs_main(v: VsIn) -> VsOut {
    var o: VsOut;
    o.clip = g.view_proj * vec4<f32>(v.pos, 1.0);
    o.world = v.pos;
    o.height_m = v.height_m;
    o.normal = v.normal;
    return o;
}

// Placeholder hypsometric palette. M3 replaces it with biome + climate colouring.
fn land_color(h: f32) -> vec3<f32> {
    let t = clamp(h / g.params.x, 0.0, 1.0);
    let lowland = vec3<f32>(0.20, 0.36, 0.17);
    let upland = vec3<f32>(0.47, 0.42, 0.26);
    let rock = vec3<f32>(0.42, 0.38, 0.35);
    let snow = vec3<f32>(0.93, 0.94, 0.96);
    var c = mix(lowland, upland, smoothstep(0.0, 0.35, t));
    c = mix(c, rock, smoothstep(0.35, 0.7, t));
    c = mix(c, snow, smoothstep(0.78, 0.9, t));
    return c;
}

fn ocean_color(h: f32) -> vec3<f32> {
    let t = clamp(-h / g.params.y, 0.0, 1.0);
    let shelf = vec3<f32>(0.11, 0.36, 0.48);
    let deep = vec3<f32>(0.02, 0.07, 0.18);
    return mix(shelf, deep, sqrt(t));
}

@fragment
fn fs_main(i: VsOut) -> @location(0) vec4<f32> {
    let up = normalize(i.world);
    let sun = normalize(g.sun_dir.xyz);
    let view = normalize(g.eye.xyz - i.world);
    var color: vec3<f32>;
    var n: vec3<f32>;
    if (i.height_m < 0.0) {
        n = up;
        color = ocean_color(i.height_m);
    } else {
        n = normalize(i.normal);
        color = land_color(i.height_m);
    }
    let diffuse = max(dot(n, sun), 0.0);
    // Soft terminator + ambient so the night side isn't pure black.
    let day = smoothstep(-0.15, 0.25, dot(up, sun));
    var lit = color * (0.08 + 0.92 * diffuse * day);
    if (i.height_m < 0.0) {
        let hlf = normalize(sun + view);
        lit += vec3<f32>(0.9, 0.85, 0.7) * pow(max(dot(n, hlf), 0.0), 120.0) * 0.6 * day;
    }
    // Thin atmosphere rim.
    let rim = pow(1.0 - max(dot(up, view), 0.0), 3.0);
    lit += vec3<f32>(0.30, 0.50, 0.90) * rim * 0.35 * (0.2 + 0.8 * day);
    return vec4<f32>(lit, 1.0);
}
