use rg_pipeline::*;
use std::time::Instant;
fn main() {
    let res: u32 = std::env::args().nth(1).unwrap().parse().unwrap();
    let w: u32 = std::env::args().nth(2).unwrap().parse().unwrap();
    let mut p = WorldParams::default();
    p.tectonics.res = res;
    let t = Instant::now();
    let tect = rg_tectonics::generate(&p.tectonics);
    let tm = t.elapsed().as_secs_f64();
    drop(tect);
    let t = Instant::now();
    let world = generate(&p);
    let total = t.elapsed().as_secs_f64();
    println!("res {res}: tect {:.2}s climate {:.2}s total {:.2}s", world.timings.tectonics_ms/1e3, world.timings.climate_ms/1e3, total);
    let _ = tm;
    let t = Instant::now();
    let eq = Equirect::new(&world, w, w / 2);
    let tn = t.elapsed().as_secs_f64();
    let t = Instant::now();
    let b = eq.bake(&world, Layer::Satellite);
    let tb = t.elapsed().as_secs_f64();
    println!("  equirect {w}x{}: new {:.2}s, bake sat {:.2}s ({} bytes)", w/2, tn, tb, b.len());
    if let Ok(st) = std::fs::read_to_string("/proc/self/status") {
        for l in st.lines() { if l.starts_with("VmHWM") { println!("  {l}"); } }
    }
}
