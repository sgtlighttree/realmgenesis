//! Decode a PNG the way a heightmap import does and print its size, timing and a grey
//! histogram (s07). Usage: `probe_import <file.png> [max_side=4096]`. Run under
//! `/usr/bin/time -l` for the peak memory.
use rg_pipeline::import;
use rg_pipeline::mask::Channel;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let path = args.get(1).expect("usage: probe_import <file.png> [max_side]");
    let max_side: u32 = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(4096);
    let bytes = std::fs::read(path).expect("read");
    let (sw, sh) = import::png_size(&bytes).expect("png");
    let step = sw.max(sh).div_ceil(max_side).max(1);
    let t = std::time::Instant::now();
    let (d, w, h, used) = import::raw_from_png(&bytes, Channel::Brightness, None, step).expect("decode");
    println!("{sw}x{sh} -> {w}x{h} (step {step}), {used:?}, {:.2} s", t.elapsed().as_secs_f64());
    let mut hist = [0u64; 256];
    for v in &d {
        hist[(v.clamp(0.0, 1.0) * 255.0).round() as usize] += 1;
    }
    let n = d.len() as f64;
    for (g, c) in hist.iter().enumerate().filter(|(_, c)| **c > 0) {
        if g % 4 == 0 {
            println!("grey {g:3}: {:6.3}%", 100.0 * *c as f64 / n);
        }
    }
}
