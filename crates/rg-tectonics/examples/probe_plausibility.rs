//! R3.1a plausibility gates (design §Gates 2), measured on the tectonic snapshot against the
//! **emergent** coast (after the sea-level bisection), not the crust edge.
//!
//! `cargo run -p rg-tectonics --release --example probe_plausibility -- [seed] [res] [active_share]`

use rg_tectonics::coast::margin;
use rg_tectonics::grid::Grid;
use rg_tectonics::plates::Kind;
use rg_tectonics::{TectonicsParams, generate_with_model};
use std::f64::consts::PI;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let seed: u64 = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(1);
    let res: u32 = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(320);
    let share: f32 = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(0.25);
    let p = TectonicsParams { seed, res, active_share: share, ..Default::default() };
    let t0 = std::time::Instant::now();
    let (snap, pm, st) = generate_with_model(&p);
    let ms = t0.elapsed().as_millis();
    let g = Grid::new(res);
    let r = p.radius_m;
    let h = |a: usize| snap.height_m.data[g.node_sample[a] as usize];
    let mg = snap.margin.as_ref().unwrap();
    let m_of = |a: usize| mg.data[g.node_sample[a] as usize];

    // Emergent coast nodes (land next to sea), with their margin class.
    let coast: Vec<usize> = (0..g.n).filter(|&a| h(a) > 0.0 && g.neighbors(a).iter().any(|&b| h(b as usize) <= 0.0)).collect();
    // Outer coasts only: land next to sea that connects to the open ocean (oceanic crust).
    // Enclosed seas inside continental crust inherit a class but aren't margins.
    let crust_at = |a: usize| snap.crust.data[g.node_sample[a] as usize];
    let mut open_sea = vec![false; g.n];
    let mut stack: Vec<usize> = (0..g.n).filter(|&a| h(a) <= 0.0 && crust_at(a) == 0).collect();
    for &a in &stack {
        open_sea[a] = true;
    }
    while let Some(a) = stack.pop() {
        for &b in g.neighbors(a) {
            let b = b as usize;
            if !open_sea[b] && h(b) <= 0.0 {
                open_sea[b] = true;
                stack.push(b);
            }
        }
    }
    let all_coast = coast.len();
    let coast: Vec<usize> = coast.into_iter().filter(|&a| g.neighbors(a).iter().any(|&b| open_sea[b as usize])).collect();
    let crust_edge: Vec<u32> = (0..g.n).filter(|&a| crust_at(a) == 1 && g.neighbors(a).iter().any(|&b| crust_at(b as usize) == 0)).map(|a| a as u32).collect();
    let ne = g.nearest_source(&crust_edge, PI as f32);
    println!("outer coast share of emergent coast nodes: {:.3}", coast.len() as f64 / all_coast.max(1) as f64);
    let sub: Vec<u32> = pm.bnodes.iter().filter(|b| b.kind == Kind::Subduction).map(|b| b.node).collect();
    let near_sub = g.nearest_source(&sub, PI as f32);
    let (mut act, mut act_ok, mut pas, mut pas_bad, mut none) = (0.0, 0.0, 0.0, 0.0, 0.0);
    let (mut act_h, mut pas_h) = (Vec::new(), Vec::new());
    let land_near = |a: usize| -> f32 {
        // Highest land within ~250 km inland (a 3–8 hop ball depending on res).
        let hops = ((250_000.0 / r) / g.spacing_rad()).ceil().max(1.0) as usize;
        let mut best = h(a);
        let mut frontier = vec![a];
        let mut seen = vec![a];
        for _ in 0..hops {
            let mut next = Vec::new();
            for &x in &frontier {
                for &y in g.neighbors(x) {
                    let y = y as usize;
                    if !seen.contains(&y) && h(y) > 0.0 {
                        seen.push(y);
                        next.push(y);
                        best = best.max(h(y));
                    }
                }
            }
            frontier = next;
        }
        best
    };
    for (i, &a) in coast.iter().enumerate() {
        let w = g.weight[a] as f64;
        let d = near_sub.dist_rad[a] as f64 * r;
        match m_of(a) {
            margin::ACTIVE => {
                act += w;
                if d < 300_000.0 {
                    act_ok += w;
                }
                if i % 7 == 0 {
                    act_h.push(land_near(a));
                }
            }
            margin::PASSIVE => {
                pas += w;
                if d < 400_000.0 {
                    pas_bad += w;
                }
                if i % 7 == 0 {
                    pas_h.push(land_near(a));
                }
            }
            _ => none += w,
        }
    }
    // Diagnostic: distance from active coast to the nearest subduction node, and to the
    // nearest crust edge (convergent or not).
    let mut dact: Vec<f32> = coast.iter().filter(|&&a| m_of(a) == margin::ACTIVE).map(|&a| (near_sub.dist_rad[a] as f64 * r / 1000.0) as f32).collect();
    dact.sort_by(f32::total_cmp);
    let q = |v: &Vec<f32>, f: f64| v.get(((v.len() as f64 - 1.0) * f) as usize).copied().unwrap_or(0.0);
    println!("active coast -> trench km: p10 {:.0} p25 {:.0} p50 {:.0} p75 {:.0} p90 {:.0}", q(&dact, 0.1), q(&dact, 0.25), q(&dact, 0.5), q(&dact, 0.75), q(&dact, 0.9));
    let mut dedge: Vec<f32> = coast.iter().filter(|&&a| m_of(a) == margin::ACTIVE).map(|&a| (ne.dist_rad[a] as f64 * r / 1000.0) as f32).collect();
    dedge.sort_by(f32::total_cmp);
    println!("active coast -> crust edge km: p10 {:.0} p50 {:.0} p90 {:.0}", q(&dedge, 0.1), q(&dedge, 0.5), q(&dedge, 0.9));
    let median = |v: &mut Vec<f32>| {
        v.sort_by(f32::total_cmp);
        v.get(v.len() / 2).copied().unwrap_or(0.0)
    };

    // Ridge symmetry: ridge nodes across a seaway axis, distance to the nearest emergent coast
    // on each plate's side.
    let k = snap.plates.len();
    let plate = |a: usize| pm.plate_of[a] as usize;
    let mut near_by_plate = Vec::with_capacity(k);
    for pl in 0..k {
        let src: Vec<u32> = coast.iter().filter(|&&a| plate(a) == pl).map(|&a| a as u32).collect();
        near_by_plate.push(if src.is_empty() { None } else { Some(g.nearest_source(&src, PI as f32)) });
    }
    let (mut sym_n, mut sym_ok) = (0usize, 0usize);
    for b in pm.bnodes.iter().filter(|b| b.kind == Kind::Ridge) {
        let a = b.node as usize;
        let own = plate(a);
        let Some(other) = g.neighbors(a).iter().map(|&x| plate(x as usize)).find(|&q| q != own) else { continue };
        let (Some(na), Some(nb)) = (&near_by_plate[own], &near_by_plate[other]) else { continue };
        let (da, db) = (na.dist_rad[a] as f64, nb.dist_rad[a] as f64);
        if da * r > 3_000_000.0 || db * r > 3_000_000.0 {
            continue; // not a seaway between two coasts
        }
        sym_n += 1;
        if (da - db).abs() / (da + db).max(1e-12) < 0.15 {
            sym_ok += 1;
        }
    }

    let counts = |kind: Kind| pm.bnodes.iter().filter(|b| b.kind == kind).count();
    println!("seed {seed} res {res} active_share {share}  tectonics {ms} ms  plates {}", k);
    println!(
        "seaways: opening {} closing {} one-sided {}   coast active share (crust edge) {:.3}",
        st.seaways[0], st.seaways[1], st.seaways[2], st.active_share
    );
    println!(
        "motion: rms residual {:.1} mm/yr, wrong-sign share {:.3}",
        st.rms_residual_mm_yr, st.wrong_sign_share
    );
    println!(
        "boundary nodes: subduction {} collision {} ridge {} rift {} transform {}",
        counts(Kind::Subduction),
        counts(Kind::Collision),
        counts(Kind::Ridge),
        counts(Kind::Rift),
        counts(Kind::Transform)
    );
    let tot = act + pas + none;
    println!(
        "emergent coast: active {:.3} passive {:.3} unclassified {:.3}",
        act / tot,
        pas / tot,
        none / tot
    );
    println!(
        "GATE active coast with trench < 300 km: {:.3} (>= 0.90)   passive coast with trench < 400 km: {:.3} (~0)",
        act_ok / act.max(1e-12),
        pas_bad / pas.max(1e-12)
    );
    println!("GATE ridge symmetry |dA-dB|/(dA+dB) < 0.15: {:.3} of {sym_n} seaway ridge nodes", sym_ok as f64 / sym_n.max(1) as f64);
    println!(
        "highest land within 250 km of coast, median: active {:.0} m, passive {:.0} m",
        median(&mut act_h),
        median(&mut pas_h)
    );
    println!("land fraction {:.3}", snap.land_fraction);

    // RG_MAP=<png>: equirect map. Grey relief; margin tint (red active, teal passive);
    // boundary kinds: subduction red, collision orange, ridge blue, rift purple, transform green.
    if let Ok(path) = std::env::var("RG_MAP") {
        let mut kind_at = vec![255u8; g.n];
        for &a in &coast {
            if m_of(a) == margin::ACTIVE && near_sub.dist_rad[a] as f64 * r >= 300_000.0 {
                kind_at[a] = 9; // failing active coast
            }
            if m_of(a) == margin::PASSIVE && near_sub.dist_rad[a] as f64 * r < 400_000.0 {
                kind_at[a] = 8; // passive coast near a trench
            }
        }
        for b in &pm.bnodes {
            kind_at[b.node as usize] = b.kind as u8;
        }
        let (w, hh) = (1440usize, 720usize);
        let mut img = vec![0u8; w * hh * 3];
        for y in 0..hh {
            let lat = std::f64::consts::FRAC_PI_2 - (y as f64 + 0.5) / hh as f64 * PI;
            for x in 0..w {
                let lon = -PI + (x as f64 + 0.5) / w as f64 * 2.0 * PI;
                let d = glam::DVec3::new(lat.cos() * lon.cos(), lat.cos() * lon.sin(), lat.sin());
                let a = g.nearest_node(d);
                let hv = h(a);
                let gr = if hv > 0.0 { 140.0 + 100.0 * (hv / 4000.0).clamp(0.0, 1.0) } else { 30.0 + 60.0 * (1.0 + hv / 6000.0).clamp(0.0, 1.0) };
                let mut c = [gr, gr, gr];
                let tint = match m_of(a) {
                    margin::ACTIVE => Some([230.0, 80.0, 50.0]),
                    margin::PASSIVE => Some([60.0, 170.0, 190.0]),
                    _ => None,
                };
                if let Some(t) = tint {
                    for i in 0..3 {
                        c[i] = 0.5 * c[i] + 0.5 * t[i];
                    }
                }
                let k = kind_at[a];
                if k != 255 {
                    c = match k {
                        0 => [255.0, 150.0, 0.0],
                        1 => [200.0, 0.0, 0.0],
                        2 => [40.0, 90.0, 255.0],
                        3 => [170.0, 60.0, 220.0],
                        9 => [255.0, 0.0, 255.0],
                        8 => [255.0, 255.0, 0.0],
                        _ => [40.0, 200.0, 60.0],
                    };
                }
                let o = (y * w + x) * 3;
                for i in 0..3 {
                    img[o + i] = c[i] as u8;
                }
            }
        }
        let file = std::io::BufWriter::new(std::fs::File::create(&path).unwrap());
        let mut enc = png::Encoder::new(file, w as u32, hh as u32);
        enc.set_color(png::ColorType::Rgb);
        enc.set_depth(png::BitDepth::Eight);
        enc.write_header().unwrap().write_image_data(&img).unwrap();
        println!("map: {path}");
    }
}
