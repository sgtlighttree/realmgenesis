//! R3.2: a small SVG reader for coastline masks (no SVG crate offline; the subset coast
//! drawings need).
//!
//! **Every filled shape is land** (even-odd fill), whatever its colour; shapes with
//! `fill="none"` are skipped. Supported: `<path d>` (M L H V C S Q T Z; arcs `A` become a
//! straight segment to their end point), `<polygon>`, `<polyline>` (filled), `<rect>`,
//! `<circle>`, `<ellipse>`, and `transform` (matrix, translate, scale, rotate) on shapes and
//! `<g>` groups. The viewBox (or width/height) is the drawing's extent.
//!
//! Placement: `WorldEquirect` stretches the viewBox over the whole sphere; for a `Patch`,
//! `km_per_px` means **km per viewBox unit**. The shapes are rasterised (exact span coverage,
//! 4 sub-rows per pixel) and then sampled like a PNG.

use rg_tectonics::mask::Placement;

/// Longest side of the rasterised mask, pixels.
const MAX_SIDE: f64 = 4096.0;

type Pt = (f64, f64);
type Mat = [f64; 6]; // a b c d e f: x' = a x + c y + e, y' = b x + d y + f

const IDENTITY: Mat = [1.0, 0.0, 0.0, 1.0, 0.0, 0.0];

fn mul(m: &Mat, n: &Mat) -> Mat {
    [
        m[0] * n[0] + m[2] * n[1],
        m[1] * n[0] + m[3] * n[1],
        m[0] * n[2] + m[2] * n[3],
        m[1] * n[2] + m[3] * n[3],
        m[0] * n[4] + m[2] * n[5] + m[4],
        m[1] * n[4] + m[3] * n[5] + m[5],
    ]
}

fn apply(m: &Mat, p: Pt) -> Pt {
    (m[0] * p.0 + m[2] * p.1 + m[4], m[1] * p.0 + m[3] * p.1 + m[5])
}

/// Numbers in an attribute or path string.
fn numbers(s: &str) -> Vec<f64> {
    let mut out = Vec::new();
    let b = s.as_bytes();
    let mut i = 0;
    while i < b.len() {
        let c = b[i];
        if c.is_ascii_digit() || c == b'-' || c == b'+' || c == b'.' {
            let start = i;
            i += 1;
            let mut seen_dot = c == b'.';
            while i < b.len() {
                let d = b[i];
                if d.is_ascii_digit() {
                    i += 1;
                } else if d == b'.' && !seen_dot {
                    seen_dot = true;
                    i += 1;
                } else if (d == b'e' || d == b'E') && i + 1 < b.len() {
                    i += 1;
                    if b[i] == b'-' || b[i] == b'+' {
                        i += 1;
                    }
                } else {
                    break;
                }
            }
            if let Ok(v) = s[start..i].parse() {
                out.push(v);
            }
        } else {
            i += 1;
        }
    }
    out
}

fn attr<'a>(tag: &'a str, name: &str) -> Option<&'a str> {
    let mut rest = tag;
    while let Some(pos) = rest.find(name) {
        let before_ok = pos == 0 || rest.as_bytes()[pos - 1].is_ascii_whitespace();
        let after = rest[pos + name.len()..].trim_start();
        if before_ok && after.starts_with('=') {
            let v = after[1..].trim_start();
            let q = v.chars().next()?;
            if q == '"' || q == '\'' {
                let end = v[1..].find(q)?;
                return Some(&v[1..1 + end]);
            }
        }
        rest = &rest[pos + name.len()..];
    }
    None
}

fn parse_transform(s: &str) -> Mat {
    let mut m = IDENTITY;
    let mut rest = s;
    while let Some(open) = rest.find('(') {
        let name = rest[..open].trim().trim_start_matches(',').trim();
        let Some(close) = rest[open..].find(')') else { break };
        let v = numbers(&rest[open + 1..open + close]);
        let t = match (name, v.as_slice()) {
            ("matrix", [a, b, c, d, e, f, ..]) => [*a, *b, *c, *d, *e, *f],
            ("translate", [x]) => [1.0, 0.0, 0.0, 1.0, *x, 0.0],
            ("translate", [x, y, ..]) => [1.0, 0.0, 0.0, 1.0, *x, *y],
            ("scale", [s]) => [*s, 0.0, 0.0, *s, 0.0, 0.0],
            ("scale", [x, y, ..]) => [*x, 0.0, 0.0, *y, 0.0, 0.0],
            ("rotate", [a, rest @ ..]) => {
                let (s, c) = a.to_radians().sin_cos();
                let r = [c, s, -s, c, 0.0, 0.0];
                if let [cx, cy, ..] = rest {
                    mul(&mul(&[1.0, 0.0, 0.0, 1.0, *cx, *cy], &r), &[1.0, 0.0, 0.0, 1.0, -cx, -cy])
                } else {
                    r
                }
            }
            _ => IDENTITY,
        };
        m = mul(&m, &t);
        rest = &rest[open + close + 1..];
    }
    m
}

/// Flatten path data into closed polygons (one per subpath).
fn path_polys(d: &str) -> Vec<Vec<Pt>> {
    // Tokenise into commands and numbers.
    let mut toks: Vec<(char, Vec<f64>)> = Vec::new();
    let mut cur = String::new();
    let mut cmd = None;
    for ch in d.chars() {
        if ch.is_ascii_alphabetic() && ch != 'e' && ch != 'E' {
            if let Some(c) = cmd {
                toks.push((c, numbers(&cur)));
            }
            cmd = Some(ch);
            cur.clear();
        } else {
            cur.push(ch);
        }
    }
    if let Some(c) = cmd {
        toks.push((c, numbers(&cur)));
    }
    let mut polys = Vec::new();
    let mut poly: Vec<Pt> = Vec::new();
    let (mut p, mut start) = ((0.0, 0.0), (0.0, 0.0));
    let mut last_ctrl: Option<Pt> = None;
    let mut flush = |poly: &mut Vec<Pt>| {
        if poly.len() >= 3 {
            polys.push(std::mem::take(poly));
        } else {
            poly.clear();
        }
    };
    let cubic = |poly: &mut Vec<Pt>, p0: Pt, c1: Pt, c2: Pt, p1: Pt| {
        for i in 1..=8 {
            let t = i as f64 / 8.0;
            let u = 1.0 - t;
            poly.push((
                u * u * u * p0.0 + 3.0 * u * u * t * c1.0 + 3.0 * u * t * t * c2.0 + t * t * t * p1.0,
                u * u * u * p0.1 + 3.0 * u * u * t * c1.1 + 3.0 * u * t * t * c2.1 + t * t * t * p1.1,
            ));
        }
    };
    for (c, v) in toks {
        let rel = c.is_ascii_lowercase();
        let off = |q: Pt, x: f64, y: f64| if rel { (q.0 + x, q.1 + y) } else { (x, y) };
        match c.to_ascii_uppercase() {
            'M' => {
                flush(&mut poly);
                for (k, xy) in v.chunks_exact(2).enumerate() {
                    p = off(p, xy[0], xy[1]);
                    if k == 0 {
                        start = p;
                    }
                    poly.push(p);
                }
                last_ctrl = None;
            }
            'L' | 'T' => {
                for xy in v.chunks_exact(2) {
                    p = off(p, xy[0], xy[1]);
                    poly.push(p);
                }
                last_ctrl = None;
            }
            'H' => {
                for &x in &v {
                    p = (if rel { p.0 + x } else { x }, p.1);
                    poly.push(p);
                }
            }
            'V' => {
                for &y in &v {
                    p = (p.0, if rel { p.1 + y } else { y });
                    poly.push(p);
                }
            }
            'C' => {
                for s in v.chunks_exact(6) {
                    let (c1, c2, e) = (off(p, s[0], s[1]), off(p, s[2], s[3]), off(p, s[4], s[5]));
                    cubic(&mut poly, p, c1, c2, e);
                    last_ctrl = Some(c2);
                    p = e;
                }
            }
            'S' => {
                for s in v.chunks_exact(4) {
                    let c1 = last_ctrl.map_or(p, |q| (2.0 * p.0 - q.0, 2.0 * p.1 - q.1));
                    let (c2, e) = (off(p, s[0], s[1]), off(p, s[2], s[3]));
                    cubic(&mut poly, p, c1, c2, e);
                    last_ctrl = Some(c2);
                    p = e;
                }
            }
            'Q' => {
                for s in v.chunks_exact(4) {
                    let (q, e) = (off(p, s[0], s[1]), off(p, s[2], s[3]));
                    let c1 = (p.0 + 2.0 / 3.0 * (q.0 - p.0), p.1 + 2.0 / 3.0 * (q.1 - p.1));
                    let c2 = (e.0 + 2.0 / 3.0 * (q.0 - e.0), e.1 + 2.0 / 3.0 * (q.1 - e.1));
                    cubic(&mut poly, p, c1, c2, e);
                    p = e;
                }
                last_ctrl = None;
            }
            'A' => {
                for s in v.chunks_exact(7) {
                    p = off(p, s[5], s[6]);
                    poly.push(p);
                }
                last_ctrl = None;
            }
            'Z' => {
                flush(&mut poly);
                p = start;
                last_ctrl = None;
            }
            _ => {}
        }
    }
    flush(&mut poly);
    polys
}

fn ellipse(cx: f64, cy: f64, rx: f64, ry: f64) -> Vec<Pt> {
    (0..64).map(|i| {
        let a = i as f64 / 64.0 * std::f64::consts::TAU;
        (cx + rx * a.cos(), cy + ry * a.sin())
    }).collect()
}

/// Shapes (in user units, transforms applied) and the viewBox `(x, y, w, h)`.
fn parse(text: &str) -> Result<(Vec<Vec<Pt>>, [f64; 4]), String> {
    let mut view: Option<[f64; 4]> = None;
    let mut stack: Vec<Mat> = vec![IDENTITY];
    let mut polys = Vec::new();
    let mut rest = text;
    while let Some(open) = rest.find('<') {
        let Some(close) = rest[open..].find('>') else { break };
        let tag = &rest[open + 1..open + close];
        rest = &rest[open + close + 1..];
        if tag.starts_with('!') || tag.starts_with('?') {
            continue;
        }
        if let Some(name) = tag.strip_prefix('/') {
            if name.trim() == "g" && stack.len() > 1 {
                stack.pop();
            }
            continue;
        }
        let name = tag.split(|c: char| c.is_whitespace() || c == '/').next().unwrap_or("");
        let self_closing = tag.trim_end().ends_with('/');
        let local = attr(tag, "transform").map_or(IDENTITY, parse_transform);
        let m = mul(stack.last().unwrap(), &local);
        if name == "svg" && view.is_none() {
            view = attr(tag, "viewBox").map(numbers).filter(|v| v.len() >= 4).map(|v| [v[0], v[1], v[2], v[3]]).or_else(|| {
                let w = attr(tag, "width").and_then(|s| numbers(s).first().copied())?;
                let h = attr(tag, "height").and_then(|s| numbers(s).first().copied())?;
                Some([0.0, 0.0, w, h])
            });
            continue;
        }
        if name == "g" {
            if !self_closing {
                stack.push(m);
            }
            continue;
        }
        if attr(tag, "fill") == Some("none") {
            continue;
        }
        let num = |k: &str| attr(tag, k).and_then(|s| numbers(s).first().copied()).unwrap_or(0.0);
        let shapes: Vec<Vec<Pt>> = match name {
            "path" => attr(tag, "d").map(path_polys).unwrap_or_default(),
            "polygon" | "polyline" => {
                let v = attr(tag, "points").map(numbers).unwrap_or_default();
                vec![v.chunks_exact(2).map(|p| (p[0], p[1])).collect()]
            }
            "rect" => {
                let (x, y, w, h) = (num("x"), num("y"), num("width"), num("height"));
                vec![vec![(x, y), (x + w, y), (x + w, y + h), (x, y + h)]]
            }
            "circle" => vec![ellipse(num("cx"), num("cy"), num("r"), num("r"))],
            "ellipse" => vec![ellipse(num("cx"), num("cy"), num("rx"), num("ry"))],
            _ => Vec::new(),
        };
        for s in shapes {
            if s.len() >= 3 {
                polys.push(s.into_iter().map(|p| apply(&m, p)).collect());
            }
        }
    }
    let view = view.ok_or("SVG: no viewBox or width/height")?;
    if view[2] <= 0.0 || view[3] <= 0.0 {
        return Err("SVG: empty viewBox".into());
    }
    if polys.is_empty() {
        return Err("SVG: no filled shapes".into());
    }
    Ok((polys, view))
}

/// Even-odd coverage of `polys` on a `w × h` raster; `to_px` maps user units to pixels.
fn rasterise(polys: &[Vec<Pt>], w: usize, h: usize, to_px: impl Fn(Pt) -> Pt) -> Vec<f32> {
    let edges: Vec<(Pt, Pt)> = polys
        .iter()
        .flat_map(|p| {
            let q: Vec<Pt> = p.iter().map(|&x| to_px(x)).collect();
            (0..q.len()).map(move |i| (q[i], q[(i + 1) % q.len()])).collect::<Vec<_>>()
        })
        .filter(|(a, b)| a.1 != b.1)
        .collect();
    const SS: usize = 4;
    let mut cov = vec![0f32; w * h];
    let mut xs: Vec<f64> = Vec::new();
    for y in 0..h {
        for s in 0..SS {
            let sy = y as f64 + (s as f64 + 0.5) / SS as f64;
            xs.clear();
            for &(a, b) in &edges {
                let (lo, hi) = if a.1 < b.1 { (a, b) } else { (b, a) };
                if sy >= lo.1 && sy < hi.1 {
                    xs.push(lo.0 + (sy - lo.1) / (hi.1 - lo.1) * (hi.0 - lo.0));
                }
            }
            xs.sort_by(f64::total_cmp);
            for pair in xs.chunks_exact(2) {
                let (x0, x1) = (pair[0].clamp(0.0, w as f64), pair[1].clamp(0.0, w as f64));
                if x1 <= x0 {
                    continue;
                }
                let (i0, i1) = (x0.floor() as usize, (x1.ceil() as usize).min(w));
                for x in i0..i1 {
                    let overlap = (x1.min(x as f64 + 1.0) - x0.max(x as f64)).max(0.0);
                    cov[y * w + x] += (overlap / SS as f64) as f32;
                }
            }
        }
    }
    for c in cov.iter_mut() {
        *c = c.min(1.0);
    }
    cov
}

/// Longest side of the drawing's viewBox, user units.
pub fn view_side(text: &str) -> Option<f64> {
    parse(text).ok().map(|(_, v)| v[2].max(v[3]))
}

/// Rasterise an SVG into coverage. Returns the raster and the placement adjusted to it (for a
/// patch, `km_per_px` per viewBox unit becomes km per raster pixel).
pub fn coverage_from_svg(bytes: &[u8], place: &Placement) -> Result<(Vec<f32>, u32, u32, Placement), String> {
    let text = std::str::from_utf8(bytes).map_err(|_| "SVG: not UTF-8 text")?;
    let (polys, v) = parse(text)?;
    let (w, h, place) = match place {
        Placement::WorldEquirect => (MAX_SIDE as usize, (MAX_SIDE / 2.0) as usize, Placement::WorldEquirect),
        Placement::Patch { center_deg, km_per_px, rotation_deg } => {
            let scale = MAX_SIDE / v[2].max(v[3]); // pixels per unit
            let (w, h) = ((v[2] * scale).round().max(1.0) as usize, (v[3] * scale).round().max(1.0) as usize);
            (w, h, Placement::Patch { center_deg: *center_deg, km_per_px: km_per_px / scale, rotation_deg: *rotation_deg })
        }
    };
    let (sx, sy) = (w as f64 / v[2], h as f64 / v[3]);
    let cov = rasterise(&polys, w, h, |p| ((p.0 - v[0]) * sx, (p.1 - v[1]) * sy));
    Ok((cov, w as u32, h as u32, place))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rect_and_path_fill() {
        let svg = br#"<svg viewBox="0 0 100 50"><rect x="0" y="0" width="50" height="50"/>
            <g transform="translate(60,10)"><path d="M0 0 h20 v20 h-20 z"/></g>
            <circle cx="90" cy="40" r="5" fill="none"/></svg>"#;
        let (cov, w, h, _) = coverage_from_svg(svg, &Placement::WorldEquirect).unwrap();
        let at = |x: f64, y: f64| cov[(y * h as f64) as usize * w as usize + (x * w as f64) as usize];
        assert!(at(0.25, 0.5) > 0.99); // inside the rect
        assert!(at(0.7, 0.4) > 0.99); // inside the translated square
        assert!(at(0.55, 0.5) < 0.01); // gap
        assert!(at(0.9, 0.8) < 0.01); // unfilled circle
    }

    #[test]
    fn even_odd_hole() {
        let svg = br#"<svg width="10" height="10"><path d="M0 0 L10 0 L10 10 L0 10 Z M3 3 L7 3 L7 7 L3 7 Z"/></svg>"#;
        let (cov, w, _, _) = coverage_from_svg(svg, &Placement::Patch { center_deg: [0.0, 0.0], km_per_px: 100.0, rotation_deg: 0.0 }).unwrap();
        let w = w as usize;
        assert!(cov[w / 10 * w + w / 10] > 0.99);
        assert!(cov[(w / 2) * w + w / 2] < 0.01);
    }
}
