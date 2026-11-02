//! Colour ramps. Colours are sRGB bytes.

pub type Rgb = [u8; 3];

/// Piecewise-linear ramp over (value, colour) stops, sorted by value.
pub fn ramp(stops: &[(f32, Rgb)], v: f32) -> Rgb {
    if v <= stops[0].0 {
        return stops[0].1;
    }
    for w in stops.windows(2) {
        let (a, ca) = w[0];
        let (b, cb) = w[1];
        if v <= b {
            let t = (v - a) / (b - a);
            return mix(ca, cb, t);
        }
    }
    stops[stops.len() - 1].1
}

pub fn mix(a: Rgb, b: Rgb, t: f32) -> Rgb {
    let t = t.clamp(0.0, 1.0);
    [
        (a[0] as f32 + (b[0] as f32 - a[0] as f32) * t).round() as u8,
        (a[1] as f32 + (b[1] as f32 - a[1] as f32) * t).round() as u8,
        (a[2] as f32 + (b[2] as f32 - a[2] as f32) * t).round() as u8,
    ]
}

/// Elevation (hypsometric tints, cartographic convention).
pub fn elevation(h: f32) -> Rgb {
    const STOPS: &[(f32, Rgb)] = &[
        (-9000.0, [8, 20, 52]),
        (-5000.0, [18, 46, 96]),
        (-3000.0, [30, 74, 132]),
        (-1000.0, [52, 112, 170]),
        (-150.0, [96, 158, 204]),
        (-1.0, [150, 196, 226]),
        (0.0, [92, 140, 84]),
        (300.0, [138, 170, 100]),
        (800.0, [200, 196, 128]),
        (1600.0, [190, 150, 98]),
        (2800.0, [150, 108, 80]),
        (4200.0, [170, 160, 156]),
        (5500.0, [245, 245, 250]),
    ];
    ramp(STOPS, h)
}

/// Temperature, °C: diverging blue–white–red.
pub fn temperature(t: f32) -> Rgb {
    const STOPS: &[(f32, Rgb)] = &[
        (-45.0, [40, 20, 90]),
        (-25.0, [52, 84, 180]),
        (-10.0, [110, 166, 222]),
        (0.0, [232, 240, 246]),
        (10.0, [250, 222, 140]),
        (20.0, [240, 150, 70]),
        (30.0, [200, 60, 44]),
        (40.0, [110, 16, 30]),
    ];
    ramp(STOPS, t)
}

/// Annual precipitation, mm: dry browns → greens → blues.
pub fn precipitation(p: f32) -> Rgb {
    const STOPS: &[(f32, Rgb)] = &[
        (0.0, [150, 110, 70]),
        (250.0, [214, 186, 120]),
        (600.0, [196, 214, 130]),
        (1200.0, [96, 172, 96]),
        (2000.0, [40, 130, 120]),
        (3500.0, [30, 80, 160]),
        (6000.0, [40, 30, 120]),
    ];
    ramp(STOPS, p)
}

/// Ocean colour for the satellite layer by depth (m, negative).
pub fn sat_ocean(h: f32) -> Rgb {
    const STOPS: &[(f32, Rgb)] = &[
        (-7000.0, [10, 26, 62]),
        (-4000.0, [16, 40, 84]),
        (-1500.0, [22, 60, 108]),
        (-200.0, [32, 96, 136]),
        (-40.0, [52, 136, 156]),
        (0.0, [78, 160, 166]),
    ];
    ramp(STOPS, h)
}

/// Plate colours: well-separated muted hues from an id.
pub fn plate(id: u16) -> Rgb {
    let h = (id as f32 * 0.618_034).fract();
    hsv(h, 0.45, 0.78)
}

pub fn hsv(h: f32, s: f32, v: f32) -> Rgb {
    let i = (h * 6.0).floor();
    let f = h * 6.0 - i;
    let (p, q, t) = (v * (1.0 - s), v * (1.0 - f * s), v * (1.0 - (1.0 - f) * s));
    let (r, g, b) = match i as i32 % 6 {
        0 => (v, t, p),
        1 => (q, v, p),
        2 => (p, v, t),
        3 => (p, q, v),
        4 => (t, p, v),
        _ => (v, p, q),
    };
    [(r * 255.0) as u8, (g * 255.0) as u8, (b * 255.0) as u8]
}
