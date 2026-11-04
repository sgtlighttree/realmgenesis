//! What the colours of a display layer mean, as plain data (no drawing). The app turns this
//! into a legend panel. Colours come from the same functions the bake uses.

use crate::layers::{self, Layer};
use crate::palette::{self, Rgb};
use rg_climate::Koppen;

#[derive(Clone, Debug)]
pub enum Legend {
    /// A continuous colour ramp, drawn as a vertical bar (`max` at the top).
    Ramp { title: &'static str, unit: &'static str, min: f32, max: f32, ticks: Vec<f32>, colour_at: fn(f32) -> Rgb },
    /// Discrete classes: a swatch and a label each.
    Swatches { title: &'static str, items: Vec<(Rgb, String)> },
    /// A layer without a meaningful scale.
    Caption(&'static str),
}

pub fn legend(layer: Layer) -> Legend {
    match layer {
        Layer::Elevation => Legend::Ramp {
            title: "Elevation",
            unit: "m",
            min: -9000.0,
            max: 5500.0,
            ticks: vec![-9000.0, -6000.0, -3000.0, 0.0, 1600.0, 2800.0, 4200.0, 5500.0],
            colour_at: palette::elevation,
        },
        Layer::TempAnnual => temperature("Temperature · year"),
        Layer::TempJanuary => temperature("Temperature · January"),
        Layer::TempJuly => temperature("Temperature · July"),
        Layer::PrecipAnnual => Legend::Ramp {
            title: "Precipitation · year",
            unit: "mm",
            min: 0.0,
            max: 6000.0,
            ticks: vec![0.0, 600.0, 1200.0, 2000.0, 3500.0, 6000.0],
            colour_at: palette::precipitation,
        },
        Layer::Koppen => {
            let mut items: Vec<(Rgb, String)> =
                Koppen::ALL.iter().map(|k| (k.color(), format!("{} {}", k.code(), koppen_name(*k)))).collect();
            items.push((layers::KOPPEN_OCEAN, "Ocean".to_string()));
            items.push((layers::KOPPEN_SEA_ICE, "Sea ice".to_string()));
            Legend::Swatches { title: "Köppen climate", items }
        }
        Layer::Satellite => {
            let land = |k: Koppen, name: &str| (layers::satellite_base(k), name.to_string());
            let items = vec![
                land(Koppen::Af, "Rainforest"),
                land(Koppen::Aw, "Savanna"),
                land(Koppen::BWh, "Hot desert"),
                land(Koppen::BWk, "Cold desert"),
                land(Koppen::BSh, "Steppe"),
                land(Koppen::Csa, "Mediterranean scrub"),
                land(Koppen::Cfb, "Temperate forest"),
                land(Koppen::Dfb, "Cool forest"),
                land(Koppen::Dfc, "Boreal forest"),
                land(Koppen::ET, "Tundra"),
                land(Koppen::EF, "Ice cap"),
                ([234, 238, 244], "Snow".to_string()),
                ([120, 112, 106], "Bare rock".to_string()),
                (palette::sat_ocean(-40.0), "Shallow sea".to_string()),
                (palette::sat_ocean(-4000.0), "Deep sea".to_string()),
            ];
            Legend::Swatches { title: "Satellite", items }
        }
        Layer::Plates => {
            Legend::Caption("Plate colours are arbitrary. Boundaries: red convergent, blue divergent, green transform.")
        }
    }
}

fn temperature(title: &'static str) -> Legend {
    Legend::Ramp {
        title,
        unit: "°C",
        min: -45.0,
        max: 40.0,
        ticks: vec![-40.0, -30.0, -20.0, -10.0, 0.0, 10.0, 20.0, 30.0, 40.0],
        colour_at: palette::temperature,
    }
}

/// Short plain-English name for each Köppen class.
fn koppen_name(k: Koppen) -> &'static str {
    use Koppen::*;
    match k {
        Af => "tropical rainforest",
        Am => "tropical monsoon",
        Aw => "tropical savanna",
        BWh => "hot desert",
        BWk => "cold desert",
        BSh => "hot steppe",
        BSk => "cold steppe",
        Csa => "hot-summer Mediterranean",
        Csb => "warm-summer Mediterranean",
        Csc => "cold-summer Mediterranean",
        Cwa => "humid subtropical, dry winter",
        Cwb => "subtropical highland",
        Cwc => "cold subtropical highland",
        Cfa => "humid subtropical",
        Cfb => "oceanic",
        Cfc => "subpolar oceanic",
        Dsa => "hot dry-summer continental",
        Dsb => "warm dry-summer continental",
        Dsc => "dry-summer subarctic",
        Dsd => "extreme dry-summer subarctic",
        Dwa => "hot continental, dry winter",
        Dwb => "warm continental, dry winter",
        Dwc => "subarctic, dry winter",
        Dwd => "extreme subarctic, dry winter",
        Dfa => "hot-summer continental",
        Dfb => "warm-summer continental",
        Dfc => "subarctic",
        Dfd => "extreme subarctic",
        ET => "tundra",
        EF => "ice cap",
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_layer_has_a_legend() {
        for l in Layer::ALL {
            match legend(l) {
                Legend::Ramp { min, max, ticks, .. } => {
                    assert!(min < max);
                    assert!(ticks.iter().all(|t| (min..=max).contains(t)));
                }
                Legend::Swatches { items, .. } => assert!(!items.is_empty()),
                Legend::Caption(c) => assert!(!c.is_empty()),
            }
        }
    }
}
