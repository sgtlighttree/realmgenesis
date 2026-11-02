//! Köppen-Geiger climate classification.
//!
//! Rules follow Peel, Finlayson & McMahon (2007), "Updated world map of the Köppen-Geiger
//! climate classification", Hydrology and Earth System Sciences 11, 1633-1644, Appendix A /
//! Table 1 (https://hess.copernicus.org/articles/11/1633/2007/). Colours follow the standard
//! map palette used by Peel 2007 and Beck et al. 2018 (Scientific Data 5:180214).
//!
//! Variables, for one location (monthly temperatures T and precipitations P):
//! - `MAT` / `MAP`: mean annual temperature / annual precipitation total
//! - `Thot` / `Tcold`: warmest / coldest month mean temperature
//! - `Tmon10`: number of months with mean temperature above 10 C
//! - `Pdry`: driest month overall; `Psdry` / `Pwdry`: driest month in summer / winter half-year;
//!   `Pswet` / `Pwwet`: wettest month in summer / winter half-year
//! - Summer half-year is April-September in the northern hemisphere and October-March in the
//!   southern hemisphere (winter is the other half).
//! - `Pthreshold = 2 MAT` if at least 70 % of MAP falls in the winter half-year,
//!   `2 MAT + 28` if at least 70 % falls in the summer half-year, else `2 MAT + 14`
//!   (never below zero).
//!
//! Decision order used here: Polar (E) first, then Arid (B), Tropical (A), Temperate (C) and
//! Cold (D). Peel's table does not fix an order between E and B; putting E first is what keeps
//! the ice caps and Tibet-like plateaus polar on the published maps. B must come before A
//! because hot deserts such as Dubai have `Tcold >= 18`.
//!
//! The "As" (dry summer) tropical subtype is not part of Peel's 30 classes, so it is not here:
//! such places are `Aw`.

/// Köppen-Geiger classes. `0` is reserved for "no class" (ocean).
#[repr(u8)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Koppen {
    Af = 1,
    Am,
    Aw,
    BWh,
    BWk,
    BSh,
    BSk,
    Csa,
    Csb,
    Csc,
    Cwa,
    Cwb,
    Cwc,
    Cfa,
    Cfb,
    Cfc,
    Dsa,
    Dsb,
    Dsc,
    Dsd,
    Dwa,
    Dwb,
    Dwc,
    Dwd,
    Dfa,
    Dfb,
    Dfc,
    Dfd,
    ET,
    EF,
}

impl Koppen {
    /// All 30 classes in `u8` order.
    pub const ALL: [Koppen; 30] = [
        Koppen::Af, Koppen::Am, Koppen::Aw, Koppen::BWh, Koppen::BWk, Koppen::BSh, Koppen::BSk,
        Koppen::Csa, Koppen::Csb, Koppen::Csc, Koppen::Cwa, Koppen::Cwb, Koppen::Cwc,
        Koppen::Cfa, Koppen::Cfb, Koppen::Cfc, Koppen::Dsa, Koppen::Dsb, Koppen::Dsc,
        Koppen::Dsd, Koppen::Dwa, Koppen::Dwb, Koppen::Dwc, Koppen::Dwd, Koppen::Dfa,
        Koppen::Dfb, Koppen::Dfc, Koppen::Dfd, Koppen::ET, Koppen::EF,
    ];

    /// `0` (ocean) and anything out of range give `None`.
    pub fn from_u8(v: u8) -> Option<Koppen> {
        if (1..=30).contains(&v) { Some(Self::ALL[v as usize - 1]) } else { None }
    }

    /// Classify one location from its twelve monthly mean temperatures (C) and monthly
    /// precipitation totals (mm). `northern_hemisphere` decides which half-year is "summer".
    pub fn classify(temp_c: &[f32; 12], precip_mm: &[f32; 12], northern_hemisphere: bool) -> Koppen {
        let t: [f64; 12] = std::array::from_fn(|i| temp_c[i] as f64);
        let p: [f64; 12] = std::array::from_fn(|i| precip_mm[i].max(0.0) as f64);

        let mat = t.iter().sum::<f64>() / 12.0;
        let map: f64 = p.iter().sum();
        let thot = t.iter().cloned().fold(f64::MIN, f64::max);
        let tcold = t.iter().cloned().fold(f64::MAX, f64::min);
        let tmon10 = t.iter().filter(|&&v| v > 10.0).count();

        // Summer = Apr..Sep (month indices 3..=8) in the north, Oct..Mar in the south.
        let is_summer = |m: usize| (3..=8).contains(&m) == northern_hemisphere;
        let (mut p_summer, mut p_winter) = (0.0, 0.0);
        let (mut psdry, mut pwdry) = (f64::MAX, f64::MAX);
        let (mut pswet, mut pwwet) = (0.0f64, 0.0f64);
        for (m, &pm) in p.iter().enumerate() {
            if is_summer(m) {
                p_summer += pm;
                psdry = psdry.min(pm);
                pswet = pswet.max(pm);
            } else {
                p_winter += pm;
                pwdry = pwdry.min(pm);
                pwwet = pwwet.max(pm);
            }
        }
        let pdry = psdry.min(pwdry);

        // E: polar.
        if thot < 10.0 {
            return if thot > 0.0 { Koppen::ET } else { Koppen::EF };
        }

        // B: arid.
        let pth = if p_winter >= 0.7 * map {
            2.0 * mat
        } else if p_summer >= 0.7 * map {
            2.0 * mat + 28.0
        } else {
            2.0 * mat + 14.0
        }
        .max(0.0);
        if map < 10.0 * pth {
            let desert = map < 5.0 * pth;
            let hot = mat >= 18.0;
            return match (desert, hot) {
                (true, true) => Koppen::BWh,
                (true, false) => Koppen::BWk,
                (false, true) => Koppen::BSh,
                (false, false) => Koppen::BSk,
            };
        }

        // A: tropical.
        if tcold >= 18.0 {
            return if pdry >= 60.0 {
                Koppen::Af
            } else if pdry >= 100.0 - map / 25.0 {
                Koppen::Am
            } else {
                Koppen::Aw
            };
        }

        // C and D share the second (precipitation) letter.
        let dry_summer = psdry < 40.0 && psdry < pwwet / 3.0;
        let dry_winter = pwdry < pswet / 10.0;
        let second = if dry_summer { 's' } else if dry_winter { 'w' } else { 'f' };

        if tcold > 0.0 {
            // C: temperate. Third letter a/b/c.
            let third = if thot >= 22.0 { 'a' } else if tmon10 >= 4 { 'b' } else { 'c' };
            match (second, third) {
                ('s', 'a') => Koppen::Csa,
                ('s', 'b') => Koppen::Csb,
                ('s', _) => Koppen::Csc,
                ('w', 'a') => Koppen::Cwa,
                ('w', 'b') => Koppen::Cwb,
                ('w', _) => Koppen::Cwc,
                (_, 'a') => Koppen::Cfa,
                (_, 'b') => Koppen::Cfb,
                _ => Koppen::Cfc,
            }
        } else {
            // D: cold. Third letter a/b/c/d (d = Tcold <= -38 C).
            let third = if thot >= 22.0 {
                'a'
            } else if tmon10 >= 4 {
                'b'
            } else if tcold > -38.0 {
                'c'
            } else {
                'd'
            };
            match (second, third) {
                ('s', 'a') => Koppen::Dsa,
                ('s', 'b') => Koppen::Dsb,
                ('s', 'c') => Koppen::Dsc,
                ('s', _) => Koppen::Dsd,
                ('w', 'a') => Koppen::Dwa,
                ('w', 'b') => Koppen::Dwb,
                ('w', 'c') => Koppen::Dwc,
                ('w', _) => Koppen::Dwd,
                (_, 'a') => Koppen::Dfa,
                (_, 'b') => Koppen::Dfb,
                (_, 'c') => Koppen::Dfc,
                _ => Koppen::Dfd,
            }
        }
    }

    /// Two-to-three letter Köppen code, e.g. `"Cfb"`.
    pub fn code(&self) -> &'static str {
        match self {
            Koppen::Af => "Af",
            Koppen::Am => "Am",
            Koppen::Aw => "Aw",
            Koppen::BWh => "BWh",
            Koppen::BWk => "BWk",
            Koppen::BSh => "BSh",
            Koppen::BSk => "BSk",
            Koppen::Csa => "Csa",
            Koppen::Csb => "Csb",
            Koppen::Csc => "Csc",
            Koppen::Cwa => "Cwa",
            Koppen::Cwb => "Cwb",
            Koppen::Cwc => "Cwc",
            Koppen::Cfa => "Cfa",
            Koppen::Cfb => "Cfb",
            Koppen::Cfc => "Cfc",
            Koppen::Dsa => "Dsa",
            Koppen::Dsb => "Dsb",
            Koppen::Dsc => "Dsc",
            Koppen::Dsd => "Dsd",
            Koppen::Dwa => "Dwa",
            Koppen::Dwb => "Dwb",
            Koppen::Dwc => "Dwc",
            Koppen::Dwd => "Dwd",
            Koppen::Dfa => "Dfa",
            Koppen::Dfb => "Dfb",
            Koppen::Dfc => "Dfc",
            Koppen::Dfd => "Dfd",
            Koppen::ET => "ET",
            Koppen::EF => "EF",
        }
    }

    /// Standard Köppen map colour (Peel 2007 / Beck 2018 legend), sRGB.
    pub fn color(&self) -> [u8; 3] {
        match self {
            Koppen::Af => [0, 0, 255],
            Koppen::Am => [0, 120, 255],
            Koppen::Aw => [70, 170, 250],
            Koppen::BWh => [255, 0, 0],
            Koppen::BWk => [255, 150, 150],
            Koppen::BSh => [245, 165, 0],
            Koppen::BSk => [255, 220, 100],
            Koppen::Csa => [255, 255, 0],
            Koppen::Csb => [200, 200, 0],
            Koppen::Csc => [150, 150, 0],
            Koppen::Cwa => [150, 255, 150],
            Koppen::Cwb => [100, 200, 100],
            Koppen::Cwc => [50, 150, 50],
            Koppen::Cfa => [200, 255, 80],
            Koppen::Cfb => [100, 255, 80],
            Koppen::Cfc => [50, 200, 0],
            Koppen::Dsa => [255, 0, 255],
            Koppen::Dsb => [200, 0, 200],
            Koppen::Dsc => [150, 50, 150],
            Koppen::Dsd => [150, 100, 150],
            Koppen::Dwa => [170, 175, 255],
            Koppen::Dwb => [90, 120, 220],
            Koppen::Dwc => [75, 80, 180],
            Koppen::Dwd => [50, 0, 135],
            Koppen::Dfa => [0, 255, 255],
            Koppen::Dfb => [56, 200, 255],
            Koppen::Dfc => [0, 125, 125],
            Koppen::Dfd => [0, 70, 95],
            Koppen::ET => [178, 178, 178],
            Koppen::EF => [102, 102, 102],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Station normals are rounded values from national weather services as tabulated in the
    // WMO 1991-2020 (or 1981-2010 / 1961-1990 where noted) climate normals. They are typed from
    // the published tables to about 0.1 C / 1 mm and are meant to be representative, not
    // authoritative: every case below sits well inside its class so rounding cannot flip it.

    /// Singapore (Changi), WMO 1991-2020 normals.
    const SINGAPORE_T: [f32; 12] = [26.8, 27.3, 27.9, 28.3, 28.5, 28.5, 28.2, 28.2, 27.9, 27.8, 27.1, 26.7];
    const SINGAPORE_P: [f32; 12] = [221.0, 105.0, 144.0, 160.0, 152.0, 135.0, 148.0, 137.0, 137.0, 146.0, 252.0, 313.0];

    /// Miami International, NOAA 1991-2020 normals.
    const MIAMI_T: [f32; 12] = [20.2, 21.2, 22.8, 24.6, 26.9, 28.2, 28.9, 29.0, 28.3, 26.5, 23.6, 21.4];
    const MIAMI_P: [f32; 12] = [51.0, 61.0, 69.0, 93.0, 154.0, 238.0, 168.0, 202.0, 232.0, 156.0, 83.0, 59.0];

    /// Cairo, NOAA / WMO 1961-1990 normals.
    const CAIRO_T: [f32; 12] = [14.0, 15.5, 18.3, 22.3, 26.0, 28.1, 29.0, 28.8, 27.0, 24.0, 19.8, 15.5];
    const CAIRO_P: [f32; 12] = [5.0, 4.0, 4.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 3.0, 5.0];

    /// Madrid (Retiro), AEMET 1981-2010 normals.
    const MADRID_T: [f32; 12] = [6.2, 7.9, 11.1, 13.1, 17.0, 22.2, 25.6, 25.1, 20.7, 14.8, 9.6, 6.7];
    const MADRID_P: [f32; 12] = [33.0, 34.0, 24.0, 46.0, 49.0, 25.0, 11.0, 12.0, 29.0, 59.0, 56.0, 52.0];

    /// Rome (Urbe), WMO 1961-1990 normals.
    const ROME_T: [f32; 12] = [7.7, 8.5, 11.0, 13.7, 18.0, 22.0, 24.5, 24.7, 20.7, 16.5, 11.9, 8.4];
    const ROME_P: [f32; 12] = [73.0, 70.0, 56.0, 72.0, 49.0, 32.0, 20.0, 29.0, 66.0, 116.0, 130.0, 90.0];

    /// London (Heathrow), Met Office 1991-2020 normals.
    const LONDON_T: [f32; 12] = [5.6, 5.8, 8.0, 10.5, 13.8, 17.0, 19.2, 18.9, 16.0, 12.4, 8.5, 5.9];
    const LONDON_P: [f32; 12] = [55.0, 40.0, 42.0, 44.0, 49.0, 45.0, 40.0, 48.0, 50.0, 70.0, 69.0, 59.0];

    /// Moscow, Roshydromet 1991-2020 normals.
    const MOSCOW_T: [f32; 12] = [-6.2, -5.0, 0.1, 8.1, 14.7, 18.7, 20.4, 18.5, 12.5, 6.0, -0.3, -4.6];
    const MOSCOW_P: [f32; 12] = [52.0, 45.0, 36.0, 40.0, 56.0, 80.0, 85.0, 80.0, 66.0, 69.0, 51.0, 52.0];

    /// Yakutsk, Roshydromet 1961-1990 normals.
    const YAKUTSK_T: [f32; 12] = [-38.6, -33.5, -20.5, -6.4, 5.9, 15.5, 19.6, 15.9, 7.0, -7.2, -26.0, -35.9];
    const YAKUTSK_P: [f32; 12] = [7.0, 5.0, 4.0, 8.0, 18.0, 38.0, 46.0, 44.0, 27.0, 17.0, 14.0, 8.0];

    /// Perth (Airport), Bureau of Meteorology 1991-2020 normals (southern hemisphere).
    const PERTH_T: [f32; 12] = [25.0, 25.5, 23.6, 20.5, 17.3, 14.9, 13.7, 14.1, 15.5, 17.5, 20.5, 23.0];
    const PERTH_P: [f32; 12] = [14.0, 13.0, 17.0, 41.0, 105.0, 157.0, 164.0, 133.0, 80.0, 55.0, 24.0, 15.0];

    /// Buenos Aires (Observatorio Central), SMN 1991-2020 normals (southern hemisphere).
    const BAIRES_T: [f32; 12] = [24.9, 23.6, 21.8, 17.6, 14.3, 11.2, 10.7, 12.1, 14.4, 17.9, 21.0, 23.5];
    const BAIRES_P: [f32; 12] = [130.0, 119.0, 129.0, 98.0, 87.0, 52.0, 53.0, 63.0, 69.0, 123.0, 99.0, 98.0];

    /// Darwin (Australia), BoM 1991-2020 (rounded): the classic Aw wet/dry tropics, southern.
    const DARWIN_T: [f32; 12] = [29.3, 29.0, 29.3, 29.8, 29.2, 27.8, 27.4, 28.4, 29.8, 30.4, 30.2, 29.9];
    const DARWIN_P: [f32; 12] = [430.0, 360.0, 320.0, 100.0, 20.0, 2.0, 1.0, 6.0, 15.0, 60.0, 140.0, 250.0];

    #[test]
    fn real_station_normals_classify_correctly() {
        let c = Koppen::classify;
        assert_eq!(c(&SINGAPORE_T, &SINGAPORE_P, true), Koppen::Af);
        assert_eq!(c(&MIAMI_T, &MIAMI_P, true), Koppen::Am);
        assert_eq!(c(&CAIRO_T, &CAIRO_P, true), Koppen::BWh);
        assert_eq!(c(&MADRID_T, &MADRID_P, true), Koppen::BSk);
        assert_eq!(c(&ROME_T, &ROME_P, true), Koppen::Csa);
        assert_eq!(c(&LONDON_T, &LONDON_P, true), Koppen::Cfb);
        assert_eq!(c(&MOSCOW_T, &MOSCOW_P, true), Koppen::Dfb);
        assert_eq!(c(&YAKUTSK_T, &YAKUTSK_P, true), Koppen::Dwd);
        assert_eq!(c(&PERTH_T, &PERTH_P, false), Koppen::Csa);
        assert_eq!(c(&BAIRES_T, &BAIRES_P, false), Koppen::Cfa);
        assert_eq!(c(&DARWIN_T, &DARWIN_P, false), Koppen::Aw);
    }

    #[test]
    fn southern_hemisphere_summer_is_october_to_march() {
        // Perth's dry summer is only recognised when Oct-Mar counts as summer.
        assert_eq!(Koppen::classify(&PERTH_T, &PERTH_P, false), Koppen::Csa);
        assert_ne!(Koppen::classify(&PERTH_T, &PERTH_P, true), Koppen::Csa);
        // Shifting a northern station by six months and flipping the flag must not change it.
        let rot = |a: &[f32; 12]| -> [f32; 12] { std::array::from_fn(|i| a[(i + 6) % 12]) };
        assert_eq!(
            Koppen::classify(&rot(&ROME_T), &rot(&ROME_P), false),
            Koppen::classify(&ROME_T, &ROME_P, true)
        );
    }

    #[test]
    fn polar_classes() {
        let tundra = [-20.0, -18.0, -12.0, -5.0, 1.0, 5.0, 8.0, 6.0, 2.0, -6.0, -14.0, -19.0];
        let ice = [-40.0, -38.0, -30.0, -25.0, -20.0, -15.0, -12.0, -15.0, -20.0, -28.0, -35.0, -39.0];
        let wet = [30.0; 12];
        assert_eq!(Koppen::classify(&tundra, &wet, true), Koppen::ET);
        assert_eq!(Koppen::classify(&ice, &wet, true), Koppen::EF);
    }

    #[test]
    fn hot_steppe_and_cold_desert() {
        // Hot semi-arid (Sahel-ish): MAT ~ 28, summer rain ~ 400 mm.
        let t = [24.0, 27.0, 31.0, 33.0, 32.0, 30.0, 28.0, 27.0, 28.0, 30.0, 28.0, 24.0];
        let p = [0.0, 0.0, 2.0, 8.0, 40.0, 80.0, 130.0, 150.0, 70.0, 15.0, 2.0, 0.0];
        assert_eq!(Koppen::classify(&t, &p, true), Koppen::BSh);
        // Cold desert (Gobi-ish).
        let t = [-15.0, -10.0, -1.0, 8.0, 15.0, 20.0, 23.0, 21.0, 14.0, 6.0, -4.0, -12.0];
        let p = [1.0, 1.0, 3.0, 5.0, 10.0, 20.0, 40.0, 30.0, 12.0, 4.0, 2.0, 1.0];
        assert_eq!(Koppen::classify(&t, &p, true), Koppen::BWk);
    }

    #[test]
    fn enum_roundtrips_and_has_distinct_colours() {
        for (i, k) in Koppen::ALL.iter().enumerate() {
            assert_eq!(*k as u8 as usize, i + 1);
            assert_eq!(Koppen::from_u8(*k as u8), Some(*k));
        }
        assert_eq!(Koppen::from_u8(0), None);
        assert_eq!(Koppen::from_u8(31), None);
        let mut cols: Vec<_> = Koppen::ALL.iter().map(|k| k.color()).collect();
        cols.sort();
        cols.dedup();
        assert_eq!(cols.len(), 30);
        assert_eq!(Koppen::Cfb.code(), "Cfb");
    }
}
