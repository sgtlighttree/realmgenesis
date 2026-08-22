import { Point, WorldParams, LandStyle } from '../types';
import { RNG, SimplexNoise } from './rng';

export interface CrustField {
  crustTypes: Uint8Array; // 0=oceanic, 1=continental
  crustThickness: Float32Array; // normalized 0-1
}

// Fractal (fBm) sample of the crust noise field at a sphere point.
// Local copy so this module has no import cycle with worldGen.ts (which is
// worldGen -> tectonicsV3 -> crust).
function fbm(
  s: SimplexNoise,
  x: number,
  y: number,
  z: number,
  octaves: number,
  persistence: number,
  lacunarity: number,
): number {
  let total = 0;
  let freq = 1;
  let amp = 1;
  let max = 0;
  for (let o = 0; o < octaves; o++) {
    total += s.noise3D(x * freq, y * freq, z * freq) * amp;
    max += amp;
    amp *= persistence;
    freq *= lacunarity;
  }
  return total / max;
}

// Per-landStyle crust-field parameters.
//
// WHY THIS EXISTS (D9, 2026-08-22): the old field sampled a SINGLE octave at
// base frequency 0.3 on the unit sphere. At that scale the whole planet sits
// inside one noise lobe, so thresholding it split the sphere into one land cap
// and one ocean cap — every default world came out as a pangea (measured: one
// connected component holding 100% of land, clump metric 0.74). The base
// frequency is now near 1.0 and the field is fractal, so continental crust
// breaks into several separate masses at continental scale.
//
// Measured component structure (10k macro points, mean of 5 seeds):
//   Continents: ~7 masses, largest ~47% of land, land ~35%   (Earth-like)
//   Pangea:     1 dominant mass ~86% + a satellite, land ~31% (intentional)
//   Archipelago: ~21 masses, largest ~21%, land ~35%
//   Islands:    ~27 masses, largest ~10%, land ~34%
//
// Higher `frequency` = smaller, more numerous masses. `threshold` sets land
// coverage. Diagnosis + tuning: tests/crustDistribution.test.ts.
interface CrustStyle {
  frequency: number;
  octaves: number;
  persistence: number;
  lacunarity: number;
  threshold: number;
}

const CRUST_STYLES: Record<Exclude<LandStyle, 'Custom'>, CrustStyle> = {
  Continents: { frequency: 1.0, octaves: 3, persistence: 0.5, lacunarity: 2.0, threshold: 0.1 },
  Pangea: { frequency: 0.4, octaves: 2, persistence: 0.4, lacunarity: 2.0, threshold: 0.22 },
  Archipelago: { frequency: 2.2, octaves: 4, persistence: 0.5, lacunarity: 2.0, threshold: 0.1 },
  Islands: { frequency: 3.5, octaves: 5, persistence: 0.5, lacunarity: 2.0, threshold: 0.12 },
};

// Seed the crust field independently of plates.
export function seedCrustField(
  points: Point[],
  params: WorldParams,
  simplex: SimplexNoise,
  crustRng: RNG,
): CrustField {
  const n = points.length;
  const crustTypes = new Uint8Array(n);
  const crustThickness = new Float32Array(n);

  // 'Custom' (any hand-tuned world) uses the Continents field.
  const style =
    params.landStyle && params.landStyle !== 'Custom'
      ? CRUST_STYLES[params.landStyle]
      : CRUST_STYLES.Continents;

  const jitterAmp = 0.08;

  for (let i = 0; i < n; i++) {
    const p = points[i];
    const f = style.frequency;
    const noise = fbm(simplex, p.x * f, p.y * f, p.z * f, style.octaves, style.persistence, style.lacunarity);
    const jitter = (crustRng.next() - 0.5) * jitterAmp;
    const isContinental = noise + jitter > style.threshold;

    crustTypes[i] = isContinental ? 1 : 0;
    const thicknessNoise = simplex.noise3D(p.x * 0.5, p.y * 0.5, p.z * 0.5);
    if (isContinental) {
      crustThickness[i] = 0.6 + (thicknessNoise * 0.5 + 0.5) * 0.4;
    } else {
      crustThickness[i] = Math.max(0, (thicknessNoise * 0.5 + 0.5) * 0.3);
    }
  }

  return { crustTypes, crustThickness };
}
