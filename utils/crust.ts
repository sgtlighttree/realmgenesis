import { Point, WorldParams } from '../types';
import { RNG, SimplexNoise } from './rng';

export interface CrustField {
  crustTypes: Uint8Array; // 0=oceanic, 1=continental
  crustThickness: Float32Array; // normalized 0-1
}

// Seed the crust field independently of plates.
// landStyle maps to a noise threshold that controls continental coverage:
//   Continents: 0.45 (default)
//   Pangea:     0.6  (large connected landmass)
//   Archipelago: 0.25 (scattered land)
//   Islands:    0.15 (tiny islands)
export function seedCrustField(
  points: Point[],
  params: WorldParams,
  simplex: SimplexNoise,
  crustRng: RNG,
): CrustField {
  const n = points.length;
  const crustTypes = new Uint8Array(n);
  const crustThickness = new Float32Array(n);

  const freq = 0.3;

  let threshold = 0.45;
  if (params.landStyle === 'Pangea') threshold = 0.6;
  else if (params.landStyle === 'Archipelago') threshold = 0.25;
  else if (params.landStyle === 'Islands') threshold = 0.15;

  const jitterAmp = 0.08;

  for (let i = 0; i < n; i++) {
    const p = points[i];
    const noise = simplex.noise3D(p.x * freq, p.y * freq, p.z * freq);
    const jitter = (crustRng.next() - 0.5) * jitterAmp;
    const isContinental = noise + jitter > threshold;

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