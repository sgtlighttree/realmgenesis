import { describe, expect, it } from 'vitest';

import { RNG } from '../utils/rng';
import { chordDistance, generateFibonacciSphere } from '../utils/spherical';
import { MacroPointGrid } from '../utils/tectonicsV3';
import { Point } from '../types';

// `MacroPointGrid` replaced a brute-force argmin in `projectTectonicsToDisplay`.
// Its whole justification is that it returns EXACTLY what that scan returned —
// a faster nearest-neighbour that disagreed even once would silently change
// generated worlds for every existing seed. So the test is not "is it close":
// it is an index-for-index comparison against the original loop.
const bruteNearest = (q: Point, pts: Point[]): number => {
  let nearest = 0;
  let minDist = Infinity;
  for (let j = 0; j < pts.length; j++) {
    const d = chordDistance(q, pts[j]);
    if (d < minDist) { minDist = d; nearest = j; }
  }
  return nearest;
};

const agree = (macroN: number, queryN: number, jitter: number) => {
  const macro = generateFibonacciSphere(macroN, new RNG('macro'), jitter);
  const queries = generateFibonacciSphere(queryN, new RNG('query'), jitter);
  const grid = new MacroPointGrid(macro);
  let mismatches = 0;
  for (const q of queries) {
    if (grid.nearestIndex(q.x, q.y, q.z) !== bruteNearest(q, macro)) mismatches++;
  }
  return mismatches;
};

describe('MacroPointGrid', () => {
  it('matches the brute-force argmin index for index', () => {
    expect(agree(1500, 4000, 0.4)).toBe(0);
  });

  // jitter 0 is a perfectly regular Fibonacci lattice — the maximum-tie case,
  // and the one where a wrong tie-break rule shows up. The original scan used a
  // strict `<`, so on an exact tie the LOWEST index wins; the grid visits
  // buckets in a different order and has to reproduce that deliberately.
  it('breaks ties toward the lower index, as the original strict < did', () => {
    expect(agree(1500, 4000, 0)).toBe(0);
  });

  // The ring walk stops once every unvisited bucket is provably farther than
  // the best distance so far. A too-eager bound is invisible at one density and
  // wrong at another, so the ratio of query points to macro points is varied.
  it('holds when queries greatly outnumber macro points', () => {
    expect(agree(400, 6000, 0.4)).toBe(0);
  });

  it('holds when macro points greatly outnumber queries', () => {
    expect(agree(6000, 400, 0.4)).toBe(0);
  });
});
