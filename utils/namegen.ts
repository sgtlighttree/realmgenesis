// Offline, seeded, char-level Markov name generator (feature A2).
// No network, no dependencies: every wordlist is embedded below. Given the
// same style and the same rng stream, output is fully deterministic, so a
// world's names are reproducible from its civSeed.

export type NameStyle = 'fantasy' | 'norse' | 'latin' | 'desert';

export const NAME_STYLES: NameStyle[] = ['fantasy', 'norse', 'latin', 'desert'];

export interface NameGenerator {
  faction(): string;
  province(): string;
  town(): string;
}

// Order-2 model boundaries. '^' pads the start of every context, '$' marks the
// end of a word so generation knows when a word may terminate.
const START = '^';
const END = '$';
const ORDER = 2;

const MIN_LEN = 4;
const MAX_LEN = 12;
const MAX_ATTEMPTS = 40;
const MAX_DEDUP_RETRIES = 8;

// Public-domain-flavored source words, hand-written per style. They are only
// training fodder for the Markov chain — generated names are novel blends, not
// verbatim picks.
const WORDLISTS: Record<NameStyle, string[]> = {
  fantasy: [
    'eldoria', 'valmoor', 'brighthold', 'ravenwood', 'silverpine', 'ashenvale',
    'greymark', 'thornfield', 'oakenshire', 'westmere', 'dawnhaven', 'stormgate',
    'frostmere', 'wyndham', 'lorendale', 'mistwood', 'highmoor', 'ironhold',
    'duskwatch', 'emberfall', 'goldcrest', 'starfall', 'northreach', 'redcliff',
    'shadowfen', 'whitewick', 'briarmoor', 'fallowmere', 'kingsreach', 'moorland',
    'ravenhall', 'silvermere', 'thistledown', 'wolfden', 'ambervale', 'blackwater',
    'cinderhold', 'dellmoor', 'elmsworth', 'faircrest', 'grimwald', 'hollowreach',
    'lanmoor', 'mournhold', 'oldbarrow', 'pinehaven', 'quellmoor', 'reedwick',
    'stonewatch', 'timberfell', 'valewood', 'winterhold', 'yarrowdale', 'brackenholt',
    'crowmoor', 'dunmere', 'everfrost', 'foxglen', 'glimmerbrook', 'harrowgate',
    'larkspire', 'meadowmere', 'nightvale', 'ravenmoor', 'sablewood', 'tanglewood',
    'underhill', 'wildemere', 'ashford', 'belmoor', 'cragmaw', 'dovewick',
  ],
  norse: [
    'bjornstad', 'valholl', 'skardvik', 'frosthelm', 'ravnfjord', 'ulfheim',
    'gunnarholm', 'thorsby', 'eiravik', 'norheim', 'skaldby', 'vestrfjord',
    'draugrborg', 'hrafnstad', 'iskeld', 'jotunheim', 'kaldrborg', 'lindholm',
    'myrkfjell', 'nordvik', 'ormstad', 'raudberg', 'sigrholm', 'trollstad',
    'vargheim', 'ysheim', 'asgardby', 'bergvik', 'dvalinheim', 'fenrisvik',
    'geirholm', 'haldorstad', 'iserfjord', 'kolbjornstad', 'leifsby', 'midgardholm',
    'njordvik', 'ostberg', 'ragnarholm', 'skjoldby', 'tyrvangr', 'vidarheim',
    'yggheim', 'brynjolf', 'egilstad', 'freyholm', 'grimsby', 'holmgard',
    'ivarheim', 'ketilstad', 'magnusfjord', 'oddvik', 'runeberg', 'steinvik',
    'thrymheim', 'valdrfjord', 'aldheim', 'bragholm', 'daegrborg', 'faldervik',
    'gjallarholm', 'helvegr', 'isgardby', 'kjartansby', 'lokstad', 'nornheim',
    'skirvik', 'tostigby', 'vanheim', 'warngard',
  ],
  latin: [
    'aurelia', 'brixia', 'castellum', 'draconis', 'emporium', 'fortunia',
    'graccia', 'hadria', 'imperia', 'juliana', 'lavinia', 'marcellum',
    'novaria', 'octavia', 'pontia', 'quirinia', 'romulia', 'severia',
    'tarquinia', 'umbria', 'valentia', 'aquilonia', 'belloria', 'caelia',
    'domitia', 'faventia', 'gallionis', 'liburnia', 'messania', 'narbonis',
    'ostia', 'placentia', 'ravennia', 'salernia', 'tibernum', 'venetia',
    'arretium', 'brundisium', 'clusium', 'firmium', 'genavium', 'herculia',
    'latium', 'mediolanum', 'neapolis', 'praeneste', 'reatinum', 'sutrium',
    'tarentum', 'veronia', 'antiquia', 'basilica', 'cornelia', 'flaminia',
    'lucania', 'maritima', 'aemilia', 'campania', 'etruria', 'ligusta',
    'sabinia', 'volscia', 'aternum', 'cremonia', 'histria', 'padua',
    'salvia', 'tuscania', 'urbinum', 'vibonia',
  ],
  desert: [
    'qasrabad', 'zahran', 'marrakah', 'sabkha', 'kharoum', 'jaziran',
    'nafudi', 'oryxa', 'sahelim', 'tademah', 'wadirum', 'yathrib',
    'zamora', 'aksara', 'basirah', 'dahnan', 'faidah', 'ghadar',
    'halwan', 'irbani', 'jerada', 'kufran', 'lahmar', 'madain',
    'najran', 'oualata', 'qatif', 'rasalem', 'sinawa', 'timbuk',
    'ubarath', 'zafran', 'amrikha', 'bahriya', 'daraqim', 'esharra',
    'fayoumi', 'gharbia', 'hejaza', 'iskandar', 'karbala', 'lakhmid',
    'mahdiya', 'nizwara', 'qusayr', 'ramlah', 'suakim', 'tabuki',
    'wahati', 'zaghwan', 'adrarim', 'bilmah', 'dongola', 'ghawar',
    'hofufa', 'juba', 'khamsa', 'meroe', 'nubara', 'oasima',
    'raziq', 'shamal', 'taghit', 'yalbis', 'zerzura', 'ainsefra',
    'darbandi', 'hamada', 'siwa', 'tozeur',
  ],
};

type Model = Map<string, string[]>;

// Lazily-built order-2 tables per style, cached across generator instances.
const MODEL_CACHE = new Map<NameStyle, Model>();

const trainModel = (style: NameStyle): Model => {
  const cached = MODEL_CACHE.get(style);
  if (cached) return cached;

  const model: Model = new Map();
  for (const raw of WORDLISTS[style]) {
    const word = raw.toLowerCase();
    const padded = START.repeat(ORDER) + word + END;
    for (let i = 0; i <= padded.length - ORDER - 1; i++) {
      const context = padded.slice(i, i + ORDER);
      const next = padded[i + ORDER];
      const bucket = model.get(context);
      if (bucket) bucket.push(next);
      else model.set(context, [next]);
    }
  }
  MODEL_CACHE.set(style, model);
  return model;
};

const capitalize = (word: string): string =>
  word.charAt(0).toUpperCase() + word.slice(1);

// One raw generation pass. Returns a lowercase word or '' if it ran off the
// end without terminating cleanly within MAX_LEN.
const generateRaw = (model: Model, rng: () => number): string => {
  let context = START.repeat(ORDER);
  let out = '';
  while (out.length <= MAX_LEN) {
    const bucket = model.get(context);
    if (!bucket || bucket.length === 0) return '';
    const next = bucket[Math.floor(rng() * bucket.length)];
    if (next === END) return out;
    out += next;
    context = (context + next).slice(-ORDER);
  }
  return '';
};

export function createNameGenerator(style: NameStyle, rng: () => number): NameGenerator {
  const model = trainModel(style);
  const used = new Set<string>();

  const makeName = (): string => {
    for (let retry = 0; retry < MAX_DEDUP_RETRIES; retry++) {
      let candidate = '';
      for (let attempt = 0; attempt < MAX_ATTEMPTS; attempt++) {
        const raw = generateRaw(model, rng);
        if (raw.length >= MIN_LEN && raw.length <= MAX_LEN) {
          candidate = raw;
          break;
        }
        // Keep a too-short fallback around in case we never hit the band.
        if (raw.length >= 3 && !candidate) candidate = raw;
      }
      // Absolute fallback: deterministic pad so we never emit an empty name.
      if (!candidate) candidate = 'val' + Math.floor(rng() * 1000);
      const name = capitalize(candidate);
      if (!used.has(name)) {
        used.add(name);
        return name;
      }
    }
    // Collisions exhausted: append a deterministic disambiguator.
    let name = capitalize(generateRaw(model, rng) || 'vale');
    let suffix = 2;
    while (used.has(name)) {
      name = `${capitalize(generateRaw(model, rng) || 'vale')}${suffix}`;
      suffix++;
    }
    used.add(name);
    return name;
  };

  return {
    faction: makeName,
    province: makeName,
    town: makeName,
  };
}
