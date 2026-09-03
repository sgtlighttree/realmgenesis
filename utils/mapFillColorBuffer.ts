// Expand the per-cell #rrggbb colour cache (utils/mapColorCache) into a
// per-vertex RGB buffer aligned to the tessellation's positions, so the GPU
// fill surface reads one colour per triangle vertex. Parsing the shared hex
// cache keeps a single source of truth for cell colour (no second getCellColor
// path that could drift from the Canvas2D fill).
const parseHex = (hex: string): [number, number, number] => {
  const n = parseInt(hex.slice(1), 16);
  return [((n >> 16) & 255) / 255, ((n >> 8) & 255) / 255, (n & 255) / 255];
};

export const buildFillColorBuffer = (
  perCellHex: string[],
  cellTriRange: Uint32Array,
  out?: Float32Array,
): Float32Array => {
  const total = cellTriRange[cellTriRange.length - 1];
  const buf = out && out.length === total * 3 ? out : new Float32Array(total * 3);
  for (let i = 0; i < perCellHex.length; i++) {
    const [r, g, b] = parseHex(perCellHex[i]);
    for (let v = cellTriRange[i]; v < cellTriRange[i + 1]; v++) {
      buf[v * 3] = r; buf[v * 3 + 1] = g; buf[v * 3 + 2] = b;
    }
  }
  return buf;
};

export const writeCellColors = (
  buf: Float32Array,
  perCellHex: string[],
  cellTriRange: Uint32Array,
  cellIds: Iterable<number>,
): void => {
  for (const i of cellIds) {
    const [r, g, b] = parseHex(perCellHex[i]);
    for (let v = cellTriRange[i]; v < cellTriRange[i + 1]; v++) {
      buf[v * 3] = r; buf[v * 3 + 1] = g; buf[v * 3 + 2] = b;
    }
  }
};
