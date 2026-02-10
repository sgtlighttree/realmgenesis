import { BiomeType, Cell, ViewMode } from '../types';
import * as THREE from 'three';

// Earth-like Natural Colors
export const BIOME_COLORS: Record<BiomeType, string> = {
  [BiomeType.DEEP_OCEAN]: '#1a237e', // Deep Blue
  [BiomeType.OCEAN]: '#0277bd',      // Standard Blue
  
  // Polar
  [BiomeType.ICE_CAP]: '#ffffff',    // White
  [BiomeType.TUNDRA]: '#cfd8dc',     // Greyish Cyan/White
  
  // Dry
  [BiomeType.HOT_DESERT]: '#e6c27e', // Sandy Orange
  [BiomeType.COLD_DESERT]: '#bcaaa4', // Greyish Brown
  [BiomeType.STEPPE]: '#c5e1a5',     // Pale dry green
  
  // Tropical
  [BiomeType.TROPICAL_RAINFOREST]: '#004d40', // Deep Jungle Green
  [BiomeType.TROPICAL_SAVANNA]: '#aed581',    // Yellowish Green
  
  // Temperate
  [BiomeType.MEDITERRANEAN]: '#8d6e63',       // Dry brownish green
  [BiomeType.TEMPERATE_FOREST]: '#2e7d32',    // Standard Forest Green
  [BiomeType.TEMPERATE_RAINFOREST]: '#1b5e20', // Darker Green
  
  // Continental
  [BiomeType.BOREAL_FOREST]: '#00695c',       // Pine Green (Blueish)
  
  // Special
  [BiomeType.BEACH]: '#fff59d',      // Sand
  [BiomeType.VOLCANIC]: '#37474f',   // Dark Grey Rock
};

const PLATE_COLORS = [
  '#ef5350', '#ab47bc', '#7e57c2', '#5c6bc0', '#42a5f5', '#29b6f6',
  '#26c6da', '#26a69a', '#66bb6a', '#9ccc65', '#d4e157', '#ffee58',
  '#ffca28', '#ffa726', '#ff7043', '#8d6e63', '#bdbdbd', '#78909c'
];

// Helper to pseudo-randomly offset color based on index with high contrast
const getProvinceVariant = (baseColorHex: string, provId: number): THREE.Color => {
    const c = new THREE.Color(baseColorHex);
    // Deterministic random based on provId
    const r = Math.sin(provId * 12.9898) * 43758.5453;
    const rnd = r - Math.floor(r); // 0..1
    const r2 = Math.cos(provId * 78.233) * 43758.5453;
    const rnd2 = r2 - Math.floor(r2); // 0..1
    
    const hsl = { h: 0, s: 0, l: 0 };
    c.getHSL(hsl);
    
    // Significantly more contrast for provinces
    // Shift lightness (+/- 25%)
    hsl.l = Math.max(0.1, Math.min(0.9, hsl.l + (rnd * 0.5 - 0.25)));
    // Shift saturation (+/- 30%)
    hsl.s = Math.max(0.1, Math.min(1.0, hsl.s + (rnd2 * 0.6 - 0.3)));
    // Shift hue slightly (+/- 8%)
    hsl.h = (hsl.h + (rnd * 0.16 - 0.08) + 1.0) % 1.0;
    
    c.setHSL(hsl.h, hsl.s, hsl.l);
    return c;
};

export const getCellColor = (cell: Cell, mode: ViewMode, seaLevel: number): THREE.Color => {
  const color = new THREE.Color();

  switch (mode) {
    case 'satellite':
      if (cell.height < seaLevel) {
         const depth = cell.height / seaLevel;
         const deep = new THREE.Color(0x051e3e);
         const shallow = new THREE.Color(0x006994);
         color.copy(deep).lerp(shallow, Math.pow(depth, 2));
      } else {
         const t = (cell.height - seaLevel) / (1 - seaLevel);
         switch(cell.biome) {
             case BiomeType.ICE_CAP: color.setHex(0xffffff); break;
             case BiomeType.TUNDRA: color.setHex(0x78766a); break;
             case BiomeType.HOT_DESERT: color.setHex(0xdabba0); break;
             case BiomeType.COLD_DESERT: color.setHex(0x9e9587); break;
             case BiomeType.TROPICAL_RAINFOREST: color.setHex(0x052e16); break;
             case BiomeType.TEMPERATE_RAINFOREST: color.setHex(0x1a4221); break;
             case BiomeType.TEMPERATE_FOREST: color.setHex(0x285020); break; 
             case BiomeType.BOREAL_FOREST: color.setHex(0x193626); break;
             case BiomeType.TROPICAL_SAVANNA: color.setHex(0x6f7d46); break;
             case BiomeType.STEPPE: color.setHex(0x8a9263); break;
             case BiomeType.MEDITERRANEAN: color.setHex(0x6b7044); break;
             case BiomeType.BEACH: color.setHex(0xe8ddc5); break;
             case BiomeType.VOLCANIC: color.setHex(0x262626); break;
             default: color.setHex(0x335533);
         }
         let snowThreshold = 0.65;
         if (cell.temperature > 20) snowThreshold = 0.85;
         if (t > 0.35 && cell.biome !== BiomeType.ICE_CAP) {
             const rockFactor = Math.min(1, (t - 0.35) * 4);
             color.lerp(new THREE.Color(0x524e49), rockFactor);
         }
         if (t > snowThreshold) {
             const snowFactor = Math.min(1, (t - snowThreshold) * 5);
             color.lerp(new THREE.Color(0xffffff), snowFactor);
         }
      }
      break;

    case 'height':
      if (cell.height < seaLevel) {
        const t = cell.height / seaLevel;
        color.setHSL(0.6, 0.7, 0.1 + t * 0.4); 
      } else {
        const t = (cell.height - seaLevel) / (1 - seaLevel);
        if (t < 0.2) {
           color.setHSL(0.25, 0.4, 0.3 + t * 0.5);
        } else if (t < 0.6) {
           color.setHSL(0.1, 0.2, 0.4 + (t - 0.2) * 0.5);
        } else {
           color.setHSL(0, 0, 0.5 + (t - 0.6) * 1.5);
        }
      }
      break;

    case 'height_bw':
      if (cell.height < seaLevel) {
         const t = cell.height / seaLevel;
         color.setScalar(t * 0.15); 
      } else {
         const t = (cell.height - seaLevel) / (1 - seaLevel);
         const val = 0.2 + Math.pow(t, 1.5) * 0.8;
         color.setScalar(val);
      }
      break;

    case 'temperature':
      const minT = -30;
      const maxT = 50;
      const tNorm = Math.max(0, Math.min(1, (cell.temperature - minT) / (maxT - minT)));
      color.setHSL(0.65 - (tNorm * 0.65), 0.8, 0.5);
      break;

    case 'moisture':
      if (cell.height < seaLevel) {
         color.setHex(0x004488);
      } else {
         color.setHSL(0.6, cell.moisture, 0.9 - (cell.moisture * 0.5));
      }
      break;

    case 'plates':
      color.set(PLATE_COLORS[cell.plateId % PLATE_COLORS.length]);
      if (cell.height < seaLevel) color.multiplyScalar(0.7);
      break;
      
    case 'political':
       if (cell.regionId !== undefined) {
          const baseColor = PLATE_COLORS[cell.regionId % PLATE_COLORS.length];
          if (cell.provinceId !== undefined) {
              color.copy(getProvinceVariant(baseColor, cell.provinceId));
          } else {
              color.set(baseColor);
          }
          // Fix for territorial waters: Blend faction color with deep blue
          if (cell.height < seaLevel) {
              color.lerp(new THREE.Color(0x1a237e), 0.65);
          }
       } else {
          // Unclaimed territory
          if (cell.height < seaLevel) {
              color.setHex(0x1a237e); 
              color.multiplyScalar(0.5 + cell.height * 0.5);
          } else {
              color.setHex(0x555555); 
          }
       }
       break;

    case 'biome':
    default:
      color.set(BIOME_COLORS[cell.biome] || '#ff00ff');
      break;
  }

  return color;
};