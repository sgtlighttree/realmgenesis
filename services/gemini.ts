import type { GoogleGenAI } from "@google/genai";
import { WorldData, LoreData } from "../types";

let ai: GoogleGenAI | null = null;
let runtimeKey: string | null = null;

export const setRuntimeApiKey = (key: string) => {
    runtimeKey = key;
    ai = null; // Reset instance to use new key
};

// The Gemini client is dynamically imported so ~200 kB of SDK stays out of
// the main bundle until the user actually generates lore.
const getAI = async (): Promise<GoogleGenAI | null> => {
    if (ai) return ai;
    const key = runtimeKey || process.env.GEMINI_API_KEY || '';
    if (!key) return null;
    const { GoogleGenAI } = await import("@google/genai");
    ai = new GoogleGenAI({ apiKey: key });
    return ai;
};

export const generateWorldLore = async (world: WorldData): Promise<LoreData> => {
  const aiInstance = await getAI();
  if (!aiInstance) {
    throw new Error("No Gemini API key configured — add one in the Sys tab (or set GEMINI_API_KEY at build time).");
  }

  const level = world.params.loreLevel;
  const civs = world.civData;
  if (!civs) throw new Error("No civilization data — generate a world first.");

  // Prepare minimal context to save tokens
  const factionSummaries = civs.factions.map(f => {
    const c = world.cells[f.capitalId];
    return {
        id: f.id,
        pop: f.totalPopulation,
        biome: c.biome,
        provinces: f.provinces.length
    };
  });

  let prompt = `
    You are a fantasy world builder. 
    World Context: 
    - ${world.params.plates} tectonic plates.
    - Dominant biome: ${getDominantBiome(world)}.
    
    Factions: ${JSON.stringify(factionSummaries)}
  `;

  if (level === 1) {
      prompt += `
      Task: Generate a Name and Description for the world, and Names for the ${factionSummaries.length} factions and their Capitals.
      Output JSON structure:
      {
        "worldName": "string",
        "description": "string",
        "factions": [ { "id": number, "name": "string", "capitalName": "string" } ]
      }
      `;
  } else if (level === 2) {
      prompt += `
      Task: Generate names for the world, factions, capitals, and ALL provinces/towns.
      Output JSON structure:
      {
        "worldName": "string",
        "description": "string",
        "factions": [ 
           { 
             "id": number, "name": "string", "capitalName": "string",
             "provinceNames": ["string", "string"...] // One for each province count
           } 
        ]
      }
      `;
  } else {
      prompt += `
      Task: Generate deep lore. World name, description, faction names, capitals, province names, AND a short backstory (50 words) for each faction.
      Output JSON structure:
      {
        "worldName": "string",
        "description": "string",
        "factions": [ 
           { 
             "id": number, "name": "string", "capitalName": "string",
             "description": "string",
             "provinceNames": ["string"...]
           } 
        ]
      }
      `;
  }

  try {
    const response = await aiInstance.models.generateContent({
      model: 'gemini-3-flash-preview',
      contents: prompt,
      config: {
        responseMimeType: "application/json",
      }
    });

    const text = response.text;
    if (!text) throw new Error("Gemini returned an empty response.");
    const json = JSON.parse(text) as {
        worldName?: unknown;
        description?: unknown;
        factions?: unknown;
    };

    // Apply names back to WorldData (mutating the object in memory).
    // The model's output is untrusted: validate every field before applying
    // so a malformed response never corrupts civData.
    if (Array.isArray(json.factions) && world.civData) {
        json.factions.forEach((fJson: { id?: unknown, name?: unknown, description?: unknown, capitalName?: unknown, provinceNames?: unknown }) => {
            if (typeof fJson?.id !== 'number') return;
            const fData = world.civData!.factions.find(f => f.id === fJson.id);
            if (!fData) return;

            if (typeof fJson.name === 'string' && fJson.name) fData.name = fJson.name;
            if (typeof fJson.description === 'string') fData.description = fJson.description;

            // Apply Capital Name (province 0 is the capital province in our gen)
            if (typeof fJson.capitalName === 'string' && fJson.capitalName &&
                fData.provinces.length > 0 && fData.provinces[0].towns.length > 0) {
                const capTown = fData.provinces[0].towns.find(t => t.isCapital);
                if (capTown) capTown.name = fJson.capitalName;
            }

            // Apply Province Names
            if (Array.isArray(fJson.provinceNames)) {
                const provinceNames = fJson.provinceNames;
                fData.provinces.forEach((p, idx) => {
                    const pName = provinceNames[idx];
                    if (typeof pName === 'string' && pName) {
                        p.name = pName;
                        // Name the main town same as province for simplicity if Level 2
                        if (p.towns.length > 0 && !p.towns[0].isCapital) {
                            p.towns[0].name = p.name + " City";
                        }
                    }
                });
            }
        });
    }

    return {
        name: typeof json.worldName === 'string' && json.worldName ? json.worldName : "Unnamed",
        description: typeof json.description === 'string' && json.description ? json.description : "No description."
    };

  } catch (error) {
    // Surface the failure to the caller (App shows it in the console log);
    // never return sentinel lore that looks like a successful result.
    console.error("Gemini Error:", error);
    throw error instanceof Error ? error : new Error(String(error));
  }
};

function getDominantBiome(world: WorldData): string {
  const counts: Record<string, number> = {};
  world.cells.forEach(c => counts[c.biome] = (counts[c.biome] || 0) + 1);
  return Object.entries(counts).sort((a,b) => b[1]-a[1])[0]?.[0] || 'Unknown';
}