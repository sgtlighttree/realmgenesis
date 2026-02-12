import { GoogleGenAI, Type } from "@google/genai";
import { WorldData, LoreData, BiomeType, CivData } from "../types";

let ai: GoogleGenAI | null = null;
let runtimeKey: string | null = null;

export const setRuntimeApiKey = (key: string) => {
    runtimeKey = key;
    ai = null; // Reset instance to use new key
};

const getAI = () => {
    if (ai) return ai;
    const key = runtimeKey || process.env.API_KEY || '';
    if (!key) return null;
    ai = new GoogleGenAI({ apiKey: key });
    return ai;
};

export const generateWorldLore = async (world: WorldData): Promise<LoreData> => {
  const aiInstance = getAI();
  if (!aiInstance) {
    return {
      name: "Lore Disabled",
      description: "Please provide a Gemini API Key in the System settings to enable AI features.",
    };
  }
  
  const level = world.params.loreLevel || 1;
  const civs = world.civData;
  if (!civs) return { name: "Wilderness", description: "No civilization data." };

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
    const aiInstance = getAI();
    if (!aiInstance) throw new Error("GoogleGenAI not initialized. API Key likely missing.");

    const response = await aiInstance.models.generateContent({
      model: 'gemini-3-flash-preview',
      contents: prompt,
      config: {
        responseMimeType: "application/json",
      }
    });

    const text = response.text;
    if (!text) throw new Error("No response text");
    const json = JSON.parse(text);

    // Apply names back to WorldData (Mutating the object in memory)
    if (json.factions && world.civData) {
        json.factions.forEach((fJson: any) => {
            const fData = world.civData!.factions.find(f => f.id === fJson.id);
            if (fData) {
                fData.name = fJson.name;
                fData.description = fJson.description;
                // Apply Capital Name
                if (fData.provinces.length > 0 && fData.provinces[0].towns.length > 0) {
                     // Usually province 0 is capital province in our gen
                     const capTown = fData.provinces[0].towns.find(t => t.isCapital);
                     if (capTown) capTown.name = fJson.capitalName;
                }
                
                // Apply Province Names
                if (fJson.provinceNames && Array.isArray(fJson.provinceNames)) {
                    fData.provinces.forEach((p, idx) => {
                        if (fJson.provinceNames[idx]) {
                            p.name = fJson.provinceNames[idx];
                            // Name the main town same as province for simplicity if Level 2
                            if (p.towns.length > 0 && !p.towns[0].isCapital) {
                                p.towns[0].name = p.name + " City";
                            }
                        }
                    });
                }
            }
        });
    }

    return {
        name: json.worldName || "Unnamed",
        description: json.description || "No description."
    };

  } catch (error) {
    console.error("Gemini Error:", error);
    return {
      name: "Error World",
      description: "Lore generation failed."
    };
  }
};

function getDominantBiome(world: WorldData): string {
  const counts: Record<string, number> = {};
  world.cells.forEach(c => counts[c.biome] = (counts[c.biome] || 0) + 1);
  return Object.entries(counts).sort((a,b) => b[1]-a[1])[0]?.[0] || 'Unknown';
}