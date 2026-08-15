import { describe, it, expect } from "vitest";
import { generateWorld } from "../utils/worldGen";
import { makeParams } from "./helpers";

describe("WorldData.currents persistence", () => {
  it("present and deterministic when currentStrength > 0", async () => {
    const p = makeParams({ seed: "currents-persist", currentStrength: 1.0 });
    const a = await generateWorld(p);
    const b = await generateWorld(p);
    expect(a.currents).toBeDefined();
    expect(a.currents!.vx.length).toBe(a.cells.length);
    expect(a.currents!.sst.length).toBe(a.cells.length);
    expect(Array.from(a.currents!.vx)).toEqual(Array.from(b.currents!.vx));
    expect(Array.from(a.currents!.sst)).toEqual(Array.from(b.currents!.sst));
  });

  it("absent when currentStrength === 0 (escape hatch)", async () => {
    const p = makeParams({ seed: "currents-persist", currentStrength: 0 });
    const w = await generateWorld(p);
    expect(w.currents).toBeUndefined();
  });
});
