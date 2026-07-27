import { describe, it, expect, vi } from 'vitest';
import { handleWorkerMessage, HandleWorkerMessageCallbacks } from '../utils/worldGenClient';
import type { WorkerResponse } from '../workers/worldGen.worker';
import type { WorldPayload } from '../utils/worldTransfer';

// Minimal payload for a single cell whose biome byte is out of range. The
// deserializer throws before touching any of the other typed arrays (it reads
// vertices, neighbors, presence/bools/biome first and throws immediately on
// an invalid biome byte), so those are the only fields that need real values;
// everything else is present but empty/zeroed to satisfy the type.
const makeInvalidBiomePayload = (): WorldPayload => ({
  cellCount: 1,
  height: new Float64Array(1),
  temperature: new Float64Array(1),
  moisture: new Float64Array(1),
  flux: new Float64Array(1),
  population: new Float64Array(1),
  plateId: new Int32Array(1),
  regionId: new Int32Array(1),
  provinceId: new Int32Array(1),
  cultureId: new Int32Array(1),
  religionId: new Int32Array(1),
  biome: new Uint8Array([255]), // out of range: no BiomeType maps to byte 255
  presence: new Uint8Array(1),
  bools: new Uint8Array(1),
  center: new Float64Array(3),
  vertOffsets: new Uint32Array([0, 0]),
  vertData: new Float64Array(0),
  nbrOffsets: new Uint32Array([0, 0]),
  nbrData: new Int32Array(0),
  ringOffsets: new Uint32Array([0, 0]),
  ringData: new Float64Array(0),
  geomPresent: new Uint8Array(1),
  propsPresent: new Uint8Array(1),
  geoSite: new Float64Array(2),
  geoSiteCoords: new Float64Array(2),
  geoNbrOffsets: new Uint32Array([0, 0]),
  geoNbrData: new Int32Array(0),
  riverOffsets: new Uint32Array([0, 0]),
  riverData: new Float64Array(0),
  hasRivers: false,
  lakeIdOffsets: new Uint32Array([0, 0]),
  lakeIdData: new Int32Array(0),
  featIdOffsets: new Uint32Array([0, 0]),
  featIdData: new Int32Array(0),
  params: {} as WorldPayload['params'],
});

describe('handleWorkerMessage', () => {
  it('rejects (does not hang) when a done message carries a payload that fails to deserialize', () => {
    // Regression test for the hang: deserializeWorld's throw happens inside
    // the message-handler callback, not inside the Promise executor's
    // synchronous frame, so without an explicit try/catch neither resolve
    // nor reject would ever be called and the awaiting caller would hang
    // forever.
    let settled = false;
    let resolved: unknown;
    let rejected: Error | undefined;
    const cb: HandleWorkerMessageCallbacks = {
      resolve: (w) => { resolved = w; },
      reject: (err) => { rejected = err; },
      isSettled: () => settled,
      finish: (fn) => {
        if (settled) return;
        settled = true;
        fn();
      },
    };

    const msg: WorkerResponse = { type: 'done', payload: makeInvalidBiomePayload() };
    handleWorkerMessage(msg, cb);

    expect(resolved).toBeUndefined();
    expect(rejected).toBeInstanceOf(Error);
    expect(rejected?.message).toMatch(/out-of-range biome byte/);
    expect(settled).toBe(true);
  });

  it('ignores a log/progress message that arrives after settlement', () => {
    let settled = true; // simulate: already settled (e.g. aborted or done)
    const onLog = vi.fn();
    const onProgress = vi.fn();
    const cb: HandleWorkerMessageCallbacks = {
      onLog,
      onProgress,
      resolve: vi.fn(),
      reject: vi.fn(),
      isSettled: () => settled,
      finish: (fn) => {
        if (settled) return;
        settled = true;
        fn();
      },
    };

    handleWorkerMessage({ type: 'progress', stage: 3, total: 10 }, cb);
    handleWorkerMessage({ type: 'log', msg: 'late message' }, cb);

    expect(onProgress).not.toHaveBeenCalled();
    expect(onLog).not.toHaveBeenCalled();
  });

  it('still delivers log/progress messages while not yet settled', () => {
    let settled = false;
    const onLog = vi.fn();
    const onProgress = vi.fn();
    const cb: HandleWorkerMessageCallbacks = {
      onLog,
      onProgress,
      resolve: vi.fn(),
      reject: vi.fn(),
      isSettled: () => settled,
      finish: (fn) => {
        if (settled) return;
        settled = true;
        fn();
      },
    };

    handleWorkerMessage({ type: 'progress', stage: 1, total: 10 }, cb);
    handleWorkerMessage({ type: 'log', msg: 'in progress' }, cb);

    expect(onProgress).toHaveBeenCalledWith(1, 10);
    expect(onLog).toHaveBeenCalledWith('in progress');
  });
});
