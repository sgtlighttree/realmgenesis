import { WorldData, WorldParams } from '../types';
import { deserializeWorld } from './worldTransfer';
import type { WorkerRequest, WorkerResponse } from '../workers/worldGen.worker';
import WorldGenWorker from '../workers/worldGen.worker?worker';

// Drop-in replacement for generateWorld with the identical signature, so the
// call sites in useWorldEngine change by one identifier.
//
// One worker per generation, terminated on abort. A long-lived worker would
// need a message-based abort, and a synchronous generation loop cannot receive
// messages — which is the whole reason the setTimeout(0) yields can go away.
// Spawn cost is ~1-5ms against a multi-second generation.
export const generateWorldInWorker = (
  params: WorldParams,
  onLog?: (msg: string) => void,
  signal?: AbortSignal,
  onProgress?: (stage: number, total: number) => void,
): Promise<WorldData> =>
  new Promise<WorldData>((resolve, reject) => {
    if (signal?.aborted) { reject(new Error('Generation Cancelled')); return; }

    const worker = new WorldGenWorker();
    let settled = false;
    const finish = (fn: () => void) => {
      if (settled) return;
      settled = true;
      signal?.removeEventListener('abort', onAbort);
      worker.terminate();
      fn();
    };
    const onAbort = () => finish(() => reject(new Error('Generation Cancelled')));
    signal?.addEventListener('abort', onAbort);

    worker.onmessage = (e: MessageEvent<WorkerResponse>) => {
      const m = e.data;
      if (m.type === 'log') onLog?.(m.msg);
      else if (m.type === 'progress') onProgress?.(m.stage, m.total);
      else if (m.type === 'done') finish(() => resolve(deserializeWorld(m.payload)));
      else if (m.type === 'error') finish(() => reject(new Error(m.message)));
    };
    worker.onerror = (e: ErrorEvent) => finish(() => reject(new Error(e.message || 'Worker failed')));

    worker.postMessage({ type: 'generate', params } satisfies WorkerRequest);
  });
