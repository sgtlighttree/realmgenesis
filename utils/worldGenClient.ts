import { WorldData, WorldParams } from '../types';
import { deserializeWorld } from './worldTransfer';
import type { WorkerRequest, WorkerResponse } from '../workers/worldGen.worker';
import WorldGenWorker from '../workers/worldGen.worker?worker';

export interface HandleWorkerMessageCallbacks {
  onLog?: (msg: string) => void;
  onProgress?: (stage: number, total: number) => void;
  resolve: (world: WorldData) => void;
  reject: (err: Error) => void;
  isSettled: () => boolean;
  finish: (fn: () => void) => void;
}

// Pure-ish message dispatcher, factored out of generateWorldInWorker so it is
// unit-testable without mocking the `?worker` import: it takes a message and
// a small set of callbacks rather than reaching into worker/promise state
// directly.
//
// Both non-terminal branches ('log', 'progress') are guarded on `isSettled()`
// so a message already queued on the main thread before an abort/settlement
// cannot fire callbacks for a generation that is already done or cancelled.
//
// The 'done' branch's `deserializeWorld` call can throw (e.g. an out-of-range
// biome byte) — that throw happens inside this handler, not inside the
// Promise executor's synchronous frame, so it would NOT auto-reject the
// promise on its own. It is wrapped in try/catch so a deserialize failure
// rejects instead of leaving the promise permanently unsettled.
export const handleWorkerMessage = (
  m: WorkerResponse,
  cb: HandleWorkerMessageCallbacks,
): void => {
  if (m.type === 'log') {
    if (cb.isSettled()) return;
    cb.onLog?.(m.msg);
  } else if (m.type === 'progress') {
    if (cb.isSettled()) return;
    cb.onProgress?.(m.stage, m.total);
  } else if (m.type === 'done') {
    cb.finish(() => {
      try {
        cb.resolve(deserializeWorld(m.payload));
      } catch (err) {
        cb.reject(err as Error);
      }
    });
  } else if (m.type === 'error') {
    cb.finish(() => cb.reject(new Error(m.message)));
  }
};

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
      // Null out the handlers before terminate(): terminate() stops the
      // worker but does not retract messages already queued on the main
      // thread's task queue, so a stale onmessage could still fire after
      // settlement without this.
      worker.onmessage = null;
      worker.onerror = null;
      worker.onmessageerror = null;
      worker.terminate();
      fn();
    };
    const onAbort = () => finish(() => reject(new Error('Generation Cancelled')));
    signal?.addEventListener('abort', onAbort);

    worker.onmessage = (e: MessageEvent<WorkerResponse>) => {
      handleWorkerMessage(e.data, {
        onLog,
        onProgress,
        resolve,
        reject,
        isSettled: () => settled,
        finish,
      });
    };
    worker.onerror = (e: ErrorEvent) => finish(() => reject(new Error(e.message || 'Worker failed')));
    // Closes the last unsettled exit path: a message that fails structured-clone
    // deserialization on the main thread fires onmessageerror, not onmessage or
    // onerror. Without this handler that message settles nothing, and the
    // promise (and isGenerating) stays pinned forever.
    worker.onmessageerror = () => finish(() => reject(new Error('Worker message failed to deserialize')));

    worker.postMessage({ type: 'generate', params } satisfies WorkerRequest);
  });
