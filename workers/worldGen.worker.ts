import { generateWorld } from '../utils/worldGen';
import { serializeWorld, WorldPayload } from '../utils/worldTransfer';
import { WorldParams } from '../types';

export type WorkerRequest = { type: 'generate'; params: WorldParams };
export type WorkerResponse =
  | { type: 'log'; msg: string }
  | { type: 'progress'; stage: number; total: number }
  | { type: 'done'; payload: WorldPayload }
  | { type: 'error'; message: string };

// This project's tsconfig includes the DOM lib, so bare `self` types as Window
// and `self.onmessage` will not accept a MessageEvent<WorkerRequest> handler.
// Adding `/// <reference lib="webworker" />` collides with DOM in the same
// program — the classic time-sink. One local alias avoids both.
const ctx = self as unknown as Worker;

// Abort is by terminate() from the client, never by message: a synchronous
// generation loop cannot drain the message queue, and draining it is exactly
// what deleting the setTimeout(0) yields gives up. See the client for why.
ctx.onmessage = async (e: MessageEvent<WorkerRequest>) => {
  if (e.data.type !== 'generate') return;
  try {
    const world = await generateWorld(
      e.data.params,
      msg => ctx.postMessage({ type: 'log', msg } satisfies WorkerResponse),
      undefined,
      (stage, total) => ctx.postMessage({ type: 'progress', stage, total } satisfies WorkerResponse),
    );
    const { payload, transfer } = serializeWorld(world);
    ctx.postMessage({ type: 'done', payload } satisfies WorkerResponse, transfer);
  } catch (err) {
    ctx.postMessage({ type: 'error', message: (err as Error).message } satisfies WorkerResponse);
  }
};
