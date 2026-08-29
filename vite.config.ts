import path from 'path';
import { defineConfig, loadEnv } from 'vite';
import react from '@vitejs/plugin-react';

export default defineConfig(({ mode }) => {
    const env = loadEnv(mode, '.', '');
    return {
      server: {
        port: 3000,
        host: '0.0.0.0',
      },
      plugins: [react()],
      worker: {
        // ES-module worker: the generation closure uses ESM imports
        // (d3-geo-voronoi et al.) and the classic-worker format cannot.
        format: 'es',
      },
      define: {
        'process.env.GEMINI_API_KEY': JSON.stringify(env.GEMINI_API_KEY || '')
      },
      resolve: {
        alias: {
          '@': path.resolve(__dirname, '.'),
        }
      },
      test: {
        // V3 terrain generation is ~9s per world on this machine and slower on
        // CI runners; several suites generate multiple worlds per test. The
        // default 5s per-test timeout is far too tight now that V3 is the live
        // path, so raise the floor suite-wide (tests with an explicit timeout
        // still win).
        testTimeout: 120000,
        // LOCALLY, LEAVE THE MACHINE USABLE.
        //
        // Vitest defaults to the `forks` pool at `availableParallelism() - 1`
        // workers — 7 separate Node processes on this 8-core M1 Air, each
        // loading the whole module graph and holding generated worlds live.
        // Measured: 992s of test time in 218s wall, about 4.5 cores pinned for
        // three and a half minutes, on a fanless machine with 16GB shared
        // between CPU and GPU. That is fine on a dedicated runner and hostile
        // on a laptop that is also running After Effects — and it is what made
        // paramLiveness blow its timeout twice.
        //
        // 3 keeps a performance core free for the app in front of you. Raise it
        // per-run when the machine is idle (`--maxWorkers=7`), or drop to 1
        // while a render is going (`VITEST_WORKERS=1`). CI keeps the default,
        // where saturating the box is the point.
        maxWorkers: process.env.CI
          ? undefined
          : Number(process.env.VITEST_WORKERS ?? 3),
      },
    };
});
