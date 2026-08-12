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
      },
    };
});
