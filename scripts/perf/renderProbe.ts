// Render-count probe with PROP ATTRIBUTION, for F4 React work.
//
// It lives under `scripts/perf/` and NOT under `hooks/` on purpose: nothing in
// the app imports it, so it is never bundled. To use it you temporarily add
// `useRenderProbe('Foo', { propA, propB })` to the components under study, take
// your measurements, and then REVERT those call sites. Committing the call
// sites is the mistake this location exists to discourage.
//
// Why attribution and not a bare counter: knowing that `Controls` rendered 40
// times during a slider drag tells you nothing you can act on. Knowing that it
// rendered 40 times and `factionColors` changed identity on 40 of them names
// the fix. `window.__rc.dump()` returns both, so a Playwright run can read it.
import { useRef } from 'react';

type Bag = Record<string, unknown>;

interface Probe {
  counts: Record<string, number>;
  keys: Record<string, Record<string, number>>;
  reset: () => void;
  dump: () => string;
}

const w = window as unknown as { __rc?: Probe };

if (!w.__rc) {
  const probe: Probe = {
    counts: {},
    keys: {},
    reset: () => {
      probe.counts = {};
      probe.keys = {};
    },
    dump: () => JSON.stringify({ counts: probe.counts, keys: probe.keys }),
  };
  w.__rc = probe;
}

export function useRenderProbe(name: string, props?: Bag): void {
  const prev = useRef<Bag | undefined>(undefined);
  const p = w.__rc!;
  p.counts[name] = (p.counts[name] || 0) + 1;
  if (props) {
    const before = prev.current;
    if (before) {
      const bucket = (p.keys[name] ||= {});
      for (const k of Object.keys(props)) {
        if (!Object.is(before[k], props[k])) bucket[k] = (bucket[k] || 0) + 1;
      }
    }
    prev.current = props;
  }
}
