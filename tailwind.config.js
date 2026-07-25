import colors from 'tailwindcss/colors';

/** @type {import('tailwindcss').Config} */
export default {
  content: [
    './index.html',
    './App.tsx',
    './components/**/*.tsx',
  ],
  theme: {
    /**
     * Sharp corners are a product decision, not an oversight: RealmGenesis is a
     * technical instrument, and square panels read as precise where rounded
     * ones read as consumer-soft. Zeroing the SCALE rather than stripping
     * `rounded-*` from every call site means it is one switch to revisit at the
     * F1b brand pass, and no future component can quietly reintroduce a radius.
     *
     * `full` is deliberately preserved: it expresses "this is a circle" (the
     * placeholder globe, any future avatar or dot), which is a different idea
     * from corner rounding.
     */
    borderRadius: {
      none: '0',
      sm: '0',
      DEFAULT: '0',
      md: '0',
      lg: '0',
      xl: '0',
      '2xl': '0',
      '3xl': '0',
      full: '9999px',
    },
    extend: {
      /**
       * SEMANTIC COLOR TOKENS — the one block F1b edits.
       *
       * Why this exists: the chrome reaches for raw palette steps (`gray-400`
       * 137 times, `gray-700` 97, `white` 81, `gray-800` 71), so a brand pass
       * without this layer is a find-and-replace across a dozen files with no
       * way to tell "this gray is a border" from "this gray is muted text".
       * Naming the ROLE lets F1b move borders without touching type.
       *
       * Every value below is the byte-identical palette step the code already
       * used, so adopting a token is provably a no-op render. That is the
       * invariant to check when converting a file: same pixels, new names.
       *
       * The ink ramp has SIX steps, which is more than a system wants. That is
       * deliberate — it is a faithful census of what the app actually uses, not
       * a proposal. Collapsing it is taste work, and taste work is F1b's job;
       * doing it here would smuggle a design change into a rename.
       *
       * Strays NOT mapped on purpose (they are single uses and each is either a
       * one-off or an accident for F1b to adjudicate): slate-800/400,
       * neutral-900, sky-500, emerald-500, blue-950.
       */
      colors: {
        // Depth ramp for panel and control fills, darkest to lightest.
        surface: {
          sunken: colors.gray[950],   // rail backing, behind panels
          DEFAULT: colors.gray[900],  // panel body — the workhorse
          raised: colors.gray[800],   // chips, buttons, inputs on a panel
          hover: colors.gray[700],    // raised control, pointer over it
        },
        // Hairlines. `edge-subtle` divides inside a panel; `edge` bounds a control.
        edge: {
          subtle: colors.gray[800],
          DEFAULT: colors.gray[700],
          strong: colors.gray[600],
        },
        // Type ramp, brightest to dimmest. See the six-step note above.
        ink: {
          strong: colors.white,       // active/selected label
          bright: colors.gray[100],   // panel titles
          DEFAULT: colors.gray[200],  // body copy
          soft: colors.gray[300],     // secondary values
          muted: colors.gray[400],    // field labels — the workhorse
          faint: colors.gray[500],    // captions, units, hints
        },
        // Single accent. Selection and state ONLY, never decoration — the
        // 4-hue "rainbow" bucket dots were removed in Session 5 for that reason.
        brand: {
          soft: colors.blue[400],
          DEFAULT: colors.blue[500],
          strong: colors.blue[600],
          deep: colors.blue[700],
        },
        // Edit mode is amber throughout; destructive confirmation is red.
        warn: { DEFAULT: colors.amber[500], soft: colors.amber[400] },
        danger: { DEFAULT: colors.red[500], soft: colors.red[400] },
      },
      /**
       * Z-SCALE — names the two-tier stack Session 6c settled on but only ever
       * wrote down in prose, plus the two layers above it that already existed
       * as bare `z-30`/`z-40`. Informal integers are how the WorldViewer pause
       * control ended up painting over the mobile sheet.
       */
      zIndex: {
        overlay: '10',  // canvas-owned: pause control, ruler readout, HUD
        chrome: '20',   // shell surfaces: rails, strips, docked panels
        sheet: '30',    // narrow-fold bottom sheets
        modal: '40',    // ConfirmDialog and anything else that takes focus
      },
    },
  },
  plugins: [],
};
