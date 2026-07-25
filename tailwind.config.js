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
    extend: {},
  },
  plugins: [],
};
