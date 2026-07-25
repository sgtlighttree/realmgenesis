---
target: components/shell (F1 redesign shell UI)
total_score: 22
p0_count: 2
p1_count: 3
timestamp: 2026-07-25T03-53-02Z
slug: components-shell
---
Method: dual-agent (A: design review · B: deterministic detector + evidence)

## Design Health Score

| # | Heuristic | Score | Key Issue |
|---|-----------|-------|-----------|
| 1 | Visibility of System Status | 3 | `rg-rise` on the Do bar is the only state-change feedback; no confirmation after paint/edit lands |
| 2 | Match System / Real World | 2 | Three widget languages for three peer decisions in one strip |
| 3 | User Control and Freedom | 2 | Undo covers paint strokes only; Generate has no undo path |
| 4 | Consistency and Standards | 2 | PANEL/FOCUS chrome centralized, form-control vocabulary is not |
| 5 | Error Prevention | 2 | Generate World discards hand edits with no gate |
| 6 | Recognition Rather Than Recall | 3 | View strip keeps mode/layer/toggles visible; Do bar discloses progressively |
| 7 | Flexibility and Efficiency | 2 | No keyboard shortcuts for view-layer or edit-mode switching |
| 8 | Aesthetic and Minimalist Design | 3 | Clean composition; Sys tab erodes into one long unstructured scroll |
| 9 | Error Recovery | 2 | Console renders errors in the same green as success lines |
| 10 | Help and Documentation | 1 | One good affordance hint ("Space + drag"); nothing else narrates the flow |
| **Total** | | **22/40** | **Acceptable — significant improvements needed** |

## Audit Health Score (technical, parallel scan)

| # | Dimension | Score | Key Finding |
|---|-----------|-------|-------------|
| 1 | Accessibility | 1 | 5 of 58 buttons carry semantic ARIA; 44 rely on `title` |
| 2 | Performance | 3 | Memo/disposal discipline sound; MiniMap gated by unmount |
| 3 | Responsive Design | 1 | ~22px touch targets vs 44px min; horizontal scrollbar in Make rail |
| 4 | Theming | 2 | No token layer: 18 hard-coded hex, opacity spread /10–/85 |
| 5 | Anti-Patterns | 3 | Detector: 9 findings, all 9 false positives |
| **Total** | | **10/20** | **Acceptable — significant work needed** |

## Anti-Patterns Verdict

**LLM assessment:** Not slop at a glance; borderline slop under 30 seconds of clicking. The static composition is competent — solid-fill rail, quiet PANEL chrome, one accent, restrained iconography. It reads as a real tool's first draft. The tells surface on interaction: three control vocabularies for peer decisions in one strip, a mobile sheet that opens mid-scroll, and a dead control sitting under a live one.

**Deterministic scan:** 9 findings, one rule (`gray-on-color`). All 9 verified FALSE POSITIVES with a shared root cause: each is a ternary `active ? 'bg-blue-600 text-white' : 'bg-gray-800 text-gray-400'`, and the rule pairs the inactive branch's gray text with the active branch's saturated background — a combination that never renders. Detector blind spot, not a codebase bug.

**Visual overlays:** not injected; browser was reserved for Assessment A's live walkthrough.

## Priority Issues

**[P0] Generate World destroys hand edits with no gate.** A user who hand-painted terrain, biomes, or borders loses all of it in one click of the always-visible Generate button. It is the highest-stakes moment in the app and has the least reassurance. Fix: gate on `undoStack.length > 0`, or relabel contextually ("Regenerate — discards N edits").

**[P0] Mobile Make sheet opens at an arbitrary scroll position.** Verified live: tapping "Make" showed World Stats / Console / Generate, not Render Mode / Seed / Resolution. Root cause: the sheet body is the same DOM node across all four tabs, so `scrollTop` persists when the content swaps. Fix: reset scroll on tab change, or key the container by `openTab`.

**[P1] Buttons rely on `title` instead of accessible names.** 58 buttons, 5 with semantic ARIA. Every icon-only control is effectively unlabeled; the 15+ biome swatches read as "button, button, button". WCAG 4.1.2.

**[P1] Touch targets far below 44x44.** Strip chips are ~22px tall; EditToolbar mode buttons similar. These are primary controls on the narrow fold. WCAG 2.5.5.

**[P1] Three control vocabularies for one decision type.** Render mode (segmented buttons), view layer (native select), layer toggles (chips) sit three inches apart. Partially addressed: the view-layer select is now the themed `Select`; 7 native selects remain.

**[P2] Do bar exposes 8 modes + up to 4 sub-controls with grouping on only one family.** Terrain has a label; Biome/Political/Edit are ungrouped peers.

**[P2] Small text stacked on low contrast.** `text-[9px] text-gray-600` on `bg-black/85` compounds two WCAG risk factors. 30 uses of `text-gray-500`, 14 of `text-gray-600`.

## Cognitive Load

3 of 8 failures (single focus, one-thing-at-a-time, minimal choices), 2 partial. Decision points above the working-memory limit: view layer 12 options, Do bar 8 modes, overlay checkboxes 7, biome swatches 15+.

## Persona Red Flags

**Alex (power user):** no keyboard shortcuts for view-layer or edit-mode switching — mouse-only except Esc and Ctrl+Z. Generate is a one-click landmine. The 12-item layer list needs a full open-scan-click cycle where a Blender/QGIS user expects number-key cycling.

**Sam (accessibility):** native checkboxes and selects had no custom focus treatment against the `FOCUS` ring everything else gets (now addressed in `index.css`). Biome swatches and faction chips carry `title` only — no `aria-label`, so a screen reader gets unlabeled buttons. The mobile sheet's close button is the one double-covered a11y win and should be the model.

## What's Working

1. **`PANEL`/`FOCUS` shared constants** — one source of truth for chrome and focus rings; every docked panel moves together. The right lever, pulled correctly.
2. **The contextual Do bar** — summoning the edit toolset only while editing, with a genuine entrance, makes mode legible as a change of shape rather than a hidden flag.
3. **Faction swatches compute their own contrast text** — luminance-based label flip so the number stays legible on any faction color. Real craft.
4. **Performance discipline** — every `useMemo` geometry has a matching disposal effect; MiniMap redraw debounced and gated by unmount, not CSS-hidden.

## Questions to Consider

1. If the Do bar's whole premise is "summon tools only when needed," why does the Sys tab disclose everything permanently?
2. `shellKit`'s stub panels are now used only by `?shell=stub`. Dead weight in the shipped bundle, or kept as the layout reference?
3. Why does exactly one complex interaction get an inline hint ("Space + drag to orbit") and no other?
