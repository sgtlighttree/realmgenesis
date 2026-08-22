# Dymaxion projection & Blender UV interop

RealmGenesis can unfold the sphere onto an icosahedron (the Fuller / Dymaxion
projection). The core lives in `utils/dymaxion.ts`; the on-screen 2D view is in
`components/Map2D.tsx`, the preview in `components/DymaxionPreview2D.tsx`, and the
PNG export in `utils/export.ts` (`exportDymaxionRaster`).

## Two net layouts

`DymaxionLayout` (`types.ts`) is `'classic' | 'blender'`. Both unfold the SAME 20
icosahedron faces; they differ only in how the faces are laid out in 2D.

| Layout | Shape | Built by | Use |
|--------|-------|----------|-----|
| `classic` | staircase / lightning-bolt | `buildDymaxionNet('classic')` — spanning-tree fold (`buildParents`) | a standalone poster-style map that fills the frame |
| `blender` | sawtooth strip | `buildBlenderNet()` — hardcoded UV table | drops straight onto a Blender default icosphere as a UV texture |

**`blender` is the default** (`hooks/useWorldEngine.ts`), so the view, preview,
and export all show the sawtooth unless the user unchecks the toggle. The default
was `classic` until 2026-08-22; that meant a user had to know to tick a box before
an export would map onto an icosphere, and a classic export can NEVER map onto the
icosphere UV — it is a different fold.

## The Blender UV net

`buildBlenderNet()` reproduces Blender's default icosphere UV unwrap exactly, so
an exported PNG is a ready-made UV texture. Key facts:

- **12 vertices, 20 faces.** UV bounds are **u ∈ [0, 1]**, **v ∈ [0, 0.472382]**
  (= 3·√3/11: three rows of equilateral triangles). The net fills the FULL width
  but only the bottom ~47% of vertical UV space — the rest of `[0,1]²` is unused,
  so a faithful export has a black upper band. That black is correct, not a bug;
  cropping it breaks the drop-on alignment.
- **Coordinate transform.** Blender is z-up; the engine is y-up. Vertices convert
  by `(x, y, z) → (x, z, y)`, so Blender `[0,0,-1]` (south pole) → engine
  `[0,-1,0]`, and `[0,0,1]` (north) → `[0,1,0]`.
- **U steps** are multiples of 1/11; **V steps** are multiples of √3/11.

The values are frozen against live-dumped Blender ground truth in
`tests/fixtures/dymaxionBlenderUV.json`, and `tests/dymaxionBlenderNet.test.ts`
asserts `buildBlenderNet` matches it to 1e-3. As of the 2026-08-22 dump it matches
to 5 decimals (max UV Δ 5e-6, max 3D Δ 6e-5).

## Why the export drops on cleanly

`exportDymaxionRaster` rasterizes the blender net at **`px = u·W`, `py = (1−v)·H`**
on a **square** canvas (`W = H = resolution`). That is the identical formula
Blender uses to sample a UV texture (`pixel = (u·imgW, (1−v)·imgH)`), and the
image aspect is irrelevant because both sides use the same mapping. So:

> Net UVs match Blender + export formula == Blender's sample formula ⇒ the PNG
> lands on the icosphere with zero scaling or moving.

The rasterizer walks each face's pixels, recovers the UV (`netPoint`), takes
barycentric weights over the face's 3D vertices, normalizes to a sphere
direction, rotates by the `lon/lat/roll` settings, and samples the
equirectangular render at that lon/lat.

## Re-extracting the ground truth from Blender

If Blender ever changes its default icosphere unwrap, re-dump and re-freeze the
fixture. In Blender: add a default Icosphere (1 subdivision), stay in **Object
Mode** (the object-mode data API reads empty in Edit Mode), then run:

```python
import bpy, json
me = bpy.data.objects['Icosphere'].data
uv = me.uv_layers.active.data
faces = []
for poly in me.polygons:
    loops = []
    for li in poly.loop_indices:
        vi = me.loops[li].vertex_index
        co = me.vertices[vi].co       # Blender z-up
        u, v = uv[li].uv
        loops.append({'v': vi, 'co': [round(co.x,5), round(co.y,5), round(co.z,5)],
                      'uv': [round(u,6), round(v,6)]})
    faces.append({'p': poly.index, 'loops': loops})
print(json.dumps({'nverts': len(me.vertices), 'nfaces': len(me.polygons), 'faces': faces}))
```

Save the output to `tests/fixtures/dymaxionBlenderUV.json`, then re-run the test.
If it fails, update `buildBlenderNet`'s UV/vertex tables to the new values.

## Related invariant

The Dymaxion pick buffer must mirror the visible rasterization — same net, same
rotation, same sizing — or clicks land on the wrong cell. See
[invariants.md](invariants.md).
