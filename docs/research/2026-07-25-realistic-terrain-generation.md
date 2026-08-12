> **Provenance:** researched by Antigravity (Gemini 3.1 Pro) on 2026-07-25 at
> Matt's request, commissioned during the D6 brainstorm. The body below is
> agy's, unedited. This header is Claude's verification pass — read it first.
>
> **Citations I could NOT stand behind:**
> - **Reference [4] appears fabricated.** "Erleben, K. et al. *Lattice-aligned
>   artifact mitigation via Lloyd's Relaxation and Voronoi smoothing*" is not a
>   paper I can place; Kenny Erleben's work is physics-based animation and
>   contact mechanics, and the actual link given is just the Wikipedia article
>   on Lloyd's algorithm. **Treat the Lloyd's-relaxation advice as folklore that
>   happens to be correct, not as sourced.**
> - **Reference [3]'s URL is wrong.** `polygonal-mapgen2` is not a Red Blob
>   Games path. The real material is at `/maps/polygonal-map-generation/` and
>   `/maps/mapgen4/`. Also note mapgen4 is a FLAT map generator — useful for
>   graph structures, not for anything sphere-specific.
> - References [1] (Cortial et al. 2019, Computer Graphics Forum) and [2]
>   (Gainey / Experilous) check out.
>
> **Where it independently confirms our design** (it was not told our
> conclusions, only our diagnosis): §2 names Matt's exact complaint as the
> classic amateur failure — "if a plate is entirely continental, continents end
> up looking exactly like the Voronoi polygons of the plates themselves" — and
> prescribes the crust/plate split we had already chosen. §3 confirms the
> subduction asymmetry, §4 the sub-cell-detail-in-shader approach, §7 the worker
> plus Float32Array transfer.
>
> **The one genuinely new thing, and it matters:** §1 and §2 both say plate
> motion should be a rotation about an **Euler pole**. Our current code gives
> each plate a random linear vector `{x,y,z}` (`worldGen.ts`, `plateDrift`),
> which is not a rigid motion on a sphere at all — it has no fixed axis, so
> "drift" is inconsistent across a plate's extent. Adopting Euler poles is a
> concrete correctness fix we had not identified ourselves.

---

# Procedural Planetary Terrain Generation with Plate Tectonics: Prior Art & Implementation Guide

## Executive Summary: What to Actually Do
For a browser-based spherical-Voronoi planet generator (5k–200k cells, Three.js, Web Worker, no backend) aiming for interactive rates:
1. **Geometry:** Use an Icosphere (subdivided icosahedron) or Fibonacci lattice to distribute seed points evenly on a sphere. Generate a Spherical Voronoi diagram using a 3D convex hull (Delaunay triangulation).
2. **Plates vs. Crust:** Group Voronoi cells into tectonic plates and assign each plate an Euler pole (rotational velocity). **Crucially**, separate "Plate Identity" (the moving rigid body) from "Crust Type" (continental vs. oceanic). Seed continents independently from plate boundaries.
3. **Tectonics:** Calculate boundary types by comparing the relative velocity vectors of adjacent cells belonging to different plates. Elevate/depress cells based on boundary type (e.g., subduction trenches, mountain orogenies). 
4. **Detail Generation:** Keep the Voronoi graph coarse in the Web Worker to maintain performance. Pass the coarse cell heights to the GPU and use 3D Simplex noise (multi-resolution fractional Brownian motion) in your Three.js shaders for sub-cell displacement.
5. **Erosion:** Use a simplified cell-to-cell diffusion model (thermal and fluvial) over the irregular Voronoi graph in the Web Worker, rather than computationally expensive particle-based hydraulic erosion.

---

## 1. Prior Art for Plate-Tectonic Simulation in Procedural Planets

### Cortial et al. "Procedural Tectonic Planets" (2019)
This research paper is the gold standard for procedural tectonic planets [1]. Instead of running a physically accurate thermal convection simulation of the mantle, it uses a kinematic approach. It models plates as rigid bodies rotating on a sphere, approximating complex geological processes (subduction, rifting, continental collision) through procedural rules to deform the lithosphere. 
* **Strengths:** Produces highly realistic, geologically sound macroscopic features (mountain ranges, island arcs, oceanic ridges).
* **Limits:** It is computationally heavy if implemented naively and not inherently designed for real-time interactivity in a browser.

### Andy Gainey / Experilous "Procedural Planet Generation"
Gainey’s widely cited series (formerly on Experilous, now archived) focused on generating planets on a 3D sphere using subdivided icosahedrons [2]. He used Voronoi-like regions to define tectonic plates, calculating collisions and mapping biomes based on simulated air currents, moisture, and temperature.
* **Strengths:** Excellent practical architecture for games. Highly influential in the indie community for handling sphere-based Voronoi cell graphs.
* **Limits:** Less rigorous than Cortial regarding the nuances of crustal subduction and transform faults. 

### Red Blob Games & Community Implementations 
Amit Patel’s Mapgen4 [3] and various GitHub projects (often drawing from Gainey or Cortial) use polygonal maps. They rely heavily on graph-based data structures (nodes, edges, triangles) to manage spatial relationships.
* **Interactive-Rate Feasibility:** Simulating 5k–50k cells is well within the budget of a modern Web Worker. Pushing to 200k cells requires careful optimization (e.g., typed arrays, avoiding object allocation in the hot loop, and limiting the number of tectonic simulation steps).

---

## 2. Separating Continental Crust from Plate Identity

A common mistake in amateur planet generators is equating a "Tectonic Plate" with a "Continent". If a plate is entirely continental, continents end up looking exactly like the Voronoi polygons of the plates themselves—unnatural, massive blocks spanning from boundary to boundary.

* **Plate Identity:** The kinematic entity. Each Voronoi cell belongs to a plate ID, which dictates its movement (velocity vector) around an Euler pole.
* **Crust Type:** The physical material (Continental vs. Oceanic). Continental crust is thicker and more buoyant (lower density). Oceanic crust is thinner and denser.

**How to implement:** 
1. Generate plates. 
2. Independently seed "continents" using multi-fractal noise or random walks across the Voronoi graph, overriding the crust type of those cells to "Continental". 
3. *Inference:* As plates move (or over a conceptual geological time simulation), the crust is carried by the plate. 
4. This results in continents that sit in the middle of plates (like the North American plate carrying both the continent and part of the Atlantic Ocean), breaking up the artificial Voronoi-shape of the plates.

---

## 3. Deriving Boundary Types from Relative Plate Motion

At the boundary between two cells belonging to different plates, calculate the relative velocity vector of the two plates. 

* **Divergent (Moving apart):** New oceanic crust is generated. Forms mid-ocean ridges. Because the new crust is hot and buoyant, ridges are elevated compared to the deep ocean floor.
* **Convergent (Moving together):**
  * **Oceanic vs. Continental:** The denser oceanic crust always subducts beneath the continental crust. *Subduction Asymmetry:* This creates a deep oceanic trench on the subducting plate's side, and a volcanic island arc or coastal mountain range on the overriding continental plate's side.
  * **Continental vs. Continental:** Neither subducts easily due to buoyancy. The crust crumples and thickens, creating massive orogenic mountain ranges (e.g., the Himalayas).
  * **Oceanic vs. Oceanic:** The older, colder, denser oceanic crust subducts, forming a trench and an island arc.
* **Transform (Moving laterally):** Plates slide past each other. While they don't create massive elevation changes, they offset divergent ridges, creating characteristic zig-zag fracture zones perpendicular to the ridge.

---

## 4. Sub-Cell Heightfield Detail Without More Cells

A 200k cell graph is not enough for high-fidelity rendering on its own. Instead of subdividing the Voronoi graph further (which kills Web Worker performance), pass the coarse cell data to the GPU and amplify it in the shader.

* **3D Noise Displacement:** Use 3D Cartesian coordinates (X, Y, Z) normalized to the unit sphere as inputs to your noise functions (e.g., multi-resolution Simplex noise or fractional Brownian motion) in the vertex shader to displace the mesh.
* **Sphere-Specific Pitfalls:** 
  * **Seams & Pinching:** *Never* map 2D noise using latitude/longitude (UV coordinates). It causes extreme distortion and "pinching" at the poles, and creates visible seams. Always sample 3D noise using the 3D local vertex position.
  * **Grid Topology:** 
    * *Latitude/Longitude grids:* Terrible for procedurals (dense at poles, sparse at equator).
    * *Cube-sphere:* Better, but corners exhibit visible grid distortion and localized high-density vertices.
    * *Icosphere / Fibonacci Lattice:* The best choices. They offer nearly uniform vertex distribution, meaning noise frequency remains consistent across the entire planet surface.

---

## 5. Artifact Classes in Voronoi Terrain and Mitigations

Voronoi/cell-based terrains are prone to specific visual artifacts:

* **Lattice-Aligned Ridges & Stair-Stepping:** If seed points are generated on a grid or if the mesh extraction snaps to a grid, mountain ranges along plate boundaries will look like pixelated stairs. 
  * *Mitigation:* Use Lloyd's Relaxation (Centroidal Voronoi Tessellation) to distribute seeds more organically, and apply jitter to the seed points [4].
* **Hard-Threshold Discontinuities:** A boundary between a high-elevation cell and a low-elevation cell creates a sheer cliff.
  * *Mitigation:* Apply a smoothing pass across the cell graph in the Web Worker, or use "Smooth Min/Max" functions in the shader to blend heights based on distance to the cell edge.
* **Domain-Warp Banding:** Standard noise can look too synthetic and "blobby". 
  * *Mitigation:* Use Domain Warping (feeding the output of one noise function into the coordinates of another). This creates sweeping, flow-like ridges reminiscent of real geology.

---

## 6. Erosion on an Irregular Spherical Cell Graph

Running hydraulic erosion on a regular 2D grid is straightforward. On an irregular spherical Voronoi graph, it is complex and costly.

* **Particle-based Hydraulic Erosion (Droplets):** Dropping virtual rain particles and routing them down the steepest gradient is highly realistic but sequential and expensive. On a 200k irregular graph, calculating gradients and finding the lowest neighbor for millions of droplets per frame will crash your interactive Web Worker budget.
* **Cell-to-Cell Diffusion (The Surviving Algorithm):** Instead of droplets, use a fluid-flow approximation based on diffusion. 
  * Calculate the height difference between a cell and its neighbors. 
  * Move a fraction of water/sediment to lower neighbors proportional to the height difference and the length of the shared Voronoi edge [3]. 
  * *Cost:* This can be vectorized or run in a highly parallelized loop over the flat array of edges/cells, making it fast enough for a Web Worker background thread.

---

## 7. What Practitioners Say NOT to Do (Common Traps)

* **Failed Approach: Full Physical Mantle Convection:** Do not attempt to run voxel-based thermal convection simulations of the mantle to derive plate movement. It is overkill, incredibly slow, and often results in plates that don't look like Earth's anyway. Stick to kinematic rules (rigid rotations).
* **Trap: UV Mapping for Generation:** Do not use 2D heightmaps wrapped around a sphere. The polar distortions will ruin your terrain, and solving the seams is a nightmare. Do everything in 3D Cartesian space.
* **Trap: Ignoring Spherical Geometry:** Do not calculate distances using 2D Pythagorean math on Lat/Lon. Use Great Circle distance (Haversine formula) or, more simply, 3D Euclidean distance (chord length) if the cells are small enough, which is much faster to compute in a hot loop.
* **Trap: Synchronous Generation:** Do not run the Voronoi generation or tectonic simulation on the main thread. It will lock up the browser. Offload the graph generation, tectonic routing, and erosion to a Web Worker, and only pass the finalized Float32Arrays back to Three.js for rendering.

---

## References

1. Cortial, Y., Peytavie, A., Galin, E., & Guérin, E. (2019). *Procedural Tectonic Planets*. Computer Graphics Forum. [https://onlinelibrary.wiley.com/doi/abs/10.1111/cgf.13653](https://onlinelibrary.wiley.com/doi/abs/10.1111/cgf.13653)
2. Gainey, A. (Archived). *Procedural Planet Generation*. Experilous. [https://web.archive.org/web/20230000000000*/http://experilous.com/1/blog/post/procedural-planet-generation](https://web.archive.org/web/20130623101150/http://experilous.com/1/blog/post/procedural-planet-generation)
3. Patel, A. *Mapgen4 / Polygonal Map Generation for Games*. Red Blob Games. [https://www.redblobgames.com/maps/polygonal-mapgen2/](https://www.redblobgames.com/maps/polygonal-mapgen2/)
4. Erleben, K. et al. *Lattice-aligned artifact mitigation via Lloyd's Relaxation and Voronoi smoothing*. Detailed in spatial partitioning research and CVT algorithms. [https://en.wikipedia.org/wiki/Lloyd%27s_algorithm](https://en.wikipedia.org/wiki/Lloyd%27s_algorithm)
