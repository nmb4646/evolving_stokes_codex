# Remesh Solver Notes

This folder contains the MATLAB membrane-evolution solver. The main production driver is [ms_multi.m](./ms_multi.m). In practice the project is usually launched from [ms_batch.ipynb](./ms_batch.ipynb), which constructs MATLAB command lines for parameter sweeps and runs them in parallel.

## Solver flow

`ms_multi.m` expects several variables from the caller workspace:

- `p`: parameter struct
- `dir`: output directory ending in `/`
- `verbose`
- `supress_outputs`

The typical time step does the following:

1. Build or reload a triangulated closed surface.
2. Construct `Geometry(M, P)` for normals, areas, curvature, and DEC operators.
3. Assemble membrane forces and the viscous surface operators.
4. Solve the inner incremental-potential minimization loop for updated positions and tension.
5. Remove rigid motion and save `geo%d.mat`.
6. Remesh when triangle quality degrades, then interpolate velocity data to the new mesh.

## Files that matter most

- [ms_multi.m](./ms_multi.m): main semipermeable membrane / hydrodynamic evolution driver.
- [ms_batch.ipynb](./ms_batch.ipynb): Python notebook used to launch parallel MATLAB parameter sweeps.
- [Geometry.m](./Geometry.m): geometry, curvature, DEC operators, bending force.
- [Mesh.m](./Mesh.m): half-edge connectivity used by `Geometry`.
- [stokeslet_SLP_triangle.m](./stokeslet_SLP_triangle.m): current triangle-quadrature Stokes single-layer potential used by `ms_multi.m`.
- [stokeslet_SLP_triangle_setup.m](./stokeslet_SLP_triangle_setup.m): topology-only cache for the triangle SLP.
- [stokeslet_SLP_advanced.m](./stokeslet_SLP_advanced.m): opt-in higher-order triangle SLP for near-panel resolution experiments.
- [stokeslet_SLP_advanced_setup.m](./stokeslet_SLP_advanced_setup.m): cache builder for the advanced SLP.
- [stokeslet_SLP.m](./stokeslet_SLP.m): older vertex-lumped reference implementation retained for comparison.
- [deformation_criterion.m](./deformation_criterion.m): remeshing trigger.
- [interpolate.m](./interpolate.m), [rm_rigid.m](./rm_rigid.m): remesh transfer and post-step cleanup.

## Hydrodynamic discretization

The old `stokeslet_SLP.m` approximates the single-layer potential as a sum of point Stokeslets at vertices with lumped vertex areas. That is simple, but it collapses each surface patch to a point and handles on-surface singular behavior poorly.

The current `stokeslet_SLP_triangle.m` instead uses a hybrid rule:

- integrates the far field with a single centroid sample per triangle
- interpolates traction linearly from vertex data to each quadrature point
- evaluates the far field in vectorized target chunks for lower MATLAB overhead
- applies a Duffy-transformed self-face correction for triangles incident to the target vertex

This is a better match to the actual surface integral than the old vertex-lumped sum while staying much cheaper than full multi-point face quadrature in MATLAB. It is still not a full high-order boundary-integral implementation. In particular, near-singular interactions from non-incident neighboring faces are still handled by the centroid rule, not a dedicated near-panel rule.

There is also an opt-in `stokeslet_SLP_advanced.m` that uses higher-order regular quadrature, adaptive refinement for near panels, and higher-order Duffy quadrature for self panels. Keep `stokeslet_SLP_triangle.m` as the default path unless use of the advanced routine is specifically requested.

## Cache behavior

`stokeslet_SLP_triangle_setup(F)` depends only on connectivity. `ms_multi.m` rebuilds the cache:

- once after initial mesh creation or reload
- again after any remeshing event

That keeps the expensive topology bookkeeping out of the inner nonlinear loop.

## External dependency

The isotropic remesher lives in [isoremesh](./isoremesh) and depends on OpenMesh. See [readme_addition.md](./readme_addition.md) and `isoremesh/README.md` for the MEX build notes.

## Practical caution

`ms_multi.m` is script-style code, not a clean MATLAB function. It relies on caller workspace variables and has some research-code conventions baked in. When changing it, check both:

- fresh-start runs with `p.start = 0`
- resumed runs that reload `geo%d.mat`
