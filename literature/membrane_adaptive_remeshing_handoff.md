# Handoff: Adaptive Surface Remeshing and Volume Control for Closed Membrane Simulations

## Purpose

This document is a technical handoff for replacing or supplementing the current remeshing backend in a MATLAB membrane-dynamics code.

The membranes are represented as closed triangulated surfaces:

- `P`: `N x 3` vertex positions
- `M`: `F x 3` triangle connectivity

The remesher is called repeatedly during time-dependent simulations.

Current backend:

- `christopherhelf/isotropicremeshing`
- https://github.com/christopherhelf/isotropicremeshing
- used through a MATLAB MEX wrapper

It has worked reasonably well, but new simulations create two requirements:

1. **Strict control of enclosed volume**
2. **Spatially adaptive mesh density**, typically determined by a scalar field such as curvature

The goal is to introduce a modern remeshing backend while keeping the MATLAB-side interface simple and while separating remeshing error from constraint-projection error.

---

# 1. Terminology: adaptive isotropic remeshing

The requirement is not necessarily directional anisotropy.

What is needed is a scalar target edge length varying over the surface,

\[
h_i = h(\mathbf{x}_i),
\]

for example,

\[
h_i = h(H_i),
\]

where \(H_i\) is a local curvature measure.

Thus the desired operation is best described as **adaptive isotropic remeshing**.

"Isotropic" means that at a given location the desired element size is the same in all tangent directions. It does not imply a globally uniform mesh density.

Typical use:

- low-curvature regions -> larger triangles
- high-curvature necks, buds, pearls -> smaller triangles

True tensor-metric anisotropy is not currently required.

---

# 2. Recommended replacement: MMGS / MMG

Repository:

https://github.com/MmgTools/mmg

Documentation:

https://mmgtools.github.io/

MMG provides several mesh-adaptation libraries. `MMGS` is the component for embedded triangular surface meshes.

Why it is attractive:

- purpose-built for mesh adaptation
- supports spatially varying scalar mesh size / metric information
- C API is suitable for a thin MEX wrapper
- mesh can be supplied and recovered in memory
- explicit controls on minimum/maximum size and metric gradation
- actively maintained
- LGPL licensed
- cleaner long-term dependency than the old OpenMesh MATLAB wrapper

Desired MATLAB interface:

```matlab
[P2, M2, info] = remesh_surface(P, M, h, opts);
```

with

```matlab
size(P) = [N,3];
size(M) = [F,3];
size(h) = [N,1];
```

and `h(i)` representing the desired local physical edge length.

The MEX wrapper should hide MMGS-specific conventions from MATLAB.

---

# 3. Other libraries considered

## Geometry Central

Documentation:

https://geometry-central.net/surface/algorithms/remeshing/

Geometry Central implements incremental remeshing using:

- edge split
- edge collapse
- edge flip
- tangential smoothing

Its remeshing supports curvature adaptation. Relevant options include concepts such as

```cpp
targetEdgeLength
curvatureAdaptation
minRelativeLength
maxIterations
```

When curvature adaptation is enabled, target lengths become smaller in high-curvature regions.

### Advantages

- clean modern C++ library
- good mesh data structures
- relatively straightforward Eigen interoperability
- access to individual mesh mutation operations
- attractive if built-in curvature adaptation is sufficient

### Disadvantage for this project

Its standard remeshing path is less naturally organized around an arbitrary MATLAB-provided per-vertex scalar `h` field than MMGS.

MMGS is therefore the first choice if MATLAB should remain responsible for constructing the sizing field.

---

## PMP Library

Documentation:

https://www.pmp-library.org/remeshing.html

PMP implements:

- edge collapse
- edge split
- edge flip
- tangential relaxation
- optional projection to the original surface

Its adaptive-remeshing routine varies element size according to local curvature.

Typical controls include:

- `min_edge_length`
- `max_edge_length`
- `approx_error`
- `iterations`
- `use_projection`

### Advantages

- simple API
- mature surface-mesh library
- built-in adaptive remeshing
- plausible MEX candidate

### Disadvantage

Its standard adaptive routine is more opinionated about how the sizing field is generated.

MMGS remains preferable if the solver should compute the scalar field itself.

---

## CGAL Polygon Mesh Processing

Relevant documentation:

https://doc.cgal.org/latest/Polygon_mesh_processing/classPMPSizingField.html

https://doc.cgal.org/latest/PMP_Remeshing/classCGAL_1_1Polygon__mesh__processing_1_1Adaptive__sizing__field.html

CGAL supports isotropic remeshing driven by a sizing-field abstraction. It also provides extensive geometry utilities:

- self-intersection checks
- curvature estimation
- mesh repair
- orientation tools
- Hausdorff-distance tools
- closed-mesh volume calculations

### Advantages

Very broad and robust geometry infrastructure.

### Disadvantages

- heavier dependency
- more difficult C++/template-heavy MEX integration
- licensing requires more care
- likely excessive if the immediate need is only adaptive remeshing

Revisit CGAL if strong mesh-repair or intersection-detection tools become necessary.

---

# 4. Replacing the remesher will not automatically solve volume conservation

General-purpose remeshers should not be assumed to preserve enclosed volume exactly.

A remesher generally performs some combination of:

- edge splitting
- edge collapsing
- edge flipping
- vertex smoothing
- reprojection
- vertex relocation

Some purely topological/local operations can preserve the represented polyhedral surface exactly, but useful mesh adaptation generally perturbs geometry.

Therefore exact volume conservation should remain the responsibility of the membrane code.

The useful question is not

> Does the remesher preserve volume exactly?

but

> How much geometric/volume error does the remesher create before correction, and how large a correction is required afterward?

This must be measured explicitly.

---

# 5. Diagnose the current pipeline before replacing anything

There are two candidate failure sources:

1. remeshing itself
2. post-remesh volume projection

Instrument the current code so these can be separated.

## State A: before remeshing

Store

\[
P_0,\qquad M_0
\]

and compute

\[
A_0=A(P_0,M_0),
\]

\[
V_0=V(P_0,M_0).
\]

Also record mesh-quality diagnostics.

## State B: raw remesher output, before any correction

Let

\[
P_r,M_r
\]

be the raw remeshed surface.

Compute

\[
A_r=A(P_r,M_r),
\]

\[
V_r=V(P_r,M_r),
\]

and

\[
\epsilon_V^{\rm raw}
=
\frac{|V_r-V_0|}{|V_0|},
\]

\[
\epsilon_A^{\rm raw}
=
\frac{|A_r-A_0|}{A_0}.
\]

Record at least:

- vertex count
- face count
- minimum triangle angle
- maximum triangle aspect ratio
- minimum/maximum edge length
- mean/median edge length
- degenerate triangles
- orientation consistency
- self-intersection status if available
- geometric distance from the old surface if practical

Useful surface-change measures include:

- nearest-point distance
- symmetric Hausdorff distance
- maximum normal displacement

## State C: after area/volume correction

Let

\[
P_c,M_r
\]

be the corrected result.

Record the same quantities again, plus:

- maximum correction displacement
- RMS correction displacement
- displacement relative to local edge length
- mesh quality after projection
- curvature change induced by projection

---

# 6. Interpretation of the diagnostic

## Case 1: raw remeshing causes the problem

If raw remeshing produces:

- large \(\epsilon_V^{\rm raw}\)
- large \(\epsilon_A^{\rm raw}\)
- strong curvature changes
- collapse of thin necks
- poor triangles
- large displacements relative to local edge length

then the remesher or its parameters are likely responsible.

Especially suspicious is raw geometric error that does not improve with input mesh refinement.

## Case 2: raw remeshing is benign but projection damages geometry

If raw remeshing is geometrically mild, but artifacts appear after projection, the projector is likely responsible.

Test this independently by perturbing a good closed mesh slightly and applying only the constraint projector.

---

# 7. Exact polyhedral volume

For a consistently oriented closed triangular mesh,

\[
V(P,M)
=
\frac16
\sum_{(i,j,k)\in M}
\mathbf p_i\cdot
\left(
\mathbf p_j\times\mathbf p_k
\right).
\]

The sign depends on orientation.

Use the same volume definition for:

- diagnostics
- constraint residual
- gradient calculation

---

# 8. Potential issue with the current minimum-displacement projection

Suppose the present projector solves the unweighted Euclidean problem

\[
\min_{\delta\mathbf x}
\frac12
\|\delta\mathbf x\|_2^2
\]

subject to the linearized volume constraint

\[
\nabla V^T\delta\mathbf x
=
V_0-V.
\]

Then

\[
\delta\mathbf x
=
\frac{V_0-V}
{\nabla V^T\nabla V}
\nabla V.
\]

This is a minimum displacement in vertex-coordinate Euclidean norm, but that norm is not invariant to mesh density.

This becomes important once the mesh is adaptive.

---

# 9. Why the unweighted projection is mesh-density dependent

Let

\[
\mathbf g_V=\nabla_P V.
\]

For a closed triangular mesh, the volume gradient at vertex \(i\) scales approximately like

\[
\mathbf g_{V,i}
\sim
A_i\mathbf n_i,
\]

where

- \(A_i\) is a lumped/associated vertex area
- \(\mathbf n_i\) is the local outward normal

Thus the unweighted correction behaves roughly like

\[
\delta\mathbf p_i
\propto
A_i\mathbf n_i.
\]

Vertices in coarse regions represent more area per vertex and therefore move more than vertices in refined regions.

This means the correction depends on the discretization.

Example:

- highly curved neck -> many small triangles
- low-curvature lobes -> fewer large triangles

An unweighted Euclidean minimum-displacement correction tends to move the coarse lobes more than the highly refined neck.

That is a discretization artifact.

---

# 10. Recommended alternative: mass-weighted projection

Instead minimize a discrete approximation to

\[
\frac12
\int_\Gamma
|\delta\mathbf x|^2\,dA.
\]

Using a lumped surface mass matrix,

\[
\min_{\delta\mathbf x}
\frac12
\delta\mathbf x^T
\mathbf W
\delta\mathbf x
\]

subject to

\[
\mathbf g_V^T\delta\mathbf x
=
V_0-V.
\]

Take

\[
\mathbf W
=
\operatorname{diag}
(A_1,A_1,A_1,
A_2,A_2,A_2,\ldots).
\]

Then

\[
\boxed{
\delta\mathbf x
=
\frac{V_0-V}
{
\mathbf g_V^T
\mathbf W^{-1}
\mathbf g_V
}
\mathbf W^{-1}\mathbf g_V
}
\]

to first order.

Since

\[
\mathbf g_{V,i}
\sim
A_i\mathbf n_i,
\]

we get approximately

\[
\mathbf W_i^{-1}\mathbf g_{V,i}
\sim
\mathbf n_i.
\]

Thus the correction behaves much more like a uniform normal offset and is substantially less sensitive to mesh density.

This is preferable for adaptive meshes.

---

# 11. Prefer simultaneous area and volume projection

If area and volume are corrected separately, one correction can partially undo the other.

Define

\[
\mathbf c(P)
=
\begin{bmatrix}
A(P)-A_0\\
V(P)-V_0
\end{bmatrix}.
\]

Let

\[
J
=
\begin{bmatrix}
\nabla_P A^T\\
\nabla_P V^T
\end{bmatrix}.
\]

Solve

\[
\min_{\delta\mathbf x}
\frac12
\delta\mathbf x^T
\mathbf W
\delta\mathbf x
\]

subject to

\[
J\delta\mathbf x=-\mathbf c.
\]

The solution is

\[
\boxed{
\delta\mathbf x
=
-
\mathbf W^{-1}J^T
\left(
J\mathbf W^{-1}J^T
\right)^{-1}
\mathbf c.
}
\]

Because there are only two constraints,

\[
J\mathbf W^{-1}J^T
\]

is only \(2\times2\).

Area and volume are nonlinear, so iterate:

1. compute \(A,V\)
2. compute \(\nabla A,\nabla V\)
3. compute correction
4. update `P`
5. repeat to tolerance

---

# 12. Suggested MATLAB constraint-projector interface

```matlab
[P2, info] = project_surface_constraints(P, M, targets, opts);
```

Example:

```matlab
targets.area   = A0;
targets.volume = V0;

opts.preserve_area   = true;
opts.preserve_volume = true;
opts.tol_area        = 1e-10;
opts.tol_volume      = 1e-10;
opts.max_iter        = 10;
```

Return diagnostics such as:

```matlab
info.area_initial
info.volume_initial
info.area_final
info.volume_final
info.rel_area_error
info.rel_volume_error
info.max_displacement
info.rms_displacement
info.iterations
info.converged
```

Do not assemble a dense `3N x 3N` matrix.

If `Ai` contains lumped vertex areas,

```matlab
winv = 1 ./ Ai;
```

and applying \(W^{-1}\) simply scales each vertex 3-vector by `winv(i)`.

---

# 13. Recommended remeshing architecture

Do not hard-code the simulation around a single backend.

Use:

```matlab
[P2, M2, info] = remesh_surface(P, M, target_h, opts);
```

with

```matlab
opts.backend = "legacy";
```

or

```matlab
opts.backend = "mmgs";
```

Later, if desired:

```matlab
opts.backend = "geometry-central";
opts.backend = "pmp";
opts.backend = "cgal";
```

The production simulation loop should not know backend-specific details.

---

# 14. Suggested remeshing options

Example:

```matlab
opts.backend = "mmgs";

opts.hmin = ...;
opts.hmax = ...;

opts.preserve_area   = true;
opts.preserve_volume = true;

opts.area_target   = A0;
opts.volume_target = V0;

opts.constraint_tol      = 1e-10;
opts.constraint_max_iter = 10;

opts.verbose = false;
```

The input

```matlab
target_h
```

should be an `N x 1` vector of desired physical edge lengths.

Backend wrappers should convert this into the library-specific metric representation internally.

---

# 15. Constructing the adaptive target-size field

The membrane code should determine the physics/geometry-based mesh density.

For example,

\[
h_i=f(H_i,K_i,\ldots).
\]

Simple starting choices include

\[
h_i
=
\operatorname{clamp}
\left[
\frac{C}
{\sqrt{|H_i|+\epsilon}},
h_{\min},
h_{\max}
\right],
\]

or

\[
h_i
=
\operatorname{clamp}
\left[
\frac{h_0}
{1+\alpha |H_i|},
h_{\min},
h_{\max}
\right].
\]

The exact formula is a design choice.

Important:

**Do not feed a noisy raw discrete-curvature field directly to the remesher.**

Consider:

- neighborhood averaging
- scalar-field Laplacian smoothing
- clipping
- explicit gradation limiting

The desired size field should vary smoothly.

MMGS also has metric-gradation controls that may be worth exposing.

---

# 16. Thin-neck considerations

Thin necks and highly curved regions are where naive remeshing is most likely to fail.

Potential failure modes:

- collapse of physically important neck edges
- excessive geometric smoothing
- curvature degradation
- local volume error
- self-intersection
- runaway vertex count if curvature-based `h` gets too small

Always impose

\[
h_{\min}
\]

and

\[
h_{\max}.
\]

Also control how rapidly target size may vary spatially.

Do not allow extreme curvature values to send

\[
h\rightarrow0.
\]

---

# 17. Proposed MMGS MEX wrapper

Preferred interface:

```matlab
[P2, M2, info] = mmgs_remesh_mex(P, M, target_h, options);
```

The MEX layer should:

1. validate `P`
2. validate `M`
3. handle MATLAB 1-based indexing
4. initialize MMGS structures
5. insert vertices
6. insert triangles
7. insert scalar sizing/metric information
8. set remeshing parameters
9. call the MMGS surface remesher
10. extract vertices
11. extract triangles
12. convert indices to MATLAB convention
13. free MMGS resources
14. return status/diagnostics

Avoid file I/O in production simulations.

The command-line executable is useful only for installation smoke tests.

---

# 18. MMGS controls worth exposing

Investigate at least:

- minimum target size
- maximum target size
- geometric/Hausdorff approximation tolerance
- metric gradation
- verbosity
- whether selected entities can be marked required/fixed

Do not expose every MMGS parameter initially.

Keep the MATLAB interface small and stable.

---

# 19. Geometry validation independent of the remesher

The wrapper or surrounding MATLAB code should check:

## Closedness

Every undirected edge belongs to exactly two triangles.

## Orientation

Neighboring triangles are consistently oriented.

Global orientation should give the expected volume sign.

## Degeneracy

Reject or warn on triangles below an area tolerance.

## Invalid data

Reject NaNs and infinities.

These checks should not depend on MMGS.

---

# 20. Field transfer after remeshing

Remeshing changes connectivity and often changes vertex count.

Do not transfer evolving fields by vertex index.

Geometry-derived quantities should generally be recomputed:

- normals
- curvature
- vertex areas
- geometric operators

If there are genuinely advected/material fields, define an interpolation from the old surface to the new one.

---

# 21. Benchmark suite

Build saved benchmark geometries.

## Test A: sphere

Checks:

- indexing
- orientation
- volume/area diagnostics
- basic mesh quality
- uniform target-size behavior

## Test B: spheroid

Checks:

- smoothly varying curvature
- adaptive-size behavior

## Test C: tubular dumbbell

Checks:

- neck refinement
- preservation of narrow high-curvature regions

## Test D: pearled tubule

Checks:

- repeated curvature scales
- neck preservation
- vertex-count control

## Test E: actual problematic simulation frames

These are the most important.

Save geometries immediately before current remeshing failures.

Run the identical input through:

1. legacy remesher
2. MMGS

Initially do this **without constraint correction** to compare raw behavior.

---

# 22. Backend comparison metrics

For every input/backend combination report:

```matlab
N_before
N_after

F_before
F_after

A_before
A_after_raw

V_before
V_after_raw

rel_area_error_raw
rel_volume_error_raw

min_edge
max_edge
mean_edge

min_triangle_angle
max_aspect_ratio

max_constraint_correction
rms_constraint_correction

rel_area_error_corrected
rel_volume_error_corrected
```

If possible:

```matlab
max_surface_distance
rms_surface_distance
```

between the pre-remesh and raw post-remesh surfaces.

---

# 23. Resolution-convergence experiment

For at least one smooth closed shape, generate several increasingly fine meshes representing the same surface.

For each:

1. remesh
2. measure raw \(\Delta A\)
3. measure raw \(\Delta V\)
4. apply constraint projection
5. measure required displacement

A healthy pipeline should generally require smaller geometric correction as resolution improves.

If raw volume error remains large under refinement, investigate remeshing or its tolerances.

If raw error is small but projection remains pathological, investigate the projector.

---

# 24. Projection unit tests independent of remeshing

Start with a good closed mesh and deliberately create a small geometry/volume error.

Then run

```matlab
[Pfix, info] = project_surface_constraints(...);
```

Check:

- target volume recovered
- target area recovered if requested
- convergence
- no triangle inversion
- correction magnitude scales with perturbation amplitude
- behavior is not strongly dependent on local mesh density

Compare:

1. current unweighted projector
2. mass-weighted projector

---

# 25. High-value test: deliberately nonuniform sphere

Construct a sphere with strongly nonuniform mesh density.

Introduce a small volume error.

Apply the current unweighted volume projection.

Measure normal displacement versus local lumped vertex area.

Expected for the unweighted method:

- larger-area/coarser vertices move more

Then apply the mass-weighted projection.

Expected:

- normal displacement becomes much closer to uniform

This is a particularly clean test because the expected continuum correction for a sphere is easy to interpret.

---

# 26. Recommended implementation order

## Phase 1: instrument current remesher

Before changing behavior:

- add raw area/volume diagnostics
- add mesh-quality diagnostics
- save problematic before/after meshes
- identify whether errors originate in remeshing or projection

## Phase 2: implement mass-weighted constraint projection

Do this independently of MMGS.

Prefer simultaneous area + volume projection.

Validate on spheres and nonuniform meshes.

## Phase 3: build MMG outside MATLAB

Compile MMG/MMGS normally.

Run one command-line smoke test.

## Phase 4: create minimal MEX gateway

Initially support

```matlab
[P2,M2] = mmgs_remesh_mex(P,M,h);
```

with sensible defaults.

## Phase 5: add key MMGS options

Expose only parameters required by the simulations.

## Phase 6: side-by-side benchmark

Run old and new remeshers on identical saved membrane geometries.

## Phase 7: production integration

Route both through

```matlab
remesh_surface(...)
```

and keep the legacy backend available until MMGS is thoroughly validated.

---

# 27. Working hypothesis

The current expectation is:

1. **MMGS is worth adopting** because it is a modern, maintained, purpose-built adaptive-remeshing backend and can work naturally with a scalar target mesh-size field.

2. **MMGS should not be relied on for strict membrane volume conservation.**
   Exact area/volume constraints belong in the membrane code.

3. **The existing minimum-displacement projection needs scrutiny.**
   If it uses an unweighted Euclidean vertex norm, it is mesh-density dependent.

4. **A lumped-area / mass-weighted projection is more appropriate**, especially with adaptive refinement.

5. **Area and volume should preferably be projected simultaneously**:

\[
\delta\mathbf x
=
-
\mathbf W^{-1}J^T
\left(
J\mathbf W^{-1}J^T
\right)^{-1}
\mathbf c.
\]

6. The decision to replace the remesher should be based on measurements of raw remeshing error **before** projection.

---

# 28. Concrete Codex tasks

Please implement this incrementally rather than replacing production behavior immediately.

## Task 1: inspect existing code

Find and document:

- current call to `isotropicremeshing`
- current MEX wrapper
- current area calculation
- current volume calculation
- current volume projection
- any later area renormalization/correction
- all fields that must survive remeshing

## Task 2: add diagnostics

Create reusable pre-remesh / raw-post-remesh / post-projection diagnostics.

Do not change numerical behavior yet.

## Task 3: inspect existing projector

Determine exactly what norm/objective it minimizes.

Check whether it is the unweighted projection described above.

## Task 4: implement mass-weighted projector

Make it a separate function first.

Support:

- volume only
- area only
- area + volume

## Task 5: projector unit tests

Include:

- uniform sphere
- deliberately nonuniform sphere
- spheroid
- tubular dumbbell

Compare weighted and unweighted behavior.

## Task 6: install/build MMG

Set up MMG so `libmmgs` can be linked by MATLAB MEX on the current Ubuntu/MATLAB environment.

Prefer a reproducible build script.

## Task 7: implement `mmgs_remesh_mex`

Start with:

```matlab
[P2,M2,status] = mmgs_remesh_mex(P,M,h);
```

Then add an options struct.

## Task 8: implement backend wrapper

Create:

```matlab
[P2,M2,info] = remesh_surface(P,M,h,opts);
```

supporting at least

```matlab
opts.backend = "legacy";
opts.backend = "mmgs";
```

## Task 9: benchmark

Run both backends on saved benchmark shapes.

Produce a comparison table/plot for:

- raw volume error
- raw area error
- vertex count
- triangle quality
- required projection displacement
- final area/volume error
- runtime

## Task 10: switch default only after validation

Do not remove the legacy backend initially.

---

# 29. References

## Current remesher

Christopher Helf's MATLAB/OpenMesh wrapper:

https://github.com/christopherhelf/isotropicremeshing

## MMG / MMGS

Official repository:

https://github.com/MmgTools/mmg

Official API/documentation:

https://mmgtools.github.io/

MMGS is the MMG component for adaptation/optimization of surface triangulations.

## Geometry Central

Surface remeshing:

https://geometry-central.net/surface/algorithms/remeshing/

## PMP Library

Adaptive remeshing:

https://www.pmp-library.org/remeshing.html

Algorithm reference:

https://www.pmp-library.org/group__algorithms.html

## CGAL

Sizing-field concept:

https://doc.cgal.org/latest/Polygon_mesh_processing/classPMPSizingField.html

Curvature-based adaptive sizing field:

https://doc.cgal.org/latest/PMP_Remeshing/classCGAL_1_1Polygon__mesh__processing_1_1Adaptive__sizing__field.html

Example:

https://doc.cgal.org/latest/Polygon_mesh_processing/Polygon_mesh_processing_2isotropic_remeshing_with_sizing_example_8cpp-example.html

---

# 30. Target production pipeline

Conceptually:

```matlab
h = compute_target_edge_length(geo, remesh_params);

[Pnew,Mnew,remesh_info] = remesh_surface(P,M,h,remesh_opts);

[Pnew,constraint_info] = project_surface_constraints( ...
    Pnew,Mnew,targets,constraint_opts);

geo = recompute_geometry(Pnew,Mnew);
```

The desired separation of responsibilities is:

- **MATLAB membrane code:** decides where mesh density is needed
- **remesher backend:** changes topology/connectivity and improves mesh quality
- **constraint projector:** enforces strict area/volume constraints
- **geometry code:** recomputes all mesh-derived quantities afterward
