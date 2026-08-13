# Adaptive Membrane Remeshing and Constraint Projection Report

## Summary

The repository now has a second, optional remeshing path based on MMGS 5.8.0,
while the original `isoremesh` MEX path remains in place and remains the
default. Selecting

```matlab
p.remesh_backend = "mmgs";
```

enables curvature-adaptive isotropic remeshing and, unless explicitly
overridden, the new mass-weighted area/volume projector. Existing runs using
the legacy backend retain the existing legacy projector and remeshing calls.

The main result is not that MMGS is universally better. It provides the
required spatially varying size field and good count control, but it can move
the represented polyhedral surface enough to require a nontrivial projection,
especially on thin tubules. The diagnostics therefore preserve a strict
separation between raw remesher error and projector error.

## Solver Interface

The streamlined production setting is:

```matlab
p.remesh_backend = "mmgs";
```

The following explicit overrides remain available:

```matlab
p.constraint_projection = "mass";    % "legacy", "mass", or "none"
p.remesh_options = struct();
p.constraint_projection_options = struct();
```

`fs_batch.ipynb` exposes:

```python
remesh_backend = "legacy"             # "legacy" or "mmgs"
constraint_projection = "auto"        # "auto", "legacy", or "mass"
```

With `auto`, MMGS selects the mass-weighted projector and legacy selects the
old projector. This preserves old behavior while requiring only one setting
to use the complete new pipeline.

All production remesh sites in `fs_multi.m`, including startup, resume,
stuck-line-search recovery, and the existing disabled periodic-remesh block,
route through the same thin wrappers. Legacy-only valence smoothing is still
applied only to the legacy backend.

## New Remeshing Architecture

### Backend abstraction

`remesh_surface(P,M,target_h,opts)` is the backend-neutral entry point.

- `opts.backend="legacy"` calls the existing `remeshing.mexa64` with one
  scalar target size.
- `opts.backend="mmgs"` passes the full per-vertex target size to MMGS.
- Input and output meshes are checked for finite values, valid indices,
  closure, orientation consistency, and degenerate faces.
- The return structure contains pre-remesh and raw-post-remesh diagnostics.

`remesh_membrane` is the solver-facing wrapper. It computes a curvature-based
size field only for MMGS and maintains the solver's existing `[M,P]` output
order.

### MMGS MEX backend

`mmgs/mmgs_remesh_mex.c` is an in-memory MEX gateway. It performs no mesh file
I/O. It:

1. Validates `P`, `M`, `h`, and the options structure.
2. Initializes MMGS mesh and metric objects.
3. Inserts MATLAB vertices and 1-based triangles directly.
4. Loads the scalar isotropic size metric.
5. Exposes `hmin`, `hmax`, `hausd`, `hgrad`, verbosity, and angle detection.
6. Executes `MMGS_mmgslib`.
7. Extracts the new vertices, faces, and MMGS return status.
8. Frees all MMGS allocations.

The reproducible build script pins:

- MMG tag: `v5.8.0`
- commit: `4d8232c8aebfed877935d75d4d4a67e850962422`
- static, position-independent `libmmgs`
- Scotch disabled

The build and installation command is:

```bash
remesh/mmgs/build_mmgs_backend.sh
```

Generated vendor, build, installation, and MEX files are ignored by Git. The
current workspace has a successfully built `mmgs_remesh_mex.mexa64`.

### Curvature sizing field

`compute_target_edge_length` defaults to the largest absolute principal
curvature. It forms

```text
H = integrated mean curvature / vertex area
K = angle defect / vertex area
k1,2 = H +/- sqrt(max(H^2-K,0))
indicator = max(abs(k1),abs(k2))
```

The indicator is neighborhood-smoothed, percentile-clipped, robustly
normalized, mapped to a target size, clamped to `[hmin,hmax]`, smoothed again,
and explicitly gradation-limited. The current conservative defaults are:

```matlab
curvature_measure = "max_abs_principal"
curvature_weight = 1.0
curvature_power = 1.0
curvature_smoothing_iterations = 12
hmin_factor = 0.35
hmax_factor = 2.0
gradation_ratio = 1.3
```

The smoothing count was chosen using the saved thin-tubule frame. Four passes
left a noisy size field that produced a minimum angle as low as 7.8 degrees in
one trial. Twelve passes produced a 21.2 degree minimum angle with a smaller
post-remesh correction.

### Vertex-budget control

Curvature refinement initially caused an unacceptable 936 to 1793 vertex
increase on the thin-tubule benchmark. Two controls now prevent this.

First, the size field is globally scaled so its metric complexity

```text
integral(1/h^2 dA)
```

matches a relative target derived from the current mesh and requested base
edge length. Consequently:

- `base_h = mean(current edges)` targets approximately the current count.
- `base_h = 0.99*mean(current edges)` requests mild refinement.
- `base_h = 1.01*mean(current edges)` requests mild coarsening.
- `base_h = 1.8*mean(current edges)` still requests strong coarsening.

Second, MMGS is rerun from the unchanged input mesh with a calibrated metric
scale until the output count is within the default 8 percent tolerance of the
predicted count, or six passes are exhausted. The closest result is retained.
This is more expensive than one MMGS call, but representative meshes usually
needed one to three passes, and the physical solver remains much more
expensive than remeshing.

Set

```matlab
p.remesh_options.preserve_vertex_budget = false;
p.remesh_options.match_target_complexity = false;
```

when a supplied `h` must be treated as an absolute physical size with no count
normalization or calibration.

### Geometric tolerance

The default MMGS geometric tolerance is

```matlab
hausd = 0.25 * min(target_h)
```

On the severe saved tubule, increasing this factor improved triangle angles
but increased geometric movement and the required constraint correction.
For example, values near `0.5` gave angles near 26 degrees but corrections
near 0.09 mean edge lengths. The selected `0.25` default gave a 21.2 degree
minimum angle and a 0.043-edge correction.

## New Constraint Projection

### Existing method

`newton_correct_volume` uses the two modes

```matlab
area_gradient = geo.lap * P;
volume_gradient = exact_polyhedral_volume_gradient(M,P);
```

and solves in their unweighted Euclidean coordinate inner product. It is a
local minimum-Euclidean-vertex-displacement correction. Because the volume
gradient at a vertex scales approximately with its represented area, coarse
vertices tend to move farther than refined vertices. The result therefore
depends on mesh density.

### Mass-weighted method

`project_surface_constraints` instead minimizes each linearized correction in
the lumped surface `L2` norm:

```text
min 1/2 deltaP' W deltaP
subject to J deltaP = -c
```

where each vertex's three coordinates receive its lumped barycentric area as
the diagonal weight. The step is

```text
deltaP = -W^-1 J' (J W^-1 J')^-1 c.
```

Only a 1-by-1 or 2-by-2 Schur matrix is formed. No dense `3N` matrix is
assembled. The implementation supports volume only, area only, and
simultaneous area plus volume.

Both area and signed-volume gradients are calculated analytically from the
same polyhedral definitions used by `Geometry`. Constraint rows are scaled by
their target magnitudes before solving. Each nonlinear iteration includes:

- Schur conditioning check and pseudoinverse fallback
- maximum step relative to mean edge length
- triangle inversion and degeneracy rejection
- residual-merit backtracking line search
- full recomputation of geometry and gradients

This is an iterative minimum-norm local projection, not a proof of the global
minimum displacement over the nonlinear constraint manifold.

### Volume-priority fallback

Area and volume gradients become nearly dependent near a sphere. In the MMGS
coarse-sphere tests, the raw remesh had a slightly different discrete reduced
volume from the input. Exact simultaneous recovery therefore required a
slow symmetry-breaking deformation and stalled before the requested
tolerance.

If simultaneous projection reaches its iteration limit, the projector now
runs a final mass-weighted volume-only correction. It reports
`info.converged=false`, records that the fallback was used, and makes volume
accurate. `correct_remesh_constraints` also warns when full area-plus-volume
convergence was not achieved. This is intentional: for an impermeable
simulation, exact volume is safer than silently retaining both small area and
volume errors.

## Projection Diagnostics

`test_project_surface_constraints` passed the following cases:

| Case | Weighting | Iterations | Final relative area error | Final relative volume error | Maximum displacement |
|---|---:|---:|---:|---:|---:|
| Uniform sphere, volume only | mass | 3 | not constrained | 0 | 1.001e-2 |
| Spheroid, area and volume | mass | 3 | 9.09e-16 | 1.03e-15 | 1.130e-2 |
| Tubular dumbbell, area and volume | mass | 3 | 4.65e-12 | 6.15e-12 | 8.794e-3 |

The deliberately nonuniform sphere had a maximum/minimum lumped vertex-area
ratio of 217.3. After a one-percent radial volume perturbation:

| Weighting | Normal-displacement coefficient of variation | Correlation with vertex area | Final relative volume error |
|---|---:|---:|---:|
| Mass | 5.58e-3 | -0.9999 | 4.37e-16 |
| Euclidean | 1.834 | 0.9999 | 4.37e-16 |

The signed correlation for the mass result reflects a very small residual
trend; the coefficient of variation is the important measure. The weighted
correction is effectively uniform, while the Euclidean correction varies by
order one and tracks local represented area almost exactly.

For the dumbbell test, the Euclidean projector needed a maximum displacement
of `3.788e-2`, compared with `8.794e-3` for mass weighting.

## Controlled MMGS Diagnostic

`test_mmgs_backend` supplies a known hemispherical target field, with smaller
sizes in the north and larger sizes in the south. It produced:

| Metric | Value |
|---|---:|
| Input vertices | 1002 |
| Predicted target vertices | 784 |
| Output vertices | 820 |
| Count error | 4.59% |
| MMGS passes | 5 |
| Median edge in fine region | 0.0937 |
| Median edge in coarse region | 0.2057 |
| Fine/coarse edge ratio | 0.455 |
| Minimum triangle angle | 32.1 degrees |
| Maximum aspect ratio | 1.66 |
| Closed and consistently oriented | yes |

This validates the MEX metric convention independently of curvature
estimation: smaller supplied scalar sizes do produce denser output.

## Backend Benchmark

`run_adaptive_remeshing_diagnostics` compares legacy and MMGS on spheres,
a spheroid, a synthetic tubular dumbbell, and two saved simulation frames.
The complete CSV and MAT outputs are in:

```text
remesh/data/adaptive_remeshing_diagnostics/
```

### Saved thin tubule

Input:

```text
Sd_1.00em06_Da_0.00ep00_gamy_p3.70em06_v_3.50em01/geo132.mat
936 vertices, 1878 faces
```

| Metric | Legacy | MMGS |
|---|---:|---:|
| Output vertices | 922 | 896 |
| Runtime | 0.033 s | 0.082 s |
| Raw relative area error | 2.88e-3 | 1.93e-3 |
| Raw relative volume error | 4.90e-3 | 1.10e-3 |
| Raw minimum angle | 29.5 deg | 21.2 deg |
| Raw maximum aspect ratio | 1.48 | 2.26 |
| Corrected relative area error | 5.65e-16 | 1.41e-16 |
| Corrected relative volume error | 3.03e-16 | 0 |
| Maximum correction / mean edge | 5.23e-3 | 4.27e-2 |
| Projection curvature RMS change | 3.56e-3 | 3.18e-2 |

MMGS has lower raw area and volume error here, and its final count is close to
target, but it is not uniformly better: its triangle quality and correction
size are worse than legacy on this already-good input. This frame should be
used as a continuing regression test before enabling MMGS for every timestep.

### Saved pearled membrane

Input:

```text
Sd_1.00em06_Da_0.00ep00_gamy_p8.00em07_v_8.63em01/geo200.mat
1196 vertices, 2388 faces
```

| Metric | Legacy | MMGS |
|---|---:|---:|
| Output vertices | 1193 | 1210 |
| Runtime | 0.043 s | 0.153 s |
| Raw relative area error | 1.27e-4 | 9.04e-5 |
| Raw relative volume error | 1.95e-4 | 1.22e-4 |
| Raw minimum angle | 41.9 deg | 33.4 deg |
| Corrected relative volume error | 3.79e-16 | 5.69e-16 |
| Maximum correction / mean edge | 7.47e-4 | 9.44e-4 |

Both pipelines are benign on this frame. MMGS preserves the vertex budget,
has slightly lower raw integral errors, and requires a very small correction.

### Smooth-shape resolution trend

For uniform MMGS remeshing of increasingly fine sphere meshes, raw relative
volume error decreased from `7.83e-3` at 92 vertices to `8.29e-4` at 252 and
`3.29e-4` at 642. Required maximum correction relative to mean edge decreased
from `4.84e-2` to `2.47e-2` and `1.05e-2`. This is the expected healthy
resolution trend.

The near-sphere simultaneous area/volume solve invoked the volume-priority
fallback. Final volume was between `5.4e-12` and `4.3e-15` relative error, but
area residuals remained between `2.0e-4` and `2.3e-5`. Spheroid and dumbbell
simultaneous projections converged without fallback.

## Diagnostics and Limitations

`mesh_surface_diagnostics` reports counts, topology, closure, orientation,
degeneracy, Euler characteristic, area, signed volume, edge statistics,
triangle angles, aspect ratio, vertex-area ratio, and optional reference
distances.

The reported symmetric surface distance is a nearest-vertex approximation,
not an exact point-to-triangle Hausdorff distance. It is useful for detecting
large changes but is sensitive to sampling redistribution and should not be
interpreted as a sharp geometric bound.

Self-intersection detection is not yet implemented. MMGS output is validated
for closedness, orientation, and degeneracy, but a strongly folded input could
still require a separate intersection library.

The curvature indicator is discrete and can be noisy. The current conservative
default deliberately smooths it strongly. On the pearled saved frame, the
default field was mild enough that local density redistribution was weak.
Increasing `curvature_power` to `2` made the response clearer; a test with
`curvature_weight=2, curvature_power=2` gave a high-curvature/low-curvature
median-edge ratio near `0.77`. This stronger setting should be tested on the
specific dynamics before becoming a global default.

No evolving material scalar field is automatically transferred. Geometry,
operators, curvature, and normals are recomputed. The solver's existing
nearest-surface `map_data` process remains responsible for velocity and
traction transfer after remeshing.

## Files Added

- `remesh_surface.m`: backend-neutral remeshing API and raw diagnostics
- `remesh_membrane.m`: solver-facing adaptive wrapper
- `compute_target_edge_length.m`: curvature sizing and count normalization
- `mesh_surface_diagnostics.m`: reusable topology/quality diagnostics
- `project_surface_constraints.m`: mass-weighted nonlinear projector
- `correct_remesh_constraints.m`: legacy/mass/off projector selector
- `remeshing_backend_available.m`: backend availability check
- `mmgs/mmgs_remesh_mex.c`: in-memory MMGS MEX gateway
- `mmgs/build_mmgs_backend.sh`: pinned reproducible build
- `mmgs/README.md`: build and use instructions
- `test_project_surface_constraints.m`: projector regression suite
- `test_mmgs_backend.m`: controlled adaptive-size regression
- `run_adaptive_remeshing_diagnostics.m`: full side-by-side benchmark

The original `isoremesh` directory and the earlier MATLAB
`remeshing_curvature.m` prototype were not modified.

## Verification Performed

- MMG 5.8.0 static library and MEX built successfully.
- MMGS command-line/native smoke test completed during installation.
- Controlled MEX size-field test passed.
- Projector regression suite passed in MATLAB R2024a.
- Legacy/MMGS benchmark completed in MATLAB R2024b with both MEX runtimes.
- All new MATLAB files passed `checkcode` with zero messages.
- `fs_multi` zero-timestep production smoke test completed using only
  `p.remesh_backend="mmgs"`; it automatically selected `mass`, remeshed the
  initial surface, rebuilt geometry/operators, saved `geo0.mat`, and returned
  a valid closed mesh with internally exact saved area and volume targets.
- `fs_batch.ipynb` passed JSON validation.
- `build_mmgs_backend.sh` passed `bash -n`.
- `git diff --check` passed.

## Recommendation

Keep `legacy` as the default while running paired short simulations with MMGS
on the actual thin-neck and pearling protocols. Enable MMGS first where
spatial adaptation is genuinely needed, monitor raw area/volume error and
`projection_max_displacement_over_edge`, and treat values above roughly 0.05
as a reason to inspect that frame rather than blindly continuing.

For stronger pearling refinement, start with:

```matlab
p.remesh_backend = "mmgs";
p.remesh_options.curvature_weight = 2;
p.remesh_options.curvature_power = 2;
```

and compare triangle quality, Willmore-energy jumps, and post-projection
displacement against the conservative defaults before using it in long runs.
