# MMGS Remeshing Backend

This directory contains the in-memory MEX adapter for the optional MMGS
adaptive surface remesher. The legacy `isoremesh` backend remains installed
and is still the default.

## Build

Requirements:

- CMake
- a C/C++ compiler
- Git
- MATLAB with a configured MEX compiler

From the repository root, run:

```bash
remesh/mmgs/build_mmgs_backend.sh
```

The script clones and pins MMG `v5.8.0`, builds a position-independent static
`libmmgs`, and creates `remesh/mmgs/mmgs_remesh_mex.mexa64`. The vendor source,
build tree, installed library, and platform-specific MEX binary are ignored by
Git. To select a non-default MATLAB executable:

```bash
MATLAB_BIN=/usr/local/MATLAB/R2024b/bin/matlab \
    remesh/mmgs/build_mmgs_backend.sh
```

## Solver Use

One setting enables both the adaptive MMGS backend and the mass-weighted
constraint projector:

```matlab
p.remesh_backend = "mmgs";
```

The legacy path is unchanged:

```matlab
p.remesh_backend = "legacy";
```

The projector can be overridden independently:

```matlab
p.constraint_projection = "mass";   % or "legacy", "none"
```

Useful adaptive sizing controls live in `p.remesh_options`:

```matlab
p.remesh_options.curvature_measure = "max_abs_principal";
p.remesh_options.curvature_weight = 1.0;
p.remesh_options.curvature_power = 1.0;
p.remesh_options.hmin_factor = 0.35;
p.remesh_options.hmax_factor = 2.0;
p.remesh_options.preserve_vertex_budget = true;
```

With `preserve_vertex_budget=true`, setting the base target to the current
mean edge length targets approximately the current vertex count. Relative
targets such as `0.99*mean_edge` and `1.8*mean_edge` still request the expected
refinement or coarsening. Set it to `false` when `h` must be interpreted as an
absolute physical target without count normalization.

MMGS-specific controls include `hausdorff_factor`, `metric_gradation`,
`match_target_complexity`, `complexity_max_passes`, and
`vertex_count_tolerance`. A stronger curvature redistribution can generally be
obtained by increasing `curvature_power`; values near `2` produced a clear
density response on the saved pearled membrane benchmark.

## Diagnostics

Run:

```matlab
test_mmgs_backend
test_project_surface_constraints
run_adaptive_remeshing_diagnostics
```

The outputs are written under
`remesh/data/adaptive_remeshing_diagnostics/` as CSV and MAT files. The full
implementation and benchmark report is in
`literature/membrane_adaptive_remeshing_implementation_report.md`.
