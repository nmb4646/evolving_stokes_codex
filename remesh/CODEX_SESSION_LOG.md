# Codex Session Log

Last updated: 2026-07-09

This file summarizes the work done in this Codex session on the evolving Stokes vesicle code. It includes implemented changes, diagnostics, analysis conclusions, experiments that were abandoned, and current caveats.

## High-Level Thread

The session centered on making `fs_multi` more robust and diagnosable for membrane simulations with bending, permeability, shear flow, and frequent remeshing. The main technical themes were:

- Newton/KKT solve conditioning and line search.
- Whether the bending Hessian matches the bending force actually used in the residual.
- Global volume conservation for impermeable vesicles.
- Parameter-sweep and plotting scripts for permeability, tilt, and deformation.
- Remeshing as a nonphysical reparameterization step that must preserve geometric invariants.

## Major Files Added

### `remesh/newton_correct_volume.m`

Added a standalone post-remesh correction helper:

```matlab
[P, geo, info] = newton_correct_volume(geo, target_area, target_volume, opts)
```

Purpose:

- Project a remeshed surface back to target area and target signed volume.
- Intended to run immediately after `remeshing`.
- For impermeable `Da == 0` runs, target volume should be the fixed initial volume.
- For permeable runs, target volume should usually be the pre-remesh volume from the current timestep, so remeshing itself adds no extra volume change.

Implementation details:

- Rebuilds `Geometry(M, P)` each Newton correction iteration.
- Uses the discrete area gradient convention already used in `fs_multi`:

  ```matlab
  area_gradient = geo.lap * P;
  ```

- Uses an exact signed polyhedral volume gradient matching `Geometry.signed_volume`.
- Forms a two-constraint correction in the span of the area and volume gradients:

  ```matlab
  modes = [area_gradient(:), volume_gradient(:)];
  J = [area_gradient(:).'; volume_gradient(:).'] * modes;
  ```

- Uses a small line search so the correction does not blow up near nearly dependent area/volume directions.
- Defaults:

  ```matlab
  max_iter = 12;
  tol_area = 1e-12;
  tol_volume = 1e-12;
  damping = 1;
  line_search = true;
  rcond_tol = 1e-10;
  ```

Verification performed:

- `checkcode('newton_correct_volume.m')` passed.
- Perturbed sphere test converged to area/volume roundoff.
- Perturbed spheroid test converged to area/volume roundoff.

## Major Files Modified

### `remesh/fs_multi.m`

#### Preconditioning Toggle

Added/used a simple boolean setting:

```matlab
p.precondition_system = true;
```

When enabled, the full KKT-like system is row/column equilibrated before backslash:

```matlab
row_scale = 1 ./ max(sum(abs(full_gamma_Jac), 2), eps);
col_scale = 1 ./ max(sum(abs(full_gamma_Jac), 1).', eps);
Dr = spdiags(sqrt(row_scale), 0, size(full_gamma_Jac, 1), size(full_gamma_Jac, 1));
Dc = spdiags(sqrt(col_scale), 0, size(full_gamma_Jac, 2), size(full_gamma_Jac, 2));
full_gamma_step = Dc * ((Dr * full_gamma_Jac * Dc) \ (Dr * full_gamma_rhs));
```

Reason:

- Condition numbers were often enormous.
- Equilibration helped but did not by itself solve nonlinear convergence failures.

#### Optional Line Search

Added an optional residual-based line search:

```matlab
p.line_search = true;
p.line_search_tau = 0.5;
p.line_search_max_iter = 12;
p.line_search_c = 0;
p.line_search_wb = 1;
p.line_search_wc = 1;
p.line_search_wd = 1;
```

Behavior:

- Computes a merit function from the force-balance residual, BIE residual, and area residual.
- Tries decreasing `alpha`.
- Accepts a step if the merit decreases.
- If no acceptable step is found but one trial improved the merit, it keeps the best trial.
- Otherwise it sets `alpha = 0` and keeps the current iterate.

Observations from run outputs:

- For a spherical membrane in shear, reducing `dt` and allowing enough iterations gave smooth convergence to near machine precision.
- For some elliptical bending-relaxation cases, line search rapidly shrank `alpha` and then got stuck, which indicated the Newton direction was no longer a descent direction for the nonlinear residual.
- Later diagnostics showed that some apparent volume drift was not from the flow solve but from remeshing.

#### Force/BIE/Area KKT System

The current solve assembles:

```matlab
full_gamma_Jac = [
    Hess,              -p.dt * M3,             force_gamma_block;
    -I_state / p.dt,    bie_df,                sparse(n_state, 1);
    area_gamma_block,   sparse(1, n_state),    0
];
```

The scalar area constraint row remains an area constraint only. There is still no volume row in the Newton system. Volume conservation is handled through the SLP projection and, now, through post-remesh correction.

#### `bending_hessian_mode`

Added/used:

```matlab
p.bending_hessian_mode = "exact";
```

Supported modes:

- `"bih"`, `"approx"`, `"cotan_bih"`: use `0.5 * bih`.
- `"exact"`, `"analytic"`: use `geo.bending_hessian(1)`.

Current helper:

```matlab
function bending_hess = bending_hessian_for_solver(geo, bih, p)
    mode = lower(string(p.bending_hessian_mode));
    if mode == "bih" || mode == "approx" || mode == "cotan_bih"
        bending_hess = 0.5 * bih;
    elseif mode == "exact" || mode == "analytic"
        bending_hess = geo.bending_hessian(1);
    else
        error("Unknown bending_hessian_mode '%s'. Use 'bih' or 'exact'.", char(mode));
    end
end
```

Important caveat:

- Fresh runs default to `"exact"` near the top of `fs_multi`.
- Older resumed states that lack this field can still fall back to `"bih"` in the resume compatibility block unless the caller overrides `p.bending_hessian_mode`.

#### Exact Discrete Bending Hessian Integration

`fs_multi` now uses `Geometry.bending_hessian(1)` when `bending_hessian_mode` is `"exact"` or `"analytic"`.

Background:

- There was concern that the solver Hessian used `bih` while the residual used `geo.bending_force(1)`.
- Replacing the bending force with a simple `bih * P` or `-0.5 * bih * P` was discussed/tried as a diagnostic, but that was not retained.
- A finite-difference exact Hessian idea was considered, but rejected as too slow and unnecessary.
- The analytical discrete Hessian from the `literature/discrete bending hessian` derivation was implemented instead.

#### `Geometry.bending_force(1)` Still Used in Residual

The force-balance residual still computes bending force through:

```matlab
fb = geo.bending_force(1);
```

So the exact Hessian option is intended to match the existing discrete bending energy/force model more closely than the simple bih approximation.

#### Translation Projection Option

Added/used:

```matlab
p.project_dP_translation = false;
```

When enabled, the Newton update direction has its area-weighted translation removed:

```matlab
dP = remove_translation_step(dP, geo.v_area);
```

Conclusion from discussion:

- Removing rigid translation can help with under-gauged directions.
- It is not a complete fix for tangential mesh modes because this solver also resolves in-plane fluid motion, so tangential sliding cannot just be projected away wholesale.

#### SLP Volume-Flux Projection

Added/used:

```matlab
p.conserve_slp_volume = true;
```

Matrix projection:

```matlab
S = project_slp_matrix_zero_flux(S, geo);
```

Residual projection:

```matlab
slpout = project_velocity_zero_flux(slpout, geo);
```

The projection subtracts a normal field chosen so the discrete net flux of the SLP output is zero:

```matlab
u = u - (flux / denom) * geo.v_normal;
```

Important conclusions:

- This preserves the ability to change volume when `Da ~= 0`, because only the SLP output is projected; the Darcy/permeation normal-slip term remains.
- This needed to be reflected in the Jacobian as well, so `project_slp_matrix_zero_flux` is applied to the matrix `S` before `bie_df` is formed.
- A sign issue in this projection was identified during testing; after fixing the sign, convergence improved substantially.

#### Remeshing Volume Correction

Before this session, the remesh path did:

```matlab
P = P * sqrt(p.area0 / geo.area);
```

That preserved area only. It did not preserve volume and could directly change reduced volume.

Current post-remesh behavior:

- Fresh initial remesh preserves pre-remesh area and volume:

  ```matlab
  geo_pre = geo;
  [M, P] = remeshing(...);
  geo = Geometry(M, P);
  [P, geo] = newton_correct_volume(geo, geo_pre.area, geo_pre.volume);
  ```

- Fresh runs save:

  ```matlab
  p.area0 = geo.area;
  p.volume0 = geo.volume;
  ```

- Older resumed states get `p.volume0` if missing:

  ```matlab
  if ~isfield(p, 'volume0')
      p.volume0 = geo.volume;
  end
  ```

- Restart remesh and per-timestep remesh call:

  ```matlab
  [P, geo] = newton_correct_volume(geo, p.area0, remesh_target_volume(geo_pre, p));
  ```

- Volume target policy:

  ```matlab
  function target_volume = remesh_target_volume(geo_pre, p)
      if isfield(p, 'volume0') && p.Da == 0
          target_volume = p.volume0;
      else
          target_volume = geo_pre.volume;
      end
  end
  ```

Interpretation:

- For `Da == 0`, remeshing is projected back to the initial volume.
- For `Da ~= 0`, the physical solve may change volume, but the remesher itself should not add another volume jump.

Current caveat:

- The per-timestep remesh condition currently reads:

  ```matlab
  if 1 %hasRemesher && deformation_criterion(geo)
  ```

  so remeshing is unconditional every timestep. This was intentional for the current tests but should be changed back to a toggle/frequency/criterion if needed.

#### Smoke Checks

Performed:

- `checkcode('fs_multi.m')`: no syntax failure; only style warnings around existing unused variables/functions.
- A cheap `T = 0`, `subdivisions = 3`, `initial_remesh = false` smoke run saved `geo0.mat` with both `p.area0` and `p.volume0`.

### `remesh/Geometry.m`

#### Analytical Exact Bending Hessian

Added/used:

```matlab
function [K, gradE] = bending_hessian(obj, Kb)
```

It computes the exact Cartesian Hessian of:

```matlab
willmore_energy(Kb)
```

using local scalar energy:

```matlab
E_i = Kb * M_i^2 / A_i
```

where:

- `A_i` is the barycentric dual area.
- `M_i` is integrated mean curvature.

The Hessian is assembled from local vertex stencils through:

```matlab
local_bending_energy_derivatives(obj, center, Kb)
```

and symmetrized:

```matlab
K = 0.5 * (K + K.');
```

History:

- A finite-difference exact Hessian was discussed and partially considered, but removed.
- The current implementation is analytical.
- It is more expensive than `0.5 * bih`, and because the bending part depends on the current `P_{n+1}`, it is built inside each nonlinear iteration for the experimental exact mode.

### `remesh/fs_batch.ipynb`

#### Parameter Sweep Support

Updated batch logic so `Da`, `Sd`, and `gamy` can each be scalar/list/array values.

The notebook builds jobs over all combinations:

```python
for Da_value, Sd_value, gamy_value in itertools.product(Da_values, Sd_values, gamy_values):
```

Example expected behavior:

```python
Da = [1, 10]
Sd = [1]
gamy = [1, 2, 3]
```

runs 6 jobs.

#### Worker Count Discussion

Discussed `max_workers`.

Conclusion:

- It controls how many MATLAB processes the batch runner launches at once.
- Practical limit depends on CPU cores, RAM, and each simulation's memory footprint.
- Start conservatively, monitor system load/RAM, and increase until MATLAB processes saturate CPU or memory pressure becomes visible.

Current notebook caveat:

- The notebook still has a single `dt = Sd[0]/gamy[0]` line in the command build. If using parameter sweeps over multiple `Sd` or `gamy`, fresh jobs may use the same `dt` unless the notebook is changed to compute `dt` per combination.
- Resumed jobs preserve saved `p.dt` because `fs_multi` currently does not override `p.dt` on resume.

### `remesh/permeation_test_graphic.m`

Created and iteratively expanded a permeability plotting script.

Current capabilities:

- Takes one or multiple `Sd` values.
- Takes any number of `Da_values`.
- Uses `gamy` to locate batch-output folders.
- Computes volume series from `geo*.mat`.
- Plots relative volume change:

  ```matlab
  (V - V0) / V0
  ```

- Optional y-axis log scaling for the volume-over-time plot.
- Optional x-axis scaling.
- Optional original time plot:

  ```matlab
  visual.show_time_plot = false;
  ```

- Computes volume-change rates in several ways:

  ```matlab
  visual.rate_source = "mean_relative_dVdt";
  ```

  Supported values include:

  - `"mean_relative_dVdt"`
  - `"net_relative_dVdt"`
  - `"mean_dVdt"`
  - `"net_dVdt"`

- Rate denominator can use physical time or frame count:

  ```matlab
  rate_time_basis = "time";
  rate_time_basis = "timestep";
  ```

- Time range can be capped:

  ```matlab
  usealltimes = true;
  maxtimestep = 10;
  ```

  If `usealltimes` is false, all rate calculations use only frames through `geo<maxtimestep>.mat`.

- Adds a log-log rate scatter plot:

  ```matlab
  Da vs volume-change rate
  ```

- Fits a power law:

  ```matlab
  Vdot = prefactor * Da^alpha
  ```

  by fitting in log space.

- Fit point range is set by first/last sorted valid Da index:

  ```matlab
  visual.rate_fit_first_index = 4;
  visual.rate_fit_last_index = 10;
  ```

- Multiple `Sd` values get different scatter colors and fit-line colors.

- Per-`Sd` fit ranges are supported:

  ```matlab
  visual.rate_fit_ranges_by_sd = [4,9; 4,8; 4,7];
  ```

- All visual settings are near the top of the script with comments, including figure size, colors, line width, marker size, fonts, grid, legends, and output files.

### `remesh/permeatio_test_graphic_dx.m`

Created a second permeability plotting script based on `permeation_test_graphic.m`.

Purpose:

- Compare three simulations by mesh spacing `dx` rather than labeling them by `Sd`.

Inputs:

- Exactly three `Sd` values, used only to locate folders.
- `Da_values`.
- `gamy`.

How `dx` is computed:

```matlab
dx = mean(geo.he_length)
```

from each run's `geo0.mat`.

Differences from `permeation_test_graphic.m`:

- Labels rate scatter by `dx`, not `Sd`.
- Does not draw fit lines by default:

  ```matlab
  visual.rate_show_fit = false;
  ```

- Same visual settings structure and time-window/rate-denominator options.

Naming note:

- The file is currently named `permeatio_test_graphic_dx.m`, missing the `n` in `permeation`.

### `remesh/diagnose_da0_volume_leak.m`

Created a diagnostic script for the `Da = 0` volume leak investigation.

Purpose:

- Load a run folder.
- Compute saved volume over time.
- Compare finite-difference `dV/dt` against surface flux diagnostics.
- Recompute BIE-related terms from saved geometry and force data.

Important conclusion from using it:

- For `Da = 0`, the BIE/velocity projection could be close to flux-conserving while the saved geometry volume still drifted.
- This meant the volume leak was not necessarily in the flow solve.
- The later remesher investigation confirmed remeshing was changing the saved polyhedral volume.

### `remesh/tilt_over_time.m`

Created and expanded a tilt/deformation plotting script.

Current capabilities:

- Selects simulations by:

  ```matlab
  Da
  Sd
  gamy
  ```

- `gamy` can be a list. If multiple values are given, curves are plotted together.
- `time_start` selects the first `geo*.mat` frame included:

  ```matlab
  time_start = 2;
  ```

- Can scale time by `gamy`:

  ```matlab
  visual.scale_time_by_gamy = true;
  ```

- Supports x/y axis scale:

  ```matlab
  visual.x_scale = "linear"; % or "log"
  visual.y_scale = "linear"; % or "log"
  ```

- Can plot tilt, deformation index `D`, or both:

  ```matlab
  visual.metrics = "tilt";
  visual.metrics = "D";
  visual.metrics = ["tilt", "D"];
  ```

- Tilt angle units:

  ```matlab
  visual.angle_units = "over_pi"; % or "radians", "degrees"
  ```

- Uses `vesicleTiltDeformation`.

### `remesh/vesicleTiltDeformation.m`

Updated the tilt/deformation metric implementation.

Key addition:

```matlab
axisLengthMethod = "ray_intersection";
```

Supported axis-length methods:

- `"projection"`: old method, axis length from max/min projected coordinate.
- `"ray_intersection"`: Zhao/Shaqfeh-style radius sum from ray intersections along the principal axes.

Default:

```matlab
"ray_intersection"
```

Reason:

- The older projection method effectively fit the deformed surface to an ellipsoid-like width estimate.
- The Shaqfeh/Zhao paper's deformation measure is closer to using physical axis lengths from intersections through the centroid/principal axes.

Behavior:

- Ray-intersection method falls back to projection if a ray miss occurs.
- Output includes the method used.

Interpretation from discussion:

- Your previous deformation index was essentially a projection-width/ellipsoid-style proxy.
- Shaqfeh's method is based on semiaxes/lengths from the vesicle shape, not simply the bounding projected width.

### `remesh/fs_plotter.m`

Updated visualization features.

#### Tilt/Deformation Metrics

- Uses `vesicleTiltDeformation(..., tilt_deformation_axis_length_method)`.
- Default axis length method:

  ```matlab
  tilt_deformation_axis_length_method = "ray_intersection";
  ```

#### Surface/Uniform Color Mode

Added color modes:

- `"surface"`
- `"uniform_surface"`
- `"uniform"`

These can hide mesh edges and show the smooth surface itself.

Current logic:

```matlab
use_uniform_color = any(color_mode == ["uniform", "surface", "uniform_surface"]);
hide_mesh_edges = any(color_mode == ["surface", "uniform_surface"]);
```

#### Lighting/Shading

Added lightweight MATLAB surface lighting options:

```matlab
lighting_options.enabled = true;
lighting_options.face_lighting = 'gouraud';
lighting_options.ambient_strength = 0.35;
lighting_options.diffuse_strength = 0.75;
lighting_options.specular_strength = 0.08;
lighting_options.specular_exponent = 12;
```

Purpose:

- Make the surface shape readable when the mesh is hidden.

#### Tangential Streamlines

Added an option for instantaneous streamlines of tangential velocity on the surface:

```matlab
show_streamlines = false;
streamline_seed_count = 40;
streamline_steps = 60;
streamline_step_size = 0.35;
streamline_surface_offset = 0.02;
streamline_projection_neighbors = 12;
streamline_line_width = 1.75;
```

Approach:

- Project velocity onto the tangent plane.
- Trace short surface paths by stepping and projecting trial points back to the mesh.

Caveat:

- This is a visualization tool, not an exact surface-flow integrator.

### `remesh/rotate_vesicle.m`

Added a utility:

```matlab
[P_rot, M_rot] = rotate_vesicle(P, M, angle_over_pi, axis)
```

Inputs:

- `P`: vertex coordinates.
- `M`: mesh connectivity.
- `angle_over_pi`: angle divided by pi, so `0.25` means `pi/4`.
- `axis`: `"x"`, `"y"`, `"z"`, or a 1x3 vector.

Outputs:

- `P_rot`: rotated vertices.
- `M_rot`: unchanged connectivity.

Use case:

- Start prolate/tilted vesicles closer to their expected steady angle to avoid wasting many timesteps waiting for startup transients.

Verification:

- A simple z-axis rotation sanity check was run.

## Analytical/Debugging Reports and Conclusions

### Iteration Blow-Up and Line Search

Observed behavior:

- Without line search or with too-large `dt`, some simulations had `eps_b`, `eps_c`, and area error explode by many orders of magnitude.
- With smaller `dt`, some shear-flow cases converged smoothly to very small residuals.
- Some relaxation cases still hit a wall even with reduced `dt`, suggesting a modeling/discretization issue rather than only a timestep issue.

Main conclusions:

- The Newton direction can stop being a descent direction because the Hessian/Jacobian can mismatch the true nonlinear residual, especially around bending, area constraint, and tangential/remeshing modes.
- A residual-based line search can prevent catastrophic blow-up but cannot fix a fundamentally poor linearization or a nonphysical postprocessing step.

### Finite-Difference Jacobian Check

Discussed as a diagnostic:

- For a residual `R(x)`, compare:

  ```matlab
  J * dx
  ```

  against:

  ```matlab
  (R(x + eps * dx) - R(x)) / eps
  ```

Purpose:

- Detect sign errors.
- Detect missing `dt` or area factors.
- Detect force/Hessian mismatch.
- Detect whether the analytical KKT blocks correspond to the residual actually being minimized/solved.

Status:

- This was discussed and recommended as a diagnostic, but not added as a permanent script in this session.

### `eps_b` and `dt`

Discussed whether `eps_b` depends on `dt`.

Conclusion:

- The force-balance residual currently contains velocity through:

  ```matlab
  u = (P - P0) / p.dt;
  fv = -2 * (KTK + p.k * DTD) * u;
  ```

- Therefore, for the same velocity field, the residual is intended to be comparable across `dt`.
- But for the same displacement `P - P0`, the residual changes with `dt`.
- The printed `eps_b` is not explicitly normalized by `dt` in the current code.

### Null Modes and Under-Gauging

Discussed whether to subtract null solutions from:

```matlab
H dP = R
```

Conclusions:

- Rigid translations and rotations can pollute an ill-conditioned solve.
- Translation projection was added as an option.
- It did not solve the deeper issue for all cases because tangential mesh sliding is physically meaningful in this formulation.
- Full removal of tangential directions is not appropriate if the solver is resolving in-plane fluid motion.

### Effect of Large `k`

Observed:

- Lowering `k` from 1000 to 100 made some simulations work much better.

Interpretation:

- High `k` behaves like a stiff penalty enforcing near-inextensibility/dilation suppression.
- Too high `k` can overconstrain the discrete problem and make the KKT system harder to solve.
- Too low `k` risks allowing unphysical dilational deformation.

Conclusion:

- `k` stiffness can be one cause of failure, but later tests showed failures could persist even with much lower `k`, so it was not the only problem.

### Local Constraint vs Global Penalty

Discussed whether local inextensibility would be "more overconstrained."

Conclusion:

- A local constraint has more constraint equations but, if formulated correctly with surface tension as a Lagrange multiplier field, it is a different saddle-point formulation rather than just "larger k."
- It can be physically more faithful, but numerically more demanding and not automatically easier.

### Bending Force/Hessian Mismatch

Concerns discussed:

- The residual used `geo.bending_force(1)`.
- Approximate Hessian used `0.5 * bih`.
- Replacing the force with `bih * P` or `-0.5 * bih * P` may introduce missing area factors, sign issues, and `dt` scaling confusion.

Conclusion:

- Keep the exact discrete bending force.
- Match it with the exact discrete analytical Hessian as an option.
- Avoid finite-difference Hessian for production use.

### Permeability Interpretation

Discussed expectations:

- Higher `Da` means easier normal permeation.
- In a relaxing spheroid, one might expect higher `Da` to allow faster volume adjustment.

Important correction:

- Sign conventions and the meaning of normal slip matter.
- The volume behavior must be interpreted through:

  ```matlab
  normal_slip = p.Gamma + dot(f, geo.v_normal, 2);
  - p.Sd * p.Da * normal_slip .* geo.v_normal
  ```

- The `Da = 0` limit should be impermeable only if the non-Darcy velocity terms and remeshing are volume-conservative.

### Da = 0 Volume Leak

Initial observation:

- `Da = 0` runs changed volume by similar amounts across mesh size and timestep.
- This was suspicious because the physical normal permeation term was off.

Diagnostic conclusion:

- The no-remesher run was volume conserving.
- With remeshing on, volume drift returned.
- Therefore the dominant leak was the remesher/post-remesh geometry update, not the BIE solve.

### Smooth Theory of SLP Volume Flux

Discussed whether smooth single-layer potential should have zero net flux.

Conclusion:

- In smooth incompressible Stokes theory, the velocity field should be divergence-free, so the net flux through a closed surface should vanish.
- Discrete SLP quadrature/interpolation can violate that identity.
- The projection option is a discretization correction enforcing the closed-surface flux identity at the discrete level.

### Mean-Flux SLP Projection

Approach:

- Compute the discrete net flux of `S[f]`.
- Subtract a normal field whose integrated flux exactly cancels that number.

This is nonunique:

- Many fields could remove net flux.
- The implemented choice subtracts a uniform-in-coefficient vertex-normal field:

  ```matlab
  u = u - (flux / denom) * geo.v_normal;
  ```

Reason for this choice:

- It is simple.
- It is local to the output velocity.
- It is consistent with the existing vertex-area flux quadrature.
- It minimally changes the velocity in a weighted normal-field sense compared with more complex corrections.

### Remesher Investigation

Read files:

- `remesh/isoremesh/README.md`
- `remesh/isoremesh/remeshing.c`
- `remesh/isoremesh/src/IsotropicRemesher.cpp`
- `remesh/isoremesh/src/IsotropicRemesher.h`

Remesher algorithm:

```cpp
splitLongEdges(mesh, high);
collapseShortEdges(mesh, low, high);
equalizeValences(mesh);
tangentialRelaxation(mesh);
projectToSurface(mesh, meshCopy, triangle_bsp);
```

Important observations:

- It is an isotropic mesh-quality remesher.
- It is not designed to preserve enclosed polyhedral volume.
- It is not designed to preserve area either.
- `projectToSurface` projects vertices to nearest points on the old triangle mesh, but new connectivity and new triangles can enclose a different volume.
- Tangential relaxation and edge collapses can change the discrete surface geometry.

Conclusion:

- Remeshing is not purely reparameterization in the current code.
- It must be followed by area/volume projection if the physics requires invariants.

## Current Known Caveats

### Remeshing Is Still Unconditional

Current `fs_multi` still remeshes every timestep because of:

```matlab
if 1 %hasRemesher && deformation_criterion(geo)
```

This may be what is desired for the current experiments, but it is not a clean final control interface.

Possible future fix:

- Add fields such as:

  ```matlab
  p.remesh_each_step = true;
  p.remesh_frequency = 1;
  p.use_deformation_criterion = false;
  ```

### `newton_correct_volume` Changes Geometry

The correction is intentionally small, but it is still a physical-space vertex displacement after remeshing.

This is better than allowing arbitrary volume leak, but the correction can:

- Alter local curvature slightly.
- Interact with very coarse meshes.
- Add a tiny nonphysical normal correction if remeshing errors are large.

Recommendation:

- Log correction magnitudes through `info.step_norm` or final relative errors during testing.

### Velocity After Remeshing

After remeshing, `velocity` and `f` are mapped by:

```matlab
[velocity, f] = map_data(geo, geo_pre, velocity, f);
```

Caveat:

- The saved/interpolated velocity after remesh is not exactly the finite difference between old saved mesh and new saved mesh.
- This can make flux diagnostics using saved `velocity` disagree with finite-difference volume changes.

Possible future fix:

- After correction and mapping, optionally re-zero the mapped velocity flux for `Da == 0`.
- Or recompute a predictor velocity based on the corrected geometry if that is more consistent with the next solve.

### Batch `dt` for Sweeps

Current `fs_batch.ipynb` can sweep `Da`, `Sd`, and `gamy`, but the visible `dt` setup is still scalar:

```python
dt = Sd[0]/gamy[0]
```

For fresh parameter sweeps, this should be changed if each `(Sd, gamy)` pair needs its own timestep.

### `permeatio_test_graphic_dx.m` Filename

The dx plotting script has a typo in the filename:

```text
permeatio_test_graphic_dx.m
```

The intended spelling is probably:

```text
permeation_test_graphic_dx.m
```

It was left as-is to avoid breaking references.

### Exact Bending Hessian Cost

`geo.bending_hessian(1)` is much more expensive than the bih approximation.

Reason:

- It assembles local Hessian blocks for each vertex each nonlinear iteration.

Current state:

- Fresh `fs_multi` defaults to exact mode.
- This may be appropriate for debugging correctness, but performance should be monitored.

## Useful Current Settings and Entry Points

### Run Solver

Main solver:

```matlab
fs_multi
```

Important parameters:

```matlab
p.precondition_system
p.line_search
p.bending_hessian_mode
p.conserve_slp_volume
p.project_dP_translation
p.Da
p.Sd
p.gamy
p.dt
p.k
p.volume0
```

### Plot Volume/Permeability

```matlab
permeation_test_graphic
```

or dx comparison:

```matlab
permeatio_test_graphic_dx
```

### Plot Tilt/Deformation

```matlab
tilt_over_time
```

### Plot Surface

```matlab
fs_plotter
```

Useful display options:

```matlab
color_mode = "surface";
lighting_options.enabled = true;
show_streamlines = true;
tilt_deformation_axis_length_method = "ray_intersection";
```

### Rotate Initial Geometry

```matlab
[P, M] = rotate_vesicle(P, M, 0.3, "y");
```

### Correct Post-Remesh Area/Volume Manually

```matlab
geo_pre = geo;
[M, P] = remeshing(int32(M), P, int32([]), r.edge_length, int32(r.n_iter));
M = cast(M, "double");
geo = Geometry(M, P);
[P, geo, info] = newton_correct_volume(geo, p.area0, p.volume0);
```

For permeable runs where volume may physically change:

```matlab
[P, geo, info] = newton_correct_volume(geo, p.area0, geo_pre.volume);
```

## Files Touched or Created During the Session

Created:

- `remesh/newton_correct_volume.m`
- `remesh/rotate_vesicle.m`
- `remesh/diagnose_da0_volume_leak.m`
- `remesh/permeation_test_graphic.m`
- `remesh/permeatio_test_graphic_dx.m`
- `remesh/tilt_over_time.m`
- `remesh/CODEX_SESSION_LOG.md`

Modified:

- `remesh/fs_multi.m`
- `remesh/Geometry.m`
- `remesh/fs_batch.ipynb`
- `remesh/fs_plotter.m`
- `remesh/vesicleTiltDeformation.m`

Generated or present as local artifacts:

- `.asv` MATLAB autosave files exist for some edited scripts.
- Some `.mat` geometry/data files were present or generated during experimentation.
- A 2011 paper PDF exists under `literature/`.

## Recommended Next Steps

1. Run a controlled `Da = 0` case with remeshing on and check whether `p.volume0` is preserved after the new correction.
2. Log `newton_correct_volume` relative area/volume errors and step norms during a short run.
3. Change unconditional remeshing to an explicit parameter once the current tests are done.
4. Decide whether fresh batch sweeps should compute `dt` per `(Sd, gamy)` pair.
5. Benchmark exact bending Hessian versus bih mode on a representative mesh.
6. If mapped velocity flux still causes predictor issues, re-zero mapped velocity flux after remeshing for `Da == 0`.
