# Force-Surface Solver Notes

This folder contains the MATLAB code for the current semipermeable membrane / evolving Stokes surface solver. The active driver is `fs_multi.m`, usually launched from `fs_batch.ipynb`. Older `ms_*` and `gs_*` drivers have been removed or are no longer the active path.

This README is intended as a working handoff note. If a future Codex session needs to resume solver work, read this file first, then inspect `fs_multi.m`, `Geometry.m`, `stokeslet_SLP_triangle.m`, and `stokeslet_SLP_triangle_matrix.m`.

## Active Files

- `fs_multi.m`: current main solver script. It expects caller-workspace variables `p`, `dir`, `verbose`, and `supress_outputs`.
- `fs_batch.ipynb`: current batch launcher. It builds MATLAB command lines for parameter sweeps and sets parameters such as `dt`, `k`, `Sd`, `Da`, `Gamma`, `gamy`, `chi`, `subdivisions`, and `remesh_size`.
- `fs_plotter.m`: plotter for saved `geo*.mat` frames.
- `Geometry.m`: discrete geometry, areas, normals, curvature, bending force, DEC operators, and evolving-surface viscosity operators.
- `Mesh.m`: half-edge connectivity used by `Geometry`.
- `stokeslet_SLP_triangle.m`: matrix-free triangle quadrature Stokes single-layer potential at surface vertices.
- `stokeslet_SLP_triangle_matrix.m`: explicit sparse matrix version of the same SLP, used in the current coupled Newton experiments.
- `stokeslet_SLP_triangle_setup.m`: topology-only cache for triangle SLP quadrature and self-panel corrections.
- `stokeslet_SLP_field.m`: off-surface SLP evaluator for plotting/field diagnostics.
- `rm_rigid_patched.m`: removes selected rigid motion after a timestep. Current `fs_multi.m` uses `"translation"`.
- `deformation_criterion.m`: remeshing trigger.
- `valence_neighbor_average_smooth.m`, `projected_neighbor_smooth.m`, `tangential_neighbor_smooth.m`: experimental smoothing/remeshing helpers.

## Model And Discrete Unknowns

The current multifield formulation solves for:

```text
P       future surface vertex positions
f       independent surface traction density stored at vertices
lambda  scalar/global tension multiplier, also referred to as gamma in notes
```

Velocity is not an independent unknown:

```matlab
u = (P(:) - P0) / p.dt;
```

The code convention is:

- `f` is traction density for the boundary integral operator.
- `traction_to_nodal(f, geo)` multiplies by `geo.v_area` before force balance.
- Bending force `fb = geo.bending_force(1)` is nodal.
- Tension force is discretized as:

```matlab
twoHn = geo.lap * P;
fs = reshape(lambda * twoHn, size(P));
```

In this codebase, `geo.lap * P` is also the discrete area gradient. Finite-difference checks showed `dA/dP` agrees with `geo.lap * P` more exactly than a separately reconstructed `M * (2H n)` expression. Therefore any gamma/lambda Jacobian block should use `geo.lap * P` for consistency with the force law and area residual.

## Current Timestep Structure

A typical `fs_multi.m` timestep:

1. Build or reload `M`, `P`, `velocity`, `lambda`, and `f`.
2. Construct `geo = Geometry(M, P)`.
3. Build frozen old-geometry operators:
   - `KTK`, `DTD` from `geo.evolving_operators()`
   - `mass0`, `M3`
   - `Hess = 2*(KTK + p.k*DTD) + bending-like biharmonic regularization + small mass regularization`
   - explicit SLP matrix `S = stokeslet_SLP_triangle_matrix(P, M, slp_cache)`
   - normal projection matrix `N`
4. Predict `P = P0 + p.dt * velocity`.
5. Iterate an inner nonlinear loop over residuals.
6. Remove translational rigid motion.
7. Remesh if `deformation_criterion(geo)` triggers.
8. Rescale area to `p.area0`, recompute `fb`, and save `geo%d.mat`.

## Residuals In The Active Code

The force balance residual is local to `fs_multi.m`:

```matlab
f_nodal = traction_to_nodal(f, geo);
fb = geo.bending_force(1);
twoHn = geo.lap * P;
u = reshape((P(:) - P0) / p.dt, size(P));
fv = reshape(-2 * (KTK + p.k * DTD) * u(:), size(P));
fs = reshape(lambda * twoHn, size(P));
f_mem = fv + fb + fs;
b_res = -f_nodal - f_mem;
```

In the inner loop this is multiplied by `p.dt`:

```matlab
b = p.dt * force_balance_residual(...);
```

The BIE residual is:

```matlab
slpout = stokeslet_SLP_triangle(P, M, f, slp_cache);
normal_slip = p.Gamma - dot(f, geo.v_normal, 2);
c_res = -u + u_background + p.Sd * slpout ...
    - p.Sd * p.Da * normal_slip .* geo.v_normal;
```

Current code evaluates the BIE residual on `P_bie = P0` / `geo_bie = Geometry(M, P0)`, while `u` comes from the current iterated `P`. This freezes geometry for the hydrodynamic part during the inner loop.

The scalar area residual is:

```matlab
d = (geo.area - p.area0) / p.area0;
```

The user has switched `eps_b` to a non-normalized residual in current experiments. Check `fs_multi.m` before relying on diagnostics.

## Iteration Schemes In `fs_multi.m`

There are currently two inner-loop solve paths arranged so the old one can be commented out and the new one tested.

### Older P/f Coupled Step

The older coupled step solves only for `[dP; df]`:

```matlab
full_Jac = [
    Hess,              -p.dt * M3;
    -I_state / p.dt,    p.Sd * (S + p.Da * N)
];
```

`lambda` is then updated separately using a constant gradient-ascent-like step:

```matlab
lambda = lambda - o.eta * darea;
```

This is not the active target, but it is retained in comments for comparison.

### Current P/f/lambda Coupled Step

The current experiment adds lambda/gamma to the Newton solve:

```matlab
area_gradient = geo.lap * P;
area_gradient = area_gradient(:);
force_gamma_block = -p.dt * area_gradient;
area_gamma_block = area_gradient.' / p.area0;

full_gamma_Jac = [
    Hess,              -p.dt * M3,             force_gamma_block;
    -I_state / p.dt,    bie_df,                sparse(n_state, 1);
    area_gamma_block,   sparse(1, n_state),    0
];
```

The RHS is:

```matlab
full_gamma_rhs = [
    -b(:);
    -c(:);
    -d
];
```

The solve currently uses row/column equilibration:

```matlab
row_scale = 1 ./ max(sum(abs(full_gamma_Jac), 2), eps);
col_scale = 1 ./ max(sum(abs(full_gamma_Jac), 1).', eps);
Dr = spdiags(sqrt(row_scale), 0, size(full_gamma_Jac, 1), size(full_gamma_Jac, 1));
Dc = spdiags(sqrt(col_scale), 0, size(full_gamma_Jac, 2), size(full_gamma_Jac, 2));
full_gamma_step = Dc * ((Dr * full_gamma_Jac * Dc) \ (Dr * full_gamma_rhs));
```

Two condition numbers are useful:

- raw: `condest(full_gamma_Jac)`
- equilibrated: `condest(Dr * full_gamma_Jac * Dc)`

The equilibrated value is the one relevant to the actual backslash solve.

## Current Convergence Behavior

Recent tests were launched from `fs_batch.ipynb` with approximately:

```text
subdivisions = 10
roughness = 0.0005
remesh_size = 1
T = 40000
Sd = 10
Da = 10
Gamma = 0
gamy = 0
chi = 0.1
dt = 0.002
k initially 100000, then tested lower such as 1000
```

Observed behavior:

- With the full gamma solve, `eps_c` often drops to near machine precision quickly.
- `eps_b` decreases for several inner iterations, then starts increasing.
- The restart logic detects monotone `eps_b` increase and reduces `p.chi`.
- The same turn-around behavior occurs even when `eps_b` is not normalized by `u_rms`.
- Reducing artificial dilational viscosity `k` from `100000` to around `1000` improved raw/equilibrated condition numbers by orders of magnitude and improved convergence behavior.
- Current raw full-gamma conditioning has been around `1e11` to `1e13` in some runs.
- Row/column equilibration has brought the effective solved-system condition estimate down roughly to `1e9` to `1e11` on early timesteps.
- The equilibration condition number appeared nearly unchanged across inner iterations in those tests. If cost becomes an issue, compute `Dr`, `Dc`, and the equilibrated condition estimate once at `j = 0` and reuse them while the frozen matrix is unchanged.

Rough interpretation:

- `condest <= 1e7`: conditioning probably not the main issue.
- `1e8` to `1e10`: acceptable but still worth watching.
- `1e11` to `1e13`: concerning; directions may be contaminated by scaling/null modes.
- `>1e13`: do not trust raw backslash directions without scaling/regularization/block solves.

## Conditioning Experiments Tried

### Area Row Scaling

The scalar area residual and area row are scaled by `p.area0`:

```matlab
d = (geo.area - p.area0) / p.area0;
area_gamma_block = area_gradient.' / p.area0;
```

This is mathematically clean, but because `area0` is order one in current runs it did not visibly change behavior.

### Lambda Unknown Scaling

Tried a blind scaling:

```matlab
lambda_scale = 1 / p.dt;
```

and solved for `dlambda_hat` with conversion back to physical `dlambda`. It had no useful effect and was removed. If revisiting lambda scaling, use a smarter norm-matching scale, not `1 / dt`.

### Row/Column Equilibration

Currently active. This improved the effective condition estimate by a few orders of magnitude:

```matlab
Dr * full_gamma_Jac * Dc
```

Keep this unless it is clearly harmful.

### Block Row Scaling

Tried scaling the force-balance, BIE, and area row blocks by average row 1-norm before row/column equilibration. It performed worse and was reverted. Do not re-add the same version blindly.

## Likely Causes Of Remaining Nonconvergence

### Frozen Jacobian Mismatch

`full_gamma_Jac` is assembled outside the inner loop from old/frozen geometry. During the inner loop, `P`, `geo`, area, `geo.lap * P`, bending force, and tension force change. The matrix is therefore a quasi-Newton/frozen-Hessian approximation, not a true Newton Jacobian.

Important next experiment: update only the cheap gamma-related Jacobian blocks each inner iteration:

```matlab
area_gradient = geo.lap * P;
area_gradient = area_gradient(:);
force_gamma_block = -p.dt * area_gradient;
area_gamma_block = area_gradient.' / p.area0;
```

This avoids rebuilding expensive `Hess`, `S`, `KTK`, and `DTD`, but keeps the area/lambda coupling consistent with the current iterated geometry. This idea was discussed and intentionally deferred.

### No Residual Line Search

A standard Armijo line search on the Rayleighian is not appropriate because the coupled system is saddle-point-like. A line search on a residual merit function is appropriate:

```text
Phi = 1/2 * (||b||^2 + ||c||^2 + |d|^2)
```

or with weights:

```text
Phi = 1/2 * (wb ||b||^2 + wc ||c||^2 + wd d^2)
```

Try `alpha = chi, chi/2, chi/4, ...` and accept when the weighted residual decreases. This idea was also discussed and intentionally deferred.

### Incomplete Jacobian

The current matrix ignores or freezes many derivatives:

- derivative of bending force wrt `P`
- derivative of `geo.lap * P` wrt `P`
- derivative of normals in the permeation term
- derivative of SLP geometry wrt `P`
- derivative of `S` and `N` as geometry changes

Some of this is intentional to keep the solve tractable. But if the frozen solve continues to turn around after conditioning and line search, the direction may simply not be a good Newton direction.

### Artificial Dilational Viscosity `k`

`p.k` multiplies `DTD` and is artificial stabilization for soft area/dilation enforcement. Large values can dominate the `P` block and badly condition the saddle solve. Lowering `k` from `100000` to `1000` improved behavior in tests. Continue treating `k` as a tuning parameter.

## Hydrodynamic Discretization

`stokeslet_SLP_triangle.m` and `stokeslet_SLP_triangle_matrix.m` implement a triangle-based Stokes single-layer potential:

- far field: one centroid sample per triangle
- traction: linearly interpolated from vertex values
- targets: evaluated in vectorized chunks
- self panels: Duffy-transformed correction for triangles incident to the target vertex

This is better than the old vertex-lumped `stokeslet_SLP.m`, but it is not a full high-order boundary-integral implementation. Near-singular interactions from non-incident nearby panels are still treated by the centroid rule in the default path.

`stokeslet_SLP_advanced.m` exists for higher-order / near-panel experiments, but it is not the default path.

## Remeshing And Restart Notes

The isotropic remesher lives in `isoremesh` and depends on OpenMesh. See `readme_addition.md` and `isoremesh/README.md` for MEX build notes.

`fs_multi.m` rebuilds the SLP cache after remeshing:

```matlab
slp_cache = stokeslet_SLP_triangle_setup(M);
```

On remesh, current code maps both `velocity` and `f` from old to new mesh using projection/interpolation. Some restart-remesh behavior has been debugged in `remesher_testing.m`.

The current code also does a final area rescale after each timestep:

```matlab
P = P * sqrt(p.area0 / geo.area);
```

Keep this in mind when interpreting area-residual convergence inside the inner solve.

## Practical Guidance For Future Work

Recommended next steps, in order:

1. Keep row/column equilibration and monitor both raw and equilibrated condition estimates.
2. Tune `p.k` downward until area/dilation behavior is acceptable but conditioning is not catastrophic.
3. Update only gamma/area Jacobian blocks inside the inner loop, leaving expensive blocks frozen.
4. Add residual-norm backtracking line search on `[b; c; d]`, not Rayleighian Armijo.
5. If conditioning remains poor, consider a bordered/Schur complement solve for the scalar lambda constraint.
6. If convergence still fails, build a tiny-mesh diagnostic driver that finite-difference checks the active Jacobian blocks against residuals.

Avoid broad refactors until the solver behavior is understood. `fs_multi.m` is intentionally script-style research code, and many comments mark sections that the user manually comments in/out while testing.
