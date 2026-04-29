# TODO: Multifield Implicit Solver

## General Idea

The current `ms_multi.m` timestep solves implicitly for the future surface positions `P`.
Velocity is not stored as an independent unknown inside the nonlinear solve; it is inferred from

```matlab
u = (P - P0) / p.dt;
```

The multifield approach keeps that implicit-in-`P` structure, but introduces the membrane traction/force field `f` as an independent unknown. The inner solve would minimize residuals over

```text
P, f, lambda
```

where `lambda` is the current scalar/global tension multiplier. The first version should not introduce an independent velocity unknown. Wherever velocity appears, use `(P - P0) / dt`.

The target residuals are

```text
u(P) = (P - P0) / dt

r1(P,f) =
    u(P) - uinf(P) + Sd S_P[f]
    - Sd Da (Gamma + f . n(P)) n(P)

r2(P,f,lambda) =
    f + L_P u(P) - fb(P) - fgamma(P,lambda)

r3(P) =
    A(P) - A0
```

Interpretation:

- `r1` enforces the permeable boundary integral / kinematic boundary condition.
- `r2` enforces the membrane force law / evolving surface Stokes balance.
- `r3` enforces global area incompressibility.

The implicit timestep is preserved because every residual is evaluated on the future geometry `P`.

## Relation To Current `ms_multi.m`

Current implementation is reduced-variable:

```matlab
u = (P - P0) / p.dt;
fb = geo.bending_force(1);
fs = reshape(lambda * twoHn, size(P));
fv = reshape((-2 * (KTK + p.k * DTD) * u(:)), size(P));
f_mem = fb + fs + fv;
slpout = stokeslet_SLP_triangle(P, M, f_mem, slp_cache);
```

So `f` is eliminated by substitution:

```text
f = fb + fgamma + fvisc
```

The multifield version does not eliminate `f`. Instead, it solves for `f` while penalizing/enforcing the residual

```text
f + L u - fb - fgamma = 0.
```

Outer code that can remain mostly unchanged:

- timestep loop over `t`
- construction/reload of geometry
- remeshing trigger and remeshing call
- saving `geo%d.mat`
- rigid motion removal, with care about how `f` is updated after `P` changes

Code that likely changes substantially:

- the inner nonlinear optimization loop
- the objective/lagrangian definition
- the meaning and storage of `f`
- residual diagnostics
- remesh interpolation of any saved independent fields

## On-Paper Next Steps

### 1. Settle Residual Signs

Check the signs in

```text
r1 = u - uinf + Sd S[f] - Sd Da (Gamma + f . n)n
r2 = f + L u - fb - fgamma
```

against the current code convention. In particular, current viscous force is

```matlab
fv = -2 * (KTK + p.k * DTD) * u(:);
```

Decide whether the discrete `L u` in `r2` corresponds to

```text
 2 * (KTK + p.k * DTD) * u
```

or the negative of that. This sign/factor is the most important convention to get right.

### 2. Decide What `f` Represents

The Stokeslet routine currently treats its `f` input as a vertex traction-like field and integrates it over triangles internally:

```matlab
S[f] = stokeslet_SLP_triangle(P, M, f, slp_cache);
```

Decide whether multifield `f` represents:

- traction density, or
- area-weighted nodal force.

The cleanest first choice is probably traction density, because that is what the boundary integral wants. Then verify whether `fb`, `fgamma`, and `L u` from `Geometry.m` / `ms_multi.m` need conversion by vertex areas before entering `r2`.

### 3. Define The Objective And Weights

A first objective could be

```text
J(P,f,lambda) =
    1/2 ||r1||_M^2
  + alpha/2 ||r2||_M^2
  + eta/2 (A(P) - A0)^2
```

where `||.||_M` is an area-weighted vertex norm.

Open choices:

- Should `r1` be weighted by `1/(Da*Sd)`?
- Should `r2` be weighted by a membrane-viscosity or force scale?
- Should area be a penalty term, a multiplier update, or both?
- Should `lambda` be optimized directly or updated by ascent as in the current code?

Bad weights will make the optimizer drive one residual while ignoring the others.

### 4. Keep Scalar Tension First

The first implementation should keep scalar/global tension:

```matlab
fgamma = reshape(lambda * (geo.lap * P), size(P));
r3 = geo.area - p.area0;
```

Do not switch to local tension/inextensibility until the scalar multifield solver works. A local tension field would require a larger saddle-point system and new remeshing/interpolation decisions.

### 5. Make A Continuum-To-Code Dictionary

Write down the exact discrete mapping before coding:

```text
S[f]        -> stokeslet_SLP_triangle(P, M, f, slp_cache)
uinf        -> shearextensionflow(P, gamy)
n           -> geo.v_normal
A           -> geo.area
fb          -> geo.bending_force(1)
fgamma      -> reshape(lambda * (geo.lap * P), size(P))
L u         -> sign/factor TBD using KTK and DTD
area weights -> geo.v_area
```

### 6. Define Success Tests

Before long runs, define short tests:

- sphere, `Gamma = 0`, `gamy = 0`: residual should converge to near zero
- compare one timestep with current `ms_multi` in a regime where both formulations should agree
- monitor `norm(r1)`, `norm(r2)`, and `abs(r3)` separately
- verify area conservation
- verify `Da` and `Sd` trends match the intended nondimensional convention

## Implementation Steps

### 1. Add A Residual Function

Create a helper, either local to `ms_multi.m` at first or in a separate file:

```matlab
function [R, parts] = multifield_residual(z, P0, M, p, lambda0, slp_cache)
```

Initial unknown packing:

```text
z = [P(:); f(:); lambda]
```

Return:

```text
R = weighted concatenation of r1, r2, r3
parts.r1
parts.r2
parts.r3
parts.P
parts.f
parts.lambda
parts.geo
```

Start simple and explicit. Optimize later.

### 2. Build A Residual Norm Objective

Add

```matlab
J = 0.5 * (R' * R);
```

Use area-weighting inside `R`, not outside, so diagnostics and optimizer behavior are clear.

### 3. Prototype The Inner Solve Separately

Before replacing the current loop, write a small driver that:

1. loads one `geo%d.mat`
2. builds `P0`, `M`, `p`, `slp_cache`
3. initializes `f` from current formula `fb + fs + fv`
4. evaluates residuals
5. attempts a few optimization iterations

This avoids disturbing the production timestep while the residual conventions are being debugged.

### 4. Choose First Optimizer

For a first implementation, use the simplest robust method available locally:

- finite-difference gradient / MATLAB optimizer if available
- Gauss-Newton with finite-difference Jacobian for tiny meshes
- hand-coded block descent as a fallback

Do not start with a full analytic Jacobian. First prove the residual formulation is correct on small meshes.

### 5. Integrate Into `ms_multi.m`

Once the standalone residual solve works:

- replace the current inner `while` force-descent loop with the multifield optimization loop
- keep the outer timestep, save, remesh, and resize logic
- save the independent `f` field if it is needed for restarts
- if `f` is saved, interpolate/map it after remeshing similarly to `velocity`

### 6. Add Diagnostics

Print or save:

```text
||r1||
||r2||
|r3|
J
area error
iteration count
```

The current single `eps_f` will not be enough. The multifield formulation is only useful if the residual components can be inspected independently.

### 7. Revisit Performance

The expensive part will be repeated applications of

```matlab
stokeslet_SLP_triangle(P, M, f, slp_cache)
```

and possibly finite-difference Jacobian evaluations. Once correctness is established:

- keep Stokeslet application matrix-free
- derive block Jacobian actions where needed
- consider matrix-free Gauss-Newton / Krylov methods
- only assemble a dense Stokeslet matrix for tiny diagnostic meshes, not production

## Open Questions

- Is `f` definitely traction density, or should current force routines be interpreted as nodal forces?
- What is the exact sign/factor mapping between continuum `L u` and code `KTK + p.k*DTD`?
- Should area be enforced by optimizing `lambda`, by lambda ascent, or by a constrained solve?
- Should the permeation term in `r1` use `Sd*Da` exactly as written, or should this be adjusted to match the newest equation `(36)` convention?
- Should `r1` be projected onto the normal direction, or should the full vector residual be retained to enforce tangential no-slip automatically?

