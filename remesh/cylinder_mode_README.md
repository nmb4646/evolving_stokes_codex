# Cylinder Mode Analysis

This subsystem analyzes the cylindrical portion of saved `geoN.mat` membrane
series without modifying the solver or any existing analysis code.

## Run an analysis

Edit the parameter block at the top of:

```matlab
cylinder_mode_analysis.m
```

Then run the script in MATLAB. Set `Sd_values`, `Da_values`, `gamy_values`,
and `v_values` to scalars or arrays. The script analyzes the Cartesian
product of those arrays and constructs the standard `fs_batch_data` folder
names automatically.

The default example analyzes:

```text
Sd_1.00em06_Da_0.00ep00_gamy_p3.70em05_v_3.50em01
```

Use `m_max = 0` for the fast axisymmetric path. Increase `m_max` to retain
bending, elliptical, helical, and other nonaxisymmetric modes.

## Validate the pipeline

Run:

```matlab
validation = cylinder_mode_synthetic_tests();
```

This creates capped-cylinder time series with known behavior and validates:

- a static translated and rotated cylinder;
- prescribed axisymmetric exponential growth;
- prescribed helical exponential decay;
- separation of opposite helical handedness;
- changing mesh connectivity;
- cylindrical-core stability under remeshing;
- paired `+epsilon/-epsilon` growth;
- cancellation of a common mode in paired runs.

## Sweep Sd at fixed capillary number

Edit and run:

```matlab
cylinder_mode_growth_vs_sd.m
```

The script pairs every requested `Sd` with `gamy = Ca*Sd`, runs the existing
mode decomposition, and plots each accepted `(m,n)` growth rate against `Sd`.
It writes the combined fit table, sweep summary, figure, and MAT file under
`remesh/data/cylinder_mode_growth_vs_sd`.

The validation table and figure are written under:

```text
remesh/data/cylinder_mode_synthetic_validation
```

## Outputs

Each analyzed series receives its own folder under
`remesh/data/cylinder_mode_analysis` by default:

```text
cylinder_frame_metrics.csv
cylinder_mode_coefficients.csv
cylinder_growth_rates.csv
cylinder_mode_analysis.mat
cylinder_resolved_config.json
cylinder_diagnostics/
```

Diagnostics include:

- area and volume evolution;
- fitted radius and core length;
- core area and axial fractions;
- axis changes;
- cylindricality and reconstruction errors;
- projection conditioning and conjugacy error;
- detected core boundaries;
- modal amplitudes and accepted exponential fits;
- initial/final spectra;
- growth rate versus mean `qR`.

Frames rejected from growth fitting remain in the output with explicit quality
flags.

## Numerical conventions

The pipeline:

1. numerically sorts `geoN.mat` files;
2. uses explicit saved time when present, otherwise reconstructs physical time
   from frame indices and `p.dt`;
3. computes barycentric vertex area weights;
4. removes global translation and tracks the principal tube axis;
5. parallel-transports the transverse basis to avoid arbitrary frame spin;
6. tracks a persistent material core while connectivity remains unchanged;
7. falls back to a broad radius-plateau core after remeshing;
8. uses a co-stretching axial coordinate and Tukey end window;
9. projects with weighted least squares by default;
10. groups only true conjugates `(m,n)` and `(-m,-n)`.

Opposite helical modes `(m,n)` and `(m,-n)` remain separate.

For paired runs with matching vertices, radial differences are formed on one
averaged cylindrical base and projected once with common weights. This avoids
subtracting coefficients produced by two different discrete projection
operators. If correspondence is unavailable, the code falls back to
independent coefficient projection and records that choice in
`response_method`.

## Interpretation

The Fourier modes are a local approximation on the central region of a finite
capped tube. Reported rates should only be used while:

- the perturbation is small;
- the core domain remains stable;
- the mode is spatially resolved;
- its amplitude is above the estimated projection noise;
- `qR` changes modestly over the fit;
- the semilog amplitude is approximately linear.

Later nonlinear pearling, necking, cap invasion, remeshing jumps, or severe
shape changes are retained for diagnostics but excluded from fitting.
