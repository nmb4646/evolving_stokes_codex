# Cylindrical Mode Analysis Pipeline for Membrane-Tubule Simulations

## Purpose

Build a reproducible pipeline that scans folders containing time-series membrane geometries, decomposes the cylindrical portion of each geometry into axial and azimuthal cylindrical modes, and estimates the short-time growth rate of each mode.

The pipeline must be robust to:

- finite rounded end caps;
- global translation and rotation;
- axial stretching of the tube;
- gradual change of the mean tube radius;
- irregular or changing surface meshes;
- small numerical noise;
- nonuniform output times;
- multiple simulations and parameter sets;
- optional unperturbed, `+epsilon`, and `-epsilon` companion runs.

The primary scientific output is a table of modal growth rates

\[
\sigma_{mn} \equiv \frac{d}{dt}\log A_{mn},
\]

reported together with the corresponding instantaneous or window-averaged dimensionless wavenumber

\[
q_n R(t).
\]

The analysis should focus on the early interval in which perturbations are small and exponential growth or decay is a reasonable approximation.

---

## High-level recommendation

Implement two analysis paths that share the same preprocessing and growth-rate code.

1. **Axisymmetric fast path (`m_max = 0`)**
   - Best for pearling and axisymmetric simulations.
   - Average the radial displacement around each cross-section and perform a one-dimensional axial decomposition.
   - This should be the default when only axisymmetric modes are scientifically relevant.

2. **General three-dimensional path (`m_max > 0`)**
   - Decompose the radial displacement into cylindrical harmonics
     \[
     e^{i(m\phi + qx)}.
     \]
   - Use this for bending, elliptical, helical, wrinkling, or symmetry-breaking modes.

The implementation should not assume that all meshes have identical connectivity, although it should use a faster material-tracking path when connectivity and vertex IDs remain consistent.

---

# 1. Expected command-line behavior

Provide a command-line entry point resembling:

```bash
analyze_cylindrical_modes \
    --root <simulation-root> \
    --config mode_analysis.yaml
```

The exact command and paths should follow the conventions of the repository.

The program should:

1. discover simulation folders;
2. identify and sort geometry frames by physical time;
3. analyze every frame;
4. save per-frame geometric quantities and modal coefficients;
5. determine valid early-time fitting windows;
6. fit growth rates;
7. generate diagnostic plots and machine-readable summaries;
8. continue processing other series if one series fails;
9. write clear logs and failure reasons.

Support processing one series, all series, or a user-selected subset.

---

# 2. Input abstraction

## 2.1 Geometry frame

Each frame must expose:

```text
time
P: N_vertices x 3 vertex coordinates
M: N_faces x 3 triangle connectivity
optional vertex IDs
optional simulation metadata
```

The I/O layer should be modular. Reuse existing project readers wherever possible.

Possible input formats may include:

- MATLAB `.mat`;
- VTK or VTU;
- PLY;
- OBJ;
- NPZ;
- HDF5;
- repository-specific formats.

Do not hard-code one format into the analysis logic.

## 2.2 Series metadata

Read simulation parameters from an existing manifest when available. Useful metadata include:

```text
reduced volume
flow strength or capillary number
bending modulus
viscosity ratio
time scale
initial perturbation amplitudes
seeded mode numbers
run type: baseline, plus, minus, ordinary
matching series ID
```

If no manifest exists, infer only what is safe to infer from folder or file names and record missing values as null rather than guessing.

---

# 3. Output structure

For each simulation series, create an analysis folder similar to:

```text
analysis/
    resolved_config.yaml
    frame_metrics.csv
    mode_coefficients.parquet
    growth_rates.csv
    summary.json
    diagnostics/
        geometry_metrics.png
        mode_amplitudes.png
        mode_growth_fits/
        spectra/
        reconstruction_error.png
    logs/
        analysis.log
```

Use formats already preferred by the repository. Parquet is recommended for long coefficient tables, with CSV as an optional export.

## 3.1 `frame_metrics`

One row per frame. Include at least:

```text
series_id
frame_index
time
n_vertices
n_faces
surface_area
volume
reduced_volume
axis_x
axis_y
axis_z
center_x
center_y
center_z
core_x_min
core_x_max
core_length
mean_core_radius
radius_std
cylindricality_error
core_area_fraction
windowed_area
reconstruction_error
axis_angle_change
quality_flags
```

## 3.2 `mode_coefficients`

Use long format, one row per `(time, m, n)` group:

```text
series_id
frame_index
time
m
n
basis_type
physical_wavenumber
dimensionless_wavenumber_qR
coefficient_real
coefficient_imag
grouped_amplitude
grouped_power
phase
noise_estimate
quality_flags
```

Store raw complex coefficients whenever possible. Constant normalization factors do not affect the fitted growth rate, but the convention must be documented and consistent.

## 3.3 `growth_rates`

One row per mode and selected fitting window:

```text
series_id
m
n
fit_start_time
fit_end_time
n_fit_points
growth_rate_sigma
growth_rate_standard_error
phase_rate_omega
r_squared
amplitude_start
amplitude_end
amplitude_ratio
mean_qR
min_qR
max_qR
mean_radius
mean_core_length
fit_method
fit_status
quality_flags
```

---

# 4. Coordinate conventions

The longitudinal tube axis is called the local \(x\)-axis.

After rigid alignment, each surface point is represented as

\[
(x_i,y_i,z_i).
\]

Define

\[
r_i = \sqrt{y_i^2+z_i^2},
\qquad
\phi_i = \operatorname{atan2}(z_i,y_i).
\]

The local radial unit vector is

\[
\mathbf e_{r,i}=(0,\cos\phi_i,\sin\phi_i).
\]

For a straight cylindrical reference surface, the radial displacement is

\[
h_i(t)=r_i(t)-R_0(t),
\]

and the dimensionless displacement is

\[
\eta_i(t)=\frac{h_i(t)}{R_0(t)}.
\]

Because the cylinder normal is radial, this equals the normal displacement to first order.

Use the area radius \(a=1\) only as the global simulation length normalization. The fitted tube radius \(R_0(t)\) remains a separate geometric quantity.

---

# 5. Frame preprocessing

## 5.1 Triangle and vertex weights

Compute triangle areas. Assign each vertex one third of the area of every incident triangle:

\[
A_i=\frac13\sum_{f\ni i}A_f.
\]

Use these vertex area weights in centering, fitting, and modal projection.

Reject or flag frames containing:

- NaN or infinite coordinates;
- zero-area triangles above a configurable fraction;
- invalid face indices;
- severe nonmanifold geometry;
- unreasonably large area or volume jumps.

## 5.2 Surface center

Use the area-weighted surface centroid as the default center:

\[
\mathbf c=
\frac{\sum_i A_i\mathbf P_i}{\sum_i A_i}.
\]

If the repository already has a reliable enclosed-volume centroid, make it an option.

## 5.3 Longitudinal axis

For the first frame:

1. subtract the weighted centroid;
2. form the area-weighted covariance matrix;
3. choose the eigenvector with the largest eigenvalue as the tube axis.

If the imposed elongational-flow direction is known, select the principal-axis sign and orientation closest to that direction.

For later frames:

- choose the sign that maximizes the dot product with the previous axis;
- optionally smooth tiny frame-to-frame orientation noise;
- flag large sudden changes.

Construct a right-handed orthonormal local frame and transform all points into it.

Do **not** independently recenter every axial cross-section if \(m=1\) bending modes are of interest. That would remove physical bending. Only remove global translation and global tilt.

---

# 6. Identify the cylindrical core

The rounded caps must not enter the modal projection.

Provide two core-selection methods.

## 6.1 Preferred fast method: persistent material core

Use when vertex IDs and connectivity remain consistent.

At the initial frame:

1. align the geometry;
2. identify vertices safely inside the initially cylindrical section;
3. exclude a configurable transition margin near both cap-cylinder junctions;
4. store the selected vertex IDs or a material-coordinate range.

Use the same material set in subsequent frames.

Advantages:

- the analysis domain does not change in response to growing modes;
- no cap vertices suddenly enter or leave the fit;
- mode amplitudes are less biased by mask motion.

This method should be preferred when available.

## 6.2 General method: data-driven core detection

Use when remeshing changes vertex IDs or topology.

For each frame:

1. bin vertices along the local \(x\)-axis;
2. calculate an area-weighted median or robust mean radius in each bin;
3. estimate the local radial slope and radius dispersion;
4. identify bins satisfying configurable cylindricality criteria:
   \[
   \left|\frac{d\bar r}{dx}\right| < s_{\max},
   \qquad
   \frac{\operatorname{std}(r)}{\bar r}<c_{\max};
   \]
5. keep the largest contiguous valid interval containing the tube center;
6. shrink that interval by a transition margin.

Avoid a core detector that aggressively follows every local pearling minimum and maximum. Smooth only enough to detect the large-scale cap transition.

Recommended defaults to expose in configuration:

```yaml
core:
  method: auto
  axial_bins: 120
  max_abs_radius_slope: 0.15
  max_relative_cross_section_std: 0.08
  transition_margin_in_radii: 1.0
  minimum_core_length_in_radii: 6.0
```

If a stable core cannot be found, flag the frame and do not silently report modal rates.

---

# 7. Time-dependent cylindrical base state

For the selected core, estimate:

\[
x_{\min}(t),\qquad x_{\max}(t),
\]

\[
L_c(t)=x_{\max}(t)-x_{\min}(t),
\]

\[
x_c(t)=\frac{x_{\max}(t)+x_{\min}(t)}{2},
\]

and the area-weighted mean core radius

\[
R_0(t)=
\frac{\sum_{i\in\mathrm{core}}A_i r_i}
     {\sum_{i\in\mathrm{core}}A_i}.
\]

Use a constant-radius base state by default. Do not remove a high-order polynomial or spline from the axial radius profile, because that can erase genuine long-wavelength unstable modes.

A low-order base correction may be implemented only as an explicitly selected option, and its effect on low modes must be documented.

Define the co-stretching axial coordinate

\[
\xi_i(t)=\frac{x_i(t)-x_c(t)}{L_c(t)},
\qquad
-\frac12\le \xi_i\le \frac12.
\]

This prevents uniform tube elongation from being misidentified as movement of spectral peaks between mode numbers.

For a mode with axial index \(n\) in a periodic Fourier convention,

\[
q_n(t)=\frac{2\pi n}{L_c(t)},
\qquad
q_nR_0(t)=\frac{2\pi nR_0(t)}{L_c(t)}.
\]

Record \(q_nR_0(t)\) at every frame.

---

# 8. End window

Apply a smooth axial window to reduce spectral leakage from the cap transitions.

Default to a Tukey window in \(\xi\):

```yaml
window:
  type: tukey
  alpha: 0.25
```

Combine the axial window with vertex area weights:

\[
W_i=A_i\,w(\xi_i).
\]

Store the effective windowed area

\[
A_{\mathrm{eff}}=\sum_iW_i.
\]

The pipeline should support a sensitivity check using at least two nearby window widths. Large variation of fitted growth rates under a modest window change should create a warning flag.

---

# 9. Modal decomposition

## 9.1 Periodic cylindrical harmonics

The default basis is

\[
\Phi_{mn}(\phi,\xi)
=
\exp\left[i\left(m\phi+2\pi n\xi\right)\right].
\]

Use:

```text
m = -m_max, ..., m_max
n = -n_max, ..., n_max
```

with the understanding that a real displacement field has conjugate symmetry.

The raw quadrature coefficient is

\[
a_{mn}(t)=
\frac{
\sum_i W_i\eta_i(t)
\exp[-i(m\phi_i+2\pi n\xi_i)]
}{
\sum_i W_i
}.
\]

This is a nonuniform discrete Fourier projection using the actual surface quadrature weights.

Implement it using vectorized matrix multiplication rather than Python loops:

```text
E_phi[i,m] = exp(-1j * m * phi[i])
E_x[i,n]   = exp(-1j * 2*pi*n * xi[i])

A_modes = E_phi^H @ ((W * eta)[:,None] * E_x) / sum(W)
```

Chunk the mode calculation if required for memory.

## 9.2 Weighted least-squares option

Provide an optional weighted least-squares projection:

\[
\min_{\{a_{mn}\}}
\sum_i W_i
\left|
\eta_i-\sum_{mn}a_{mn}\Phi_{mn}(\phi_i,\xi_i)
\right|^2.
\]

This is slower but reduces leakage when sampling is irregular or the window makes the basis noticeably nonorthogonal.

Recommended configuration:

```yaml
projection:
  method: quadrature_nudft
  alternative_method: weighted_lstsq
  regularization: 1.0e-12
```

Use QR, SVD, or a stable Hermitian solve. Avoid explicitly inverting matrices.

## 9.3 Axisymmetric fast path

When `m_max: 0`, avoid constructing the full azimuthal basis.

Either:

1. directly evaluate
   \[
   a_{0n}=
   \frac{\sum_iW_i\eta_ie^{-i2\pi n\xi_i}}
        {\sum_iW_i},
   \]
   or

2. robustly bin and average \(\eta\) over cross-sections, then perform a one-dimensional transform.

The direct weighted formula is preferable because it avoids an interpolation step.

## 9.4 Alternative finite-interval basis

Support, but do not default to, a cosine basis:

\[
\cos\left(\frac{n\pi(x-x_{\min})}{L_c}\right).
\]

Use this only when the seeded perturbations and theoretical comparison explicitly use finite-interval cosine modes.

The configuration must make the wavenumber convention unambiguous:

```yaml
modes:
  axial_basis: fourier_periodic
  m_max: 0
  n_max: 24
```

For `fourier_periodic`, use \(q_n=2\pi n/L_c\).

For `cosine_neumann`, use \(q_n=n\pi/L_c\).

---

# 10. Grouped amplitudes and mode symmetry

Always store the raw complex coefficients.

For growth-rate fitting, use a phase-invariant grouped amplitude. Group symmetry-related coefficients with the same absolute mode indices:

\[
\mathcal G_{mn}
=
\{(s_m m,s_n n):s_m,s_n\in\{-1,+1\}\},
\]

removing duplicate entries when \(m=0\) or \(n=0\).

Define

\[
A_{mn}(t)
=
\left[
\sum_{(p,q)\in\mathcal G_{mn}}
|a_{pq}(t)|^2
\right]^{1/2}.
\]

Any constant normalization factor changes the reported amplitude but not

\[
\frac{d}{dt}\log A_{mn}.
\]

Document the exact convention in `resolved_config.yaml`.

Exclude the base mode `(m,n)=(0,0)` from instability growth-rate reporting. Retain it for diagnostics.

Typical interpretation:

```text
m = 0: axisymmetric varicose or pearling modes
m = 1: bending or sinuous modes
m = 2: elliptical cross-section modes
m > 2: higher azimuthal wrinkles or fluting
```

---

# 11. Reference, plus/minus, and ordinary series

The pipeline should support four analysis modes.

## 11.1 Ordinary single series

Analyze the coefficients directly:

\[
a_{mn}^{\mathrm{signal}}(t)=a_{mn}(t).
\]

This is the least clean option because base-state relaxation and numerical transients can contaminate low modes.

## 11.2 Perturbed minus unperturbed reference

For a perturbation and matched unperturbed simulation:

\[
\delta a_{mn}(t)
=
a_{mn}^{\mathrm{perturbed}}(t)
-
a_{mn}^{\mathrm{reference}}(t).
\]

Match frames by physical time, interpolating complex coefficients in time only when necessary.

## 11.3 Paired `+epsilon` and `-epsilon`

For matched simulations initialized with opposite perturbations:

\[
\delta a_{mn}(t)
=
\frac{
a_{mn}^{+}(t)-a_{mn}^{-}(t)
}{2\epsilon}.
\]

This cancels the common evolving base state and leading even-order nonlinear terms.

It does not require pointwise mesh correspondence if modal coefficients are computed separately before subtraction.

## 11.4 Plus/minus with reference

If all three are available, use the plus/minus signal for the linear response and the reference run for additional diagnostics of base-state evolution and numerical noise.

---

# 12. Efficient simulation design

For the most efficient parameter sweep, support a low-amplitude multisine initial perturbation rather than requiring one simulation per mode.

For example,

\[
\eta(\phi,\xi,0)
=
\epsilon
\sum_{(m,n)\in\mathcal S}
c_{mn}
\cos(m\phi+2\pi n\xi+\psi_{mn}),
\]

with random fixed phases \(\psi_{mn}\) and coefficients normalized so that the total RMS perturbation is small.

Recommended experiment per physical parameter set:

```text
one +multisine simulation
one -multisine simulation
optional unperturbed reference simulation
```

This can recover the short-time rates of many modes from two or three simulations, provided:

- every seeded amplitude is above the numerical noise floor;
- the total perturbation remains in the linear regime;
- output is frequent enough to resolve growth or decay;
- the selected modes are spatially resolved;
- modal cross-coupling remains small over the fitting interval.

As a starting point, keep total radial perturbation RMS below roughly one percent of the mean tube radius and verify amplitude independence with at least one smaller-amplitude test.

---

# 13. Spatial-resolution checks

For each axial mode, the wavelength is

\[
\lambda_n=\frac{2\pi}{q_n}.
\]

Estimate the effective axial point spacing in the core. Flag modes with fewer than a configurable number of samples per wavelength.

Recommended default:

```yaml
quality:
  minimum_axial_samples_per_wavelength: 10
  minimum_azimuthal_samples_per_wavelength: 8
```

For azimuthal index \(m\), require sufficient sampling around the circumference.

Also flag modes too close to the mesh Nyquist limit or modes for which the projection changes strongly under mesh refinement.

---

# 14. Growth-rate extraction

## 14.1 Primary fit

For each grouped amplitude \(A_{mn}(t)\), fit

\[
\log A_{mn}(t)=b+\sigma_{mn}t
\]

over an early-time linear window.

Use ordinary least squares by default and robust regression as an option.

Return:

- slope \(\sigma_{mn}\);
- standard error;
- confidence interval when available;
- \(R^2\);
- fit residuals;
- number of points;
- amplitude ratio over the window.

The output times need not be uniform.

## 14.2 Optional phase rate

For a single complex coefficient with coherent phase, unwrap its phase and fit

\[
\arg a_{mn}(t)=\phi_0+\omega_{mn}t.
\]

Report \(\omega_{mn}\) separately from the real growth rate.

## 14.3 Instantaneous diagnostic rate

Optionally compute

\[
\sigma_{\mathrm{inst}}(t)
=
\frac{d}{dt}\log A_{mn}(t)
\]

using a Savitzky-Golay derivative or a local polynomial fit.

This is diagnostic only. The primary reported value should be the windowed exponential-fit slope.

## 14.4 Automatic fitting-window selection

Allow the user to specify a fixed fitting interval. Also implement an automatic mode.

A candidate frame is eligible when:

1. its geometry passed quality checks;
2. the mode amplitude is above the estimated noise floor;
3. the dimensionless amplitude is below a nonlinear threshold;
4. the mode is spatially resolved;
5. the core detector is stable;
6. no topological event or severe shape transition has occurred.

Recommended configurable thresholds:

```yaml
growth_fit:
  mode: auto
  minimum_points: 6
  minimum_amplitude_ratio: 1.25
  maximum_dimensionless_amplitude: 0.05
  minimum_r_squared: 0.90
  maximum_relative_qR_change: 0.15
  robust_regression: false
```

Automatic-window algorithm:

1. enumerate contiguous early-time windows satisfying the minimum number of points;
2. reject windows violating amplitude, geometry, or wavenumber criteria;
3. fit `log(amplitude)` versus time;
4. score each window using:
   - early start time;
   - larger number of points;
   - high \(R^2\);
   - stable slope under removal of one endpoint;
   - sufficient amplitude change;
5. choose the earliest stable window, not simply the window with the highest \(R^2\).

Do not report a confident rate when the amplitude remains at the numerical noise floor or changes too little.

## 14.5 Changing base state

Because the cylinder stretches,

\[
q_nR_0
\]

may change during the fit.

Record the minimum, maximum, and mean \(qR\) in every fit window.

If \(qR\) changes substantially, either:

- shorten the fit window;
- report a local or piecewise rate;
- or mark the rate as an average over a changing base state.

Do not present it as a rate at one sharply defined wavenumber without a warning.

---

# 15. Noise-floor estimation

Provide at least one of the following methods:

1. use coefficients from an unperturbed reference run;
2. use high-wavenumber modes outside the seeded range;
3. use the variation among repeated simulations;
4. use the earliest frames before measurable growth;
5. use projection residuals.

Store the estimated noise amplitude per mode.

A mode should not be fit until its amplitude exceeds a configurable multiple of the noise floor, for example:

```yaml
growth_fit:
  minimum_signal_to_noise: 5.0
```

For strongly decaying modes, the fit may end when the amplitude approaches the noise floor.

---

# 16. Diagnostics and quality flags

Create explicit quality flags rather than silently dropping questionable results.

Suggested flags:

```text
INVALID_GEOMETRY
AREA_DRIFT
VOLUME_DRIFT
AXIS_JUMP
CORE_NOT_FOUND
CORE_TOO_SHORT
CAP_CONTAMINATION
LARGE_CYLINDRICALITY_ERROR
MODE_UNRESOLVED
AMPLITUDE_BELOW_NOISE
AMPLITUDE_NONLINEAR
INSUFFICIENT_TIME_POINTS
POOR_EXPONENTIAL_FIT
QR_CHANGES_TOO_MUCH
WINDOW_SENSITIVE
LARGE_RECONSTRUCTION_ERROR
TOPOLOGY_CHANGED
```

## 16.1 Cylindricality error

Compute

\[
E_{\mathrm{cyl}}
=
\left[
\frac{\sum_i W_i(r_i-R_0)^2}
     {\sum_i W_i}
\right]^{1/2}
\frac{1}{R_0}.
\]

This measures total radial variation, including physical modes, so interpret it as a linear-regime diagnostic rather than an error that must always be tiny.

## 16.2 Reconstruction error

Reconstruct the displacement using retained modes:

\[
\eta_i^{\mathrm{rec}}
=
\sum_{mn}a_{mn}\Phi_{mn}(\phi_i,\xi_i).
\]

Compute the weighted relative residual:

\[
E_{\mathrm{rec}}
=
\frac{
\left[\sum_iW_i|\eta_i-\eta_i^{\mathrm{rec}}|^2\right]^{1/2}
}{
\left[\sum_iW_i|\eta_i|^2\right]^{1/2}
}.
\]

Use the same projection convention used to obtain the coefficients.

## 16.3 Sensitivity checks

Optionally repeat the analysis with:

- a slightly smaller core;
- a different Tukey-window width;
- one higher and one lower `n_max`;
- an alternative axis fit;
- the weighted least-squares projection.

Store differences in a sensitivity table.

---

# 17. Suggested plots

Generate, at minimum:

1. area, volume, reduced volume, \(R_0(t)\), and \(L_c(t)\);
2. \(q_nR_0(t)\) for analyzed modes;
3. semilog modal amplitude versus time;
4. selected linear fits overlaid on amplitudes;
5. growth rate versus mean \(qR\);
6. initial and final spectra;
7. reconstruction error versus time;
8. core boundaries overlaid on representative meridian profiles;
9. window-sensitivity comparison when enabled.

Use one plot per figure unless the repository has a standard dashboard style.

---

# 18. Configuration example

```yaml
discovery:
  series_glob: "*"
  frame_glob: "*"
  recursive: true

io:
  reader: auto
  time_source: metadata_or_filename

geometry:
  center: area_centroid
  axis_method: weighted_pca
  use_known_flow_axis: true
  enforce_axis_sign_continuity: true

core:
  method: auto
  axial_bins: 120
  max_abs_radius_slope: 0.15
  max_relative_cross_section_std: 0.08
  transition_margin_in_radii: 1.0
  minimum_core_length_in_radii: 6.0

window:
  type: tukey
  alpha: 0.25

modes:
  axial_basis: fourier_periodic
  m_max: 0
  n_min: 1
  n_max: 24
  normalize_displacement_by_radius: true
  coordinate: co_stretching

projection:
  method: quadrature_nudft
  regularization: 1.0e-12
  chunk_size: 32

pairing:
  enabled: true
  prefer_plus_minus: true
  allow_reference_subtraction: true
  time_match_tolerance: 1.0e-8

growth_fit:
  mode: auto
  fixed_start_time: null
  fixed_end_time: null
  minimum_points: 6
  minimum_amplitude_ratio: 1.25
  minimum_signal_to_noise: 5.0
  maximum_dimensionless_amplitude: 0.05
  minimum_r_squared: 0.90
  maximum_relative_qR_change: 0.15
  robust_regression: false
  compute_instantaneous_rate: true

quality:
  minimum_axial_samples_per_wavelength: 10
  minimum_azimuthal_samples_per_wavelength: 8
  maximum_area_relative_drift: 0.01
  maximum_volume_relative_drift: 0.01
  maximum_axis_angle_jump_degrees: 5.0
  maximum_reconstruction_error: 0.25

output:
  save_csv: true
  save_parquet: true
  save_mat: false
  save_plots: true
  overwrite: false

parallel:
  enabled: true
  workers: auto
  parallelize_over_series: true
```

The final implementation should validate the configuration and save the fully resolved version for every run.

---

# 19. Recommended software architecture

Follow the repository language and conventions. A possible modular structure is:

```text
mode_analysis/
    cli
    config
    discovery
    io
    mesh_geometry
    alignment
    core_detection
    cylindrical_projection
    series_pairing
    growth_fitting
    diagnostics
    plotting
    output
```

Key logical interfaces:

```text
discover_series(root, config) -> list[SeriesDescriptor]

load_series(descriptor, config) -> list[GeometryFrame]

preprocess_frame(frame, previous_state, config) -> ProcessedFrame

detect_core(processed_frame, initial_state, config) -> CoreDefinition

project_modes(processed_frame, core, config) -> ModeFrameResult

combine_related_series(mode_results, pairing_metadata, config)
    -> LinearResponseSeries

fit_growth_rates(linear_response_series, frame_metrics, config)
    -> GrowthRateTable

write_outputs(...)
```

Parallelize over simulation series, not over individual frames within one series, unless memory constraints require otherwise. Frame order matters for axis continuity, core tracking, and diagnostics.

---

# 20. Pseudocode

```text
config = load_and_validate_config()

series_list = discover_series(root, config)

for series in parallel(series_list):

    frames = load_series(series)
    frames = sort_by_physical_time(frames)

    previous_state = None
    initial_state = None
    frame_metrics = []
    coefficients = []

    for frame in frames:

        processed = preprocess_frame(
            frame,
            previous_state=previous_state,
            config=config
        )

        if initial_state is None:
            initial_state = initialize_core_tracking(processed, config)

        core = detect_or_track_core(
            processed,
            initial_state=initial_state,
            config=config
        )

        base = fit_cylindrical_base(processed, core, config)

        mode_result = project_modes(
            processed,
            core,
            base,
            config
        )

        diagnostics = compute_frame_diagnostics(
            processed,
            core,
            base,
            mode_result,
            config
        )

        frame_metrics.append(diagnostics.metrics)
        coefficients.append(mode_result)

        previous_state = processed.state

    save_raw_frame_and_mode_results(series, frame_metrics, coefficients)

group matched baseline, plus, and minus series

for response_group:

    response = construct_linear_response(response_group, config)

    noise = estimate_noise(response_group, response, config)

    rates = fit_growth_rates(
        response,
        frame_metrics,
        noise,
        config
    )

    create_plots(...)
    save_growth_rate_outputs(...)
```

---

# 21. Unit and regression tests

Implement synthetic tests before applying the pipeline to production simulations.

## 21.1 Exact static cylinder

Generate vertices on a perfect cylinder with a known core and verify:

```text
all nonzero modes are near machine or discretization noise
R0 and Lc are recovered correctly
```

## 21.2 Single imposed axisymmetric mode

Use

\[
\eta=A\cos(2\pi n\xi+\psi).
\]

Verify:

```text
the correct n is recovered
amplitude is independent of phase
neighboring modes are small
```

## 21.3 Single azimuthal mode

Use

\[
\eta=A\cos(m\phi+\psi).
\]

Verify recovery of \(m\).

## 21.4 Combined mode

Use

\[
\eta=A\cos(m\phi+2\pi n\xi+\psi).
\]

Verify grouped amplitude and conjugate symmetry.

## 21.5 Known exponential growth

Generate synthetic frames with

\[
A(t)=A_0e^{\sigma t}
\]

and verify that the recovered growth rate matches the prescribed \(\sigma\).

Test:

- uniform and nonuniform time spacing;
- small additive noise;
- decaying modes;
- changing \(R_0(t)\);
- affine axial stretching;
- irregular vertex sampling;
- altered mesh connectivity.

## 21.6 Cap contamination

Create a cylinder with hemispherical caps. Verify that:

- full-surface projection produces leakage;
- core masking and windowing suppress it;
- recovered central mode amplitudes are accurate.

## 21.7 Plus/minus cancellation

Create synthetic base relaxation common to two series and opposite linear perturbations. Verify that plus/minus subtraction removes the common base contribution.

---

# 22. Acceptance criteria

The implementation is complete when it can:

1. analyze all series under a root folder without manual per-series editing;
2. read the repository's actual geometry format;
3. consistently align each tube;
4. identify or track a cylindrical core;
5. compute modal coefficients for every valid frame;
6. correctly account for changing tube length and radius;
7. combine plus/minus or reference runs when available;
8. fit and report early-time modal growth rates;
9. flag unresolved, noisy, nonlinear, or poorly fit modes;
10. produce tests that recover prescribed synthetic growth rates;
11. save reproducible configuration, tables, logs, and diagnostic plots.

---

# 23. Scientific cautions to preserve in code and documentation

- A finite capped tube is not an infinite translationally invariant cylinder. The extracted Fourier modes are a controlled local approximation on the central cylindrical region.
- Uniform stretching changes the physical wavenumber even when the normalized mode index \(n\) is fixed.
- Subtracting an overly flexible smooth base profile can remove genuine low-wavenumber modes.
- Cross-sectional recentering can erase physical \(m=1\) bending.
- The earliest frames may contain cap-junction or discretization transients.
- A high \(R^2\) alone does not prove a valid linear regime.
- Modes below the noise floor or near the spatial-resolution limit must not be assigned confident rates.
- If substantial off-diagonal mode transfer appears immediately, individual Fourier growth rates may not correspond to independent finite-vesicle eigenmodes.
- Report the fitting interval and the range of \(qR\) with every rate.

---

# 24. Recommended first production configuration

For the first implementation and validation:

```text
analyze only m = 0
use n = 1 through approximately 20
use the persistent initial core if vertex IDs are preserved
use a Tukey window with alpha = 0.25
use direct area-weighted nonuniform Fourier projection
fit log amplitude over a manually inspected early-time interval
compare against synthetic tests before enabling automatic windows
then add plus/minus pairing and general m > 0 support
```

This staged approach minimizes implementation risk while immediately producing the axisymmetric pearling growth rates needed for comparison with a cylindrical linear-stability analysis.
