# Specification: Extract Early-Time Membrane Energy Relaxation Timescales

## Objective

Process a collection of membrane-energy time series indexed by the nondimensional Saffman–Delbrück number, `Sd`.

Each simulation contains exactly **300 timesteps**. The physically relevant relaxation occurs near the beginning of the run, typically within approximately the first **20–30 timesteps**. Later-time oscillations are believed to be spurious and must not influence the extracted relaxation timescale.

The analysis must:

1. Estimate the early-time relaxation timescale for every simulation.
2. Reject or ignore the later oscillatory portion of each trace.
3. quantify uncertainty associated with the fitting-window choice.
4. Test the expected scaling with `Sd`.
5. Produce diagnostic plots and machine-readable summary tables.
6. Avoid requiring dimensional values of surface viscosity, bulk viscosity, bending modulus, or radius.

The codebase and file layout are already known. Integrate this analysis into the existing project without redesigning the input layer.

---

## Physical Background

Let the stored nondimensional time be

\[
\hat t = \frac{t}{T_s},
\]

where `T_s` is the surface-viscous bending timescale.

For a characteristic system length \(R\),

\[
T_s \sim \frac{\mu_s R^2}{\kappa},
\qquad
T_b \sim \frac{\eta_b R^3}{\kappa},
\]

and

\[
Sd = \frac{\mu_s}{\eta_b R},
\]

up to any fixed convention-dependent prefactor already present in the governing nondimensional equations.

Therefore,

\[
\frac{T_s}{T_b} \propto Sd.
\]

If the code's nondimensionalization contains no additional factor, the bulk-scaled nondimensional time is

\[
t_b^* = Sd\,\hat t.
\]

If inspection of the nondimensional equations shows that the bulk dissipation term appears as

\[
\frac{C}{Sd}\mathcal D_b,
\]

then instead use

\[
t_b^* = \frac{Sd}{C}\hat t.
\]

Keep this factor configurable. Call it `bulk_time_factor`, with default value `1.0`, so that

\[
t_b^* = bulk\_time\_factor \times Sd \times \hat t.
\]

The expected relaxation-time scaling in surface-time units is

\[
\hat\tau_E \sim
\begin{cases}
\text{constant}, & Sd \gg 1,\\
Sd^{-1}, & Sd \ll 1.
\end{cases}
\]

Equivalently,

\[
Sd\,\hat\tau_E
\]

should approach a constant at low `Sd`, again up to the fixed convention factor.

---

## Important Interpretation of the Energy Timescale

Fit the early energy relaxation using

\[
E(t) = E_\infty + A\exp\left[-\frac{t-t_0}{\tau_E}\right].
\]

The fitted \(\tau_E\) is the **energy-decay timescale**.

Near equilibrium, if a shape amplitude \(x\) relaxes according to

\[
x(t) \propto e^{-t/\tau_x}
\]

and the excess energy is quadratic,

\[
E-E_\infty \propto x^2,
\]

then

\[
\tau_x = 2\tau_E.
\]

Store both:

- `tau_energy = tau_E`
- `tau_shape_equivalent = 2 * tau_E`

Do not silently replace one with the other.

---

## Required Inputs

For every simulation, obtain:

- `Sd`
- nondimensional time array `t`
- energy array `E`

Assume:

- each run has 300 samples;
- time may or may not be uniformly spaced;
- the initial relaxation generally finishes within the first 20–30 samples;
- later samples may contain oscillations of varying amplitude.

Validate that:

- `len(t) == len(E)`;
- all values are finite after any project-specific preprocessing;
- time is strictly increasing;
- at least 12 usable early-time points remain.

If a trace fails validation, record a failed status rather than crashing the entire batch.

---

## General Analysis Strategy

Use a two-stage strategy:

1. **Detect a conservative upper bound for the clean early-time region.**
2. **Fit the exponential repeatedly over several candidate endpoints and retain a stable window plateau.**

Do not choose one fitting endpoint based only on the nominal statement that relaxation occurs by timestep 20–30. Instead, use that knowledge to constrain the search.

The final fit must use the original unsmoothed energy values. Smoothing is permitted only for diagnostics and automatic interval detection.

---

## Default Configuration

Expose the following settings in a configuration object or clearly grouped constants:

```python
N_TOTAL = 300

START_INDEX_MIN = 0
START_INDEX_MAX = 5

MIN_FIT_POINTS = 8
MAX_EARLY_INDEX = 40

CANDIDATE_END_MIN = 10
CANDIDATE_END_MAX = 35
CANDIDATE_END_STEP = 1

SMOOTH_WINDOW = 5
SMOOTH_POLYORDER = 2

PERSISTENT_RISE_POINTS = 3
DERIVATIVE_NOISE_MULTIPLIER = 3.0

ROBUST_LOSS = "soft_l1"
F_SCALE = 1.0

MIN_R2 = 0.98
MAX_RELATIVE_TAU_SE = 0.5

PLATEAU_RELATIVE_TOLERANCE = 0.15
MIN_PLATEAU_WINDOWS = 4

BULK_TIME_FACTOR = 1.0
```

These defaults may be adjusted slightly to match the existing numerical stack, but preserve the intent.

The primary candidate fit endpoints should be between approximately samples 10 and 35. Never allow the automatic initial-relaxation fit to extend beyond sample 40 unless explicitly overridden.

---

## Step 1: Preprocess Only for Detection

Create a diagnostic smoothed energy trace, for example with a Savitzky–Golay filter.

Requirements:

- Use an odd smoothing window.
- Default to 5 points.
- Do not smooth across missing values.
- Do not use the smoothed trace in the final nonlinear fit.
- If there are too few points for Savitzky–Golay, use a short centered moving average.

Denote the smoothed trace by `E_smooth`.

Compute a numerical derivative using actual time spacing:

```python
dE_dt = np.gradient(E_smooth, t)
```

Estimate a derivative noise scale from a short part of the early trace after the steepest initial drop, or with a robust statistic such as

```python
sigma_d = 1.4826 * median(abs(dE_dt - median(dE_dt)))
```

Use this only as a heuristic detector.

---

## Step 2: Choose the Fit Start

The default fit start is the first valid sample, but avoid obvious startup artifacts.

Test candidate start indices from 0 through 5.

For each candidate start:

- require that the energy initially decreases overall;
- reject a start point if it is an isolated discontinuity or extreme one-step jump inconsistent with the next several samples;
- prefer the earliest start that produces stable fits.

Do not move the start later merely to improve the fit unless the first samples are demonstrably numerical startup artifacts.

Store the selected start index as `fit_start_index`.

---

## Step 3: Detect the Onset of Oscillations or Reversal

Define a conservative oscillation-onset index, `oscillation_start_index`.

Search only after the first several samples, for example after index 5.

A candidate onset occurs when either of the following is true:

### Criterion A: Persistent Increase

The smoothed derivative is significantly positive for at least `PERSISTENT_RISE_POINTS` consecutive samples:

\[
\frac{dE}{dt} >
DERIVATIVE\_NOISE\_MULTIPLIER \times \sigma_d.
\]

### Criterion B: Repeated Derivative Sign Changes

Within a short rolling window, the derivative changes sign more than once after the main decay has substantially slowed.

This criterion should not fire during the steep initial transient due to small derivative noise.

### Conservative Default

If no clear onset is detected, set

```python
oscillation_start_index = min(40, len(E) - 1)
```

The fit-end search must satisfy

```python
fit_end_index < oscillation_start_index
```

unless no valid fit is otherwise possible.

Because the physical relaxation is expected within the first 20–30 samples, do not use late data to improve the equilibrium estimate.

---

## Step 4: Fit the Early-Time Exponential

For every candidate fitting interval

```python
[start_index : end_index + 1]
```

fit

\[
E(t)=E_\infty+A\exp[-(t-t_0)/\tau_E].
\]

Use nonlinear least squares with a robust loss.

### Parameter Bounds

For a decaying trace, enforce:

```python
A > 0
tau_E > 0
```

Constrain `E_inf` to a physically reasonable range based on the early data. A useful default is:

```python
E_min_early - margin <= E_inf <= E_end_candidate + margin
```

where the margin is a modest fraction of the early energy range.

Do not determine `E_inf` from the final samples of the 300-step simulation because those samples contain spurious oscillations.

### Initial Guesses

Use:

```python
E_inf_guess = median(E over the last 3 points of the candidate interval)
A_guess = max(E[start] - E_inf_guess, small_positive_value)
tau_guess = 0.25 * (t[end] - t[start])
```

A log-linear fit may be used only to generate an initial guess after selecting a provisional `E_inf`. It must not be the final estimator.

### Weighting

Use equal weights by default.

If the numerical energy noise is demonstrably heteroscedastic, permit user-supplied weights, but do not infer aggressive weighting from one trace.

---

## Step 5: Candidate Fit Endpoints

Construct candidate endpoints using:

```python
lower_end = max(
    fit_start_index + MIN_FIT_POINTS - 1,
    CANDIDATE_END_MIN
)

upper_end = min(
    CANDIDATE_END_MAX,
    oscillation_start_index - 1,
    len(E) - 1
)
```

Test every integer endpoint in this interval.

If fewer than four candidate endpoints are available, relax `CANDIDATE_END_MAX` only up to `MAX_EARLY_INDEX = 40`.

For each candidate fit, store:

- start index;
- end index;
- fitted `E_inf`;
- fitted `A`;
- fitted `tau_energy`;
- `tau_shape_equivalent`;
- parameter standard errors if available;
- residual sum of squares;
- RMSE;
- \(R^2\);
- Akaike information criterion if convenient;
- convergence status;
- boundary-hit flags.

---

## Step 6: Reject Bad Fits

Reject a candidate fit if any of the following applies:

- optimizer did not converge;
- `tau_energy <= 0`;
- `A <= 0`;
- fitted parameters are on imposed bounds;
- fewer than `MIN_FIT_POINTS` are present;
- fitted `E_inf` is clearly incompatible with a decaying approach;
- normalized RMSE is excessive;
- \(R^2 < MIN_R2\), unless the trace is flagged as non-single-exponential;
- relative standard error of `tau_energy` exceeds `MAX_RELATIVE_TAU_SE`;
- residuals show a strong systematic oscillation within the fit interval;
- the fitted timescale is much smaller than one timestep or much larger than the full candidate interval without a correspondingly good fit.

Do not reject solely because `R^2` is slightly below the threshold if the energy range is extremely small. Record a lower-confidence status instead.

---

## Step 7: Select a Stable Fitting-Window Plateau

The preferred timescale is not the result from one arbitrary endpoint.

Sort valid fits by endpoint and inspect the sequence

\[
\tau_E(k),
\]

where \(k\) is the fit-end index.

Find the longest contiguous group of at least `MIN_PLATEAU_WINDOWS` candidate endpoints for which the timescales remain stable.

A practical criterion is:

```python
abs(tau_i - median_tau) / median_tau <= PLATEAU_RELATIVE_TOLERANCE
```

with default tolerance 15%.

Prefer a plateau that:

1. lies before the detected oscillation onset;
2. includes endpoints near the expected relaxation completion, approximately samples 20–30;
3. has good residual metrics;
4. is not formed only from the shortest available windows.

Select the final timescale as:

```python
tau_energy = median(tau values in plateau)
```

Select final `E_inf` and `A` either from:

- the fit whose `tau` is closest to the plateau median, or
- a refit using the median plateau endpoint.

Use the first approach by default to avoid inventing a new endpoint after selection.

### Window-Selection Uncertainty

Calculate:

```python
tau_window_low = percentile(plateau_tau, 16)
tau_window_high = percentile(plateau_tau, 84)
tau_window_std = robust_std(plateau_tau)
```

Combine this with the nonlinear-fit parameter uncertainty if available.

A simple conservative total uncertainty is

\[
\sigma_{\tau,\mathrm{total}}
=
\sqrt{
\sigma_{\tau,\mathrm{fit}}^2
+
\sigma_{\tau,\mathrm{window}}^2
}.
\]

Record both components separately.

---

## Step 8: Local Decay-Rate Diagnostic

For the selected fit, define

\[
\Delta E(t)=E(t)-E_\infty.
\]

Only evaluate points where \(\Delta E>0\) by a safe numerical margin.

Compute

\[
r_E(t)
=
-\frac{d}{dt}\log[\Delta E(t)]
\]

using a smoothed diagnostic signal.

For a single exponential,

\[
r_E(t)\approx \frac{1}{\tau_E}.
\]

Generate a plot of `r_E(t)` over the first 40 samples and draw a horizontal line at `1 / tau_energy`.

This is diagnostic only. Do not replace the nonlinear fit with noisy pointwise derivative estimates.

Flag a trace as `multimode_or_nonlinear` if the local rate drifts systematically rather than showing an approximate plateau.

---

## Step 9: Threshold-Based Cross-Check

Using the fitted `E_inf`, define

\[
y(t)
=
\frac{E(t)-E_\infty}
     {E(t_0)-E_\infty}.
\]

Interpolate the first time at which

\[
y(t)=e^{-1}.
\]

Define

```python
tau_1e = t_at_y_equals_1_over_e - t0
```

For an ideal single exponential,

```python
tau_1e ≈ tau_energy
```

Record the relative discrepancy:

```python
tau_crosscheck_error = abs(tau_1e - tau_energy) / tau_energy
```

Use this as a quality diagnostic.

Do not use thresholds very close to equilibrium because the later oscillations contaminate them.

---

## Step 10: Optional Competing Models

The primary reported model must remain the single exponential.

However, if the semilog curve is visibly curved or the local decay rate drifts, optionally fit a double exponential over the same early interval:

\[
E(t)
=
E_\infty
+
A_1e^{-(t-t_0)/\tau_1}
+
A_2e^{-(t-t_0)/\tau_2}.
\]

Only retain this as a diagnostic if:

- both amplitudes are positive;
- timescales are well separated enough to be identifiable;
- the fit materially improves AIC or BIC;
- parameters are not at bounds.

Do not automatically report the slower double-exponential timescale as the physical answer. The main batch comparison should use the consistently defined early-time single-exponential effective timescale.

---

## Step 11: Per-Trace Diagnostic Plots

For every simulation, create one diagnostic figure containing clearly separated panels or separate image files showing:

### A. Full Energy Trace

Plot all 300 samples.

Overlay:

- selected fit start;
- selected fit end;
- detected oscillation onset;
- fitted early-time exponential extended only modestly beyond the fit interval;
- fitted equilibrium energy as a horizontal line.

The full trace is useful for verifying that late oscillations are excluded.

### B. Early-Time Energy Trace

Plot only samples 0–40.

Show the raw data and fitted curve.

### C. Semilog Excess Energy

Plot

\[
\log(E-E_\infty)
\]

over the selected fitting interval.

A single exponential should look approximately linear.

### D. Timescale Versus Candidate Endpoint

Plot fitted `tau_energy` against candidate end index.

Highlight:

- accepted fits;
- rejected fits;
- selected plateau;
- final median timescale.

This is the most important fit-window diagnostic.

### E. Local Decay Rate

Plot `r_E(t)` and `1/tau_energy`.

Save figures with filenames that include the `Sd` value or the existing simulation identifier.

---

## Step 12: Batch Scaling Analysis

After extracting one timescale per trace, create a summary dataframe sorted by `Sd`.

Required derived columns:

```python
tau_surface = tau_energy
tau_shape_surface = 2 * tau_energy

tau_bulk_scaled = BULK_TIME_FACTOR * Sd * tau_energy
tau_shape_bulk_scaled = BULK_TIME_FACTOR * Sd * 2 * tau_energy

decay_rate_surface = 1 / tau_energy
decay_rate_bulk_scaled = 1 / tau_bulk_scaled
```

The name `tau_surface` means the fitted nondimensional time measured in the simulation's existing surface-time units.

---

## Step 13: Fit the Sd Dependence

Fit the crossover model

\[
\hat\tau_E(Sd)=a+\frac{b}{Sd}.
\]

Use positive bounds:

```python
a > 0
b > 0
```

Weight by the extracted timescale uncertainty when reliable uncertainties are available. Otherwise use equal weights and report that fact.

Store:

- `a`;
- `b`;
- confidence intervals;
- weighted RMSE;
- \(R^2\) in linear space;
- residuals versus `Sd`.

Treat this as a motivated interpolation, not a guaranteed exact law.

---

## Step 14: Estimate Asymptotic Log-Log Slopes

On log-log axes, estimate

\[
p = \frac{d\log\tau}{d\log Sd}.
\]

Do this separately in low- and high-`Sd` subsets.

Do not hard-code the low/high split unless the project already has known ranges. By default:

- sort by `Sd`;
- use approximately the lowest third for the low-`Sd` slope;
- use approximately the highest third for the high-`Sd` slope;
- require at least three points in each fit.

Expected values:

```text
low Sd:   slope ≈ -1
high Sd:  slope ≈ 0
```

Report uncertainty on each slope.

If there are too few simulations, omit numerical slope fitting and show only the plots.

---

## Step 15: Required Batch Plots

Create the following plots.

### 1. Normalized Energy in Surface Time

For each run, plot

\[
y(t)
=
\frac{E(t)-E_\infty}
     {E(t_0)-E_\infty}
\]

against the original nondimensional time \(\hat t\).

Restrict the x-axis to the early-time region, preferably samples 0–40 or the corresponding time range.

High-`Sd` curves should collapse if the relaxation follows the surface-viscous timescale.

### 2. Normalized Energy in Bulk Time

Plot the same normalized energy against

\[
t_b^* =
BULK\_TIME\_FACTOR \times Sd \times \hat t.
\]

Low-`Sd` curves should collapse if the relaxation follows the bulk-viscous timescale.

### 3. Timescale Versus Sd

Plot

\[
\tau_E
\]

versus `Sd` on log-log axes.

Overlay:

- the fit \(a+b/Sd\);
- a reference slope of `-1` in the low-`Sd` region;
- a horizontal reference in the high-`Sd` region if helpful.

### 4. Bulk-Scaled Timescale Versus Sd

Plot

\[
BULK\_TIME\_FACTOR \times Sd \times \tau_E
\]

versus `Sd`.

Look for a low-`Sd` plateau.

### 5. Decay Rate Versus Sd

Plot

\[
1/\tau_E
\]

versus `Sd` on log-log axes.

Expected behavior:

```text
low Sd:   rate proportional to Sd
high Sd:  rate approaches a constant
```

Do not claim a scaling regime unless there are enough points spanning a meaningful range.

---

## Step 16: Output Table

Write a CSV or parquet file with one row per simulation.

Include at least:

```text
simulation_id
Sd
status
quality_flag

n_total
fit_start_index
fit_end_index
oscillation_start_index
n_fit_points

E_initial
E_inf
A

tau_energy
tau_energy_fit_se
tau_energy_window_std
tau_energy_total_se
tau_energy_ci_low
tau_energy_ci_high

tau_shape_equivalent
tau_1e
tau_crosscheck_error

tau_bulk_scaled
tau_shape_bulk_scaled

r_squared
rmse
normalized_rmse

plateau_start_end_index
plateau_end_end_index
n_plateau_windows

multimode_or_nonlinear
fit_hit_bounds
notes
```

Use clear column names and document units as nondimensional.

---

## Step 17: Quality Flags

Assign one of these flags:

### `good`

Use when:

- a stable endpoint plateau exists;
- the fit converges;
- fit quality is high;
- uncertainty is moderate;
- `tau_1e` agrees reasonably with `tau_energy`;
- no oscillation intrudes into the fitting interval.

### `acceptable`

Use when:

- a fit is usable but the plateau is short;
- local rate shows mild drift;
- uncertainty is larger;
- the trace may be weakly multimode.

### `poor`

Use when:

- no stable plateau exists;
- fewer than the minimum number of clean points are available;
- oscillations begin almost immediately;
- the fitted equilibrium is poorly constrained;
- timescale uncertainty is comparable to the estimate.

### `failed`

Use when:

- data are invalid;
- optimization repeatedly fails;
- no physically admissible fit is obtained.

Do not omit poor or failed traces from the summary table. Exclude them from the primary scaling regression by default, but show them visibly in diagnostic plots.

---

## Step 18: Robustness Checks

Implement or make easy to run the following sensitivity checks:

1. Repeat with candidate endpoint limits `10–30`, `10–35`, and `10–40`.
2. Repeat with smoothing windows 3, 5, and 7.
3. Repeat with start indices 0 through 5.
4. Compare robust and ordinary least-squares fits.
5. Compare nonlinear exponential fit with the `1/e` threshold estimate.
6. Check whether the fitted scaling changes when `acceptable` traces are excluded.
7. Check timestep-convergence results if multiple temporal resolutions exist.

Summarize how much the final `Sd` scaling conclusions change under these choices.

---

## Step 19: Interpretation Rules

Use the following interpretation.

### Evidence for Surface-Viscous Scaling

Conclude that the high-`Sd` region is surface-time controlled when:

- normalized traces collapse in the original surface-scaled time;
- `tau_energy` approaches an `Sd`-independent plateau;
- the high-`Sd` log-log slope is statistically consistent with zero.

### Evidence for Bulk-Viscous Scaling

Conclude that the low-`Sd` region is bulk-time controlled when:

- normalized traces collapse under `Sd * t`, including any configured convention factor;
- `tau_energy` is proportional to `Sd**(-1)`;
- `Sd * tau_energy` approaches a plateau;
- the low-`Sd` log-log slope is statistically consistent with `-1`.

### Evidence Against a Clean Single-Timescale Picture

Flag this when:

- semilog excess-energy curves are strongly curved;
- the local decay rate has no plateau;
- the extracted timescale depends strongly on the fitting endpoint;
- changing `Sd` also changes the dominant spatial deformation mode.

In this case, still report the consistently defined early-time effective timescale, but do not describe it as a unique eigenmode relaxation time.

---

## Step 20: Treatment of Late Oscillations

The later oscillations are not part of the timescale fit.

Never:

- estimate `E_inf` from the final 300-step tail by default;
- fit one exponential over the entire trace;
- let a low residual over the late tail determine the early relaxation timescale;
- use thresholds near the final oscillatory energy level.

Do:

- retain the full trace in diagnostic plots;
- record the detected oscillation-onset index;
- optionally compute late-time oscillation amplitude as a numerical-quality metric;
- keep the fitting and oscillation-analysis modules separate.

An optional oscillation metric is:

```python
late_amplitude = percentile(E_late, 95) - percentile(E_late, 5)
```

where `E_late` begins after the detected oscillation onset or after sample 40.

This metric is diagnostic only and must not modify the extracted early relaxation timescale.

---

## Suggested Function Structure

Adapt names to the existing codebase, but aim for a structure similar to:

```python
def validate_trace(t, E) -> ValidationResult:
    ...

def smooth_for_detection(t, E, config) -> np.ndarray:
    ...

def choose_fit_start(t, E, E_smooth, config) -> int:
    ...

def detect_oscillation_start(t, E_smooth, config) -> int:
    ...

def fit_single_exponential(t, E, start_idx, end_idx, config) -> FitResult:
    ...

def scan_fit_windows(t, E, start_idx, oscillation_idx, config) -> list[FitResult]:
    ...

def select_timescale_plateau(fits, config) -> PlateauResult:
    ...

def compute_local_decay_rate(t, E, E_inf, config):
    ...

def compute_threshold_timescale(t, E, E_inf, start_idx):
    ...

def analyze_trace(simulation_id, Sd, t, E, config) -> TraceAnalysis:
    ...

def fit_sd_crossover(summary_df, config) -> ScalingFit:
    ...

def make_trace_diagnostics(analysis, output_dir):
    ...

def make_batch_scaling_plots(summary_df, scaling_fit, output_dir, config):
    ...

def write_summary(summary_df, output_dir):
    ...
```

Use dataclasses or typed dictionaries for fit results so rejected fits still carry a reason.

---

## Testing Requirements

Add synthetic tests.

### Test 1: Clean Single Exponential

Generate:

\[
E(t)=E_\infty+A e^{-t/\tau}
\]

with small noise.

Confirm that the recovered `tau_energy` is within approximately 5–10% of the truth.

### Test 2: Exponential Plus Late Oscillations

Generate a clean exponential for the first 25 samples and add oscillations after sample 30.

Confirm that:

- the fit endpoint remains before the oscillations;
- the recovered timescale remains accurate;
- fitting all 300 points would fail or be biased, demonstrating why the early-window method is necessary.

### Test 3: Double Exponential

Generate a two-mode trace.

Confirm that:

- the single-exponential result is stable only over a limited interval or is flagged as multimode;
- the code still returns an effective timescale and quality warning.

### Test 4: Expected Sd Scaling

Generate synthetic timescales from

\[
\tau(Sd)=a+\frac{b}{Sd}
\]

and synthetic energy traces using those timescales.

Confirm that the batch fit recovers `a` and `b`, the low-`Sd` slope approaches `-1`, and the high-`Sd` slope approaches `0`.

### Test 5: Irregular Time Grid

Generate the same trace on a nonuniform increasing time array.

Confirm correct recovery using actual time values rather than sample indices.

---

## Final Deliverables

Produce:

1. A reusable analysis module integrated with the known data-loading code.
2. A configuration section containing all fitting and scaling choices.
3. One summary CSV or parquet file.
4. One diagnostic figure per simulation.
5. Batch collapse and scaling figures.
6. A short text or Markdown report containing:
   - number of good, acceptable, poor, and failed traces;
   - fitted crossover parameters `a` and `b`;
   - low- and high-`Sd` log-log slopes;
   - whether surface-time collapse is observed at high `Sd`;
   - whether bulk-time collapse is observed at low `Sd`;
   - important sensitivity-analysis results;
   - a list of traces requiring manual inspection.

---

## Minimum Acceptable Implementation

At minimum, the implementation must:

- fit only the early portion of each 300-step trace;
- fit `E_inf`, `A`, and `tau_energy` simultaneously;
- scan candidate endpoints primarily between samples 10 and 35;
- select a stable endpoint plateau rather than one arbitrary fit;
- report fitting-window uncertainty;
- compare original surface time against `Sd`-rescaled bulk time;
- plot `tau_energy` versus `Sd`;
- plot `Sd * tau_energy` versus `Sd`;
- fit or assess the expected asymptotes:
  - high `Sd`: constant `tau_energy`;
  - low `Sd`: `tau_energy ~ 1/Sd`;
- save all extracted values and quality flags in a table.

The analysis must be reproducible, parameterized, and conservative about claiming a unique physical timescale when the early trace is not well described by a single exponential.
