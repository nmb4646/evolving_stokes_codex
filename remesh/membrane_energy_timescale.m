function [summary, scaling, analyses] = membrane_energy_timescale(mode, config_overrides)
%MEMBRANE_ENERGY_TIMESCALE Extract early-time energy relaxation timescales.
%
%   membrane_energy_timescale
%       Processes the fs_batch runs selected in default_config(), writes a
%       summary CSV and Markdown report, and saves diagnostic figures.
%
%   membrane_energy_timescale("test")
%       Runs the built-in synthetic tests without reading simulation data.
%
%   membrane_energy_timescale("analyze", overrides)
%       Applies fields from the scalar struct overrides after default_config().
%
% The fitted model is
%
%   E(t) = E_inf + A * exp(-(t - t0) / tau_energy).
%
% Fits always use the original energy values. Smoothing is used only for
% oscillation detection and local-rate diagnostics.

if nargin < 1
    mode = "analyze";
end
if nargin < 2
    config_overrides = struct();
end

cfg = default_config();
cfg = apply_config_overrides(cfg, config_overrides);
mode = lower(string(mode));

if mode == "test"
    run_synthetic_tests(cfg);
    summary = table();
    scaling = empty_scaling_result();
    analyses = {};
    return
elseif mode ~= "analyze"
    error("Unknown mode '%s'. Use 'analyze' or 'test'.", mode);
end

script_dir = fileparts(mfilename("fullpath"));
remesh_dir = find_remesh_dir(script_dir);
addpath(remesh_dir);

if strlength(cfg.data_dir) == 0
    data_dir = find_data_dir(script_dir, remesh_dir);
else
    data_dir = char(cfg.data_dir);
end

if strlength(cfg.output_dir) == 0
    output_dir = fullfile(remesh_dir, "data", "membrane_energy_timescale");
else
    output_dir = char(cfg.output_dir);
end
diagnostic_dir = fullfile(output_dir, "trace_diagnostics");
if ~isfolder(output_dir)
    mkdir(output_dir);
end
if cfg.save_trace_diagnostics && ~isfolder(diagnostic_dir)
    mkdir(diagnostic_dir);
end

Sd_values = cfg.Sd_values(:).';
analyses = cell(numel(Sd_values), 1);

fprintf("Membrane energy timescale analysis: %d Sd values\n", numel(Sd_values));
for i = 1:numel(Sd_values)
    Sd = Sd_values(i);
    simulation_id = make_run_tag(Sd, cfg.Da, cfg.gamy, cfg.v);
    folder = resolve_run_folder(data_dir, simulation_id, Sd, cfg);
    fprintf("[%d/%d] %s\n", i, numel(Sd_values), simulation_id);

    if strlength(folder) == 0
        analyses{i} = failed_analysis(simulation_id, Sd, ...
            "Run folder was not found under " + string(data_dir));
        continue
    end

    try
        trace = load_energy_trace(char(folder), cfg);
        analyses{i} = analyze_trace(simulation_id, Sd, trace.t, trace.E, ...
            trace.frame_ids, trace.dt, cfg);
    catch ME
        analyses{i} = failed_analysis(simulation_id, Sd, string(ME.message));
    end

    a = analyses{i};
    if a.status == "ok"
        fprintf("  tau_E = %.6g, quality = %s, fit samples %d:%d\n", ...
            a.tau_energy, a.quality_flag, a.fit_start_index, a.fit_end_index);
        if cfg.save_trace_diagnostics
            make_trace_diagnostic(a, diagnostic_dir, cfg);
        end
    else
        fprintf("  failed: %s\n", a.notes);
    end
end

summary = make_summary_table(analyses, cfg);
summary = sortrows(summary, "Sd");
scaling = fit_sd_scaling(summary, cfg);

summary_file = fullfile(output_dir, "membrane_energy_timescale_summary.csv");
writetable(summary, summary_file);
fprintf("Saved %s\n", summary_file);

if cfg.save_batch_plots
    make_batch_plots(analyses, summary, scaling, output_dir, cfg);
end

sensitivity = table();
if cfg.run_sensitivity_analysis
    sensitivity = run_sensitivity_analysis(analyses, summary, cfg);
    sensitivity_file = fullfile(output_dir, "membrane_energy_timescale_sensitivity.csv");
    writetable(sensitivity, sensitivity_file);
    fprintf("Saved %s\n", sensitivity_file);
end

write_report(summary, scaling, sensitivity, output_dir, cfg);
disp(summary(:, ["Sd", "status", "quality_flag", "tau_energy", ...
    "tau_energy_total_se", "tau_bulk_scaled", "fit_start_index", ...
    "fit_end_index", "oscillation_start_index"]));
end

function cfg = apply_config_overrides(cfg, overrides)
    if isempty(overrides)
        return
    end
    if ~isstruct(overrides) || ~isscalar(overrides)
        error("config_overrides must be a scalar struct.");
    end
    names = fieldnames(overrides);
    for i = 1:numel(names)
        name = names{i};
        if ~isfield(cfg, name)
            error("Unknown configuration field '%s'.", name);
        end
        cfg.(name) = overrides.(name);
    end
end

function cfg = default_config()
    %%% Simulation selection
    cfg.Sd_values = [3e-3,1e-2, 3e-2, 1e-1, 3e-1, 1, 3, 1e1, 3e1, 1e2];
    %cfg.Sd_values = [1, 3, 1e1, 3e1, 1e2];
    cfg.Da = 0;
    cfg.gamy = 0;
    cfg.v = 0.97;
    cfg.energy_Kb = 1;
    cfg.data_dir = "";   % Empty finds remesh/data/fs_batch_data automatically.
    cfg.output_dir = ""; % Empty writes remesh/data/membrane_energy_timescale.
    cfg.allow_legacy_folder_without_v = true;

    %%% Input validation and early-time search
    cfg.expected_n_samples = 300;
    cfg.allow_one_extra_initial_sample = true; % Allows geo0 through geo300.
    cfg.require_expected_n_samples = false;
    cfg.start_index_min = 0;
    cfg.start_index_max = 5;
    cfg.min_usable_points = 12;
    cfg.min_fit_points = 8;
    cfg.max_early_index = 200;
    cfg.candidate_end_min = 10;
    cfg.candidate_end_max = 35;
    cfg.candidate_end_step = 1;

    %%% Detection settings
    cfg.smooth_window = 5;
    cfg.smooth_polyorder = 2;
    cfg.persistent_rise_points = 3;
    cfg.derivative_noise_multiplier = 3.0;
    cfg.sign_change_window = 7;
    cfg.sign_changes_required = 2;
    cfg.slow_decay_fraction = 0.20;
    cfg.startup_jump_multiplier = 8.0;

    %%% Nonlinear fit and candidate rejection
    cfg.robust_loss = "soft_l1"; % "soft_l1" or "linear".
    cfg.f_scale = 1.0;
    cfg.energy_bound_margin_fraction = 0.50;
    cfg.tau_lower_step_fraction = 0.05;
    cfg.tau_upper_interval_factor = 100;
    cfg.min_r_squared = 0.98;
    cfg.max_relative_tau_se = 0.50;
    cfg.max_normalized_rmse = 0.10;
    cfg.boundary_tolerance = 2e-3;
    cfg.max_tau_interval_ratio = 10;
    cfg.min_tau_step_ratio = 0.10;
    cfg.optimizer_max_iterations = 2000;
    cfg.optimizer_max_evaluations = 8000;

    %%% Endpoint plateau and quality settings
    cfg.plateau_relative_tolerance = 0.15;
    cfg.min_plateau_windows = 4;
    cfg.crosscheck_good_tolerance = 0.30;
    cfg.local_rate_drift_tolerance = 0.50;
    cfg.double_exponential_aic_improvement = 6.0;
    cfg.double_exponential_min_amplitude_fraction = 0.05;
    cfg.double_exponential_min_tau_ratio = 1.8;

    %%% Scaling analysis
    cfg.bulk_time_factor = 1.0;
    cfg.include_acceptable_in_scaling = true;
    cfg.minimum_scaling_points = 3;

    %%% Outputs and visual settings
    cfg.save_trace_diagnostics = true;
    cfg.save_batch_plots = true;
    cfg.batch_timescale_only = false; % True saves only the tau_E-versus-Sd scaling plot.
    cfg.figure_visible = "off"; % Use "on" to display figures while processing.
    cfg.figure_position = [100, 100, 1250, 780];
    cfg.line_width = 1.7;
    cfg.marker_size = 38;
    cfg.axes_font_size = 11;
    cfg.output_resolution = 220;

    %%% Optional robustness checks
    % This reuses the already-loaded energy traces but performs many extra fits.
    cfg.run_sensitivity_analysis = false;
    cfg.sensitivity_end_max = [30, 35, 40];
    cfg.sensitivity_smooth_windows = [3, 5, 7];
    cfg.sensitivity_losses = ["soft_l1", "linear"];
end

function trace = load_energy_trace(folder, cfg)
    files = dir(fullfile(folder, "geo*.mat"));
    if isempty(files)
        error("No geo*.mat files found in %s.", folder);
    end

    frame_ids = NaN(numel(files), 1);
    for i = 1:numel(files)
        token = regexp(files(i).name, "^geo(\d+)\.mat$", "tokens", "once");
        if ~isempty(token)
            frame_ids(i) = str2double(token{1});
        end
    end
    valid = isfinite(frame_ids);
    files = files(valid);
    frame_ids = frame_ids(valid);
    [frame_ids, order] = sort(frame_ids);
    files = files(order);

    if isempty(frame_ids)
        error("No files matching geo<number>.mat were found in %s.", folder);
    end

    E = NaN(numel(files), 1);
    saved_time = NaN(numel(files), 1);
    dt = NaN;

    for i = 1:numel(files)
        data = load(fullfile(folder, files(i).name), "M", "P", "p");
        if ~isfield(data, "M") || ~isfield(data, "P")
            error("%s does not contain both M and P.", files(i).name);
        end
        if isfield(data, "p")
            if isnan(dt) && isfield(data.p, "dt") && isscalar(data.p.dt)
                dt = double(data.p.dt);
            end
            saved_time(i) = extract_saved_time(data.p);
        end
        geo = Geometry(data.M, data.P);
        E(i) = geo.willmore_energy(cfg.energy_Kb);
    end

    if all(isfinite(saved_time)) && all(diff(saved_time) > 0)
        t = saved_time;
    else
        if ~isfinite(dt) || dt <= 0
            error("A positive p.dt was not found in %s.", folder);
        end
        t = frame_ids * dt;
    end

    trace = struct("t", t(:), "E", E(:), "frame_ids", frame_ids(:), "dt", dt);
end

function value = extract_saved_time(p)
    value = NaN;
    names = ["current_time", "simulation_time", "time_current"];
    for name = names
        if isfield(p, name) && isscalar(p.(name)) && isfinite(p.(name))
            value = double(p.(name));
            return
        end
    end
end

function analysis = analyze_trace(simulation_id, Sd, t, E, frame_ids, dt, cfg)
    analysis = failed_analysis(simulation_id, Sd, "Analysis did not complete.");
    analysis.t = t(:);
    analysis.E = E(:);
    analysis.frame_ids = frame_ids(:);
    analysis.dt = dt;

    [valid, validation_notes] = validate_trace(t, E, cfg);
    if ~valid
        analysis.notes = validation_notes;
        return
    end

    E_smooth = smooth_for_detection(t, E, cfg);
    oscillation_pos = detect_oscillation_start(t, E_smooth, cfg);
    [start_pos, fits, plateau, start_notes] = choose_fit_start( ...
        t, E, E_smooth, oscillation_pos, cfg);

    if isempty(fits)
        analysis.notes = join_notes(validation_notes, start_notes, ...
            "No candidate fitting windows were available.");
        return
    end

    admissible = [fits.admissible];
    if ~any(admissible)
        reasons = unique(string({fits.rejection_reason}));
        analysis.notes = join_notes(validation_notes, start_notes, ...
            "No physically admissible exponential fit. " + strjoin(reasons, "; "));
        analysis.candidate_fits = fits;
        analysis.E_smooth = E_smooth;
        return
    end

    if plateau.exists
        plateau_fits = fits(plateau.fit_indices);
        plateau_tau = [plateau_fits.tau];
        tau_energy = median(plateau_tau);
        [~, nearest] = min(abs(plateau_tau - tau_energy));
        selected_fit_index = plateau.fit_indices(nearest);
    else
        usable_indices = find(admissible);
        [~, best_local] = max([fits(usable_indices).selection_score]);
        selected_fit_index = usable_indices(best_local);
        plateau_tau = fits(selected_fit_index).tau;
        tau_energy = fits(selected_fit_index).tau;
    end

    selected = fits(selected_fit_index);
    tau_window_low = local_percentile(plateau_tau, 16);
    tau_window_high = local_percentile(plateau_tau, 84);
    tau_window_std = robust_std(plateau_tau);
    if isscalar(plateau_tau)
        tau_window_std = NaN;
    end
    fit_se = selected.tau_se;
    total_se = combine_uncertainties(fit_se, tau_window_std);

    [local_time, local_rate, local_rate_drift, rate_multimode] = ...
        compute_local_decay_rate(t, E_smooth, selected.E_inf, ...
        start_pos, selected.end_pos, tau_energy, cfg);
    [double_support, double_delta_aic] = double_exponential_diagnostic( ...
        t, E, selected, cfg);
    multimode = rate_multimode || double_support;
    [tau_1e, crosscheck_error] = compute_threshold_timescale( ...
        t, E, selected.E_inf, start_pos, oscillation_pos, tau_energy);

    quality = assign_quality_flag(selected, plateau, total_se, ...
        crosscheck_error, multimode, oscillation_pos, cfg);

    late_start = min(max(oscillation_pos, cfg.max_early_index + 1), numel(E));
    E_late = E(late_start:end);
    late_amplitude = local_percentile(E_late, 95) - local_percentile(E_late, 5);

    analysis.status = "ok";
    analysis.quality_flag = quality;
    analysis.notes = join_notes(validation_notes, start_notes, selected.rejection_reason);
    analysis.n_total = numel(E);
    analysis.fit_start_index = start_pos - 1;
    analysis.fit_end_index = selected.end_pos - 1;
    analysis.oscillation_start_index = oscillation_pos - 1;
    analysis.n_fit_points = selected.n_points;
    analysis.E_initial = E(start_pos);
    analysis.E_inf = selected.E_inf;
    analysis.A = selected.A;
    analysis.tau_energy = tau_energy;
    analysis.tau_energy_fit_se = fit_se;
    analysis.tau_energy_window_std = tau_window_std;
    analysis.tau_energy_total_se = total_se;
    analysis.tau_energy_ci_low = max(0, tau_energy - 1.96 * total_se);
    analysis.tau_energy_ci_high = tau_energy + 1.96 * total_se;
    if ~isfinite(total_se)
        analysis.tau_energy_ci_low = tau_window_low;
        analysis.tau_energy_ci_high = tau_window_high;
    end
    analysis.tau_window_low = tau_window_low;
    analysis.tau_window_high = tau_window_high;
    analysis.tau_shape_equivalent = 2 * tau_energy;
    analysis.tau_1e = tau_1e;
    analysis.tau_crosscheck_error = crosscheck_error;
    analysis.r_squared = selected.r_squared;
    analysis.rmse = selected.rmse;
    analysis.normalized_rmse = selected.normalized_rmse;
    analysis.aic = selected.aic;
    analysis.plateau_start_end_index = plateau.start_end_index;
    analysis.plateau_end_end_index = plateau.end_end_index;
    analysis.n_plateau_windows = numel(plateau.fit_indices);
    analysis.multimode_or_nonlinear = multimode;
    analysis.fit_hit_bounds = selected.hit_bounds;
    analysis.late_oscillation_amplitude = late_amplitude;
    analysis.E_smooth = E_smooth;
    analysis.candidate_fits = fits;
    analysis.plateau_fit_indices = plateau.fit_indices;
    analysis.selected_fit_index = selected_fit_index;
    analysis.local_rate_time = local_time;
    analysis.local_decay_rate = local_rate;
    analysis.local_rate_drift = local_rate_drift;
    analysis.double_exponential_delta_aic = double_delta_aic;
end

function [valid, notes] = validate_trace(t, E, cfg)
    valid = false;
    notes = "";
    if numel(t) ~= numel(E)
        notes = "Time and energy arrays have different lengths.";
        return
    end
    if numel(t) < cfg.min_usable_points
        notes = sprintf("Only %d samples are available; at least %d are required.", ...
            numel(t), cfg.min_usable_points);
        return
    end
    if any(~isfinite(t)) || any(~isfinite(E))
        notes = "Time or energy contains a nonfinite value.";
        return
    end
    if any(diff(t) <= 0)
        notes = "Time is not strictly increasing.";
        return
    end
    if max(E) - min(E) <= 10 * eps(max(abs(E)))
        notes = "Energy range is numerically zero.";
        return
    end

    expected = cfg.expected_n_samples;
    expected_ok = numel(E) == expected;
    if cfg.allow_one_extra_initial_sample
        expected_ok = expected_ok || numel(E) == expected + 1;
    end
    if ~expected_ok
        notes = sprintf("Expected %d samples%s but found %d.", expected, ...
            extra_sample_text(cfg.allow_one_extra_initial_sample), numel(E));
        if cfg.require_expected_n_samples
            return
        end
    end
    valid = true;
end

function text = extra_sample_text(allow_extra)
    if allow_extra
        text = " (or one extra initial snapshot)";
    else
        text = "";
    end
end

function E_smooth = smooth_for_detection(t, E, cfg)
    n = numel(E);
    window = min(cfg.smooth_window, n);
    if mod(window, 2) == 0
        window = window - 1;
    end
    if window < 3
        E_smooth = E;
        return
    end

    polyorder = min(cfg.smooth_polyorder, window - 1);
    half = floor(window / 2);
    E_smooth = zeros(size(E));
    for i = 1:n
        first = max(1, i - half);
        last = min(n, i + half);
        idx = first:last;
        if numel(idx) <= polyorder
            E_smooth(i) = mean(E(idx));
            continue
        end
        tc = t(idx) - t(i);
        coefficients = polyfit(tc, E(idx), polyorder);
        E_smooth(i) = polyval(coefficients, 0);
    end
end

function onset = detect_oscillation_start(t, E_smooth, cfg)
    n = numel(E_smooth);
    default_pos = min(cfg.max_early_index + 1, n);
    dE_dt = gradient(E_smooth, t);
    early_last = min(cfg.max_early_index + 1, n);
    d_early = dE_dt(1:early_last);
    sigma_d = 1.4826 * median(abs(d_early - median(d_early)));
    sigma_d = max(sigma_d, eps(max(abs(d_early))));
    threshold = cfg.derivative_noise_multiplier * sigma_d;
    search_first = min(cfg.start_index_max + 2, n);

    onset_a = inf;
    run_length = cfg.persistent_rise_points;
    for i = search_first:max(search_first, early_last - run_length + 1)
        idx = i:min(i + run_length - 1, n);
        if numel(idx) == run_length && all(dE_dt(idx) > threshold)
            onset_a = i;
            break
        end
    end

    onset_b = inf;
    max_initial_drop = max(abs(dE_dt(1:min(6, n))));
    rolling = cfg.sign_change_window;
    for i = max(search_first, 7):max(search_first, early_last - rolling + 1)
        idx = i:min(i + rolling - 1, n);
        values = dE_dt(idx);
        if numel(values) < rolling
            continue
        end
        slowed = median(abs(values)) <= cfg.slow_decay_fraction * max_initial_drop;
        signs = sign(values);
        signs(abs(values) <= sigma_d) = 0;
        signs = signs(signs ~= 0);
        changes = 0;
        if numel(signs) >= 3
            changes = nnz(diff(signs) ~= 0);
        end
        if slowed && changes >= cfg.sign_changes_required
            onset_b = i;
            break
        end
    end

    onset = min([default_pos, onset_a, onset_b]);
    onset = max(onset, cfg.min_fit_points + 1);
end

function [start_pos, fits, plateau, notes] = choose_fit_start(t, E, E_smooth, oscillation_pos, cfg)
    start_positions = (cfg.start_index_min:cfg.start_index_max) + 1;
    start_positions = start_positions(start_positions <= numel(E) - cfg.min_fit_points + 1);
    best = struct("start_pos", NaN, "fits", [], "plateau", empty_plateau(), ...
        "score", -inf);
    notes = "";

    for start_pos = start_positions
        if ~initially_decreasing(E_smooth, start_pos)
            continue
        end
        if is_startup_artifact(E_smooth, start_pos, cfg)
            continue
        end

        candidate_fits = scan_fit_windows(t, E, start_pos, oscillation_pos, cfg);
        candidate_plateau = select_timescale_plateau(candidate_fits, cfg);
        score = start_selection_score(candidate_fits, candidate_plateau, start_pos);
        if score > best.score
            best = struct("start_pos", start_pos, "fits", candidate_fits, ...
                "plateau", candidate_plateau, "score", score);
        end
        if candidate_plateau.exists
            break
        end
    end

    if ~isfinite(best.score)
        for start_pos = start_positions
            candidate_fits = scan_fit_windows(t, E, start_pos, oscillation_pos, cfg);
            candidate_plateau = select_timescale_plateau(candidate_fits, cfg);
            score = start_selection_score(candidate_fits, candidate_plateau, start_pos);
            if score > best.score
                best = struct("start_pos", start_pos, "fits", candidate_fits, ...
                    "plateau", candidate_plateau, "score", score);
            end
        end
        notes = "No start candidate passed the startup checks; used the best admissible fallback.";
    elseif best.start_pos > cfg.start_index_min + 1
        notes = sprintf("Skipped %d startup sample(s).", best.start_pos - 1);
    end

    start_pos = best.start_pos;
    fits = best.fits;
    plateau = best.plateau;
end

function tf = initially_decreasing(E, start_pos)
    last = min(numel(E), start_pos + 6);
    tf = E(start_pos) > median(E(max(start_pos + 2, start_pos):last));
end

function tf = is_startup_artifact(E, start_pos, cfg)
    tf = false;
    if start_pos + 5 > numel(E)
        return
    end
    jumps = abs(diff(E(start_pos:start_pos + 5)));
    reference = median(jumps(2:end));
    if reference <= eps(max(abs(E)))
        return
    end
    first_jump_isolated = jumps(1) > cfg.startup_jump_multiplier * reference;
    reverses_immediately = sign(E(start_pos + 1) - E(start_pos)) ~= ...
        sign(median(diff(E(start_pos + 1:start_pos + 5))));
    tf = first_jump_isolated && reverses_immediately;
end

function score = start_selection_score(fits, plateau, start_pos)
    if isempty(fits)
        score = -inf;
        return
    end
    score = 1000 * plateau.exists + 20 * numel(plateau.fit_indices) ...
        + 2 * nnz([fits.valid]) + nnz([fits.admissible]) - 0.01 * start_pos;
end

function fits = scan_fit_windows(t, E, start_pos, oscillation_pos, cfg)
    start_index = start_pos - 1;
    lower_end = max(start_index + cfg.min_fit_points - 1, cfg.candidate_end_min);
    upper_end = min([cfg.candidate_end_max, oscillation_pos - 2, numel(E) - 1]);
    endpoints = lower_end:cfg.candidate_end_step:upper_end;

    if numel(endpoints) < cfg.min_plateau_windows
        upper_relaxed = min(cfg.max_early_index, numel(E) - 1);
        endpoints = lower_end:cfg.candidate_end_step:upper_relaxed;
    end
    if isempty(endpoints)
        fits = repmat(empty_fit(), 0, 1);
        return
    end

    fits = repmat(empty_fit(), numel(endpoints), 1);
    for i = 1:numel(endpoints)
        fits(i) = fit_single_exponential(t, E, start_pos, endpoints(i) + 1, cfg);
        fits(i).extends_past_oscillation = endpoints(i) >= oscillation_pos - 1;
    end
end

function fit = fit_single_exponential(t, E, start_pos, end_pos, cfg)
    fit = empty_fit();
    fit.start_pos = start_pos;
    fit.end_pos = end_pos;
    fit.start_index = start_pos - 1;
    fit.end_index = end_pos - 1;
    fit.n_points = end_pos - start_pos + 1;

    if fit.n_points < cfg.min_fit_points
        fit.rejection_reason = "too few points";
        return
    end

    tw = t(start_pos:end_pos);
    Ew = E(start_pos:end_pos);
    t0 = tw(1);
    elapsed = tw - t0;
    duration = elapsed(end);
    step = median(diff(tw));
    energy_range = max(Ew) - min(Ew);
    if duration <= 0 || step <= 0 || energy_range <= eps(max(abs(Ew)))
        fit.rejection_reason = "degenerate fit interval";
        return
    end

    margin = max(cfg.energy_bound_margin_fraction * energy_range, ...
        100 * eps(max(abs(Ew))));
    E_lower = min(Ew) - margin;
    E_upper = min(Ew(end) + margin, Ew(1) - 10 * eps(max(abs(Ew))));
    if E_upper <= E_lower
        E_upper = E_lower + max(energy_range, margin);
    end
    A_lower = max(100 * eps(max(abs(Ew))), 1e-14);
    A_upper = max(10 * energy_range + abs(Ew(1) - E_lower), 100 * A_lower);
    tau_lower = max(cfg.tau_lower_step_fraction * min(diff(tw)), eps(duration));
    tau_upper = cfg.tau_upper_interval_factor * duration;
    bounds = [E_lower, E_upper; A_lower, A_upper; tau_lower, tau_upper];

    E_inf_guess = median(Ew(max(1, end - 2):end));
    E_inf_guess = clamp_inside(E_inf_guess, E_lower, E_upper);
    A_guess = clamp_inside(max(Ew(1) - E_inf_guess, A_lower), A_lower, A_upper);
    tau_guesses = duration * [0.15, 0.25, 0.50];
    tau_guesses = min(max(tau_guesses, tau_lower), tau_upper);
    scale = max(energy_range, sqrt(eps) * max(abs(Ew)));

    options = optimset("Display", "off", ...
        "MaxIter", cfg.optimizer_max_iterations, ...
        "MaxFunEvals", cfg.optimizer_max_evaluations, ...
        "TolX", 1e-10, "TolFun", 1e-12);
    best_objective = inf;
    best_parameters = [NaN, NaN, NaN];
    best_exitflag = -1;

    for tau_guess = tau_guesses
        initial = [E_inf_guess, A_guess, tau_guess];
        z0 = physical_to_unconstrained(initial, bounds);
        objective = @(z) exponential_objective(z, bounds, elapsed, Ew, scale, cfg);
        [z, objective_value, exitflag] = fminsearch(objective, z0, options);
        parameters = unconstrained_to_physical(z, bounds);
        if isfinite(objective_value) && objective_value < best_objective
            best_objective = objective_value;
            best_parameters = parameters;
            best_exitflag = exitflag;
        end
    end

    fit.exitflag = best_exitflag;
    fit.converged = best_exitflag > 0 && all(isfinite(best_parameters));
    if ~fit.converged
        fit.rejection_reason = "optimizer did not converge";
        return
    end

    fit.E_inf = best_parameters(1);
    fit.A = best_parameters(2);
    fit.tau = best_parameters(3);
    prediction = fit.E_inf + fit.A * exp(-elapsed / fit.tau);
    residual = Ew - prediction;
    fit.rss = sum(residual .^ 2);
    fit.rmse = sqrt(mean(residual .^ 2));
    fit.normalized_rmse = fit.rmse / max(energy_range, eps(max(abs(Ew))));
    total_sum_squares = sum((Ew - mean(Ew)) .^ 2);
    fit.r_squared = 1 - fit.rss / max(total_sum_squares, eps(total_sum_squares));
    fit.aic = fit.n_points * log(max(fit.rss / fit.n_points, realmin)) + 6;

    exponential = exp(-elapsed / fit.tau);
    J = [ones(fit.n_points, 1), exponential, ...
        fit.A * exponential .* elapsed / fit.tau ^ 2];
    dof = fit.n_points - 3;
    if dof > 0 && rank(J) == 3
        covariance = (fit.rss / dof) * pinv(J' * J);
        parameter_se = sqrt(max(diag(covariance), 0));
        fit.E_inf_se = parameter_se(1);
        fit.A_se = parameter_se(2);
        fit.tau_se = parameter_se(3);
    end

    spans = bounds(:, 2) - bounds(:, 1);
    relative_to_lower = (best_parameters(:) - bounds(:, 1)) ./ spans;
    relative_to_upper = (bounds(:, 2) - best_parameters(:)) ./ spans;
    fit.hit_bounds = any(relative_to_lower < cfg.boundary_tolerance ...
        | relative_to_upper < cfg.boundary_tolerance);

    relative_tau_se = fit.tau_se / fit.tau;
    if ~isfinite(relative_tau_se)
        relative_tau_se = inf;
    end
    residual_correlation = lag_one_correlation(residual);
    fit.systematic_residual = isfinite(residual_correlation) ...
        && abs(residual_correlation) > 0.90;

    reasons = strings(0, 1);
    if fit.A <= 0 || fit.tau <= 0
        reasons(end + 1) = "nonpositive amplitude or timescale";
    end
    if fit.hit_bounds
        reasons(end + 1) = "parameter hit a bound";
    end
    if fit.E_inf >= Ew(1)
        reasons(end + 1) = "equilibrium energy is incompatible with decay";
    end
    if fit.r_squared < cfg.min_r_squared
        reasons(end + 1) = "R-squared below threshold";
    end
    if relative_tau_se > cfg.max_relative_tau_se
        reasons(end + 1) = "relative tau standard error too large";
    end
    if fit.normalized_rmse > cfg.max_normalized_rmse
        reasons(end + 1) = "normalized RMSE too large";
    end
    if fit.tau < cfg.min_tau_step_ratio * step
        reasons(end + 1) = "timescale is below one-tenth of a timestep";
    end
    if fit.tau > cfg.max_tau_interval_ratio * duration && fit.r_squared < 0.995
        reasons(end + 1) = "timescale is much longer than the fit interval";
    end
    if fit.systematic_residual
        reasons(end + 1) = "strongly correlated residuals";
    end

    fit.admissible = fit.converged && fit.A > 0 && fit.tau > 0 ...
        && fit.E_inf < Ew(1) && ~fit.hit_bounds;
    fit.valid = fit.admissible && isempty(reasons);
    fit.rejection_reason = strjoin(reasons, "; ");
    fit.selection_score = fit.r_squared - fit.normalized_rmse ...
        - 0.05 * min(relative_tau_se, 10) - 0.25 * fit.hit_bounds;
end

function objective = exponential_objective(z, bounds, elapsed, E, scale, cfg)
    parameters = unconstrained_to_physical(z, bounds);
    prediction = parameters(1) + parameters(2) * exp(-elapsed / parameters(3));
    residual = (E - prediction) / scale;
    if any(~isfinite(residual))
        objective = realmax / 100;
        return
    end
    if lower(cfg.robust_loss) == "soft_l1"
        f_scale = max(cfg.f_scale, eps);
        objective = sum(2 * f_scale ^ 2 * (sqrt(1 + (residual / f_scale) .^ 2) - 1));
    elseif any(lower(cfg.robust_loss) == ["linear", "least_squares", "ols"])
        objective = sum(residual .^ 2);
    else
        error("Unknown robust_loss '%s'.", cfg.robust_loss);
    end
end

function [supported, delta_aic] = double_exponential_diagnostic(t, E, single, cfg)
    supported = false;
    delta_aic = NaN;
    positions = single.start_pos:single.end_pos;
    tw = t(positions);
    Ew = E(positions);
    elapsed = tw - tw(1);
    duration = elapsed(end);
    step = median(diff(tw));
    energy_range = max(Ew) - min(Ew);
    if numel(Ew) < 10 || duration <= 0 || energy_range <= eps(max(abs(Ew)))
        return
    end

    margin = max(cfg.energy_bound_margin_fraction * energy_range, ...
        100 * eps(max(abs(Ew))));
    E_lower = min(Ew) - margin;
    E_upper = min(Ew(end) + margin, Ew(1) - 10 * eps(max(abs(Ew))));
    amplitude_lower = max(100 * eps(max(abs(Ew))), 1e-14);
    amplitude_upper = max(10 * energy_range, 100 * amplitude_lower);
    tau_lower = max(cfg.tau_lower_step_fraction * min(diff(tw)), eps(duration));
    tau_upper = cfg.tau_upper_interval_factor * duration;
    bounds = [E_lower, E_upper; ...
        amplitude_lower, amplitude_upper; tau_lower, tau_upper; ...
        amplitude_lower, amplitude_upper; tau_lower, tau_upper];

    scale = max(energy_range, sqrt(eps) * max(abs(Ew)));
    options = optimset("Display", "off", "MaxIter", cfg.optimizer_max_iterations, ...
        "MaxFunEvals", cfg.optimizer_max_evaluations, "TolX", 1e-9, "TolFun", 1e-11);
    tau_pairs = [max(tau_lower, single.tau / 3), min(tau_upper, single.tau * 3); ...
        max(tau_lower, single.tau / 2), min(tau_upper, single.tau * 2)];
    best_rss = inf;
    best_parameters = NaN(1, 5);
    for pair = 1:size(tau_pairs, 1)
        E_inf_guess = clamp_inside(single.E_inf, E_lower, E_upper);
        amplitude_guess = clamp_inside(single.A / 2, amplitude_lower, amplitude_upper);
        initial = [E_inf_guess, amplitude_guess, tau_pairs(pair, 1), ...
            amplitude_guess, tau_pairs(pair, 2)];
        z0 = physical_to_unconstrained(initial, bounds);
        objective = @(z) double_exponential_objective(z, bounds, elapsed, Ew, scale, cfg);
        [z, ~, exitflag] = fminsearch(objective, z0, options);
        if exitflag <= 0
            continue
        end
        parameters = unconstrained_to_physical(z, bounds);
        prediction = parameters(1) ...
            + parameters(2) * exp(-elapsed / parameters(3)) ...
            + parameters(4) * exp(-elapsed / parameters(5));
        rss = sum((Ew - prediction) .^ 2);
        if rss < best_rss
            best_rss = rss;
            best_parameters = parameters;
        end
    end
    if ~all(isfinite(best_parameters))
        return
    end

    double_aic = numel(Ew) * log(max(best_rss / numel(Ew), realmin)) + 10;
    delta_aic = single.aic - double_aic;
    amplitude_fraction = [best_parameters(2), best_parameters(4)] / energy_range;
    tau_ratio = max(best_parameters([3, 5])) / min(best_parameters([3, 5]));
    spans = bounds(:, 2) - bounds(:, 1);
    lower_fraction = (best_parameters(:) - bounds(:, 1)) ./ spans;
    upper_fraction = (bounds(:, 2) - best_parameters(:)) ./ spans;
    hit_bounds = any(lower_fraction < cfg.boundary_tolerance ...
        | upper_fraction < cfg.boundary_tolerance);
    supported = delta_aic >= cfg.double_exponential_aic_improvement ...
        && all(amplitude_fraction >= cfg.double_exponential_min_amplitude_fraction) ...
        && tau_ratio >= cfg.double_exponential_min_tau_ratio && ~hit_bounds ...
        && min(best_parameters([3, 5])) >= cfg.min_tau_step_ratio * step;
end

function objective = double_exponential_objective(z, bounds, elapsed, E, scale, cfg)
    parameters = unconstrained_to_physical(z, bounds);
    prediction = parameters(1) ...
        + parameters(2) * exp(-elapsed / parameters(3)) ...
        + parameters(4) * exp(-elapsed / parameters(5));
    residual = (E - prediction) / scale;
    if any(~isfinite(residual))
        objective = realmax / 100;
    elseif lower(cfg.robust_loss) == "soft_l1"
        f_scale = max(cfg.f_scale, eps);
        objective = sum(2 * f_scale ^ 2 * (sqrt(1 + (residual / f_scale) .^ 2) - 1));
    else
        objective = sum(residual .^ 2);
    end
end

function z = physical_to_unconstrained(parameters, bounds)
    fraction = (parameters(:) - bounds(:, 1)) ./ (bounds(:, 2) - bounds(:, 1));
    fraction = min(max(fraction, 1e-8), 1 - 1e-8);
    z = log(fraction ./ (1 - fraction));
end

function parameters = unconstrained_to_physical(z, bounds)
    z = min(max(z(:), -50), 50);
    fraction = 1 ./ (1 + exp(-z));
    parameters = (bounds(:, 1) + fraction .* (bounds(:, 2) - bounds(:, 1))).';
end

function value = clamp_inside(value, lower_bound, upper_bound)
    margin = 1e-6 * (upper_bound - lower_bound);
    value = min(max(value, lower_bound + margin), upper_bound - margin);
end

function value = lag_one_correlation(residual)
    value = NaN;
    if numel(residual) < 5 || std(residual(1:end-1)) == 0 || std(residual(2:end)) == 0
        return
    end
    C = corrcoef(residual(1:end-1), residual(2:end));
    value = C(1, 2);
end

function plateau = select_timescale_plateau(fits, cfg)
    plateau = empty_plateau();
    if isempty(fits)
        return
    end
    valid_indices = find([fits.valid]);
    if numel(valid_indices) < cfg.min_plateau_windows
        return
    end

    best_length = 0;
    best_bonus = -inf;
    best_indices = [];
    for first = 1:numel(valid_indices)
        for last = first:numel(valid_indices)
            indices = valid_indices(first:last);
            endpoints = [fits(indices).end_index];
            if any(diff(endpoints) ~= cfg.candidate_end_step)
                continue
            end
            tau = [fits(indices).tau];
            median_tau = median(tau);
            stable = all(abs(tau - median_tau) / median_tau ...
                <= cfg.plateau_relative_tolerance);
            if ~stable || numel(indices) < cfg.min_plateau_windows
                continue
            end
            expected_overlap = nnz(endpoints >= 20 & endpoints <= 30);
            bonus = 10 * expected_overlap + mean(endpoints) / 100;
            if numel(indices) > best_length ...
                    || (numel(indices) == best_length && bonus > best_bonus)
                best_length = numel(indices);
                best_bonus = bonus;
                best_indices = indices;
            end
        end
    end

    if ~isempty(best_indices)
        plateau.exists = true;
        plateau.fit_indices = best_indices;
        plateau.start_end_index = fits(best_indices(1)).end_index;
        plateau.end_end_index = fits(best_indices(end)).end_index;
    end
end

function [rate_time, rate, drift, multimode] = compute_local_decay_rate( ...
        t, E_smooth, E_inf, start_pos, end_pos, tau, cfg)
    early_end = min([numel(t), cfg.max_early_index + 1, end_pos + 5]);
    idx = start_pos:early_end;
    excess = E_smooth(idx) - E_inf;
    margin = max(1e-10 * max(abs(E_smooth(idx))), eps(max(abs(E_smooth(idx)))));
    keep = excess > margin;
    idx = idx(keep);
    idx = idx(:);
    excess = excess(keep);
    if numel(idx) < 4
        rate_time = [];
        rate = [];
        drift = NaN;
        multimode = true;
        return
    end
    rate_time = t(idx);
    rate = -gradient(log(excess), rate_time);
    fit_keep = idx <= end_pos & isfinite(rate) & rate > 0;
    rate_for_fit = rate(fit_keep);
    time_for_fit = rate_time(fit_keep);
    if numel(rate_for_fit) < 4
        drift = NaN;
        multimode = true;
        return
    end
    coefficients = polyfit(time_for_fit, rate_for_fit, 1);
    drift = abs(coefficients(1)) * (time_for_fit(end) - time_for_fit(1)) / (1 / tau);
    multimode = drift > cfg.local_rate_drift_tolerance;
end

function [tau_1e, error_relative] = compute_threshold_timescale( ...
        t, E, E_inf, start_pos, oscillation_pos, tau)
    tau_1e = NaN;
    error_relative = NaN;
    denominator = E(start_pos) - E_inf;
    if denominator <= 0
        return
    end
    last = min(oscillation_pos, numel(E));
    y = (E(start_pos:last) - E_inf) / denominator;
    crossing = find(y <= exp(-1), 1, "first");
    if isempty(crossing) || crossing == 1
        return
    end
    positions = start_pos:last;
    i1 = positions(crossing - 1);
    i2 = positions(crossing);
    if y(crossing) == y(crossing - 1)
        crossing_time = t(i2);
    else
        fraction = (exp(-1) - y(crossing - 1)) / (y(crossing) - y(crossing - 1));
        crossing_time = t(i1) + fraction * (t(i2) - t(i1));
    end
    tau_1e = crossing_time - t(start_pos);
    error_relative = abs(tau_1e - tau) / tau;
end

function quality = assign_quality_flag(selected, plateau, total_se, ...
        crosscheck_error, multimode, oscillation_pos, cfg)
    relative_total_se = total_se / selected.tau;
    if ~isfinite(relative_total_se)
        relative_total_se = inf;
    end
    crosscheck_good = isfinite(crosscheck_error) ...
        && crosscheck_error <= cfg.crosscheck_good_tolerance;
    before_oscillation = selected.end_pos < oscillation_pos;

    if plateau.exists && selected.valid && relative_total_se <= 0.35 ...
            && crosscheck_good && ~multimode && before_oscillation
        quality = "good";
    elseif selected.admissible && (plateau.exists || selected.r_squared >= cfg.min_r_squared) ...
            && relative_total_se <= 1.0
        quality = "acceptable";
    else
        quality = "poor";
    end
end

function value = robust_std(values)
    values = values(isfinite(values));
    if isempty(values)
        value = NaN;
    else
        value = 1.4826 * median(abs(values - median(values)));
    end
end

function value = combine_uncertainties(fit_se, window_std)
    components = [fit_se, window_std];
    components = components(isfinite(components));
    if isempty(components)
        value = NaN;
    else
        value = sqrt(sum(components .^ 2));
    end
end

function value = local_percentile(values, percentile)
    values = sort(values(isfinite(values)));
    if isempty(values)
        value = NaN;
        return
    end
    if isscalar(values)
        value = values;
        return
    end
    position = 1 + (numel(values) - 1) * percentile / 100;
    lower = floor(position);
    upper = ceil(position);
    fraction = position - lower;
    value = values(lower) * (1 - fraction) + values(upper) * fraction;
end

function summary = make_summary_table(analyses, cfg)
    rows = repmat(summary_row_template(), numel(analyses), 1);
    for i = 1:numel(analyses)
        rows(i) = analysis_to_row(analyses{i}, cfg);
    end
    summary = struct2table(rows, "AsArray", true);
end

function row = analysis_to_row(a, cfg)
    row = summary_row_template();
    fields = fieldnames(row);
    for i = 1:numel(fields)
        name = fields{i};
        if isfield(a, name)
            row.(name) = a.(name);
        end
    end
    if a.status == "ok"
        row.tau_surface = a.tau_energy;
        row.tau_shape_surface = 2 * a.tau_energy;
        row.tau_bulk_scaled = cfg.bulk_time_factor * a.Sd * a.tau_energy;
        row.tau_shape_bulk_scaled = cfg.bulk_time_factor * a.Sd * 2 * a.tau_energy;
        row.decay_rate_surface = 1 / a.tau_energy;
        row.decay_rate_bulk_scaled = 1 / row.tau_bulk_scaled;
    end
end

function row = summary_row_template()
    row = struct( ...
        "simulation_id", "", "Sd", NaN, "status", "failed", "quality_flag", "failed", ...
        "n_total", 0, "dt", NaN, "fit_start_index", NaN, "fit_end_index", NaN, ...
        "oscillation_start_index", NaN, "n_fit_points", 0, ...
        "E_initial", NaN, "E_inf", NaN, "A", NaN, ...
        "tau_energy", NaN, "tau_energy_fit_se", NaN, ...
        "tau_energy_window_std", NaN, "tau_energy_total_se", NaN, ...
        "tau_energy_ci_low", NaN, "tau_energy_ci_high", NaN, ...
        "tau_shape_equivalent", NaN, "tau_1e", NaN, ...
        "tau_crosscheck_error", NaN, "tau_surface", NaN, ...
        "tau_shape_surface", NaN, "tau_bulk_scaled", NaN, ...
        "tau_shape_bulk_scaled", NaN, "decay_rate_surface", NaN, ...
        "decay_rate_bulk_scaled", NaN, "r_squared", NaN, "rmse", NaN, ...
        "normalized_rmse", NaN, "aic", NaN, ...
        "plateau_start_end_index", NaN, "plateau_end_end_index", NaN, ...
        "n_plateau_windows", 0, "multimode_or_nonlinear", false, ...
        "fit_hit_bounds", false, "double_exponential_delta_aic", NaN, ...
        "late_oscillation_amplitude", NaN, "notes", "");
end

function scaling = fit_sd_scaling(summary, cfg, include_acceptable)
    if nargin < 3
        include_acceptable = cfg.include_acceptable_in_scaling;
    end
    scaling = empty_scaling_result();
    included_quality = summary.quality_flag == "good";
    if include_acceptable
        included_quality = included_quality | summary.quality_flag == "acceptable";
    end
    use = summary.status == "ok" & included_quality ...
        & isfinite(summary.Sd) & summary.Sd > 0 ...
        & isfinite(summary.tau_energy) & summary.tau_energy > 0;
    scaling.used_rows = find(use);
    scaling.n_used = nnz(use);
    if scaling.n_used < cfg.minimum_scaling_points
        scaling.notes = "Too few good or acceptable traces for a scaling fit.";
        return
    end

    Sd = summary.Sd(use);
    tau = summary.tau_energy(use);
    se = summary.tau_energy_total_se(use);
    reliable_uncertainty = isfinite(se) & se > 0 & se ./ tau < 1;
    if nnz(reliable_uncertainty) == numel(se)
        weights = 1 ./ se .^ 2;
        scaling.weighting = "inverse variance";
    else
        weights = ones(size(tau));
        scaling.weighting = "equal";
    end
    weights = weights / mean(weights);

    objective = @(z) sum(weights .* (tau - exp(z(1)) - exp(z(2)) ./ Sd) .^ 2);
    initial_a = max(median(tau(Sd >= median(Sd))), eps);
    initial_b = max(median(tau(Sd <= median(Sd)) .* Sd(Sd <= median(Sd))), eps);
    [z, ~, exitflag] = fminsearch(objective, log([initial_a, initial_b]), ...
        optimset("Display", "off", "MaxIter", 2000, "MaxFunEvals", 8000));
    if exitflag <= 0
        scaling.notes = "Crossover optimizer did not converge.";
        return
    end

    scaling.a = exp(z(1));
    scaling.b = exp(z(2));
    prediction = scaling.a + scaling.b ./ Sd;
    residual = tau - prediction;
    scaling.weighted_rmse = sqrt(sum(weights .* residual .^ 2) / sum(weights));
    scaling.r_squared = 1 - sum(residual .^ 2) / sum((tau - mean(tau)) .^ 2);
    scaling.residuals = residual;
    scaling.Sd = Sd;
    scaling.prediction = prediction;

    J = [ones(size(Sd)), 1 ./ Sd];
    dof = numel(Sd) - 2;
    if dof > 0 && rank(J) == 2
        W = diag(weights);
        covariance = (sum(weights .* residual .^ 2) / dof) * pinv(J' * W * J);
        parameter_se = sqrt(max(diag(covariance), 0));
        scaling.a_se = parameter_se(1);
        scaling.b_se = parameter_se(2);
        scaling.a_ci = [max(0, scaling.a - 1.96 * scaling.a_se), ...
            scaling.a + 1.96 * scaling.a_se];
        scaling.b_ci = [max(0, scaling.b - 1.96 * scaling.b_se), ...
            scaling.b + 1.96 * scaling.b_se];
    end

    [Sd_sorted, order] = sort(Sd);
    tau_sorted = tau(order);
    subset_count = max(3, ceil(numel(Sd_sorted) / 3));
    if numel(Sd_sorted) >= 6
        [scaling.low_slope, scaling.low_slope_se] = loglog_slope( ...
            Sd_sorted(1:subset_count), tau_sorted(1:subset_count));
        [scaling.high_slope, scaling.high_slope_se] = loglog_slope( ...
            Sd_sorted(end - subset_count + 1:end), tau_sorted(end - subset_count + 1:end));
    end
    scaling.success = true;
end

function [slope, slope_se] = loglog_slope(x, y)
    X = [ones(numel(x), 1), log(x(:))];
    coefficients = X \ log(y(:));
    slope = coefficients(2);
    residual = log(y(:)) - X * coefficients;
    dof = numel(x) - 2;
    if dof > 0
        covariance = (sum(residual .^ 2) / dof) * pinv(X' * X);
        slope_se = sqrt(max(covariance(2, 2), 0));
    else
        slope_se = NaN;
    end
end

function make_trace_diagnostic(a, output_dir, cfg)
    fig = figure("Visible", cfg.figure_visible, "Color", "w", ...
        "Position", cfg.figure_position);
    layout = tiledlayout(fig, 2, 3, "TileSpacing", "compact", "Padding", "compact");
    selected = a.candidate_fits(a.selected_fit_index);
    fit_last = min([numel(a.t), selected.end_pos + 5, a.oscillation_start_index + 1]);
    fit_positions = selected.start_pos:fit_last;
    fit_curve = selected.E_inf + selected.A * exp( ...
        -(a.t(fit_positions) - a.t(selected.start_pos)) / selected.tau);

    ax = nexttile(layout, 1);
    plot(ax, a.t, a.E, "k-", "LineWidth", 1);
    hold(ax, "on");
    plot(ax, a.t(fit_positions), fit_curve, "r--", "LineWidth", cfg.line_width);
    yline(ax, selected.E_inf, ":", "E_\infty");
    xline(ax, a.t(selected.start_pos), "--", "start");
    xline(ax, a.t(selected.end_pos), "--", "end");
    xline(ax, a.t(a.oscillation_start_index + 1), ":", "oscillation");
    xlabel(ax, "Surface-scaled time"); ylabel(ax, "Energy"); title(ax, "Full trace");

    ax = nexttile(layout, 2);
    early_last = min(numel(a.t), cfg.max_early_index + 1);
    plot(ax, a.t(1:early_last), a.E(1:early_last), "ko-", ...
        "MarkerSize", 3, "LineWidth", 1);
    hold(ax, "on");
    plot(ax, a.t(fit_positions), fit_curve, "r--", "LineWidth", cfg.line_width);
    xlabel(ax, "Surface-scaled time"); ylabel(ax, "Energy"); title(ax, "Early trace");

    ax = nexttile(layout, 3);
    selected_positions = selected.start_pos:selected.end_pos;
    excess = a.E(selected_positions) - selected.E_inf;
    keep = excess > 0;
    semilogy(ax, a.t(selected_positions(keep)), excess(keep), "ko", "MarkerSize", 4);
    hold(ax, "on");
    fitted_excess = selected.A * exp( ...
        -(a.t(selected_positions) - a.t(selected.start_pos)) / selected.tau);
    semilogy(ax, a.t(selected_positions), fitted_excess, "r--", "LineWidth", cfg.line_width);
    xlabel(ax, "Surface-scaled time"); ylabel(ax, "E - E_\infty");
    title(ax, "Semilog excess energy");

    ax = nexttile(layout, 4);
    endpoints = [a.candidate_fits.end_index];
    tau = [a.candidate_fits.tau];
    accepted = [a.candidate_fits.valid];
    scatter(ax, endpoints(~accepted), tau(~accepted), cfg.marker_size, [0.65, 0.65, 0.65], "x");
    hold(ax, "on");
    scatter(ax, endpoints(accepted), tau(accepted), cfg.marker_size, [0.1, 0.35, 0.8], "filled");
    if ~isempty(a.plateau_fit_indices)
        plateau_end = [a.candidate_fits(a.plateau_fit_indices).end_index];
        plateau_tau = [a.candidate_fits(a.plateau_fit_indices).tau];
        scatter(ax, plateau_end, plateau_tau, cfg.marker_size + 18, [0.85, 0.25, 0.1], "o", ...
            "LineWidth", 1.5);
    end
    yline(ax, a.tau_energy, "r--", "median \tau_E");
    xlabel(ax, "Candidate end index"); ylabel(ax, "\tau_E");
    title(ax, "Endpoint stability");

    ax = nexttile(layout, 5);
    if ~isempty(a.local_rate_time)
        plot(ax, a.local_rate_time, a.local_decay_rate, "k.-", "LineWidth", 1);
        hold(ax, "on");
        yline(ax, 1 / a.tau_energy, "r--", "1/\tau_E");
    end
    xlabel(ax, "Surface-scaled time"); ylabel(ax, "-d log(E-E_\infty)/dt");
    title(ax, "Local decay rate");

    title_text = sprintf("%s: %s, tau_E = %.4g", a.simulation_id, a.quality_flag, a.tau_energy);
    title(layout, strrep(title_text, "_", "\_"), "Interpreter", "tex");
    set(findall(fig, "Type", "axes"), "FontSize", cfg.axes_font_size);
    filename = fullfile(output_dir, sanitize_filename(a.simulation_id) + "_diagnostic.png");
    exportgraphics(fig, filename, "Resolution", cfg.output_resolution);
    close(fig);
end

function make_batch_plots(analyses, summary, scaling, output_dir, cfg)
    ok = summary.status == "ok" & summary.Sd > 0 & summary.tau_energy > 0;
    accepted = ok & (summary.quality_flag == "good" | summary.quality_flag == "acceptable");
    poor = ok & summary.quality_flag == "poor";

    if cfg.batch_timescale_only
        fig = figure("Visible", cfg.figure_visible, "Color", "w", ...
            "Position", cfg.figure_position);
        ax_tau = axes(fig);
        hold(ax_tau, "on");
        draw_energy_timescale_plot(ax_tau, summary, scaling, accepted, poor, cfg);
        title(ax_tau, sprintf("Energy relaxation timescale: Da = %.3g, gamy = %.3g, v = %.3g", ...
            cfg.Da, cfg.gamy, cfg.v));
        set(ax_tau, "FontSize", cfg.axes_font_size, "Box", "on");
        grid(ax_tau, "on");
        filename = fullfile(output_dir, "membrane_energy_timescale_batch.png");
        exportgraphics(fig, filename, "Resolution", cfg.output_resolution);
        close(fig);
        fprintf("Saved %s\n", filename);
        return
    end

    colors = turbo(max(numel(analyses), 1));
    fig = figure("Visible", cfg.figure_visible, "Color", "w", ...
        "Position", cfg.figure_position);
    layout = tiledlayout(fig, 2, 3, "TileSpacing", "compact", "Padding", "compact");

    ax_surface = nexttile(layout, 1); hold(ax_surface, "on");
    ax_bulk = nexttile(layout, 2); hold(ax_bulk, "on");
    for i = 1:numel(analyses)
        a = analyses{i};
        if a.status ~= "ok"
            continue
        end
        denominator = a.E(a.fit_start_index + 1) - a.E_inf;
        if denominator <= 0
            continue
        end
        last = min(numel(a.t), cfg.max_early_index + 1);
        idx = (a.fit_start_index + 1):last;
        normalized = (a.E(idx) - a.E_inf) / denominator;
        plot(ax_surface, a.t(idx) - a.t(idx(1)), normalized, ...
            "Color", colors(i, :), "LineWidth", cfg.line_width, ...
            "DisplayName", sprintf("Sd = %.3g", a.Sd));
        bulk_time = cfg.bulk_time_factor * a.Sd * (a.t(idx) - a.t(idx(1)));
        plot(ax_bulk, bulk_time, normalized, "Color", colors(i, :), ...
            "LineWidth", cfg.line_width, "DisplayName", sprintf("Sd = %.3g", a.Sd));
    end
    xlabel(ax_surface, "Surface-scaled time"); ylabel(ax_surface, "Normalized excess energy");
    title(ax_surface, "Surface-time collapse");
    legend(ax_surface, "Location", "best", "FontSize", 8);
    xlabel(ax_bulk, "Bulk time factor \times Sd \times t"); ylabel(ax_bulk, "Normalized excess energy");
    title(ax_bulk, "Bulk-time collapse");
    legend(ax_bulk, "Location", "best", "FontSize", 8);

    ax_tau = nexttile(layout, 3); hold(ax_tau, "on");
    draw_energy_timescale_plot(ax_tau, summary, scaling, accepted, poor, cfg);
    title(ax_tau, "Energy timescale");

    ax_bulk_tau = nexttile(layout, 4); hold(ax_bulk_tau, "on");
    loglog(ax_bulk_tau, summary.Sd(accepted), summary.tau_bulk_scaled(accepted), "o", ...
        "MarkerFaceColor", [0.1, 0.35, 0.8], "MarkerEdgeColor", "k", ...
        "MarkerSize", 6);
    if any(poor)
        loglog(ax_bulk_tau, summary.Sd(poor), summary.tau_bulk_scaled(poor), "x", ...
            "Color", [0.7, 0.2, 0.2], "MarkerSize", 7);
    end
    set(ax_bulk_tau, "XScale", "log", "YScale", "log");
    xlabel(ax_bulk_tau, "Sd"); ylabel(ax_bulk_tau, "Bulk time factor \times Sd \times \tau_E");
    title(ax_bulk_tau, "Bulk-scaled timescale");

    ax_rate = nexttile(layout, 5); hold(ax_rate, "on");
    loglog(ax_rate, summary.Sd(accepted), summary.decay_rate_surface(accepted), "o", ...
        "MarkerFaceColor", [0.1, 0.35, 0.8], "MarkerEdgeColor", "k", ...
        "MarkerSize", 6);
    if scaling.success
        xfit = logspace(log10(min(summary.Sd(accepted))), log10(max(summary.Sd(accepted))), 300);
        loglog(ax_rate, xfit, 1 ./ (scaling.a + scaling.b ./ xfit), "k-", ...
            "LineWidth", cfg.line_width);
    end
    set(ax_rate, "XScale", "log", "YScale", "log");
    xlabel(ax_rate, "Sd"); ylabel(ax_rate, "1/\tau_E"); title(ax_rate, "Decay rate");

    ax_residual = nexttile(layout, 6); hold(ax_residual, "on");
    if scaling.success
        semilogx(ax_residual, scaling.Sd, scaling.residuals, "ko", "MarkerFaceColor", [0.1, 0.35, 0.8]);
        yline(ax_residual, 0, "k:");
    end
    set(ax_residual, "XScale", "log");
    xlabel(ax_residual, "Sd"); ylabel(ax_residual, "\tau_E - (a+b/Sd)");
    title(ax_residual, "Crossover residuals");

    title(layout, sprintf("Early membrane-energy relaxation: Da = %.3g, gamy = %.3g, v = %.3g", ...
        cfg.Da, cfg.gamy, cfg.v));
    axes_list = findall(fig, "Type", "axes");
    set(axes_list, "FontSize", cfg.axes_font_size, "Box", "on");
    grid(axes_list, "on");
    filename = fullfile(output_dir, "membrane_energy_timescale_batch.png");
    exportgraphics(fig, filename, "Resolution", cfg.output_resolution);
    close(fig);
    fprintf("Saved %s\n", filename);
end

function draw_energy_timescale_plot(ax, summary, scaling, accepted, poor, cfg)
    plot_timescale_points(ax, summary, accepted, poor, "tau_energy", cfg);
    if scaling.success
        xfit = logspace(log10(min(summary.Sd(accepted))), ...
            log10(max(summary.Sd(accepted))), 300);
        loglog(ax, xfit, scaling.a + scaling.b ./ xfit, "k-", ...
            "LineWidth", cfg.line_width, "DisplayName", "a + b/Sd");
        low_anchor = scaling.a + scaling.b / xfit(1);
        reference = low_anchor * (xfit / xfit(1)) .^ (-1);
        loglog(ax, xfit, reference, "k:", "LineWidth", 1.2, ...
            "DisplayName", "slope -1");
    end
    xlabel(ax, "Sd");
    ylabel(ax, "\tau_E");
    legend(ax, "Location", "best");
end

function plot_timescale_points(ax, summary, accepted, poor, field, ~)
    x = summary.Sd(accepted);
    y = summary.(field)(accepted);
    error_value = summary.tau_energy_total_se(accepted);
    finite_error = isfinite(error_value) & error_value >= 0;
    if any(finite_error)
        errorbar(ax, x(finite_error), y(finite_error), error_value(finite_error), ...
            "o", "LineStyle", "none", "Color", [0.1, 0.35, 0.8], ...
            "MarkerFaceColor", [0.1, 0.35, 0.8], "MarkerEdgeColor", "k", ...
            "LineWidth", 1, "DisplayName", "good/acceptable");
    end
    if any(~finite_error)
        loglog(ax, x(~finite_error), y(~finite_error), "o", ...
            "MarkerFaceColor", [0.1, 0.35, 0.8], "MarkerEdgeColor", "k", ...
            "DisplayName", "good/acceptable");
    end
    if any(poor)
        loglog(ax, summary.Sd(poor), summary.(field)(poor), "x", ...
            "Color", [0.75, 0.2, 0.2], "MarkerSize", 7, ...
            "LineWidth", 1.2, "DisplayName", "poor");
    end
    set(ax, "XScale", "log", "YScale", "log");
    if isempty(x)
        text(ax, 0.5, 0.5, "No accepted traces", "Units", "normalized", ...
            "HorizontalAlignment", "center");
    end
end

function sensitivity = run_sensitivity_analysis(analyses, primary_summary, cfg)
    successful = analyses(cellfun(@(a) a.status == "ok", analyses));
    rows = repmat(sensitivity_row_template(), 0, 1);
    row_index = 0;
    for end_max = cfg.sensitivity_end_max
        for smooth_window = cfg.sensitivity_smooth_windows
            for loss = cfg.sensitivity_losses
                variant_cfg = cfg;
                variant_cfg.candidate_end_max = end_max;
                variant_cfg.smooth_window = smooth_window;
                variant_cfg.robust_loss = loss;
                variant_analyses = cell(numel(successful), 1);
                for i = 1:numel(successful)
                    base = successful{i};
                    variant_analyses{i} = analyze_trace(base.simulation_id, base.Sd, ...
                        base.t, base.E, base.frame_ids, base.dt, variant_cfg);
                end
                variant_summary = make_summary_table(variant_analyses, variant_cfg);
                for include_acceptable = [true, false]
                    fit = fit_sd_scaling(variant_summary, variant_cfg, include_acceptable);
                    row_index = row_index + 1;
                    rows(row_index) = sensitivity_row_template();
                    rows(row_index).candidate_end_max = end_max;
                    rows(row_index).smooth_window = smooth_window;
                    rows(row_index).robust_loss = loss;
                    rows(row_index).include_acceptable = include_acceptable;
                    rows(row_index).n_used = fit.n_used;
                    rows(row_index).a = fit.a;
                    rows(row_index).b = fit.b;
                    rows(row_index).low_slope = fit.low_slope;
                    rows(row_index).high_slope = fit.high_slope;
                    rows(row_index).median_relative_tau_change = compare_tau_summaries( ...
                        primary_summary, variant_summary);
                end
            end
        end
    end
    sensitivity = struct2table(rows, "AsArray", true);
end

function row = sensitivity_row_template()
    row = struct("candidate_end_max", NaN, "smooth_window", NaN, ...
        "robust_loss", "", "include_acceptable", false, "n_used", 0, ...
        "a", NaN, "b", NaN, "low_slope", NaN, "high_slope", NaN, ...
        "median_relative_tau_change", NaN);
end

function value = compare_tau_summaries(primary, variant)
    value = NaN;
    relative_changes = [];
    for i = 1:height(variant)
        match = find(primary.simulation_id == variant.simulation_id(i), 1);
        if isempty(match) || ~isfinite(primary.tau_energy(match)) ...
                || ~isfinite(variant.tau_energy(i)) || primary.tau_energy(match) <= 0
            continue
        end
        relative_changes(end + 1) = abs(variant.tau_energy(i) ...
            - primary.tau_energy(match)) / primary.tau_energy(match); %#ok<AGROW>
    end
    if ~isempty(relative_changes)
        value = median(relative_changes);
    end
end

function write_report(summary, scaling, sensitivity, output_dir, cfg)
    report_file = fullfile(output_dir, "membrane_energy_timescale_report.md");
    file_id = fopen(report_file, "w");
    if file_id < 0
        warning("Could not open %s for writing.", report_file);
        return
    end
    cleanup = onCleanup(@() fclose(file_id));

    fprintf(file_id, "# Membrane Energy Timescale Analysis\n\n");
    fprintf(file_id, "All reported times are nondimensional simulation time. ");
    fprintf(file_id, "The bulk-rescaled time uses `bulk_time_factor = %.8g`.\n\n", ...
        cfg.bulk_time_factor);
    flags = ["good", "acceptable", "poor", "failed"];
    for flag = flags
        if flag == "failed"
            count = nnz(summary.status == "failed");
        else
            count = nnz(summary.quality_flag == flag);
        end
        fprintf(file_id, "- %s: %d\n", flag, count);
    end

    fprintf(file_id, "\n## Scaling fit\n\n");
    if scaling.success
        fprintf(file_id, "The crossover interpolation `tau_E = a + b/Sd` used %d traces ", scaling.n_used);
        fprintf(file_id, "with %s weighting.\n\n", scaling.weighting);
        fprintf(file_id, "- `a = %.8g` (95%% CI %.8g to %.8g)\n", ...
            scaling.a, scaling.a_ci(1), scaling.a_ci(2));
        fprintf(file_id, "- `b = %.8g` (95%% CI %.8g to %.8g)\n", ...
            scaling.b, scaling.b_ci(1), scaling.b_ci(2));
        fprintf(file_id, "- linear-space R-squared: %.6g\n", scaling.r_squared);
        fprintf(file_id, "- weighted RMSE: %.6g\n", scaling.weighted_rmse);
        if isfinite(scaling.low_slope)
            fprintf(file_id, "- low-Sd log-log slope: %.6g +/- %.6g\n", ...
                scaling.low_slope, scaling.low_slope_se);
            fprintf(file_id, "- high-Sd log-log slope: %.6g +/- %.6g\n", ...
                scaling.high_slope, scaling.high_slope_se);
        else
            fprintf(file_id, "- Too few traces were available for separate low/high-Sd slopes.\n");
        end
        fprintf(file_id, "\nSurface-time control at high Sd is supported when the high-Sd slope ");
        fprintf(file_id, "is consistent with zero and the surface-time curves collapse.\n\n");
        fprintf(file_id, "Bulk-time control at low Sd is supported when the low-Sd slope is ");
        fprintf(file_id, "consistent with -1 and `bulk_time_factor * Sd * tau_E` plateaus. ");
        fprintf(file_id, "The collapse panels must be inspected before making either claim.\n");
    else
        fprintf(file_id, "%s\n", scaling.notes);
    end

    fprintf(file_id, "\n## Sensitivity\n\n");
    if isempty(sensitivity)
        fprintf(file_id, "Sensitivity analysis was disabled. Set ");
        fprintf(file_id, "`cfg.run_sensitivity_analysis = true` to scan endpoint maxima, ");
        fprintf(file_id, "smoothing windows, robust/ordinary fits, and quality inclusion.\n");
    else
        finite_change = sensitivity.median_relative_tau_change( ...
            isfinite(sensitivity.median_relative_tau_change));
        if isempty(finite_change)
            fprintf(file_id, "Sensitivity variants did not produce comparable timescales.\n");
        else
            fprintf(file_id, "Across configured variants, the median per-trace relative ");
            fprintf(file_id, "timescale change ranged from %.4g to %.4g.\n", ...
                min(finite_change), max(finite_change));
        end
    end

    inspect = summary.simulation_id(summary.quality_flag == "poor" | summary.status == "failed");
    fprintf(file_id, "\n## Manual inspection\n\n");
    if isempty(inspect)
        fprintf(file_id, "No traces were flagged poor or failed.\n");
    else
        for id = inspect.'
            fprintf(file_id, "- `%s`\n", id);
        end
    end
    fprintf("Saved %s\n", report_file);
end

function run_synthetic_tests(cfg)
    fprintf("Running membrane_energy_timescale synthetic tests...\n");
    rng(7);
    cfg.save_trace_diagnostics = false;
    cfg.save_batch_plots = false;
    cfg.run_sensitivity_analysis = false;
    cfg.require_expected_n_samples = false;
    frame_ids = (0:299).';

    t = linspace(0, 24, 300).';
    truth = 1.7;
    E = 18 + 6 * exp(-t / truth) + 0.002 * randn(size(t));
    a = analyze_trace("synthetic_clean", 1, t, E, frame_ids, median(diff(t)), cfg);
    assert(a.status == "ok" && abs(a.tau_energy - truth) / truth < 0.10, ...
        "Clean exponential test failed: recovered tau %.6g, expected %.6g.", a.tau_energy, truth);
    fprintf("  clean exponential: tau %.5g (truth %.5g)\n", a.tau_energy, truth);

    E_osc = 18 + 6 * exp(-t / truth) + 0.002 * randn(size(t));
    onset = 31;
    E_osc(onset:end) = E_osc(onset:end) + 0.35 * sin((0:numel(E_osc) - onset).' * 0.9);
    a_osc = analyze_trace("synthetic_oscillatory", 1, t, E_osc, frame_ids, median(diff(t)), cfg);
    full_fit = fit_single_exponential(t, E_osc, 1, numel(t), cfg);
    early_error = abs(a_osc.tau_energy - truth) / truth;
    full_error = abs(full_fit.tau - truth) / truth;
    assert(a_osc.status == "ok" && early_error < 0.10, ...
        "Late-oscillation test failed to recover the early timescale.");
    assert(a_osc.fit_end_index < onset || full_error > early_error, ...
        "Late-oscillation test did not exclude or improve upon the contaminated fit.");
    fprintf("  late oscillations: early error %.3g, full-trace error %.3g\n", early_error, full_error);

    E_double = 12 + 4 * exp(-t / 0.7) + 3 * exp(-t / 3.5) + 0.003 * randn(size(t));
    a_double = analyze_trace("synthetic_double", 1, t, E_double, frame_ids, median(diff(t)), cfg);
    assert(a_double.status == "ok" && isfinite(a_double.tau_energy), ...
        "Double-exponential test did not return an effective timescale.");
    fprintf("  double exponential diagnostic: drift %.5g, quality %s\n", ...
        a_double.local_rate_drift, a_double.quality_flag);
    assert(a_double.multimode_or_nonlinear || a_double.quality_flag ~= "good", ...
        "Double-exponential test was not flagged as limited-confidence or multimode.");
    fprintf("  double exponential: tau_eff %.5g, quality %s\n", ...
        a_double.tau_energy, a_double.quality_flag);

    Sd_values = logspace(-2, 2, 9).';
    truth_a = 0.8;
    truth_b = 0.25;
    scaling_analyses = cell(numel(Sd_values), 1);
    for i = 1:numel(Sd_values)
        tau = truth_a + truth_b / Sd_values(i);
        ti = linspace(0, 25 * tau, 300).';
        Ei = 9 + 2.5 * exp(-ti / tau) + 2e-4 * randn(size(ti));
        scaling_analyses{i} = analyze_trace("synthetic_scaling_" + i, Sd_values(i), ...
            ti, Ei, frame_ids, median(diff(ti)), cfg);
    end
    scaling_summary = make_summary_table(scaling_analyses, cfg);
    fit = fit_sd_scaling(scaling_summary, cfg, true);
    assert(fit.success && abs(fit.a - truth_a) / truth_a < 0.20 ...
        && abs(fit.b - truth_b) / truth_b < 0.20, ...
        "Sd crossover test failed: a %.6g, b %.6g.", fit.a, fit.b);
    fprintf("  Sd crossover: a %.5g (truth %.5g), b %.5g (truth %.5g)\n", ...
        fit.a, truth_a, fit.b, truth_b);

    increments = 0.04 + 0.08 * rand(299, 1);
    t_irregular = [0; cumsum(increments)];
    irregular_truth = 1.3;
    E_irregular = 5 + 2 * exp(-t_irregular / irregular_truth) + 5e-4 * randn(size(t_irregular));
    a_irregular = analyze_trace("synthetic_irregular", 1, t_irregular, E_irregular, ...
        frame_ids, NaN, cfg);
    assert(a_irregular.status == "ok" ...
        && abs(a_irregular.tau_energy - irregular_truth) / irregular_truth < 0.10, ...
        "Irregular-time-grid test failed.");
    fprintf("  irregular time: tau %.5g (truth %.5g)\n", ...
        a_irregular.tau_energy, irregular_truth);
    fprintf("All synthetic tests passed.\n");
end

function analysis = failed_analysis(simulation_id, Sd, notes)
    analysis = struct( ...
        "simulation_id", string(simulation_id), "Sd", Sd, "status", "failed", ...
        "quality_flag", "failed", "notes", string(notes), "n_total", 0, ...
        "dt", NaN, "fit_start_index", NaN, "fit_end_index", NaN, ...
        "oscillation_start_index", NaN, "n_fit_points", 0, ...
        "E_initial", NaN, "E_inf", NaN, "A", NaN, "tau_energy", NaN, ...
        "tau_energy_fit_se", NaN, "tau_energy_window_std", NaN, ...
        "tau_energy_total_se", NaN, "tau_energy_ci_low", NaN, ...
        "tau_energy_ci_high", NaN, "tau_shape_equivalent", NaN, ...
        "tau_1e", NaN, "tau_crosscheck_error", NaN, "r_squared", NaN, ...
        "rmse", NaN, "normalized_rmse", NaN, "aic", NaN, ...
        "plateau_start_end_index", NaN, "plateau_end_end_index", NaN, ...
        "n_plateau_windows", 0, "multimode_or_nonlinear", false, ...
        "fit_hit_bounds", false, "late_oscillation_amplitude", NaN, ...
        "t", [], "E", [], "frame_ids", [], "E_smooth", [], ...
        "candidate_fits", repmat(empty_fit(), 0, 1), ...
        "plateau_fit_indices", [], "selected_fit_index", NaN, ...
        "local_rate_time", [], "local_decay_rate", [], "local_rate_drift", NaN, ...
        "double_exponential_delta_aic", NaN);
end

function fit = empty_fit()
    fit = struct("start_pos", NaN, "end_pos", NaN, "start_index", NaN, ...
        "end_index", NaN, "n_points", 0, "E_inf", NaN, "A", NaN, ...
        "tau", NaN, "E_inf_se", NaN, "A_se", NaN, "tau_se", NaN, ...
        "rss", NaN, "rmse", NaN, "normalized_rmse", NaN, ...
        "r_squared", NaN, "aic", NaN, "exitflag", -1, "converged", false, ...
        "admissible", false, "valid", false, "hit_bounds", false, ...
        "systematic_residual", false, "extends_past_oscillation", false, ...
        "selection_score", -inf, "rejection_reason", "");
end

function plateau = empty_plateau()
    plateau = struct("exists", false, "fit_indices", [], ...
        "start_end_index", NaN, "end_end_index", NaN);
end

function scaling = empty_scaling_result()
    scaling = struct("success", false, "n_used", 0, "used_rows", [], ...
        "a", NaN, "b", NaN, "a_se", NaN, "b_se", NaN, ...
        "a_ci", [NaN, NaN], "b_ci", [NaN, NaN], ...
        "weighted_rmse", NaN, "r_squared", NaN, "weighting", "none", ...
        "low_slope", NaN, "low_slope_se", NaN, ...
        "high_slope", NaN, "high_slope_se", NaN, ...
        "Sd", [], "prediction", [], "residuals", [], "notes", "");
end

function notes = join_notes(varargin)
    parts = strings(0, 1);
    for i = 1:nargin
        value = strtrim(string(varargin{i}));
        if strlength(value) > 0
            parts(end + 1) = value; %#ok<AGROW>
        end
    end
    notes = strjoin(unique(parts, "stable"), " ");
end

function folder = resolve_run_folder(data_dir, simulation_id, Sd, cfg)
    folder = "";
    exact = fullfile(data_dir, simulation_id);
    if isfolder(exact)
        folder = string(exact);
        return
    end
    if cfg.allow_legacy_folder_without_v
        legacy = make_legacy_run_tag(Sd, cfg.Da, cfg.gamy);
        legacy_folder = fullfile(data_dir, legacy);
        if isfolder(legacy_folder)
            folder = string(legacy_folder);
        end
    end
end

function run_tag = make_run_tag(Sd, Da, gamy, v)
    run_tag = sprintf("Sd_%.2e_Da_%.2e_gamy_%+.2e_v_%.2e", Sd, Da, gamy, v);
    run_tag = strrep(run_tag, "+", "p");
    run_tag = strrep(run_tag, "-", "m");
end

function run_tag = make_legacy_run_tag(Sd, Da, gamy)
    run_tag = sprintf("Sd_%.2e_Da_%.2e_gamy_%+.2e", Sd, Da, gamy);
    run_tag = strrep(run_tag, "+", "p");
    run_tag = strrep(run_tag, "-", "m");
end

function name = sanitize_filename(name)
    name = regexprep(string(name), "[^A-Za-z0-9_.-]", "_");
end

function remesh_dir = find_remesh_dir(script_dir)
    candidates = string({script_dir, pwd, fullfile(pwd, "remesh")});
    geometry_path = which("Geometry");
    if strlength(geometry_path) > 0
        candidates(end + 1) = string(fileparts(geometry_path));
    end
    for candidate = candidates
        if isfile(fullfile(candidate, "Geometry.m"))
            remesh_dir = char(candidate);
            return
        end
    end
    error("Could not locate remesh directory.");
end

function data_dir = find_data_dir(script_dir, remesh_dir)
    candidates = string({ ...
        fullfile(remesh_dir, "data", "fs_batch_data"), ...
        fullfile(script_dir, "data", "fs_batch_data"), ...
        fullfile(pwd, "data", "fs_batch_data"), ...
        fullfile(pwd, "remesh", "data", "fs_batch_data")});
    for candidate = candidates
        if isfolder(candidate)
            data_dir = char(candidate);
            return
        end
    end
    error("Could not locate data/fs_batch_data. Checked:%s", ...
        sprintf("\n  %s", candidates));
end
