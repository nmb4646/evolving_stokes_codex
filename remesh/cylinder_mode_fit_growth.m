function [growth_rates, fit_details] = cylinder_mode_fit_growth( ...
    time, frame_index, amplitude, phase, noise, qR, group_m, group_n, ...
    frame_valid, cfg)
%CYLINDER_MODE_FIT_GROWTH Fit early exponential rates for grouped modes.

time = time(:);
frame_index = frame_index(:);
group_m = group_m(:);
group_n = group_n(:);
frame_valid = logical(frame_valid(:));

n_groups = numel(group_m);
fit_details = repmat(empty_fit_detail(), n_groups, 1);

rows = cell(n_groups, 20);
for g = 1:n_groups
    A = amplitude(:, g);
    N = noise(:, min(g, size(noise, 2)));
    q = qR(:, g);

    eligible = frame_valid & isfinite(time) & isfinite(A) & A > 0 ...
        & isfinite(N) & A >= cfg.growth.minimum_signal_to_noise .* max(N, cfg.growth.absolute_noise_floor) ...
        & A <= cfg.growth.maximum_dimensionless_amplitude ...
        & isfinite(q);

    if cfg.growth.mode == "fixed"
        eligible = eligible & time >= cfg.growth.fixed_start_time ...
            & time <= cfg.growth.fixed_end_time;
        candidate = fit_selected_indices(find(eligible), time, A, q, cfg);
        candidates = candidate;
    else
        candidates = enumerate_early_candidates(eligible, time, A, q, cfg);
    end

    [selected, status, flags] = select_candidate(candidates, cfg);
    if selected.valid
        idx = selected.indices;
        phase_rate = fit_phase_rate(time(idx), phase(idx, g));
        fit_details(g) = selected;
        fit_details(g).phase_rate = phase_rate;
        rows(g, :) = {group_m(g), group_n(g), time(idx(1)), time(idx(end)), ...
            frame_index(idx(1)), frame_index(idx(end)), numel(idx), ...
            selected.slope, selected.slope_se, phase_rate, selected.r_squared, ...
            A(idx(1)), A(idx(end)), selected.amplitude_ratio, ...
            mean(q(idx), "omitnan"), min(q(idx)), max(q(idx)), ...
            selected.normalized_rmse, status, join_flags(flags)};
    else
        rows(g, :) = {group_m(g), group_n(g), NaN, NaN, NaN, NaN, 0, ...
            NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, ...
            status, join_flags(flags)};
    end
end

growth_rates = cell2table(rows, "VariableNames", [ ...
    "m", "n", "fit_start_time", "fit_end_time", ...
    "fit_start_frame", "fit_end_frame", "n_fit_points", ...
    "growth_rate_sigma", "growth_rate_standard_error", "phase_rate_omega", ...
    "r_squared", "amplitude_start", "amplitude_end", "amplitude_ratio", ...
    "mean_qR", "min_qR", "max_qR", "normalized_rmse", ...
    "fit_status", "quality_flags"]);
end

function candidates = enumerate_early_candidates(eligible, time, A, q, cfg)
    last = min(numel(time), cfg.growth.maximum_early_frames);
    eligible((last + 1):end) = false;
    runs = contiguous_runs(eligible);
    candidates = repmat(empty_fit_detail(), 0, 1);

    for r = 1:size(runs, 1)
        run_idx = runs(r, 1):runs(r, 2);
        if numel(run_idx) < cfg.growth.minimum_points
            continue
        end
        for first_local = 1:(numel(run_idx) - cfg.growth.minimum_points + 1)
            maximum_last = min(numel(run_idx), ...
                first_local + cfg.growth.maximum_points - 1);
            for last_local = (first_local + cfg.growth.minimum_points - 1):maximum_last
                idx = run_idx(first_local:last_local);
                candidates(end + 1, 1) = fit_selected_indices(idx, time, A, q, cfg); %#ok<AGROW>
            end
        end
    end
end

function candidate = fit_selected_indices(idx, time, A, q, cfg)
    candidate = empty_fit_detail();
    if numel(idx) < cfg.growth.minimum_points
        return
    end

    t = time(idx);
    y = log(A(idx));
    if any(~isfinite(y)) || max(t) <= min(t)
        return
    end

    [slope, intercept, slope_se, r_squared, nrmse, yfit] = linear_fit(t, y);
    amplitude_ratio = max(A(idx)) / max(min(A(idx)), realmin);
    q_mean = mean(q(idx), "omitnan");
    q_change = (max(q(idx)) - min(q(idx))) / max(abs(q_mean), eps);
    endpoint_change = endpoint_slope_change(t, y, slope);

    candidate.valid = true;
    candidate.indices = idx(:);
    candidate.slope = slope;
    candidate.intercept = intercept;
    candidate.slope_se = slope_se;
    candidate.r_squared = r_squared;
    candidate.normalized_rmse = nrmse;
    candidate.amplitude_ratio = amplitude_ratio;
    candidate.relative_qR_change = q_change;
    candidate.endpoint_slope_change = endpoint_change;
    candidate.yfit = yfit;
    candidate.passes = r_squared >= cfg.growth.minimum_r_squared ...
        && amplitude_ratio >= cfg.growth.minimum_amplitude_ratio ...
        && q_change <= cfg.growth.maximum_relative_qR_change ...
        && endpoint_change <= cfg.growth.maximum_endpoint_slope_change;
end

function [selected, status, flags] = select_candidate(candidates, cfg)
    selected = empty_fit_detail();
    flags = strings(0, 1);

    if isempty(candidates) || ~any([candidates.valid])
        status = "no_fit";
        flags(end + 1) = "INSUFFICIENT_VALID_POINTS";
        return
    end

    valid_candidates = candidates([candidates.valid]);
    passing = valid_candidates([valid_candidates.passes]);
    if ~isempty(passing)
        starts = arrayfun(@(x) x.indices(1), passing);
        earliest = min(starts);
        same_start = passing(starts == earliest);
        lengths = arrayfun(@(x) numel(x.indices), same_start);
        r2 = [same_start.r_squared];
        [~, order] = sortrows([-lengths(:), -r2(:)], [1, 2]);
        selected = same_start(order(1));
        status = "good";
        return
    end

    if ~cfg.growth.accept_poor_fits
        status = "no_fit";
        flags(end + 1) = "POOR_EXPONENTIAL_FIT";
        return
    end

    score = arrayfun(@candidate_score, valid_candidates);
    [~, best] = max(score);
    selected = valid_candidates(best);
    status = "poor";

    if selected.r_squared < cfg.growth.minimum_r_squared
        flags(end + 1) = "POOR_EXPONENTIAL_FIT";
    end
    if selected.amplitude_ratio < cfg.growth.minimum_amplitude_ratio
        flags(end + 1) = "INSUFFICIENT_AMPLITUDE_CHANGE";
    end
    if selected.relative_qR_change > cfg.growth.maximum_relative_qR_change
        flags(end + 1) = "QR_CHANGES_TOO_MUCH";
    end
    if selected.endpoint_slope_change > cfg.growth.maximum_endpoint_slope_change
        flags(end + 1) = "ENDPOINT_SENSITIVE";
    end
end

function score = candidate_score(candidate)
    score = candidate.r_squared ...
        + 0.05 * min(log(max(candidate.amplitude_ratio, 1)), 2) ...
        - 0.15 * min(candidate.relative_qR_change, 2) ...
        - 0.10 * min(candidate.endpoint_slope_change, 2) ...
        - 1e-4 * candidate.indices(1) ...
        + 1e-5 * numel(candidate.indices);
end

function [slope, intercept, slope_se, r_squared, nrmse, yfit] = linear_fit(t, y)
    tc = t - mean(t);
    X = [ones(size(tc)), tc];
    beta = X \ y;
    yfit = X * beta;
    residual = y - yfit;
    slope = beta(2);
    intercept = beta(1) - slope * mean(t);

    sse = sum(residual .^ 2);
    sst = sum((y - mean(y)) .^ 2);
    r_squared = 1 - sse / max(sst, eps(max(abs(y))) ^ 2);
    nrmse = sqrt(mean(residual .^ 2)) / max(max(y) - min(y), eps);

    if numel(y) > 2 && sum(tc .^ 2) > 0
        slope_se = sqrt((sse / (numel(y) - 2)) / sum(tc .^ 2));
    else
        slope_se = NaN;
    end
end

function change = endpoint_slope_change(t, y, slope)
    if numel(t) < 2 + 3
        change = inf;
        return
    end
    slope_left = simple_slope(t(2:end), y(2:end));
    slope_right = simple_slope(t(1:end-1), y(1:end-1));
    scale = max(abs(slope), 1 / max(t(end) - t(1), eps));
    change = max(abs([slope_left - slope, slope_right - slope])) / scale;
end

function slope = simple_slope(t, y)
    tc = t - mean(t);
    slope = sum(tc .* (y - mean(y))) / max(sum(tc .^ 2), eps);
end

function omega = fit_phase_rate(t, phase)
    valid = isfinite(t) & isfinite(phase);
    if nnz(valid) < 3
        omega = NaN;
        return
    end
    unwrapped = unwrap(phase(valid));
    omega = simple_slope(t(valid), unwrapped);
end

function runs = contiguous_runs(mask)
    mask = mask(:);
    edges = diff([false; mask; false]);
    runs = [find(edges == 1), find(edges == -1) - 1];
end

function detail = empty_fit_detail()
    detail = struct( ...
        "valid", false, ...
        "passes", false, ...
        "indices", zeros(0, 1), ...
        "slope", NaN, ...
        "intercept", NaN, ...
        "slope_se", NaN, ...
        "phase_rate", NaN, ...
        "r_squared", NaN, ...
        "normalized_rmse", NaN, ...
        "amplitude_ratio", NaN, ...
        "relative_qR_change", NaN, ...
        "endpoint_slope_change", NaN, ...
        "yfit", zeros(0, 1));
end

function value = join_flags(flags)
    flags = unique(flags(strlength(flags) > 0), "stable");
    if isempty(flags)
        value = "";
    else
        value = join(flags, ";");
    end
end
