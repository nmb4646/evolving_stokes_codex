function responses = cylinder_mode_pair_responses(series_results, cfg)
%CYLINDER_MODE_PAIR_RESPONSES Build plus/minus or reference-subtracted signals.

n_series = numel(series_results);
responses = {};
if isempty(cfg.pairing.roles)
    warning("CylinderMode:MissingPairRoles", ...
        "Pairing is enabled, but cfg.pairing.roles is empty.");
    return
end
if numel(cfg.pairing.roles) ~= n_series
    error("CylinderMode:PairingConfig", ...
        "cfg.pairing.roles must have one entry per series.");
end

roles = string(cfg.pairing.roles(:));
if isempty(cfg.pairing.group_ids)
    group_ids = repmat("default", n_series, 1);
else
    group_ids = string(cfg.pairing.group_ids(:));
    if numel(group_ids) ~= n_series
        error("CylinderMode:PairingConfig", ...
            "cfg.pairing.group_ids must have one entry per series.");
    end
end

successful = cellfun(@(r) isfield(r, "status") && r.status == "ok", series_results);
for group_id = unique(group_ids(successful), "stable").'
    members = find(successful & group_ids == group_id);
    member_roles = roles(members);
    plus = members(member_roles == "plus");
    minus = members(member_roles == "minus");
    reference = members(member_roles == "reference");

    if ~isempty(plus) && ~isempty(minus)
        epsilon = pairing_epsilon(cfg, plus(1), n_series);
        response_id = "cylinder_response_" + sanitize_name(group_id) + "_plus_minus";
        response = combine_series( ...
            series_results{plus(1)}, series_results{minus(1)}, ...
            2 * epsilon, response_id, "plus_minus", cfg);
        responses{end + 1, 1} = response; %#ok<AGROW>
        continue
    end

    perturbed = members(member_roles == "ordinary" | member_roles == "plus");
    if ~isempty(reference)
        for index = perturbed(:).'
            response_id = "cylinder_response_" + sanitize_name(group_id) ...
                + "_" + sanitize_name(series_results{index}.series_id) ...
                + "_minus_reference";
            response = combine_series( ...
                series_results{index}, series_results{reference(1)}, ...
                1, response_id, "reference_subtraction", cfg);
            responses{end + 1, 1} = response; %#ok<AGROW>
        end
    else
        warning("CylinderMode:IncompletePair", ...
            "Pairing group '%s' has no complete plus/minus or reference pair.", group_id);
    end
end
end

function response = combine_series(a, b, denominator, response_id, method, cfg)
    verify_compatible_modes(a, b);
    overlap_start = max(min(a.time), min(b.time));
    overlap_end = min(max(a.time), max(b.time));
    selected_a = a.time >= overlap_start & a.time <= overlap_end;
    time = a.time(selected_a);
    if numel(time) < cfg.growth.minimum_points
        error("CylinderMode:PairingTime", ...
            "Related series do not have enough overlapping physical times.");
    end

    [joint_available, b_state_indices] = matching_projection_states(a, b, selected_a, time, cfg);
    if joint_available
        [coefficient, group_qR, group_noise] = joint_pair_projection( ...
            a, b, find(selected_a), b_state_indices, denominator, cfg);
        method = string(method) + "_shared_projection";
    else
        coefficient_b = interpolate_matrix(b.time, b.coefficient, time);
        coefficient = (a.coefficient(selected_a, :) - coefficient_b) / denominator;
        qR_b = interpolate_matrix(b.time, b.group_qR, time);
        group_qR = 0.5 * (a.group_qR(selected_a, :) + qR_b);
        noise_b = interpolate_matrix(b.time, b.group_noise, time);
        group_noise = hypot(a.group_noise(selected_a, :), noise_b) / abs(denominator);
        method = string(method) + "_independent_projection";
    end
    resolved_b = interpolate_matrix(b.time, double(b.group_resolved), time) >= 0.5;
    group_resolved = a.group_resolved(selected_a, :) & resolved_b;

    n_groups = numel(a.group_m);
    group_amplitude = nan(numel(time), n_groups);
    group_phase = nan(numel(time), n_groups);
    for g = 1:n_groups
        primary = find(a.mode_m == a.group_m(g) & a.mode_n == a.group_n(g), 1);
        conjugate = find(a.mode_m == -a.group_m(g) & a.mode_n == -a.group_n(g), 1);
        c1 = coefficient(:, primary);
        if primary == conjugate
            group_amplitude(:, g) = abs(c1);
        else
            group_amplitude(:, g) = hypot(abs(c1), abs(coefficient(:, conjugate)));
        end
        group_phase(:, g) = angle(c1);
    end

    frame_metrics = a.frame_metrics(selected_a, :);
    frame_metrics.series_id(:) = response_id;
    frame_metrics.time = time;
    valid_b = interpolate_matrix(b.time, double(b.frame_metrics.frame_valid), time) >= 0.5;
    frame_valid = frame_metrics.frame_valid & valid_b;

    noise_for_fit = group_noise;
    noise_for_fit(~group_resolved) = inf;
    [growth_rates, fit_details] = cylinder_mode_fit_growth( ...
        time, frame_metrics.frame_index, group_amplitude, group_phase, ...
        noise_for_fit, group_qR, a.group_m, a.group_n, frame_valid, cfg);

    growth_rates.series_id = repmat(response_id, height(growth_rates), 1);
    growth_rates.mean_radius = nan(height(growth_rates), 1);
    growth_rates.mean_core_length = nan(height(growth_rates), 1);
    growth_rates.resolved_fraction = mean(group_resolved, 1).';
    for g = 1:n_groups
        idx = fit_details(g).indices;
        if ~isempty(idx)
            growth_rates.mean_radius(g) = mean(frame_metrics.mean_core_radius(idx), "omitnan");
            growth_rates.mean_core_length(g) = mean(frame_metrics.core_length(idx), "omitnan");
        end
        if a.group_m(g) == 0 && a.group_n(g) == 0
            growth_rates.fit_status(g) = "excluded_base_mode";
            growth_rates.growth_rate_sigma(g) = NaN;
            growth_rates.growth_rate_standard_error(g) = NaN;
            growth_rates.phase_rate_omega(g) = NaN;
            growth_rates.quality_flags(g) = append_flag( ...
                growth_rates.quality_flags(g), "BASE_MODE");
        end
    end
    growth_rates = movevars(growth_rates, "series_id", "Before", 1);

    response = struct();
    response.series_id = response_id;
    response.series_folder = "";
    response.response_method = string(method);
    response.source_series = [a.series_id, b.series_id];
    response.frame_metrics = frame_metrics;
    response.mode_coefficients = table();
    response.growth_rates = growth_rates;
    response.fit_details = fit_details;
    response.time = time;
    response.frame_index = frame_metrics.frame_index;
    response.mode_m = a.mode_m;
    response.mode_n = a.mode_n;
    response.coefficient = coefficient;
    response.group_m = a.group_m;
    response.group_n = a.group_n;
    response.group_amplitude = group_amplitude;
    response.group_phase = group_phase;
    response.group_noise = group_noise;
    response.group_qR = group_qR;
    response.group_resolved = group_resolved;
    response.profiles = {};
    response.status = "ok";
    response.error_message = "";
end

function [available, b_indices] = matching_projection_states(a, b, selected_a, time, cfg)
    available = isfield(a, "projection_states") && isfield(b, "projection_states");
    b_indices = zeros(numel(time), 1);
    if ~available
        return
    end
    a_indices = find(selected_a);
    for k = 1:numel(time)
        [distance, b_index] = min(abs(b.time - time(k)));
        scale = max(1, abs(time(k)));
        if distance > cfg.pairing.time_tolerance * scale
            available = false;
            return
        end
        state_a = a.projection_states{a_indices(k)};
        state_b = b.projection_states{b_index};
        if isempty(state_a) || isempty(state_b) ...
                || ~isequal(state_a.vertex_ids, state_b.vertex_ids)
            available = false;
            return
        end
        b_indices(k) = b_index;
    end
end

function [coefficient, group_qR, group_noise] = joint_pair_projection( ...
        a, b, a_indices, b_indices, denominator, cfg)
    n_frames = numel(a_indices);
    n_modes = numel(a.mode_m);
    n_groups = numel(a.group_m);
    coefficient = complex(nan(n_frames, n_modes));
    group_qR = nan(n_frames, n_groups);
    group_noise = nan(n_frames, n_groups);

    for k = 1:n_frames
        state_a = a.projection_states{a_indices(k)};
        state_b = b.projection_states{b_indices(k)};
        common_x = 0.5 * (state_a.x + state_b.x);
        common_phi = angle(exp(1i * state_a.phi) + exp(1i * state_b.phi));
        common_area = 0.5 * (state_a.vertex_area + state_b.vertex_area);
        x_min = 0.5 * (state_a.x_min + state_b.x_min);
        x_max = 0.5 * (state_a.x_max + state_b.x_max);
        length_core = x_max - x_min;
        xi = (common_x - 0.5 * (x_min + x_max)) / length_core;
        W = common_area .* paired_axial_window(xi, cfg);
        radius = 0.5 * (state_a.radius + state_b.radius);
        delta_eta = (state_a.r - state_b.r) / (denominator * radius);

        [coeff, reconstructed] = paired_project( ...
            common_phi, xi, delta_eta, W, a.mode_m, a.mode_n, cfg);
        coefficient(k, :) = coeff.';
        residual = delta_eta - reconstructed;
        effective_vertices = sum(W) ^ 2 / max(sum(W .^ 2), eps);
        residual_rms = sqrt(sum(W .* abs(residual) .^ 2) / max(sum(W), eps));
        coefficient_noise = max(residual_rms / sqrt(max(effective_vertices, 1)), ...
            cfg.growth.absolute_noise_floor);

        for g = 1:n_groups
            primary = find(a.mode_m == a.group_m(g) & a.mode_n == a.group_n(g), 1);
            conjugate = find(a.mode_m == -a.group_m(g) & a.mode_n == -a.group_n(g), 1);
            pair_count = 1 + (primary ~= conjugate);
            group_noise(k, g) = sqrt(pair_count) * coefficient_noise;
            group_qR(k, g) = 2 * pi * abs(a.group_n(g)) * radius / length_core;
        end
    end
end

function [coeff, reconstructed] = paired_project( ...
        phi, xi, eta, W, mode_m, mode_n, cfg)
    Phi = exp(1i * (phi * mode_m.' + 2 * pi * xi * mode_n.'));
    if cfg.modes.projection == "quadrature_nudft"
        coeff = Phi' * (W .* eta) / sum(W);
    else
        sqrt_W = sqrt(W);
        weighted_basis = Phi .* sqrt_W;
        weighted_eta = eta .* sqrt_W;
        [U, S, V] = svd(weighted_basis, "econ");
        singular_values = diag(S);
        lambda = cfg.modes.regularization * singular_values(1) ^ 2;
        filter = singular_values ./ (singular_values .^ 2 + lambda);
        coeff = V * (filter .* (U' * weighted_eta));
    end
    reconstructed = real(Phi * coeff);
end

function window = paired_axial_window(xi, cfg)
    if cfg.window.type == "none" || cfg.window.alpha <= 0
        window = ones(size(xi));
        return
    end
    u = xi + 0.5;
    alpha = min(max(cfg.window.alpha, 0), 1);
    window = ones(size(u));
    outside = u < 0 | u > 1;
    left = u >= 0 & u < alpha / 2;
    right = u > 1 - alpha / 2 & u <= 1;
    window(left) = 0.5 * (1 + cos(pi * (2 * u(left) / alpha - 1)));
    window(right) = 0.5 * (1 + cos(pi * (2 * u(right) / alpha - 2 / alpha + 1)));
    window(outside) = 0;
end

function verify_compatible_modes(a, b)
    if ~isequal(a.mode_m, b.mode_m) || ~isequal(a.mode_n, b.mode_n) ...
            || ~isequal(a.group_m, b.group_m) || ~isequal(a.group_n, b.group_n)
        error("CylinderMode:PairingModes", ...
            "Related series were not analyzed with identical mode sets.");
    end
end

function values = interpolate_matrix(source_time, source_values, query_time)
    values = interp1(source_time, source_values, query_time, "linear");
end

function epsilon = pairing_epsilon(cfg, series_index, n_series)
    values = cfg.pairing.epsilon;
    if isempty(values)
        error("CylinderMode:PairingEpsilon", ...
            "Set cfg.pairing.epsilon for plus/minus analysis.");
    elseif isscalar(values)
        epsilon = values;
    elseif numel(values) == n_series
        epsilon = values(series_index);
    else
        error("CylinderMode:PairingEpsilon", ...
            "Pairing epsilon must be scalar or have one value per series.");
    end
    if ~isfinite(epsilon) || epsilon <= 0
        error("CylinderMode:PairingEpsilon", "Pairing epsilon must be positive.");
    end
end

function value = sanitize_name(value)
    value = regexprep(string(value), "[^A-Za-z0-9_+-]", "_");
end

function value = append_flag(value, flag)
    if strlength(value) == 0
        value = string(flag);
    elseif ~contains(";" + value + ";", ";" + flag + ";")
        value = value + ";" + flag;
    end
end
