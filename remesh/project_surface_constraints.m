function [P, geo, info] = project_surface_constraints(P, M, targets, opts)
%PROJECT_SURFACE_CONSTRAINTS Mass-weighted area/volume projection.
%
%   [P2, geo2, info] = project_surface_constraints(P, M, targets, opts)
%
% targets.area and targets.volume specify the desired polyhedral values.
% The default correction minimizes a lumped-area approximation to
% 0.5*integral(|delta x|^2 dA), simultaneously for all enabled constraints.
% Set opts.weighting="euclidean" only for comparison with the old projector.

    if nargin < 3 || isempty(targets)
        targets = struct();
    end
    if nargin < 4 || isempty(opts)
        opts = struct();
    end
    opts = complete_options(opts, targets);

    M = double(M);
    P = double(P);
    validateattributes(P, {'double'}, {'2d', 'ncols', 3, 'finite', 'real'});
    validateattributes(M, {'numeric'}, {'2d', 'ncols', 3, 'integer', ...
        '>=', 1, '<=', size(P, 1)});

    geo = Geometry(M, P);
    if ~isfield(targets, "area") || isempty(targets.area)
        targets.area = geo.area;
    end
    if ~isfield(targets, "volume") || isempty(targets.volume)
        targets.volume = geo.volume;
    end
    validateattributes(targets.area, {'numeric'}, {'scalar', 'finite', 'positive'});
    validateattributes(targets.volume, {'numeric'}, {'scalar', 'finite', 'nonzero'});

    P_initial = P;
    geo_initial = geo;
    max_states = opts.max_iter + 1;
    history_area = nan(max_states, 1);
    history_volume = nan(max_states, 1);
    history_rel_area = nan(max_states, 1);
    history_rel_volume = nan(max_states, 1);
    history_step = nan(opts.max_iter, 1);
    history_rcond = nan(opts.max_iter, 1);
    history_alpha = nan(opts.max_iter, 1);
    converged = false;
    reason = "MAX_ITERATIONS";
    n_steps = 0;

    for iter = 0:opts.max_iter
        geo = Geometry(M, P);
        rel_area = abs(geo.area - targets.area) / abs(targets.area);
        rel_volume = abs(geo.volume - targets.volume) / abs(targets.volume);
        history_area(iter + 1) = geo.area;
        history_volume(iter + 1) = geo.volume;
        history_rel_area(iter + 1) = rel_area;
        history_rel_volume(iter + 1) = rel_volume;

        area_ok = ~opts.preserve_area || rel_area <= opts.tol_area;
        volume_ok = ~opts.preserve_volume || rel_volume <= opts.tol_volume;
        if area_ok && volume_ok
            converged = true;
            reason = "CONVERGED";
            break
        end
        if iter == opts.max_iter
            break
        end

        [gradients, rhs, labels] = active_constraint_system(geo, M, P, targets, opts);
        vertex_weight = projection_weight(geo, opts);
        n_constraints = numel(gradients);
        weighted_gradients = cell(n_constraints, 1);
        schur = zeros(n_constraints);
        for constraint = 1:n_constraints
            weighted_gradients{constraint} = gradients{constraint} ./ vertex_weight;
        end
        for row = 1:n_constraints
            for column = row:n_constraints
                value = sum(dot(gradients{row}, weighted_gradients{column}, 2));
                schur(row, column) = value;
                schur(column, row) = value;
            end
        end

        schur_rcond = rcond(schur);
        history_rcond(iter + 1) = schur_rcond;
        if ~isfinite(schur_rcond) || schur_rcond < opts.rcond_tol
            coefficient = pinv(schur, opts.pinv_tolerance) * rhs;
        else
            coefficient = schur \ rhs;
        end
        if any(~isfinite(coefficient))
            reason = "NONFINITE_STEP";
            break
        end

        dP = zeros(size(P));
        for constraint = 1:n_constraints
            dP = dP + coefficient(constraint) * weighted_gradients{constraint};
        end
        dP = limit_step(dP, geo, opts);
        [alpha, accepted] = choose_step(P, M, dP, geo, targets, opts);
        history_alpha(iter + 1) = alpha;
        history_step(iter + 1) = max(vecnorm(alpha * dP, 2, 2));
        if ~accepted
            reason = "LINE_SEARCH_FAILED";
            break
        end
        P = P + alpha * dP;
        n_steps = n_steps + 1;

        if opts.verbose
            fprintf("Constraint projection iter %d (%s): relA %.3e, relV %.3e, alpha %.3g\n", ...
                iter + 1, strjoin(labels, "+"), rel_area, rel_volume, alpha);
        end
    end

    geo = Geometry(M, P);
    fallback_info = struct("applied", false, "converged", false, ...
        "iterations", 0);
    if ~converged && opts.volume_priority_on_failure ...
            && opts.preserve_area && opts.preserve_volume
        fallback_options = opts;
        fallback_options.preserve_area = false;
        fallback_options.preserve_volume = true;
        fallback_options.volume_priority_on_failure = false;
        fallback_options.max_iter = opts.volume_fallback_max_iter;
        [P, geo, volume_info] = project_surface_constraints( ...
            P, M, struct("volume", targets.volume), fallback_options);
        fallback_info.applied = true;
        fallback_info.converged = volume_info.converged;
        fallback_info.iterations = volume_info.iterations;
        reason = reason + ";VOLUME_PRIORITY_FALLBACK";
    end
    displacement = P - P_initial;
    displacement_norm = vecnorm(displacement, 2, 2);
    initial_weight = max(geo_initial.v_area(:), opts.mass_floor_absolute);
    info = struct();
    info.method = string(opts.weighting);
    info.preserve_area = opts.preserve_area;
    info.preserve_volume = opts.preserve_volume;
    info.target_area = targets.area;
    info.target_volume = targets.volume;
    info.area_initial = geo_initial.area;
    info.volume_initial = geo_initial.volume;
    info.area_final = geo.area;
    info.volume_final = geo.volume;
    info.rel_area_error_initial = abs(geo_initial.area - targets.area) / abs(targets.area);
    info.rel_volume_error_initial = abs(geo_initial.volume - targets.volume) / abs(targets.volume);
    info.rel_area_error = abs(geo.area - targets.area) / abs(targets.area);
    info.rel_volume_error = abs(geo.volume - targets.volume) / abs(targets.volume);
    info.max_displacement = max(displacement_norm);
    info.rms_displacement = sqrt(mean(displacement_norm.^2));
    info.area_weighted_rms_displacement = sqrt( ...
        sum(initial_weight .* displacement_norm.^2) / sum(initial_weight));
    info.displacement_objective = 0.5 * sum(initial_weight .* displacement_norm.^2);
    info.iterations = n_steps;
    info.converged = converged;
    info.reason = reason;
    info.volume_priority_fallback = fallback_info.applied;
    info.volume_fallback_converged = fallback_info.converged;
    info.volume_fallback_iterations = fallback_info.iterations;
    info.history = table((0:n_steps).', history_area(1:n_steps + 1), ...
        history_volume(1:n_steps + 1), history_rel_area(1:n_steps + 1), ...
        history_rel_volume(1:n_steps + 1), ...
        'VariableNames', {'iteration', 'area', 'volume', ...
        'rel_area_error', 'rel_volume_error'});
    info.step_max_displacement = history_step(1:n_steps);
    info.step_alpha = history_alpha(1:n_steps);
    info.schur_rcond = history_rcond(1:n_steps);
end

function opts = complete_options(opts, targets)
    opts = with_default(opts, "preserve_area", isfield(targets, "area"));
    opts = with_default(opts, "preserve_volume", isfield(targets, "volume"));
    opts = with_default(opts, "weighting", "mass");
    opts = with_default(opts, "tol_area", 1e-12);
    opts = with_default(opts, "tol_volume", 1e-12);
    opts = with_default(opts, "max_iter", 15);
    opts = with_default(opts, "line_search", true);
    opts = with_default(opts, "line_search_tau", 0.5);
    opts = with_default(opts, "line_search_max_iter", 16);
    opts = with_default(opts, "sufficient_decrease", 1e-4);
    opts = with_default(opts, "max_step_edge_fraction", 0.25);
    opts = with_default(opts, "mass_floor_relative", 1e-12);
    opts = with_default(opts, "mass_floor_absolute", realmin);
    opts = with_default(opts, "rcond_tol", 1e-12);
    opts = with_default(opts, "pinv_tolerance", 1e-14);
    opts = with_default(opts, "volume_priority_on_failure", true);
    opts = with_default(opts, "volume_fallback_max_iter", 8);
    opts = with_default(opts, "verbose", false);
    opts.weighting = lower(string(opts.weighting));
    if ~ismember(opts.weighting, ["mass", "euclidean"])
        error("opts.weighting must be 'mass' or 'euclidean'.");
    end
    if ~opts.preserve_area && ~opts.preserve_volume
        error("At least one of preserve_area or preserve_volume must be true.");
    end
end

function [gradients, rhs, labels] = active_constraint_system(geo, M, P, targets, opts)
    gradients = cell(0, 1);
    rhs = zeros(0, 1);
    labels = strings(0, 1);
    if opts.preserve_area
        scale = abs(targets.area);
        gradients{end + 1, 1} = surface_area_gradient(M, P) / scale;
        rhs(end + 1, 1) = (targets.area - geo.area) / scale;
        labels(end + 1, 1) = "area";
    end
    if opts.preserve_volume
        scale = abs(targets.volume);
        gradients{end + 1, 1} = signed_volume_gradient(M, P) / scale;
        rhs(end + 1, 1) = (targets.volume - geo.volume) / scale;
        labels(end + 1, 1) = "volume";
    end
end

function weight = projection_weight(geo, opts)
    if opts.weighting == "mass"
        area = geo.v_area(:);
        floor_value = max(opts.mass_floor_absolute, ...
            opts.mass_floor_relative * mean(area));
        weight = max(area, floor_value);
    else
        weight = ones(size(geo.V, 1), 1);
    end
end

function gradient = surface_area_gradient(M, P)
    p1 = P(M(:, 1), :);
    p2 = P(M(:, 2), :);
    p3 = P(M(:, 3), :);
    normal_vector = cross(p2 - p1, p3 - p1, 2);
    normal = normal_vector ./ max(vecnorm(normal_vector, 2, 2), realmin);
    g1 = 0.5 * cross(p2 - p3, normal, 2);
    g2 = 0.5 * cross(p3 - p1, normal, 2);
    g3 = 0.5 * cross(p1 - p2, normal, 2);
    gradient = accumulate_face_vectors(M, [g1; g2; g3], size(P, 1));
end

function gradient = signed_volume_gradient(M, P)
    p1 = P(M(:, 1), :);
    p2 = P(M(:, 2), :);
    p3 = P(M(:, 3), :);
    g1 = cross(p2, p3, 2) / 6;
    g2 = cross(p3, p1, 2) / 6;
    g3 = cross(p1, p2, 2) / 6;
    gradient = accumulate_face_vectors(M, [g1; g2; g3], size(P, 1));
end

function result = accumulate_face_vectors(M, values, n_vertices)
    vertices = M(:);
    result = zeros(n_vertices, 3);
    for dimension = 1:3
        result(:, dimension) = accumarray(vertices, values(:, dimension), ...
            [n_vertices, 1], @sum, 0);
    end
end

function dP = limit_step(dP, geo, opts)
    if isempty(opts.max_step_edge_fraction) || ~isfinite(opts.max_step_edge_fraction)
        return
    end
    maximum = opts.max_step_edge_fraction * mean(geo.he_length);
    actual = max(vecnorm(dP, 2, 2));
    if actual > maximum && maximum > 0
        dP = dP * (maximum / actual);
    end
end

function [alpha, accepted] = choose_step(P, M, dP, geo, targets, opts)
    alpha = 1;
    accepted = false;
    if ~opts.line_search
        accepted = trial_is_valid(P, M, dP);
        return
    end
    merit0 = constraint_merit(geo, targets, opts);
    for iteration = 1:opts.line_search_max_iter
        trial_step = alpha * dP;
        if trial_is_valid(P, M, trial_step)
            geo_trial = Geometry(M, P + trial_step);
            merit = constraint_merit(geo_trial, targets, opts);
            if merit <= (1 - opts.sufficient_decrease * alpha) * merit0
                accepted = true;
                return
            end
        end
        alpha = alpha * opts.line_search_tau;
    end
    alpha = 0;
end

function merit = constraint_merit(geo, targets, opts)
    residual = zeros(0, 1);
    if opts.preserve_area
        residual(end + 1, 1) = (geo.area - targets.area) / abs(targets.area);
    end
    if opts.preserve_volume
        residual(end + 1, 1) = (geo.volume - targets.volume) / abs(targets.volume);
    end
    merit = 0.5 * sum(residual.^2);
end

function valid = trial_is_valid(P, M, dP)
    trial = P + dP;
    old_normal = cross(P(M(:, 2), :) - P(M(:, 1), :), ...
        P(M(:, 3), :) - P(M(:, 1), :), 2);
    new_normal = cross(trial(M(:, 2), :) - trial(M(:, 1), :), ...
        trial(M(:, 3), :) - trial(M(:, 1), :), 2);
    valid = all(isfinite(trial), "all") ...
        && all(vecnorm(new_normal, 2, 2) > realmin) ...
        && all(dot(old_normal, new_normal, 2) > 0);
end

function opts = with_default(opts, field, value)
    if ~isfield(opts, field) || isempty(opts.(field))
        opts.(field) = value;
    end
end
