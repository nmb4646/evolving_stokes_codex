function [target_h, info] = compute_target_edge_length(geo, base_h, opts)
%COMPUTE_TARGET_EDGE_LENGTH Smooth curvature-adaptive isotropic size field.
%
%   h = compute_target_edge_length(geo, base_h, opts)
%
% By default the indicator is the largest absolute principal curvature.
% The field is smoothed, robustly normalized, clipped, and gradation-limited
% before being passed to a remeshing backend.

    if nargin < 2 || isempty(base_h)
        base_h = median(geo.he_length);
    end
    if nargin < 3 || isempty(opts)
        opts = struct();
    end
    opts = complete_options(opts);
    validateattributes(base_h, {'numeric'}, {'scalar', 'finite', 'positive'});

    if ~opts.adaptive
        target_h = repmat(double(base_h), geo.mesh.n_v, 1);
        info = make_info(zeros(geo.mesh.n_v, 1), target_h, base_h, opts, NaN);
        return
    end

    curvature = curvature_indicator(geo, opts);
    curvature = smooth_vertex_field(geo.F, curvature, ...
        opts.curvature_smoothing_iterations, opts.smoothing_alpha);
    curvature = max(curvature, 0);
    cap = local_percentile(curvature, opts.curvature_clip_percentile);
    if isfinite(cap) && cap > 0
        curvature = min(curvature, cap);
    end

    positive = curvature(curvature > opts.curvature_floor);
    if isempty(positive)
        reference = 1;
        ratio = zeros(size(curvature));
    else
        reference = median(positive);
        ratio = curvature / max(reference, eps);
    end
    response = ratio .^ opts.curvature_power;
    target_h = base_h * (1 + opts.curvature_weight) ...
        ./ (1 + opts.curvature_weight * response);

    hmin = opts.hmin_factor * base_h;
    hmax = opts.hmax_factor * base_h;
    target_h = min(max(target_h, hmin), hmax);
    target_h = smooth_vertex_field(geo.F, target_h, ...
        opts.target_smoothing_iterations, opts.smoothing_alpha);
    target_h = min(max(target_h, hmin), hmax);
    target_h = limit_gradation(geo.F, target_h, opts.gradation_ratio, ...
        opts.gradation_iterations);
    target_vertex_count = max(4, round(2 + (geo.mesh.n_v - 2) ...
        * (mean(geo.he_length) / base_h)^2));
    target_complexity = sqrt(3) / 2 * (target_vertex_count - 2);
    complexity_before_normalization = metric_complexity_ratio( ...
        target_h, geo.v_area, target_complexity);
    if opts.preserve_vertex_budget
        target_h = normalize_metric_complexity(target_h, geo.v_area, ...
            target_complexity, hmin, hmax);
        target_h = limit_gradation(geo.F, target_h, opts.gradation_ratio, ...
            opts.gradation_iterations);
        target_h = normalize_metric_complexity(target_h, geo.v_area, ...
            target_complexity, hmin, hmax);
    end

    info = make_info(curvature, target_h, base_h, opts, reference);
    info.target_vertex_count = target_vertex_count;
    info.metric_complexity_before_normalization = complexity_before_normalization;
    info.metric_complexity_ratio = metric_complexity_ratio( ...
        target_h, geo.v_area, target_complexity);
end

function opts = complete_options(opts)
    opts = with_default(opts, "adaptive", true);
    opts = with_default(opts, "curvature_measure", "max_abs_principal");
    opts = with_default(opts, "curvature", []);
    opts = with_default(opts, "curvature_weight", 1.0);
    opts = with_default(opts, "curvature_power", 1.0);
    opts = with_default(opts, "curvature_floor", 1e-14);
    opts = with_default(opts, "curvature_clip_percentile", 98);
    opts = with_default(opts, "curvature_smoothing_iterations", 12);
    opts = with_default(opts, "target_smoothing_iterations", 2);
    opts = with_default(opts, "smoothing_alpha", 0.5);
    opts = with_default(opts, "hmin_factor", 0.35);
    opts = with_default(opts, "hmax_factor", 2.0);
    opts = with_default(opts, "gradation_ratio", 1.3);
    opts = with_default(opts, "gradation_iterations", 20);
    opts = with_default(opts, "preserve_vertex_budget", true);
    if opts.hmin_factor <= 0 || opts.hmax_factor < opts.hmin_factor
        error("Require 0 < hmin_factor <= hmax_factor.");
    end
    if opts.gradation_ratio < 1
        error("gradation_ratio must be at least one.");
    end
end

function ratio = metric_complexity_ratio(h, vertex_area, target_complexity)
% Element count scales approximately with integral(1/h^2 dA).
    vertex_area = max(vertex_area(:), 0);
    ratio = sum(vertex_area ./ h(:).^2) / max(target_complexity, realmin);
end

function h = normalize_metric_complexity(h, vertex_area, target, hmin, hmax)
    lower_scale = hmin / max(max(h), realmin);
    upper_scale = hmax / max(min(h), realmin);
    for iteration = 1:60
        scale = sqrt(lower_scale * upper_scale);
        candidate = min(max(scale * h, hmin), hmax);
        complexity = sum(vertex_area ./ candidate.^2);
        if complexity > target
            lower_scale = scale;
        else
            upper_scale = scale;
        end
    end
    h = min(max(sqrt(lower_scale * upper_scale) * h, hmin), hmax);
end

function curvature = curvature_indicator(geo, opts)
    if ~isempty(opts.curvature)
        curvature = double(opts.curvature);
        if ismatrix(curvature) && size(curvature, 2) > 1
            curvature = vecnorm(curvature, 2, 2);
        end
        curvature = curvature(:);
        if numel(curvature) ~= geo.mesh.n_v
            error("opts.curvature must contain one value per vertex.");
        end
        curvature = abs(curvature);
        return
    end

    vertex_area = max(geo.v_area(:), realmin);
    H = geo.v_mean_curvature(:) ./ vertex_area;
    K = geo.v_gaussian_curvature(:) ./ vertex_area;
    discriminant = sqrt(max(H.^2 - K, 0));
    k1 = H + discriminant;
    k2 = H - discriminant;
    measure = lower(string(opts.curvature_measure));
    if measure == "max_abs_principal" || measure == "principal"
        curvature = max(abs(k1), abs(k2));
    elseif measure == "mean"
        curvature = abs(H);
    elseif measure == "gaussian"
        curvature = sqrt(abs(K));
    else
        error("Unknown curvature_measure '%s'.", measure);
    end
end

function values = smooth_vertex_field(M, values, iterations, alpha)
    if iterations <= 0 || alpha <= 0
        return
    end
    n_vertices = numel(values);
    edges = unique(sort([M(:, [1, 2]); M(:, [2, 3]); M(:, [3, 1])], 2), "rows");
    adjacency = sparse([edges(:, 1); edges(:, 2)], ...
        [edges(:, 2); edges(:, 1)], 1, n_vertices, n_vertices);
    degree = full(sum(adjacency, 2));
    for iteration = 1:iterations
        neighbor_mean = (adjacency * values) ./ max(degree, 1);
        values = (1 - alpha) * values + alpha * neighbor_mean;
    end
end

function h = limit_gradation(M, h, ratio, iterations)
    if ratio <= 1 || iterations <= 0
        return
    end
    edges = unique(sort([M(:, [1, 2]); M(:, [2, 3]); M(:, [3, 1])], 2), "rows");
    for iteration = 1:iterations
        old_h = h;
        for direction = 1:2
            if direction == 1
                source = edges(:, 1);
                destination = edges(:, 2);
            else
                source = edges(:, 2);
                destination = edges(:, 1);
            end
            allowed = ratio * h(source);
            proposed = accumarray(destination, allowed, [numel(h), 1], @min, inf);
            h = min(h, proposed);
        end
        if max(abs(log(h ./ old_h))) < 1e-12
            break
        end
    end
end

function info = make_info(curvature, target_h, base_h, opts, reference)
    info = struct();
    info.base_h = base_h;
    info.adaptive = logical(opts.adaptive);
    info.curvature_measure = string(opts.curvature_measure);
    info.curvature_reference = reference;
    info.curvature_min = min(curvature);
    info.curvature_median = median(curvature);
    info.curvature_max = max(curvature);
    info.h_min = min(target_h);
    info.h_median = median(target_h);
    info.h_mean = mean(target_h);
    info.h_max = max(target_h);
    info.hmin_bound = opts.hmin_factor * base_h;
    info.hmax_bound = opts.hmax_factor * base_h;
end

function value = local_percentile(values, percentile)
    values = sort(values(isfinite(values)));
    if isempty(values)
        value = NaN;
        return
    end
    percentile = min(max(percentile, 0), 100);
    position = 1 + (numel(values) - 1) * percentile / 100;
    lower_index = floor(position);
    upper_index = ceil(position);
    fraction = position - lower_index;
    value = (1 - fraction) * values(lower_index) + fraction * values(upper_index);
end

function opts = with_default(opts, field, value)
    if ~isfield(opts, field) || isempty(opts.(field))
        opts.(field) = value;
    end
end
