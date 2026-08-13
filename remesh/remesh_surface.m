function [P2, M2, info] = remesh_surface(P, M, target_h, opts)
%REMESH_SURFACE Remesh a surface through the legacy or MMGS backend.
%
%   [P2,M2,info] = remesh_surface(P,M,target_h,opts)
%
% opts.backend is "legacy" (default) or "mmgs". The legacy backend uses
% a single target length; MMGS consumes the full per-vertex target_h field.

    if nargin < 3 || isempty(target_h)
        geo = Geometry(M, P);
        target_h = mean(geo.he_length);
    end
    if nargin < 4 || isempty(opts)
        opts = struct();
    end
    opts = complete_options(opts);

    P = double(P);
    M = double(M);
    if isscalar(target_h)
        target_h = repmat(double(target_h), size(P, 1), 1);
    else
        target_h = double(target_h(:));
    end
    validateattributes(P, {'double'}, {'2d', 'ncols', 3, 'finite', 'real'});
    validateattributes(M, {'double'}, {'2d', 'ncols', 3, 'integer', ...
        '>=', 1, '<=', size(P, 1)});
    validateattributes(target_h, {'double'}, {'column', 'numel', size(P, 1), ...
        'finite', 'positive'});

    before = mesh_surface_diagnostics(M, P);
    if ~before.valid
        error("remesh_surface:InvalidInput", ...
            "Input mesh failed validation: %s", before.quality_flags);
    end

    timer = tic;
    backend = lower(string(opts.backend));
    if backend == "legacy"
        require_legacy_backend();
        scalar_h = opts.legacy_target_statistic(target_h);
        [M2, P2] = remeshing(int32(M), P, int32([]), scalar_h, ...
            int32(opts.legacy_iterations));
        M2 = double(M2);
        P2 = double(P2);
        status = struct("success", true, "return_code", 0);
    elseif backend == "mmgs"
        mex_file = require_mmgs_backend();
        addpath(fileparts(mex_file));
        [P2, M2, status] = run_mmgs_with_complexity_control( ...
            P, M, target_h, before.area, opts);
        P2 = double(P2);
        M2 = double(M2);
    else
        error("remesh_surface:UnknownBackend", ...
            "Unknown backend '%s'. Use 'legacy' or 'mmgs'.", backend);
    end
    runtime = toc(timer);

    after = mesh_surface_diagnostics(M2, P2, struct("P", P, "M", M));
    if ~after.valid
        error("remesh_surface:InvalidOutput", ...
            "%s returned an invalid mesh: %s", upper(backend), after.quality_flags);
    end
    target_h_output = nearest_old_vertex_field(P, target_h, P2);
    [edge_ratio_min, edge_ratio_mean, edge_ratio_max] = ...
        target_edge_ratio(M2, P2, target_h_output);

    info = struct();
    info.backend = backend;
    info.runtime_seconds = runtime;
    info.status = status;
    info.before = before;
    info.raw_after = after;
    info.rel_area_error_raw = abs(after.area - before.area) / before.area;
    info.rel_volume_error_raw = abs(after.signed_volume - before.signed_volume) ...
        / max(abs(before.signed_volume), eps);
    info.target_h_min = min(target_h);
    info.target_h_mean = mean(target_h);
    info.target_h_max = max(target_h);
    info.edge_target_ratio_min = edge_ratio_min;
    info.edge_target_ratio_mean = edge_ratio_mean;
    info.edge_target_ratio_max = edge_ratio_max;
end

function opts = complete_options(opts)
    opts = with_default(opts, "backend", "legacy");
    opts = with_default(opts, "legacy_iterations", 20);
    opts = with_default(opts, "legacy_target_statistic", @mean);
    opts = with_default(opts, "metric_gradation", 1.3);
    opts = with_default(opts, "hausdorff_factor", 0.25);
    opts = with_default(opts, "mmgs_size_scale", 1.5);
    opts = with_default(opts, "match_target_complexity", true);
    opts = with_default(opts, "complexity_max_passes", 6);
    opts = with_default(opts, "vertex_count_tolerance", 0.08);
    opts = with_default(opts, "mmgs_size_scale_min", 0.4);
    opts = with_default(opts, "mmgs_size_scale_max", 4.0);
    opts = with_default(opts, "mmgs", struct());
    validateattributes(opts.mmgs_size_scale, {'numeric'}, ...
        {'scalar', 'finite', 'positive'});
    validateattributes(opts.complexity_max_passes, {'numeric'}, ...
        {'scalar', 'integer', '>=', 1});
    validateattributes(opts.vertex_count_tolerance, {'numeric'}, ...
        {'scalar', 'finite', 'nonnegative'});
    validateattributes(opts.mmgs_size_scale_min, {'numeric'}, ...
        {'scalar', 'finite', 'positive'});
    validateattributes(opts.mmgs_size_scale_max, {'numeric'}, ...
        {'scalar', 'finite', '>=', opts.mmgs_size_scale_min});
end

function [P_best, M_best, status_best] = run_mmgs_with_complexity_control( ...
        P, M, target_h, surface_area, opts)
    target_count = max(4, round(2 * surface_area / sqrt(3) ...
        * area_weighted_inverse_square(P, M, target_h) + 2));
    scale = opts.mmgs_size_scale;
    best_error = inf;
    passes = 0;
    lower_scale = NaN;
    upper_scale = NaN;

    for pass = 1:opts.complexity_max_passes
        passes = pass;
        metric_h = scale * target_h;
        mmgs_options = make_mmgs_options(metric_h, target_h, opts);
        [P_candidate, M_candidate, status_candidate] = ...
            mmgs_remesh_mex(P, M, metric_h, mmgs_options);
        count = size(P_candidate, 1);
        relative_error = abs(count - target_count) / target_count;
        if relative_error < best_error
            best_error = relative_error;
            P_best = P_candidate;
            M_best = M_candidate;
            status_best = status_candidate;
            best_scale = scale;
            best_count = count;
        end
        if ~opts.match_target_complexity ...
                || relative_error <= opts.vertex_count_tolerance
            break
        end
        if count > target_count
            lower_scale = scale;
        else
            upper_scale = scale;
        end
        if isfinite(lower_scale) && isfinite(upper_scale)
            scale = sqrt(lower_scale * upper_scale);
        else
            % A full log-space count correction converges faster than the
            % idealized square-root law when geometric tolerances dominate.
            scale = scale * count / target_count;
        end
        scale = min(max(scale, opts.mmgs_size_scale_min), ...
            opts.mmgs_size_scale_max);
    end

    status_best.complexity_passes = passes;
    status_best.target_vertex_count = target_count;
    status_best.output_vertex_count = best_count;
    status_best.vertex_count_relative_error = best_error;
    status_best.mmgs_size_scale = best_scale;
end

function options = make_mmgs_options(metric_h, target_h, opts)
    options = opts.mmgs;
    default_hmin = min(metric_h);
    default_hmax = max(metric_h);
    if default_hmax <= default_hmin * (1 + 1e-12)
        default_hmin = 0.99 * default_hmin;
        default_hmax = 1.01 * default_hmax;
    end
    options.hmin = with_default_value(options, "hmin", default_hmin);
    options.hmax = with_default_value(options, "hmax", default_hmax);
    options.hgrad = with_default_value(options, "hgrad", opts.metric_gradation);
    options.hausd = with_default_value(options, "hausd", ...
        opts.hausdorff_factor * min(target_h));
end

function mean_inverse_square = area_weighted_inverse_square(P, M, h)
    p1 = P(M(:, 1), :);
    p2 = P(M(:, 2), :);
    p3 = P(M(:, 3), :);
    face_area = 0.5 * vecnorm(cross(p2 - p1, p3 - p1, 2), 2, 2);
    vertex_area = accumarray(M(:), repmat(face_area / 3, 3, 1), ...
        [size(P, 1), 1], @sum, 0);
    mean_inverse_square = sum(vertex_area ./ h.^2) / sum(vertex_area);
end

function require_legacy_backend()
    if exist("remeshing", "file") == 0
        legacy_dir = fullfile(fileparts(mfilename("fullpath")), "isoremesh");
        addpath(legacy_dir);
    end
    if exist("remeshing", "file") == 0
        error("remesh_surface:LegacyUnavailable", ...
            "Legacy remeshing.mexa64 was not found.");
    end
end

function mex_file = require_mmgs_backend()
    backend_dir = fullfile(fileparts(mfilename("fullpath")), "mmgs");
    mex_file = fullfile(backend_dir, "mmgs_remesh_mex.mexa64");
    if ~isfile(mex_file)
        error("remesh_surface:MMGSUnavailable", ...
            "MMGS MEX backend is missing. Run remesh/mmgs/build_mmgs_backend.sh.");
    end
end

function value = with_default_value(options, field, fallback)
    if isfield(options, field) && ~isempty(options.(field))
        value = options.(field);
    else
        value = fallback;
    end
end

function values = nearest_old_vertex_field(P_old, values_old, P_new)
    if exist("knnsearch", "file") ~= 0
        index = knnsearch(P_old, P_new, "K", 1);
    else
        index = zeros(size(P_new, 1), 1);
        for i = 1:size(P_new, 1)
            [~, index(i)] = min(sum((P_old - P_new(i, :)).^2, 2));
        end
    end
    values = values_old(index);
end

function [minimum, average, maximum] = target_edge_ratio(M, P, target_h)
    edges = unique(sort([M(:, [1, 2]); M(:, [2, 3]); M(:, [3, 1])], 2), "rows");
    length_value = vecnorm(P(edges(:, 2), :) - P(edges(:, 1), :), 2, 2);
    edge_target = 0.5 * (target_h(edges(:, 1)) + target_h(edges(:, 2)));
    ratio = length_value ./ edge_target;
    minimum = min(ratio);
    average = mean(ratio);
    maximum = max(ratio);
end

function opts = with_default(opts, field, value)
    if ~isfield(opts, field) || isempty(opts.(field))
        opts.(field) = value;
    end
end
