function [M, P, info] = remeshing_curvature(M, P, curvature, target_edge_length, n_iter, opts)
%REMESHING_CURVATURE Curvature-weighted remeshing in pure MATLAB.
%
%   [M2, P2] = remeshing_curvature(M, P, curvature, target_edge_length, n_iter)
%   remeshes a closed triangular surface using the same broad loop as the
%   isotropic OpenMesh remesher:
%
%       split long edges
%       collapse short edges
%       equalize valences by edge flips
%       tangential relaxation
%       projection back to the input surface
%
%   The difference is that the target edge length is spatially varying.
%   High curvature gives smaller local target lengths. This file is
%   intentionally independent from isoremesh/remeshing.mexa64.
%
%   curvature may be vertex-located (#V), face-located (#F), or half-edge
%   located (3#F). Vector data is converted to row-wise norm.

    if nargin < 3
        curvature = [];
    end
    if nargin < 4
        target_edge_length = [];
    end
    if nargin < 5 || isempty(n_iter)
        n_iter = 8;
    end
    if nargin < 6 || isempty(opts)
        opts = struct();
    end
    if isstruct(n_iter)
        opts = n_iter;
        n_iter = default_opt(opts, "n_iter", 8);
    end

    M = double(M);
    P = double(P);
    n_iter = max(0, double(n_iter));

    opts = complete_options(opts, size(P, 1));
    if isempty(target_edge_length)
        target_edge_length = median_edge_length(M, P);
    end
    target_edge_length = double(target_edge_length);

    ref_M = M;
    ref_P = P;
    ref_projector = make_reference_projector(ref_M, ref_P, opts);

    curv_v = curvature_to_vertex(M, P, curvature, opts);
    curv_v = smooth_scalar_on_mesh(M, curv_v, ...
        opts.curvature_smooth_iterations, opts.curvature_smooth_alpha);
    target_h = target_lengths_from_curvature(M, P, curv_v, ...
        target_edge_length, opts);

    info.initial_vertex_count = size(P, 1);
    info.initial_face_count = size(M, 1);
    info.target_vertex_count = opts.target_vertex_count;
    info.target_edge_length = target_edge_length;
    info.total_splits = 0;
    info.total_collapses = 0;
    info.total_flips = 0;
    info.total_relaxations = 0;
    info.projection_failed = false;
    info.iteration = repmat(empty_iteration_info(), max(1, n_iter), 1);

    for iter = 1:n_iter
        n_before = size(P, 1);

        [M, P, curv_v, target_h, n_splits] = split_long_edges( ...
            M, P, curv_v, target_h, opts);

        [M, P, curv_v, target_h, n_collapses] = collapse_short_edges( ...
            M, P, curv_v, target_h, opts);

        n_flips = 0;
        for pass = 1:opts.flip_passes_per_iteration
            [M, pass_flips] = equalize_valences(M, P, target_h, opts);
            n_flips = n_flips + pass_flips;
            if pass_flips == 0
                break
            end
        end
        M = orient_faces_consistently(M);

        if opts.relaxation_iterations_per_iteration > 0
            P = tangential_relaxation(M, P, ...
                opts.relaxation_iterations_per_iteration, opts.relaxation_alpha);
        end

        if opts.project_to_input
            [P, projection_failed] = project_to_reference(P, ref_M, ref_P, ref_projector, opts);
            info.projection_failed = info.projection_failed || projection_failed;
        end

        [~, edge_length, h_edge, ratio] = edge_stats(M, P, target_h);
        info.iteration(iter).n_vertices = size(P, 1);
        info.iteration(iter).n_faces = size(M, 1);
        info.iteration(iter).n_splits = n_splits;
        info.iteration(iter).n_collapses = n_collapses;
        info.iteration(iter).n_flips = n_flips;
        info.iteration(iter).mean_edge_ratio = mean(edge_length ./ max(h_edge, eps));
        info.iteration(iter).max_edge_ratio = max(ratio);
        info.iteration(iter).min_edge_ratio = min(ratio);
        info.total_splits = info.total_splits + n_splits;
        info.total_collapses = info.total_collapses + n_collapses;
        info.total_flips = info.total_flips + n_flips;
        info.total_relaxations = info.total_relaxations + opts.relaxation_iterations_per_iteration;

        if opts.verbose
            fprintf("remeshing_curvature iter %d: V %d -> %d, split %d, collapse %d, flip %d, edge ratio [%g, %g]\n", ...
                iter, n_before, size(P, 1), n_splits, n_collapses, n_flips, ...
                info.iteration(iter).min_edge_ratio, info.iteration(iter).max_edge_ratio);
        end

        count_ok = abs(size(P, 1) - opts.target_vertex_count) <= opts.count_tolerance;
        length_ok = max(ratio) <= opts.split_factor && min(ratio) >= 0.5 * opts.collapse_factor;
        if count_ok && length_ok && n_splits == 0 && n_collapses == 0 && n_flips == 0
            info.iteration = info.iteration(1:iter);
            break
        end
    end

    for polish = 1:opts.polish_iterations
        n_flips = 0;
        for pass = 1:opts.flip_passes_per_iteration
            [M, pass_flips] = equalize_valences(M, P, target_h, opts);
            n_flips = n_flips + pass_flips;
            if pass_flips == 0
                break
            end
        end
        M = orient_faces_consistently(M);
        if opts.relaxation_iterations_per_iteration > 0
            P = tangential_relaxation(M, P, ...
                opts.relaxation_iterations_per_iteration, opts.relaxation_alpha);
        end
        if opts.project_to_input
            [P, projection_failed] = project_to_reference(P, ref_M, ref_P, ref_projector, opts);
            info.projection_failed = info.projection_failed || projection_failed;
        end
        info.total_flips = info.total_flips + n_flips;
    end

    M = orient_faces_consistently(M);
    [M, P, ~, target_h] = compact_vertices(M, P, curv_v, target_h);
    info.final_vertex_count = size(P, 1);
    info.final_face_count = size(M, 1);
    info.final_target_h_min = min(target_h);
    info.final_target_h_max = max(target_h);
    info.final_bad_edge_count = count_bad_edges(M);
end

function s = empty_iteration_info()
    s = struct( ...
        "n_vertices", 0, ...
        "n_faces", 0, ...
        "n_splits", 0, ...
        "n_collapses", 0, ...
        "n_flips", 0, ...
        "mean_edge_ratio", 0, ...
        "max_edge_ratio", 0, ...
        "min_edge_ratio", 0);
end

function opts = complete_options(opts, n_v)
    opts.curvature_location = string(default_opt(opts, "curvature_location", "auto"));
    opts.curvature_abs = default_opt(opts, "curvature_abs", true);
    opts.curvature_low_percentile = default_opt(opts, "curvature_low_percentile", 5);
    opts.curvature_high_percentile = default_opt(opts, "curvature_high_percentile", 95);
    opts.curvature_weight = default_opt(opts, "curvature_weight", 6);
    opts.curvature_power = default_opt(opts, "curvature_power", 1.25);
    opts.curvature_smooth_iterations = default_opt(opts, "curvature_smooth_iterations", 3);
    opts.curvature_smooth_alpha = default_opt(opts, "curvature_smooth_alpha", 0.35);
    opts.target_smooth_iterations = default_opt(opts, "target_smooth_iterations", 2);
    opts.target_smooth_alpha = default_opt(opts, "target_smooth_alpha", 0.35);
    opts.min_edge_length_factor = default_opt(opts, "min_edge_length_factor", 0.25);
    opts.max_edge_length_factor = default_opt(opts, "max_edge_length_factor", 2.5);
    opts.split_factor = default_opt(opts, "split_factor", 4 / 3);
    opts.collapse_factor = default_opt(opts, "collapse_factor", 4 / 5);
    opts.forced_collapse_factor = default_opt(opts, "forced_collapse_factor", 1.0);
    opts.collapse_max_edge_factor = default_opt(opts, "collapse_max_edge_factor", 4 / 3);
    opts.flip_max_edge_factor = default_opt(opts, "flip_max_edge_factor", 1.55);
    opts.target_vertex_count = default_opt(opts, "target_vertex_count", n_v);
    opts.count_tolerance = default_opt(opts, "count_tolerance", max(6, ceil(0.04 * n_v)));
    opts.max_vertex_count_factor = default_opt(opts, "max_vertex_count_factor", 1.35);
    opts.max_splits_per_iteration = default_opt(opts, "max_splits_per_iteration", max(50, ceil(0.25 * n_v)));
    opts.max_collapses_per_iteration = default_opt(opts, "max_collapses_per_iteration", max(50, ceil(0.35 * n_v)));
    opts.high_curvature_collapse_penalty = default_opt(opts, "high_curvature_collapse_penalty", 0.8);
    opts.min_normal_dot = default_opt(opts, "min_normal_dot", 0.02);
    opts.min_triangle_quality = default_opt(opts, "min_triangle_quality", 0.04);
    opts.collapse_quality_factor = default_opt(opts, "collapse_quality_factor", 0.25);
    opts.flip_quality_factor = default_opt(opts, "flip_quality_factor", 0.6);
    opts.flip_accept_equal_valence_if_quality_improves = default_opt(opts, "flip_accept_equal_valence_if_quality_improves", true);
    opts.flip_passes_per_iteration = default_opt(opts, "flip_passes_per_iteration", 6);
    opts.relaxation_iterations_per_iteration = default_opt(opts, "relaxation_iterations_per_iteration", 1);
    opts.relaxation_alpha = default_opt(opts, "relaxation_alpha", 0.85);
    opts.polish_iterations = default_opt(opts, "polish_iterations", 2);
    opts.project_to_input = default_opt(opts, "project_to_input", true);
    opts.projection_neighbors = default_opt(opts, "projection_neighbors", 20);
    opts.projection_max_neighbors = default_opt(opts, "projection_max_neighbors", 80);
    opts.verbose = default_opt(opts, "verbose", false);
end

function value = default_opt(opts, name, default_value)
    if isfield(opts, name) && ~isempty(opts.(name))
        value = opts.(name);
    else
        value = default_value;
    end
end

function curv_v = curvature_to_vertex(M, P, curvature, opts)
    n_v = size(P, 1);
    n_f = size(M, 1);
    n_he = 3 * n_f;

    if isempty(curvature)
        geo = Geometry(M, P);
        curv_v = geo.v_mean_curvature ./ max(geo.v_area, eps);
        return
    end

    if size(curvature, 2) > 1
        curv = vecnorm(double(curvature), 2, 2);
    else
        curv = double(curvature(:));
    end
    curv(~isfinite(curv)) = 0;
    if opts.curvature_abs
        curv = abs(curv);
    end

    location = lower(opts.curvature_location);
    if location == "auto"
        if numel(curv) == n_v
            location = "vertex";
        elseif numel(curv) == n_he
            location = "halfedge";
        elseif numel(curv) == n_f
            location = "face";
        else
            error("Could not infer curvature location from %d values for %d vertices and %d faces.", ...
                numel(curv), n_v, n_f);
        end
    end

    mesh = Mesh(M, P);
    if location == "vertex"
        if numel(curv) ~= n_v
            error("Vertex curvature must have one value per vertex.");
        end
        curv_v = curv;
    elseif location == "halfedge"
        if numel(curv) ~= n_he
            error("Half-edge curvature must have 3*#faces values.");
        end
        [sum_curv, n_neighbor] = mesh.halfedge_to_vertex(curv);
        curv_v = sum_curv ./ max(n_neighbor, 1);
    elseif location == "face"
        if numel(curv) ~= n_f
            error("Face curvature must have one value per face.");
        end
        [sum_curv, n_neighbor] = mesh.face_to_vertex(curv);
        curv_v = sum_curv ./ max(n_neighbor, 1);
    else
        error("Unknown curvature_location '%s'.", location);
    end
end

function target_h = target_lengths_from_curvature(M, P, curv_v, base_h, opts)
    q_lo = local_percentile(curv_v, opts.curvature_low_percentile);
    q_hi = local_percentile(curv_v, opts.curvature_high_percentile);
    if q_hi <= q_lo
        q_norm = zeros(size(curv_v));
    else
        q_norm = (curv_v - q_lo) ./ (q_hi - q_lo);
        q_norm = min(max(q_norm, 0), 1);
    end

    density = 1 + opts.curvature_weight * q_norm .^ opts.curvature_power;
    target_h = base_h ./ sqrt(density);

    [edges, ~, ~, ~] = edge_stats(M, P, target_h);
    edge_h = 0.5 * (target_h(edges(:, 1)) + target_h(edges(:, 2)));
    if mean(edge_h) > 0
        target_h = target_h * (base_h / mean(edge_h));
    end

    target_h = smooth_scalar_on_mesh(M, target_h, ...
        opts.target_smooth_iterations, opts.target_smooth_alpha);
    target_h = min(max(target_h, opts.min_edge_length_factor * base_h), ...
        opts.max_edge_length_factor * base_h);
end

function [M, P, curv_v, target_h, n_splits] = split_long_edges(M, P, curv_v, target_h, opts)
    n_splits = 0;
    max_vertices = ceil(opts.max_vertex_count_factor * opts.target_vertex_count);
    max_splits = opts.max_splits_per_iteration;

    [edges, ~, ~, ratio, occ] = edge_stats(M, P, target_h);
    candidates = find(occ == 2 & ratio > opts.split_factor);
    if isempty(candidates)
        return
    end

    [~, order] = sort(ratio(candidates), "descend");
    candidates = candidates(order);

    for k = 1:numel(candidates)
        if n_splits >= max_splits || size(P, 1) >= max_vertices
            break
        end
        edge = edges(candidates(k), :);
        [M, P, curv_v, target_h, did_split] = split_one_edge( ...
            M, P, curv_v, target_h, edge(1), edge(2));
        if did_split
            n_splits = n_splits + 1;
        end
    end
end

function [M, P, curv_v, target_h, n_collapses] = collapse_short_edges(M, P, curv_v, target_h, opts)
    n_collapses = 0;
    max_collapses = opts.max_collapses_per_iteration;

    while n_collapses < max_collapses
        [edges, ~, h_edge, ratio, occ] = edge_stats(M, P, target_h);
        count_too_high = size(P, 1) > opts.target_vertex_count + opts.count_tolerance;
        if count_too_high
            threshold = opts.forced_collapse_factor;
        else
            threshold = opts.collapse_factor;
        end

        candidates = find(occ == 2 & ratio < threshold);
        if isempty(candidates)
            break
        end

        target_density = 1 ./ max(h_edge, eps) .^ 2;
        target_density = target_density / max(target_density);
        score = ratio + opts.high_curvature_collapse_penalty * target_density;
        [~, order] = sort(score(candidates), "ascend");
        candidates = candidates(order);

        did_any = false;
        for k = 1:numel(candidates)
            edge = edges(candidates(k), :);
            [M_trial, P_trial, curv_trial, h_trial, did_collapse] = collapse_one_edge( ...
                M, P, curv_v, target_h, edge(1), edge(2), opts);
            if did_collapse
                M = M_trial;
                P = P_trial;
                curv_v = curv_trial;
                target_h = h_trial;
                did_any = true;
                break
            end
        end

        if ~did_any
            break
        end
        n_collapses = n_collapses + 1;
    end
end

function [M, P, curv_v, target_h, did_split] = split_one_edge(M, P, curv_v, target_h, a, b)
    did_split = false;
    face_ids = find(sum(M == a | M == b, 2) == 2);
    if numel(face_ids) ~= 2
        return
    end

    new_vertex = size(P, 1) + 1;
    P(new_vertex, :) = 0.5 * (P(a, :) + P(b, :));
    curv_v(new_vertex, 1) = 0.5 * (curv_v(a) + curv_v(b));
    target_h(new_vertex, 1) = 0.5 * (target_h(a) + target_h(b));

    new_faces = zeros(4, 3);
    out = 1;
    for i = 1:numel(face_ids)
        tri = M(face_ids(i), :);
        [f1, f2, ok] = split_face_by_edge(tri, a, b, new_vertex);
        if ~ok
            continue
        end
        new_faces(out, :) = f1;
        new_faces(out + 1, :) = f2;
        out = out + 2;
    end
    if out ~= 5
        P(new_vertex, :) = [];
        curv_v(new_vertex) = [];
        target_h(new_vertex) = [];
        return
    end

    M(face_ids, :) = [];
    M = [M; new_faces];
    did_split = true;
end

function [f1, f2, ok] = split_face_by_edge(tri, a, b, new_vertex)
    ok = true;
    if edge_matches(tri(1), tri(2), a, b)
        x = tri(1); y = tri(2); z = tri(3);
    elseif edge_matches(tri(2), tri(3), a, b)
        x = tri(2); y = tri(3); z = tri(1);
    elseif edge_matches(tri(3), tri(1), a, b)
        x = tri(3); y = tri(1); z = tri(2);
    else
        ok = false;
        f1 = [0, 0, 0];
        f2 = [0, 0, 0];
        return
    end
    f1 = [x, new_vertex, z];
    f2 = [new_vertex, y, z];
end

function tf = edge_matches(x, y, a, b)
    tf = (x == a && y == b) || (x == b && y == a);
end

function [M, P, curv_v, target_h, did_collapse] = collapse_one_edge(M, P, curv_v, target_h, a, b, opts)
    did_collapse = false;
    n_v = size(P, 1);
    if ~passes_link_condition(M, n_v, a, b)
        return
    end

    A = vertex_adjacency(M, n_v);
    valence = full(sum(A, 2));
    if target_h(a) < target_h(b) || (target_h(a) == target_h(b) && valence(a) >= valence(b))
        keep = a;
        drop = b;
    else
        keep = b;
        drop = a;
    end

    neighbors = find(A(drop, :));
    neighbors = setdiff(neighbors, keep);
    for i = 1:numel(neighbors)
        n = neighbors(i);
        h_new = 0.5 * (target_h(keep) + target_h(n));
        if norm(P(keep, :) - P(n, :)) > opts.collapse_max_edge_factor * h_new
            return
        end
    end

    edge_faces = find(sum(M == keep | M == drop, 2) == 2);
    local_face_mask = any(M == keep | M == drop, 2) & ~ismember((1:size(M, 1))', edge_faces);
    old_tri = M(local_face_mask, :);
    old_cross = triangle_cross(P, old_tri);
    old_quality = triangle_quality(P, old_tri);

    P_trial = P;
    M_trial = M;
    M_trial(M_trial == drop) = keep;

    degenerate = M_trial(:, 1) == M_trial(:, 2) ...
        | M_trial(:, 2) == M_trial(:, 3) ...
        | M_trial(:, 3) == M_trial(:, 1);
    M_trial(degenerate, :) = [];
    if isempty(M_trial)
        return
    end

    local_after = old_tri;
    local_after(local_after == drop) = keep;
    local_degenerate = local_after(:, 1) == local_after(:, 2) ...
        | local_after(:, 2) == local_after(:, 3) ...
        | local_after(:, 3) == local_after(:, 1);
    local_after(local_degenerate, :) = [];
    old_cross(local_degenerate, :) = [];
    old_quality(local_degenerate) = [];

    if ~isempty(local_after)
        new_cross = triangle_cross(P_trial, local_after);
        old_norm = vecnorm(old_cross, 2, 2);
        new_norm = vecnorm(new_cross, 2, 2);
        normal_dot = dot(old_cross, new_cross, 2) ./ max(old_norm .* new_norm, eps);
        new_quality = triangle_quality(P_trial, local_after);
        if any(new_norm <= eps) || any(normal_dot < opts.min_normal_dot)
            return
        end
        min_old_quality = max(min(old_quality), opts.min_triangle_quality);
        if min(new_quality) < opts.collapse_quality_factor * min_old_quality ...
                || min(new_quality) < opts.min_triangle_quality
            return
        end
    end

    [~, unique_idx] = unique(sort(M_trial, 2), "rows", "stable");
    if numel(unique_idx) < size(M_trial, 1)
        return
    end

    curv_v(keep) = max(curv_v(keep), curv_v(drop));
    target_h(keep) = min(target_h(keep), target_h(drop));
    [M, P, curv_v, target_h] = compact_vertices(M_trial, P_trial, curv_v, target_h);
    did_collapse = true;
end

function [M, n_flips] = equalize_valences(M, P, target_h, opts)
    n_flips = 0;
    [edges, ~, ~, ratio, occ] = edge_stats(M, P, target_h);
    candidates = find(occ == 2 & ratio < opts.flip_max_edge_factor);
    if isempty(candidates)
        return
    end

    [~, order] = sort(abs(ratio(candidates) - 1), "descend");
    candidates = candidates(order);

    for k = 1:numel(candidates)
        edge = edges(candidates(k), :);
        [M_trial, did_flip] = flip_one_edge(M, P, target_h, edge(1), edge(2), opts);
        if did_flip
            M = M_trial;
            n_flips = n_flips + 1;
        end
    end
end

function [M, did_flip] = flip_one_edge(M, P, target_h, a, b, opts)
    did_flip = false;
    face_ids = find(sum(M == a | M == b, 2) == 2);
    if numel(face_ids) ~= 2
        return
    end

    tri1 = M(face_ids(1), :);
    tri2 = M(face_ids(2), :);
    c = tri1(~ismember(tri1, [a, b]));
    d = tri2(~ismember(tri2, [a, b]));
    if numel(c) ~= 1 || numel(d) ~= 1 || c == d
        return
    end

    if any(sum(M == c | M == d, 2) == 2)
        return
    end

    new_h = 0.5 * (target_h(c) + target_h(d));
    if norm(P(c, :) - P(d, :)) > opts.flip_max_edge_factor * new_h
        return
    end

    A = vertex_adjacency(M, size(P, 1));
    valence = full(sum(A, 2));
    verts = [a; b; c; d];
    deviation_pre = sum(abs(valence(verts) - 6));
    valence_after = valence;
    valence_after(a) = valence_after(a) - 1;
    valence_after(b) = valence_after(b) - 1;
    valence_after(c) = valence_after(c) + 1;
    valence_after(d) = valence_after(d) + 1;
    deviation_post = sum(abs(valence_after(verts) - 6));

    old_quality = triangle_quality(P, [tri1; tri2]);
    old_patch_normal = sum(triangle_cross(P, [tri1; tri2]), 1);
    if norm(old_patch_normal) <= eps
        return
    end

    [new_tri1, new_tri2] = orient_triangle_pair_to_normal( ...
        [c, d, b], [d, c, a], P, old_patch_normal);
    new_quality = triangle_quality(P, [new_tri1; new_tri2]);
    if any(new_quality < opts.min_triangle_quality)
        return
    end

    quality_ok = min(new_quality) >= opts.flip_quality_factor * max(min(old_quality), opts.min_triangle_quality);
    valence_improves = deviation_post < deviation_pre;
    quality_improves = min(new_quality) > min(old_quality) * 1.02;
    if ~quality_ok
        return
    end
    if ~valence_improves
        if ~(opts.flip_accept_equal_valence_if_quality_improves ...
                && deviation_post == deviation_pre && quality_improves)
            return
        end
    end

    M(face_ids(1), :) = new_tri1;
    M(face_ids(2), :) = new_tri2;
    did_flip = true;
end

function [tri1, tri2] = orient_triangle_pair_to_normal(tri1, tri2, P, normal_ref)
    c = sum(triangle_cross(P, [tri1; tri2]), 1);
    if dot(c, normal_ref) < 0
        tri1 = tri1([1, 3, 2]);
        tri2 = tri2([1, 3, 2]);
    end
end

function P = tangential_relaxation(M, P, n_steps, alpha)
    A = vertex_adjacency(M, size(P, 1));
    deg = max(sum(A, 2), 1);
    for step = 1:n_steps
        geo = Geometry(M, P);
        P_avg = (A * P) ./ deg;
        dP = P_avg - P;
        dP = dP - dot(dP, geo.v_normal, 2) .* geo.v_normal;
        P = P + alpha * dP;
    end
end

function ref_projector = make_reference_projector(ref_M, ref_P, opts)
    ref_projector = struct();
    if ~opts.project_to_input
        return
    end
    geo_ref = Geometry(ref_M, ref_P);
    ref_projector.kdtree = KDTreeSearcher(geo_ref.f_center);
end

function [P, projection_failed] = project_to_reference(P, ref_M, ref_P, ref_projector, opts)
    projection_failed = false;
    n_neighbor = opts.projection_neighbors;
    while n_neighbor <= opts.projection_max_neighbors
        [face, uv, ~, fail] = project(ref_P, ref_M, P, ref_projector.kdtree, n_neighbor);
        if ~fail
            P = interpolate(ref_M, face, uv, ref_P);
            return
        end
        n_neighbor = 2 * n_neighbor;
    end

    projection_failed = true;
    warning("remeshing_curvature:ProjectionFailed", ...
        "Projection failed up to %d neighboring faces; leaving relaxed points unprojected.", ...
        opts.projection_max_neighbors);
end

function ok = passes_link_condition(M, n_v, a, b)
    face_ids = find(sum(M == a | M == b, 2) == 2);
    if numel(face_ids) ~= 2
        ok = false;
        return
    end

    A = vertex_adjacency(M, n_v);
    na = find(A(a, :));
    nb = find(A(b, :));
    common_neighbors = intersect(na, nb);
    opposite = M(face_ids, :);
    opposite = opposite(opposite ~= a & opposite ~= b);
    ok = numel(common_neighbors) == 2 ...
        && numel(opposite) == 2 ...
        && isequal(sort(common_neighbors(:)), sort(opposite(:)));
end

function [edges, edge_length, h_edge, ratio, occ] = edge_stats(M, P, target_h)
    directed_edges = [M(:, [1, 2]); M(:, [2, 3]); M(:, [3, 1])];
    sorted_edges = sort(directed_edges, 2);
    [edges, ~, ic] = unique(sorted_edges, "rows");
    occ = accumarray(ic, 1, [size(edges, 1), 1]);
    dP = P(edges(:, 1), :) - P(edges(:, 2), :);
    edge_length = vecnorm(dP, 2, 2);
    h_edge = 0.5 * (target_h(edges(:, 1)) + target_h(edges(:, 2)));
    ratio = edge_length ./ max(h_edge, eps);
end

function x = smooth_scalar_on_mesh(M, x, n_steps, alpha)
    x = x(:);
    if n_steps <= 0 || alpha <= 0
        return
    end
    A = vertex_adjacency(M, numel(x));
    deg = max(sum(A, 2), 1);
    for step = 1:n_steps
        x_avg = (A * x) ./ deg;
        x = (1 - alpha) * x + alpha * x_avg;
    end
end

function A = vertex_adjacency(M, n_v)
    rows = [M(:, 1); M(:, 2); M(:, 3); M(:, 2); M(:, 3); M(:, 1)];
    cols = [M(:, 2); M(:, 3); M(:, 1); M(:, 1); M(:, 2); M(:, 3)];
    A = sparse(rows, cols, true, n_v, n_v);
end

function M = orient_faces_consistently(M)
    n_f = size(M, 1);
    directed_edges = [M(:, [1, 2]); M(:, [2, 3]); M(:, [3, 1])];
    face_id = repmat((1:n_f)', 3, 1);
    edge_dir = ones(size(directed_edges, 1), 1);
    edge_dir(directed_edges(:, 1) > directed_edges(:, 2)) = -1;

    [~, ~, edge_id] = unique(sort(directed_edges, 2), "rows");
    edge_faces = accumarray(edge_id, face_id, [], @(x) {x});
    edge_dirs = accumarray(edge_id, edge_dir, [], @(x) {x});

    adj_i = zeros(2 * numel(edge_faces), 1);
    adj_j = zeros(2 * numel(edge_faces), 1);
    adj_relation = zeros(2 * numel(edge_faces), 1);
    n_adj = 0;
    for e = 1:numel(edge_faces)
        faces = edge_faces{e};
        dirs = edge_dirs{e};
        if numel(faces) ~= 2
            continue
        end
        relation = -dirs(1) * dirs(2);
        n_adj = n_adj + 1;
        adj_i(n_adj) = faces(1);
        adj_j(n_adj) = faces(2);
        adj_relation(n_adj) = relation;
        n_adj = n_adj + 1;
        adj_i(n_adj) = faces(2);
        adj_j(n_adj) = faces(1);
        adj_relation(n_adj) = relation;
    end
    adj_i = adj_i(1:n_adj);
    adj_j = adj_j(1:n_adj);
    adj_relation = adj_relation(1:n_adj);

    incident = cell(n_f, 1);
    for k = 1:n_adj
        incident{adj_i(k)}(end + 1) = k;
    end

    face_sign = zeros(n_f, 1);
    for start = 1:n_f
        if face_sign(start) ~= 0
            continue
        end
        face_sign(start) = 1;
        queue = zeros(n_f, 1);
        queue(1) = start;
        head = 1;
        tail = 1;
        while head <= tail
            f = queue(head);
            head = head + 1;
            edges = incident{f};
            for k = edges
                g = adj_j(k);
                required = adj_relation(k) * face_sign(f);
                if face_sign(g) == 0
                    face_sign(g) = required;
                    tail = tail + 1;
                    queue(tail) = g;
                end
            end
        end
    end

    flip_faces = face_sign < 0;
    M(flip_faces, :) = M(flip_faces, [1, 3, 2]);
end

function c = triangle_cross(P, T)
    if isempty(T)
        c = zeros(0, 3);
        return
    end
    c = cross(P(T(:, 2), :) - P(T(:, 1), :), ...
        P(T(:, 3), :) - P(T(:, 1), :), 2);
end

function q = triangle_quality(P, T)
    if isempty(T)
        q = zeros(0, 1);
        return
    end
    e1 = P(T(:, 2), :) - P(T(:, 1), :);
    e2 = P(T(:, 3), :) - P(T(:, 2), :);
    e3 = P(T(:, 1), :) - P(T(:, 3), :);
    a2 = sum(e1 .^ 2, 2);
    b2 = sum(e2 .^ 2, 2);
    c2 = sum(e3 .^ 2, 2);
    area2 = vecnorm(cross(e1, -e3, 2), 2, 2);
    q = 2 * sqrt(3) * area2 ./ max(a2 + b2 + c2, eps);
end

function [M, P, curv_v, target_h] = compact_vertices(M, P, curv_v, target_h)
    used = unique(M(:));
    map = zeros(size(P, 1), 1);
    map(used) = 1:numel(used);
    M = map(M);
    P = P(used, :);
    curv_v = curv_v(used);
    target_h = target_h(used);
end

function h = median_edge_length(M, P)
    directed_edges = [M(:, [1, 2]); M(:, [2, 3]); M(:, [3, 1])];
    edges = unique(sort(directed_edges, 2), "rows");
    h = median(vecnorm(P(edges(:, 1), :) - P(edges(:, 2), :), 2, 2));
end

function n_bad = count_bad_edges(M)
    edges = [M(:, [1, 2]); M(:, [2, 3]); M(:, [3, 1])];
    [~, ~, edge_id] = unique(sort(edges, 2), "rows");
    edge_count = accumarray(edge_id, 1);
    n_bad = nnz(edge_count ~= 2);
end

function q = local_percentile(x, percent)
    x = sort(x(isfinite(x)));
    if isempty(x)
        q = 0;
        return
    end
    percent = min(max(percent, 0), 100);
    t = 1 + (numel(x) - 1) * percent / 100;
    lo = floor(t);
    hi = ceil(t);
    if lo == hi
        q = x(lo);
    else
        w = t - lo;
        q = (1 - w) * x(lo) + w * x(hi);
    end
end
