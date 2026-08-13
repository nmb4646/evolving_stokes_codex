function info = mesh_surface_diagnostics(M, P, reference)
%MESH_SURFACE_DIAGNOSTICS Geometry and quality metrics for a triangle mesh.
%
%   info = mesh_surface_diagnostics(M, P)
%   info = mesh_surface_diagnostics(M, P, struct("M", M0, "P", P0))
%
% The optional reference adds corresponding-vertex displacement when the
% connectivity sizes agree and symmetric nearest-vertex distances otherwise.

    if nargin < 3
        reference = [];
    end

    M = double(M);
    P = double(P);
    info = empty_info();
    info.n_vertices = size(P, 1);
    info.n_faces = size(M, 1);
    info.finite_vertices = size(P, 2) == 3 && all(isfinite(P), "all");
    info.integer_faces = size(M, 2) == 3 && all(M == round(M), "all");
    info.indices_in_range = info.integer_faces && ~isempty(M) ...
        && min(M, [], "all") >= 1 && max(M, [], "all") <= size(P, 1);

    if ~info.finite_vertices || ~info.indices_in_range
        info.valid = false;
        info.quality_flags = "INVALID_INPUT";
        return
    end

    p1 = P(M(:, 1), :);
    p2 = P(M(:, 2), :);
    p3 = P(M(:, 3), :);
    e12 = p2 - p1;
    e23 = p3 - p2;
    e31 = p1 - p3;
    l12 = vecnorm(e12, 2, 2);
    l23 = vecnorm(e23, 2, 2);
    l31 = vecnorm(e31, 2, 2);
    twice_area_vector = cross(e12, p3 - p1, 2);
    twice_area = vecnorm(twice_area_vector, 2, 2);
    face_area = 0.5 * twice_area;

    scale = max(max(P, [], 1) - min(P, [], 1));
    area_tolerance = max(eps(max(scale, 1)^2) * 100, realmin);
    info.n_degenerate_faces = nnz(face_area <= area_tolerance);
    info.min_face_area = min(face_area);
    info.max_face_area = max(face_area);
    info.area = sum(face_area);
    info.signed_volume = sum(dot(p1, cross(p2, p3, 2), 2)) / 6;
    info.volume = abs(info.signed_volume);

    directed_edges = [M(:, [1, 2]); M(:, [2, 3]); M(:, [3, 1])];
    sorted_edges = sort(directed_edges, 2);
    [edges, ~, edge_group] = unique(sorted_edges, "rows");
    incidence = accumarray(edge_group, 1);
    direction = 2 * (directed_edges(:, 1) < directed_edges(:, 2)) - 1;
    direction_balance = accumarray(edge_group, direction);
    info.n_edges = size(edges, 1);
    info.n_boundary_edges = nnz(incidence == 1);
    info.n_nonmanifold_edges = nnz(incidence > 2);
    info.closed = all(incidence == 2);
    info.orientation_consistent = info.closed && all(direction_balance == 0);
    info.euler_characteristic = info.n_vertices - info.n_edges + info.n_faces;

    edge_vector = P(edges(:, 2), :) - P(edges(:, 1), :);
    edge_length = vecnorm(edge_vector, 2, 2);
    info.min_edge_length = min(edge_length);
    info.max_edge_length = max(edge_length);
    info.mean_edge_length = mean(edge_length);
    info.median_edge_length = median(edge_length);
    info.edge_length_cv = std(edge_length) / max(info.mean_edge_length, eps);

    angle1 = triangle_angle(l12, l31, l23);
    angle2 = triangle_angle(l12, l23, l31);
    angle3 = triangle_angle(l23, l31, l12);
    all_angles = [angle1; angle2; angle3];
    info.min_triangle_angle_deg = min(all_angles) * 180 / pi;
    info.max_triangle_angle_deg = max(all_angles) * 180 / pi;

    perimeter = l12 + l23 + l31;
    longest = max([l12, l23, l31], [], 2);
    aspect = longest .* perimeter ./ max(4 * sqrt(3) * face_area, realmin);
    info.max_aspect_ratio = max(aspect);
    info.mean_aspect_ratio = mean(aspect);

    vertex_area = accumarray(M(:), repmat(face_area / 3, 3, 1), ...
        [size(P, 1), 1], @sum, 0);
    info.min_vertex_area = min(vertex_area);
    info.max_vertex_area = max(vertex_area);
    info.vertex_area_ratio = info.max_vertex_area / max(info.min_vertex_area, realmin);

    flags = strings(0, 1);
    if ~info.closed
        flags(end + 1) = "NOT_CLOSED";
    end
    if ~info.orientation_consistent
        flags(end + 1) = "INCONSISTENT_ORIENTATION";
    end
    if info.n_degenerate_faces > 0
        flags(end + 1) = "DEGENERATE_FACES";
    end
    info.quality_flags = strjoin(flags, ";");
    info.valid = info.closed && info.orientation_consistent ...
        && info.n_degenerate_faces == 0;

    if ~isempty(reference)
        info = add_reference_distances(info, M, P, reference, vertex_area);
    end
end

function angle = triangle_angle(side_a, side_b, opposite)
    denominator = 2 * side_a .* side_b;
    cosine = (side_a.^2 + side_b.^2 - opposite.^2) ...
        ./ max(denominator, realmin);
    cosine = min(max(cosine, -1), 1);
    angle = acos(cosine);
end

function info = add_reference_distances(info, M, P, reference, vertex_area)
    if ~isstruct(reference) || ~isfield(reference, "P")
        error("reference must be a struct containing P and optionally M.");
    end
    P0 = double(reference.P);
    if size(P0, 2) ~= 3 || any(~isfinite(P0), "all")
        error("reference.P must be a finite N-by-3 array.");
    end

    if size(P0, 1) == size(P, 1) ...
            && (~isfield(reference, "M") || isequal(double(reference.M), M))
        displacement = P - P0;
        displacement_norm = vecnorm(displacement, 2, 2);
        info.max_corresponding_displacement = max(displacement_norm);
        info.rms_corresponding_displacement = sqrt(mean(displacement_norm.^2));
        info.area_weighted_rms_displacement = sqrt( ...
            sum(vertex_area .* displacement_norm.^2) / max(sum(vertex_area), eps));
        info.max_displacement_over_mean_edge = ...
            info.max_corresponding_displacement / max(info.mean_edge_length, eps);
    end

    distance_to_reference = nearest_vertex_distance(P, P0);
    distance_from_reference = nearest_vertex_distance(P0, P);
    info.approx_symmetric_hausdorff = max( ...
        max(distance_to_reference), max(distance_from_reference));
    info.approx_symmetric_rms_distance = sqrt(0.5 * ( ...
        mean(distance_to_reference.^2) + mean(distance_from_reference.^2)));
end

function distance = nearest_vertex_distance(query, reference)
    if exist("knnsearch", "file") ~= 0
        [~, distance] = knnsearch(reference, query, "K", 1);
        return
    end
    distance = zeros(size(query, 1), 1);
    chunk_size = 250;
    for first = 1:chunk_size:size(query, 1)
        last = min(first + chunk_size - 1, size(query, 1));
        difference = permute(query(first:last, :), [1, 3, 2]) ...
            - permute(reference, [3, 1, 2]);
        distance_squared = sum(difference.^2, 3);
        distance(first:last) = sqrt(min(distance_squared, [], 2));
    end
end

function info = empty_info()
    info = struct( ...
        "valid", false, "quality_flags", "", ...
        "n_vertices", NaN, "n_faces", NaN, "n_edges", NaN, ...
        "finite_vertices", false, "integer_faces", false, ...
        "indices_in_range", false, "closed", false, ...
        "orientation_consistent", false, "n_boundary_edges", NaN, ...
        "n_nonmanifold_edges", NaN, "n_degenerate_faces", NaN, ...
        "euler_characteristic", NaN, "area", NaN, ...
        "signed_volume", NaN, "volume", NaN, ...
        "min_face_area", NaN, "max_face_area", NaN, ...
        "min_edge_length", NaN, "max_edge_length", NaN, ...
        "mean_edge_length", NaN, "median_edge_length", NaN, ...
        "edge_length_cv", NaN, "min_triangle_angle_deg", NaN, ...
        "max_triangle_angle_deg", NaN, "max_aspect_ratio", NaN, ...
        "mean_aspect_ratio", NaN, "min_vertex_area", NaN, ...
        "max_vertex_area", NaN, "vertex_area_ratio", NaN, ...
        "max_corresponding_displacement", NaN, ...
        "rms_corresponding_displacement", NaN, ...
        "area_weighted_rms_displacement", NaN, ...
        "max_displacement_over_mean_edge", NaN, ...
        "approx_symmetric_hausdorff", NaN, ...
        "approx_symmetric_rms_distance", NaN);
end
