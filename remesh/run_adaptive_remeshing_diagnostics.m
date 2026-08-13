function results = run_adaptive_remeshing_diagnostics(output_dir)
%RUN_ADAPTIVE_REMESHING_DIAGNOSTICS Compare legacy and MMGS remeshing.
%
%   results = run_adaptive_remeshing_diagnostics()
%
% Each row separates the raw remesher error from the subsequent
% mass-weighted simultaneous area/volume correction. Results are saved as
% CSV and MAT files under data/adaptive_remeshing_diagnostics by default.

    root = fileparts(mfilename("fullpath"));
    if nargin < 1 || isempty(output_dir)
        output_dir = fullfile(root, "data", "adaptive_remeshing_diagnostics");
    end
    if ~isfolder(output_dir)
        mkdir(output_dir);
    end

    cases = benchmark_cases(root);
    backends = ["legacy", "mmgs"];
    rows = repmat(empty_row(), 0, 1);
    remeshed_surfaces = struct();

    for case_index = 1:numel(cases)
        test_case = cases(case_index);
        geo_before = Geometry(test_case.M, test_case.P);
        size_options = test_case.size_options;
        [target_h, target_info] = compute_target_edge_length( ...
            geo_before, test_case.base_h, size_options);

        for backend = backends
            options = size_options;
            options.backend = backend;
            options.legacy_iterations = 20;
            options.match_target_complexity = true;
            try
                [P_raw, M_raw, remesh_info] = remesh_surface( ...
                    test_case.P, test_case.M, target_h, options);
            catch exception
                warning("Benchmark %s/%s failed: %s", ...
                    test_case.name, backend, exception.message);
                row = empty_row();
                row.case_name = test_case.name;
                row.source = test_case.source;
                row.backend = backend;
                row.success = false;
                row.failure = string(exception.message);
                rows(end + 1, 1) = row; %#ok<AGROW>
                continue
            end

            geo_raw = Geometry(M_raw, P_raw);
            targets = struct("area", geo_before.area, ...
                "volume", geo_before.volume);
            projection_options = struct("tol_area", 1e-10, ...
                "tol_volume", 1e-10, "max_iter", 20, ...
                "weighting", "mass");
            [P_corrected, geo_corrected, projection_info] = ...
                project_surface_constraints(P_raw, M_raw, targets, ...
                projection_options);
            corrected_diagnostics = mesh_surface_diagnostics(M_raw, P_corrected);
            density = adaptive_density_diagnostic( ...
                test_case.P, target_h, P_raw, M_raw);
            raw_curvature = geo_raw.v_mean_curvature ./ geo_raw.v_area;
            corrected_curvature = geo_corrected.v_mean_curvature ...
                ./ geo_corrected.v_area;

            row = empty_row();
            row.case_name = test_case.name;
            row.source = test_case.source;
            row.backend = backend;
            row.success = true;
            row.adaptive = target_info.adaptive;
            row.n_vertices_before = size(test_case.P, 1);
            row.n_faces_before = size(test_case.M, 1);
            row.n_vertices_target = predicted_vertex_count( ...
                test_case.M, test_case.P, target_h);
            row.n_vertices_raw = size(P_raw, 1);
            row.n_faces_raw = size(M_raw, 1);
            row.vertex_count_relative_error = abs(row.n_vertices_raw ...
                - row.n_vertices_target) / row.n_vertices_target;
            row.runtime_seconds = remesh_info.runtime_seconds;
            row.rel_area_error_raw = remesh_info.rel_area_error_raw;
            row.rel_volume_error_raw = remesh_info.rel_volume_error_raw;
            row.rel_area_error_corrected = projection_info.rel_area_error;
            row.rel_volume_error_corrected = projection_info.rel_volume_error;
            row.min_angle_before_deg = remesh_info.before.min_triangle_angle_deg;
            row.min_angle_raw_deg = remesh_info.raw_after.min_triangle_angle_deg;
            row.min_angle_corrected_deg = corrected_diagnostics.min_triangle_angle_deg;
            row.max_aspect_before = remesh_info.before.max_aspect_ratio;
            row.max_aspect_raw = remesh_info.raw_after.max_aspect_ratio;
            row.max_aspect_corrected = corrected_diagnostics.max_aspect_ratio;
            row.min_edge_raw = remesh_info.raw_after.min_edge_length;
            row.mean_edge_raw = remesh_info.raw_after.mean_edge_length;
            row.max_edge_raw = remesh_info.raw_after.max_edge_length;
            row.approx_surface_hausdorff = ...
                remesh_info.raw_after.approx_symmetric_hausdorff;
            row.approx_surface_rms = ...
                remesh_info.raw_after.approx_symmetric_rms_distance;
            row.projection_converged = projection_info.converged;
            row.projection_iterations = projection_info.iterations;
            row.projection_max_displacement = projection_info.max_displacement;
            row.projection_rms_displacement = projection_info.rms_displacement;
            row.projection_max_displacement_over_edge = ...
                projection_info.max_displacement / max(row.mean_edge_raw, eps);
            row.projection_curvature_relative_rms_change = rms_value( ...
                corrected_curvature - raw_curvature) ...
                / max(rms_value(raw_curvature), eps);
            row.density_edge_target_correlation = density.correlation;
            row.high_to_low_resolution_edge_ratio = density.high_to_low_ratio;
            row.valid_raw_mesh = remesh_info.raw_after.valid;
            row.valid_corrected_mesh = corrected_diagnostics.valid;
            if backend == "mmgs"
                row.mmgs_passes = remesh_info.status.complexity_passes;
                row.mmgs_size_scale = remesh_info.status.mmgs_size_scale;
            end
            rows(end + 1, 1) = row; %#ok<AGROW>

            key = matlab.lang.makeValidName(test_case.name + "_" + backend);
            remeshed_surfaces.(key) = struct( ...
                "M_before", test_case.M, "P_before", test_case.P, ...
                "target_h", target_h, "M_raw", M_raw, "P_raw", P_raw, ...
                "P_corrected", P_corrected);
        end
    end

    results = struct2table(rows);
    writetable(results, fullfile(output_dir, "remeshing_backend_comparison.csv"));
    save(fullfile(output_dir, "remeshing_backend_comparison.mat"), ...
        "results", "remeshed_surfaces");

    fprintf("\nAdaptive remeshing benchmark\n");
    disp(results(:, ["case_name", "backend", "n_vertices_before", ...
        "n_vertices_target", "n_vertices_raw", "rel_area_error_raw", ...
        "rel_volume_error_raw", "rel_area_error_corrected", ...
        "rel_volume_error_corrected", "min_angle_raw_deg", ...
        "projection_max_displacement_over_edge"]));
end

function cases = benchmark_cases(root)
    cases = repmat(make_case("", "", zeros(0, 3), zeros(0, 3), ...
        NaN, struct()), 0, 1);
    sphere_resolutions = [3, 5, 8];
    for resolution = sphere_resolutions
        [P, M] = subdivided_sphere(resolution);
        P = P ./ vecnorm(P, 2, 2);
        geo = Geometry(M, P);
        cases(end + 1) = make_case( ...
            "sphere_N" + string(resolution), "synthetic", M, P, ...
            mean(geo.he_length), struct("adaptive", false)); %#ok<AGROW>
    end

    [P, M] = subdivided_sphere(7);
    P = P ./ vecnorm(P, 2, 2) .* [1.6, 0.8, 0.8];
    geo = Geometry(M, P);
    cases(end + 1) = make_case("spheroid", "synthetic", M, P, ...
        mean(geo.he_length), default_adaptive_options());

    [M, P] = initial_dumbbell(0.35, [24, 48]);
    geo = Geometry(M, P);
    cases(end + 1) = make_case("tubular_dumbbell", "synthetic", M, P, ...
        mean(geo.he_length), default_adaptive_options());

    saved = [ ...
        struct("name", "saved_thin_tubule", "relative_path", fullfile( ...
            "data", "fs_batch_data", ...
            "Sd_1.00em06_Da_0.00ep00_gamy_p3.70em06_v_3.50em01", ...
            "geo132.mat")), ...
        struct("name", "saved_pearled_membrane", "relative_path", fullfile( ...
            "data", "fs_batch_data", ...
            "Sd_1.00em06_Da_0.00ep00_gamy_p8.00em07_v_8.63em01", ...
            "geo200.mat"))];
    for index = 1:numel(saved)
        path = fullfile(root, saved(index).relative_path);
        if ~isfile(path)
            warning("Saved benchmark frame not found: %s", path);
            continue
        end
        data = load(path, "M", "P");
        geo = Geometry(data.M, data.P);
        cases(end + 1) = make_case(saved(index).name, ...
            string(saved(index).relative_path), data.M, data.P, ...
            mean(geo.he_length), default_adaptive_options()); %#ok<AGROW>
    end
end

function options = default_adaptive_options()
    options = struct( ...
        "adaptive", true, ...
        "curvature_measure", "max_abs_principal", ...
        "curvature_weight", 1.0, ...
        "curvature_power", 1.0, ...
        "hmin_factor", 0.35, ...
        "hmax_factor", 2.0, ...
        "preserve_vertex_budget", true);
end

function test_case = make_case(name, source, M, P, base_h, size_options)
    test_case = struct("name", string(name), "source", string(source), ...
        "M", double(M), "P", double(P), "base_h", base_h, ...
        "size_options", size_options);
end

function density = adaptive_density_diagnostic(P_old, h_old, P_new, M_new)
    index = nearest_vertex_index(P_new, P_old);
    h_new = h_old(index);
    edges = unique(sort([M_new(:, [1, 2]); M_new(:, [2, 3]); ...
        M_new(:, [3, 1])], 2), "rows");
    edge_length = vecnorm(P_new(edges(:, 2), :) - P_new(edges(:, 1), :), 2, 2);
    edge_h = 0.5 * (h_new(edges(:, 1)) + h_new(edges(:, 2)));
    density.correlation = correlation_coefficient(log(edge_h), log(edge_length));
    lower = local_percentile(edge_h, 25);
    upper = local_percentile(edge_h, 75);
    high_resolution = edge_h <= lower;
    low_resolution = edge_h >= upper;
    density.high_to_low_ratio = median(edge_length(high_resolution)) ...
        / median(edge_length(low_resolution));
end

function count = predicted_vertex_count(M, P, h)
    p1 = P(M(:, 1), :);
    p2 = P(M(:, 2), :);
    p3 = P(M(:, 3), :);
    face_area = 0.5 * vecnorm(cross(p2 - p1, p3 - p1, 2), 2, 2);
    vertex_area = accumarray(M(:), repmat(face_area / 3, 3, 1), ...
        [size(P, 1), 1], @sum, 0);
    count = max(4, round(2 / sqrt(3) * sum(vertex_area ./ h.^2) + 2));
end

function index = nearest_vertex_index(query, reference)
    if exist("knnsearch", "file") ~= 0
        index = knnsearch(reference, query, "K", 1);
        return
    end
    index = zeros(size(query, 1), 1);
    chunk_size = 250;
    for first = 1:chunk_size:size(query, 1)
        last = min(first + chunk_size - 1, size(query, 1));
        difference = permute(query(first:last, :), [1, 3, 2]) ...
            - permute(reference, [3, 1, 2]);
        distance_squared = sum(difference.^2, 3);
        [~, index(first:last)] = min(distance_squared, [], 2);
    end
end

function value = correlation_coefficient(x, y)
    x = x(:) - mean(x);
    y = y(:) - mean(y);
    value = dot(x, y) / max(norm(x) * norm(y), realmin);
end

function value = local_percentile(values, percentile)
    values = sort(values(:));
    position = 1 + (numel(values) - 1) * percentile / 100;
    lower = floor(position);
    upper = ceil(position);
    fraction = position - lower;
    value = (1 - fraction) * values(lower) + fraction * values(upper);
end

function value = rms_value(values)
    value = sqrt(mean(values(:).^2));
end

function row = empty_row()
    row = struct( ...
        "case_name", "", "source", "", "backend", "", ...
        "success", false, "failure", "", "adaptive", false, ...
        "n_vertices_before", NaN, "n_faces_before", NaN, ...
        "n_vertices_target", NaN, "n_vertices_raw", NaN, ...
        "n_faces_raw", NaN, "vertex_count_relative_error", NaN, ...
        "runtime_seconds", NaN, "rel_area_error_raw", NaN, ...
        "rel_volume_error_raw", NaN, "rel_area_error_corrected", NaN, ...
        "rel_volume_error_corrected", NaN, "min_angle_before_deg", NaN, ...
        "min_angle_raw_deg", NaN, "min_angle_corrected_deg", NaN, ...
        "max_aspect_before", NaN, "max_aspect_raw", NaN, ...
        "max_aspect_corrected", NaN, "min_edge_raw", NaN, ...
        "mean_edge_raw", NaN, "max_edge_raw", NaN, ...
        "approx_surface_hausdorff", NaN, "approx_surface_rms", NaN, ...
        "projection_converged", false, "projection_iterations", NaN, ...
        "projection_max_displacement", NaN, ...
        "projection_rms_displacement", NaN, ...
        "projection_max_displacement_over_edge", NaN, ...
        "projection_curvature_relative_rms_change", NaN, ...
        "density_edge_target_correlation", NaN, ...
        "high_to_low_resolution_edge_ratio", NaN, ...
        "valid_raw_mesh", false, "valid_corrected_mesh", false, ...
        "mmgs_passes", NaN, "mmgs_size_scale", NaN);
end
