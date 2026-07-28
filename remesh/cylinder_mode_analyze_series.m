function result = cylinder_mode_analyze_series(series_folder, cfg, shared_basis)
%CYLINDER_MODE_ANALYZE_SERIES Analyze cylindrical modes in one geo*.mat series.

if nargin < 3
    shared_basis = [];
end

series_folder = char(series_folder);
if ~isfolder(series_folder)
    error("CylinderMode:MissingSeries", "Series folder not found: %s", series_folder);
end

[frame_files, frame_index] = discover_frames(series_folder, cfg);
n_frames = numel(frame_files);
if n_frames < 2
    error("CylinderMode:TooFewFrames", ...
        "Need at least two selected geo*.mat files in %s.", series_folder);
end

series_id = string(get_folder_name(series_folder));
[mode_m, mode_n] = signed_mode_list(cfg);
[group_m, group_n, group_raw_primary, group_raw_conjugate] = ...
    grouped_mode_list(mode_m, mode_n, cfg);
n_modes = numel(mode_m);
n_groups = numel(group_m);

coefficient = complex(nan(n_frames, n_modes));
group_amplitude = nan(n_frames, n_groups);
group_phase = nan(n_frames, n_groups);
group_noise = nan(n_frames, n_groups);
group_qR = nan(n_frames, n_groups);
group_resolved = false(n_frames, n_groups);

metric = initialize_metric_arrays(n_frames);
profiles = cell(3, 1);
profile_positions = unique([1, round((n_frames + 1) / 2), n_frames]);
projection_states = cell(n_frames, 1);

initial_connectivity = [];
persistent_core_ids = [];
previous_axis = [];
previous_transverse = [];
initial_area = NaN;
initial_volume = NaN;
initial_core_axial_fraction = NaN;
previous_area = NaN;
previous_volume = NaN;

for k = 1:n_frames
    data = load_frame(frame_files(k), cfg);
    validate_loaded_frame(data, frame_files(k));
    P = double(data.P);
    M = double(data.M);

    geometry = basic_geometry(P, M);
    frame_flags = geometry.flags;

    [axis_vector, transverse_vector, center, axis_jump] = align_frame( ...
        P, geometry.vertex_area, previous_axis, previous_transverse, ...
        shared_basis, cfg);
    binormal_vector = cross(axis_vector, transverse_vector);
    binormal_vector = binormal_vector / norm(binormal_vector);
    basis = [axis_vector, transverse_vector, binormal_vector];
    local_P = (P - center) * basis;
    x = local_P(:, 1);
    r = hypot(local_P(:, 2), local_P(:, 3));
    phi = atan2(local_P(:, 3), local_P(:, 2));

    same_connectivity = k == 1 || ...
        (size(M, 1) == size(initial_connectivity, 1) && isequal(M, initial_connectivity));
    use_persistent = same_connectivity && ~isempty(persistent_core_ids) ...
        && (cfg.core.method == "persistent" || cfg.core.method == "auto");

    if k == 1
        [core_mask, core_info, core_flags] = detect_core( ...
            x, r, geometry.vertex_area, cfg);
        initial_connectivity = M;
        persistent_core_ids = find(core_mask);
    elseif use_persistent && all(persistent_core_ids <= size(P, 1))
        core_mask = false(size(P, 1), 1);
        core_mask(persistent_core_ids) = true;
        core_info = core_info_from_mask(x, r, geometry.vertex_area, core_mask, cfg);
        core_flags = strings(0, 1);
    else
        [core_mask, core_info, core_flags] = detect_core( ...
            x, r, geometry.vertex_area, cfg);
        if ~same_connectivity
            core_flags(end + 1) = "TOPOLOGY_CHANGED";
        end
    end
    frame_flags = [frame_flags; core_flags(:)]; %#ok<AGROW>

    [base, W, eta, base_flags] = cylindrical_base( ...
        x, r, geometry.vertex_area, core_mask, core_info, cfg);
    frame_flags = [frame_flags; base_flags(:)]; %#ok<AGROW>
    full_x_min = weighted_quantile(x, geometry.vertex_area, ...
        cfg.core.robust_endpoint_quantile);
    full_x_max = weighted_quantile(x, geometry.vertex_area, ...
        1 - cfg.core.robust_endpoint_quantile);
    core_axial_fraction = base.length / (full_x_max - full_x_min);
    if k == 1
        initial_core_axial_fraction = core_axial_fraction;
    elseif isfinite(core_axial_fraction) && isfinite(initial_core_axial_fraction) ...
            && abs(core_axial_fraction / initial_core_axial_fraction - 1) ...
            > cfg.core.maximum_relative_axial_fraction_change
        frame_flags(end + 1) = "CORE_DOMAIN_DRIFT";
    end

    if base.valid
        [coeff, eta_reconstructed, projection_condition] = project_modes( ...
            phi(core_mask), base.xi, eta, W, mode_m, mode_n, cfg);
        coefficient(k, :) = coeff.';

        residual = eta - eta_reconstructed;
        reconstruction_error = weighted_relative_norm(residual, eta, W);
        effective_vertices = sum(W) ^ 2 / max(sum(W .^ 2), eps);
        residual_rms = sqrt(sum(W .* abs(residual) .^ 2) / max(sum(W), eps));
        coefficient_noise = max(residual_rms / sqrt(max(effective_vertices, 1)), ...
            cfg.growth.absolute_noise_floor);

        if cfg.pairing.enabled
            projection_states{k} = struct( ...
                "vertex_ids", find(core_mask), ...
                "x", x(core_mask), ...
                "r", r(core_mask), ...
                "phi", phi(core_mask), ...
                "vertex_area", geometry.vertex_area(core_mask), ...
                "x_min", base.x_min, ...
                "x_max", base.x_max, ...
                "radius", base.radius, ...
                "mean_edge_length", geometry.mean_edge_length);
        end

        for g = 1:n_groups
            primary = coeff(group_raw_primary(g));
            conjugate_index = group_raw_conjugate(g);
            if conjugate_index == group_raw_primary(g)
                amp = abs(primary);
                pair_count = 1;
            else
                amp = hypot(abs(primary), abs(coeff(conjugate_index)));
                pair_count = 2;
            end
            group_amplitude(k, g) = amp;
            group_phase(k, g) = angle(primary);
            group_noise(k, g) = sqrt(pair_count) * coefficient_noise;
            group_qR(k, g) = 2 * pi * abs(group_n(g)) * base.radius / base.length;
            group_resolved(k, g) = mode_is_resolved( ...
                group_m(g), group_n(g), base.radius, base.length, ...
                geometry.mean_edge_length, cfg);
        end

        conjugacy_error = coefficient_conjugacy_error( ...
            coeff, mode_m, mode_n);
        if projection_condition > cfg.modes.maximum_condition_number
            frame_flags(end + 1) = "ILL_CONDITIONED_PROJECTION";
        end
        if reconstruction_error > cfg.quality.maximum_reconstruction_error
            frame_flags(end + 1) = "LARGE_RECONSTRUCTION_ERROR";
        end
        if effective_vertices < cfg.quality.minimum_effective_vertices
            frame_flags(end + 1) = "TOO_FEW_EFFECTIVE_VERTICES";
        end
    else
        projection_condition = NaN;
        reconstruction_error = NaN;
        effective_vertices = NaN;
        conjugacy_error = NaN;
    end

    if k == 1
        initial_area = geometry.area;
        initial_volume = geometry.volume;
    else
        if abs(geometry.area / previous_area - 1) > cfg.quality.maximum_area_relative_jump
            frame_flags(end + 1) = "AREA_JUMP";
        end
        if previous_volume > 0 && ...
                abs(geometry.volume / previous_volume - 1) > cfg.quality.maximum_volume_relative_jump
            frame_flags(end + 1) = "VOLUME_JUMP";
        end
        if abs(geometry.area / initial_area - 1) > cfg.quality.maximum_area_relative_drift
            frame_flags(end + 1) = "AREA_DRIFT";
        end
        if initial_volume > 0 && ...
                abs(geometry.volume / initial_volume - 1) > cfg.quality.maximum_volume_relative_drift
            frame_flags(end + 1) = "VOLUME_DRIFT";
        end
    end
    if axis_jump > cfg.alignment.maximum_axis_jump_degrees
        frame_flags(end + 1) = "AXIS_JUMP";
    end

    frame_flags = unique(frame_flags(strlength(frame_flags) > 0), "stable");
    hard_failure = any(ismember(frame_flags, [ ...
        "INVALID_GEOMETRY", "CORE_NOT_FOUND", "CORE_TOO_SHORT", ...
        "CORE_DOMAIN_DRIFT", "AREA_JUMP", "VOLUME_JUMP", ...
        "ILL_CONDITIONED_PROJECTION"]));

    metric.explicit_time(k) = extract_explicit_time(data, cfg);
    metric.dt(k) = extract_dt(data);
    metric.n_vertices(k) = size(P, 1);
    metric.n_faces(k) = size(M, 1);
    metric.surface_area(k) = geometry.area;
    metric.volume(k) = geometry.volume;
    metric.reduced_volume(k) = geometry.reduced_volume;
    metric.center(k, :) = center;
    metric.axis(k, :) = axis_vector;
    metric.axis_angle_change(k) = axis_jump;
    metric.core_x_min(k) = base.x_min;
    metric.core_x_max(k) = base.x_max;
    metric.core_length(k) = base.length;
    metric.mean_core_radius(k) = base.radius;
    metric.radius_std(k) = base.radius_std;
    metric.cylindricality_error(k) = base.cylindricality_error;
    metric.core_area_fraction(k) = sum(geometry.vertex_area(core_mask)) / geometry.area;
    metric.core_axial_fraction(k) = core_axial_fraction;
    metric.windowed_area(k) = sum(W);
    metric.effective_vertices(k) = effective_vertices;
    metric.reconstruction_error(k) = reconstruction_error;
    metric.projection_condition(k) = projection_condition;
    metric.conjugacy_error(k) = conjugacy_error;
    metric.mean_edge_length(k) = geometry.mean_edge_length;
    metric.frame_valid(k) = ~hard_failure && base.valid;
    metric.quality_flags(k) = join_flags(frame_flags);

    profile_slot = find(profile_positions == k, 1);
    if ~isempty(profile_slot)
        profiles{profile_slot} = struct( ...
            "frame_index", frame_index(k), "x", x, "r", r, ...
            "core_mask", core_mask, "x_min", base.x_min, "x_max", base.x_max);
    end

    previous_axis = axis_vector;
    previous_transverse = transverse_vector;
    previous_area = geometry.area;
    previous_volume = geometry.volume;

    if cfg.verbose && (k == 1 || k == n_frames || mod(k, max(1, floor(n_frames / 10))) == 0)
        fprintf("  %s: frame %d (%d/%d), R = %.4g, L = %.4g, Erec = %.3g\n", ...
            series_id, frame_index(k), k, n_frames, base.radius, base.length, ...
            reconstruction_error);
    end
end

time = resolve_physical_times(frame_index, metric.explicit_time, metric.dt);
frame_metrics = make_frame_metrics_table( ...
    series_id, frame_index, time, metric);

group_noise_for_fit = group_noise;
group_noise_for_fit(~group_resolved) = inf;
[growth_rates, fit_details] = cylinder_mode_fit_growth( ...
    time, frame_index, group_amplitude, group_phase, group_noise_for_fit, ...
    group_qR, group_m, group_n, metric.frame_valid, cfg);

growth_rates.series_id = repmat(series_id, height(growth_rates), 1);
growth_rates.mean_radius = nan(height(growth_rates), 1);
growth_rates.mean_core_length = nan(height(growth_rates), 1);
growth_rates.resolved_fraction = mean(group_resolved, 1).';
for g = 1:n_groups
    idx = fit_details(g).indices;
    if ~isempty(idx)
        growth_rates.mean_radius(g) = mean(metric.mean_core_radius(idx), "omitnan");
        growth_rates.mean_core_length(g) = mean(metric.core_length(idx), "omitnan");
    end
    if ~any(group_resolved(:, g))
        growth_rates.quality_flags(g) = append_flag( ...
            growth_rates.quality_flags(g), "MODE_UNRESOLVED");
        growth_rates.fit_status(g) = "no_fit";
    end
    if group_m(g) == 0 && group_n(g) == 0
        growth_rates.fit_status(g) = "excluded_base_mode";
        growth_rates.quality_flags(g) = append_flag( ...
            growth_rates.quality_flags(g), "BASE_MODE");
        growth_rates.growth_rate_sigma(g) = NaN;
        growth_rates.growth_rate_standard_error(g) = NaN;
        growth_rates.phase_rate_omega(g) = NaN;
    end
end
growth_rates = movevars(growth_rates, "series_id", "Before", 1);

mode_coefficients = make_mode_table( ...
    series_id, frame_index, time, coefficient, mode_m, mode_n, ...
    group_m, group_n, group_amplitude, group_noise, group_qR, ...
    metric.core_length, metric.mean_core_radius, metric.quality_flags);

result = struct();
result.series_id = series_id;
result.series_folder = string(series_folder);
result.frame_metrics = frame_metrics;
result.mode_coefficients = mode_coefficients;
result.growth_rates = growth_rates;
result.fit_details = fit_details;
result.time = time;
result.frame_index = frame_index;
result.mode_m = mode_m;
result.mode_n = mode_n;
result.coefficient = coefficient;
result.group_m = group_m;
result.group_n = group_n;
result.group_amplitude = group_amplitude;
result.group_phase = group_phase;
result.group_noise = group_noise;
result.group_qR = group_qR;
result.group_resolved = group_resolved;
result.profiles = profiles;
result.projection_states = projection_states;
result.status = "ok";
result.error_message = "";
end

function [files, indices] = discover_frames(folder, cfg)
    listing = dir(fullfile(folder, "geo*.mat"));
    names = string({listing.name}).';
    indices = nan(numel(names), 1);
    keep = false(numel(names), 1);
    for i = 1:numel(names)
        token = regexp(names(i), "^geo(\d+)\.mat$", "tokens", "once");
        if ~isempty(token)
            indices(i) = str2double(token{1});
            keep(i) = true;
        end
    end
    listing = listing(keep);
    indices = indices(keep);
    [indices, order] = sort(indices);
    listing = listing(order);

    selected = indices >= cfg.frame.first_index ...
        & indices <= cfg.frame.last_index;
    selected_positions = find(selected);
    selected_positions = selected_positions(1:cfg.frame.stride:end);
    if isfinite(cfg.frame.max_frames)
        selected_positions = selected_positions(1:min(numel(selected_positions), cfg.frame.max_frames));
    end
    listing = listing(selected_positions);
    indices = indices(selected_positions);
    files = strings(numel(listing), 1);
    for i = 1:numel(listing)
        files(i) = fullfile(listing(i).folder, listing(i).name);
    end
end

function validate_loaded_frame(data, filename)
    if ~isfield(data, "P") || ~isfield(data, "M")
        error("CylinderMode:InvalidFrame", "%s does not contain P and M.", filename);
    end
    if size(data.P, 2) ~= 3 || size(data.M, 2) ~= 3
        error("CylinderMode:InvalidFrame", "%s does not contain triangular 3D geometry.", filename);
    end
end

function data = load_frame(filename, cfg)
    available = string(who("-file", filename));
    requested = ["P", "M", "p", cfg.frame.time_field_names];
    selected = requested(ismember(requested, available));
    names = cellstr(selected);
    data = load(filename, names{:});
end

function geometry = basic_geometry(P, M)
    flags = strings(0, 1);
    if any(~isfinite(P), "all") || any(~isfinite(M), "all") ...
            || any(M(:) < 1) || any(M(:) > size(P, 1)) ...
            || any(M(:) ~= round(M(:)))
        flags(end + 1) = "INVALID_GEOMETRY";
    end

    v1 = P(M(:, 1), :);
    v2 = P(M(:, 2), :);
    v3 = P(M(:, 3), :);
    cross_product = cross(v2 - v1, v3 - v1, 2);
    face_area = 0.5 * vecnorm(cross_product, 2, 2);
    if nnz(face_area <= eps(max(face_area))) > max(1, 0.001 * numel(face_area))
        flags(end + 1) = "DEGENERATE_TRIANGLES";
    end

    vertex_area = accumarray(M(:), repmat(face_area / 3, 3, 1), ...
        [size(P, 1), 1], @sum, 0);
    area = sum(face_area);
    signed_volume = sum(dot(v1, cross(v2, v3, 2), 2)) / 6;
    volume = abs(signed_volume);
    reduced_volume = 6 * sqrt(pi) * volume / max(area ^ 1.5, eps);

    edges = [M(:, [1, 2]); M(:, [2, 3]); M(:, [3, 1])];
    edges = sort(edges, 2);
    edges = unique(edges, "rows");
    edge_length = vecnorm(P(edges(:, 1), :) - P(edges(:, 2), :), 2, 2);

    geometry = struct( ...
        "face_area", face_area, ...
        "vertex_area", vertex_area, ...
        "area", area, ...
        "volume", volume, ...
        "signed_volume", signed_volume, ...
        "reduced_volume", reduced_volume, ...
        "mean_edge_length", mean(edge_length), ...
        "flags", flags);
end

function [axis_vector, transverse, center, angle_change] = align_frame( ...
        P, vertex_area, previous_axis, previous_transverse, shared_basis, cfg)
    center = sum(P .* vertex_area, 1) / sum(vertex_area);
    centered = P - center;

    if ~isempty(shared_basis)
        initial_axis = shared_basis(:, 1);
        initial_transverse = shared_basis(:, 2);
    else
        initial_axis = [];
        initial_transverse = [];
    end

    if cfg.alignment.axis_mode == "known"
        axis_raw = cfg.alignment.known_axis(:);
        if numel(axis_raw) ~= 3 || norm(axis_raw) == 0
            error("CylinderMode:InvalidKnownAxis", ...
                "cfg.alignment.known_axis must be a nonzero three-vector.");
        end
    elseif cfg.alignment.axis_mode == "fixed_initial" && ~isempty(previous_axis)
        axis_raw = previous_axis;
    elseif isempty(previous_axis) && ~isempty(initial_axis)
        axis_raw = initial_axis;
    else
        covariance = centered.' * (centered .* vertex_area) / sum(vertex_area);
        [vectors, values] = eig((covariance + covariance.') / 2, "vector");
        [~, largest] = max(values);
        axis_raw = vectors(:, largest);
    end
    axis_raw = axis_raw / norm(axis_raw);

    if ~isempty(previous_axis)
        if dot(axis_raw, previous_axis) < 0
            axis_raw = -axis_raw;
        end
        smooth = cfg.alignment.smooth_axis_fraction;
        axis_vector = (1 - smooth) * axis_raw + smooth * previous_axis;
        axis_vector = axis_vector / norm(axis_vector);
        angle_change = acosd(max(-1, min(1, dot(axis_vector, previous_axis))));
    else
        axis_vector = axis_raw;
        angle_change = 0;
    end

    if isempty(previous_transverse)
        if ~isempty(initial_transverse)
            candidate = initial_transverse;
        else
            coordinate_axes = eye(3);
            [~, least_parallel] = min(abs(coordinate_axes.' * axis_vector));
            candidate = coordinate_axes(:, least_parallel);
        end
    else
        candidate = previous_transverse;
    end
    transverse = candidate - axis_vector * dot(axis_vector, candidate);
    if norm(transverse) < 1e-12
        coordinate_axes = eye(3);
        [~, least_parallel] = min(abs(coordinate_axes.' * axis_vector));
        candidate = coordinate_axes(:, least_parallel);
        transverse = candidate - axis_vector * dot(axis_vector, candidate);
    end
    transverse = transverse / norm(transverse);
end

function [mask, info, flags] = detect_core(x, r, vertex_area, cfg)
    flags = strings(0, 1);
    n_bins = min(cfg.core.axial_bins, max(20, floor(numel(x) / 8)));
    x_min = min(x);
    x_max = max(x);
    edges = linspace(x_min, x_max, n_bins + 1);
    centers = 0.5 * (edges(1:end-1) + edges(2:end));
    bin = discretize(x, edges);
    mean_radius = nan(1, n_bins);
    radius_std = nan(1, n_bins);
    bin_weight = zeros(1, n_bins);

    for b = 1:n_bins
        selected = bin == b;
        if nnz(selected) < 3
            continue
        end
        w = vertex_area(selected);
        values = r(selected);
        bin_weight(b) = sum(w);
        mean_radius(b) = sum(w .* values) / sum(w);
        radius_std(b) = sqrt(sum(w .* (values - mean_radius(b)) .^ 2) / sum(w));
    end

    valid_bins = isfinite(mean_radius);
    if nnz(valid_bins) < 5
        [mask, info] = fallback_core(x, r, vertex_area, cfg);
        flags(end + 1) = "CORE_FALLBACK";
        if ~core_is_usable(mask, info, cfg)
            flags(end + 1) = "CORE_NOT_FOUND";
        end
        return
    end
    mean_radius = fillmissing(mean_radius, "linear", "EndValues", "nearest");
    radius_std = fillmissing(radius_std, "linear", "EndValues", "nearest");
    smooth_radius = movmean(mean_radius, cfg.core.profile_smoothing_bins);
    smooth_std = movmean(radius_std, cfg.core.profile_smoothing_bins);
    slope = gradient(smooth_radius, centers);

    radius_reference = max(smooth_radius);
    cylindrical = smooth_radius >= cfg.core.minimum_radius_fraction * radius_reference;

    runs = contiguous_runs(cylindrical(:));
    if isempty(runs)
        [mask, info] = fallback_core(x, r, vertex_area, cfg);
        flags(end + 1) = "CORE_FALLBACK";
        if ~core_is_usable(mask, info, cfg)
            flags(end + 1) = "CORE_NOT_FOUND";
        end
        return
    end

    zero_bin = max(1, min(n_bins, discretize(0, edges)));
    contains_zero = find(runs(:, 1) <= zero_bin & runs(:, 2) >= zero_bin, 1);
    if ~isempty(contains_zero)
        chosen = runs(contains_zero, :);
    else
        lengths = runs(:, 2) - runs(:, 1) + 1;
        run_centers = 0.5 * (runs(:, 1) + runs(:, 2));
        score = lengths - 0.1 * abs(run_centers - zero_bin);
        [~, best] = max(score);
        chosen = runs(best, :);
    end

    raw_min = edges(chosen(1));
    raw_max = edges(chosen(2) + 1);
    selected_bins = chosen(1):chosen(2);
    if max(abs(slope(selected_bins))) > cfg.core.maximum_abs_radius_slope
        flags(end + 1) = "CORE_PROFILE_STEEP";
    end
    if max(smooth_std(selected_bins) ./ max(smooth_radius(selected_bins), eps)) ...
            > cfg.core.maximum_relative_cross_section_std
        flags(end + 1) = "CORE_CROSS_SECTION_VARIATION";
    end
    radius_estimate = sum(smooth_radius(selected_bins) .* bin_weight(selected_bins)) ...
        / max(sum(bin_weight(selected_bins)), eps);
    margin = cfg.core.transition_margin_in_radii * radius_estimate;
    core_min = raw_min + margin;
    core_max = raw_max - margin;
    mask = x >= core_min & x <= core_max;

    info = core_info_from_mask(x, r, vertex_area, mask, cfg);
    if ~core_is_usable(mask, info, cfg)
        [mask, info] = fallback_core(x, r, vertex_area, cfg);
        flags(end + 1) = "CORE_FALLBACK";
        if ~core_is_usable(mask, info, cfg)
            flags(end + 1) = "CORE_TOO_SHORT";
        end
    end
end

function usable = core_is_usable(mask, info, cfg)
    usable = nnz(mask) >= cfg.core.minimum_vertices ...
        && isfinite(info.length) && isfinite(info.radius) && info.radius > 0 ...
        && info.length >= cfg.core.minimum_core_length_in_radii * info.radius;
end

function [mask, info] = fallback_core(x, r, vertex_area, cfg)
    center = weighted_quantile(x, vertex_area, 0.5);
    half_length = cfg.core.fallback_half_length_fraction * (max(x) - min(x));
    mask = abs(x - center) <= half_length;
    info = core_info_from_mask(x, r, vertex_area, mask, cfg);
end

function info = core_info_from_mask(x, r, vertex_area, mask, cfg)
    if nnz(mask) < 2
        info = struct("x_min", NaN, "x_max", NaN, "length", NaN, ...
            "radius", NaN);
        return
    end
    q = cfg.core.robust_endpoint_quantile;
    x_min = weighted_quantile(x(mask), vertex_area(mask), q);
    x_max = weighted_quantile(x(mask), vertex_area(mask), 1 - q);
    radius = sum(vertex_area(mask) .* r(mask)) / sum(vertex_area(mask));
    info = struct("x_min", x_min, "x_max", x_max, ...
        "length", x_max - x_min, "radius", radius);
end

function [base, W, eta, flags] = cylindrical_base( ...
        x, r, vertex_area, core_mask, core_info, cfg)
    flags = strings(0, 1);
    base = empty_base();
    W = zeros(0, 1);
    eta = zeros(0, 1);
    if nnz(core_mask) < cfg.core.minimum_vertices || ...
            ~isfinite(core_info.length) || core_info.length <= 0
        flags(end + 1) = "CORE_NOT_FOUND";
        return
    end

    x_core = x(core_mask);
    r_core = r(core_mask);
    area_core = vertex_area(core_mask);
    x_min = core_info.x_min;
    x_max = core_info.x_max;
    length_core = x_max - x_min;
    x_center = 0.5 * (x_min + x_max);
    xi = (x_core - x_center) / length_core;
    window = axial_window(xi, cfg);
    W = area_core .* window;
    if sum(W) <= 0
        flags(end + 1) = "INVALID_WINDOW";
        return
    end

    radius_area = sum(area_core .* r_core) / sum(area_core);
    radius = sum(W .* r_core) / sum(W);
    eta = (r_core - radius) / radius;
    eta = eta - sum(W .* eta) / sum(W);
    radius_std = sqrt(sum(W .* (r_core - radius) .^ 2) / sum(W));

    base = struct( ...
        "valid", true, ...
        "x_min", x_min, ...
        "x_max", x_max, ...
        "x_center", x_center, ...
        "length", length_core, ...
        "radius", radius, ...
        "radius_area", radius_area, ...
        "radius_std", radius_std, ...
        "cylindricality_error", radius_std / radius, ...
        "xi", xi);
end

function window = axial_window(xi, cfg)
    if cfg.window.type == "none" || cfg.window.alpha <= 0
        window = ones(size(xi));
        return
    end
    if cfg.window.type ~= "tukey"
        error("CylinderMode:InvalidWindow", "Unknown window type: %s", cfg.window.type);
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

function [coeff, reconstructed, condition_number] = project_modes( ...
        phi, xi, eta, W, mode_m, mode_n, cfg)
    phase = phi * mode_m.' + 2 * pi * xi * mode_n.';
    Phi = exp(1i * phase);

    if cfg.modes.projection == "quadrature_nudft"
        coeff = Phi' * (W .* eta) / sum(W);
        gram = Phi' * (Phi .* W);
        condition_number = cond(gram);
    elseif cfg.modes.projection == "weighted_lstsq"
        sqrt_W = sqrt(W);
        weighted_basis = Phi .* sqrt_W;
        weighted_eta = eta .* sqrt_W;
        [U, S, V] = svd(weighted_basis, "econ");
        singular_values = diag(S);
        if isempty(singular_values) || singular_values(1) == 0
            coeff = complex(nan(size(mode_m)));
            reconstructed = nan(size(eta));
            condition_number = inf;
            return
        end
        lambda = cfg.modes.regularization * singular_values(1) ^ 2;
        filter = singular_values ./ (singular_values .^ 2 + lambda);
        coeff = V * (filter .* (U' * weighted_eta));
        condition_number = singular_values(1) / max(singular_values(end), eps);
    else
        error("CylinderMode:InvalidProjection", ...
            "Unknown projection method: %s", cfg.modes.projection);
    end
    reconstructed = real(Phi * coeff);
end

function resolved = mode_is_resolved(m, n, radius, length_core, mean_edge, cfg)
    if n == 0
        axial_samples = inf;
    else
        axial_samples = (length_core / abs(n)) / mean_edge;
    end
    if m == 0
        azimuthal_samples = inf;
    else
        azimuthal_samples = (2 * pi * radius / abs(m)) / mean_edge;
    end
    resolved = axial_samples >= cfg.quality.minimum_axial_samples_per_wavelength ...
        && azimuthal_samples >= cfg.quality.minimum_azimuthal_samples_per_wavelength;
end

function error_value = coefficient_conjugacy_error(coeff, mode_m, mode_n)
    errors = zeros(numel(coeff), 1);
    for k = 1:numel(coeff)
        partner = find(mode_m == -mode_m(k) & mode_n == -mode_n(k), 1);
        errors(k) = abs(coeff(k) - conj(coeff(partner)));
    end
    error_value = max(errors) / max(max(abs(coeff)), eps);
end

function value = weighted_relative_norm(residual, signal, W)
    numerator = sqrt(sum(W .* abs(residual) .^ 2));
    denominator = sqrt(sum(W .* abs(signal) .^ 2));
    value = numerator / max(denominator, eps);
end

function [mode_m, mode_n] = signed_mode_list(cfg)
    [m_grid, n_grid] = ndgrid(-cfg.modes.m_max:cfg.modes.m_max, ...
        -cfg.modes.n_max:cfg.modes.n_max);
    keep = abs(n_grid) >= cfg.modes.n_min | (m_grid == 0 & n_grid == 0);
    if cfg.modes.include_n_zero_nonaxisymmetric
        keep = keep | n_grid == 0;
    end
    mode_m = m_grid(keep);
    mode_n = n_grid(keep);
    mode_m = mode_m(:);
    mode_n = mode_n(:);
end

function [group_m, group_n, primary, conjugate] = grouped_mode_list( ...
        mode_m, mode_n, cfg)
    group_m = zeros(0, 1);
    group_n = zeros(0, 1);

    group_m(end + 1, 1) = 0;
    group_n(end + 1, 1) = 0;
    for n = max(1, cfg.modes.n_min):cfg.modes.n_max
        group_m(end + 1, 1) = 0; %#ok<AGROW>
        group_n(end + 1, 1) = n; %#ok<AGROW>
    end
    for m = 1:cfg.modes.m_max
        for n = -cfg.modes.n_max:cfg.modes.n_max
            include = abs(n) >= cfg.modes.n_min;
            if n == 0 && cfg.modes.include_n_zero_nonaxisymmetric
                include = true;
            end
            if include
                group_m(end + 1, 1) = m; %#ok<AGROW>
                group_n(end + 1, 1) = n; %#ok<AGROW>
            end
        end
    end

    primary = zeros(numel(group_m), 1);
    conjugate = zeros(numel(group_m), 1);
    for g = 1:numel(group_m)
        primary(g) = find(mode_m == group_m(g) & mode_n == group_n(g), 1);
        conjugate(g) = find(mode_m == -group_m(g) & mode_n == -group_n(g), 1);
    end
end

function table_out = make_frame_metrics_table(series_id, frame_index, time, metric)
    table_out = table( ...
        repmat(series_id, numel(frame_index), 1), frame_index, time, ...
        metric.dt, metric.n_vertices, metric.n_faces, ...
        metric.surface_area, metric.volume, metric.reduced_volume, ...
        metric.axis(:, 1), metric.axis(:, 2), metric.axis(:, 3), ...
        metric.center(:, 1), metric.center(:, 2), metric.center(:, 3), ...
        metric.core_x_min, metric.core_x_max, metric.core_length, ...
        metric.mean_core_radius, metric.radius_std, ...
        metric.cylindricality_error, metric.core_area_fraction, ...
        metric.core_axial_fraction, ...
        metric.windowed_area, metric.effective_vertices, ...
        metric.reconstruction_error, metric.projection_condition, ...
        metric.conjugacy_error, metric.mean_edge_length, ...
        metric.axis_angle_change, metric.frame_valid, metric.quality_flags, ...
        'VariableNames', [ ...
        "series_id", "frame_index", "time", "dt", "n_vertices", "n_faces", ...
        "surface_area", "volume", "reduced_volume", ...
        "axis_x", "axis_y", "axis_z", "center_x", "center_y", "center_z", ...
        "core_x_min", "core_x_max", "core_length", "mean_core_radius", ...
        "radius_std", "cylindricality_error", "core_area_fraction", ...
        "core_axial_fraction", ...
        "windowed_area", "effective_vertices", "reconstruction_error", ...
        "projection_condition", "conjugacy_error", "mean_edge_length", ...
        "axis_angle_change", "frame_valid", "quality_flags"]);
end

function table_out = make_mode_table( ...
        series_id, frame_index, time, coefficient, mode_m, mode_n, ...
        group_m, group_n, group_amplitude, group_noise, group_qR, ...
        core_length, core_radius, frame_flags)
    n_frames = numel(frame_index);
    n_modes = numel(mode_m);
    raw_group_index = zeros(n_modes, 1);
    for k = 1:n_modes
        [canonical_m, canonical_n] = canonical_mode(mode_m(k), mode_n(k));
        raw_group_index(k) = find(group_m == canonical_m & group_n == canonical_n, 1);
    end

    grouped_amp_raw = group_amplitude(:, raw_group_index);
    grouped_noise_raw = group_noise(:, raw_group_index);
    grouped_qR_raw = group_qR(:, raw_group_index);
    group_m_raw = group_m(raw_group_index);
    group_n_raw = group_n(raw_group_index);
    physical_q = 2 * pi * mode_n.' ./ core_length;

    table_out = table( ...
        repmat(series_id, n_frames * n_modes, 1), ...
        repelem(frame_index, n_modes), repelem(time, n_modes), ...
        repmat(mode_m, n_frames, 1), repmat(mode_n, n_frames, 1), ...
        repmat(group_m_raw, n_frames, 1), repmat(group_n_raw, n_frames, 1), ...
        reshape(physical_q.', [], 1), reshape(grouped_qR_raw.', [], 1), ...
        reshape(real(coefficient).', [], 1), ...
        reshape(imag(coefficient).', [], 1), ...
        reshape(grouped_amp_raw.', [], 1), ...
        reshape(grouped_amp_raw.'.^2, [], 1), ...
        reshape(angle(coefficient).', [], 1), ...
        reshape(grouped_noise_raw.', [], 1), ...
        repelem(core_radius, n_modes), ...
        repelem(frame_flags, n_modes), ...
        'VariableNames', [ ...
        "series_id", "frame_index", "time", "m", "n", "group_m", "group_n", ...
        "physical_wavenumber", "dimensionless_wavenumber_qR", ...
        "coefficient_real", "coefficient_imag", "grouped_amplitude", ...
        "grouped_power", "phase", "noise_estimate", "mean_core_radius", ...
        "quality_flags"]);
end

function [canonical_m, canonical_n] = canonical_mode(m, n)
    if m < 0
        canonical_m = -m;
        canonical_n = -n;
    elseif m == 0 && n < 0
        canonical_m = 0;
        canonical_n = -n;
    else
        canonical_m = m;
        canonical_n = n;
    end
end

function metric = initialize_metric_arrays(n)
    scalar_fields = [ ...
        "explicit_time", "dt", "n_vertices", "n_faces", "surface_area", ...
        "volume", "reduced_volume", "axis_angle_change", "core_x_min", ...
        "core_x_max", "core_length", "mean_core_radius", "radius_std", ...
        "cylindricality_error", "core_area_fraction", "windowed_area", ...
        "core_axial_fraction", ...
        "effective_vertices", "reconstruction_error", "projection_condition", ...
        "conjugacy_error", "mean_edge_length"];
    metric = struct();
    for field = scalar_fields
        metric.(field) = nan(n, 1);
    end
    metric.center = nan(n, 3);
    metric.axis = nan(n, 3);
    metric.frame_valid = false(n, 1);
    metric.quality_flags = strings(n, 1);
end

function time = resolve_physical_times(frame_index, explicit_time, dt)
    if all(isfinite(explicit_time))
        time = explicit_time;
        return
    end
    if any(~isfinite(dt))
        error("CylinderMode:MissingTime", ...
            "Frames lack both explicit physical time and p.dt.");
    end

    time = nan(size(frame_index));
    time(1) = frame_index(1) * dt(1);
    for k = 2:numel(time)
        time(k) = time(k - 1) + (frame_index(k) - frame_index(k - 1)) * dt(k);
    end
    explicit = isfinite(explicit_time);
    if any(explicit)
        offset = median(explicit_time(explicit) - time(explicit));
        time = time + offset;
        time(explicit) = explicit_time(explicit);
    end
end

function value = extract_explicit_time(data, cfg)
    value = NaN;
    for name = cfg.frame.time_field_names
        field = char(name);
        if isfield(data, field) && isscalar(data.(field)) && isfinite(data.(field))
            value = double(data.(field));
            return
        end
    end
end

function value = extract_dt(data)
    value = NaN;
    if isfield(data, "p") && isstruct(data.p) ...
            && isfield(data.p, "dt") && isscalar(data.p.dt)
        value = double(data.p.dt);
    end
end

function value = weighted_quantile(values, weights, probability)
    [values, order] = sort(values(:));
    weights = weights(order);
    cumulative = cumsum(weights) / sum(weights);
    index = find(cumulative >= probability, 1);
    value = values(index);
end

function runs = contiguous_runs(mask)
    edges = diff([false; mask(:); false]);
    runs = [find(edges == 1), find(edges == -1) - 1];
end

function base = empty_base()
    base = struct( ...
        "valid", false, "x_min", NaN, "x_max", NaN, "x_center", NaN, ...
        "length", NaN, "radius", NaN, "radius_area", NaN, ...
        "radius_std", NaN, "cylindricality_error", NaN, ...
        "xi", zeros(0, 1));
end

function name = get_folder_name(folder)
    pieces = split(string(folder), filesep);
    pieces = pieces(strlength(pieces) > 0);
    name = pieces(end);
end

function value = join_flags(flags)
    flags = unique(flags(strlength(flags) > 0), "stable");
    if isempty(flags)
        value = "";
    else
        value = join(flags, ";");
    end
end

function value = append_flag(value, flag)
    if strlength(value) == 0
        value = string(flag);
    elseif ~contains(";" + value + ";", ";" + flag + ";")
        value = value + ";" + flag;
    end
end
