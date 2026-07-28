function pipeline = cylinder_mode_pipeline(series_folders, cfg)
%CYLINDER_MODE_PIPELINE Run cylindrical mode analysis over simulation folders.

if nargin < 2 || isempty(cfg)
    cfg = cylinder_mode_defaults();
end
cfg = validate_config(cfg);
series_folders = resolve_series_folders(series_folders, cfg.data_root);
n_series = numel(series_folders);

if ~isfolder(cfg.output_root)
    mkdir(cfg.output_root);
end

shared_basis = [];
if cfg.pairing.enabled && cfg.modes.m_max > 0
    shared_basis = initial_basis_from_series(series_folders(1), cfg);
end

series_results = cell(n_series, 1);
summary_rows = cell(n_series, 6);
for k = 1:n_series
    folder = series_folders(k);
    series_id = folder_name(folder);
    output_dir = fullfile(cfg.output_root, series_id);
    if ~isfolder(output_dir)
        mkdir(output_dir);
    end

    if cfg.verbose
        fprintf("Cylinder mode analysis [%d/%d]: %s\n", k, n_series, series_id);
    end

    try
        cfg_series = cfg;
        if cfg.pairing.enabled && cfg.modes.m_max > 0
            cfg_series.alignment.axis_mode = "fixed_initial";
        end
        result = cylinder_mode_analyze_series(folder, cfg_series, shared_basis);
        write_series_outputs(result, output_dir, cfg_series);
        cylinder_mode_diagnostics(result, fullfile(output_dir, "cylinder_diagnostics"), cfg_series);
        series_results{k} = result;

        n_good = nnz(result.growth_rates.fit_status == "good");
        n_poor = nnz(result.growth_rates.fit_status == "poor");
        summary_rows(k, :) = {result.series_id, string(folder), "ok", ...
            height(result.frame_metrics), n_good, n_poor};
    catch exception
        summary_rows(k, :) = {string(series_id), string(folder), "failed", 0, 0, 0};
        series_results{k} = struct( ...
            "series_id", string(series_id), ...
            "series_folder", string(folder), ...
            "status", "failed", ...
            "error_message", string(getReport(exception, "extended", "hyperlinks", "off")));
        write_failure_log(output_dir, exception);
        if cfg.verbose
            fprintf(2, "  Failed: %s\n", exception.message);
        end
        if ~cfg.continue_on_error
            rethrow(exception)
        end
    end
end

summary = cell2table(summary_rows, 'VariableNames', ...
    ["series_id", "series_folder", "status", "n_frames", "n_good_fits", "n_poor_fits"]);
writetable(summary, fullfile(cfg.output_root, "cylinder_pipeline_summary.csv"));

response_results = {};
if cfg.pairing.enabled
    response_results = cylinder_mode_pair_responses(series_results, cfg);
    for k = 1:numel(response_results)
        response = response_results{k};
        output_dir = fullfile(cfg.output_root, char(response.series_id));
        if ~isfolder(output_dir)
            mkdir(output_dir);
        end
        write_series_outputs(response, output_dir, cfg);
        cylinder_mode_diagnostics(response, ...
            fullfile(output_dir, "cylinder_diagnostics"), cfg);
    end
end

pipeline = struct();
pipeline.config = cfg;
pipeline.series_folders = series_folders;
pipeline.series_results = series_results;
pipeline.response_results = response_results;
pipeline.summary = summary;
end

function cfg = validate_config(cfg)
    if cfg.modes.m_max < 0 || cfg.modes.m_max ~= round(cfg.modes.m_max)
        error("CylinderMode:InvalidConfig", "m_max must be a nonnegative integer.");
    end
    if cfg.modes.n_min < 1 || cfg.modes.n_max < cfg.modes.n_min ...
            || cfg.modes.n_max ~= round(cfg.modes.n_max)
        error("CylinderMode:InvalidConfig", ...
            "Require integer mode indices with 1 <= n_min <= n_max.");
    end
    if ~ismember(string(cfg.modes.projection), ...
            ["weighted_lstsq", "quadrature_nudft"])
        error("CylinderMode:InvalidConfig", "Unknown projection method.");
    end
    if ~ismember(string(cfg.core.method), ...
            ["auto", "persistent", "detect_each_frame"])
        error("CylinderMode:InvalidConfig", "Unknown core method.");
    end
    if ~ismember(string(cfg.alignment.axis_mode), ...
            ["tracked_pca", "fixed_initial", "known"])
        error("CylinderMode:InvalidConfig", "Unknown alignment axis mode.");
    end
    if cfg.frame.stride < 1 || cfg.frame.stride ~= round(cfg.frame.stride)
        error("CylinderMode:InvalidConfig", "Frame stride must be a positive integer.");
    end
end

function folders = resolve_series_folders(series_folders, data_root)
    series_folders = string(series_folders);
    series_folders = series_folders(:);
    if isempty(series_folders)
        error("CylinderMode:NoSeries", ...
            "Provide at least one simulation folder name or full path.");
    end
    folders = strings(size(series_folders));
    for k = 1:numel(series_folders)
        candidate = series_folders(k);
        if isfolder(candidate)
            folders(k) = string(java.io.File(char(candidate)).getCanonicalPath());
        else
            candidate = fullfile(data_root, candidate);
            if ~isfolder(candidate)
                error("CylinderMode:MissingSeries", ...
                    "Series folder not found: %s", series_folders(k));
            end
            folders(k) = string(java.io.File(char(candidate)).getCanonicalPath());
        end
    end
end

function write_series_outputs(result, output_dir, cfg)
    if cfg.output.save_frame_csv && isfield(result, "frame_metrics")
        writetable(result.frame_metrics, ...
            fullfile(output_dir, "cylinder_frame_metrics.csv"));
    end
    if cfg.output.save_mode_csv && isfield(result, "mode_coefficients")
        writetable(result.mode_coefficients, ...
            fullfile(output_dir, "cylinder_mode_coefficients.csv"));
    end
    if cfg.output.save_growth_csv && isfield(result, "growth_rates")
        writetable(result.growth_rates, ...
            fullfile(output_dir, "cylinder_growth_rates.csv"));
    end
    if cfg.output.save_mat
        save(fullfile(output_dir, "cylinder_mode_analysis.mat"), ...
            "result", "cfg", "-v7.3");
    end
    if cfg.output.save_resolved_json
        json = jsonencode(cfg, "PrettyPrint", true);
        file = fopen(fullfile(output_dir, "cylinder_resolved_config.json"), "w");
        if file < 0
            error("CylinderMode:OutputFailure", ...
                "Could not create resolved configuration file.");
        end
        cleanup = onCleanup(@() fclose(file));
        fwrite(file, json, "char");
        clear cleanup
    end
end

function write_failure_log(output_dir, exception)
    file = fopen(fullfile(output_dir, "cylinder_analysis_failure.txt"), "w");
    if file < 0
        return
    end
    cleanup = onCleanup(@() fclose(file));
    fprintf(file, "%s\n", getReport(exception, "extended", "hyperlinks", "off"));
    clear cleanup
end

function basis = initial_basis_from_series(folder, cfg)
    listing = dir(fullfile(folder, "geo*.mat"));
    indices = nan(numel(listing), 1);
    keep = false(numel(listing), 1);
    for k = 1:numel(listing)
        token = regexp(listing(k).name, "^geo(\d+)\.mat$", "tokens", "once");
        if ~isempty(token)
            indices(k) = str2double(token{1});
            keep(k) = true;
        end
    end
    listing = listing(keep);
    indices = indices(keep);
    [~, first] = min(indices);
    data = load(fullfile(listing(first).folder, listing(first).name), "P", "M");
    P = double(data.P);
    M = double(data.M);

    v1 = P(M(:, 1), :);
    v2 = P(M(:, 2), :);
    v3 = P(M(:, 3), :);
    area = 0.5 * vecnorm(cross(v2 - v1, v3 - v1, 2), 2, 2);
    vertex_area = accumarray(M(:), repmat(area / 3, 3, 1), ...
        [size(P, 1), 1], @sum, 0);
    center = sum(P .* vertex_area, 1) / sum(vertex_area);
    centered = P - center;

    if cfg.alignment.axis_mode == "known" && ~isempty(cfg.alignment.known_axis)
        axis_vector = cfg.alignment.known_axis(:);
    else
        covariance = centered.' * (centered .* vertex_area) / sum(vertex_area);
        [vectors, values] = eig((covariance + covariance.') / 2, "vector");
        [~, largest] = max(values);
        axis_vector = vectors(:, largest);
    end
    axis_vector = axis_vector / norm(axis_vector);
    coordinate_axes = eye(3);
    [~, least_parallel] = min(abs(coordinate_axes.' * axis_vector));
    transverse = coordinate_axes(:, least_parallel);
    transverse = transverse - axis_vector * dot(axis_vector, transverse);
    transverse = transverse / norm(transverse);
    binormal = cross(axis_vector, transverse);
    binormal = binormal / norm(binormal);
    basis = [axis_vector, transverse, binormal];
end

function name = folder_name(folder)
    pieces = split(string(folder), filesep);
    pieces = pieces(strlength(pieces) > 0);
    name = char(pieces(end));
end
