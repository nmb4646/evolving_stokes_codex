close all; clc;

% Comparison to Shaqfeh and Zhao (2011) Figure 4, steady-state tilt angle
% and deformation index of vesicles in shear flow



%%% Simulation selection
% Edit these values to select runs from data/fs_batch_data.
Da = 0;
Sd = [1e-6];
gamy = [8e-7];

% First geo timestep to include. Use 0 to include geo0.mat.
time_start = 1;
% Last geo timestep to include. Use inf to include all later frames.
time_final = inf;

%%% Tilt calculation
% For shear_flow.m, u_x = gamy * z, so flow_dim = 1 and grad_dim = 3.
flow_dim = 1;
grad_dim = 3;
axis_length_method = "ray_intersection"; % Options: "ray_intersection" or "projection".

% If true, remove +/- pi branch jumps using the axis-periodic angle.
% Leave false to match the wrapped psi/pi values printed by fs_plotter.
tilt_options.unwrap = false;

%%% Visual settings
% Options: "tilt", "D", or ["tilt", "D"].
visual.metrics = ["tilt", "D"];
% Options: "over_pi" matches fs_plotter, "radians", or "degrees".
visual.angle_units = "over_pi";
visual.line_width = 2.0;
visual.marker = "none";
visual.marker_size = 5;
visual.figure_position = [140, 140, 900, 560];
visual.scale_time_by_gamy = true; % If true, plot t_scaled = t * gamy_current.
visual.time_scale_by_Sd = true; % If true, plot t_scaled = t / Sd_current.
visual.x_scale = "linear"; % Options: "linear" or "log".
visual.y_scale = "linear"; % Options: "linear" or "log".
visual.show_grid = true;
visual.save_figure = false;
visual.output_file = fullfile("data", "tilt_over_time.png");
visual.output_resolution = 300;
visual.D_label = "Deformation index D";

script_dir = fileparts(mfilename("fullpath"));
remesh_dir = find_remesh_dir(script_dir);
addpath(remesh_dir);

data_dir = find_data_dir(script_dir, remesh_dir);
sd_values = Sd(:).';
gamy_values = gamy(:).';
multiple_sd = numel(sd_values) > 1;
multiple_gamy = numel(gamy_values) > 1;
run_count = numel(sd_values) * numel(gamy_values);
multiple_runs = run_count > 1;
metrics = normalize_metrics(visual.metrics);
plot_count = numel(metrics);

fig = figure("Color", "w", "Position", visual.figure_position);
if plot_count == 1
    axes_list = axes(fig);
else
    tiled = tiledlayout(fig, plot_count, 1, "TileSpacing", "compact", "Padding", "compact");
    axes_list = gobjects(plot_count, 1);
    for metric_idx = 1:plot_count
        axes_list(metric_idx) = nexttile(tiled);
    end
end
for metric_idx = 1:plot_count
    hold(axes_list(metric_idx), "on");
end

legend_entries = strings(run_count, 1);
x_label = "";
metric_labels = strings(plot_count, 1);
run_idx = 0;

for sd_idx = 1:numel(sd_values)
    Sd_current = sd_values(sd_idx);

    for gamy_idx = 1:numel(gamy_values)
        gamy_current = gamy_values(gamy_idx);
        run_idx = run_idx + 1;
        run_tag = make_run_tag(Sd_current, Da, gamy_current);
        folder = fullfile(data_dir, run_tag);

        fprintf("Loading %s\n", folder);

        [x, x_label_current, tilt, D] = load_tilt_deformation_series(folder, flow_dim, grad_dim, ...
            axis_length_method, tilt_options, time_start, time_final);
        [x, x_label_current] = apply_time_scaling(x, x_label_current, ...
            gamy_current, Sd_current, visual);

        if strlength(x_label) == 0
            x_label = x_label_current;
        elseif x_label ~= x_label_current
            x_label = "Time / timestep";
        end

        for metric_idx = 1:plot_count
            [y, metric_labels(metric_idx)] = metric_to_plot_values( ...
                metrics(metric_idx), tilt, D, visual);
            plot(axes_list(metric_idx), x, y, ...
                "LineWidth", visual.line_width, ...
                "Marker", visual.marker, ...
                "MarkerSize", visual.marker_size);
        end

        legend_entries(run_idx) = make_legend_entry(Sd_current, gamy_current, ...
            multiple_sd, multiple_gamy);

        fprintf("Sd = %.8g, gamy = %.8g, initial tilt: %.8g psi/pi, %.8g rad, %.8g deg\n", ...
            Sd_current, gamy_current, tilt(1) / pi, tilt(1), rad2deg(tilt(1)));
        fprintf("Sd = %.8g, gamy = %.8g, final tilt:   %.8g psi/pi, %.8g rad, %.8g deg\n", ...
            Sd_current, gamy_current, tilt(end) / pi, tilt(end), rad2deg(tilt(end)));
        fprintf("Sd = %.8g, gamy = %.8g, initial D: %.8g, final D: %.8g\n", ...
            Sd_current, gamy_current, D(1), D(end));
    end
end

for metric_idx = 1:plot_count
    ax = axes_list(metric_idx);
    ylabel(ax, metric_labels(metric_idx));
    set(ax, "XScale", visual.x_scale);
    set(ax, "YScale", visual.y_scale);

    if visual.show_grid
        grid(ax, "on");
    end
    box(ax, "on");
end
xlabel(axes_list(end), x_label);

if multiple_runs
    title(axes_list(1), sprintf("%s over time: %s", ...
        metrics_title(metrics), comparison_title_suffix(Da, sd_values, gamy_values)));
    legend(axes_list(1), legend_entries, "Location", "best");
else
    title(axes_list(1), sprintf("%s over time: Sd = %.4g, Da = %.4g, gamy = %.4g", ...
        metrics_title(metrics), sd_values, Da, gamy_values));
end

if visual.save_figure
    exportgraphics(fig, visual.output_file, "Resolution", visual.output_resolution);
    fprintf("Saved %s\n", visual.output_file);
end

function [x, x_label, tilt, D] = load_tilt_deformation_series(folder, flow_dim, grad_dim, ...
        axis_length_method, tilt_options, time_start, time_final)
    files = dir(fullfile(folder, "geo*.mat"));
    if isempty(files)
        error("No geo*.mat files found in %s", folder);
    end

    frame_ids = frame_ids_from_files(files);
    [frame_ids, order] = sort(frame_ids);
    files = files(order);

    valid = ~isnan(frame_ids);
    frame_ids = frame_ids(valid);
    files = files(valid);

    if time_final < time_start
        error("time_final must be greater than or equal to time_start.");
    end

    keep = frame_ids >= time_start & frame_ids <= time_final;
    frame_ids = frame_ids(keep);
    files = files(keep);

    if isempty(frame_ids)
        error("No valid geo*.mat frame numbers from time_start = %.17g through time_final = %.17g found in %s", ...
            time_start, time_final, folder);
    end

    tilt = zeros(numel(frame_ids), 1);
    D = zeros(numel(frame_ids), 1);
    dt = NaN;

    for i = 1:numel(frame_ids)
        frame_path = fullfile(folder, files(i).name);
        data = load(frame_path, "M", "P", "p");

        if i == 1 && isfield(data, "p") && isfield(data.p, "dt")
            dt = data.p.dt;
        end

        tilt_out = vesicleTiltDeformation(data.P, data.M, flow_dim, grad_dim, axis_length_method);
        tilt(i) = tilt_out.psi;
        D(i) = tilt_out.D;
    end

    if tilt_options.unwrap
        tilt = 0.5 * unwrap(2 * tilt);
    end

    if isnan(dt)
        x = frame_ids;
        x_label = "Timestep";
    else
        x = frame_ids * dt;
        x_label = "Time";
    end
end

function [x, x_label] = apply_time_scaling(x, base_label, gamy, Sd, visual)
    scale_by_gamy = visual.scale_time_by_gamy;
    scale_by_Sd = visual.time_scale_by_Sd;

    if ~scale_by_gamy && ~scale_by_Sd
        x_label = base_label;
        return
    end

    if scale_by_Sd && Sd == 0
        error("Cannot scale time by Sd because Sd is zero.");
    end

    if scale_by_gamy
        x = x * gamy;
    end

    if scale_by_Sd
        x = x / Sd;
    end

    if base_label == "Time"
        quantity = "time";
        symbol = "t";
    elseif base_label == "Timestep"
        quantity = "timestep";
        symbol = "n";
    else
        quantity = "time / timestep";
        symbol = "x";
    end

    if scale_by_gamy && scale_by_Sd
        expression = "gamy " + symbol + " / Sd";
    elseif scale_by_gamy
        expression = "gamy " + symbol;
    else
        expression = symbol + " / Sd";
    end

    x_label = "Scaled " + quantity + " (" + expression + ")";
end

function metrics = normalize_metrics(metrics)
    metrics = lower(string(metrics));
    if isempty(metrics)
        error('visual.metrics must include "tilt", "D", or both.');
    end

    for i = 1:numel(metrics)
        switch metrics(i)
            case {"tilt", "psi", "angle"}
                metrics(i) = "tilt";
            case {"d", "deformation", "deformation_index"}
                metrics(i) = "D";
            otherwise
                error('Unknown visual metric "%s". Use "tilt", "D", or ["tilt", "D"].', metrics(i));
        end
    end
    metrics = unique(metrics, "stable");
end

function [y, y_label] = metric_to_plot_values(metric, tilt, D, visual)
    switch metric
        case "tilt"
            [y, y_label] = tilt_to_plot_values(tilt, visual.angle_units);
        case "D"
            y = D;
            y_label = visual.D_label;
    end
end

function [y, y_label] = tilt_to_plot_values(tilt, angle_units)
    switch lower(string(angle_units))
        case {"over_pi", "pi", "psi_over_pi"}
            y = tilt / pi;
            y_label = "Tilt angle (\psi/\pi)";
        case {"radians", "rad"}
            y = tilt;
            y_label = "Tilt angle (rad)";
        case {"degrees", "deg"}
            y = rad2deg(tilt);
            y_label = "Tilt angle (deg)";
        otherwise
            error("Unknown visual.angle_units '%s'. Use 'over_pi', 'radians', or 'degrees'.", angle_units);
    end
end

function title_text = metrics_title(metrics)
    if isequal(metrics, "tilt")
        title_text = "Tilt";
    elseif isequal(metrics, "D")
        title_text = "Deformation";
    else
        title_text = "Tilt and deformation";
    end
end

function run_tag = make_run_tag(Sd, Da, gamy)
    run_tag = sprintf("Sd_%.2e_Da_%.2e_gamy_%+.2e", Sd, Da, gamy);
    run_tag = strrep(run_tag, "+", "p");
    run_tag = strrep(run_tag, "-", "m");
end

function label = make_legend_entry(Sd, gamy, show_sd, show_gamy)
    parts = strings(0, 1);

    if show_sd
        parts(end + 1) = sprintf("Sd = %.4g", Sd);
    end

    if show_gamy
        parts(end + 1) = sprintf("gamy = %.4g", gamy);
    end

    if isempty(parts)
        label = sprintf("Sd = %.4g, gamy = %.4g", Sd, gamy);
    else
        label = strjoin(parts, ", ");
    end
end

function suffix = comparison_title_suffix(Da, sd_values, gamy_values)
    parts = "Da = " + sprintf("%.4g", Da);

    if isscalar(sd_values)
        parts(end + 1) = "Sd = " + sprintf("%.4g", sd_values);
    end

    if isscalar(gamy_values)
        parts(end + 1) = "gamy = " + sprintf("%.4g", gamy_values);
    end

    suffix = strjoin(parts, ", ");
end

function frame_ids = frame_ids_from_files(files)
    frame_ids = NaN(numel(files), 1);

    for i = 1:numel(files)
        token = regexp(files(i).name, "geo(\d+)\.mat", "tokens", "once");
        if ~isempty(token)
            frame_ids(i) = str2double(token{1});
        end
    end
end

function remesh_dir = find_remesh_dir(script_dir)
    candidates = string({script_dir, pwd, fullfile(pwd, "remesh")});
    geometry_path = which("Geometry");
    if strlength(geometry_path) > 0
        candidates(end + 1) = string(fileparts(geometry_path));
    end

    for i = 1:numel(candidates)
        candidate = candidates(i);
        if isfile(fullfile(candidate, "Geometry.m"))
            remesh_dir = char(candidate);
            return
        end
    end

    error("Could not locate remesh directory. Run from the repo root or add remesh/ to the MATLAB path.");
end

function data_dir = find_data_dir(script_dir, remesh_dir)
    candidates = string({
        fullfile(remesh_dir, "data", "fs_batch_data"), ...
        fullfile(script_dir, "data", "fs_batch_data"), ...
        fullfile(pwd, "data", "fs_batch_data"), ...
        fullfile(pwd, "remesh", "data", "fs_batch_data")});

    for i = 1:numel(candidates)
        candidate = candidates(i);
        if isfolder(candidate)
            data_dir = char(candidate);
            return
        end
    end

    error("Could not locate data/fs_batch_data. Checked:%s", sprintf("\n  %s", candidates));
end
