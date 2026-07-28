%CYLINDER_ONEMODE_OVERTIME Plot one cylindrical mode amplitude across Da.
%
% Edit the parameters below, then run this script. It analyzes the standard
% fs_batch_data folder for each Da and plots the selected m = 0 mode against
% physical time.

clearvars;
close all;
clc;

%%% Simulation and mode parameters

Sd = 1e-6;
Da_values = [1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1];
gamy = 4e-6;
v = 0.35;

mode_n = 1;                    % Positive axial mode number to plot.
analysis_n_max = 6;           % Retain this many modes in the joint projection.

%%% Frame range

first_frame = 0;
last_frame = inf;
frame_stride = 1;
maximum_frames = inf;

%%% Cylindrical decomposition

projection_method = "weighted_lstsq"; % "weighted_lstsq" or "quadrature_nudft".
core_method = "auto";                 % "auto", "persistent", or "detect_each_frame".
axis_mode = "tracked_pca";            % "tracked_pca", "fixed_initial", or "known".
known_axis = [];                      % Example: [1, 0, 0].
tukey_alpha = .055;

%%% Reliability masks

mask_invalid_frames = false;   % True removes frames rejected by core/geometry checks.
mask_unresolved_frames = false; % True removes frames lacking spatial resolution for mode n.

%%% Plot appearance

figure_visible = "on";
figure_position = [100, 100, 1050, 720];
figure_resolution = 250;
figure_background_color = [1, 1, 1];

x_axis_scale = "log";
y_axis_scale = "log";
x_axis_limits = [];            % Empty selects limits automatically.
y_axis_limits = [];            % Empty selects limits automatically.
x_axis_label = "Physical time";
y_axis_label = "Grouped mode amplitude";
axes_font_size = 14;
label_font_size = 17;
title_font_size = 17;
axes_line_width = 1.1;
grid_alpha = 0.18;

data_line_width = 2;
marker_symbol = "none";
marker_size = 6;
Da_color_order = [];           % Empty uses MATLAB's lines map; otherwise N-by-3 RGB.

show_legend = true;
legend_location = "eastoutside";
legend_font_size = 12;
legend_box = "off";

%%% Saved outputs

save_figure = true;
save_matlab_figure = true;
save_combined_mat = true;
save_amplitude_csv = true;
save_per_series_diagnostics = false;

%%% Resolve requested simulation folders

script_dir = fileparts(mfilename("fullpath"));
addpath(script_dir);
cfg = cylinder_mode_defaults();

validateattributes(Sd, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'}, mfilename, "Sd");
validateattributes(Da_values, {'numeric'}, ...
    {'vector', 'real', 'finite', 'nonnegative', 'nonempty'}, ...
    mfilename, "Da_values");
validateattributes(gamy, {'numeric'}, ...
    {'scalar', 'real', 'finite'}, mfilename, "gamy");
validateattributes(v, {'numeric'}, ...
    {'scalar', 'real', 'finite'}, mfilename, "v");
validateattributes(mode_n, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'integer', 'positive'}, mfilename, "mode_n");
validateattributes(analysis_n_max, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'integer', 'positive'}, ...
    mfilename, "analysis_n_max");

analysis_n_max = max(analysis_n_max, mode_n);
Da_requested = unique(double(Da_values(:)), "stable");
series_folders = strings(numel(Da_requested), 1);
for simulation_index = 1:numel(Da_requested)
    series_folders(simulation_index) = cylinder_parameter_tag( ...
        Sd, Da_requested(simulation_index), gamy, v);
end

folder_exists = false(size(series_folders));
for simulation_index = 1:numel(series_folders)
    folder_exists(simulation_index) = isfolder( ...
        fullfile(cfg.data_root, series_folders(simulation_index)));
    if ~folder_exists(simulation_index)
        warning("CylinderMode:MissingOneModeSeries", ...
            "Skipping missing series folder: %s", series_folders(simulation_index));
    end
end
if ~any(folder_exists)
    error("CylinderMode:NoOneModeSeries", ...
        "None of the requested simulation folders exist under %s.", cfg.data_root);
end

Da_used = Da_requested(folder_exists);
series_folders = series_folders(folder_exists);

output_tag = cylinder_onemode_tag(Sd, gamy, v, mode_n);
cfg.output_root = fullfile(script_dir, "data", ...
    "cylinder_onemode_overtime", output_tag);

%%% Configure and run the established mode-analysis pipeline

cfg.frame.first_index = first_frame;
cfg.frame.last_index = last_frame;
cfg.frame.stride = frame_stride;
cfg.frame.max_frames = maximum_frames;

cfg.modes.m_max = 0;
cfg.modes.n_min = 1;
cfg.modes.n_max = analysis_n_max;
cfg.modes.projection = projection_method;

cfg.core.method = core_method;
cfg.alignment.axis_mode = axis_mode;
cfg.alignment.known_axis = known_axis;
cfg.window.alpha = tukey_alpha;

cfg.output.save_frame_csv = false;
cfg.output.save_mode_csv = false;
cfg.output.save_growth_csv = false;
cfg.output.save_mat = false;
cfg.output.save_resolved_json = false;
cfg.diagnostics.enabled = save_per_series_diagnostics;
cfg.diagnostics.save_individual_fit_plots = false;
cfg.diagnostics.figure_visible = "off";

cfg.core.minimum_radius_fraction = 0.80;
cfg.core.transition_margin_in_radii = 3.0;
cfg.core.fallback_half_length_fraction = 0.20;
cfg.core.profile_smoothing_bins = 3;
cfg.core.method = core_method;

cylinder_results = cylinder_mode_pipeline(series_folders, cfg);

%%% Extract the requested amplitude from each successful series

all_mode_amplitudes = table();
plot_series = cell(numel(Da_used), 1);
series_successful = false(numel(Da_used), 1);

for simulation_index = 1:numel(cylinder_results.series_results)
    result = cylinder_results.series_results{simulation_index};
    if ~isfield(result, "group_amplitude")
        warning("CylinderMode:FailedOneModeSeries", ...
            "No mode amplitudes were produced for Da = %.6g.", Da_used(simulation_index));
        continue
    end

    group_index = find(result.group_m == 0 & result.group_n == mode_n, 1);
    if isempty(group_index)
        warning("CylinderMode:MissingRequestedMode", ...
            "Mode m = 0, n = %d was not produced for Da = %.6g.", ...
            mode_n, Da_used(simulation_index));
        continue
    end

    amplitude = result.group_amplitude(:, group_index);
    noise = result.group_noise(:, group_index);
    resolved = result.group_resolved(:, group_index);
    valid = result.frame_metrics.frame_valid;
    n_frames = numel(result.time);

    series_table = table( ...
        repmat(Sd, n_frames, 1), ...
        repmat(Da_used(simulation_index), n_frames, 1), ...
        repmat(gamy, n_frames, 1), ...
        repmat(v, n_frames, 1), ...
        zeros(n_frames, 1), ...
        repmat(mode_n, n_frames, 1), ...
        result.frame_index(:), ...
        result.time(:), ...
        amplitude(:), ...
        noise(:), ...
        valid(:), ...
        resolved(:), ...
        result.frame_metrics.quality_flags, ...
        'VariableNames', [ ...
        "Sd", "Da", "gamy", "v", "m", "n", "frame_index", "time", ...
        "amplitude", "noise_estimate", "frame_valid", "mode_resolved", ...
        "quality_flags"]);

    if isempty(all_mode_amplitudes)
        all_mode_amplitudes = series_table;
    else
        all_mode_amplitudes = [all_mode_amplitudes; series_table]; %#ok<AGROW>
    end

    plot_series{simulation_index} = struct( ...
        "Da", Da_used(simulation_index), ...
        "time", result.time(:), ...
        "amplitude", amplitude(:), ...
        "valid", valid(:), ...
        "resolved", resolved(:));
    series_successful(simulation_index) = true;

    fprintf("Da = %.6g: extracted %d frames; valid = %d, resolved = %d.\n", ...
        Da_used(simulation_index), n_frames, nnz(valid), nnz(resolved));
end

if ~any(series_successful)
    error("CylinderMode:NoOneModeAmplitudes", ...
        "Every available simulation failed before producing the requested mode.");
end

if ~isfolder(cfg.output_root)
    mkdir(cfg.output_root);
end
if save_amplitude_csv
    writetable(all_mode_amplitudes, ...
        fullfile(cfg.output_root, "cylinder_onemode_amplitudes.csv"));
end

%%% Plot one amplitude history per successful Da

Da_plotted = Da_used(series_successful);
plot_series = plot_series(series_successful);

if isempty(Da_color_order)
    plot_colors = lines(numel(Da_plotted));
else
    validateattributes(Da_color_order, {'numeric'}, ...
        {'2d', 'ncols', 3, '>=', 0, '<=', 1}, ...
        mfilename, "Da_color_order");
    repeats = ceil(numel(Da_plotted) / size(Da_color_order, 1));
    plot_colors = repmat(Da_color_order, repeats, 1);
    plot_colors = plot_colors(1:numel(Da_plotted), :);
end

fig = figure( ...
    "Color", figure_background_color, ...
    "Visible", figure_visible, ...
    "Position", figure_position);
ax = axes(fig);
hold(ax, "on");

for simulation_index = 1:numel(plot_series)
    series = plot_series{simulation_index};
    amplitude = series.amplitude;
    if mask_invalid_frames
        amplitude(~series.valid) = NaN;
    end
    if mask_unresolved_frames
        amplitude(~series.resolved) = NaN;
    end
    if y_axis_scale == "log"
        amplitude(amplitude <= 0) = NaN;
    end

    plot(ax, series.time, amplitude, ...
        "LineStyle", "-", ...
        "LineWidth", data_line_width, ...
        "Marker", marker_symbol, ...
        "MarkerSize", marker_size, ...
        "Color", plot_colors(simulation_index, :), ...
        "DisplayName", sprintf("D_a = %.4g", series.Da));
end

set(ax, ...
    "XScale", x_axis_scale, ...
    "YScale", y_axis_scale, ...
    "FontSize", axes_font_size, ...
    "LineWidth", axes_line_width, ...
    "Box", "on", ...
    "XMinorGrid", "on", ...
    "YMinorGrid", "on", ...
    "GridAlpha", grid_alpha, ...
    "MinorGridAlpha", grid_alpha);
grid(ax, "on");

if ~isempty(x_axis_limits)
    xlim(ax, x_axis_limits);
end
if ~isempty(y_axis_limits)
    ylim(ax, y_axis_limits);
end

xlabel(ax, x_axis_label, "FontSize", label_font_size);
ylabel(ax, y_axis_label, "FontSize", label_font_size);
title(ax, sprintf("Cylindrical mode amplitude, m = 0, n = %d", mode_n), ...
    "FontSize", title_font_size, ...
    "FontWeight", "normal");

if show_legend
    legend_handle = legend(ax, ...
        "Location", legend_location, ...
        "FontSize", legend_font_size, ...
        "Box", legend_box);
    legend_handle.AutoUpdate = "off";
end

if save_figure
    exportgraphics(fig, fullfile(cfg.output_root, ...
        "cylinder_onemode_overtime.png"), ...
        "Resolution", figure_resolution);
end
if save_matlab_figure
    savefig(fig, fullfile(cfg.output_root, ...
        "cylinder_onemode_overtime.fig"));
end
if save_combined_mat
    save(fullfile(cfg.output_root, "cylinder_onemode_overtime.mat"), ...
        "all_mode_amplitudes", "plot_series", "Da_used", "Da_plotted", ...
        "Sd", "gamy", "v", "mode_n", "cfg");
end

fprintf("Outputs: %s\n", cfg.output_root);

function tag = cylinder_parameter_tag(Sd, Da, gamy, v)
    tag = sprintf("Sd_%.2e_Da_%.2e_gamy_%+.2e_v_%.2e", ...
        Sd, Da, gamy, v);
    tag = strrep(tag, "+", "p");
    tag = strrep(tag, "-", "m");
end

function tag = cylinder_onemode_tag(Sd, gamy, v, mode_n)
    tag = sprintf("Sd_%.2e_gamy_%+.2e_v_%.2e_m_0_n_%d", ...
        Sd, gamy, v, mode_n);
    tag = strrep(tag, "+", "p");
    tag = strrep(tag, "-", "m");
end
