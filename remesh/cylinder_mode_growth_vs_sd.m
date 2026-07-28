%CYLINDER_MODE_GROWTH_VS_SD Plot cylindrical-mode growth rates versus Sd.
%
% For one fixed capillary number Ca = gamy/Sd, this script pairs every Sd
% with gamy = Ca*Sd, analyzes the corresponding simulation, and plots one
% growth-rate curve for each retained cylindrical mode (m,n).

clearvars;
close all;
clc;

%%% Physical sweep

Sd_values = 10.^(-8:-4);     % Sd values to analyze, in any order.
Ca = 3.9;                    % Fixed capillary number: Ca = gamy/Sd.
Da = 0;
v = 0.35;

%%% Modes and frame range

m_max = 0;                   % 0 selects only axisymmetric modes.
n_min = 1;
n_max = 16;
modes_to_plot = [];          % Empty plots all fits; or use [m1,n1; m2,n2].

first_frame = 0;
last_frame = inf;
frame_stride = 1;
maximum_frames = inf;

%%% Cylindrical decomposition

projection_method = "weighted_lstsq"; % "weighted_lstsq" or "quadrature_nudft".
core_method = "auto";                 % "auto", "persistent", or "detect_each_frame".
axis_mode = "tracked_pca";            % "tracked_pca", "fixed_initial", or "known".
known_axis = [];                      % Example: [1, 0, 0].
tukey_alpha = 0.25;

%%% Growth-rate fitting

growth_fit_mode = "auto";     % "auto" or "fixed".
fixed_fit_start_time = -inf;
fixed_fit_end_time = inf;
minimum_fit_points = 6;
maximum_fit_points = 40;
maximum_early_frames = 100;
minimum_amplitude_ratio = 1.20;
minimum_signal_to_noise = 3.0;
maximum_dimensionless_amplitude = 0.05;
minimum_r_squared = 0.90;
include_poor_fits = false;    % False plots only fits classified as "good".
minimum_Sd_points_per_mode = 2;

%%% Plot appearance

figure_visible = "on";
figure_position = [100, 100, 1050, 720];
figure_resolution = 250;
figure_background_color = [1, 1, 1];

x_axis_scale = "log";
y_axis_scale = "log";        % "log" displays only positive growth rates.
x_axis_label = "S_d";
y_axis_label = "Growth rate, \sigma";
axes_font_size = 14;
label_font_size = 17;
title_font_size = 17;
axes_line_width = 1.1;
grid_alpha = 0.18;

connect_points = true;
data_line_width = 1.8;
marker_symbol = "o";
marker_size = 7;
mode_color_order = [];        % Empty uses MATLAB's lines map; otherwise N-by-3 RGB.
show_zero_line = true;

show_legend = true;
legend_location = "eastoutside";
legend_font_size = 11;
legend_num_columns = 1;
legend_box = "off";

%%% Saved outputs

save_figure = true;
save_matlab_figure = true;
save_combined_mat = true;
save_per_series_diagnostics = false;
save_individual_fit_plots = false;
save_mode_coefficient_csv = false;
save_full_series_mat = false;

%%% Build the one-to-one Sd/gamy sweep

script_dir = fileparts(mfilename("fullpath"));
addpath(script_dir);
cfg = cylinder_mode_defaults();

validateattributes(Sd_values, {'numeric'}, ...
    {'vector', 'real', 'finite', 'positive', 'nonempty'}, mfilename, "Sd_values");
validateattributes(Ca, {'numeric'}, ...
    {'scalar', 'real', 'finite'}, mfilename, "Ca");
validateattributes(Da, {'numeric'}, ...
    {'scalar', 'real', 'finite'}, mfilename, "Da");
validateattributes(v, {'numeric'}, ...
    {'scalar', 'real', 'finite'}, mfilename, "v");

Sd_requested = double(Sd_values(:));
Sd_requested = sort(Sd_requested);
gamy_requested = Ca .* Sd_requested;

series_folders = strings(numel(Sd_requested), 1);
for simulation_index = 1:numel(Sd_requested)
    series_folders(simulation_index) = cylinder_parameter_tag( ...
        Sd_requested(simulation_index), Da, ...
        gamy_requested(simulation_index), v);
end

folder_exists = false(size(series_folders));
for simulation_index = 1:numel(series_folders)
    folder_exists(simulation_index) = isfolder( ...
        fullfile(cfg.data_root, series_folders(simulation_index)));
    if ~folder_exists(simulation_index)
        warning("CylinderMode:MissingSweepSeries", ...
            "Skipping missing series folder: %s", series_folders(simulation_index));
    end
end
if ~any(folder_exists)
    error("CylinderMode:NoSweepSeries", ...
        "None of the requested fixed-Ca simulation folders exist under %s.", ...
        cfg.data_root);
end

Sd_used = Sd_requested(folder_exists);
gamy_used = gamy_requested(folder_exists);
series_folders = series_folders(folder_exists);

sweep_tag = cylinder_sweep_tag(Ca, Da, v);
cfg.output_root = fullfile(script_dir, "data", ...
    "cylinder_mode_growth_vs_sd", sweep_tag);

cfg.frame.first_index = first_frame;
cfg.frame.last_index = last_frame;
cfg.frame.stride = frame_stride;
cfg.frame.max_frames = maximum_frames;

cfg.modes.m_max = m_max;
cfg.modes.n_min = n_min;
cfg.modes.n_max = n_max;
cfg.modes.projection = projection_method;

cfg.core.method = core_method;
cfg.alignment.axis_mode = axis_mode;
cfg.alignment.known_axis = known_axis;
cfg.window.alpha = tukey_alpha;

cfg.growth.mode = growth_fit_mode;
cfg.growth.fixed_start_time = fixed_fit_start_time;
cfg.growth.fixed_end_time = fixed_fit_end_time;
cfg.growth.minimum_points = minimum_fit_points;
cfg.growth.maximum_points = maximum_fit_points;
cfg.growth.maximum_early_frames = maximum_early_frames;
cfg.growth.minimum_amplitude_ratio = minimum_amplitude_ratio;
cfg.growth.minimum_signal_to_noise = minimum_signal_to_noise;
cfg.growth.maximum_dimensionless_amplitude = maximum_dimensionless_amplitude;
cfg.growth.minimum_r_squared = minimum_r_squared;

cfg.output.save_frame_csv = true;
cfg.output.save_mode_csv = save_mode_coefficient_csv;
cfg.output.save_growth_csv = true;
cfg.output.save_mat = save_full_series_mat;
cfg.diagnostics.enabled = save_per_series_diagnostics;
cfg.diagnostics.save_individual_fit_plots = save_individual_fit_plots;
cfg.diagnostics.figure_visible = "off";

cylinder_results = cylinder_mode_pipeline(series_folders, cfg);

%%% Collect all mode fits into one long table

all_growth_rates = table();
for simulation_index = 1:numel(cylinder_results.series_results)
    result = cylinder_results.series_results{simulation_index};
    if ~isfield(result, "growth_rates")
        continue
    end

    rates = result.growth_rates;
    row_count = height(rates);
    rates.Sd = repmat(Sd_used(simulation_index), row_count, 1);
    rates.Da = repmat(Da, row_count, 1);
    rates.gamy = repmat(gamy_used(simulation_index), row_count, 1);
    rates.Ca = repmat(Ca, row_count, 1);
    rates = movevars(rates, ["Sd", "Da", "gamy", "Ca"], "Before", 1);

    if isempty(all_growth_rates)
        all_growth_rates = rates;
    else
        all_growth_rates = [all_growth_rates; rates]; %#ok<AGROW>
    end
end

if isempty(all_growth_rates)
    error("CylinderMode:NoGrowthRates", ...
        "Every available simulation failed before producing growth-rate data.");
end

sweep_summary = cylinder_results.summary;
sweep_summary.Sd = Sd_used;
sweep_summary.Da = repmat(Da, height(sweep_summary), 1);
sweep_summary.gamy = gamy_used;
sweep_summary.Ca = repmat(Ca, height(sweep_summary), 1);
sweep_summary = movevars(sweep_summary, ...
    ["Sd", "Da", "gamy", "Ca"], "Before", 1);

if ~isfolder(cfg.output_root)
    mkdir(cfg.output_root);
end
writetable(all_growth_rates, ...
    fullfile(cfg.output_root, "cylinder_growth_rates_vs_sd.csv"));
writetable(sweep_summary, ...
    fullfile(cfg.output_root, "cylinder_growth_sweep_summary.csv"));

%%% Select accepted fits and plot one curve per mode

if include_poor_fits
    accepted_status = ismember(all_growth_rates.fit_status, ["good", "poor"]);
else
    accepted_status = all_growth_rates.fit_status == "good";
end
plot_rows = accepted_status ...
    & isfinite(all_growth_rates.growth_rate_sigma) ...
    & isfinite(all_growth_rates.Sd);
if y_axis_scale == "log"
    plot_rows = plot_rows & all_growth_rates.growth_rate_sigma > 0;
end

if ~isempty(modes_to_plot)
    validateattributes(modes_to_plot, {'numeric'}, ...
        {'2d', 'ncols', 2, 'integer', 'finite'}, mfilename, "modes_to_plot");
    requested_mode = false(height(all_growth_rates), 1);
    for mode_index = 1:size(modes_to_plot, 1)
        requested_mode = requested_mode ...
            | (all_growth_rates.m == modes_to_plot(mode_index, 1) ...
            & all_growth_rates.n == modes_to_plot(mode_index, 2));
    end
    plot_rows = plot_rows & requested_mode;
end

plot_rates = all_growth_rates(plot_rows, :);
if isempty(plot_rates)
    error("CylinderMode:NoAcceptedGrowthRates", ...
        "No growth-rate fits satisfy the current status and mode filters.");
end

mode_pairs = unique([plot_rates.m, plot_rates.n], "rows", "sorted");
mode_point_count = zeros(size(mode_pairs, 1), 1);
for mode_index = 1:size(mode_pairs, 1)
    mode_point_count(mode_index) = nnz( ...
        plot_rates.m == mode_pairs(mode_index, 1) ...
        & plot_rates.n == mode_pairs(mode_index, 2));
end
mode_pairs = mode_pairs(mode_point_count >= minimum_Sd_points_per_mode, :);
if isempty(mode_pairs)
    error("CylinderMode:InsufficientSweepCoverage", ...
        "No mode has accepted fits at the required number of Sd values.");
end

if isempty(mode_color_order)
    plot_colors = lines(size(mode_pairs, 1));
else
    validateattributes(mode_color_order, {'numeric'}, ...
        {'2d', 'ncols', 3, '>=', 0, '<=', 1}, ...
        mfilename, "mode_color_order");
    repeats = ceil(size(mode_pairs, 1) / size(mode_color_order, 1));
    plot_colors = repmat(mode_color_order, repeats, 1);
    plot_colors = plot_colors(1:size(mode_pairs, 1), :);
end

fig = figure( ...
    "Color", figure_background_color, ...
    "Visible", figure_visible, ...
    "Position", figure_position);
ax = axes(fig);
hold(ax, "on");

if connect_points
    line_style = "-";
else
    line_style = "none";
end

for mode_index = 1:size(mode_pairs, 1)
    m = mode_pairs(mode_index, 1);
    n = mode_pairs(mode_index, 2);
    mode_rows = plot_rates.m == m & plot_rates.n == n;
    mode_rates = sortrows(plot_rates(mode_rows, :), "Sd");

    plot(ax, mode_rates.Sd, mode_rates.growth_rate_sigma, ...
        "LineStyle", line_style, ...
        "LineWidth", data_line_width, ...
        "Marker", marker_symbol, ...
        "MarkerSize", marker_size, ...
        "MarkerFaceColor", plot_colors(mode_index, :), ...
        "MarkerEdgeColor", plot_colors(mode_index, :), ...
        "Color", plot_colors(mode_index, :), ...
        "DisplayName", sprintf("m = %d, n = %d", m, n));
end

if show_zero_line && y_axis_scale == "linear"
    yline(ax, 0, ":", "Color", [0.25, 0.25, 0.25], ...
        "LineWidth", 1, "HandleVisibility", "off");
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
xlabel(ax, x_axis_label, "FontSize", label_font_size);
ylabel(ax, y_axis_label, "FontSize", label_font_size);
title(ax, sprintf("Cylindrical mode growth rates, Ca = %.4g", Ca), ...
    "FontSize", title_font_size, "FontWeight", "normal");

if show_legend
    legend_handle = legend(ax, ...
        "Location", legend_location, ...
        "FontSize", legend_font_size, ...
        "NumColumns", legend_num_columns, ...
        "Box", legend_box);
    legend_handle.AutoUpdate = "off";
end

figure_path = fullfile(cfg.output_root, ...
    "cylinder_mode_growth_rates_vs_sd.png");
if save_figure
    exportgraphics(fig, figure_path, "Resolution", figure_resolution);
end
if save_matlab_figure
    savefig(fig, fullfile(cfg.output_root, ...
        "cylinder_mode_growth_rates_vs_sd.fig"));
end
if save_combined_mat
    save(fullfile(cfg.output_root, ...
        "cylinder_mode_growth_rates_vs_sd.mat"), ...
        "all_growth_rates", "sweep_summary", "mode_pairs", ...
        "Sd_used", "gamy_used", "Ca", "Da", "v", "cfg");
end

fprintf("Analyzed %d fixed-Ca simulations and plotted %d modes.\n", ...
    height(sweep_summary), size(mode_pairs, 1));
fprintf("Ca = %.6g, with gamy = Ca*Sd.\n", Ca);
fprintf("Outputs: %s\n", cfg.output_root);

function tag = cylinder_parameter_tag(Sd, Da, gamy, v)
    tag = sprintf("Sd_%.2e_Da_%.2e_gamy_%+.2e_v_%.2e", ...
        Sd, Da, gamy, v);
    tag = strrep(tag, "+", "p");
    tag = strrep(tag, "-", "m");
end

function tag = cylinder_sweep_tag(Ca, Da, v)
    tag = sprintf("Ca_%.2e_Da_%.2e_v_%.2e", Ca, Da, v);
    tag = strrep(tag, "+", "p");
    tag = strrep(tag, "-", "m");
end
