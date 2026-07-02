clear; clc; close all;

%%% Simulation selection
% Sd value used in the fs_batch output folder names.
Sd = 1200;

% Background shear/extension parameter used in the fs_batch output folder names.
gamy = 0;

% Darcy numbers to load and plot on the same axes.
Da_values = [0, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-2];

%%% Visual settings
% Figure position in pixels: [left, bottom, width, height].
visual.figure_position = [100, 100, 1450, 720];

% Figure and axes background colors.
visual.figure_color = "w";
visual.axes_color = "w";

% Color map used to assign one line color per Da. Examples: "lines", "turbo", "parula".
visual.color_map = "lines";

% Optional explicit color order as an n x 3 RGB matrix. Leave [] to use visual.color_map.
visual.color_order = [];

% Line appearance.
visual.line_width = 2.0;
visual.line_styles = ["-", "--", ":", "-."]; % Cycled through if there are more Da values than styles.
visual.marker = "none";                      % Examples: "none", "o", "s", ".", "^".
visual.marker_size = 6;

% Axes text sizes.
visual.axes_font_size = 13;
visual.label_font_size = 15;
visual.title_font_size = 16;
visual.legend_font_size = 11;

% Axes labels and title.
visual.x_label = "Time";
visual.y_label = "(V - V0) / V0";
visual.title = sprintf("Relative volume change over time, Sd = %.4g, gamy = %.4g", Sd, gamy);

% Axes scaling and limits. Use [] for automatic limits.
visual.x_scale = "linear";                   % "linear" or "log".
visual.y_axis_log = false;                   % true gives a log y-axis; nonpositive values are omitted.
visual.x_limits = [];
visual.y_limits = [];

% Grid, box, and legend controls.
visual.show_grid = true;
visual.show_box = true;
visual.legend_location = "best";
visual.show_mean_dVdt_in_legend = true;

% Optional output file. Set save_figure = true to export after plotting.
visual.save_figure = false;
visual.output_file = fullfile("data", "permeation_volume_plot.png");
visual.output_resolution = 300;

script_dir = fileparts(mfilename("fullpath"));
remesh_dir = find_remesh_dir(script_dir);
addpath(remesh_dir);
data_dir = find_data_dir(script_dir, remesh_dir);

figure("Color", visual.figure_color, "Position", visual.figure_position);
ax = axes("Color", visual.axes_color);
hold on;
if visual.show_grid
    grid on;
else
    grid off;
end
if visual.show_box
    box on;
else
    box off;
end

colors = visual.color_order;
if isempty(colors)
    colors = feval(char(visual.color_map), max(numel(Da_values), 1));
end

summary = table('Size', [0, 8], ...
    'VariableTypes', ["double", "double", "double", "double", "double", "double", "double", "double"], ...
    'VariableNames', ["Da", "initial_volume", "final_volume", "final_relative_volume", ...
        "net_dVdt", "mean_dVdt", "net_relative_dVdt", "mean_relative_dVdt"]);

for k = 1:numel(Da_values)
    Da = Da_values(k);
    folder = fullfile(data_dir, make_run_tag(Sd, Da, gamy));
    frames = list_geo_frames(folder);

    if isempty(frames)
        warning("No geo*.mat files found for Da = %.4g in %s", Da, folder);
        continue
    end

    times = zeros(numel(frames), 1);
    volumes = zeros(numel(frames), 1);

    for j = 1:numel(frames)
        frame_id = frames(j);
        data = load(fullfile(folder, sprintf("geo%d.mat", frame_id)));
        geo = Geometry(data.M, data.P);

        if isfield(data, "p") && isfield(data.p, "dt")
            dt = data.p.dt;
        else
            dt = 1;
        end

        times(j) = frame_id * dt;
        volumes(j) = geo.volume;
    end

    [times, order] = sort(times);
    volumes = volumes(order);
    volume0 = volumes(1);
    relative_volumes = (volumes - volume0) / volume0;

    if numel(times) > 1
        dVdt = diff(volumes) ./ diff(times);
        relative_dVdt = diff(relative_volumes) ./ diff(times);
        mean_dVdt = mean(dVdt, "omitnan");
        mean_relative_dVdt = mean(relative_dVdt, "omitnan");
        net_dVdt = (volumes(end) - volumes(1)) / (times(end) - times(1));
        net_relative_dVdt = (relative_volumes(end) - relative_volumes(1)) / (times(end) - times(1));
    else
        mean_dVdt = NaN;
        mean_relative_dVdt = NaN;
        net_dVdt = NaN;
        net_relative_dVdt = NaN;
    end

    if visual.show_mean_dVdt_in_legend
        display_name = sprintf("Da = %.3g, mean dVrel/dt = %.3e", Da, mean_relative_dVdt);
    else
        display_name = sprintf("Da = %.3g", Da);
    end

    plot_values = relative_volumes;
    if visual.y_axis_log && any(plot_values <= 0)
        warning("Log y-axis requested: omitting nonpositive relative volume values for Da = %.4g.", Da);
        plot_values(plot_values <= 0) = NaN;
    end

    color_idx = mod(k - 1, size(colors, 1)) + 1;
    style_idx = mod(k - 1, numel(visual.line_styles)) + 1;
    plot(times, plot_values, ...
        "Color", colors(color_idx, :), ...
        "LineStyle", visual.line_styles(style_idx), ...
        "LineWidth", visual.line_width, ...
        "Marker", visual.marker, ...
        "MarkerSize", visual.marker_size, ...
        "DisplayName", display_name);

    summary(end + 1, :) = {Da, volumes(1), volumes(end), relative_volumes(end), ...
        net_dVdt, mean_dVdt, net_relative_dVdt, mean_relative_dVdt}; %#ok<SAGROW>
end

xlabel(visual.x_label, "FontSize", visual.label_font_size);
ylabel(visual.y_label, "FontSize", visual.label_font_size);
title(visual.title, "FontSize", visual.title_font_size);
y_scale = "linear";
if visual.y_axis_log
    y_scale = "log";
end
set(ax, ...
    "FontSize", visual.axes_font_size, ...
    "XScale", visual.x_scale, ...
    "YScale", y_scale);
if ~isempty(visual.x_limits)
    xlim(visual.x_limits);
end
if ~isempty(visual.y_limits)
    ylim(visual.y_limits);
end
legend("Location", visual.legend_location, "FontSize", visual.legend_font_size);

if visual.save_figure
    exportgraphics(gcf, visual.output_file, "Resolution", visual.output_resolution);
end

disp(summary);

function run_tag = make_run_tag(Sd, Da, gamy)
    run_tag = sprintf('Sd_%.2e_Da_%.2e_gamy_%+.2e', Sd, Da, gamy);
    run_tag = strrep(run_tag, '+', 'p');
    run_tag = strrep(run_tag, '-', 'm');
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

function frames = list_geo_frames(folder)
    files = dir(fullfile(folder, "geo*.mat"));
    frames = zeros(numel(files), 1);
    keep = false(numel(files), 1);

    for i = 1:numel(files)
        token = regexp(files(i).name, '^geo(\d+)\.mat$', 'tokens', 'once');
        if ~isempty(token)
            frames(i) = str2double(token{1});
            keep(i) = true;
        end
    end

    frames = sort(frames(keep));
end
