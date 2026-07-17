close all;clear;clc;
%%% Simulation selection
% Sd value used in the fs_batch output folder names.
Sd = [1e2]%,1e3,1.01e4];

% Background shear/extension parameter used in the fs_batch output folder names.
gamy = 0;

% Darcy numbers to load and plot on the same axes.
%Da_values = [1e-12,1e-11,1e-10,1e-9, 1e-8, 4e-8,6e-8, 1e-7, 2e-7,4e-7,6e-7,8e-7];%, 1e-6,2e-6,4e-6,6e-6,8e-6,1e-5, 1e-4, 1e-2];
Da_values = 1

% Time-step range. If true, use every available geo*.mat frame for each run.
% If false, use only geo0.mat through geo<maxtimestep>.mat.
usealltimes = true;
maxtimestep = 10;

% Rate denominator. "time" divides by physical time; "timestep" divides by geo frame index.
rate_time_basis = "timestep";

Sd_values = Sd(:).';
multiple_sd = numel(Sd_values) > 1;
if multiple_sd
    sd_title_text = sprintf("%d Sd values", numel(Sd_values));
else
    sd_title_text = sprintf("Sd = %.4g", Sd_values);
end
[rate_time_basis, rate_denominator_label, rate_legend_label] = rate_basis_labels(rate_time_basis);

%%% Visual settings
% First figure: relative volume change over time.
visual.show_time_plot = false;

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
visual.title = sprintf("Relative volume change over time, %s", sd_title_text);

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

% Second figure: log-log scatter of Da against a volume-change rate.
visual.show_rate_scatter = true;
visual.rate_figure_position = [140, 140, 780, 620];
visual.rate_source = "mean_relative_dVdt";   % Options: "mean_relative_dVdt", "net_relative_dVdt", "mean_dVdt", "net_dVdt".
visual.rate_use_absolute_value = false;      % true plots abs(rate), useful if rates are negative.
visual.rate_marker_size = 70;
visual.rate_marker = "o";
visual.rate_marker_face_color = [0.10, 0.35, 0.80];
visual.rate_marker_edge_color = "k";
visual.rate_line_width = 1.0;
visual.rate_x_label = "Da";
visual.rate_y_label = sprintf("(V-V_0/V_0)/%s", rate_denominator_label);
visual.rate_title = sprintf("Volume-change rate vs permeability, %s", ...
     sd_title_text);
visual.rate_show_grid = true;
visual.rate_show_box = true;
visual.rate_x_limits = [];
visual.rate_y_limits = [];
visual.rate_save_figure = false;
visual.rate_output_file = fullfile("data", "permeation_rate_scatter.png");
visual.rate_sd_color_map = "lines";         % Used to color each Sd group when multiple Sd values are given.
visual.rate_sd_color_order = [];            % Optional nSd x 3 RGB colors. Leave [] to use visual.rate_sd_color_map.

% Power-law fit on the rate scatter: Vdot = prefactor * Da^alpha.
% Fit indices refer to the sorted valid positive Da values used in the scatter.
visual.rate_show_fit = true;
visual.rate_fit_first_index = 4;             % First sorted valid Da point included in the fit.
visual.rate_fit_last_index = 10;              % Last sorted valid Da point included in the fit.
visual.rate_fit_ranges_by_sd = [4,9;4,8;4,7];           % Optional nSd x 2 matrix: row i is [first, last] for Sd_values(i).
visual.rate_fit_line_color = [0.85, 0.10, 0.10];
visual.rate_fit_line_style = "-";
visual.rate_fit_line_width = 2.5;
visual.rate_fit_n_plot_points = 150;
visual.rate_show_legend = true;
visual.rate_legend_location = "best";

script_dir = fileparts(mfilename("fullpath"));
remesh_dir = find_remesh_dir(script_dir);
addpath(remesh_dir);
data_dir = find_data_dir(script_dir, remesh_dir);

if visual.show_time_plot
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
end

colors = visual.color_order;
if isempty(colors)
    colors = feval(char(visual.color_map), max(numel(Da_values), 1));
end
sd_colors = visual.rate_sd_color_order;
if isempty(sd_colors)
    sd_colors = feval(char(visual.rate_sd_color_map), max(numel(Sd_values), 1));
end

summary = table('Size', [0, 9], ...
    'VariableTypes', ["double", "double", "double", "double", "double", "double", "double", "double", "double"], ...
    'VariableNames', ["Sd", "Da", "initial_volume", "final_volume", "final_relative_volume", ...
        "net_dVdt", "mean_dVdt", "net_relative_dVdt", "mean_relative_dVdt"]);

for sd_idx = 1:numel(Sd_values)
    Sd_current = Sd_values(sd_idx);
    for k = 1:numel(Da_values)
    Da = Da_values(k);
    folder = fullfile(data_dir, make_run_tag(Sd_current, Da, gamy));
    frames = list_geo_frames(folder);
    [frames, time_range] = apply_time_range(frames, usealltimes, maxtimestep);

    if isempty(frames)
        warning("%s", no_frames_message(Da, folder, time_range));
        continue
    end
    if time_range.is_capped && frames(end) < time_range.maxtimestep
        warning("Sd = %.4g, Da = %.4g only reaches geo%d.mat before maxtimestep = %d; rates use the shorter available interval.", ...
            Sd_current, Da, frames(end), time_range.maxtimestep);
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
    frame_ids = frames(order);
    volume0 = volumes(1);
    relative_volumes = (volumes - volume0) / volume0;

    if numel(times) > 1
        [rate_denominator, net_denominator] = rate_denominators(times, frame_ids, rate_time_basis);
        dVdt = diff(volumes) ./ rate_denominator;
        relative_dVdt = diff(relative_volumes) ./ rate_denominator;
        mean_dVdt = mean(dVdt, "omitnan");
        mean_relative_dVdt = mean(relative_dVdt, "omitnan");
        net_dVdt = (volumes(end) - volumes(1)) / net_denominator;
        net_relative_dVdt = (relative_volumes(end) - relative_volumes(1)) / net_denominator;
    else
        mean_dVdt = NaN;
        mean_relative_dVdt = NaN;
        net_dVdt = NaN;
        net_relative_dVdt = NaN;
    end

    if visual.show_mean_dVdt_in_legend
        if multiple_sd
            display_name = sprintf("Sd = %.3g, Da = %.3g, %s = %.3e", ...
                Sd_current, Da, rate_legend_label, mean_relative_dVdt);
        else
            display_name = sprintf("Da = %.3g, %s = %.3e", Da, rate_legend_label, mean_relative_dVdt);
        end
    else
        if multiple_sd
            display_name = sprintf("Sd = %.3g, Da = %.3g", Sd_current, Da);
        else
            display_name = sprintf("Da = %.3g", Da);
        end
    end

    plot_values = relative_volumes;
    if visual.y_axis_log && any(plot_values <= 0)
        warning("Log y-axis requested: omitting nonpositive relative volume values for Da = %.4g.", Da);
        plot_values(plot_values <= 0) = NaN;
    end

    if multiple_sd
        color_idx = mod(sd_idx - 1, size(sd_colors, 1)) + 1;
        line_color = sd_colors(color_idx, :);
    else
        color_idx = mod(k - 1, size(colors, 1)) + 1;
        line_color = colors(color_idx, :);
    end
    if visual.show_time_plot
        style_idx = mod(k - 1, numel(visual.line_styles)) + 1;
        plot(times, plot_values, ...
            "Color", line_color, ...
            "LineStyle", visual.line_styles(style_idx), ...
            "LineWidth", visual.line_width, ...
            "Marker", visual.marker, ...
            "MarkerSize", visual.marker_size, ...
            "DisplayName", display_name);
    end

    summary(end + 1, :) = {Sd_current, Da, volumes(1), volumes(end), relative_volumes(end), ...
        net_dVdt, mean_dVdt, net_relative_dVdt, mean_relative_dVdt}; %#ok<SAGROW>
    end
end

if visual.show_time_plot
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
end

if visual.show_time_plot && visual.save_figure
    exportgraphics(gcf, visual.output_file, "Resolution", visual.output_resolution);
end

if visual.show_rate_scatter
    if ~ismember(visual.rate_source, string(summary.Properties.VariableNames))
        error("Unknown visual.rate_source '%s'. Use one of: %s", ...
            visual.rate_source, strjoin(string(summary.Properties.VariableNames), ", "));
    end
    rate_values = summary.(visual.rate_source);
    if visual.rate_use_absolute_value
        rate_values = abs(rate_values);
    end

    valid_rate_all = summary.Da > 0 & rate_values > 0 & isfinite(rate_values);
    if any(~valid_rate_all)
        warning("Omitting %d point(s) from rate scatter because log-log axes require Da > 0 and rate > 0.", ...
            nnz(~valid_rate_all));
    end

    figure("Color", visual.figure_color, "Position", visual.rate_figure_position);
    rate_ax = axes("Color", visual.axes_color);
    hold on;
    if visual.rate_show_grid
        grid on;
    else
        grid off;
    end
    if visual.rate_show_box
        box on;
    else
        box off;
    end

    for sd_idx = 1:numel(Sd_values)
        Sd_current = Sd_values(sd_idx);
        sd_group = summary.Sd == Sd_current;
        valid_rate = valid_rate_all & sd_group;
        if ~any(valid_rate)
            warning("Skipping rate scatter for Sd = %.4g because no positive log-log points are available.", ...
                Sd_current);
            continue
        end

        if multiple_sd
            sd_color_idx = mod(sd_idx - 1, size(sd_colors, 1)) + 1;
            scatter_color = sd_colors(sd_color_idx, :);
            fit_color = scatter_color;
            scatter_name = sprintf("Sd = %.1g data", Sd_current);
        else
            scatter_color = visual.rate_marker_face_color;
            fit_color = visual.rate_fit_line_color;
            scatter_name = "data";
        end

        scatter(summary.Da(valid_rate), rate_values(valid_rate), ...
            visual.rate_marker_size, ...
            "Marker", visual.rate_marker, ...
            "MarkerFaceColor", scatter_color, ...
            "MarkerEdgeColor", visual.rate_marker_edge_color, ...
            "LineWidth", visual.rate_line_width, ...
            "DisplayName", scatter_name);

        if visual.rate_show_fit
            fit_Da_all = summary.Da(valid_rate);
            fit_rate_all = rate_values(valid_rate);
            [fit_Da_all, fit_order] = sort(fit_Da_all);
            fit_rate_all = fit_rate_all(fit_order);
            [first_fit, last_fit] = fit_index_range_for_sd(visual, sd_idx, numel(fit_Da_all));
            n_fit = last_fit - first_fit + 1;

            if n_fit >= 2
                fit_Da = fit_Da_all(first_fit:last_fit);
                fit_rate = fit_rate_all(first_fit:last_fit);
                fit_coeff = polyfit(log(fit_Da), log(fit_rate), 1);
                alpha = fit_coeff(1);
                prefactor = exp(fit_coeff(2));

                log_fit_rate = polyval(fit_coeff, log(fit_Da));
                ss_res = sum((log(fit_rate) - log_fit_rate) .^ 2);
                ss_tot = sum((log(fit_rate) - mean(log(fit_rate))) .^ 2);
                r_squared = 1 - ss_res / ss_tot;

                x_fit = logspace(log10(min(fit_Da)), log10(max(fit_Da)), ...
                    visual.rate_fit_n_plot_points);
                y_fit = prefactor * x_fit .^ alpha;
                if multiple_sd
                    fit_name = sprintf("Sd = %.1g fit: \\alpha = %.3g", Sd_current, alpha);
                else
                    fit_name = sprintf("fit: \\alpha = %.3g", alpha);
                end
                plot(x_fit, y_fit, ...
                    "Color", fit_color, ...
                    "LineStyle", visual.rate_fit_line_style, ...
                    "LineWidth", visual.rate_fit_line_width, ...
                    "DisplayName", fit_name);

                if multiple_sd
                    fprintf("Sd = %.6g power-law fit using sorted Da indices %d:%d: Vdot = %.6g * Da^%.6g, R^2 = %.6g\n", ...
                        Sd_current, first_fit, last_fit, prefactor, alpha, r_squared);
                else
                    fprintf("Power-law fit using sorted Da indices %d:%d: Vdot = %.6g * Da^%.6g, R^2 = %.6g\n", ...
                        first_fit, last_fit, prefactor, alpha, r_squared);
                end
            else
                warning("Skipping power-law fit for Sd = %.4g because the selected sorted Da index range contains fewer than 2 points.", ...
                    Sd_current);
            end
        end
    end

    set(rate_ax, ...
        "FontSize", visual.axes_font_size, ...
        "XScale", "log", ...
        "YScale", "log");
    xlabel(visual.rate_x_label, "FontSize", visual.label_font_size);
    ylabel(visual.rate_y_label, "FontSize", visual.label_font_size);
    title(visual.rate_title, "FontSize", visual.title_font_size);
    if ~isempty(visual.rate_x_limits)
        xlim(visual.rate_x_limits);
    end
    if ~isempty(visual.rate_y_limits)
        ylim(visual.rate_y_limits);
    end
    if visual.rate_show_legend
        legend("Location", visual.rate_legend_location, "FontSize", visual.legend_font_size);
    end

    if visual.rate_save_figure
        exportgraphics(gcf, visual.rate_output_file, "Resolution", visual.output_resolution);
    end
end

if multiple_sd
    disp(summary);
else
    disp(removevars(summary, "Sd"));
end

function [rate_time_basis, denominator_label, legend_label] = rate_basis_labels(rate_time_basis)
    rate_time_basis = lower(string(rate_time_basis));
    if rate_time_basis == "time"
        denominator_label = "dt";
        legend_label = "mean dVrel/dt";
    elseif rate_time_basis == "timestep" || rate_time_basis == "step" || rate_time_basis == "frame"
        rate_time_basis = "timestep";
        denominator_label = "timestep";
        legend_label = "mean dVrel/step";
    else
        error('Unknown rate_time_basis "%s". Use "time" or "timestep".', rate_time_basis);
    end
end

function [rate_denominator, net_denominator] = rate_denominators(times, frame_ids, rate_time_basis)
    if rate_time_basis == "time"
        rate_denominator = diff(times);
        net_denominator = times(end) - times(1);
    else
        rate_denominator = diff(frame_ids);
        net_denominator = frame_ids(end) - frame_ids(1);
    end

    if any(rate_denominator == 0) || net_denominator == 0
        error("Rate denominator contains zero. Check for duplicate frame indices or times.");
    end
end

function [frames, time_range] = apply_time_range(frames, usealltimes, maxtimestep)
    time_range.is_capped = ~usealltimes;
    time_range.maxtimestep = maxtimestep;

    if time_range.is_capped
        frames = frames(frames <= maxtimestep);
    end
end

function [first_fit, last_fit] = fit_index_range_for_sd(visual, sd_idx, n_points)
    if isempty(visual.rate_fit_ranges_by_sd)
        first_fit = visual.rate_fit_first_index;
        last_fit = visual.rate_fit_last_index;
    else
        if size(visual.rate_fit_ranges_by_sd, 2) ~= 2
            error("visual.rate_fit_ranges_by_sd must be an nSd x 2 matrix of [first, last] fit indices.");
        end
        if sd_idx > size(visual.rate_fit_ranges_by_sd, 1)
            error("visual.rate_fit_ranges_by_sd has fewer rows than the number of Sd values.");
        end
        first_fit = visual.rate_fit_ranges_by_sd(sd_idx, 1);
        last_fit = visual.rate_fit_ranges_by_sd(sd_idx, 2);
    end

    first_fit = max(1, first_fit);
    last_fit = min(last_fit, n_points);
end

function message = no_frames_message(Da, folder, time_range)
    if time_range.is_capped
        message = sprintf("No geo*.mat files found at or before maxtimestep = %d for Da = %.4g in %s", ...
            time_range.maxtimestep, Da, folder);
    else
        message = sprintf("No geo*.mat files found for Da = %.4g in %s", Da, folder);
    end
end

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
