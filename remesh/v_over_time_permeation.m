%V_OVER_TIME_PERMEATION Plot reduced volume for gamy/v parameter combinations.

clearvars;
close all;
clc;

%%% Simulation selection

Sd = 1e-6;                         % One fixed Saffman-Delbruck number.
Da = 1e-1;                         % One fixed permeability.
gamy_values = [1.01e-6,1.01e-5];       % Every gamy is combined with every v.
v_values = [0.5, 0.7, 0.9];       % Initial reduced-volume parameters.

time_start = 0;                    % First geoN.mat index to include.
time_final = inf;                  % Last geoN.mat index; inf uses all.
frame_stride = 1;                  % Analyze every Nth available frame.

%%% Visual settings

visual.figure_position = [100, 100, 1100, 700];
visual.figure_color = "w";
visual.axes_color = "w";

visual.x_scale = "log";         % "linear" or "log" time axis.
visual.x_limits = [];              % Empty uses automatic limits.
visual.y_limits = [];              % Empty uses automatic limits.
visual.prepend_initial_point = true; % Begin each plotted curve at the nominal v.
visual.initial_point_time = 1;   % Near-zero time used for that initial point.
visual.exponential_initial_extension = true; % Smoothly extend the first real point left to v.
visual.initial_extension_point_count = 60;
visual.initial_extension_decay = 5; % Number of exponential e-folds across the extension.

visual.line_width = 4.0;
visual.line_styles = ["-", ":", ":", "-."]; % Assigned by gamy, then cycled.
visual.marker = "none";
visual.marker_size = 5;

% Values of v are distinguished by different shades of blue.
visual.light_blue = [0.42, 0.70, 0.93];
visual.dark_blue = [0.02, 0.18, 0.48];

visual.axes_font_size = 13;
visual.label_font_size = 16;
visual.title_font_size = 16;
visual.legend_font_size = 16;
visual.legend_location = "northwest"; % Interior axes location.
visual.legend_num_columns = 1;
visual.legend_box = "on";

visual.show_grid = true;
visual.grid_alpha = 0.18;
visual.show_box = true;
visual.x_label = "Surface-scaled time, t = n \Delta t";
visual.y_label = "Reduced volume";
visual.title = sprintf("Reduced volume over time: Sd = %.4g, Da = %.4g", Sd, Da);

visual.save_figure = true;
visual.output_file = fullfile("data", "./figures/figure3.png");
visual.output_resolution = 300;
visual.save_summary = false;
visual.summary_file = fullfile("data", "v_over_time_permeation.csv");

%%% Locate data and validate settings

script_dir = fileparts(mfilename("fullpath"));
data_dir = find_data_dir(script_dir);

validateattributes(Sd, {'numeric'}, ...
    {'scalar', 'real', 'finite'}, mfilename, "Sd");
validateattributes(Da, {'numeric'}, ...
    {'scalar', 'real', 'finite'}, mfilename, "Da");
validateattributes(gamy_values, {'numeric'}, ...
    {'vector', 'real', 'finite', 'nonempty'}, mfilename, "gamy_values");
validateattributes(v_values, {'numeric'}, ...
    {'vector', 'real', 'finite', 'nonempty'}, mfilename, "v_values");
validateattributes(frame_stride, {'numeric'}, ...
    {'scalar', 'integer', 'positive'}, mfilename, "frame_stride");

if time_final < time_start
    error("time_final must be greater than or equal to time_start.");
end
if ~ismember(lower(string(visual.x_scale)), ["linear", "log"])
    error("visual.x_scale must be either 'linear' or 'log'.");
end
if Sd == 0
    error("Sd must be nonzero to calculate the capillary number chi = gamy/Sd.");
end
if isempty(visual.line_styles)
    error("visual.line_styles must contain at least one line style.");
end
validateattributes(visual.initial_point_time, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'}, ...
    mfilename, "visual.initial_point_time");
validateattributes(visual.initial_extension_point_count, {'numeric'}, ...
    {'scalar', 'integer', '>=', 2}, ...
    mfilename, "visual.initial_extension_point_count");
validateattributes(visual.initial_extension_decay, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'}, ...
    mfilename, "visual.initial_extension_decay");

gamy_values = gamy_values(:).';
v_values = v_values(:).';
v_colors = blue_shades(numel(v_values), visual.light_blue, visual.dark_blue);

%%% Load and plot every gamy/v combination

fig = figure( ...
    "Color", visual.figure_color, ...
    "Position", visual.figure_position);
ax = axes(fig, "Color", visual.axes_color);
hold(ax, "on");

summary = table('Size', [0, 8], ...
    'VariableTypes', repmat("double", 1, 8), ...
    'VariableNames', ["Sd", "Da", "gamy", "v_parameter", ...
        "n_frames", "initial_reduced_volume", ...
        "final_reduced_volume", "net_reduced_volume_change"]);

plotted_count = 0;
for gamy_index = 1:numel(gamy_values)
    gamy = gamy_values(gamy_index);
    style_index = mod(gamy_index - 1, numel(visual.line_styles)) + 1;
    line_style = visual.line_styles(style_index);

    for v_index = 1:numel(v_values)
        v_parameter = v_values(v_index);
        run_tag = make_run_tag(Sd, Da, gamy, v_parameter);
        folder = fullfile(data_dir, run_tag);

        if ~isfolder(folder)
            warning("Run folder not found: %s", folder);
            continue
        end

        fprintf("Loading %s\n", folder);
        try
            [time, reduced_volume_values, frame_ids] = load_reduced_volume_series( ...
                folder, time_start, time_final, frame_stride);
        catch exception
            warning("Could not process %s: %s", run_tag, exception.message);
            continue
        end

        plot_time = time;
        plot_reduced_volume = reduced_volume_values;
        if visual.prepend_initial_point
            keep_plot_point = plot_time > visual.initial_point_time;
            plot_time = plot_time(keep_plot_point);
            plot_reduced_volume = plot_reduced_volume(keep_plot_point);
            if isempty(plot_time)
                warning("No data occur after visual.initial_point_time for %s.", run_tag);
                continue
            end
            if visual.exponential_initial_extension
                [extension_time, extension_volume] = exponential_initial_extension( ...
                    visual.initial_point_time, v_parameter, ...
                    plot_time(1), plot_reduced_volume(1), ...
                    visual.x_scale, visual.initial_extension_point_count, ...
                    visual.initial_extension_decay);
                plot_time = [extension_time; plot_time(2:end)];
                plot_reduced_volume = [extension_volume; ...
                    plot_reduced_volume(2:end)];
            else
                plot_time = [visual.initial_point_time; plot_time];
                plot_reduced_volume = [v_parameter; plot_reduced_volume];
            end
        elseif lower(string(visual.x_scale)) == "log"
            keep_plot_point = plot_time > 0;
            plot_time = plot_time(keep_plot_point);
            plot_reduced_volume = plot_reduced_volume(keep_plot_point);
        end

        plot(ax, plot_time, plot_reduced_volume, ...
            "Color", v_colors(v_index, :), ...
            "LineStyle", line_style, ...
            "LineWidth", visual.line_width, ...
            "Marker", visual.marker, ...
            "MarkerSize", visual.marker_size, ...
            "HandleVisibility", "off");

        plotted_count = plotted_count + 1;
        summary(end + 1, :) = { ...
            Sd, Da, gamy, v_parameter, numel(frame_ids), ...
            reduced_volume_values(1), reduced_volume_values(end), ...
            reduced_volume_values(end) - reduced_volume_values(1)}; %#ok<SAGROW>

        fprintf("  gamy = %.8g, v = %.8g: %.9g -> %.9g over %d frames\n", ...
            gamy, v_parameter, reduced_volume_values(1), ...
            reduced_volume_values(end), numel(frame_ids));
    end
end

if plotted_count == 0
    close(fig);
    error("No requested gamy/v simulation combinations could be plotted.");
end

%%% Finish the figure

set(ax, ...
    "XScale", visual.x_scale, ...
    "FontSize", visual.axes_font_size, ...
    "GridAlpha", visual.grid_alpha, ...
    "MinorGridAlpha", visual.grid_alpha);
xlabel(ax, visual.x_label, "FontSize", visual.label_font_size);
ylabel(ax, visual.y_label, "FontSize", visual.label_font_size);
title(ax, visual.title, ...
    "FontSize", visual.title_font_size, ...
    "FontWeight", "normal");

if visual.show_grid
    grid(ax, "on");
    ax.XMinorGrid = "on";
    ax.YMinorGrid = "on";
end
if visual.show_box
    box(ax, "on");
else
    box(ax, "off");
end
if ~isempty(visual.x_limits)
    xlim(ax, visual.x_limits);
end
if ~isempty(visual.y_limits)
    ylim(ax, visual.y_limits);
end


xlim([1e1,4e6])
ylim([.4,1])

legend_handles = gobjects(numel(gamy_values), 1);
for gamy_index = 1:numel(gamy_values)
    style_index = mod(gamy_index - 1, numel(visual.line_styles)) + 1;
    chi = gamy_values(gamy_index) / Sd;
    display_name = sprintf("%cchi = %.0f", char(92), chi);
    legend_handles(gamy_index) = plot(ax, NaN, NaN, ...
        "Color", visual.dark_blue, ...
        "LineStyle", visual.line_styles(style_index), ...
        "LineWidth", visual.line_width, ...
        "DisplayName", display_name);
end

legend_handle = legend(ax, legend_handles, ...
    "Location", visual.legend_location, ...
    "FontSize", visual.legend_font_size, ...
    "NumColumns", visual.legend_num_columns, ...
    "Box", visual.legend_box);
legend_handle.AutoUpdate = "off";

if visual.save_figure
    ensure_parent_folder(visual.output_file);
    exportgraphics(fig, visual.output_file, ...
        "Resolution", visual.output_resolution);
    fprintf("Saved figure: %s\n", visual.output_file);
end
if visual.save_summary
    ensure_parent_folder(visual.summary_file);
    writetable(summary, visual.summary_file);
    fprintf("Saved summary: %s\n", visual.summary_file);
end

function [time, reduced_volume_values, frame_ids] = ...
        load_reduced_volume_series(folder, time_start, time_final, frame_stride)
    files = dir(fullfile(folder, "geo*.mat"));
    frame_ids = frame_ids_from_files(files);
    valid = isfinite(frame_ids);
    files = files(valid);
    frame_ids = frame_ids(valid);

    [frame_ids, order] = sort(frame_ids);
    files = files(order);
    keep = frame_ids >= time_start & frame_ids <= time_final;
    frame_ids = frame_ids(keep);
    files = files(keep);
    frame_ids = frame_ids(1:frame_stride:end);
    files = files(1:frame_stride:end);

    if isempty(files)
        error("No geoN.mat frames fall within the requested frame range.");
    end

    frame_count = numel(files);
    reduced_volume_values = nan(frame_count, 1);
    dt = nan(frame_count, 1);

    for frame_index = 1:frame_count
        data = load(fullfile(files(frame_index).folder, files(frame_index).name), ...
            "M", "P", "p");
        reduced_volume_values(frame_index) = mesh_reduced_volume(data.M, data.P);
        if isfield(data, "p") && isstruct(data.p) ...
                && isfield(data.p, "dt") && isscalar(data.p.dt)
            dt(frame_index) = double(data.p.dt);
        end
    end

    if any(~isfinite(dt))
        error("At least one selected frame does not contain a finite p.dt.");
    end

    time = nan(frame_count, 1);
    time(1) = frame_ids(1) * dt(1);
    for frame_index = 2:frame_count
        frame_step = frame_ids(frame_index) - frame_ids(frame_index - 1);
        time(frame_index) = time(frame_index - 1) + frame_step * dt(frame_index);
    end
end

function value = mesh_reduced_volume(M, P)
    M = double(M);
    P = double(P);
    p1 = P(M(:, 1), :);
    p2 = P(M(:, 2), :);
    p3 = P(M(:, 3), :);

    face_cross = cross(p2 - p1, p3 - p1, 2);
    area = 0.5 * sum(vecnorm(face_cross, 2, 2));
    signed_volume = sum(dot(p1, cross(p2, p3, 2), 2)) / 6;
    value = 6 * sqrt(pi) * signed_volume / area^(3 / 2);
end

function frame_ids = frame_ids_from_files(files)
    frame_ids = nan(numel(files), 1);
    for file_index = 1:numel(files)
        token = regexp(files(file_index).name, ...
            "^geo(\d+)\.mat$", "tokens", "once");
        if ~isempty(token)
            frame_ids(file_index) = str2double(token{1});
        end
    end
end

function colors = blue_shades(count, light_blue, dark_blue)
    if count == 1
        colors = 0.5 * (light_blue + dark_blue);
        return
    end
    fraction = linspace(0, 1, count).';
    colors = light_blue .* (1 - fraction) + dark_blue .* fraction;
end

function [extension_time, extension_volume] = exponential_initial_extension( ...
        initial_time, initial_volume, first_time, first_volume, ...
        x_scale, point_count, decay)
    if first_time <= initial_time
        error("The first real time must exceed visual.initial_point_time.");
    end

    s = linspace(0, 1, point_count).';
    if lower(string(x_scale)) == "log"
        extension_time = logspace( ...
            log10(initial_time), log10(first_time), point_count).';
    else
        extension_time = initial_time + s * (first_time - initial_time);
    end

    shape = expm1(decay * s) / expm1(decay);
    extension_volume = initial_volume ...
        + (first_volume - initial_volume) * shape;
end

function run_tag = make_run_tag(Sd, Da, gamy, v)
    run_tag = sprintf( ...
        "Sd_%.2e_Da_%.2e_gamy_%+.2e_v_%.2e", Sd, Da, gamy, v);
    run_tag = strrep(run_tag, "+", "p");
    run_tag = strrep(run_tag, "-", "m");
end

function data_dir = find_data_dir(script_dir)
    candidates = string({ ...
        fullfile(script_dir, "data", "fs_batch_data"), ...
        fullfile(pwd, "data", "fs_batch_data"), ...
        fullfile(pwd, "remesh", "data", "fs_batch_data")});

    for candidate_index = 1:numel(candidates)
        if isfolder(candidates(candidate_index))
            data_dir = char(candidates(candidate_index));
            return
        end
    end
    error("Could not locate data/fs_batch_data. Checked:%s", ...
        sprintf("\n  %s", candidates));
end

function ensure_parent_folder(filename)
    parent = fileparts(filename);
    if strlength(parent) > 0 && ~isfolder(parent)
        mkdir(parent);
    end
end
