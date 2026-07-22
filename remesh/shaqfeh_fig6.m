clear; close; clc;

%%% Parameters
figure_position = [600, 200, 1000, 700];
figure_color = 'w';

% Current-work simulation data to extract.
broussinos_Sd = 1e-6;
broussinos_Da = 0;
broussinos_gamy_values = [8e-7,8e-6,8e-5]; % Scalar or array, e.g. [4e-7, 8e-7, 1.2e-6].
broussinos_v_values = [0.988940, 0.97734, 0.96597, 0.954817, 0.9438, ...
    0.93314, 0.9226, 0.91227, 0.90213, 0.89217, 0.88240, 0.87280];
broussinos_tf = "last"; % "last" or a numeric frame id, e.g. 200.

% Axis data options.
x_axis_mode = "Lam"; % Options: "excess_area" or "Lam".
tilt_scaled_by_pi = true; % true plots psi/pi; false plots psi in radians.

x_limits = [0.05, 1.5]; % Excess-area limits; auto-converted if x_axis_mode = "Lam".
y_limits = [0.09, 0.25]; % psi/pi limits; auto-converted if tilt_scaled_by_pi = false.
axis_line_width = 1.25;
axis_font_size = 18;
axis_font_name = 'Helvetica';
axis_font_weight = 'normal';

curve_line_width = 1.5;
scatter_marker_size = 100;
scatter_line_width = 2;

shaqfeh_color = 'k';
misbah_color = 'k';
kantsler_edge_color = 'k';
broussinos_edge_color = 'r';
broussinos_color_map = "lines";
broussinos_color_order = [.8 .1 .1;.1 .8 .1; .1 .1 .8]; % Optional nGamy x 3 RGB matrix. Leave [] to use broussinos_color_map.

title_text = "Membrane shape vs tilt angle, \chi = 8";
title_font_size = 30;
title_font_weight = 'normal';
x_label_text = "";
y_label_text = "";
label_font_size = 22;
label_font_weight = 'normal';

legend_font_size = 20;
legend_location = 'northeast';

shaqfeh_data = [0.09142212189616253, 0.21728650137741048
0.1574492099322799, 0.20695592286501377
0.21839729119638826, 0.1993801652892562
0.3081264108352144, 0.19042699724517906
0.3876975169300226, 0.18353994490358128
0.48250564334085777, 0.17665289256198347
0.577313769751693, 0.17079889807162535
0.6992099322799097, 0.16391184573002754
0.8278781038374717, 0.15736914600550964
0.9633182844243792, 0.15151515151515152
1.1241534988713318, 0.14462809917355374
1.327313769751693, 0.13705234159779614];

kantsler_data=[0.12358916478555304, 0.22176308539944906
0.18623024830699775, 0.20695592286501377
0.25225733634311515, 0.19765840220385675
0.33860045146726864, 0.1825068870523416
0.43510158013544015, 0.17389807162534435
0.5468397291196389, 0.16460055096418733
0.6670428893905191, 0.1553030303030303
0.7990970654627539, 0.14738292011019283
0.9446952595936794, 0.1391184573002755
1.0987584650112867, 0.13292011019283748
1.2714446952595937, 0.12637741046831957];

misbah_data=[0.08803611738148984, 0.21694214876033058
0.1879232505643341, 0.20075757575757577
0.27257336343115124, 0.19008264462809918
0.4621896162528217, 0.17011019283746556
0.6992099322799097, 0.14910468319559228
0.9260722347629796, 0.13085399449035812
1.1580135440180586, 0.11225895316804407
1.3205417607223475, 0.09951790633608816];

[broussinos_data, broussinos_gamy_by_point] = extract_broussinos_data( ...
    broussinos_Sd, broussinos_Da, broussinos_gamy_values, ...
    broussinos_v_values, broussinos_tf);

[shaqfeh_data, x_label_text, y_label_text, x_limits, y_limits] = apply_axis_options( ...
    shaqfeh_data, x_axis_mode, tilt_scaled_by_pi, x_label_text, y_label_text, ...
    x_limits, y_limits);
[kantsler_data, ~, ~, ~, ~] = apply_axis_options( ...
    kantsler_data, x_axis_mode, tilt_scaled_by_pi, x_label_text, y_label_text, ...
    x_limits, y_limits);
[misbah_data, ~, ~, ~, ~] = apply_axis_options( ...
    misbah_data, x_axis_mode, tilt_scaled_by_pi, x_label_text, y_label_text, ...
    x_limits, y_limits);
[broussinos_data, ~, ~, ~, ~] = apply_axis_options( ...
    broussinos_data, x_axis_mode, tilt_scaled_by_pi, x_label_text, y_label_text, ...
    x_limits, y_limits);

fig = figure("Color", figure_color, "Position", figure_position);
ax = axes(fig);
hold(ax, "on");

plot(ax, shaqfeh_data(:,1), shaqfeh_data(:,2), ...
    Color=shaqfeh_color, LineWidth=curve_line_width, ...
    DisplayName="Zhao & Shaqfeh simulations");
plot(ax, misbah_data(:,1), misbah_data(:,2), ...
    Color=misbah_color, LineStyle='--', LineWidth=curve_line_width, ...
    DisplayName="Misbah perturbation theory");
scatter(ax, kantsler_data(:,1), kantsler_data(:,2), scatter_marker_size, ...
    LineWidth=scatter_line_width, MarkerEdgeColor=kantsler_edge_color, ...
    DisplayName="Kantsler & Steinberg measurements");
plot_broussinos_scatter(ax, broussinos_data, broussinos_gamy_by_point, ...
    broussinos_gamy_values, broussinos_color_order, broussinos_color_map, ...
    broussinos_edge_color, scatter_marker_size, scatter_line_width);

title(ax, title_text, Interpreter='tex', FontSize=title_font_size, ...
    FontWeight=title_font_weight);
xlabel(ax, x_label_text, FontSize=label_font_size, FontWeight=label_font_weight);
ylabel(ax, y_label_text, Interpreter='tex', FontSize=label_font_size, ...
    FontWeight=label_font_weight);
xlim(ax, x_limits);
ylim(ax, y_limits);
set(ax, LineWidth=axis_line_width, FontSize=axis_font_size, ...
    FontName=axis_font_name, FontWeight=axis_font_weight);
box on;
legend(ax, FontSize=legend_font_size, Location=legend_location);

function [data, x_label_text, y_label_text, x_limits, y_limits] = apply_axis_options( ...
        data, x_axis_mode, tilt_scaled_by_pi, x_label_text, y_label_text, ...
        x_limits, y_limits)
    x_axis_mode = lower(string(x_axis_mode));
    tilt_scaled_by_pi = logical(tilt_scaled_by_pi);

    switch x_axis_mode
        case {"excess_area", "excess", "ea"}
            default_x_label = "Membrane excess area";
        case {"lam", "lambda"}
            data(:, 1) = Lam_from_ea(data(:, 1));
            if ~isempty(x_limits)
                x_limits = sort(Lam_from_ea(x_limits));
            end
            default_x_label = "\Lambda";
        otherwise
            error('Unknown x_axis_mode "%s". Use "excess_area" or "Lam".', x_axis_mode);
    end

    if tilt_scaled_by_pi
        default_y_label = "Tilt angle \psi/\pi";
    else
        data(:, 2) = pi * data(:, 2);
        if ~isempty(y_limits)
            y_limits = pi * y_limits;
        end
        default_y_label = "Tilt angle \psi";
    end

    if strlength(string(x_label_text)) == 0
        x_label_text = default_x_label;
    end
    if strlength(string(y_label_text)) == 0
        y_label_text = default_y_label;
    end
end

function [broussinos_data, gamy_by_point] = extract_broussinos_data(Sd, Da, gamy_values, v_values, tf)
    script_dir = fileparts(mfilename("fullpath"));
    if strlength(script_dir) == 0
        script_dir = pwd;
    end
    remesh_dir = find_remesh_dir(script_dir);
    addpath(remesh_dir);
    data_dir = find_data_dir(script_dir, remesh_dir);

    gamy_values = gamy_values(:).';
    v_values = v_values(:).';
    broussinos_data = NaN(numel(gamy_values) * numel(v_values), 2);
    gamy_by_point = NaN(numel(gamy_values) * numel(v_values), 1);

    row = 0;
    for g = 1:numel(gamy_values)
        gamy = gamy_values(g);
        for i = 1:numel(v_values)
            row = row + 1;
            v = v_values(i);
            run_tag = make_run_tag(Sd, Da, gamy, v);
            folder = fullfile(data_dir, run_tag);
            frames = list_geo_frames(folder);
            if isempty(frames)
                warning("No geo*.mat files found for gamy = %.6g, v = %.6g in %s", ...
                    gamy, v, folder);
                continue
            end

            frame_id = select_frame(frames, tf);
            if isempty(frame_id)
                warning("Requested frame %s not found for gamy = %.6g, v = %.6g in %s", ...
                    string(tf), gamy, v, folder);
                continue
            end

            data = load(fullfile(folder, sprintf("geo%d.mat", frame_id)), "M", "P");
            geo = Geometry(data.M, data.P);
            tilt = vesicleTiltDeformation(data.P, data.M, 1, 3, "ray_intersection");
            broussinos_data(row, :) = [excess_area(geo), tilt.psi_over_pi];
            gamy_by_point(row) = gamy;
        end
    end

    valid = all(isfinite(broussinos_data), 2) & isfinite(gamy_by_point);
    broussinos_data = broussinos_data(valid, :);
    gamy_by_point = gamy_by_point(valid);
end

function plot_broussinos_scatter(ax, data, gamy_by_point, gamy_values, color_order, color_map, ...
        single_edge_color, marker_size, line_width)
    if isempty(data)
        warning("No Broussinos/current-work data was extracted.");
        return
    end

    gamy_values = gamy_values(:).';
    if isempty(color_order)
        if isscalar(gamy_values)
            color_order = single_edge_color;
        else
            color_order = feval(char(color_map), numel(gamy_values));
        end
    end

    for g = 1:numel(gamy_values)
        gamy = gamy_values(g);
        mask = abs(gamy_by_point - gamy) <= max(eps(abs(gamy)), eps);
        if ~any(mask)
            continue
        end

        if ischar(color_order) || isstring(color_order)
            marker_color = color_order;
        else
            color_idx = mod(g - 1, size(color_order, 1)) + 1;
            marker_color = color_order(color_idx, :);
        end

        if isscalar(gamy_values)
            display_name = "Current work";
        else
            display_name = sprintf("Current work, \\gamma = %.3g", gamy);
        end

        scatter(ax, data(mask, 1), data(mask, 2), marker_size, ...
            LineWidth=line_width, MarkerEdgeColor=marker_color, ...
            DisplayName=display_name);
    end
end

function frame_id = select_frame(frames, tf)
    frame_id = [];
    if isstring(tf) || ischar(tf)
        tf = string(tf);
        if lower(tf) == "last"
            frame_id = max(frames);
            return
        end
        requested = str2double(tf);
    else
        requested = double(tf);
    end

    if ismember(requested, frames)
        frame_id = requested;
    end
end

function frames = list_geo_frames(folder)
    files = dir(fullfile(folder, "geo*.mat"));
    frames = zeros(numel(files), 1);
    keep = false(numel(files), 1);
    for i = 1:numel(files)
        token = regexp(files(i).name, "^geo(\d+)\.mat$", "tokens", "once");
        if ~isempty(token)
            frames(i) = str2double(token{1});
            keep(i) = true;
        end
    end
    frames = sort(frames(keep));
end

function run_tag = make_run_tag(Sd, Da, gamy, v)
    run_tag = sprintf("Sd_%.2e_Da_%.2e_gamy_%+.2e_v_%.2e", Sd, Da, gamy, v);
    run_tag = strrep(run_tag, "+", "p");
    run_tag = strrep(run_tag, "-", "m");
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
