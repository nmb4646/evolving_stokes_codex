close all; clc;

%%% Simulation selection
% Edit these three values to select a run from data/fs_batch_data.
Da = 0;
Sd = 1e-3;
gamy = 1e-3;

%%% Tilt calculation
% For shear_flow.m, u_x = gamy * z, so flow_dim = 1 and grad_dim = 3.
flow_dim = 1;
grad_dim = 3;

% If true, remove +/- pi branch jumps using the axis-pex   riodic angle.
% Leave false to match the wrapped psi/pi values printed by fs_plotter.
tilt_options.unwrap = false;

%%% Visual settings
% Options: "over_pi" matches fs_plotter, "radians", or "degrees".
visual.angle_units = "over_pi";
visual.line_width = 2.0;
visual.marker = "none";
visual.marker_size = 5;
visual.figure_position = [140, 140, 900, 560];
visual.x_scale = "linear"; % Options: "linear" or "log".
visual.y_scale = "linear"; % Options: "linear" or "log".
visual.show_grid = true;
visual.save_figure = false;
visual.output_file = fullfile("data", "tilt_over_time.png");
visual.output_resolution = 300;

script_dir = fileparts(mfilename("fullpath"));
remesh_dir = find_remesh_dir(script_dir);
addpath(remesh_dir);

data_dir = find_data_dir(script_dir, remesh_dir);
run_tag = make_run_tag(Sd, Da, gamy);
folder = fullfile(data_dir, run_tag);

fprintf("Loading %s\n", folder);

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

if isempty(frame_ids)
    error("No valid geo*.mat frame numbers found in %s", folder);
end

tilt = zeros(numel(frame_ids), 1);
dt = NaN;

for i = 1:numel(frame_ids)
    frame_path = fullfile(folder, files(i).name);
    data = load(frame_path, "M", "P", "p");

    if i == 1 && isfield(data, "p") && isfield(data.p, "dt")
        dt = data.p.dt;
    end

    tilt_out = vesicleTiltDeformation(data.P, data.M, flow_dim, grad_dim);
    tilt(i) = tilt_out.psi;
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

switch lower(string(visual.angle_units))
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
        error("Unknown visual.angle_units '%s'. Use 'over_pi', 'radians', or 'degrees'.", visual.angle_units);
end

fig = figure("Color", "w", "Position", visual.figure_position);
ax = axes(fig);
plot(ax, x, y, ...
    "LineWidth", visual.line_width, ...
    "Marker", visual.marker, ...
    "MarkerSize", visual.marker_size);

xlabel(ax, x_label);
ylabel(ax, y_label);
title(ax, sprintf("Tilt over time: Sd = %.4g, Da = %.4g, gamy = %.4g", Sd, Da, gamy));
set(ax, "XScale", visual.x_scale);
set(ax, "YScale", visual.y_scale);

if visual.show_grid
    grid(ax, "on");
end
box(ax, "on");

fprintf("Initial tilt: %.8g psi/pi, %.8g rad, %.8g deg\n", ...
    tilt(1) / pi, tilt(1), rad2deg(tilt(1)));
fprintf("Final tilt:   %.8g psi/pi, %.8g rad, %.8g deg\n", ...
    tilt(end) / pi, tilt(end), rad2deg(tilt(end)));

if visual.save_figure
    exportgraphics(fig, visual.output_file, "Resolution", visual.output_resolution);
    fprintf("Saved %s\n", visual.output_file);
end

function run_tag = make_run_tag(Sd, Da, gamy)
    run_tag = sprintf("Sd_%.2e_Da_%.2e_gamy_%+.2e", Sd, Da, gamy);
    run_tag = strrep(run_tag, "+", "p");
    run_tag = strrep(run_tag, "-", "m");
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
