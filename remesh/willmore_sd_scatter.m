%WILLMORE_SD_SCATTER Scatter inverse normalized Willmore drop over Sd.
%
% For each Sd in Sd_values, this script loads geo0.mat and geo<tf>.mat from
% the corresponding fs_batch output folder and plots
%
%   y = ((W0 - Wtf) / (Sd * dt))^(-1)
%
% on log-log axes against Sd, where W = geo.willmore_energy(1).

clearvars;
close all;
clc;

%%% Parameters
Da = 0;
v = 0.98;
gamy = 0;
Sd_values = 10 .^ (-2:1:2);

tf = "last";               % "last" or a numeric frame id, e.g. 50
initial_frame = 0;
time_denominator_mode = "dt"; % "dt" uses Sd*dt; "total_time" uses Sd*(tf-initial_frame)*dt.

%%% Visual settings
visual.figure_position = [120, 120, 760, 620];
visual.figure_color = "w";
visual.marker = "o";
visual.marker_size = 80;
visual.marker_face_color = [0.10, 0.35, 0.80];
visual.marker_edge_color = "k";
visual.line_width = 1.0;
visual.axes_font_size = 13;
visual.label_font_size = 15;
visual.title_font_size = 16;
visual.show_grid = true;
visual.show_box = true;
visual.x_limits = [];
visual.y_limits = [];
visual.save_figure = false;
visual.output_file = fullfile("data", "willmore_sd_scatter.png");
visual.output_resolution = 300;

%%% Locate data
script_dir = fileparts(mfilename("fullpath"));
if strlength(script_dir) == 0
    script_dir = pwd;
end
remesh_dir = find_remesh_dir(script_dir);
addpath(remesh_dir);
data_dir = find_data_dir(script_dir, remesh_dir);

summary = table('Size', [0, 8], ...
    'VariableTypes', ["double", "double", "double", "double", "double", "double", "double", "double"], ...
    'VariableNames', ["Sd", "frame_final", "dt", "W0", "Wtf", "energy_drop", ...
        "normalized_drop_rate", "inverse_normalized_drop_rate"]);

for i = 1:numel(Sd_values)
    Sd = Sd_values(i);
    run_tag = make_run_tag(Sd, Da, gamy, v);
    folder = fullfile(data_dir, run_tag);
    frames = list_geo_frames(folder);
    if isempty(frames)
        warning("No geo*.mat files found for Sd = %.4g in %s", Sd, folder);
        continue
    end

    if ~ismember(initial_frame, frames)
        warning("Skipping Sd = %.4g because geo%d.mat was not found in %s.", ...
            Sd, initial_frame, folder);
        continue
    end

    frame_final = select_final_frame(frames, tf);
    if isempty(frame_final)
        warning("Skipping Sd = %.4g because requested tf = %s was not found in %s.", ...
            Sd, string(tf), folder);
        continue
    end

    data0 = load(fullfile(folder, sprintf("geo%d.mat", initial_frame)), "M", "P", "p");
    dataf = load(fullfile(folder, sprintf("geo%d.mat", frame_final)), "M", "P", "p");

    geo0 = Geometry(data0.M, data0.P);
    geof = Geometry(dataf.M, dataf.P);
    W0 = geo0.willmore_energy(1);
    Wtf = geof.willmore_energy(1);
    dt = extract_dt(dataf, data0);

    time_factor = dt;
    if time_denominator_mode == "total_time"
        time_factor = max(frame_final - initial_frame, 1) * dt;
    end

    energy_drop = W0 - Wtf;
    normalized_drop_rate = energy_drop / (Sd * time_factor);
    inverse_normalized_drop_rate = 1 / normalized_drop_rate;

    summary(end + 1, :) = {Sd, frame_final, dt, W0, Wtf, energy_drop, ...
        normalized_drop_rate, inverse_normalized_drop_rate}; %#ok<SAGROW>
end

disp(summary);

valid = isfinite(summary.Sd) ...
    & isfinite(summary.inverse_normalized_drop_rate) ...
    & summary.Sd > 0 ...
    & summary.inverse_normalized_drop_rate > 0;

if nnz(valid) < height(summary)
    warning("Omitting %d nonpositive or nonfinite points from the log-log scatter.", ...
        height(summary) - nnz(valid));
end

figure("Color", visual.figure_color, "Position", visual.figure_position);
ax = axes;
hold on;
scatter(summary.Sd(valid), summary.inverse_normalized_drop_rate(valid), ...
    visual.marker_size, ...
    "Marker", visual.marker, ...
    "MarkerFaceColor", visual.marker_face_color, ...
    "MarkerEdgeColor", visual.marker_edge_color, ...
    "LineWidth", visual.line_width);

set(ax, ...
    "XScale", "log", ...
    "YScale", "log", ...
    "FontSize", visual.axes_font_size);
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
if ~isempty(visual.x_limits)
    xlim(visual.x_limits);
end
if ~isempty(visual.y_limits)
    ylim(visual.y_limits);
end

xlabel("Sd", "FontSize", visual.label_font_size);
if time_denominator_mode == "total_time"
    ylabel("((W_0-W_{tf})/(Sd \Delta t_{total}))^{-1}", ...
        "FontSize", visual.label_font_size);
else
    ylabel("((W_0-W_{tf})/(Sd dt))^{-1}", ...
        "FontSize", visual.label_font_size);
end
title(sprintf("Inverse normalized Willmore drop, Da = %.3g, v = %.3g, gamy = %.3g", ...
    Da, v, gamy), "FontSize", visual.title_font_size);

if visual.save_figure
    exportgraphics(gcf, visual.output_file, "Resolution", visual.output_resolution);
end

function frame_final = select_final_frame(frames, tf)
    frame_final = [];
    if isstring(tf) || ischar(tf)
        tf = string(tf);
        if lower(tf) == "last"
            frame_final = max(frames);
        else
            requested = str2double(tf);
            if ismember(requested, frames)
                frame_final = requested;
            end
        end
    else
        requested = double(tf);
        if ismember(requested, frames)
            frame_final = requested;
        end
    end
end

function dt = extract_dt(dataf, data0)
    if isfield(dataf, "p") && isfield(dataf.p, "dt")
        dt = dataf.p.dt;
    elseif isfield(data0, "p") && isfield(data0.p, "dt")
        dt = data0.p.dt;
    else
        warning("No p.dt found; using dt = 1.");
        dt = 1;
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
