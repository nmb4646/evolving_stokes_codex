%CURVATURE_REMESH_TESTER Visual test for curvature-weighted remeshing.
%
% This script loads one saved geo*.mat file, applies remeshing_curvature,
% and compares the triangulated surface before/after.

clearvars;
close all;

%%% Parameters
run_folder_name = "Sd_1.00em06_Da_0.00ep00_gamy_p3.70em05_v_3.50em01";
geo_timestep = "last";        % "last" or a numeric timestep, e.g. 25
curvature_source = "max_abs_principal";  % "max_abs_principal", "v_mean", "he_mean", "v_gaussian", "f_gaussian"

target_edge_length_factor = 1.0;  % Multiplies mean pre-remesh half-edge length
remesh_iterations = 10;

remesh_opts = struct();
remesh_opts.curvature_weight = 8;
remesh_opts.curvature_power = 2;
remesh_opts.min_edge_length_factor = 0.1;
remesh_opts.max_edge_length_factor = 10;
remesh_opts.curvature_smooth_iterations = 3;
remesh_opts.target_smooth_iterations = 2;

remesh_opts.split_factor = 4 / 3;
remesh_opts.max_vertex_count_factor = 1.1;   % or 1.10
remesh_opts.count_tolerance = 20;             % optional explicit tolerance
remesh_opts.max_collapses_per_iteration = 600;
remesh_opts.forced_collapse_factor = 1.15;
remesh_opts.collapse_factor = 0.90;

remesh_opts.flip_passes_per_iteration = 8;
remesh_opts.relaxation_iterations_per_iteration = 2;
remesh_opts.relaxation_alpha = 0.85;
remesh_opts.polish_iterations = 3;
remesh_opts.project_to_input = true;
remesh_opts.projection_neighbors = 20;
remesh_opts.verbose = true;

visual = struct();
visual.figure_position = [120, 60, 760, 1100];
visual.color_by_curvature = true;
visual.shared_color_limits = true;
visual.edge_color = [0.08, 0.08, 0.08];
visual.edge_alpha = 0.55;
visual.edge_width = 0.25;
visual.face_color = [0.78, 0.84, 0.91];
visual.face_alpha = 1;
visual.view_azimuth = 35;
visual.view_elevation = 20;
visual.font_size = 13;
visual.title_font_size = 14;
visual.colormap = turbo(256);

%%% Load geometry
script_dir = fileparts(mfilename("fullpath"));
if strlength(script_dir) == 0
    script_dir = pwd;
end
remesh_dir = find_remesh_dir(script_dir);
addpath(remesh_dir);

data_dir = find_data_dir(script_dir, remesh_dir);
run_dir = fullfile(data_dir, run_folder_name);
if ~isfolder(run_dir)
    error("Run folder not found: %s", run_dir);
end

geo_file = select_geo_file(run_dir, geo_timestep);
[~, geo_file_base] = fileparts(geo_file);
loaded = load(geo_file, "M", "P");
M_before = double(loaded.M);
P_before = double(loaded.P);

geo_before = Geometry(M_before, P_before);
curvature = select_curvature(geo_before, curvature_source);
target_edge_length = target_edge_length_factor * mean(geo_before.he_length);

fprintf("Loaded %s\n", geo_file);
print_geo_stats("before", M_before, P_before, geo_before);

%%% Remesh
[M_after, P_after, info] = remeshing_curvature(int32(M_before), P_before, ...
    curvature, target_edge_length, int32(remesh_iterations), remesh_opts);
M_after = double(M_after);
geo_after = Geometry(M_after, P_after);

print_geo_stats("after", M_after, P_after, geo_after);
fprintf("remeshing_curvature: splits = %d, collapses = %d, projection_failed = %d\n", ...
    info.total_splits, info.total_collapses, info.projection_failed);
fprintf("bad undirected edge count after remesh = %d\n", count_bad_edges(M_after));

%%% Plot
c_before = abs(select_curvature(geo_before, curvature_source));
c_after = abs(select_curvature(geo_after, curvature_source));
if visual.shared_color_limits
    c_limits = [min([c_before; c_after]), max([c_before; c_after])];
else
    c_limits = [];
end

figure("Color", "w", "Position", visual.figure_position);
subplot(2, 1, 1);
plot_trisurf_mesh(M_before, P_before, c_before, visual, c_limits);
title(sprintf("Before: V=%d, F=%d", size(P_before, 1), size(M_before, 1)), ...
    "FontSize", visual.title_font_size);

subplot(2, 1, 2);
plot_trisurf_mesh(M_after, P_after, c_after, visual, c_limits);
title(sprintf("After: V=%d, F=%d", size(P_after, 1), size(M_after, 1)), ...
    "FontSize", visual.title_font_size);

sgtitle(sprintf("curvature remesh test: %s, %s", run_folder_name, geo_file_base), ...
    "Interpreter", "none", "FontSize", visual.title_font_size);

function geo_file = select_geo_file(run_dir, geo_timestep)
    files = dir(fullfile(run_dir, "geo*.mat"));
    if isempty(files)
        error("No geo*.mat files found in %s", run_dir);
    end

    names = string({files.name});
    tokens = regexp(names, "^geo(\d+)\.mat$", "tokens", "once");
    valid = ~cellfun(@isempty, tokens);
    names = names(valid);
    tokens = tokens(valid);
    indices = cellfun(@(token) str2double(token{1}), tokens);

    if isempty(indices)
        error("No numbered geo*.mat files found in %s", run_dir);
    end

    if isstring(geo_timestep) || ischar(geo_timestep)
        geo_timestep = string(geo_timestep);
        if lower(geo_timestep) == "last"
            [~, idx] = max(indices);
        else
            requested = str2double(geo_timestep);
            idx = find(indices == requested, 1);
        end
    else
        idx = find(indices == double(geo_timestep), 1);
    end

    if isempty(idx) || isnan(idx)
        error("Requested geo timestep %s not found in %s", string(geo_timestep), run_dir);
    end
    geo_file = fullfile(run_dir, names(idx));
end

function curvature = select_curvature(geo, source)
    source = lower(string(source));
    if source == "v_mean"
        curvature = geo.v_mean_curvature;
    elseif source == "he_mean"
        curvature = geo.he_mean_curvature;
    elseif source == "v_gaussian"
        curvature = geo.v_gaussian_curvature;
    elseif source == "f_gaussian"
        curvature = geo.f_gaussian_curvature;
    elseif source == "max_abs_principal"
        curvature = max_abs_principal_curvature(geo);
    else
        error("Unknown curvature_source '%s'.", source);
    end
end

function curvature = max_abs_principal_curvature(geo)
    H = geo.v_mean_curvature ./ max(geo.v_area, eps);
    K = geo.v_gaussian_curvature ./ max(geo.v_area, eps);
    discriminant = max(H .^ 2 - K, 0);
    root = sqrt(discriminant);
    k1 = H + root;
    k2 = H - root;
    curvature = max(abs(k1), abs(k2));
end

function plot_trisurf_mesh(M, P, color_data, visual, c_limits)
    if visual.color_by_curvature
        h = trisurf(M, P(:, 1), P(:, 2), P(:, 3), color_data, ...
            "FaceColor", "interp", ...
            "EdgeColor", visual.edge_color, ...
            "EdgeAlpha", visual.edge_alpha, ...
            "LineWidth", visual.edge_width, ...
            "FaceAlpha", visual.face_alpha);
        colormap(gca, visual.colormap);
        if ~isempty(c_limits) && c_limits(2) > c_limits(1)
            clim(c_limits);
        end
        colorbar;
    else
        h = trisurf(M, P(:, 1), P(:, 2), P(:, 3), ...
            "FaceColor", visual.face_color, ...
            "EdgeColor", visual.edge_color, ...
            "EdgeAlpha", visual.edge_alpha, ...
            "LineWidth", visual.edge_width, ...
            "FaceAlpha", visual.face_alpha);
    end

    h.SpecularStrength = 0.15;
    h.DiffuseStrength = 0.8;
    axis equal;
    axis off;
    view(visual.view_azimuth, visual.view_elevation);
    camlight headlight;
    lighting gouraud;
    set(gca, "FontSize", visual.font_size);
end

function print_geo_stats(label, M, P, geo)
    fprintf("%s: V = %d, F = %d, area = %.12g, volume = %.12g, mean he = %.6g\n", ...
        label, size(P, 1), size(M, 1), geo.area, geo.volume, mean(geo.he_length));
end

function n_bad = count_bad_edges(M)
    edges = [M(:, [1, 2]); M(:, [2, 3]); M(:, [3, 1])];
    [~, ~, edge_id] = unique(sort(edges, 2), "rows");
    edge_count = accumarray(edge_id, 1);
    n_bad = nnz(edge_count ~= 2);
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
