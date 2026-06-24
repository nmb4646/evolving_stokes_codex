close all; clc;

% Reproduce the fs_multi restart-remesh path for the gamy = 4 run.
% This intentionally mirrors the start ~= 0, p.remesh_size ~= 0 branch.

run_folder = "./data/fs_batch_data/Sd_1.00ep01_Da_1.00ep02_gamy_p4.00ep00/";
frame_id = 79;
remesh_size = 1;
smooth_steps = 2;
smooth_alpha = .25;
smooth_target_valences = [5, 7];
output_dir = "./data/remesher_testing/";
if ~exist(output_dir, "dir")
    mkdir(output_dir);
end

addpath('./isoremesh')

loaded = load(run_folder + sprintf("geo%d.mat", frame_id), ...
    "M", "P", "velocity", "lambda", "f", "p", "o", "r");

M_pre = loaded.M;
P_pre = loaded.P;
velocity_pre = loaded.velocity;
f_pre = loaded.f;
p = loaded.p;
o = loaded.o;
r = loaded.r;
lambda = loaded.lambda;

geo_pre = Geometry(M_pre, P_pre);
fprintf("Loaded geo%d from %s\n", frame_id, run_folder);
fprintf("Old vertices: %d, faces: %d, area: %.16g\n", ...
    size(P_pre, 1), size(M_pre, 1), geo_pre.area);

r.edge_length = remesh_size * mean(geo_pre.he_length);
fprintf("Restart remesh target edge length: %.16g\n", r.edge_length);

[M, P] = remeshing(int32(M_pre), P_pre, int32([]), r.edge_length, int32(100));


M = cast(M, "double");
geo = Geometry(M, P);
fprintf("After remeshing, before area rescale: vertices: %d, faces: %d, area: %.16g\n", ...
    size(P, 1), size(M, 1), geo.area);

if smooth_steps > 0 && smooth_alpha ~= 0
    [P, smooth_info] = valence_neighbor_average_smooth(M, P, ...
        smooth_target_valences, smooth_steps, smooth_alpha);
    geo = Geometry(M, P);
    fprintf("After valence neighbor averaging (%d steps, alpha %.4g), before area rescale: area: %.16g\n", ...
        smooth_steps, smooth_alpha, geo.area);
    fprintf("Smoothed %d vertices with valence in [%s]; max displacement: %.16g, mean displacement: %.16g\n", ...
        smooth_info.n_target, sprintf("%d ", smooth_target_valences), ...
        max(smooth_info.max_displacement), max(smooth_info.mean_displacement));
end

P = P * sqrt(p.area0 / geo.area);
geo = Geometry(M, P);
fprintf("After area rescale: vertices: %d, faces: %d, area: %.16g\n", ...
    size(P, 1), size(M, 1), geo.area);



[velocity, f, transfer] = map_data_debug(geo, geo_pre, velocity_pre, f_pre);

velocity_pre_mat = reshape(velocity_pre, [], 3);
velocity_mat = reshape(velocity, [], 3);
fprintf("Projection max distance: %.16g\n", max(transfer.distance));
fprintf("Projection mean distance: %.16g\n", mean(transfer.distance));
fprintf("Projection max bary overshoot: %.16g\n", max(transfer.bary_overshoot));
fprintf("Old velocity RMS: %.16g\n", rms_vec(velocity_pre_mat));
fprintf("Mapped velocity RMS: %.16g\n", rms_vec(velocity_mat));
fprintf("Old f RMS: %.16g\n", rms_vec(f_pre));
fprintf("Mapped f RMS: %.16g\n", rms_vec(f));

save(run_folder + "remesher_testing_geo79_remeshed.mat", ...
    "M", "P", "velocity", "lambda", "f", "p", "o", "r", ...
    "M_pre", "P_pre", "velocity_pre", "f_pre", "transfer");
fprintf("Saved %s\n", run_folder + "remesher_testing_geo79_remeshed.mat");

delete_if_exists(output_dir + "geo79_remesh_curvature.png");
delete_if_exists(output_dir + "geo79_remesh_velocity.png");
delete_if_exists(output_dir + "geo79_remesh_projection_distance.png");

view_azi = 179;
view_ele = 1;
alph = .9;
edge_color = [.3, .3, .3];
velocity_stride = 1;
velocity_scale = 4;
velocity_color = 'k';

C_pre = geo_pre.v_mean_curvature ./ geo_pre.v_area;
C = geo.v_mean_curvature ./ geo.v_area;
[color_min, color_max] = percentile_bounds([C_pre; C], 2, 98);
if color_max <= color_min
    color_max = color_min + eps;
end

plot_limits = padded_limits([P_pre; P], 0.04);

fig = figure("Color", "w", "Visible", "off", "Position", [100, 100, 1800, 850], ...
    "InvertHardcopy", "off");
ax = axes(fig, "Position", [0.01, 0.035, 0.49, 0.93], "Color", [.9, .9, .9]);
plot_curvature_velocity_panel(ax, M_pre, P_pre, C_pre, velocity_pre_mat, ...
    "geo79 before remesh", color_min, color_max, plot_limits, view_azi, view_ele, ...
    alph, edge_color, velocity_stride, velocity_scale, velocity_color);

ax = axes(fig, "Position", [0.50, 0.035, 0.49, 0.93], "Color", [.9, .9, .9]);
plot_curvature_velocity_panel(ax, M, P, C, velocity_mat, ...
    "geo79 after restart remesh + map", color_min, color_max, plot_limits, view_azi, view_ele, ...
    alph, edge_color, velocity_stride, velocity_scale, velocity_color);

output_file = output_dir + "geo79_remesh_curvature_velocity.png";
exportgraphics(fig, output_file, "Resolution", 200);
fprintf("Saved %s\n", output_file);

function [velocity, f, transfer] = map_data_debug(geo, geo_pre, velocity_pre, f_pre)
    kdtree = KDTreeSearcher(geo_pre.f_center);
    [face, uv, count, fail] = project(geo_pre.V, geo_pre.F, geo.V, kdtree, 6);
    if fail
        error("projection failed.");
    end

    velocity = interpolate(geo_pre.F, face, uv, reshape(velocity_pre, [], 3));
    velocity = velocity(:);
    f = interpolate(geo_pre.F, face, uv, f_pre);

    projected = interpolate(geo_pre.F, face, uv, geo_pre.V);
    bary3 = 1 - uv(:, 1) - uv(:, 2);
    bary = [uv, bary3];

    transfer.face = face;
    transfer.uv = uv;
    transfer.count = count;
    transfer.projected = projected;
    transfer.distance = vecnorm(geo.V - projected, 2, 2);
    transfer.bary = bary;
    transfer.bary_overshoot = max(max(-bary, 0), [], 2) + max(max(bary - 1, 0), [], 2);
end

function value = rms_vec(x)
    value = norm(x(:)) / sqrt(numel(x));
end

function plot_curvature_velocity_panel(ax, M, P, C, velocity, plot_title, color_min, color_max, plot_limits, ...
        view_azi, view_ele, alph, edge_color, velocity_stride, velocity_scale, velocity_color)
    trisurf(M, P(:, 1), P(:, 2), P(:, 3), ...
        "Parent", ax, ...
        "FaceColor", "interp", ...
        "FaceVertexCData", C, ...
        "CDataMapping", "scaled", ...
        "EdgeColor", edge_color, ...
        "FaceAlpha", alph);
    hold(ax, "on");
    colormap(ax, turbo(256));
    clim(ax, [color_min, color_max]);

    velocity_ids = 1:velocity_stride:size(P, 1);
    quiver3(ax, P(velocity_ids, 1), P(velocity_ids, 2), P(velocity_ids, 3), ...
        velocity(velocity_ids, 1), velocity(velocity_ids, 2), velocity(velocity_ids, 3), ...
        velocity_scale, velocity_color);

    axis(ax, "manual");
    axis(ax, "equal");
    axis(ax, "vis3d");
    axis(ax, "off");
    daspect(ax, [1, 1, 1]);
    xlim(ax, plot_limits.x);
    ylim(ax, plot_limits.y);
    zlim(ax, plot_limits.z);
    view(ax, view_azi, view_ele);
    title(ax, plot_title);
end

function [lo, hi] = percentile_bounds(x, lo_pct, hi_pct)
    x = sort(x(:));
    x = x(isfinite(x));
    if isempty(x)
        lo = 0;
        hi = 1;
        return
    end
    n = numel(x);
    lo_idx = max(1, min(n, round(1 + (n - 1) * lo_pct / 100)));
    hi_idx = max(1, min(n, round(1 + (n - 1) * hi_pct / 100)));
    lo = x(lo_idx);
    hi = x(hi_idx);
end

function delete_if_exists(path)
    if isfile(path)
        delete(path);
    end
end

function limits = padded_limits(P, pad_fraction)
    lo = min(P, [], 1);
    hi = max(P, [], 1);
    span = hi - lo;
    max_span = max(span);
    if max_span == 0
        max_span = 1;
    end
    span(span == 0) = max_span;
    pad = pad_fraction * span;
    limits.x = [lo(1) - pad(1), hi(1) + pad(1)];
    limits.y = [lo(2) - pad(2), hi(2) + pad(2)];
    limits.z = [lo(3) - pad(3), hi(3) + pad(3)];
end
