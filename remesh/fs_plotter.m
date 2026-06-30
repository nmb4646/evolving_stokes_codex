close all; clc; clear;

% Barebones plotter for fs_multi output.
% Override run_tag/directory/name before running if needed.

dt = 1e-1;
k = 1000;
Sd = 1000;
Da = 1e-2;
Gamma = 0;
gamy = 0;

run_tag = sprintf('Sd_%.2e_Da_%.2e_gamy_%+.2e', Sd, Da, gamy);
run_tag = strrep(run_tag, '+', 'p');
run_tag = strrep(run_tag, '-', 'm');

if ~exist('directory', 'var')
    directory = "./data/fs_batch_data/";
end
folder = directory + run_tag;
name = directory + run_tag + ".gif";
fprintf("Plotting %s\n", folder);
fprintf("Writing %s\n", name);

view_azi = 179;
view_ele = 1;
view_rotation_angle_deg = 0;
view_rotation_axis = [0, 1, 0];
view_rotation_center = [0, 0, 0];
alph = .9;
edge_color = [.3, .3, .3];
gif_delay = 1 / 14;
size_x = 1000;
size_y = 900;
plot_stride = 1;

if ~exist('show_velocity', 'var')
    show_velocity = true;
end
if ~exist('velocity_stride', 'var')
    velocity_stride = 1;
end
if ~exist('velocity_scale', 'var')
    velocity_scale = 4;
end
if ~exist('velocity_color', 'var')
    velocity_color = 'k';
end
if ~exist('color_mode', 'var')
    color_mode = "permeation";
end
if ~exist('uniform_color', 'var')
    uniform_color = [0.75, 0.78, 0.82];
end

files = dir(fullfile(folder, 'geo*.mat'));
if isempty(files)
    error("No geo*.mat files found in %s", folder);
end

frame_ids = zeros(numel(files), 1);
for i = 1:numel(files)
    token = regexp(files(i).name, 'geo(\d+)\.mat', 'tokens', 'once');
    frame_ids(i) = str2double(token{1});
end
frame_ids = sort(frame_ids);
frame_ids = frame_ids(~isnan(frame_ids));

tf = frame_ids(end);
ts = plot_stride; 


plot_ids = frame_ids(1:ts:end).';
if plot_ids(end) ~= tf
    plot_ids(end + 1) = tf;
end
plot_ids = unique(plot_ids, 'stable');

geoj = load(fullfile(folder, sprintf("geo%d.mat", 0)));
geo = Geometry(geoj.M, geoj.P);
volume0 = geo.volume;
geoj = load(fullfile(folder, sprintf("geo%d.mat", tf)));
geo = Geometry(geoj.M, geoj.P);
view_rotation = make_view_rotation(view_rotation_angle_deg, view_rotation_axis, view_rotation_center);
P_view = apply_view_rotation(geoj.P, view_rotation);
use_uniform_color = strcmp(string(color_mode), "uniform");
if use_uniform_color
    C = [];
else
    C = frame_color_data(geoj, geo, color_mode, gamy);
    [color_min, color_max] = percentile_bounds(C, 2, 98);
    if color_max <= color_min
        color_max = color_min + eps;
    end
end
edge = 1.08 * max(abs(P_view), [], "all");
if edge == 0
    edge = 1;
end

if isfile(name)
    delete(name);
end

fig = figure('Color', 'w', 'Position', [200, 200, size_x, size_y], 'InvertHardcopy', 'off');
ax = axes(fig);
set(ax, 'Units', 'normalized', 'Position', [0.03, 0.03, 0.94, 0.94], ...
    'ActivePositionProperty', 'position', 'Color', [.9, .9, .9]);
colormap(ax, turbo(256));
if ~use_uniform_color
    clim(ax, [color_min, color_max]);
end
hold(ax, 'on');

if use_uniform_color
    h = trisurf(geoj.M, P_view(:,1), P_view(:,2), P_view(:,3), ...
        'Parent', ax, ...
        'FaceColor', uniform_color, ...
        'EdgeColor', edge_color, ...
        'FaceAlpha', alph);
else
    h = trisurf(geoj.M, P_view(:,1), P_view(:,2), P_view(:,3), ...
        'Parent', ax, ...
        'FaceColor', 'interp', ...
        'FaceVertexCData', C, ...
        'CDataMapping', 'scaled', ...
        'EdgeColor', edge_color, ...
        'FaceAlpha', alph);
end
if show_velocity
    velocity = reshape(geoj.velocity, [], 3);
    velocity_view = rotate_view_vectors(velocity, view_rotation);
    %velocity = -shear_flow(geoj.P,1);
    %velocity = dot(reshape(geoj.velocity, [], 3),geo.v_normal,2).*geo.v_normal;
    %velocity = reshape(geoj.velocity, [], 3) - dot(reshape(geoj.velocity, [], 3),geo.v_normal,2).*geo.v_normal;
    velocity_ids = 1:velocity_stride:size(geoj.P, 1);
    qv = quiver3(ax, P_view(velocity_ids,1), P_view(velocity_ids,2), P_view(velocity_ids,3), ...
        velocity_view(velocity_ids,1), velocity_view(velocity_ids,2), velocity_view(velocity_ids,3), ...
        velocity_scale, velocity_color);
    %"off",velocity_color);
        
else
    qv = [];
    velocity_ids = [];
end
% lighting(ax, 'gouraud');
% camlight(ax, 'headlight');
axis(ax, 'manual');
axis(ax, 'equal');
axis(ax, 'vis3d');
axis(ax, 'off');
pbaspect(ax, [1 2 1]);
xlim(ax, [-edge, edge]);
ylim(ax, [-edge, edge]);
zlim(ax, [-edge, edge]);
view(ax, view_azi, view_ele);

gif_frame_count = 0;

%%% OVERRIDE

%plot_ids = 0:1:2;


for frame_idx = 1:numel(plot_ids)
    n = plot_ids(frame_idx);
    geoj = load(fullfile(folder, sprintf("geo%d.mat", n)));
    geo = Geometry(geoj.M, geoj.P);
    P_view = apply_view_rotation(geoj.P, view_rotation);

    %disp((geo.volume - 4/3*pi)/(4/3*pi));
    fprintf("Volume: %.9f\n", geo.volume/volume0);

    if ~use_uniform_color
        C = frame_color_data(geoj, geo, color_mode, gamy);
    end

    if use_uniform_color
        set(h, 'Faces', geoj.M, ...
            'Vertices', P_view, ...
            'FaceColor', uniform_color);
    else
        set(h, 'Faces', geoj.M, ...
            'Vertices', P_view, ...
            'FaceVertexCData', C);
    end
    if show_velocity
        velocity = reshape(geoj.velocity, [], 3);
        velocity_view = rotate_view_vectors(velocity, view_rotation);
        %velocity = -shear_flow(geoj.P,1);
        %velocity = dot(reshape(geoj.velocity, [], 3),geo.v_normal,2).*geo.v_normal;
        %velocity = velocity - dot(velocity, geo.v_normal, 2) .* geo.v_normal;
        velocity_ids = 1:velocity_stride:size(geoj.P, 1);
        set(qv, ...
            'XData', P_view(velocity_ids,1), ...
            'YData', P_view(velocity_ids,2), ...
            'ZData', P_view(velocity_ids,3), ...
            'UData', velocity_view(velocity_ids,1), ...
            'VData', velocity_view(velocity_ids,2), ...
            'WData', velocity_view(velocity_ids,3));
    end
    if ~use_uniform_color
        colorbar;%clim([0 2])
    end
    drawnow;

    frame = getframe(fig);
    frame_rgb = frame2im(frame);
    [frame_ind, map] = rgb2ind(frame_rgb, 256);
    if gif_frame_count == 0
        imwrite(frame_ind, map, name, 'gif', ...
            'LoopCount', inf, ...
            'DelayTime', gif_delay, ...
            'DisposalMethod', 'restoreBG');
    else
        imwrite(frame_ind, map, name, 'gif', ...
            'WriteMode', 'append', ...
            'DelayTime', gif_delay, ...
            'DisposalMethod', 'restoreBG');
    end
    gif_frame_count = gif_frame_count + 1;
    fprintf("Wrote frame %d/%d from geo%d.mat\n", frame_idx, numel(plot_ids), n);
end

fprintf("Saved %s with %d frames.\n", name, gif_frame_count);

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

function view_rotation = make_view_rotation(angle_deg, axis, center)
    if nargin < 1 || isempty(angle_deg)
        angle_deg = 0;
    end
    if nargin < 2 || isempty(axis)
        axis = [0, 0, 1];
    end
    if nargin < 3 || isempty(center)
        center = [0, 0, 0];
    end

    axis = parse_rotation_axis(axis);
    theta = deg2rad(angle_deg);
    K = [0, -axis(3), axis(2); ...
         axis(3), 0, -axis(1); ...
         -axis(2), axis(1), 0];
    R = eye(3) + sin(theta) * K + (1 - cos(theta)) * (K * K);

    view_rotation.R = R;
    view_rotation.center = center(:).';
end

function axis = parse_rotation_axis(axis)
    if isstring(axis) || ischar(axis)
        switch lower(string(axis))
            case "x"
                axis = [1, 0, 0];
            case "y"
                axis = [0, 1, 0];
            case "z"
                axis = [0, 0, 1];
            otherwise
                error("Unknown view_rotation_axis '%s'. Use 'x', 'y', 'z', or a 1x3 vector.", axis);
        end
    end
    axis = double(axis(:).');
    if numel(axis) ~= 3 || any(~isfinite(axis)) || norm(axis) == 0
        error("view_rotation_axis must be 'x', 'y', 'z', or a finite nonzero 1x3 vector.");
    end
    axis = axis / norm(axis);
end

function P_view = apply_view_rotation(P, view_rotation)
    P_view = (P - view_rotation.center) * view_rotation.R.' + view_rotation.center;
end

function V_view = rotate_view_vectors(V, view_rotation)
    V_view = V * view_rotation.R.';
end

function C = frame_color_data(geoj, geo, color_mode, ~)
    switch string(color_mode)
        case "uniform"
            C = [];
        case "curvature"
            C = geo.v_mean_curvature ./ geo.v_area;
        case {"permeation", "permeation_velocity", "permeation_norm"}
            if ~isfield(geoj, 'f')
                error("Cannot color by permeation velocity because this frame does not contain f.");
            end
            if isfield(geoj, 'p') && isfield(geoj.p, 'Gamma')
                Gamma = geoj.p.Gamma;
            else
                Gamma = 0;
            end
            C = abs(Gamma + dot(geoj.f, geo.v_normal, 2));
        otherwise
            error("Unknown color_mode '%s'. Use 'curvature' or 'permeation_velocity'.", color_mode);
    end
end
