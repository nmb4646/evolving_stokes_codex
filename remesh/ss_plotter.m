close all; clc;

% Barebones plotter for ss_multi output.
% Override run_tag/directory/name before running if needed.

dt = 1e-2;
k = 10000;
run_tag = sprintf('dt_%.2e_k_%.2e', dt, k);
run_tag = strrep(run_tag, '+', 'p');
run_tag = strrep(run_tag, '-', 'm');

if ~exist('directory', 'var')
    directory = "./data/ss_batch_data/";
end
folder = directory + run_tag;
name = directory + run_tag + ".gif";
fprintf("Plotting %s\n", folder);
fprintf("Writing %s\n", name);

view_azi = -60;
view_ele = -40;
alph = .9;
edge_color = [.3, .3, .3];
gif_delay = 1 / 14;
size_x = 1000;
size_y = 900;
plot_stride = 1;

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
ts=3;

plot_ids = frame_ids(1:ts:end).';
if plot_ids(end) ~= tf
    plot_ids(end + 1) = tf;
end
plot_ids = unique(plot_ids, 'stable');

geoj = load(fullfile(folder, sprintf("geo%d.mat", tf)));
geo = Geometry(geoj.M, geoj.P);
H = geo.v_mean_curvature ./ geo.v_area;
[curv_min, curv_max] = percentile_bounds(H, 2, 98);
if curv_max <= curv_min
    curv_max = curv_min + eps;
end
edge = 1.08 * max(abs(geoj.P), [], "all");
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
clim(ax, [curv_min, curv_max]);
hold(ax, 'on');

h = trisurf(geoj.M, geoj.P(:,1), geoj.P(:,2), geoj.P(:,3), ...
    'Parent', ax, ...
    'FaceColor', 'interp', ...
    'FaceVertexCData', H, ...
    'CDataMapping', 'scaled', ...
    'EdgeColor', edge_color, ...
    'FaceAlpha', alph);
% lighting(ax, 'gouraud');
% camlight(ax, 'headlight');
axis(ax, 'manual');
axis(ax, 'equal');
axis(ax, 'vis3d');
axis(ax, 'off');
pbaspect(ax, [1 1 1]);
xlim(ax, [-edge, edge]);
ylim(ax, [-edge, edge]);
zlim(ax, [-edge, edge]);
view(ax, view_azi, view_ele);

gif_frame_count = 0;



for frame_idx = 1:numel(plot_ids)
    n = plot_ids(frame_idx);
    geoj = load(fullfile(folder, sprintf("geo%d.mat", n)));
    geo = Geometry(geoj.M, geoj.P);
    H = geo.v_mean_curvature ./ geo.v_area;

    set(h, 'Faces', geoj.M, ...
        'Vertices', geoj.P, ...
        'FaceVertexCData', H);
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
