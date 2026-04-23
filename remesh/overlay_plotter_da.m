close all; clc;

% Animated overlay of x-z cross-sections at y = 0 for two Da runs.

gamy = 1;
Sd = 1;
Da_values = [0, 1];

directory = "./data/ms_batch_data/";
name = "./data/overlay_Da";

edge = 1.5;
plane_tol = 1e-10;
time_stride = 30;
line_width = 2.0;
show_filled = false;
show_frame_text = true;

curve_colors = [
    0.10, 0.35, 0.85;
    0.85, 0.25, 0.15
];

run_tags = strings(size(Da_values));
for k = 1:numel(Da_values)
    run_tags(k) = make_run_tag(gamy, Da_values(k), Sd);
end

folders = directory + run_tags;
last_frames = zeros(size(Da_values));
for k = 1:numel(Da_values)
    files = dir(fullfile(folders(k), "geo*.mat"));
    last_frames(k) = numel(files);
end

tf = min(last_frames);
frame_ids = 1:max(1, time_stride):tf;

delete(name)
v = VideoWriter(name);
v.FrameRate = 14;
open(v)

figure("Color", "w");
size_x = 700;
size_y = 700;
set(gcf, "Position", [100, 100, size_x, size_y]);
set(gcf, "Resize", "off");

for j = 1:numel(frame_ids)
    n = frame_ids(j);
    clf;
    ax = axes(gcf);
    hold(ax, "on");

    legend_handles = gobjects(1, numel(Da_values));

    for k = 1:numel(Da_values)
        geoj = load(folders(k) + sprintf("/geo%d.mat", n));
        [loop_xz, loop_area] = extract_cross_section_xz(geoj.P, geoj.M, plane_tol);

        if isempty(loop_xz)
            continue
        end

        if show_filled
            patch(ax, ...
                "XData", loop_xz(:, 1), ...
                "YData", loop_xz(:, 2), ...
                "FaceColor", curve_colors(k, :), ...
                "FaceAlpha", 0.12, ...
                "EdgeColor", "none");
        end

        legend_handles(k) = plot(ax, loop_xz(:, 1), loop_xz(:, 2), ...
            "Color", curve_colors(k, :), ...
            "LineWidth", line_width, ...
            "DisplayName", sprintf("Da = %.2f, A_{xz} = %.4f", Da_values(k), loop_area));
    end

    axis(ax, "equal");
    xlim(ax, [-edge, edge]);
    ylim(ax, [-edge, edge]);
    xlabel(ax, "x");
    ylabel(ax, "z");
    box(ax, "on");
    grid(ax, "on");
    title(ax, sprintf("y = 0 cross-section overlay, frame %d", n));

    valid_handles = legend_handles(isgraphics(legend_handles));
    if ~isempty(valid_handles)
        legend(ax, valid_handles, "Location", "southoutside");
    end

    if show_frame_text
        text(ax, 0.02, 0.98, sprintf("frame %d / %d", n, tf), ...
            "Units", "normalized", ...
            "HorizontalAlignment", "left", ...
            "VerticalAlignment", "top");
    end

    drawnow;
    frame = getframe(gcf);
    frame_rgb = imresize(frame.cdata, [size_y, size_x]);
    writeVideo(v, frame_rgb);
end

close(v)


function run_tag = make_run_tag(gamy, Da, Sd)
run_tag = sprintf("gamy_%+.2f_Da_%+.2f_Sd_%+.2f", gamy, Da, Sd);
run_tag = strrep(run_tag, "+", "p");
run_tag = strrep(run_tag, "-", "m");
run_tag = string(run_tag);
end


function [loop_xz, loop_area] = extract_cross_section_xz(P, F, plane_tol)
segments = zeros(0, 2, 2);

for face_idx = 1:size(F, 1)
    tri = P(F(face_idx, :), :);
    points = zeros(0, 3);
    edges = [1 2; 2 3; 3 1];

    for edge_idx = 1:3
        a = tri(edges(edge_idx, 1), :);
        b = tri(edges(edge_idx, 2), :);
        ya = a(2);
        yb = b(2);

        if abs(ya) <= plane_tol && abs(yb) <= plane_tol
            points = append_unique_point(points, a, plane_tol);
            points = append_unique_point(points, b, plane_tol);
        elseif abs(ya) <= plane_tol
            points = append_unique_point(points, a, plane_tol);
        elseif abs(yb) <= plane_tol
            points = append_unique_point(points, b, plane_tol);
        elseif (ya < 0 && yb > 0) || (ya > 0 && yb < 0)
            t = -ya / (yb - ya);
            p = a + t * (b - a);
            points = append_unique_point(points, p, plane_tol);
        end
    end

    if size(points, 1) >= 2
        if size(points, 1) > 2
            points = unique(round(points / plane_tol) * plane_tol, "rows", "stable");
        end
        if size(points, 1) >= 2
            segments(end + 1, :, :) = points(1:2, [1, 3]); %#ok<AGROW>
        end
    end
end

if isempty(segments)
    loop_xz = [];
    loop_area = NaN;
    return
end

all_points = reshape(segments, [], 2);
all_points = unique(round(all_points / plane_tol) * plane_tol, "rows", "stable");

if size(all_points, 1) < 3
    loop_xz = [];
    loop_area = NaN;
    return
end

center = mean(all_points, 1);
angles = atan2(all_points(:, 2) - center(2), all_points(:, 1) - center(1));
[~, order] = sort(angles);
loop_xz = all_points(order, :);
loop_xz(end + 1, :) = loop_xz(1, :);
loop_area = polyarea(loop_xz(:, 1), loop_xz(:, 2));
end


function points = append_unique_point(points, point, plane_tol)
if isempty(points)
    points = point;
    return
end

if all(vecnorm(points - point, 2, 2) > 10 * plane_tol)
    points(end + 1, :) = point;
end
end
