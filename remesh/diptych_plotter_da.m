close all; clc;

% Diptych movie for two Da values with fixed gamy and Sd.

gamy = 1;
Sd = 1;
Da_values = [0, 1];

view_azi = 154;
view_ele = 21;
alph = 1;
edge = 1.5;
color_mode = "uniform"; % options: "curvature", "permeation_flux", "uniform"
show_velocity = false;
show_background_shear_xz = true;
reverse_background_shear_xz = true;
background_shear_samples = 13;
background_shear_scale = 0.9;
background_shear_color = [0.1, 0.1, 0.1];

directory = "./data/ms_batch_data/";
run_tags = strings(size(Da_values));
for k = 1:numel(Da_values)
    run_tags(k) = make_run_tag(gamy, Da_values(k), Sd);
end

name = "./data/diptych_Da";
delete(name)
v = VideoWriter(name);
v.FrameRate = 14;
open(v)

figure('Color', 'w');
size_x = 1040;
size_y = 520;
set(gcf, 'Position', [100, 100, size_x, size_y]);
set(gcf, 'Resize', 'off');
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

folders = directory + run_tags;
last_frames = zeros(size(Da_values));
for k = 1:numel(Da_values)
    files = dir(fullfile(folders(k), 'geo*.mat'));
    last_frames(k) = length(files);
end

tf = min(last_frames);
ti = 1;
ts = max(1, ceil(tf / 1000));
ts = 10;

curvature_min = inf(1, numel(Da_values));
curvature_max = -inf(1, numel(Da_values));
flux_abs_max = ones(1, numel(Da_values));

for k = 1:numel(Da_values)
    geoj = load(folders(k) + sprintf("/geo%d.mat", tf));
    geoj_Geom = Geometry(geoj.M, geoj.P);
    mcs = geoj_Geom.v_mean_curvature;
    A0 = mean(mcs(geoj_Geom.F), 2);
    curvature_min(k) = min(A0);
    curvature_max(k) = max(A0);

    if color_mode == "permeation_flux"
        [~, ~, ~, ~, KTK, DTD] = geoj_Geom.evolving_operators();
        twoHn = geoj_Geom.lap * geoj.P;
        fb = geoj_Geom.bending_force(1);
        fs = reshape(geoj.lambda * twoHn, size(geoj.P));
        fv = reshape((- 2 * (KTK + geoj.p.k * DTD) * geoj.velocity / geoj.p.dt), size(geoj.P));
        qperm = geoj.p.pnd + dot(fb + fs + fv, geoj_Geom.v_normal, 2);
        face_flux = geoj_Geom.f_area .* mean(qperm(geoj_Geom.F), 2);
        flux_abs_max(k) = max(1, max(abs(face_flux)));
    end
end

for n = ti:ts:tf
    for k = 1:numel(Da_values)
        nexttile(k);
        cla;

        geoj = load(folders(k) + sprintf("/geo%d.mat", n));
        geoj_Geom = Geometry(geoj.M, geoj.P);
        disp(sprintf("Da = %.2f, volume = ", Da_values(k)) + geoj_Geom.volume)
        fi = trisurf(geoj.M, geoj.P(:,1), geoj.P(:,2), geoj.P(:,3), ...
            'FaceColor', [.6,.6,.6], 'EdgeColor', [.3,.3,.3], 'FaceAlpha', alph);
        hold on;

        if color_mode == "curvature"
            mcs = geoj_Geom.v_mean_curvature;
            A0 = mean(mcs(geoj_Geom.F), 2);
            denom = max(curvature_max(k) - curvature_min(k), eps);
            set(fi, 'FaceColor', 'flat', ...
                'FaceVertexCData', (A0 - curvature_min(k)) / denom * 255, ...
                'CDataMapping', 'direct');
        elseif color_mode == "permeation_flux"
            [~, ~, ~, ~, KTK, DTD] = geoj_Geom.evolving_operators();
            twoHn = geoj_Geom.lap * geoj.P;
            fb = geoj_Geom.bending_force(1);
            fs = reshape(geoj.lambda * twoHn, size(geoj.P));
            fv = reshape((- 2 * (KTK + geoj.p.k * DTD) * geoj.velocity / geoj.p.dt), size(geoj.P));
            qperm = geoj.p.pnd + dot(fb + fs + fv, geoj_Geom.v_normal, 2);
            face_flux = geoj_Geom.f_area .* mean(qperm(geoj_Geom.F), 2);
            set(fi, 'FaceColor', 'flat', ...
                'FaceVertexCData', face_flux, ...
                'CDataMapping', 'scaled');
            clim([-flux_abs_max(k), flux_abs_max(k)]);
            colormap(turbo(256));
        end

        if show_velocity
            L = size(geoj.P, 1);
            quiver3(geoj.P(:,1), geoj.P(:,2), geoj.P(:,3), ...
                geoj.velocity(1:L), geoj.velocity(L+1:2*L), geoj.velocity(2*L+1:end), ...
                4, 'k');
        end

        if show_background_shear_xz
            x_vals = linspace(-edge, edge, background_shear_samples);
            z_vals = linspace(-edge, edge, background_shear_samples);
            [X_bg, Z_bg] = meshgrid(x_vals, z_vals);
            Y_bg = zeros(size(X_bg));
            bg_points = [X_bg(:), Y_bg(:), Z_bg(:)];
            U_bg = shear_flow(bg_points, gamy);
            if reverse_background_shear_xz
                U_bg = -U_bg;
            end
            quiver3(X_bg, Y_bg, Z_bg, ...
                reshape(U_bg(:,1), size(X_bg)), ...
                reshape(U_bg(:,2), size(X_bg)), ...
                reshape(U_bg(:,3), size(X_bg)), ...
                background_shear_scale, ...
                'Color', background_shear_color, ...
                'LineWidth', 0.75, ...
                'MaxHeadSize', 0.5);
        end

        axis equal;
        xlim([-edge, edge]); ylim([-edge, edge]); zlim([-edge, edge]);
        view(view_azi, view_ele);
        axis off;
        grid off;
        title(sprintf('Da = %.2f', Da_values(k)));
    end

    drawnow;
    frame = getframe(gcf);
    frame_rgb = imresize(frame.cdata, [size_y, size_x]);
    writeVideo(v, frame_rgb);
end

close(v)

function run_tag = make_run_tag(gamy, Da, Sd)
run_tag = sprintf('gamy_%+.2f_Da_%+.2f_Sd_%+.2f', gamy, Da, Sd);
run_tag = strrep(run_tag, '+', 'p');
run_tag = strrep(run_tag, '-', 'm');
run_tag = string(run_tag);
end
