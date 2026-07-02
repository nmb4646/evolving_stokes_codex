clear; clc; close all;

run_tag = "Sd_1.10ep03_Da_0.00ep00_gamy_p0.00ep00";

script_dir = fileparts(mfilename("fullpath"));
remesh_dir = find_remesh_dir(script_dir);
addpath(remesh_dir);
data_dir = find_data_dir(script_dir, remesh_dir);
folder = fullfile(data_dir, run_tag);

frames = list_geo_frames(folder);
if numel(frames) < 2
    error("Need at least geo0.mat and geo1.mat in %s", folder);
end

rows = [];

for idx = 2:numel(frames)
    frame_prev = frames(idx - 1);
    frame_curr = frames(idx);
    prev = load(fullfile(folder, sprintf("geo%d.mat", frame_prev)));
    curr = load(fullfile(folder, sprintf("geo%d.mat", frame_curr)));

    if ~isequal(prev.M, curr.M)
        warning("Skipping frame %d because mesh connectivity changed.", frame_curr);
        continue
    end

    geo_prev = Geometry(prev.M, prev.P);
    geo_curr = Geometry(curr.M, curr.P);
    p = curr.p;
    dt = p.dt;

    f = as_vertex_vector(curr.f, size(prev.P, 1));
    u_saved = as_vertex_vector(curr.velocity, size(prev.P, 1));

    slp_cache = stokeslet_SLP_triangle_setup(prev.M);
    slpout = stokeslet_SLP_triangle(prev.P, prev.M, f, slp_cache);
    u_background = background_velocity(prev.P, p);

    normal_slip = p.Gamma - dot(f, geo_prev.v_normal, 2);
    u_bie = u_background - p.Sd * slpout - p.Sd * p.Da * normal_slip .* geo_prev.v_normal;
    c_res = -u_saved + u_bie;

    dVdt_fd = (geo_curr.volume - geo_prev.volume) / dt;
    flux_u_prev = surface_flux(u_saved, geo_prev);
    flux_u_curr = surface_flux(u_saved, geo_curr);
    flux_slp = surface_flux(slpout, geo_prev);
    flux_background = surface_flux(u_background, geo_prev);
    flux_darcy = surface_flux(-p.Sd * p.Da * normal_slip .* geo_prev.v_normal, geo_prev);
    flux_bie = surface_flux(u_bie, geo_prev);
    c_rms = norm(c_res(:)) / sqrt(numel(c_res));

    rows = [rows; frame_curr, frame_curr * dt, geo_prev.volume, geo_curr.volume, ...
        dVdt_fd, flux_u_prev, flux_u_curr, flux_background, flux_slp, ...
        -p.Sd * flux_slp, flux_darcy, flux_bie, c_rms]; %#ok<AGROW>
end

diagnostic = array2table(rows, "VariableNames", ...
    ["frame", "time", "V_prev", "V_curr", "dVdt_fd", ...
     "flux_u_prev_geom", "flux_u_curr_geom", "flux_background", ...
     "flux_slp", "minus_Sd_flux_slp", "flux_darcy", "flux_bie", "c_rms"]);

disp("Volume-leak diagnostic for " + run_tag);
disp(diagnostic);

fprintf("\nMeans over analyzed steps:\n");
fprintf("  mean dV/dt finite difference       = %.12e\n", mean(diagnostic.dVdt_fd));
fprintf("  mean int u_saved dot n dA prev geom = %.12e\n", mean(diagnostic.flux_u_prev_geom));
fprintf("  mean -Sd int S[f] dot n dA          = %.12e\n", mean(diagnostic.minus_Sd_flux_slp));
fprintf("  mean BIE-predicted flux             = %.12e\n", mean(diagnostic.flux_bie));
fprintf("  mean BIE residual RMS               = %.12e\n", mean(diagnostic.c_rms));

figure("Color", "w");
tiledlayout(3, 1, "TileSpacing", "compact");

nexttile;
plot(diagnostic.time, diagnostic.V_curr, "o-", "LineWidth", 1.5);
grid on;
xlabel("time");
ylabel("volume");
title("Saved volume");

nexttile;
plot(diagnostic.time, diagnostic.dVdt_fd, "o-", "LineWidth", 1.5, ...
    "DisplayName", "finite-difference dV/dt");
hold on;
plot(diagnostic.time, diagnostic.flux_u_prev_geom, "s-", "LineWidth", 1.5, ...
    "DisplayName", "int u_{saved} dot n dA, previous geometry");
plot(diagnostic.time, diagnostic.flux_bie, ".-", "LineWidth", 1.5, ...
    "DisplayName", "int u_{BIE} dot n dA, previous geometry");
grid on;
xlabel("time");
ylabel("volume rate");
legend("Location", "best");

nexttile;
plot(diagnostic.time, diagnostic.minus_Sd_flux_slp, "o-", "LineWidth", 1.5, ...
    "DisplayName", "-Sd int S[f] dot n dA");
hold on;
plot(diagnostic.time, diagnostic.c_rms, "s-", "LineWidth", 1.5, ...
    "DisplayName", "BIE residual RMS");
grid on;
xlabel("time");
ylabel("diagnostic");
legend("Location", "best");

function flux = surface_flux(u, geo)
    flux = sum(dot(u, geo.v_normal, 2) .* geo.v_area);
end

function u = as_vertex_vector(value, n_v)
    if isempty(value)
        u = zeros(n_v, 3);
        return
    end
    if isequal(size(value), [n_v, 3])
        u = value;
    else
        u = reshape(value, [], 3);
    end
end

function u_background = background_velocity(P, p)
    if ~isfield(p, "gamy") || p.gamy == 0
        u_background = zeros(size(P));
    elseif exist("shearextensionflow", "file") == 2
        u_background = shearextensionflow(P, p.gamy);
    elseif exist("shear_flow", "file") == 2
        u_background = shear_flow(P, p.gamy);
    else
        error("No recognized background-flow function is on the MATLAB path.");
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

function frames = list_geo_frames(folder)
    files = dir(fullfile(folder, "geo*.mat"));
    frames = zeros(numel(files), 1);
    keep = false(numel(files), 1);

    for i = 1:numel(files)
        token = regexp(files(i).name, '^geo(\d+)\.mat$', 'tokens', 'once');
        if ~isempty(token)
            frames(i) = str2double(token{1});
            keep(i) = true;
        end
    end

    frames = sort(frames(keep));
end
