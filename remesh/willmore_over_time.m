close all; clc; clear;

% Plot Helfrich/Willmore bending energy over time for fs_multi output.

%%% Simulation selection
% Edit these values to select runs from data/fs_batch_data.
Da = 1e-1;
Sd = [1e-6];
gamy = 1.01e-5;
v = .5;

% First geo timestep to include. Use 0 to include geo0.mat.
time_start = 0;
% Last geo timestep to include. Use inf to include all later frames.
time_final = inf;
% Options: "index" plots n; "time" plots t = n * dt from saved p.dt.
time_axis_mode = "time";

%%% Energy calculation
energy_Kb = 1; % Bending modulus passed to Geometry.willmore_energy.
subtract_initial_energy = false; % If true, plot E(t) - E(0) for each run.
normalize_by_initial_energy = false; % If true, plot E(t) / E(0) for each run.

%%% Visual settings
visual.line_width = 2.0;
visual.marker = "none";
visual.marker_size = 5;
visual.figure_position = [140, 140, 900, 560];
visual.scale_time_by_gamy = false; % If true, plot t_scaled = t * gamy_current.
visual.time_scale_by_Sd = false; % If true, plot t_scaled = t / Sd_current.
visual.x_scale = "linear"; % Options: "linear" or "log".
visual.y_scale = "linear"; % Options: "linear" or "log".
visual.show_grid = false;
visual.save_figure = false;
visual.output_file = fullfile("data", "willmore_over_time.png");
visual.output_resolution = 300;

script_dir = fileparts(mfilename("fullpath"));
remesh_dir = find_remesh_dir(script_dir);
addpath(remesh_dir);

data_dir = find_data_dir(script_dir, remesh_dir);
sd_values = Sd(:).';
gamy_values = gamy(:).';
v_values = v(:).';
multiple_sd = numel(sd_values) > 1;
multiple_gamy = numel(gamy_values) > 1;
multiple_v = numel(v_values) > 1;
run_count = numel(sd_values) * numel(gamy_values) * numel(v_values);
multiple_runs = run_count > 1;

fig = figure("Color", "w", "Position", visual.figure_position);
ax = axes(fig);
hold(ax, "on");

legend_entries = strings(run_count, 1);
x_label = "";
run_idx = 0;

for sd_idx = 1:numel(sd_values)
    Sd_current = sd_values(sd_idx);

    for gamy_idx = 1:numel(gamy_values)
        gamy_current = gamy_values(gamy_idx);

        for v_idx = 1:numel(v_values)
            v_current = v_values(v_idx);
            run_idx = run_idx + 1;
            run_tag = make_run_tag(Sd_current, Da, gamy_current, v_current);
            folder = fullfile(data_dir, run_tag);

            fprintf("Loading %s\n", folder);

            [x, x_label_current, energy] = load_willmore_series( ...
                folder, energy_Kb, time_start, time_final, time_axis_mode);
            [x, x_label_current] = apply_time_scaling(x, x_label_current, ...
                gamy_current, Sd_current, visual);
            energy = apply_energy_scaling(energy, subtract_initial_energy, ...
                normalize_by_initial_energy);

            if strlength(x_label) == 0
                x_label = x_label_current;
            elseif x_label ~= x_label_current
                x_label = "Time / timestep";
            end

            plot(ax, x, energy, ...
                "LineWidth", visual.line_width, ...
                "Marker", visual.marker, ...
                "MarkerSize", visual.marker_size);

            legend_entries(run_idx) = make_legend_entry(Sd_current, gamy_current, v_current, ...
                multiple_sd, multiple_gamy, multiple_v);

            fprintf("Sd = %.8g, gamy = %.8g, v = %.8g, initial energy: %.12g, final energy: %.12g\n", ...
                Sd_current, gamy_current, v_current, energy(1), energy(end));
        end
    end
end

xlabel(ax, x_label);
ylabel(ax, energy_label(subtract_initial_energy, normalize_by_initial_energy));
set(ax, "XScale", visual.x_scale);
set(ax, "YScale", visual.y_scale);

if visual.show_grid
    grid(ax, "on");
end
box(ax, "on");

if multiple_runs
    title(ax, sprintf("Helfrich energy over time: %s", ...
        comparison_title_suffix(Da, sd_values, gamy_values, v_values)));
    legend(ax, legend_entries, "Location", "best");
else
    title(ax, sprintf("Helfrich energy over time: Sd = %.4g, Da = %.4g, gamy = %.4g, v = %.4g", ...
        sd_values, Da, gamy_values, v_values));
end

if visual.save_figure
    exportgraphics(fig, visual.output_file, "Resolution", visual.output_resolution);
    fprintf("Saved %s\n", visual.output_file);
end

function [x, x_label, energy] = load_willmore_series(folder, energy_Kb, time_start, time_final, time_axis_mode)
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

    if time_final < time_start
        error("time_final must be greater than or equal to time_start.");
    end

    keep = frame_ids >= time_start & frame_ids <= time_final;
    frame_ids = frame_ids(keep);
    files = files(keep);

    if isempty(frame_ids)
        error("No valid geo*.mat frame numbers from time_start = %.17g through time_final = %.17g found in %s", ...
            time_start, time_final, folder);
    end

    energy = zeros(numel(frame_ids), 1);
    dt = NaN;

    for i = 1:numel(frame_ids)
        frame_path = fullfile(folder, files(i).name);
        data = load(frame_path, "M", "P", "p");

        if i == 1 && isfield(data, "p") && isfield(data.p, "dt")
            dt = data.p.dt;
        end

        geo = Geometry(data.M, data.P);
        energy(i) = geo.willmore_energy(energy_Kb);
    end

    time_axis_mode = lower(string(time_axis_mode));
    switch time_axis_mode
        case {"index", "n", "step", "timestep"}
            x = frame_ids;
            x_label = "Timestep n";
        case {"time", "t", "physical_time"}
            if isnan(dt)
                error("Cannot use time_axis_mode = 'time' because p.dt was not found in %s.", folder);
            end
            x = frame_ids * dt;
            x_label = "Time t = n dt";
        otherwise
            error("Unknown time_axis_mode '%s'. Use 'index' or 'time'.", time_axis_mode);
    end
end

function energy = apply_energy_scaling(energy, subtract_initial, normalize_by_initial)
    if subtract_initial && normalize_by_initial
        error("Use only one of subtract_initial_energy or normalize_by_initial_energy.");
    end

    if subtract_initial
        energy = energy - energy(1);
    elseif normalize_by_initial
        if energy(1) == 0
            error("Cannot normalize by initial energy because E(0) is zero.");
        end
        energy = energy / energy(1);
    end
end

function label = energy_label(subtract_initial, normalize_by_initial)
    if subtract_initial
        label = "Helfrich energy E - E_0";
    elseif normalize_by_initial
        label = "Helfrich energy E / E_0";
    else
        label = "Helfrich energy";
    end
end

function [x, x_label] = apply_time_scaling(x, base_label, gamy, Sd, visual)
    scale_by_gamy = visual.scale_time_by_gamy;
    scale_by_Sd = visual.time_scale_by_Sd;

    if ~scale_by_gamy && ~scale_by_Sd
        x_label = base_label;
        return
    end

    if scale_by_Sd && Sd == 0
        error("Cannot scale time by Sd because Sd is zero.");
    end

    if scale_by_gamy
        x = x * gamy;
    end

    if scale_by_Sd
        x = x / Sd;
    end

    if startsWith(base_label, "Time")
        quantity = "time";
        symbol = "t";
    elseif startsWith(base_label, "Timestep")
        quantity = "timestep";
        symbol = "n";
    else
        quantity = "time / timestep";
        symbol = "x";
    end

    if scale_by_gamy && scale_by_Sd
        expression = "gamy " + symbol + " / Sd";
    elseif scale_by_gamy
        expression = "gamy " + symbol;
    else
        expression = symbol + " / Sd";
    end

    x_label = "Scaled " + quantity + " (" + expression + ")";
end

function run_tag = make_run_tag(Sd, Da, gamy, v)
    run_tag = sprintf("Sd_%.2e_Da_%.2e_gamy_%+.2e_v_%.2e", Sd, Da, gamy, v);
    run_tag = strrep(run_tag, "+", "p");
    run_tag = strrep(run_tag, "-", "m");
end

function label = make_legend_entry(Sd, gamy, v, show_sd, show_gamy, show_v)
    parts = strings(0, 1);

    if show_sd
        parts(end + 1) = sprintf("Sd = %.4g", Sd);
    end

    if show_gamy
        parts(end + 1) = sprintf("gamy = %.4g", gamy);
    end

    if show_v
        parts(end + 1) = sprintf("v = %.4g", v);
    end

    if isempty(parts)
        label = sprintf("Sd = %.4g, gamy = %.4g, v = %.4g", Sd, gamy, v);
    else
        label = strjoin(parts, ", ");
    end
end

function suffix = comparison_title_suffix(Da, sd_values, gamy_values, v_values)
    parts = "Da = " + sprintf("%.4g", Da);

    if isscalar(sd_values)
        parts(end + 1) = "Sd = " + sprintf("%.4g", sd_values);
    end

    if isscalar(gamy_values)
        parts(end + 1) = "gamy = " + sprintf("%.4g", gamy_values);
    end

    if isscalar(v_values)
        parts(end + 1) = "v = " + sprintf("%.4g", v_values);
    end

    suffix = strjoin(parts, ", ");
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
