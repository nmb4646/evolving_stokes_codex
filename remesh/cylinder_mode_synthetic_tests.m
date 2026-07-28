function validation = cylinder_mode_synthetic_tests()
%CYLINDER_MODE_SYNTHETIC_TESTS End-to-end validation of the cylinder pipeline.

this_dir = fileparts(mfilename("fullpath"));
addpath(this_dir);
root = fullfile(this_dir, "data", "cylinder_mode_synthetic_validation");
series_root = fullfile(root, "series");
analysis_root = fullfile(root, "analysis");
if isfolder(root)
    rmdir(root, "s");
end
mkdir(series_root);
mkdir(analysis_root);

times_growth = [0, 0.13, 0.29, 0.48, 0.70, 0.95, ...
    1.23, 1.55, 1.90, 2.30, 2.75, 3.20].';
times_decay = [0, 0.16, 0.35, 0.57, 0.82, 1.10, ...
    1.42, 1.78, 2.18, 2.62, 3.10, 3.62, 4.18].';

static_spec = synthetic_spec("cylinder_static", times_growth(1:6));
static_spec.modes = struct([]);
static_spec.rotation_rate = 0.12;
static_spec.translation_scale = 0.3;
create_synthetic_series(series_root, static_spec);

axis_spec = synthetic_spec("cylinder_axisymmetric_growth", times_growth);
axis_spec.modes = mode_spec(0, 3, 1.5e-3, 0.35, 0.41);
axis_spec.axial_stretch_rate = 0.018;
axis_spec.rotation_rate = 0.08;
axis_spec.translation_scale = 0.25;
create_synthetic_series(series_root, axis_spec);

helical_spec = synthetic_spec("cylinder_helical_decay", times_decay);
helical_spec.modes = mode_spec(2, 2, 7e-3, -0.28, -0.63);
helical_spec.rotation_rate = 0.11;
helical_spec.translation_scale = 0.2;
create_synthetic_series(series_root, helical_spec);

remesh_spec = synthetic_spec("cylinder_remeshed_growth", times_growth);
remesh_spec.modes = mode_spec(0, 2, 1.8e-3, 0.22, 0.27);
remesh_spec.change_connectivity = true;
remesh_spec.axial_stretch_rate = 0.012;
create_synthetic_series(series_root, remesh_spec);

plus_spec = synthetic_spec("cylinder_pair_plus", times_growth);
plus_spec.modes = [ ...
    mode_spec(2, 1, 3e-3, -0.10, 0.17), ...
    mode_spec(0, 3, 1.2e-3, 0.31, -0.22)];
minus_spec = plus_spec;
minus_spec.name = "cylinder_pair_minus";
minus_spec.modes(2).amplitude0 = -minus_spec.modes(2).amplitude0;
create_synthetic_series(series_root, plus_spec);
create_synthetic_series(series_root, minus_spec);

common = test_config(series_root, analysis_root);

cfg = common;
cfg.output_root = fullfile(analysis_root, "cylinder_static");
cfg.modes.m_max = 0;
cfg.modes.n_max = 6;
static_run = cylinder_mode_pipeline("cylinder_static", cfg);
static_result = static_run.series_results{1};

cfg = common;
cfg.output_root = fullfile(analysis_root, "cylinder_axisymmetric_growth");
cfg.modes.m_max = 0;
cfg.modes.n_max = 6;
axis_run = cylinder_mode_pipeline("cylinder_axisymmetric_growth", cfg);
axis_result = axis_run.series_results{1};

cfg = common;
cfg.output_root = fullfile(analysis_root, "cylinder_helical_decay");
cfg.modes.m_max = 2;
cfg.modes.n_max = 4;
helical_run = cylinder_mode_pipeline("cylinder_helical_decay", cfg);
helical_result = helical_run.series_results{1};

cfg = common;
cfg.output_root = fullfile(analysis_root, "cylinder_remeshed_growth");
cfg.modes.m_max = 0;
cfg.modes.n_max = 5;
remesh_run = cylinder_mode_pipeline("cylinder_remeshed_growth", cfg);
remesh_result = remesh_run.series_results{1};

cfg = common;
cfg.output_root = fullfile(analysis_root, "cylinder_pair_response");
cfg.modes.m_max = 2;
cfg.modes.n_max = 5;
cfg.pairing.enabled = true;
cfg.pairing.roles = ["plus"; "minus"];
cfg.pairing.group_ids = ["validation"; "validation"];
cfg.pairing.epsilon = 1;
pair_run = cylinder_mode_pipeline( ...
    ["cylinder_pair_plus"; "cylinder_pair_minus"], cfg);
pair_result = pair_run.response_results{1};

static_amplitude = maximum_nonbase_amplitude(static_result);
axis_rate = extract_rate(axis_result, 0, 3);
helical_rate = extract_rate(helical_result, 2, 2);
remesh_rate = extract_rate(remesh_result, 0, 2);
pair_rate = extract_rate(pair_result, 0, 3);

opposite = mode_column(helical_result, 2, -2);
target = mode_column(helical_result, 2, 2);
helicity_leakage = max(helical_result.group_amplitude(:, opposite)) ...
    / max(helical_result.group_amplitude(:, target));
core_length_variation = range(remesh_result.frame_metrics.core_length) ...
    / mean(remesh_result.frame_metrics.core_length);
pair_common = mode_column(pair_result, 2, 1);
plus_common = mode_column(pair_run.series_results{1}, 2, 1);
pair_cancellation = max(pair_result.group_amplitude(:, pair_common)) ...
    / max(pair_run.series_results{1}.group_amplitude(:, plus_common));

test_name = [ ...
    "Static cylinder noise";
    "Axisymmetric growth rate";
    "Helical decay rate";
    "Opposite-helicity separation";
    "Remeshed growth rate";
    "Remeshed core stability";
    "Plus/minus growth rate";
    "Plus/minus common-mode cancellation"];
expected = [0; 0.35; -0.28; 0; 0.22; 0; 0.31; 0];
measured = [static_amplitude; axis_rate; helical_rate; helicity_leakage; ...
    remesh_rate; core_length_variation; pair_rate; pair_cancellation];
tolerance = [5e-8; 0.03; 0.03; 0.08; 0.04; 0.08; 0.04; 0.08];
absolute_error = abs(measured - expected);
passed = absolute_error <= tolerance;

validation = table(test_name, expected, measured, tolerance, ...
    absolute_error, passed);
writetable(validation, fullfile(root, "cylinder_synthetic_validation.csv"));
save(fullfile(root, "cylinder_synthetic_validation.mat"), "validation");
save_validation_figure(validation, root);

disp(validation);
if ~all(passed)
    failed = join(test_name(~passed), ", ");
    error("CylinderMode:SyntheticValidation", ...
        "Synthetic cylinder validation failed: %s", failed);
end
fprintf("All synthetic cylinder-mode validations passed.\n");
end

function cfg = test_config(series_root, analysis_root)
    cfg = cylinder_mode_defaults();
    cfg.data_root = series_root;
    cfg.output_root = analysis_root;
    cfg.verbose = false;
    cfg.continue_on_error = false;
    cfg.modes.projection = "weighted_lstsq";
    cfg.core.method = "auto";
    cfg.growth.minimum_points = 6;
    cfg.growth.maximum_points = 20;
    cfg.growth.maximum_early_frames = 30;
    cfg.growth.minimum_amplitude_ratio = 1.15;
    cfg.growth.minimum_signal_to_noise = 2.0;
    cfg.growth.maximum_dimensionless_amplitude = 0.10;
    cfg.growth.minimum_r_squared = 0.97;
    cfg.growth.maximum_relative_qR_change = 0.25;
    cfg.growth.maximum_endpoint_slope_change = 0.40;
    cfg.quality.maximum_reconstruction_error = 0.5;
    cfg.diagnostics.figure_visible = "off";
    cfg.diagnostics.resolution = 150;
    cfg.diagnostics.save_individual_fit_plots = false;
    cfg.output.save_mode_csv = false;
end

function spec = synthetic_spec(name, times)
    spec = struct();
    spec.name = string(name);
    spec.times = times(:);
    spec.radius = 0.5;
    spec.cylinder_length = 8.0;
    spec.n_axial = 48;
    spec.n_azimuthal = 32;
    spec.n_cap = 10;
    spec.modes = struct([]);
    spec.axial_stretch_rate = 0;
    spec.rotation_rate = 0;
    spec.translation_scale = 0;
    spec.change_connectivity = false;
end

function mode = mode_spec(m, n, amplitude0, growth_rate, phase)
    mode = struct( ...
        "m", m, ...
        "n", n, ...
        "amplitude0", amplitude0, ...
        "growth_rate", growth_rate, ...
        "phase", phase);
end

function create_synthetic_series(root, spec)
    folder = fullfile(root, spec.name);
    if isfolder(folder)
        rmdir(folder, "s");
    end
    mkdir(folder);

    for k = 1:numel(spec.times)
        time = spec.times(k); %#ok<NASGU>
        if spec.change_connectivity && mod(k, 2) == 0
            n_axial = spec.n_axial + 3;
            n_azimuthal = spec.n_azimuthal + 4;
        else
            n_axial = spec.n_axial;
            n_azimuthal = spec.n_azimuthal;
        end
        current_length = spec.cylinder_length * (1 + spec.axial_stretch_rate * time);
        [P, M] = capped_cylinder_mesh( ...
            spec.radius, current_length, n_axial, n_azimuthal, ...
            spec.n_cap, spec.modes, time);

        if spec.rotation_rate ~= 0
            rotation_axis = [0.4, -0.7, 0.58];
            rotation_axis = rotation_axis / norm(rotation_axis);
            rotation = axis_angle_rotation(rotation_axis, spec.rotation_rate * time);
            P = P * rotation.';
        end
        if spec.translation_scale ~= 0
            translation = spec.translation_scale * ...
                [sin(0.7 * time), cos(0.4 * time), sin(0.3 * time + 0.2)];
            P = P + translation;
        end

        p = struct(); %#ok<NASGU>
        if numel(spec.times) > 1
            p.dt = median(diff(spec.times));
        else
            p.dt = 1;
        end
        save(fullfile(folder, sprintf("geo%d.mat", k - 1)), ...
            "P", "M", "p", "time");
    end
end

function [P, M] = capped_cylinder_mesh( ...
        radius, length_cylinder, n_axial, n_phi, n_cap, modes, time)
    left_angles = (1:(n_cap - 1)) / n_cap * (pi / 2);
    left_x = -length_cylinder / 2 - radius * cos(left_angles);
    left_r = radius * sin(left_angles);

    center_x = linspace(-length_cylinder / 2, length_cylinder / 2, n_axial);
    center_r = radius * ones(size(center_x));

    right_angles = ((n_cap - 1):-1:1) / n_cap * (pi / 2);
    right_x = length_cylinder / 2 + radius * cos(right_angles);
    right_r = radius * sin(right_angles);

    ring_x = [left_x, center_x, right_x];
    ring_r = [left_r, center_r, right_r];
    n_rings = numel(ring_x);
    P = zeros(2 + n_rings * n_phi, 3);
    P(1, :) = [-length_cylinder / 2 - radius, 0, 0];

    for ring = 1:n_rings
        phi = 2 * pi * (0:(n_phi - 1)) / n_phi;
        phi = phi + 0.08 * sin(0.73 * ring) / n_phi;
        perturbation = zeros(size(phi));
        if abs(ring_x(ring)) <= length_cylinder / 2 + 10 * eps
            xi = ring_x(ring) / length_cylinder;
            for mode = modes
                amplitude = mode.amplitude0 * exp(mode.growth_rate * time);
                perturbation = perturbation + amplitude * cos( ...
                    mode.m * phi + 2 * pi * mode.n * xi + mode.phase);
            end
        end
        local_radius = ring_r(ring) .* (1 + perturbation);
        indices = 1 + (ring - 1) * n_phi + (1:n_phi);
        P(indices, 1) = ring_x(ring);
        P(indices, 2) = local_radius .* cos(phi);
        P(indices, 3) = local_radius .* sin(phi);
    end
    right_tip = size(P, 1);
    P(right_tip, :) = [length_cylinder / 2 + radius, 0, 0];

    M = zeros(2 * n_phi + 2 * n_phi * (n_rings - 1), 3);
    face = 0;
    first_ring = 1 + (1:n_phi);
    for j = 1:n_phi
        next = mod(j, n_phi) + 1;
        face = face + 1;
        M(face, :) = [1, first_ring(next), first_ring(j)];
    end

    for ring = 1:(n_rings - 1)
        lower = 1 + (ring - 1) * n_phi + (1:n_phi);
        upper = 1 + ring * n_phi + (1:n_phi);
        for j = 1:n_phi
            next = mod(j, n_phi) + 1;
            face = face + 1;
            M(face, :) = [lower(j), lower(next), upper(next)];
            face = face + 1;
            M(face, :) = [lower(j), upper(next), upper(j)];
        end
    end

    last_ring = 1 + (n_rings - 1) * n_phi + (1:n_phi);
    for j = 1:n_phi
        next = mod(j, n_phi) + 1;
        face = face + 1;
        M(face, :) = [last_ring(j), last_ring(next), right_tip];
    end
    M = M(1:face, :);
end

function rotation = axis_angle_rotation(axis, angle)
    K = [ ...
        0, -axis(3), axis(2);
        axis(3), 0, -axis(1);
        -axis(2), axis(1), 0];
    rotation = eye(3) + sin(angle) * K + (1 - cos(angle)) * (K * K);
end

function rate = extract_rate(result, m, n)
    row = result.growth_rates.m == m & result.growth_rates.n == n;
    if nnz(row) ~= 1
        error("CylinderMode:ValidationMode", ...
            "Could not identify growth-rate row m=%d, n=%d.", m, n);
    end
    rate = result.growth_rates.growth_rate_sigma(row);
end

function column = mode_column(result, m, n)
    column = find(result.group_m == m & result.group_n == n, 1);
    if isempty(column)
        error("CylinderMode:ValidationMode", ...
            "Could not identify mode m=%d, n=%d.", m, n);
    end
end

function value = maximum_nonbase_amplitude(result)
    nonbase = ~(result.group_m == 0 & result.group_n == 0);
    value = max(result.group_amplitude(:, nonbase), [], "all");
end

function save_validation_figure(validation, root)
    fig = figure("Visible", "off", "Color", "w", ...
        "Position", [100, 100, 1200, 620]);
    ax = axes(fig);
    relative_scale = max(abs(validation.expected), validation.tolerance);
    normalized_error = validation.absolute_error ./ relative_scale;
    colors = repmat([0.12, 0.55, 0.28], height(validation), 1);
    colors(~validation.passed, :) = repmat([0.82, 0.18, 0.14], ...
        nnz(~validation.passed), 1);
    bars = bar(ax, normalized_error, "FaceColor", "flat");
    bars.CData = colors;
    yline(ax, 1, "k--", "Tolerance");
    set(ax, ...
        "XTick", 1:height(validation), ...
        "XTickLabel", validation.test_name, ...
        "XTickLabelRotation", 25, ...
        "FontSize", 12, ...
        "Box", "on");
    ylabel(ax, "Error / tolerance");
    title(ax, "Synthetic cylindrical-mode validation");
    grid(ax, "on");
    exportgraphics(fig, ...
        fullfile(root, "cylinder_synthetic_validation.png"), ...
        "Resolution", 200);
    close(fig);
end
