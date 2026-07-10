close all; clc;

%%% Inputs
subdivisions = 11;

% Reduced-volume values used to numerically build excess_area(v).
% This range must cover the target excess-area values below.
v_scan = linspace(0.70, 0.9999, 300);

% Target excess-area values. Edit this to whatever range you want.
target_excess_area = 0.1:0.1:1.3;

interpolation_method = "pchip"; % Options: "linear", "pchip", "makima", etc.
make_plot = true;

%%% Numerically form excess area as a function of reduced volume.
[P0, M] = subdivided_sphere(subdivisions);

v_scan = v_scan(:);
stretch_factor = zeros(size(v_scan));
numeric_reduced_volume = zeros(size(v_scan));
numeric_excess_area = zeros(size(v_scan));

for i = 1:numel(v_scan)
    stretch_factor(i) = a_from_v(v_scan(i));

    P = P0;
    P(:, 3) = stretch_factor(i) * P(:, 3);
    geo = Geometry(M, P);

    numeric_reduced_volume(i) = reduced_volume(geo);
    numeric_excess_area(i) = excess_area(geo);
end

%%% Interpolate reduced volume and stretch factor for the requested excess areas.
[excess_sorted, order] = sort(numeric_excess_area);
v_sorted = v_scan(order);
stretch_sorted = stretch_factor(order);
numeric_v_sorted = numeric_reduced_volume(order);

[excess_unique, unique_idx] = unique(excess_sorted, "stable");
v_unique = v_sorted(unique_idx);
stretch_unique = stretch_sorted(unique_idx);
numeric_v_unique = numeric_v_sorted(unique_idx);

target_excess_area = target_excess_area(:);
outside = target_excess_area < min(excess_unique) | target_excess_area > max(excess_unique);
if any(outside)
    warning("Some target excess areas are outside the scanned range [%.6g, %.6g]; those rows will be NaN.", ...
        min(excess_unique), max(excess_unique));
end

target_v = interp1(excess_unique, v_unique, target_excess_area, interpolation_method, NaN);
target_stretch_factor = interp1(excess_unique, stretch_unique, target_excess_area, interpolation_method, NaN);
target_numeric_reduced_volume = interp1(excess_unique, numeric_v_unique, target_excess_area, interpolation_method, NaN);

target_table = table(target_excess_area, target_v, target_numeric_reduced_volume, target_stretch_factor, ...
    'VariableNames', ["target_excess_area", "input_v", "mesh_reduced_volume", "stretch_factor"]);

scan_table = table(v_scan, numeric_reduced_volume, numeric_excess_area, stretch_factor, ...
    'VariableNames', ["input_v", "mesh_reduced_volume", "mesh_excess_area", "stretch_factor"]);

fprintf("Scanned excess-area range: %.8g to %.8g\n", min(numeric_excess_area), max(numeric_excess_area));
disp(target_table);

if make_plot
    figure("Color", "w");
    plot(v_scan, numeric_excess_area, "LineWidth", 1.5);
    hold on;
    scatter(target_v, target_excess_area, 45, "filled");
    grid on;
    box on;
    xlabel("Input reduced volume v");
    ylabel("Numerical excess area");
    title("Numerical excess area vs reduced volume");
    legend("scan", "interpolated targets", "Location", "best");
end

function a = a_from_v(v)
% a_from_v  Prolate spheroid aspect ratio from reduced volume.
%
%   a = a_from_v(v)
%
% Returns the aspect ratio a = c/b >= 1 for a prolate spheroid whose
% reduced volume is v, using
%
%   v = V / ((4*pi/3) * R_A^3),
%   R_A = sqrt(A/(4*pi)).
%
% Valid input:
%   0 < v <= 1
%
% For v = 1, the result is a = 1, the sphere.

    if v <= 0 || v > 1
        error('v must satisfy 0 < v <= 1.');
    end

    if abs(v - 1) < 1e-12
        a = 1;
        return
    end

    f = @(a) reduced_volume_prolate(a) - v;

    % Find an upper bracket where reduced_volume(a) < v
    lo = 1;
    hi = 2;

    while f(hi) > 0
        hi = 2 * hi;
    end

    % Solve by bisection/fzero with bracket
    a = fzero(f, [lo hi]);
end


function v = reduced_volume_prolate(a)
% reduced_volume  Reduced volume of prolate spheroid with aspect ratio a.

    if abs(a - 1) < 1e-12
        v = 1;
        return
    end

    e = sqrt(1 - 1/a^2);

    F = 1 + (a/e) * asin(e);

    v = 2 * sqrt(2) * a / F^(3/2);
end
