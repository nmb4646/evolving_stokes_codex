[P, M] = subdivided_sphere(8);

v = .8635;
stretch_factor = a_from_v(v);
P(:,3) = stretch_factor*P(:,3);
geo = Geometry(M,P);

excess_area(geo)

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