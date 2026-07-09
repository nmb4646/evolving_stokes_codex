function out = vesicleTiltDeformation(P, M, flowDim, gradDim, axisLengthMethod)
%VESICLETILTDEFORMATION Compute vesicle tilt angle psi and deformation D.
%
% Inputs
%   P       N x 3 vertex coordinates
%   M       T x 3 triangle connectivity
%   flowDim flow direction dimension, default 1 for x
%   gradDim shear-gradient dimension, default 3 for z
%   axisLengthMethod
%           "projection"        use max/min projection width along each axis
%           "ray_intersection"  use radius(+axis) + radius(-axis), matching
%                               the Zhao/Shaqfeh-style definition (default)
%
% For the paper's shear flow u = gamma_dot*z*e_x:
%   flowDim = 1;
%   gradDim = 3;
%
% If your simulation uses u = gamma_dot*y*e_x, use:
%   flowDim = 1;
%   gradDim = 2;
%
% Outputs in struct out:
%   out.psi          tilt angle in radians, wrapped to [-pi/2, pi/2)
%   out.psi_over_pi  psi/pi
%   out.D            deformation factor
%   out.L1           longer axis length in shear plane
%   out.L2           shorter axis length in shear plane
%   out.eLong        longer principal direction
%   out.eShort       shorter principal direction
%   out.center       centroid used for centering
%   out.I            surface inertia tensor
%   out.axis_length_method
%                    requested axis length method
%   out.axis_length_used
%                    actual method used for the two shear-plane axes

    if nargin < 3 || isempty(flowDim)
        flowDim = 1;   % x
    end
    if nargin < 4 || isempty(gradDim)
        gradDim = 3;   % z
    end
    if nargin < 5 || isempty(axisLengthMethod)
        axisLengthMethod = "ray_intersection";
    end
    axisLengthMethod = normalizeAxisLengthMethod(axisLengthMethod);

    P = double(P);
    M = double(M);

    % ---- Center the vesicle ----
    center = surfaceCentroid(P, M);
    Pc = P - center;

    % ---- Surface inertia tensor ----
    I = surfaceInertiaTensor(Pc, M);

    % Eigenvectors are principal axes.
    [V, Lambda] = eig(I);
    eigvals = diag(Lambda);

    % Store axes as columns.
    axes = V;

    % Choose the two principal axes most aligned with the shear plane.
    planeScore = axes(flowDim,:).^2 + axes(gradDim,:).^2;
    [~, idx] = sort(planeScore, 'descend');
    iA = idx(1);
    iB = idx(2);

    eA = axes(:, iA);
    eB = axes(:, iB);

    % Axis lengths. The ray-intersection method follows the radius-sum
    % definition used by Zhao/Shaqfeh; projection width is the older robust
    % approximation.
    [LA, usedA] = axisLength(Pc, M, eA, axisLengthMethod);
    [LB, usedB] = axisLength(Pc, M, eB, axisLengthMethod);

    % Pick long and short axes in the shear plane.
    if LA >= LB
        eLong = eA;
        eShort = eB;
        L1 = LA;
        L2 = LB;
    else
        eLong = eB;
        eShort = eA;
        L1 = LB;
        L2 = LA;
    end

    % Fix sign ambiguity of eigenvector.
    if eLong(flowDim) < 0
        eLong = -eLong;
    end

    % Tilt angle from flow axis toward gradient axis.
    psi = atan2(eLong(gradDim), eLong(flowDim));

    % An axis is a line, not an arrow, so wrap to [-pi/2, pi/2).
    if psi >= pi/2
        psi = psi - pi;
    elseif psi < -pi/2
        psi = psi + pi;
    end

    D = (L1 - L2) / (L1 + L2);

    out.psi = psi;
    out.psi_over_pi = psi / pi;
    out.D = D;
    out.L1 = L1;
    out.L2 = L2;
    out.eLong = eLong;
    out.eShort = eShort;
    out.center = center;
    out.I = I;
    out.eigvals = eigvals;
    out.axes = axes;
    out.axis_length_method = axisLengthMethod;
    out.axis_length_used = [usedA, usedB];
end


function center = surfaceCentroid(P, M)
% Area-weighted surface centroid.

    tri = P(M, :);
    tri = reshape(tri, [], 3, 3);

    A = squeeze(tri(:,1,:));
    B = squeeze(tri(:,2,:));
    C = squeeze(tri(:,3,:));

    normals = cross(B - A, C - A, 2);
    area = 0.5 * vecnorm(normals, 2, 2);

    faceCenter = (A + B + C) / 3;

    center = sum(faceCenter .* area, 1) / sum(area);
end


function I = surfaceInertiaTensor(P, M)
% Computes I = integral_M (|x|^2 I_3 - x x^T) dA
% using exact quadratic integration over each triangle.

    tri = P(M, :);
    tri = reshape(tri, [], 3, 3);

    A = squeeze(tri(:,1,:));
    B = squeeze(tri(:,2,:));
    C = squeeze(tri(:,3,:));

    normals = cross(B - A, C - A, 2);
    area = 0.5 * vecnorm(normals, 2, 2);

    I = zeros(3,3);

    for k = 1:size(M,1)
        a = A(k,:).';
        b = B(k,:).';
        c = C(k,:).';
        s = a + b + c;
        Ak = area(k);

        % Q = integral_T x x^T dA
        Q = Ak/12 * (a*a.' + b*b.' + c*c.' + s*s.');

        I = I + trace(Q)*eye(3) - Q;
    end
end


function L = axisLengthProjection(Pc, e)
% Width of the surface along direction e.
% For ellipsoid-like shapes this is a good estimate of r_plus + r_minus.

    e = e(:) / norm(e);
    q = Pc * e;
    L = max(q) - min(q);
end


function method = normalizeAxisLengthMethod(method)
    method = lower(string(method));
    switch method
        case {"projection", "projected", "width", "projection_width"}
            method = "projection";
        case {"ray", "ray_intersection", "radius_sum", "shaqfeh", "zhao_shaqfeh"}
            method = "ray_intersection";
        otherwise
            error("Unknown axisLengthMethod '%s'. Use 'projection' or 'ray_intersection'.", method);
    end
end


function [L, usedMethod] = axisLength(Pc, M, e, method)
    switch method
        case "projection"
            L = axisLengthProjection(Pc, e);
            usedMethod = "projection";
        case "ray_intersection"
            [rPlus, hitPlus] = raySurfaceRadius(Pc, M, e);
            [rMinus, hitMinus] = raySurfaceRadius(Pc, M, -e);
            if hitPlus && hitMinus
                L = rPlus + rMinus;
                usedMethod = "ray_intersection";
            else
                L = axisLengthProjection(Pc, e);
                usedMethod = "projection_fallback";
            end
    end
end


function [radius, hit] = raySurfaceRadius(Pc, M, e)
%RAYSURFACERADIUS First positive ray/triangle intersection distance.
% The ray starts at the centered origin and points along unit direction e.

    e = e(:).' / norm(e);

    A = Pc(M(:,1), :);
    B = Pc(M(:,2), :);
    C = Pc(M(:,3), :);
    edge1 = B - A;
    edge2 = C - A;
    n_tri = size(M, 1);

    ray = repmat(e, n_tri, 1);
    h = cross(ray, edge2, 2);
    det = dot(edge1, h, 2);

    scale = max(1, max(abs(Pc), [], "all"));
    detTol = 1e-14 * scale^2;
    baryTol = 1e-10;
    tTol = 1e-12 * scale;

    notParallel = abs(det) > detTol;
    invDet = zeros(n_tri, 1);
    invDet(notParallel) = 1 ./ det(notParallel);

    s = -A;
    u = invDet .* dot(s, h, 2);
    q = cross(s, edge1, 2);
    v = invDet .* dot(ray, q, 2);
    t = invDet .* dot(edge2, q, 2);

    valid = notParallel ...
        & u >= -baryTol ...
        & v >= -baryTol ...
        & (u + v) <= 1 + baryTol ...
        & t > tTol ...
        & isfinite(t);

    if any(valid)
        radius = min(t(valid));
        hit = true;
    else
        radius = NaN;
        hit = false;
    end
end
