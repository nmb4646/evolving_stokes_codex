function out = vesicleTiltDeformation(P, M, flowDim, gradDim)
%VESICLETILTDEFORMATION Compute vesicle tilt angle psi and deformation D.
%
% Inputs
%   P       N x 3 vertex coordinates
%   M       T x 3 triangle connectivity
%   flowDim flow direction dimension, default 1 for x
%   gradDim shear-gradient dimension, default 3 for z
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

    if nargin < 3 || isempty(flowDim)
        flowDim = 1;   % x
    end
    if nargin < 4 || isempty(gradDim)
        gradDim = 3;   % z
    end

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

    % Axis lengths. For convex vesicles, projection width is robust and
    % usually agrees with the paper's radius-sum definition.
    LA = axisLengthProjection(Pc, eA);
    LB = axisLengthProjection(Pc, eB);

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