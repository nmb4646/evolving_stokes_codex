function [M, P, info] = initial_capsule(v, res)
%INITIAL_CAPSULE Long cylindrical vesicle with hemispherical end caps.
%
%   [M, P] = initial_capsule(v)
%   [M, P, info] = initial_capsule(v, res)
%
%   v   : target reduced volume, 0 < v <= 1
%   res : scalar or [nMeridian nPhi]
%
%   The final mesh is scaled to area A = 4*pi, so the area radius is a = 1.
%
%   The shape consists of:
%       - a cylinder of radius R and length ell
%       - two hemispherical end caps of radius R

    if nargin < 2 || isempty(res)
        nMeridian = 160;
        nPhi      = 40;
    elseif isscalar(res)
        nMeridian = max(40, round(res));
        nPhi      = max(24, round(res/4));
    else
        nMeridian = max(40, round(res(1)));
        nPhi      = max(24, round(res(2)));
    end

    if ~(isscalar(v) && isfinite(v) && v > 0 && v <= 1)
        error('Reduced volume v must satisfy 0 < v <= 1.');
    end

    % For an exact cylinder plus hemispherical caps with area 4*pi:
    %
    %     v = (3*R - R^3)/2
    %
    % The physical root lies in 0 < R <= 1.
    radiusEquation = @(R) 0.5*(3*R - R.^3) - v;

    if abs(v - 1) < 1e-12
        R = 1;
    else
        R = fzero(radiusEquation, [eps, 1]);
    end

    % Cylindrical-section length from area constraint:
    %
    %     2*pi*R*ell + 4*pi*R^2 = 4*pi
    %
    ell = 2/R - 2*R;

    [M, P] = makeCapsuleMesh(R, ell, nMeridian, nPhi);

    % Correct the discrete mesh area to exactly 4*pi.
    [A, ~, ~] = meshAreaVolumeReduced(P, M);
    scale = sqrt(4*pi/A);
    P = scale*P;

    [A, V, vred] = meshAreaVolumeReduced(P, M);

    if nargout > 2
        info = struct();
        info.targetReducedVolume = v;
        info.reducedVolume       = vred;
        info.area                = A;
        info.volume              = V;
        info.areaRadius          = sqrt(A/(4*pi));
        info.radius              = scale*R;
        info.cylinderLength      = scale*ell;
        info.totalLength         = scale*(ell + 2*R);
        info.nMeridian           = nMeridian;
        info.nPhi                = nPhi;
    end
end


function [M, P] = makeCapsuleMesh(R, ell, nMeridian, nPhi)
%MAKECAPSULEMESH Construct an axisymmetric triangulated capsule.

    % Arc length along the generating curve from the +x pole to -x pole.
    %
    % Top hemisphere:    pi*R/2
    % Cylinder:          ell
    % Bottom hemisphere: pi*R/2
    totalArcLength = pi*R + ell;

    s = linspace(0, totalArcLength, nMeridian).';

    topCapEnd = pi*R/2;
    cylinderEnd = topCapEnd + ell;

    x   = zeros(nMeridian, 1);
    rho = zeros(nMeridian, 1);

    % Top hemisphere: +x pole to the cylindrical section.
    I = s <= topCapEnd;
    alpha = s(I)/R;

    x(I)   = ell/2 + R*cos(alpha);
    rho(I) = R*sin(alpha);

    % Cylindrical section.
    I = s > topCapEnd & s < cylinderEnd;
    sc = s(I) - topCapEnd;

    x(I)   = ell/2 - sc;
    rho(I) = R;

    % Bottom hemisphere: cylindrical section to -x pole.
    I = s >= cylinderEnd;
    beta = (s(I) - cylinderEnd)/R;

    x(I)   = -ell/2 - R*sin(beta);
    rho(I) = R*cos(beta);

    rho(1)   = 0;
    rho(end) = 0;

    phi = linspace(0, 2*pi, nPhi + 1);
    phi(end) = [];

    nInterior = nMeridian - 2;
    nVerts = 2 + nInterior*nPhi;

    P = zeros(nVerts, 3);

    % +x pole.
    P(1,:) = [x(1), 0, 0];

    idx = @(i,j) 1 + (i - 2)*nPhi + j;

    % Interior rings.
    for i = 2:nMeridian-1
        for j = 1:nPhi
            P(idx(i,j),:) = [ ...
                x(i), ...
                rho(i)*cos(phi(j)), ...
                rho(i)*sin(phi(j)) ...
            ];
        end
    end

    % -x pole.
    bottom = nVerts;
    P(bottom,:) = [x(end), 0, 0];

    nFaces = 2*nPhi + 2*(nInterior - 1)*nPhi;
    M = zeros(nFaces, 3);

    wrap = @(j) mod(j - 1, nPhi) + 1;
    f = 0;

    % Top cap.
    for j = 1:nPhi
        jp = wrap(j + 1);

        f = f + 1;
        M(f,:) = [1, idx(2,j), idx(2,jp)];
    end

    % Interior ring connections.
    for i = 2:nMeridian-2
        for j = 1:nPhi
            jp = wrap(j + 1);

            v00 = idx(i,   j);
            v01 = idx(i,   jp);
            v10 = idx(i+1, j);
            v11 = idx(i+1, jp);

            f = f + 1;
            M(f,:) = [v00, v10, v11];

            f = f + 1;
            M(f,:) = [v00, v11, v01];
        end
    end

    % Bottom cap.
    for j = 1:nPhi
        jp = wrap(j + 1);

        f = f + 1;
        M(f,:) = [idx(nMeridian-1,j), bottom, idx(nMeridian-1,jp)];
    end

    M = M(1:f,:);
end


function [A, V, vred] = meshAreaVolumeReduced(P, M)

    p1 = P(M(:,1),:);
    p2 = P(M(:,2),:);
    p3 = P(M(:,3),:);

    cr = cross(p2 - p1, p3 - p1, 2);

    triArea = 0.5*sqrt(sum(cr.^2, 2));
    A = sum(triArea);

    Vsigned = sum(dot(p1, cross(p2, p3, 2), 2))/6;
    V = abs(Vsigned);

    a = sqrt(A/(4*pi));
    vred = 3*V/(4*pi*a^3);
end