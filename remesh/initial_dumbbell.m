function [M, P, info] = initial_dumbbell(v, res)
%INITIAL_DUMBBELL  Tubular dumbbell vesicle mesh with rounded end bulbs.
%
%   [M, P] = initial_dumbbell(v)
%   [M, P, info] = initial_dumbbell(v, res)
%
%   v    : target reduced volume, 0 < v <= 1
%   res  : scalar or [nTheta nPhi]
%
%   Output:
%       M : nFaces x 3 triangle connectivity
%       P : nVerts x 3 point matrix
%
%   The surface is scaled to area A = 4*pi, so the area-radius is a = 1.
%
%   Shape idea:
%       x(theta)   = L cos(theta)
%       rho(theta) = sin(theta) * [rTube + (rBulb-rTube)*|cos(theta)|^p]
%
%   Large p gives a more cylindrical middle section and gradual end bulges.

    if nargin < 2 || isempty(res)
        nTheta = 60; %50, 50 gave good results before
        nPhi   = 60;
    elseif isscalar(res)
        nTheta = max(24, round(res));
        nPhi   = max(48, 2*nTheta);
    else
        nTheta = max(24, round(res(1)));
        nPhi   = max(48, round(res(2)));
    end

    if ~(isscalar(v) && isfinite(v) && v > 0 && v <= 1)
        error('Reduced volume v must satisfy 0 < v <= 1.');
    end

    % q controls deflation/tubularity. Larger q -> longer tube, smaller v.
    qlo = 0;
    qhi = 1;

    vlo = reducedVolumeForQ(qlo, nTheta, nPhi);
    if v >= vlo - 1e-5
        q = 0;
    else
        vhi = reducedVolumeForQ(qhi, nTheta, nPhi);

        while vhi > v
            qhi = 2*qhi;
            vhi = reducedVolumeForQ(qhi, nTheta, nPhi);

            if qhi > 1e4
                error('Could not reach target reduced volume with this shape family.');
            end
        end

        for k = 1:80
            qmid = 0.5*(qlo + qhi);
            vmid = reducedVolumeForQ(qmid, nTheta, nPhi);

            if vmid > v
                qlo = qmid;
            else
                qhi = qmid;
            end
        end

        q = 0.5*(qlo + qhi);
    end

    [M, P] = makeTubularDumbbellMesh(q, nTheta, nPhi);

    % Scale to area A = 4*pi.
    [A, ~, ~] = meshAreaVolumeReduced(P, M);
    P = P * sqrt(4*pi/A);

    [A, V, vred] = meshAreaVolumeReduced(P, M);

    if nargout > 2
        info = struct();
        info.targetReducedVolume = v;
        info.reducedVolume       = vred;
        info.area                = A;
        info.volume              = V;
        info.areaRadius          = sqrt(A/(4*pi));
        info.shapeParameterQ     = q;
        info.nTheta              = nTheta;
        info.nPhi                = nPhi;
    end
end


function vred = reducedVolumeForQ(q, nTheta, nPhi)
    [M, P] = makeTubularDumbbellMesh(q, nTheta, nPhi);
    [~, ~, vred] = meshAreaVolumeReduced(P, M);
end


function [M, P] = makeTubularDumbbellMesh(q, nTheta, nPhi)

    % Tunable shape family.
    %
    % q = 0 gives sphere-like.
    % increasing q gives:
    %   longer vesicle,
    %   thinner central tube,
    %   rounder/larger end bulbs.
    %
    % Adjust these constants if you want more/less bulbous ends.
    L     = 1 + 4.5*q;
    rTube = 1 ./ (1 + 1.45*q);
    rBulb = 1 + 0.20*q ./ (1 + 0.25*q);

    % Higher p -> longer nearly cylindrical middle.
    p = 6;

    theta = linspace(0, pi, nTheta).';
    phi   = linspace(0, 2*pi, nPhi+1);
    phi(end) = [];

    c = cos(theta);
    s = sin(theta);

    x = L * c;

    % Smooth radius: central tube radius rTube, gradually bulged ends.
    shapeRadius = rTube + (rBulb - rTube) * abs(c).^p;
    rho = s .* shapeRadius;

    % Slight extra smoothing toward poles to avoid visually flat caps.
    % This keeps rho = 0 at the poles but makes the bulbs rounder.
    poleRound = 1 + 0.10*q .* abs(c).^2 .* (1 - abs(c)).^2;
    rho = rho .* poleRound;

    nInterior = nTheta - 2;
    nVerts = 2 + nInterior*nPhi;
    P = zeros(nVerts, 3);

    % +x pole
    P(1,:) = [x(1), 0, 0];

    idx = @(i,j) 1 + (i-2)*nPhi + j;

    for i = 2:nTheta-1
        for j = 1:nPhi
            P(idx(i,j),:) = [ ...
                x(i), ...
                rho(i)*cos(phi(j)), ...
                rho(i)*sin(phi(j)) ...
            ];
        end
    end

    % -x pole
    bottom = nVerts;
    P(bottom,:) = [x(end), 0, 0];

    nFaces = 2*nPhi + 2*(nInterior-1)*nPhi;
    M = zeros(nFaces, 3);
    f = 0;

    wrap = @(j) mod(j-1, nPhi) + 1;

    % Top cap.
    for j = 1:nPhi
        jp = wrap(j+1);
        f = f + 1;
        M(f,:) = [1, idx(2,j), idx(2,jp)];
    end

    % Interior.
    for i = 2:nTheta-2
        for j = 1:nPhi
            jp = wrap(j+1);

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
        jp = wrap(j+1);
        f = f + 1;
        M(f,:) = [idx(nTheta-1,j), bottom, idx(nTheta-1,jp)];
    end

    M = M(1:f,:);
end


function [A, V, vred] = meshAreaVolumeReduced(P, M)
    p1 = P(M(:,1),:);
    p2 = P(M(:,2),:);
    p3 = P(M(:,3),:);

    cr = cross(p2 - p1, p3 - p1, 2);
    triArea = 0.5 * sqrt(sum(cr.^2, 2));
    A = sum(triArea);

    Vsigned = sum(dot(p1, cross(p2, p3, 2), 2)) / 6;
    V = abs(Vsigned);

    a = sqrt(A/(4*pi));
    vred = 3*V/(4*pi*a^3);
end