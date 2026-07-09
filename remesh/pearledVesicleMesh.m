function [P,M,rz] = pearledVesicleMesh(q, nu, nt)
% budMeshSmooth(q) makes two nearly spherical bulbs smoothly joined by a neck.
%
% q  = upper bulb radius / lower bulb radius
% P  = vertices, #V x 3
% M  = triangles, #F x 3
% rz = meridian curve [r z]
%
% Example:
%   [P,M] = budMeshSmooth(0.55);
%   trisurf(M,P(:,1),P(:,2),P(:,3), ...
%       'EdgeColor','none','FaceColor',[0.55 0.8 0.45]);
%   axis equal off
%   camlight
%   lighting gouraud

if nargin < 1 || isempty(q),  q  = 0.55; end
if nargin < 2 || isempty(nu), nu = 220;  end
if nargin < 3 || isempty(nt), nt = 180;  end

if q <= 0
    error('q must be positive.');
end

% Only exposed shape parameter is q.
R1 = 1.0;
R2 = q;

% Internal automatic geometry.
joinFrac   = 0.68;   % where we leave each sphere
neckFrac   = 0.24;   % neck radius relative to smaller sphere
blendScale = 1.0;   % larger = smoother, longer neck transition

rA = joinFrac * R1;
rB = joinFrac * R2;
rN = neckFrac * min(R1,R2);

% Lower sphere centered at z = 0.
zA = sqrt(R1^2 - rA^2);
mA = -zA / rA;             % dr/dz
kA = -R1^2 / rA^3;         % d2r/dz2

% Upper sphere center is chosen automatically.
dzB = sqrt(R2^2 - rB^2);
mB  = dzB / rB;
kB  = -R2^2 / rB^3;

% Blend lengths chosen automatically.
hA = max(0.25*R1, blendScale*(rA-rN)/abs(mA));
hB = max(0.25*R2, blendScale*(rB-rN)/abs(mB));

zN = zA + hA;      % neck center
zB = zN + hB;      % upper sphere join
c2 = zB + dzB;     % upper sphere center

% Positive second derivative at the neck: makes r(z) have a smooth minimum.
kN = 1.2 * min((rA-rN)/hA^2, (rB-rN)/hB^2);

% Dense meridian construction first.
nd = max(1200, 8*nu);

n1 = round(0.30*nd);
n2 = round(0.20*nd);
n3 = round(0.20*nd);
n4 = nd - n1 - n2 - n3;

% Lower exact spherical bulb.
phiA = acos(zA/R1);
phi1 = linspace(pi, phiA, n1).';
r1 = R1*sin(phi1);
z1 = R1*cos(phi1);

% C2 blend: lower sphere -> neck.
z2 = linspace(zA, zN, n2).';
r2 = quinticHermite1d(z2, zA, zN, rA, rN, mA, 0, kA, kN);

% C2 blend: neck -> upper sphere.
z3 = linspace(zN, zB, n3).';
r3 = quinticHermite1d(z3, zN, zB, rN, rB, 0, mB, kN, kB);

% Upper exact spherical bulb.
phiB = acos((zB-c2)/R2);
phi4 = linspace(phiB, 0, n4).';
r4 = R2*sin(phi4);
z4 = c2 + R2*cos(phi4);

% Merge, removing duplicate endpoints.
r0 = [r1; r2(2:end); r3(2:end); r4(2:end)];
z0 = [z1; z2(2:end); z3(2:end); z4(2:end)];

% Resample by meridian arclength to avoid uneven/jagged triangle bands.
s = [0; cumsum(sqrt(diff(r0).^2 + diff(z0).^2))];
[s,ia] = unique(s, 'stable');
r0 = r0(ia);
z0 = z0(ia);

sq = linspace(0, s(end), nu).';
r = interp1(s, r0, sq, 'pchip');
z = interp1(s, z0, sq, 'pchip');

r(1) = 0;
r(end) = 0;
r = max(r, 0);

rz = [r z];

% Revolve meridian into a triangle mesh.
theta = 2*pi*(0:nt-1)'/nt;
nz = numel(z);

P = zeros(2 + (nz-2)*nt, 3);
P(1,:) = [0 0 z(1)];

for i = 2:nz-1
    ids = 2 + (i-2)*nt + (0:nt-1);
    P(ids,1) = r(i)*cos(theta);
    P(ids,2) = r(i)*sin(theta);
    P(ids,3) = z(i);
end

top = size(P,1);
P(top,:) = [0 0 z(end)];

ring = @(i,k) 2 + (i-2)*nt + mod(k-1,nt);

M = zeros(2*nt*(nz-2), 3);
f = 0;

% Bottom cap.
for k = 1:nt
    kp = k + 1;
    f = f + 1;
    M(f,:) = [1, ring(2,kp), ring(2,k)];
end

% Side bands.
for i = 2:nz-2
    for k = 1:nt
        kp = k + 1;

        a = ring(i,k);
        b = ring(i,kp);
        c = ring(i+1,k);
        d = ring(i+1,kp);

        f = f + 1;
        M(f,:) = [a b c];

        f = f + 1;
        M(f,:) = [b d c];
    end
end

% Top cap.
for k = 1:nt
    kp = k + 1;
    f = f + 1;
    M(f,:) = [top, ring(nz-1,k), ring(nz-1,kp)];
end

M = M(1:f,:);

end

function y = quinticHermite1d(x,x0,x1,y0,y1,m0,m1,k0,k1)
% Quintic Hermite interpolant matching:
% y, dy/dx, and d2y/dx2 at both endpoints.

L = x1 - x0;
t = (x - x0) / L;

H00 = 1 - 10*t.^3 + 15*t.^4 - 6*t.^5;
H10 = t - 6*t.^3 + 8*t.^4 - 3*t.^5;
H20 = 0.5*(t.^2 - 3*t.^3 + 3*t.^4 - t.^5);

H01 = 10*t.^3 - 15*t.^4 + 6*t.^5;
H11 = -4*t.^3 + 7*t.^4 - 3*t.^5;
H21 = 0.5*(t.^3 - 2*t.^4 + t.^5);

y = y0*H00 + L*m0*H10 + L^2*k0*H20 + ...
    y1*H01 + L*m1*H11 + L^2*k1*H21;
end