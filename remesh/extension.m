function u_out = extension(P, gamy)
%EXTENSION Classic incompressible uniaxial extensional flow.
%
%   u_out = extension(P, gamy)
%
%   For P = [x, y, z], returns
%       u = gamy * [-x / 2, -y / 2, z].

if size(P, 2) ~= 3
    error('P must be an Nx3 matrix.');
end

u_out = zeros(size(P));
u_out(:, 1) = -0.5 * gamy * P(:, 1);
u_out(:, 2) = -0.5 * gamy * P(:, 2);
u_out(:, 3) = gamy * P(:, 3);
end
