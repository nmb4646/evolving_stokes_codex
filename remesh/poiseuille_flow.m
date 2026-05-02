function u_out = poiseuille_flow(P, gamy, channel_width)
%POISEUILLE_FLOW Parabolic channel flow at 3D points.
%
%   u_out = poiseuille_flow(P, gamy, channel_width)
%
%   The channel is centered at z = 0 with walls at
%   z = +/- channel_width / 2. The flow is in the x direction and gamy is
%   the centerline speed.

if channel_width <= 0
    error('channel_width must be positive.');
end

if size(P, 2) ~= 3
    error('P must be an Nx3 matrix.');
end

z = P(:, 3);
u_out = zeros(size(P));
u_out(:, 1) = gamy * (1 - (2 * z / channel_width).^2);
end
