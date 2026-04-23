function V = shearextensionflow(P,gamy)
%VESICLEBACKGROUNDVELOCITY Background flow velocity at 3D points.
%
%   V = vesicleBackgroundVelocity(P)
%
%   Input:
%       P : Nx3 matrix of point coordinates [x, y, z]
%
%   Output:
%       V : Nx3 matrix of velocities [vx, vy, vz]
%
%   The flow is given in cylindrical coordinates by
%       U = (5*r*z) e_theta + (-0.5*r) e_r + z e_z
%
%   Converted to Cartesian coordinates:
%       vx = -0.5*x - 5*y*z
%       vy = -0.5*y + 5*x*z
%       vz = z

    if size(P,2) ~= 3
        error('Input P must be an Nx3 matrix.');
    end

    x = P(:,1);
    y = P(:,2);
    z = P(:,3);

    vx = -0.5 .* x - 5 .* y .* z;
    vy = -0.5 .* y + 5 .* x .* z;
    vz = z;

    V = gamy*[vx, vy, vz];
end