function [P_rot, M_rot] = rotate_vesicle(P, M, angle_over_pi, axis)
%ROTATE_VESICLE Rotate vesicle vertices around a coordinate or vector axis.
%
% Inputs
%   P             N x 3 vertex coordinates
%   M             T x 3 triangle connectivity
%   angle_over_pi rotation angle divided by pi. Example: 0.25 rotates pi/4.
%   axis          "x", "y", "z", or a finite nonzero 1 x 3 vector
%
% Outputs
%   P_rot         rotated vertex coordinates
%   M_rot         unchanged triangle connectivity

    if nargin < 4 || isempty(axis)
        axis = "y";
    end

    P = double(P);
    M_rot = M;

    axis = parse_rotation_axis(axis);
    theta = double(angle_over_pi) * pi;

    K = [0, -axis(3), axis(2); ...
         axis(3), 0, -axis(1); ...
         -axis(2), axis(1), 0];
    R = eye(3) + sin(theta) * K + (1 - cos(theta)) * (K * K);

    P_rot = P * R.';
end

function axis = parse_rotation_axis(axis)
    if isstring(axis) || ischar(axis)
        switch lower(string(axis))
            case "x"
                axis = [1, 0, 0];
            case "y"
                axis = [0, 1, 0];
            case "z"
                axis = [0, 0, 1];
            otherwise
                error('Unknown rotation axis "%s". Use "x", "y", "z", or a 1 x 3 vector.', axis);
        end
    end

    axis = double(axis(:).');
    if numel(axis) ~= 3 || any(~isfinite(axis)) || norm(axis) == 0
        error('Rotation axis must be "x", "y", "z", or a finite nonzero 1 x 3 vector.');
    end
    axis = axis / norm(axis);
end
