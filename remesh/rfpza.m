function [P_rot, R] = rfpza(P)
%ROTATEFIRSTPOINTTOZAXIS Rotate points so the first point lies on the z-axis.
%
% Inputs:
%   P     - Nx3 array of 3D points
%
% Outputs:
%   P_rot - Nx3 rotated points
%   R     - 3x3 rotation matrix
%
% The rotation is the minimal rotation that maps P(1,:) onto the z-axis.

    % Validate input
    if size(P,2) ~= 3
        error('Input P must be an Nx3 array.');
    end

    v = P(1,:).';   % first point as column vector

    if norm(v) == 0
        error('The first point is the zero vector, so its direction is undefined.');
    end

    v = v / norm(v);              % normalize
    z = [0; 0; 1];                % target direction

    % Case 1: already on +z axis
    if norm(v - z) < 1e-12
        R = eye(3);

    % Case 2: on -z axis
    elseif norm(v + z) < 1e-12
        % 180-degree rotation about x-axis (could also use y-axis)
        R = [1  0  0;
             0 -1  0;
             0  0 -1];

    % General case: Rodrigues' rotation formula
    else
        k = cross(v, z);          % rotation axis
        s = norm(k);
        c = dot(v, z);

        K = [   0   -k(3)  k(2);
              k(3)    0   -k(1);
             -k(2)  k(1)    0   ];

        R = eye(3) + K + K*K * ((1 - c) / s^2);
    end

    % Rotate all points
    P_rot = (R * P.').';
end