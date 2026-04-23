function [M, P] = cube_triangulation()
%CUBE_TRIANGULATION Triangulation of a unit cube centered at the origin.
%   [M, P] = cube_triangulation()
%
% Outputs:
%   P : 8x3 array of vertex coordinates
%   M : 12x3 array of triangular faces (indices into P)
%
% The cube has side length 1 and is centered at (0,0,0), so its
% coordinates range from -0.5 to 0.5 in x, y, and z.

    % Vertex coordinates
    P = [
        -0.5, -0.5, -0.5;  % 1
         0.5, -0.5, -0.5;  % 2
         0.5,  0.5, -0.5;  % 3
        -0.5,  0.5, -0.5;  % 4
        -0.5, -0.5,  0.5;  % 5
         0.5, -0.5,  0.5;  % 6
         0.5,  0.5,  0.5;  % 7
        -0.5,  0.5,  0.5   % 8
    ];

    % Triangular faces (2 triangles per cube face, 12 total)
    M = [
        1, 2, 3;  1, 3, 4;  % bottom  z = -0.5
        5, 7, 6;  5, 8, 7;  % top     z =  0.5
        1, 5, 6;  1, 6, 2;  % front   y = -0.5
        2, 6, 7;  2, 7, 3;  % right   x =  0.5
        3, 7, 8;  3, 8, 4;  % back    y =  0.5
        4, 8, 5;  4, 5, 1   % left    x = -0.5
    ];
end