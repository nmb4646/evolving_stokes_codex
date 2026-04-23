function [VV,FF] = subdivided_cube(N)
%SUBDIVIDED_CUBE Triangulated cube surface by recursive-style face subdivision
%
%   [VV,FF] = subdivided_cube(N)
%
% Generates a simplicial triangulation of the surface of the cube
% of side length 1 centered at the origin.
%
% Input:
%   N  - subdivision parameter (positive integer)
%
% Output:
%   VV - vertices
%   FF - triangular faces
%
% The construction mirrors subdivided_sphere(N):
% each seed triangle is subdivided into N^2 smaller triangles
% using barycentric interpolation, then duplicate vertices are merged.

if N < 1 || round(N) ~= N
    error('N must be a positive integer.');
end

% Cube vertices: side length 1, centered at origin
V = [
   -0.5 -0.5 -0.5;  % 1
    0.5 -0.5 -0.5;  % 2
    0.5  0.5 -0.5;  % 3
   -0.5  0.5 -0.5;  % 4
   -0.5 -0.5  0.5;  % 5
    0.5 -0.5  0.5;  % 6
    0.5  0.5  0.5;  % 7
   -0.5  0.5  0.5   % 8
];

% 12 triangles (2 per square face), consistently oriented outward
F = [
    1 3 2; 1 4 3;   % bottom  z = -0.5
    5 6 7; 5 7 8;   % top     z =  0.5
    1 2 6; 1 6 5;   % front   y = -0.5
    2 3 7; 2 7 6;   % right   x =  0.5
    3 4 8; 3 8 7;   % back    y =  0.5
    4 1 5; 4 5 8    % left    x = -0.5
];

VV = [];
FF = [];

% Subdivide each face triangle exactly as in subdivided_sphere
for i = 1:size(F,1)

    v1 = V(F(i,1),:);
    v2 = V(F(i,2),:);
    v3 = V(F(i,3),:);

    % Vertices generated on this face triangle
    c = 0;
    v123 = [];

    for j = 0:N
        v12 = (v1*(N-j) + v2*j)/N;
        v13 = (v1*(N-j) + v3*j)/N;

        for k = 0:j
            c = c + 1;
            if j == 0
                v123(c,:) = v12;
            else
                v123(c,:) = (v12*(j-k) + v13*k)/j;
            end
        end
    end

    % Connectivity of subdivided small triangles
    c = 0;
    FT = [];

    for j = 1:N
        for k = 1:j
            c = c + 1;
            n = j*(j-1)/2 + k;
            FT(c,:) = [n, n+j, n+j+1];
        end
    end

    for j = 2:N
        for k = 1:j-1
            c = c + 1;
            n = j*(j-1)/2 + k;
            FT(c,:) = [n, n+j+1, n+1];
        end
    end

    % Merge into global arrays
    m  = size(VV,1);
    VV = [VV; v123];
    FF = [FF; FT + m];
end

% Remove duplicate vertices
tol = 1e-10;
[VV,~,ix] = uniquetol(VV, tol, 'ByRows', true);
FF = ix(FF);

end