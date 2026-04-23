function [uv, distance] = bary_coord(V, F, point, face)
    % project the point onto face and return its barycentric coordinate
    % 
    % Inputs:
    %   V: vertex matrix 
    %   F: face matrix
    %   point: N by M, point coordinate
    %   face: N by M, face index
    % Outputs:
    %   uv: N by M by 2, point uv coordinate
    %   distance: N by M, distance

    %%% linearize
    s = size(face);
    point = reshape(point, [], 3);
    face = face(:);
    %%% uv 
    v1 = V(F(face,1),:);
    v2 = V(F(face,2),:);
    v3 = V(F(face,3),:);
    n = cross(v2 - v1, v3 - v2, 2);
    f_area = 0.5 * sqrt(sum(n.^2, 2));
    n = 0.5 * n ./ f_area;
    distance = dot(point - v1, n, 2);
    point = point - n .* distance;
    u = dot(cross(point - v3, v2 - v3, 2), n, 2) ./ 2./ f_area;
    v = dot(cross(v1 - v3, point - v3, 2), n, 2) ./ 2./ f_area;
    %%% matrix
    uv = reshape([u, v], s(1), s(2), 2);
    distance = reshape(abs(distance), s);
end