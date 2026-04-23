function [face, uv, count, fail] = project(V, F, points, kdtree, n_neighbor)
    % project the points onto the mesh and return the face index and uv
    % inputs:
    %   V: vertex matrix
    %   F: face matrix
    %   points: N by 3, point coordinate
    %   kdtree: KDTreeSearcher
    %   n_neighbor: int, number of closest faces to search
    % outputs:
    %   face: N by 1, face index
    %   uv: N by 2, point uv coordinate
    %   count: N by 1, number of faces that contain the projected point
    %   fail: bool, if the projection failed
    
    f_list =knnsearch(kdtree, points, 'K', n_neighbor);
    [bary, distance] = bary_coord(V, F, permute(points(:, :, ones(n_neighbor, 1)), [1, 3, 2]), f_list);
    %%% triangles that contain projected point
    tol = 1e-1; up = 1. + tol; down = 0. - tol; % adjust the tolerance for uv
    is_inside = (bary(:, :, 1) < up) & (bary(:, :, 2) < up) & ...
                ((1 - bary(:, :, 1) - bary(:, :, 2)) < up) & ...
                (bary(:, :, 1) > down) & (bary(:, :, 2) > down) & ...
                ((1 - bary(:, :, 1) - bary(:, :, 2)) > down);
    
    count = sum(is_inside, 2);
    face = zeros(size(points, 1), 1);
    uv = zeros(size(points, 1), 2);
    if all(count > 0)
        %%% the closest one from is_inside triangles
        distance(~is_inside) = Inf;
        [~, min_ind] = min(distance, [], 2); 
        min_ind = sub2ind(size(distance), 1:size(distance, 1), min_ind');
        face = f_list(min_ind);
        uv = reshape(bary, [], 2);
        uv = uv(min_ind, :);
        fail = false;
    else      
        fail = true;
        disp("Closest points not exist. " + ...
            "Increase search range K. Check variable 'count'. ")
    end    
end


% %%
% close all; clear all;
% load("bob.mat")
% M = [M(:, 3), M(:, 2), M(:, 1)];
% geo = Geometry(M, P);
% kdtree = KDTreeSearcher(geo.f_center);
% a = linspace(-0.3, 0.3, 100)';
% v = geo.f_center(1,:) - dot(geo.f_center(1,:), geo.f_normal(1, :), 2) * geo.f_normal(1, :);
% p = repmat(geo.f_center(1,:), size(a, 1), 1) + a .* repmat(v, size(a, 1), 1);
% [face, uv, count, fail] = project(P, M, p, kdtree, 20);
% if fail
%     error("projection failed.");
% end
% projected_p = interpolate(M, face, uv, geo.V);
% IO.scatter(projected_p + 0 * geo.f_normal(face, :));