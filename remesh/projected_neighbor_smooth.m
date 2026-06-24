function [P, info] = projected_neighbor_smooth(M, P, M_ref, P_ref, n_steps, alpha, n_neighbor)
%PROJECTED_NEIGHBOR_SMOOTH Smooth positions and project back to a reference mesh.
%
%   [P, info] = projected_neighbor_smooth(M, P, M_ref, P_ref, n_steps,
%   alpha, n_neighbor) repeatedly moves vertices toward their one-ring
%   neighbor averages, then projects the trial points back to the reference
%   mesh. This can improve local vertex layout without allowing free normal
%   shrinkage away from the reference surface.

    if nargin < 7 || isempty(n_neighbor)
        n_neighbor = 12;
    end

    n_v = size(P, 1);
    rows = [M(:, 1); M(:, 2); M(:, 3); M(:, 2); M(:, 3); M(:, 1)];
    cols = [M(:, 2); M(:, 3); M(:, 1); M(:, 1); M(:, 2); M(:, 3)];
    A = sparse(rows, cols, 1, n_v, n_v);
    A = double(A > 0);
    deg = sum(A, 2);

    geo_ref = Geometry(M_ref, P_ref);
    kdtree = KDTreeSearcher(geo_ref.f_center);

    info.fail = false;
    info.max_projection_distance = zeros(n_steps, 1);
    info.mean_projection_distance = zeros(n_steps, 1);

    for step = 1:n_steps
        P_avg = (A * P) ./ max(deg, 1);
        P_trial = P + alpha * (P_avg - P);

        [face, uv, ~, fail] = project(P_ref, M_ref, P_trial, kdtree, n_neighbor);
        if fail
            error("projected_neighbor_smooth: projection failed at smoothing step %d.", step);
        end

        P_projected = interpolate(M_ref, face, uv, P_ref);
        projection_distance = vecnorm(P_trial - P_projected, 2, 2);
        info.max_projection_distance(step) = max(projection_distance);
        info.mean_projection_distance(step) = mean(projection_distance);
        P = P_projected;
    end
end
