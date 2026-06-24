function P = neighbor_position_smooth(M, P, n_steps, alpha)
%NEIGHBOR_POSITION_SMOOTH Smooth vertex positions with one-ring averaging.
%
%   Unlike tangential_neighbor_smooth, this keeps the full displacement
%   toward the neighbor average, including the normal component. It can
%   change the shape and area, so callers should use small alpha and
%   re-normalize constraints as needed.

    n_v = size(P, 1);
    rows = [M(:, 1); M(:, 2); M(:, 3); M(:, 2); M(:, 3); M(:, 1)];
    cols = [M(:, 2); M(:, 3); M(:, 1); M(:, 1); M(:, 2); M(:, 3)];
    A = sparse(rows, cols, 1, n_v, n_v);
    A = double(A > 0);
    deg = sum(A, 2);

    for step = 1:n_steps
        P_avg = (A * P) ./ max(deg, 1);
        P = P + alpha * (P_avg - P);
    end
end
