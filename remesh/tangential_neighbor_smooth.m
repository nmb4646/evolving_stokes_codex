function P = tangential_neighbor_smooth(M, P, n_steps, alpha)
%TANGENTIAL_NEIGHBOR_SMOOTH Smooth vertex positions mostly tangentially.
%
%   P = tangential_neighbor_smooth(M, P, n_steps, alpha) repeatedly moves
%   vertices toward their one-ring neighbor average, with the normal
%   component removed at each step. This mainly redistributes vertices on
%   the surface rather than intentionally changing the shape.

    n_v = size(P, 1);
    rows = [M(:, 1); M(:, 2); M(:, 3); M(:, 2); M(:, 3); M(:, 1)];
    cols = [M(:, 2); M(:, 3); M(:, 1); M(:, 1); M(:, 2); M(:, 3)];
    A = sparse(rows, cols, 1, n_v, n_v);
    A = double(A > 0);
    deg = sum(A, 2);

    for step = 1:n_steps
        geo = Geometry(M, P);
        P_avg = (A * P) ./ max(deg, 1);
        dP = P_avg - P;
        dP = dP - dot(dP, geo.v_normal, 2) .* geo.v_normal;
        P = P + alpha * dP;
    end
end
