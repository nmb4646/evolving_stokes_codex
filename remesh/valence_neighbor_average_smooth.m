function [P, info] = valence_neighbor_average_smooth(M, P, target_valences, n_steps, alpha)
%VALENCE_NEIGHBOR_AVERAGE_SMOOTH Average selected-valence vertices.
%
%   P = valence_neighbor_average_smooth(M, P, [5 7], n_steps, alpha)
%   replaces vertices with target valence by their one-ring neighbor average.
%   The update is simultaneous within each pass. alpha = 1 is exact
%   replacement; smaller alpha damps the move.

    if nargin < 3 || isempty(target_valences)
        target_valences = [5, 7];
    end
    if nargin < 4 || isempty(n_steps)
        n_steps = 1;
    end
    if nargin < 5 || isempty(alpha)
        alpha = 1;
    end

    n_v = size(P, 1);
    rows = [M(:, 1); M(:, 2); M(:, 3); M(:, 2); M(:, 3); M(:, 1)];
    cols = [M(:, 2); M(:, 3); M(:, 1); M(:, 1); M(:, 2); M(:, 3)];
    A = sparse(rows, cols, 1, n_v, n_v);
    A = double(A > 0);
    valence = full(sum(A, 2));
    target = ismember(valence, target_valences(:));

    info.valence = valence;
    info.target = target;
    info.n_target = nnz(target);
    info.max_displacement = zeros(n_steps, 1);
    info.mean_displacement = zeros(n_steps, 1);

    deg = max(valence, 1);
    for step = 1:n_steps
        P_avg = (A * P) ./ deg;
        dP = zeros(size(P));
        dP(target, :) = P_avg(target, :) - P(target, :);
        P = P + alpha * dP;

        step_disp = vecnorm(alpha * dP(target, :), 2, 2);
        if isempty(step_disp)
            info.max_displacement(step) = 0;
            info.mean_displacement(step) = 0;
        else
            info.max_displacement(step) = max(step_disp);
            info.mean_displacement(step) = mean(step_disp);
        end
    end
end
