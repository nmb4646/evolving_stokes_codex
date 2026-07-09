function [P, geo, info] = newton_correct_volume(geo, target_area, target_volume, opts)
%NEWTON_CORRECT_VOLUME Project a remeshed surface back to target area/volume.
%
%   [P, geo, info] = newton_correct_volume(geo, target_area, target_volume)
%   applies small normal corrections to geo.V until the Geometry area and
%   signed volume match the requested targets.
%
%   This is intended for use immediately after remeshing. For impermeable
%   runs, target_volume should usually be the initial volume. For permeable
%   runs, target_volume should usually be the pre-remesh volume from the
%   current timestep, so remeshing itself adds no extra volume change.

    if nargin < 2 || isempty(target_area)
        target_area = geo.area;
    end
    if nargin < 3 || isempty(target_volume)
        target_volume = geo.volume;
    end
    if nargin < 4 || isempty(opts)
        opts = struct();
    end

    opts = with_default(opts, "max_iter", 12);
    opts = with_default(opts, "tol_area", 1e-12);
    opts = with_default(opts, "tol_volume", 1e-12);
    opts = with_default(opts, "damping", 1);
    opts = with_default(opts, "verbose", false);
    opts = with_default(opts, "line_search", true);
    opts = with_default(opts, "line_search_tau", 0.5);
    opts = with_default(opts, "line_search_max_iter", 12);
    opts = with_default(opts, "rcond_tol", 1e-10);

    M = geo.F;
    P = geo.V;

    info.area = zeros(opts.max_iter + 1, 1);
    info.volume = zeros(opts.max_iter + 1, 1);
    info.rel_area_error = zeros(opts.max_iter + 1, 1);
    info.rel_volume_error = zeros(opts.max_iter + 1, 1);
    info.step = zeros(opts.max_iter, 2);
    info.step_norm = zeros(opts.max_iter, 1);
    info.converged = false;
    info.iter = 0;

    for iter = 0:opts.max_iter
        geo = Geometry(M, P);
        area_error = target_area - geo.area;
        volume_error = target_volume - geo.volume;
        rel_area_error = abs(area_error) / max(abs(target_area), eps);
        rel_volume_error = abs(volume_error) / max(abs(target_volume), eps);

        info.area(iter + 1) = geo.area;
        info.volume(iter + 1) = geo.volume;
        info.rel_area_error(iter + 1) = rel_area_error;
        info.rel_volume_error(iter + 1) = rel_volume_error;

        if opts.verbose
            fprintf("newton_correct_volume iter %d: rel area err %.4e, rel volume err %.4e\n", ...
                iter, rel_area_error, rel_volume_error);
        end

        if rel_area_error <= opts.tol_area && rel_volume_error <= opts.tol_volume
            info.converged = true;
            info.iter = iter;
            trim_info();
            return
        end

        if iter == opts.max_iter
            info.iter = iter;
            trim_info();
            return
        end

        area_gradient = geo.lap * P;
        volume_gradient = signed_volume_gradient(M, P);
        modes = [area_gradient(:), volume_gradient(:)];
        J = [area_gradient(:).'; volume_gradient(:).'] * modes;
        rhs = [area_error; volume_error];

        if rcond(full(J)) < opts.rcond_tol
            warning("newton_correct_volume:IllConditionedCorrection", ...
                "Area/volume correction Jacobian is ill-conditioned; using pseudoinverse step.");
            coeff = pinv(full(J)) * rhs;
        else
            coeff = J \ rhs;
        end

        coeff = opts.damping * coeff;
        dP = reshape(modes * coeff, size(P));
        alpha = choose_step(P, M, dP, target_area, target_volume, ...
            rel_area_error, rel_volume_error, opts);
        P = P + alpha * dP;

        info.step(iter + 1, :) = (alpha * coeff).';
        info.step_norm(iter + 1) = norm(alpha * dP(:));
    end

    function trim_info()
        keep_state = 1:(info.iter + 1);
        keep_step = 1:max(info.iter, 1);
        info.area = info.area(keep_state);
        info.volume = info.volume(keep_state);
        info.rel_area_error = info.rel_area_error(keep_state);
        info.rel_volume_error = info.rel_volume_error(keep_state);
        if info.iter == 0
            info.step = zeros(0, 2);
            info.step_norm = zeros(0, 1);
        else
            info.step = info.step(keep_step, :);
            info.step_norm = info.step_norm(keep_step);
        end
    end
end

function opts = with_default(opts, field, value)
    field = char(field);
    if ~isfield(opts, field) || isempty(opts.(field))
        opts.(field) = value;
    end
end

function alpha = choose_step(P, M, dP, target_area, target_volume, ...
        rel_area_error, rel_volume_error, opts)
    alpha = 1;
    if ~opts.line_search
        return
    end

    phi0 = correction_merit(rel_area_error, rel_volume_error);
    best_alpha = 0;
    best_phi = Inf;
    for ls_iter = 1:opts.line_search_max_iter
        geo_trial = Geometry(M, P + alpha * dP);
        area_error_trial = abs(target_area - geo_trial.area) / max(abs(target_area), eps);
        volume_error_trial = abs(target_volume - geo_trial.volume) / max(abs(target_volume), eps);
        phi_trial = correction_merit(area_error_trial, volume_error_trial);
        if phi_trial < best_phi
            best_phi = phi_trial;
            best_alpha = alpha;
        end
        if phi_trial < phi0
            return
        end
        alpha = alpha * opts.line_search_tau;
    end

    alpha = best_alpha;
end

function phi = correction_merit(rel_area_error, rel_volume_error)
    phi = 0.5 * (rel_area_error^2 + rel_volume_error^2);
end

function grad = signed_volume_gradient(M, P)
    grad = zeros(size(P));
    i = M(:, 1);
    j = M(:, 2);
    k = M(:, 3);

    gi = cross(P(j, :), P(k, :), 2) / 6;
    gj = cross(P(k, :), P(i, :), 2) / 6;
    gk = cross(P(i, :), P(j, :), 2) / 6;

    grad = add_rows(grad, i, gi);
    grad = add_rows(grad, j, gj);
    grad = add_rows(grad, k, gk);
end

function A = add_rows(A, rows, values)
    for dim = 1:size(A, 2)
        A(:, dim) = A(:, dim) + accumarray(rows, values(:, dim), [size(A, 1), 1], @sum, 0);
    end
end
