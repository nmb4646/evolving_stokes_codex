% Nondimensional coupled Newton-Krylov force-surface solver.
%
% Expected caller variables:
%   p       parameter struct
%   dir     output directory ending in '/'
%   verbose logical/scalar
%
% This is the coupled multifield update:
%   P, f, and lambda are solved simultaneously with a matrix-free
%   Newton-Krylov iteration on the force-balance, BIE, and area residuals.



if ~exist('verbose', 'var')
    verbose = false;
end
if ~exist('supress_outputs', 'var')
    supress_outputs = false;
end
if ~isfield(p, 'start')
    p.start = 0;
end
if ~isfield(p, 'subdivisions')
    p.subdivisions = 4;
end
if ~isfield(p, 'roughness')
    p.roughness = 0;
end
if ~isfield(p, 'dt')
    p.dt = 0.05;
end
if ~isfield(p, 'T')
    p.T = 100;
end
if ~isfield(p, 'k')
    p.k = 10;
end
if ~isfield(p, 'remesh_size')
    p.remesh_size = 0;
end
if ~isfield(p, 'initial_remesh')
    p.initial_remesh = false;
end
if ~isfield(p, 'Sd')
    p.Sd = 1;
end
if ~isfield(p, 'Da')
    p.Da = 1;
end
if ~isfield(p, 'Gamma')
    p.Gamma = 0;
end
if ~isfield(p, 'gamy')
    p.gamy = 0;
end
if ~isfield(p, 'chi')
    p.chi = 0.1;
end
if ~isfield(p, 'nk_tol')
    p.nk_tol = 1e-2;
end
if ~isfield(p, 'nk_gmres_tol')
    p.nk_gmres_tol = 0.2;
end
if ~isfield(p, 'nk_gmres_restart')
    p.nk_gmres_restart = 100;
end
if ~isfield(p, 'nk_gmres_maxit')
    p.nk_gmres_maxit = 50;
end
if ~isfield(p, 'nk_fd_relstep')
    p.nk_fd_relstep = 1e-6;
end
if ~isfield(p, 'nk_max_step')
    p.nk_max_step = 0.25;
end
if ~isfield(p, 'nk_area_weight')
    p.nk_area_weight = 1;
end
if ~isfield(p, 'nk_precondition')
    p.nk_precondition = false;
end
if ~isfield(p, 'nk_frozen_jvp')
    p.nk_frozen_jvp = true;
end
if ~isfield(p, 'nk_frozen_switch')
    p.nk_frozen_switch = 0;
end
if ~isfield(p, 'nk_split_force_scale_ratio')
    p.nk_split_force_scale_ratio = 10;
end
if ~isfield(p, 'nk_lm')
    p.nk_lm = 0;
end
if ~isfield(p, 'nk_schur_fd')
    p.nk_schur_fd = true;
end
if ~isfield(p, 'nk_schur_fd_switch')
    p.nk_schur_fd_switch = 0.2;
end

start = p.start;
[status, msg, msgID] = mkdir(dir);

addpath('./isoremesh')

hasRemesher = (exist('remeshing', 'file') ~= 0);
if ~hasRemesher
    warning('Remesher is not installed and has been disabled.');
end

if start == 0
    p.total_time = 0;

    if isfield(p, 'initial_geometry') && strlength(string(p.initial_geometry)) > 0
        loaded = load(p.initial_geometry);
        P = loaded.P;
        M = loaded.M;
    else
        [P, M] = subdivided_sphere(p.subdivisions);
        %P(:,3) = P(:,3)*2;
    end

    geo = Geometry(M, P);
    if p.initial_remesh && hasRemesher
        [M, P] = remeshing(int32(M), P, int32([]), mean(geo.he_length), int32(100));
        M = cast(M, "double");
        geo = Geometry(M, P);
    end
    if p.roughness ~= 0
        geo = perturb_rand(geo, p.roughness);
        P = geo.V;
        M = geo.F;
    end

    p.area0 = geo.area;

    o.h = 1;
    o.eta = 100;
    o.tol_b = .05;
    o.tol_c = .05;
    o.tol_d = .01;
    o.tol_f = o.tol_b;
    o.max_iter = 10000;
    if isfield(p, 'tol_b')
        o.tol_b = p.tol_b;
    end
    if isfield(p, 'tol_c')
        o.tol_c = p.tol_c;
    end
    if isfield(p, 'tol_d')
        o.tol_d = p.tol_d;
    end
    if isfield(p, 'max_iter')
        o.max_iter = p.max_iter;
    end
    if isfield(p, 'h')
        o.h = p.h;
    end
    if isfield(p, 'eta')
        o.eta = p.eta;
    end

    lsr.c = 1e-4;
    lsr.tau = 0.5;
    lsr.max_iter = 20;
    lsr.f_tau = 0.5;
    lsr.f_max_iter = 12;

    r.edge_length = mean(geo.he_length);
    r.n_iter = 50;

    velocity = zeros(size(P, 1) * 3, 1);
    lambda = 0;
    fb = geo.bending_force(1);
    [~, ~, ~, ~, KTK0, DTD0] = geo.evolving_operators();
    twoHn0 = geo.lap * P;
    fv0 = reshape(-2 * (KTK0 + p.k * DTD0) * velocity, size(P));
    fs0 = reshape(lambda * twoHn0, size(P));
    f = -nodal_to_traction(fv0 + fb + fs0, geo);
    slp_cache = stokeslet_SLP_triangle_setup(M);

    save(dir + "geo0.mat", "M", "P", "velocity", "lambda", "f", "fb", "p", "o", "r", "lsr");
    if ~supress_outputs
        fprintf("Save geo0.mat \n");
    end
else
    p_input = p;
    dt_override = p.dt;
    remesh_size = p.remesh_size;
    load(sprintf("%sgeo%d.mat", dir, start), "M", "P", "velocity", "lambda", "p", "o", "r");
    if exist(sprintf("%sgeo%d.mat", dir, start), "file")
        loaded = load(sprintf("%sgeo%d.mat", dir, start), "f");
        if isfield(loaded, "f")
            f = loaded.f;
        else
            f = zeros(size(P));
        end
    end
    p.dt = dt_override;
    p.remesh_size = remesh_size;
    override_fields = ["Sd", "Da", "Gamma", "gamy", "chi", ...
        "nk_tol", "nk_gmres_tol", "nk_gmres_restart", "nk_gmres_maxit", ...
        "nk_fd_relstep", "nk_max_step", "nk_area_weight", ...
        "nk_precondition", "nk_frozen_jvp", "nk_frozen_switch", ...
        "nk_split_force_scale_ratio", "nk_lm", ...
        "nk_schur_fd", "nk_schur_fd_switch", ...
        "tol_b", "tol_c", "tol_d", "max_iter", "h", "eta"];
    for override_idx = 1:numel(override_fields)
        field = override_fields(override_idx);
        if isfield(p_input, field)
            p.(field) = p_input.(field);
        end
    end
    geo = Geometry(M, P);
    if ~isfield(p, 'total_time')
        p.total_time = 0;
    end

    lsr.c = 1e-4;
    lsr.tau = 0.5;
    lsr.max_iter = 20;
    lsr.f_tau = 0.5;
    lsr.f_max_iter = 12;
    
    p.f_gmres_tol = 1e-4;
    p.f_gmres_restart = 50;
    p.f_gmres_maxit = 50;
    p.f_relax = 0.5;

    if ~isfield(o, 'tol_b')
        o.tol_b = o.tol_f;
    end
    if ~isfield(o, 'tol_c')
        o.tol_c = 1e-4;
    end
    if isfield(p, 'tol_b')
        o.tol_b = p.tol_b;
    end
    if isfield(p, 'tol_c')
        o.tol_c = p.tol_c;
    end
    if isfield(p, 'tol_d')
        o.tol_d = p.tol_d;
    end
    if isfield(p, 'max_iter')
        o.max_iter = p.max_iter;
    end
    if isfield(p, 'h')
        o.h = p.h;
    end
    if isfield(p, 'eta')
        o.eta = p.eta;
    end
    slp_cache = stokeslet_SLP_triangle_setup(M);

    if p.remesh_size ~= 0 && hasRemesher
        r.edge_length = p.remesh_size * mean(geo.he_length);
        geo_pre = geo;
        [M, P] = remeshing(int32(M), P, int32([]), r.edge_length, int32(20));
        M = cast(M, "double");
        geo = Geometry(M, P);
        P = P * sqrt(p.area0 / geo.area);
        geo = Geometry(M, P);
        [velocity, f] = map_data(geo, geo_pre, velocity, f);
        slp_cache = stokeslet_SLP_triangle_setup(M);
    end

    %%% OVERRRIDES
    o.tol_b = .05;
    o.tol_c = .05;
    o.tol_d = .01;

    %.h = o.h*3;
    %o.eta = o.eta*.5;
    




end

for t = (start + 1):p.T
    tic;

    [~, ~, ~, ~, KTK, DTD] = geo.evolving_operators();
    P0 = P(:);

    P(:) = P0 + p.dt * velocity;
    geo = Geometry(M, P);

    z = pack_state(P, f, lambda);
    [~, initial_parts] = gs_residual(z, P0, M, KTK, DTD, slp_cache, p, []);
    scales = gs_residual_scales(initial_parts, p);

    j = 0;
    eps_b = Inf;
    eps_c = Inf;
    eps_d = Inf;
    phi = Inf;
    solve_failed = false;
    while ((eps_b > o.tol_b) || (eps_c > o.tol_c) || (eps_d > o.tol_d)) && (j < o.max_iter)
        [R, parts] = gs_residual(z, P0, M, KTK, DTD, slp_cache, p, scales);
        phi = norm(R) / sqrt(numel(R));
        eps_b = residual_rms(parts.r2) / scales.f;
        eps_c = residual_rms(parts.r1) / scales.u;
        eps_d = abs(parts.r3) / scales.area;

        if eps_b <= o.tol_b && eps_c <= o.tol_c && eps_d <= o.tol_d
            break
        end

        if eps_b <= o.tol_b && eps_c > o.tol_c ...
                && scales.u / max(scales.f, 1e-14) <= p.nk_split_force_scale_ratio
            [z_block, block_info] = gs_force_block_step(z, parts, P0, M, KTK, DTD, slp_cache, p, scales, o);
            if block_info.accepted
                z = z_block;
                if verbose
                    fprintf("t = %d, j = %d, f-block GMRES flag = %d, relres = %0.4g, iter = [%d %d], alpha = %0.4g, eps_b = %0.4g, eps_c = %0.4g\n", ...
                        t, j, block_info.flag, block_info.relres, block_info.iter(1), block_info.iter(2), block_info.alpha, block_info.eps_b, block_info.eps_c);
                end
                j = j + 1;
                continue
            end
        end

        schur_scale = gs_schur_state_scale(parts, scales);
        if verbose && j == 0
            fprintf("    NK scales: residual u = %.3e, f = %.3e, area = %.3e; Schur state P = %.3e, lambda = %.3e\n", ...
                scales.u, scales.f, scales.area, schur_scale(1), schur_scale(end));
        end

        p_schur = p;
        p_schur.nk_schur_fd_active = p.nk_schur_fd && (phi < p.nk_schur_fd_switch);
        Schur = @(y) gs_schur_apply_scaled(y, schur_scale, parts, P0, M, KTK, DTD, slp_cache, p_schur, scales);
        rhs = gs_schur_rhs(parts, M, slp_cache, p, scales);
        gmres_tol = min(0.5, max(p.nk_gmres_tol, 0.9 * phi));
        [dy, gmres_flag, gmres_relres, gmres_iter] = gmres(Schur, rhs, ...
            p.nk_gmres_restart, gmres_tol, p.nk_gmres_maxit);
        [dP, dlambda] = unpack_schur_step(schur_scale .* dy, size(parts.P, 1));
        df = gs_recover_force_step(dP, dlambda, parts, KTK, DTD, p);
        dz = pack_state(dP, df, dlambda);

        if gmres_flag ~= 0 || any(~isfinite(dz))
            warning("Schur Newton-Krylov GMRES did not converge at t = %d, j = %d with flag = %d, relres = %0.4g. Timestep not saved.", ...
                t, j, gmres_flag, gmres_relres);
            solve_failed = true;
            break;
        end

        dz_norm = norm(dz);
        if dz_norm > p.nk_max_step * max(1, norm(z))
            dz = dz * (p.nk_max_step * max(1, norm(z)) / dz_norm);
        end

        accepted = false;
        alpha = 1;
        for ls_iter = 1:lsr.max_iter
            z_trial = z + alpha * dz;
            [R_trial, parts_trial] = gs_residual(z_trial, P0, M, KTK, DTD, slp_cache, p, scales);
            phi_trial = norm(R_trial) / sqrt(numel(R_trial));
            if phi_trial <= (1 - lsr.c * alpha) * phi
                z = z_trial;
                parts = parts_trial;
                phi = phi_trial;
                accepted = true;
                break
            end
            alpha = lsr.tau * alpha;
        end

        if ~accepted
            warning("Newton-Krylov line search failed at t = %d, j = %d after converged GMRES relres = %0.4g. Timestep not saved.", ...
                t, j, gmres_relres);
            solve_failed = true;
            break;
        end

        eps_b = residual_rms(parts.r2) / scales.f;
        eps_c = residual_rms(parts.r1) / scales.u;
        eps_d = abs(parts.r3) / scales.area;

        if verbose
            fprintf("t = %d, j = %d, Schur NK phi = %0.4g, eps_b = %0.4g, eps_c = %0.4g, eps_d = %0.4g, GMRES flag = %d, relres = %0.4g, tol = %0.4g, iter = [%d %d], alpha = %0.4g\n", ...
                t, j, phi, eps_b, eps_c, eps_d, gmres_flag, gmres_relres, gmres_tol, gmres_iter(1), gmres_iter(2), alpha);
        end

        j = j + 1;
    end

    if solve_failed
        break;
    end

    [P, f, lambda] = unpack_state(z, size(P, 1));
    geo = Geometry(M, P);

    if j >= o.max_iter
        warning("Terminating at t = %d because j reached o.max_iter = %d. NK phi = %0.4g, eps_b = %0.4g, eps_c = %0.4g, eps_d = %0.4g", ...
            t, o.max_iter, phi, eps_b, eps_c, eps_d);
        break;
    end

    p.total_time = p.total_time + toc;
    [P, velocity] = rm_rigid(P, (P(:) - P0) / p.dt, geo.v_area);
    geo = Geometry(M, P);

    if hasRemesher && 0%deformation_criterion(geo)
        geo_pre = geo;
        if ~supress_outputs
            fprintf("Remeshing. t = %d \n", t);
        end
        [M, P] = remeshing(int32(M), P, int32([]), r.edge_length, int32(r.n_iter));
        M = cast(M, "double");
        geo = Geometry(M, P);
        P = P * sqrt(p.area0 / geo.area);
        geo = Geometry(M, P);
        [velocity, f] = map_data(geo, geo_pre, velocity, f);
        slp_cache = stokeslet_SLP_triangle_setup(M);
    end

    geo = Geometry(M, P);
    P = P * sqrt(p.area0 / geo.area);
    geo = Geometry(M, P);
    fb = geo.bending_force(1);

    save(dir + sprintf("geo%d.mat", t), "M", "P", "velocity", "lambda", "f", "fb", "p", "o", "r", "lsr");
    if ~supress_outputs
        fprintf("Save geo%d.mat at j = %d, eps_b = %0.4g, eps_c = %0.4g, eps_d = %0.4g, total time: %0.4f\n", ...
            t, j, eps_b, eps_c, eps_d, p.total_time);
    end
end

function [velocity, f] = map_data(geo, geo_pre, velocity_pre, f_pre)
    kdtree = KDTreeSearcher(geo_pre.f_center);
    [face, uv, ~, fail] = project(geo_pre.V, geo_pre.F, geo.V, kdtree, 6);
    if fail
        error("projection failed.");
    end
    velocity = interpolate(geo_pre.F, face, uv, reshape(velocity_pre, [], 3));
    velocity = velocity(:);
    f = interpolate(geo_pre.F, face, uv, f_pre);
end

function z = pack_state(P, f, lambda)
    z = [P(:); f(:); lambda];
end

function [P, f, lambda] = unpack_state(z, n_v)
    n_p = 3 * n_v;
    P = reshape(z(1:n_p), [], 3);
    f = reshape(z((n_p + 1):(2 * n_p)), [], 3);
    lambda = z(end);
end

function [R, parts] = gs_residual(z, P0, M, KTK, DTD, slp_cache, p, scales)
    n_v = size(P0, 1) / 3;
    [P, f, lambda] = unpack_state(z, n_v);
    geo = Geometry(M, P);

    u = reshape((P(:) - P0) / p.dt, size(P));
    u_background = shear_flow(P, p.gamy);
    r1 = bie_residual(P, M, f, geo, u, u_background, slp_cache, p);

    fb = geo.bending_force(1);
    twoHn = geo.lap * P;
    fv = reshape(-2 * (KTK + p.k * DTD) * u(:), size(P));
    fs = reshape(lambda * twoHn, size(P));
    r2 = f + nodal_to_traction(fv + fb + fs, geo);

    r3 = geo.area - p.area0;

    parts.P = P;
    parts.f = f;
    parts.lambda = lambda;
    parts.geo = geo;
    parts.u = u;
    parts.r1 = r1;
    parts.r2 = r2;
    parts.r3 = r3;

    if isempty(scales)
        R = [];
        return
    end

    R = [
        r1(:) / scales.u;
        r2(:) / scales.f;
        p.nk_area_weight * r3 / scales.area
    ];
end

function scales = gs_residual_scales(parts, p)
    scales.u = max([residual_rms(parts.r1), residual_rms(parts.u), abs(p.Gamma), 1e-8]);
    scales.f = max([residual_rms(parts.r2), residual_rms(parts.f), 1e-8]);
    scales.area = max(abs(p.area0), 1e-8);
end

function x_scale = gs_schur_state_scale(parts, scales)
    P_scale = max([residual_rms(parts.P), residual_rms(parts.geo.he_length), 1e-8]);
    lambda_shape = gs_G_apply(1, parts);
    lambda_scale = scales.f / max(residual_rms(lambda_shape), 1e-8);
    lambda_scale = max([lambda_scale, abs(parts.lambda), 1e-8]);
    x_scale = [
        P_scale * ones(numel(parts.P), 1);
        lambda_scale
    ];
end

function rhs = gs_schur_rhs(parts, M, slp_cache, p, scales)
    Br2 = gs_B_apply(parts.r2, parts.P, M, parts.geo, slp_cache, p);
    rhs = [
        (-parts.r1(:) + Br2(:)) / scales.u;
        -p.nk_area_weight * parts.r3 / scales.area
    ];
end

function y = gs_schur_apply_scaled(y_scaled, x_scale, parts, P0, M, KTK, DTD, slp_cache, p, scales)
    [dP, dlambda] = unpack_schur_step(x_scale .* y_scaled, size(parts.P, 1));
    if isfield(p, 'nk_schur_fd_active') && p.nk_schur_fd_active
        [A_dP, C_dP, h_dP] = gs_AC_area_apply_fd(dP, parts, P0, M, KTK, DTD, slp_cache, p);
    else
        A_dP = gs_A_apply(dP, parts, p);
        C_dP = gs_C_apply(dP, parts, KTK, DTD, p);
        h_dP = gs_area_apply(dP, parts);
    end
    G_dl = gs_G_apply(dlambda, parts);
    B_CG = gs_B_apply(C_dP + G_dl, parts.P, M, parts.geo, slp_cache, p);
    y = [
        (A_dP(:) - B_CG(:)) / scales.u;
        p.nk_area_weight * h_dP / scales.area
    ];
end

function [A_dP, C_dP, h_dP] = gs_AC_area_apply_fd(dP, parts, P0, M, KTK, DTD, slp_cache, p)
    dP_norm = norm(dP(:));
    if dP_norm == 0
        A_dP = zeros(size(parts.P));
        C_dP = zeros(size(parts.P));
        h_dP = 0;
        return
    end

    fd_step = p.nk_fd_relstep * max(1, norm(parts.P(:))) / dP_norm;
    z_plus = pack_state(parts.P + fd_step * dP, parts.f, parts.lambda);
    z_minus = pack_state(parts.P - fd_step * dP, parts.f, parts.lambda);
    [~, parts_plus] = gs_residual(z_plus, P0, M, KTK, DTD, slp_cache, p, []);
    [~, parts_minus] = gs_residual(z_minus, P0, M, KTK, DTD, slp_cache, p, []);

    A_dP = (parts_plus.r1 - parts_minus.r1) / (2 * fd_step);
    C_dP = (parts_plus.r2 - parts_minus.r2) / (2 * fd_step);
    h_dP = (parts_plus.r3 - parts_minus.r3) / (2 * fd_step);
end

function [dP, dlambda] = unpack_schur_step(x, n_v)
    n_p = 3 * n_v;
    dP = reshape(x(1:n_p), [], 3);
    dlambda = x(end);
end

function df = gs_recover_force_step(dP, dlambda, parts, KTK, DTD, p)
    C_dP = gs_C_apply(dP, parts, KTK, DTD, p);
    G_dl = gs_G_apply(dlambda, parts);
    df = -parts.r2 - C_dP - G_dl;
end

function out = gs_A_apply(dP, parts, p)
    du_background = zeros(size(dP));
    du_background(:, 1) = p.gamy * dP(:, 3);
    out = dP / p.dt - du_background;
end

function out = gs_B_apply(df, P, M, geo, slp_cache, p)
    slp_df = stokeslet_SLP_triangle(P, M, df, slp_cache);
    normal_df = dot(df, geo.v_normal, 2);
    out = p.Sd * slp_df - p.Sd * p.Da * normal_df .* geo.v_normal;
end

function out = gs_C_apply(dP, parts, KTK, DTD, p)
    du = dP / p.dt;
    visc_step = reshape(-2 * (KTK + p.k * DTD) * du(:), [], 3);
    out = nodal_to_traction(visc_step, parts.geo);
end

function out = gs_G_apply(dlambda, parts)
    out = nodal_to_traction(dlambda * (parts.geo.lap * parts.P), parts.geo);
end

function out = gs_area_apply(dP, parts)
    twoHn = parts.geo.lap * parts.P;
    out = -sum(dot(twoHn, dP, 2) .* parts.geo.v_area);
end

function z_scale = gs_state_scale(parts, scales, p)
    P_scale = max([residual_rms(parts.P), residual_rms(parts.geo.he_length), 1e-8]);
    f_scale = max(scales.f, 1e-8);

    lambda_shape = nodal_to_traction(parts.geo.lap * parts.P, parts.geo);
    lambda_scale = f_scale / max(residual_rms(lambda_shape), 1e-8);
    lambda_scale = max([lambda_scale, abs(parts.lambda), 1e-8]);

    z_scale = [
        P_scale * ones(numel(parts.P), 1);
        f_scale * ones(numel(parts.f), 1);
        lambda_scale
    ];
end

function dz = gs_scale_step(y, parts, scales, p)
    n_v = size(parts.P, 1);
    n_p = 3 * n_v;

    P_scale = max([residual_rms(parts.P), residual_rms(parts.geo.he_length), 1e-8]);
    f_tangent_scale = max(scales.f, 1e-8);
    if ~isfield(p, 'nk_use_split_force_scale') || ~p.nk_use_split_force_scale ...
            || scales.u / f_tangent_scale > p.nk_split_force_scale_ratio
        f_normal_scale = f_tangent_scale;
    else
        f_normal_scale = min(f_tangent_scale, max(scales.u / max(p.Sd * p.Da, 1e-14), 1e-8));
    end

    yP = reshape(y(1:n_p), [], 3);
    yf = reshape(y((n_p + 1):(2 * n_p)), [], 3);

    n = parts.geo.v_normal;
    yf_normal = dot(yf, n, 2) .* n;
    yf_tangent = yf - yf_normal;

    lambda_shape = nodal_to_traction(parts.geo.lap * parts.P, parts.geo);
    lambda_scale = f_tangent_scale / max(residual_rms(lambda_shape), 1e-8);
    lambda_scale = max([lambda_scale, abs(parts.lambda), 1e-8]);

    df = f_tangent_scale * yf_tangent + f_normal_scale * yf_normal;
    dz = [
        P_scale * yP(:);
        df(:);
        lambda_scale * y(end)
    ];
end

function dy0 = gs_initial_krylov_guess(parts, scales, z_scale, p)
    n_p = numel(parts.P);
    dy0 = zeros(2 * n_p + 1, 1);
    if scales.u / max(scales.f, 1e-14) > p.nk_split_force_scale_ratio
        return
    end

    normal_residual = dot(parts.r1, parts.geo.v_normal, 2);
    df0 = (normal_residual / max(p.Sd * p.Da, 1e-14)) .* parts.geo.v_normal;
    dz0 = [
        zeros(n_p, 1);
        df0(:);
        0
    ];
    dy0 = dz0 ./ z_scale;
    dy0(~isfinite(dy0)) = 0;
end

function [z_new, info] = gs_force_block_step(z, parts, P0, M, KTK, DTD, slp_cache, p, scales, o)
    info.accepted = false;
    info.flag = NaN;
    info.relres = Inf;
    info.iter = [0, 0];
    info.alpha = 0;
    info.eps_b = Inf;
    info.eps_c = Inf;
    z_new = z;

    rhs = -parts.r1(:) / scales.u;
    Af = @(df_vec) gs_force_block_apply(df_vec, parts.P, M, parts.geo, slp_cache, p, scales);
    [df_vec, flag, relres, iter] = gmres(Af, rhs, p.nk_gmres_restart, ...
        max(p.nk_gmres_tol, 0.2), p.nk_gmres_maxit);

    info.flag = flag;
    info.relres = relres;
    info.iter = iter;
    if flag ~= 0 || any(~isfinite(df_vec))
        return
    end

    n_p = numel(parts.P);
    phi = norm(gs_residual(z, P0, M, KTK, DTD, slp_cache, p, scales)) / sqrt(2 * n_p + 1);
    for ls_iter = 1:20
        z_trial = z;
        z_trial((n_p + 1):(2 * n_p)) = z_trial((n_p + 1):(2 * n_p)) + df_vec;
        z_trial((n_p + 1):(2 * n_p)) = z((n_p + 1):(2 * n_p)) + (0.5 ^ (ls_iter - 1)) * df_vec;
        R_trial = gs_residual(z_trial, P0, M, KTK, DTD, slp_cache, p, scales);
        [~, parts_trial] = gs_residual(z_trial, P0, M, KTK, DTD, slp_cache, p, []);
        eps_b_trial = residual_rms(parts_trial.r2) / scales.f;
        eps_c_trial = residual_rms(parts_trial.r1) / scales.u;
        phi_trial = norm(R_trial) / sqrt(numel(R_trial));
        eps_b_limit = max(o.tol_b, 1.2 * residual_rms(parts.r2) / scales.f);
        if phi_trial < phi && eps_b_trial <= eps_b_limit
            z_new = z_trial;
            info.accepted = true;
            info.alpha = 0.5 ^ (ls_iter - 1);
            info.eps_b = eps_b_trial;
            info.eps_c = eps_c_trial;
            return
        end
    end
end

function out = gs_force_block_apply(df_vec, P, M, geo, slp_cache, p, scales)
    df = reshape(df_vec, [], 3);
    slp_df = stokeslet_SLP_triangle(P, M, df, slp_cache);
    normal_df = dot(df, geo.v_normal, 2);
    out = (p.Sd * slp_df - p.Sd * p.Da * normal_df .* geo.v_normal) / scales.u;
    out = out(:);
end

function Jv = gs_scaled_jacobian_vector(v, z_scale, z, R, P0, M, KTK, DTD, slp_cache, p, scales)
    dz = z_scale .* v;
    Jv = gs_jacobian_vector(dz, z, R, P0, M, KTK, DTD, slp_cache, p, scales);
end

function dz = gs_precondition_step(y, parts, scales, KTK, DTD, p)
    n_v = size(parts.P, 1);
    n_p = 3 * n_v;

    y1 = reshape(y(1:n_p), [], 3) * scales.u;
    y2 = reshape(y((n_p + 1):(2 * n_p)), [], 3) * scales.f;
    y3 = y(end) * scales.area / max(p.nk_area_weight, 1e-14);

    n = parts.geo.v_normal;
    y1n_scalar = dot(y1, n, 2);
    y1n = y1n_scalar .* n;
    y1t = y1 - y1n;

    y2n_scalar = dot(y2, n, 2);
    y2t = y2 - y2n_scalar .* n;

    slip_coef = max(p.Sd * p.Da, 1e-14);
    df = y2t - (y1n_scalar / slip_coef) .* n;
    dP = p.dt * y1t;

    area_grad = parts.geo.lap * parts.P;
    area_dir = -area_grad;
    area_slope = sum(dot(area_grad, area_dir, 2) .* parts.geo.v_area);
    if abs(area_slope) > 1e-14
        dP = dP + (y3 / area_slope) * area_dir;
    end

    u_step = dP(:) / p.dt;
    visc_step = reshape(-2 * (KTK + p.k * DTD) * u_step, [], 3);
    df = df - nodal_to_traction(visc_step, parts.geo);

    lambda_shape = nodal_to_traction(parts.geo.lap * parts.P, parts.geo);
    lambda_norm2 = sum(lambda_shape(:) .^ 2);
    if lambda_norm2 > 1e-14
        r2_after = df + nodal_to_traction(visc_step, parts.geo);
        dlambda = -sum(r2_after(:) .* lambda_shape(:)) / lambda_norm2;
    else
        dlambda = 0;
    end

    dz = pack_state(dP, df, dlambda);
end

function Jv = gs_jacobian_vector(v, z, R, P0, M, KTK, DTD, slp_cache, p, scales)
    n_v = size(P0, 1) / 3;
    [P, ~, ~] = unpack_state(z, n_v);
    [dP, df, dlambda] = unpack_state(v, n_v);

    if p.nk_frozen_jvp
        Jv = gs_frozen_jacobian_vector(P, dP, df, dlambda, M, KTK, DTD, slp_cache, p, scales);
        return
    end

    JR = zeros(size(R));

    dP_norm = norm(dP(:));
    if dP_norm > 0
        fd_step = p.nk_fd_relstep * max(1, norm(P(:))) / dP_norm;
        z_plus = z;
        z_minus = z;
        z_plus(1:numel(P)) = z_plus(1:numel(P)) + fd_step * dP(:);
        z_minus(1:numel(P)) = z_minus(1:numel(P)) - fd_step * dP(:);
        R_plus = gs_residual(z_plus, P0, M, KTK, DTD, slp_cache, p, scales);
        R_minus = gs_residual(z_minus, P0, M, KTK, DTD, slp_cache, p, scales);
        JR = JR + (R_plus - R_minus) / (2 * fd_step);
    end

    geo = Geometry(M, P);
    df_norm = norm(df(:));
    if df_norm > 0
        slp_df = stokeslet_SLP_triangle(P, M, df, slp_cache);
        normal_df = dot(df, geo.v_normal, 2);
        Jr1_df = p.Sd * slp_df - p.Sd * p.Da * normal_df .* geo.v_normal;
        Jr2_df = df;
        JR = JR + [
            Jr1_df(:) / scales.u;
            Jr2_df(:) / scales.f;
            0
        ];
    end

    if dlambda ~= 0
        Jr2_lambda = nodal_to_traction(dlambda * (geo.lap * P), geo);
        JR = JR + [
            zeros(numel(P), 1);
            Jr2_lambda(:) / scales.f;
            0
        ];
    end

    Jv = JR;
end

function Jv = gs_frozen_jacobian_vector(P, dP, df, dlambda, M, KTK, DTD, slp_cache, p, scales)
    geo = Geometry(M, P);
    n_p = numel(P);

    du = dP / p.dt;
    du_background = zeros(size(P));
    du_background(:, 1) = p.gamy * dP(:, 3);

    slp_df = stokeslet_SLP_triangle(P, M, df, slp_cache);
    normal_df = dot(df, geo.v_normal, 2);
    Jr1 = du - du_background ...
        + p.Sd * slp_df ...
        - p.Sd * p.Da * normal_df .* geo.v_normal;

    visc_step = reshape(-2 * (KTK + p.k * DTD) * du(:), [], 3);
    tension_step = dlambda * (geo.lap * P);
    Jr2 = df + nodal_to_traction(visc_step + tension_step, geo);

    twoHn = geo.lap * P;
    Jr3 = -sum(dot(twoHn, dP, 2) .* geo.v_area);

    Jv = [
        Jr1(:) / scales.u;
        Jr2(:) / scales.f;
        p.nk_area_weight * Jr3 / scales.area
    ];

    if numel(Jv) ~= 2 * n_p + 1
        error("Frozen JVP returned unexpected size.");
    end
end

function out = residual_rms(x)
    out = norm(x(:)) / sqrt(numel(x));
end

function c = bie_residual(P, M, f, geo, u, u_background, slp_cache, p)
    slpout = stokeslet_SLP_triangle(P, M, f, slp_cache);
    normal_slip = p.Gamma + dot(f, geo.v_normal, 2);
    c = u - u_background + p.Sd * slpout ...
        - p.Sd * p.Da * normal_slip .* geo.v_normal;
end

function eps_b = force_balance_residual(P, P0, M, f, lambda, KTK, DTD, p)
    geo = Geometry(M, P);
    f_nodal = traction_to_nodal(f, geo);
    fb = geo.bending_force(1);
    twoHn = geo.lap * P;
    u = reshape((P(:) - P0) / p.dt, size(P));
    fv = reshape(-2 * (KTK + p.k * DTD) * u(:), size(P));
    fs = reshape(lambda * twoHn, size(P));
    f_mem = fv + fb + fs;
    res = f_nodal + f_mem;
    res_rms = norm(res(:)) / sqrt(numel(res));
    u_rms = norm(u(:)) / sqrt(numel(u));
    eps_b = res_rms / max(u_rms, 1e-14);
end

function nodal = traction_to_nodal(traction, geo)
    nodal = traction .* geo.v_area;
end

function traction = nodal_to_traction(nodal, geo)
    traction = nodal ./ geo.v_area;
end
