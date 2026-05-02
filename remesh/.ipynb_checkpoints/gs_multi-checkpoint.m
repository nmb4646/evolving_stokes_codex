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
    p.nk_tol = 1e-3;
end
if ~isfield(p, 'nk_gmres_tol')
    p.nk_gmres_tol = 1e-2;
end
if ~isfield(p, 'nk_gmres_restart')
    p.nk_gmres_restart = 50;
end
if ~isfield(p, 'nk_gmres_maxit')
    p.nk_gmres_maxit = 20;
end
if ~isfield(p, 'nk_fd_relstep')
    p.nk_fd_relstep = sqrt(eps);
end
if ~isfield(p, 'nk_max_step')
    p.nk_max_step = 1;
end
if ~isfield(p, 'nk_area_weight')
    p.nk_area_weight = 1;
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
    o.tol_b = .01;
    o.tol_c = .01;
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
    while ((eps_b > o.tol_b) || (eps_c > o.tol_c) || (eps_d > o.tol_d) || (phi > p.nk_tol)) && (j < o.max_iter)
        [R, parts] = gs_residual(z, P0, M, KTK, DTD, slp_cache, p, scales);
        phi = norm(R) / sqrt(numel(R));
        eps_b = residual_rms(parts.r2) / scales.f;
        eps_c = residual_rms(parts.r1) / scales.u;
        eps_d = abs(parts.r3) / scales.area;

        if phi <= p.nk_tol && eps_b <= o.tol_b && eps_c <= o.tol_c && eps_d <= o.tol_d
            break
        end

        Jv = @(v) gs_jacobian_vector(v, z, R, P0, M, KTK, DTD, slp_cache, p, scales);
        [dz, gmres_flag, gmres_relres, gmres_iter] = gmres(Jv, -R, ...
            p.nk_gmres_restart, p.nk_gmres_tol, p.nk_gmres_maxit);

        if gmres_flag > 1 || any(~isfinite(dz))
            warning("Newton-Krylov GMRES failed at t = %d, j = %d with flag = %d, relres = %0.4g.", ...
                t, j, gmres_flag, gmres_relres);
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
            warning("Newton-Krylov line search failed at t = %d, j = %d.", t, j);
            break;
        end

        eps_b = residual_rms(parts.r2) / scales.f;
        eps_c = residual_rms(parts.r1) / scales.u;
        eps_d = abs(parts.r3) / scales.area;

        if verbose
            fprintf("t = %d, j = %d, NK phi = %0.4g, eps_b = %0.4g, eps_c = %0.4g, eps_d = %0.4g, GMRES flag = %d, relres = %0.4g, iter = [%d %d], alpha = %0.4g\n", ...
                t, j, phi, eps_b, eps_c, eps_d, gmres_flag, gmres_relres, gmres_iter(1), gmres_iter(2), alpha);
        end

        j = j + 1;
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

function Jv = gs_jacobian_vector(v, z, R, P0, M, KTK, DTD, slp_cache, p, scales)
    v_norm = norm(v);
    if v_norm == 0
        Jv = zeros(size(R));
        return
    end
    fd_step = p.nk_fd_relstep * max(1, norm(z)) / v_norm;
    R_perturbed = gs_residual(z + fd_step * v, P0, M, KTK, DTD, slp_cache, p, scales);
    Jv = (R_perturbed - R) / fd_step;
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
