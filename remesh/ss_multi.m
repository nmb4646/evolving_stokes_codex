% Nondimensional scalar-surface relaxation solver.
%
% Expected caller variables:
%   p       parameter struct
%   dir     output directory ending in '/'
%   verbose logical/scalar
%
% This is the barebones Willmore/nematic NCGS flow with all nematic,
% activity, hydrodynamic BIE, and local pressure fields removed.

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
        load("genus6_smooth.mat");
    end

    geo = Geometry(M, P);
    if p.initial_remesh && hasRemesher
        [M, P] = remeshing(int32(M), P, int32([]), mean(geo.he_length), int32(10));
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
    o.tol_f = 1e-4;
    o.tol_d = 5e-4;
    o.max_iter = 10000;

    lsr.c = 1e-4;
    lsr.tau = 0.5;
    lsr.max_iter = 20;

    r.edge_length = mean(geo.he_length);
    r.n_iter = 50;

    velocity = zeros(size(P, 1) * 3, 1);
    lambda = 0;
    fb = geo.bending_force(1);

    save(dir + "geo0.mat", "M", "P", "velocity", "lambda", "fb", "p", "o", "r", "lsr");
    if ~supress_outputs
        fprintf("Save geo0.mat \n");
    end
else
    dt_override = p.dt;
    remesh_size = p.remesh_size;
    load(sprintf("%sgeo%d.mat", dir, start), "M", "P", "velocity", "lambda", "p", "o", "r");
    p.dt = dt_override;
    p.remesh_size = remesh_size;
    geo = Geometry(M, P);
    if ~isfield(p, 'total_time')
        p.total_time = 0;
    end

    lsr.c = 1e-4;
    lsr.tau = 0.5;
    lsr.max_iter = 20;

    if p.remesh_size ~= 0 && hasRemesher
        r.edge_length = p.remesh_size * mean(geo.he_length);
        geo_pre = geo;
        [M, P] = remeshing(int32(M), P, int32([]), r.edge_length, int32(20));
        M = cast(M, "double");
        geo = Geometry(M, P);
        P = P * sqrt(p.area0 / geo.area);
        geo = Geometry(M, P);
        velocity = map_data(geo, geo_pre, velocity);
    end
end

for t = (start + 1):p.T
    tic;

    [~, ~, ~, ~, KTK, DTD] = geo.evolving_operators();
    P0 = P(:);

    mass0 = spdiags(geo.v_area, 0, geo.mesh.n_v, geo.mesh.n_v);
    mass0_inv = spdiags(1 ./ geo.v_area, 0, geo.mesh.n_v, geo.mesh.n_v);
    bih = geo.lap * mass0_inv * geo.lap;
    bih = blkdiag(bih, bih, bih);
    Hess = 2 * (KTK + p.k * DTD) ...
        + 0.5 * p.dt * bih ...
        + 1e-3 * r.edge_length^(-2) * blkdiag(mass0, mass0, mass0);

    lagrangian = @(P_, willmore_, lambda_, area_) ...
        (P_ - P0)' * (KTK + p.k * DTD) * (P_ - P0) ...
        + p.dt * willmore_ ...
        - p.dt * (area_ - p.area0) * lambda_;

    P(:) = P0 + p.dt * velocity;
    geo = Geometry(M, P);

    eps_f = Inf;
    eps_d = Inf;
    j = 0;
    alpha_mem = o.h;

    while ((eps_f > o.tol_f) || (eps_d > o.tol_d)) && (j < o.max_iter)
        E = lagrangian(P(:), geo.willmore_energy(1), lambda, geo.area);

        fb = geo.bending_force(1);
        twoHn = geo.lap * P;
        b = -2 * (KTK + p.k * DTD) * (P(:) - P0) ...
            + p.dt * fb(:) ...
            + p.dt * lambda * twoHn(:);

        dP = Hess \ b;

        quad = dP' * b;
        if quad <= 0
            error('quad = %0.4g should be positive definite', quad);
        end
        eps_f = sqrt(quad);

        accepted = false;
        alpha = min(alpha_mem, 1 / max(1, eps_f));
        for ls_iter = 1:lsr.max_iter
            P_new = reshape(P(:) + alpha * dP, [], 3);
            geo_new = Geometry(M, P_new);
            dE = E - lagrangian(P_new(:), geo_new.willmore_energy(1), lambda, geo_new.area);
            pred = alpha * quad;
            if dE >= lsr.c * pred
                P = P_new;
                geo = geo_new;
                accepted = true;
                alpha_mem = min(alpha / lsr.tau, o.h);
                break;
            end
            alpha = lsr.tau * alpha;
        end
        if ~accepted
            error("backtracking fails");
        end

        darea = geo.area - p.area0;
        eps_d = abs(darea) / p.area0;
        lambda = lambda - o.eta * darea;

        if verbose
            fprintf("t = %d, j = %d, ls_iter = %d, E = %0.10g, eps_f = %0.4g, eps_d = %0.4g, darea = %0.4g \n", ...
                t, j, ls_iter, E, eps_f, eps_d, darea);
        end

        j = j + 1;
    end

    if j >= o.max_iter
        warning("Terminating at t = %d because j reached o.max_iter = %d. eps_f = %0.4g, eps_d = %0.4g", ...
            t, o.max_iter, eps_f, eps_d);
        break;
    end

    p.total_time = p.total_time + toc;
    [P, velocity] = rm_rigid(P, (P(:) - P0) / p.dt, geo.v_area);
    geo = Geometry(M, P);

    if hasRemesher && deformation_criterion(geo)
        geo_pre = geo;
        if ~supress_outputs
            fprintf("Remeshing. t = %d \n", t);
        end
        [M, P] = remeshing(int32(M), P, int32([]), r.edge_length, int32(r.n_iter));
        M = cast(M, "double");
        geo = Geometry(M, P);
        P = P * sqrt(p.area0 / geo.area);
        geo = Geometry(M, P);
        velocity = map_data(geo, geo_pre, velocity);
    end

    geo = Geometry(M, P);
    P = P * sqrt(p.area0 / geo.area);
    geo = Geometry(M, P);
    fb = geo.bending_force(1);

    save(dir + sprintf("geo%d.mat", t), "M", "P", "velocity", "lambda", "fb", "p", "o", "r", "lsr");
    if ~supress_outputs
        fprintf("Save geo%d.mat at j = %d, eps_f = %0.4g, eps_d = %0.4g, total time: %0.4f\n", ...
            t, j, eps_f, eps_d, p.total_time);
    end
end

function velocity = map_data(geo, geo_pre, velocity_pre)
    kdtree = KDTreeSearcher(geo_pre.f_center);
    [face, uv, ~, fail] = project(geo_pre.V, geo_pre.F, geo.V, kdtree, 6);
    if fail
        error("projection failed.");
    end
    velocity = interpolate(geo_pre.F, face, uv, reshape(velocity_pre, [], 3));
    velocity = velocity(:);
end
