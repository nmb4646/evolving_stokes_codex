% Nondimensional force-surface relaxation solver.
%
% Expected caller variables:
%   p       parameter struct
%   dir     output directory ending in '/'
%   verbose logical/scalar
%
% This is the first multifield force update:
%   P, lambda are advanced by the scalar surface relaxation step with f fixed.
%   f is then updated from the permeable BIE residual.



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
        P(:,3) = P(:,3)*2;
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

    o.h = .01;
    o.eta = 1;
    o.tol_b = .001;
    o.tol_c = .001;
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
    f = -0*nodal_to_traction(fv0 + fb + fs0, geo);
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
    override_fields = ["Sd", "Da", "Gamma", "gamy", "chi", "tol_b", "tol_c", "tol_d", "max_iter", "h", "eta"];
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
        disp("Remeshing on start")
        r.edge_length = p.remesh_size * mean(geo.he_length);
        geo_pre = geo;
        [M, P] = remeshing(int32(M), P, int32([]), r.edge_length, int32(20));
        M = cast(M, "double");
        geo = Geometry(M, P);

        [P, smooth_info] = valence_neighbor_average_smooth(M, P, ...
        [5, 7], 1, .25);

        P = P * sqrt(p.area0 / geo.area);
        geo = Geometry(M, P);
        [velocity, f] = map_data(geo, geo_pre, velocity, f);
        slp_cache = stokeslet_SLP_triangle_setup(M);
        disp("Old: " + length(geo_pre.V))
        disp("New: " + length(geo.V))
    end

    %%% OVERRRIDES
    o.tol_b = .01;
    % o.tol_c = .03;
    % o.tol_d = .01;

    %.h = o.h*3;
    %o.eta = 0;
    




end

for t = (start + 1):p.T
    tic;

    [~, ~, ~, ~, KTK, DTD] = geo.evolving_operators();
    P0 = P(:);
    lambda_restart = lambda;
    f_restart = f;

    mass0 = spdiags(geo.v_area, 0, geo.mesh.n_v, geo.mesh.n_v);
    mass0_inv = spdiags(1 ./ geo.v_area, 0, geo.mesh.n_v, geo.mesh.n_v);
    bih = geo.lap * mass0_inv * geo.lap;
    bih = blkdiag(bih, bih, bih);
    Hess = 2 * (KTK + p.k * DTD) ...
        + 0.5 * p.dt * bih ...
        + 1e-3 * r.edge_length^(-2) * blkdiag(mass0, mass0, mass0);



    lagrangian = @(P_, willmore_, lambda_, area_, f_nodal_) ...
        (P_ - P0)' * (KTK + p.k * DTD) * (P_ - P0) ...
        + p.dt * willmore_ ...
        - p.dt * (area_ - p.area0) * lambda_...
        - p.dt * f_nodal_(:)' * (P_ - P0);





    S = stokeslet_SLP_triangle_matrix(P, M, slp_cache);
    N = normal_projection_matrix(geo.v_normal);
    NN = blkdiag(mass0, mass0, mass0) * N;

    P(:) = P0 + p.dt * velocity;
    geo = Geometry(M, P);

    eps_b = Inf;
    eps_c = Inf;
    eps_d = Inf;
    eps_b_prev = Inf;
    eps_b_rise_count = 0;
    j = 0;
    alpha_mem = 1;


    % I will comment out from here ...
    n_state = numel(P);
    M3 = blkdiag(mass0, mass0, mass0);
    I_state = speye(n_state);
    bie_df = p.Sd * (S + p.Da * N);
    % full_Jac = [
    %     Hess,              -p.dt * M3;
    %     -I_state / p.dt,    bie_df
    % ];
    % fprintf('R_{uu} condition number: %.4e \n', condest(Hess))
    % fprintf('Hessian condition number: %.4e \n', condest(full_Jac))
    % to here



    %Codex: add new hessian from here

    o.eta = 0;



    % to here

    while ((eps_b > o.tol_b) || (eps_c > o.tol_c) || (eps_d > o.tol_d)) && (j < o.max_iter)

        f_nodal = traction_to_nodal(f, geo);
        fb = geo.bending_force(1);
        twoHn = geo.lap * P;
        u = reshape((P(:) - P0) / p.dt, size(P));
        fv = reshape(-2 * (KTK + p.k * DTD) * u(:), size(P));
        fs = reshape(lambda * twoHn, size(P));
        b = p.dt * force_balance_residual(P, P0, M, f, lambda, KTK, DTD, p);
        P_bie = reshape(P0, [], 3);
        geo_bie = Geometry(M, P_bie);
        u_background = shearextensionflow(P_bie, -p.gamy);
        c = bie_residual(P_bie, M, f, geo_bie, u, u_background, slp_cache, p);


        % I will comment out from here ...
        % full_rhs = [
        %     -b(:);
        %     -c(:);
        % ];
        % full_step = full_Jac \ full_rhs;
        % 
        % 
        % alpha = p.chi;
        % dP = full_step(1:n_state);
        % df = reshape(full_step((n_state + 1):end), size(P));
        % P = reshape(P(:) + alpha * dP, [], 3);
        % f = f + alpha * df;
        % geo = Geometry(M, P);
        % ls_iter = 1;
        %... to here.

        % Codex: only add stuff in these bounds...

        d = (geo.area - p.area0) / p.area0;
        area_gradient = geo.lap * P;
        area_gradient = area_gradient(:);
        force_gamma_block = -p.dt * area_gradient;
        area_gamma_block = area_gradient.' / p.area0;
        full_gamma_Jac = [
            Hess,              -p.dt * M3,             force_gamma_block;
            -I_state / p.dt,    bie_df,                sparse(n_state, 1);
            area_gamma_block,   sparse(1, n_state),    0
        ];
        if j==0 
        fprintf('Full gamma Hessian condition number: %.4e \n', condest(full_gamma_Jac))
        end
        full_gamma_rhs = [
            -b(:);
            -c(:);
            -d
        ];
        row_scale = 1 ./ max(sum(abs(full_gamma_Jac), 2), eps);
        col_scale = 1 ./ max(sum(abs(full_gamma_Jac), 1).', eps);
        Dr = spdiags(sqrt(row_scale), 0, size(full_gamma_Jac, 1), size(full_gamma_Jac, 1));
        Dc = spdiags(sqrt(col_scale), 0, size(full_gamma_Jac, 2), size(full_gamma_Jac, 2));
        if j==0
            fprintf('Equilibrated full gamma Hessian condition number: %.4e \n', ...
                condest(Dr * full_gamma_Jac * Dc))
        end
        full_gamma_step = Dc * ((Dr * full_gamma_Jac * Dc) \ (Dr * full_gamma_rhs));

        alpha = p.chi;
        dP = full_gamma_step(1:n_state);
        df = reshape(full_gamma_step((n_state + 1):(2 * n_state)), size(P));
        dlambda = full_gamma_step(end);
        P = reshape(P(:) + alpha * dP, [], 3);
        f = f + alpha * df;
        lambda = lambda + alpha * dlambda;
        geo = Geometry(M, P);
        ls_iter = 1;






        % ... to here

        
        u = reshape((P(:) - P0) / p.dt, size(P));
        c = bie_residual(P_bie, M, f, geo_bie, u, u_background, slp_cache, p);
        eps_c_raw = norm(c(:)) / sqrt(numel(c));
        u_rms = norm(u(:)) / sqrt(numel(u));
        eps_c = eps_c_raw/ max(u_rms, 1e-14);

        darea = geo.area - p.area0;
        eps_d = abs(darea) / p.area0;
        lambda = lambda - o.eta * darea;
        b_res = force_balance_residual(P, P0, M, f, lambda, KTK, DTD, p);
        eps_b_raw = norm(b_res(:)) / sqrt(numel(b_res));
        eps_b = eps_b_raw;%/ max(u_rms, 1e-14);




        if verbose
            fprintf("t = %d, j = %d, ls_iter = %d, eps_b = %0.4g, eps_c = %0.4g, eps_d = %0.4g, darea = %0.4g \n", ...
                t, j, ls_iter, eps_b, eps_c, eps_d, darea);
        end

        if eps_b > eps_b_prev
            eps_b_rise_count = eps_b_rise_count + 1;
            if eps_b_rise_count > 2
                if ~supress_outputs
                    fprintf("eps_b monotonically increasing; halving o.h and restarting t = %d\n", t);
                end
                p.chi = 0.8 * p.chi;
                lambda = lambda_restart;
                f = f_restart;
                P = reshape(P0 + p.dt * velocity, [], 3);
                geo = Geometry(M, P);
                eps_b = Inf;
                eps_c = Inf;
                eps_d = Inf;
                eps_b_prev = Inf;
                eps_b_rise_count = 0;
                j = 0;
                alpha_mem = o.h;
                continue
            end
        else
            eps_b_rise_count = 0;
        end
        eps_b_prev = eps_b;

        j = j + 1;
    end

    if j >= o.max_iter
        warning("Terminating at t = %d because j reached o.max_iter = %d. eps_b = %0.4g, eps_c = %0.4g, eps_d = %0.4g", ...
            t, o.max_iter, eps_b, eps_c, eps_d);
        break;
    end

    p.total_time = p.total_time + toc;
    %[P, velocity] = rm_rigid(P, (P(:) - P0) / p.dt, geo.v_area);
    [P, velocity] = rm_rigid_patched(P, (P(:) - P0) / p.dt, geo.v_area,"translation");
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

function c_res = bie_residual(P, M, f, geo, u, u_background, slp_cache, p) %takes in  traction f, returns traction c
    slpout = stokeslet_SLP_triangle(P, M, f, slp_cache);
    normal_slip = p.Gamma - dot(f, geo.v_normal, 2);
    c_res = - u + u_background + p.Sd * slpout - p.Sd*p.Da*normal_slip.*geo.v_normal;
end

function b_res = force_balance_residual(P, P0, M, f, lambda, KTK, DTD, p) % takes in traction f, returns nodal b_res
    geo = Geometry(M, P);
    f_nodal = traction_to_nodal(f, geo);
    fb = geo.bending_force(1);
    twoHn = geo.lap * P;
    u = reshape((P(:) - P0) / p.dt, size(P));
    fv = reshape(-2 * (KTK + p.k * DTD) * u(:), size(P));
    fs = reshape(lambda * twoHn, size(P));
    f_mem = fv + fb + fs;
    b_res = -f_nodal - f_mem;
end

function nodal = traction_to_nodal(traction, geo)
    nodal = traction .* geo.v_area;
end

function traction = nodal_to_traction(nodal, geo)
    traction = nodal ./ geo.v_area;
end

function N = normal_projection_matrix(v_normal)
    n_v = size(v_normal, 1);
    rows = zeros(9 * n_v, 1);
    cols = zeros(9 * n_v, 1);
    vals = zeros(9 * n_v, 1);
    cursor = 0;

    for i = 1:n_v
        ni = v_normal(i, :).';
        block = ni * ni.';
        for row_dim = 1:3
            row = i + (row_dim - 1) * n_v;
            for col_dim = 1:3
                cursor = cursor + 1;
                rows(cursor) = row;
                cols(cursor) = i + (col_dim - 1) * n_v;
                vals(cursor) = block(row_dim, col_dim);
            end
        end
    end

    N = sparse(rows, cols, vals, 3 * n_v, 3 * n_v);
end
