verbose = true;
time_readouts = false;
%%% directory
dir = "./data/willmore/genus6_smooth/"; 
[status, msg, msgID] = mkdir(dir);
run('./isoremesh/create_mex.m')


%% system init
%%% 0 for new simulation, otherwise continue from geo"start".mat
start = 0;

if start == 0
    %%% geometry
    %[P, M] = subdivided_sphere(4);
    load("./genus6_smooth.mat")
    %P(:, 1) = P(:, 1) * 2;
    geo = Geometry(M, P);
    %%% parameters
    p.dt = 1e-2; % time
    p.kappa = 5e-2; % bending
    p.T = 6000; % total time (frames)
    %%% optimizer
    o.h = 1; o.eta = 1e2; o.k = 0; o.metric = "lap";
    o.tol_f = 5e-4; o.tol_d = 5e-4; o.max_iter = 10000;
    %%% remesh
    r.edge_length = mean(geo.he_length); % target edge length
    r.n_iter = 50; % remeshing iterations
    r.smooth = 0; % smoothing factor
    %%% initialize
    velocity = zeros(size(P, 1) * 3, 1);
    pressure = zeros(size(M, 1), 1);
else
    load(dir + sprintf("geo%d.mat", start), ...
     "M", "P", "velocity", "pressure", "p", "o", "r"); 
    geo = Geometry(M, P);
end
%% REMOVE LATER!
%r.n_iter = 1250;
p.T = 6000;
%% main loop (not supposed to be modified)
for t = (start + 1):p.T
    tic;
    %%% incremental potential minimization
    [~, K, ~, div, KTK, DTD] = geo.evolving_operators();
    P0 = P(:); P(:) = P0 + p.dt * velocity;
    eps_f = Inf; eps_d = Inf; j = 0;
    while ((eps_f > o.tol_f) || (eps_d > o.tol_d)) && (j < o.max_iter)
        %%% update bending force
        fb = geo.bending_force(p.kappa);
        b = - 2 * KTK * (P(:) - P0) - p.dt * (- fb(:) - div' * pressure);
        expan = div * (P(:) - P0);
        %%% measure residual
        eps_f = norm_f(b, geo.v_area); eps_d = norm_d(expan, geo.f_area);
        if verbose fprintf("t = %d, j = %d, eps_f = %0.4g, eps_d = %0.4g \n", t, j, eps_f, eps_d); end
        %%% gradient descent/ascent
        P(:) = P(:) + o.h * ((preconditioner(geo, o.metric) + 1e-5 * sparse(eye(3*geo.mesh.n_v))) ...
                            \ (b - o.k * DTD * (P(:) - P0)));
        pressure = pressure - o.eta * o.h * expan ./ geo.f_area;
        %%% update geometry
        geo = Geometry(M, P);
        j = j + 1;
    end
    %%% save data
    [P, velocity] = rm_rigid(P, (P(:) - P0) / p.dt, geo.v_area);
    save(dir + sprintf("geo%d.mat", t), "M", "P", "velocity", "pressure", "fb", "p", "o", "r");
    if time_readouts
        fprintf('t= %i \n',t);
        fprintf('Iteration time: %0.4f \n',toc)
    else
        fprintf("Save geo%d.mat at j =%d, eps_f = %0.4g, eps_d = %0.4g \n", t, j, eps_f, eps_d);
    end
    %%% update geometry, remesh if needed
    geo = Geometry(M, P);
    if ~geo.is_delaunay(0)
        geo_pre = geo;
        if ~time_readouts
            fprintf("Remeshing. t = %d \n", t);
        end
        tic;
        %[P, M] = IO.remesh(P, M, r, dir, t); %Old Version using houdini
        [M, P] = remeshing(int32(M), P, int32([]), r.edge_length, int32(r.n_iter));
        geo = Geometry(M, P);
        [velocity, pressure] = map_data(geo, geo_pre, velocity, pressure);
        if time_readouts
            fprintf('Remeshing time: %0.4f \n',toc)
            fprintf('Newmesh delaunay? %i \n',geo.is_delaunay(0))
        end
    end
end
%% helper functions
function H = preconditioner(geo, metric)
    % preconditioner for incremental potential minimization
    switch  metric 
        case "bih"
            mass0 = spdiags(geo.v_area, 0, geo.mesh.n_v, geo.mesh.n_v);
            mass0_inv = spdiags(1./geo.v_area, 0, geo.mesh.n_v, geo.mesh.n_v);
            H = mass0 * geo.lap * mass0_inv * geo.lap;
            H = blkdiag(H, H, H);
        case "lap"
            H = geo.lap;
            H = blkdiag(H, H, H);
        case "stks"
            H = KTK;
        case "l2"
            H = sparse(eye(3*geo.mesh.n_v));
    end
end

function e = norm_f(b, area)
    % normalized norm of vector field
    e = sqrt(sum(b.^2 ./ [area; area; area])) / sum(area);
end

function e = norm_d(expan, area)
    % normalized norm of scalar field
    e = sqrt(sum(expan.^2 ./ area)) / sum(area);
end

function [velocity, pressure] = map_data(geo, geo_pre, velocity_pre, pressure_pre)
    % interpolate data from previous geometry to current geometry
    kdtree = KDTreeSearcher(geo_pre.f_center);
    %%% interpolate vertex data - velocity
    [face, uv, count, fail] = project(geo_pre.V, geo_pre.F, geo.V, kdtree, 6);
    if fail
        error("projection failed.");
    end
    velocity = interpolate(geo_pre.F, face, uv, reshape(velocity_pre, [], 3));
    velocity = velocity(:);
    %%% interpolate face data - pressure
    [face, uv, count, fail] = project(geo_pre.V, geo_pre.F, geo.f_center, kdtree, 6);
    if fail
        error("projection failed.");
    end
    [pressure_pre_v, neighbor] = geo_pre.mesh.face_to_vertex(pressure_pre);
    pressure_pre_v = pressure_pre_v ./ neighbor;
    pressure = interpolate(geo_pre.F, face, uv, pressure_pre_v);
end

function [velocity, pressure] = aug_lag_guess(face, vertex, parameter)
    % instantaneous Stokes solution using augmented lagrangian
    geo = Geometry(face, vertex);
    id = sparse(eye(3*size(vertex, 1)));
    force = geo.bending_force(parameter.kappa);
    pressure = zeros(size(face, 1), 1);
    [~, ~, ~, div, KTK, DTD] = geo.evolving_operators();
    for i=1:30
        b = (force(:) + div' * pressure);
        velocity = (KTK +  1e-5 * id + 1 * DTD) \ b;
        % [velocity, fail] = cgs(KTK + rho * DTD +  1e-2 * id , b, 1e-4, 10000, [], [], velocity);
        pressure = pressure - 1 * (div * velocity) ./ geo.f_area;
    end
    [~, velocity] = rm_rigid(vertex, velocity, geo.v_area);
end
