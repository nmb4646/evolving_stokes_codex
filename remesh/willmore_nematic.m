verbose = false;
%%% directory
addpath("evolve_nematic/")
dir = "./data/rod8/"; 
[status, msg, msgID] = mkdir(dir); 
%%% check remeshing 
hasRemesher = (exist('remeshing', 'file') ~= 0);
if ~hasRemesher
    warning('Remesher is not installed and has been disabled.');
end

%% system init
%%% 0 for new simulation, otherwise continue from geo"start".mat
start = 0;
if start == 0
    %%% geometry
    % load("../assets/torus.mat")
    [P, M] = subdivided_sphere(4);
    % [P, M] = loop(P, M, 1);
    % settings = struct('edge_length', 0.05, 'n_iter', 20, 'smooth', 0.1, ...
    %             'radius', 0.2, 'height', 4, 'n_col', 60);
    % [P, M] = IO.rod(settings, dir, false);
    fprintf("mesh obtained");
    geo = Geometry(M, P);
    p.rank = 2; % nematic init based on eigen rank
    %%% parameters
    p.dt = 0.1; % time
    p.kappa_inv = 1./100.; % viscosity
    p.alpha = 100; % activity: extensile (+), contractile (-)
    p.T = 5000; % total time (frames)
    %%% optimizer
    o.h = 500; o.eta = 500; o.k = 20; o.e = 1e-3; o.metric = "bih";
    o.tol_f = 1e-6; o.tol_d = 1e-5; o.max_iter = 10000;
    ls.c       = 0.5;   % small required decrease in E (monotone E decrease)
    ls.tau     = 0.5;     % backtracking factor
    ls.maxiter = 20;      % max backtracking iterations
    %%% remesh
    r.edge_length = mean(geo.he_length); % target edge length
    r.n_iter = 50; % remeshing iterations
    %%% initialize
    velocity = zeros(size(P, 1) * 3, 1);
    pressure = zeros(size(M, 1), 1);
    %%% initialize nematics
    nematic = initialize_nematic(M, P, p.rank);
    % for i = 1:30
    %     nematic = diffuse(geo, nematic, 0.05);
    % end
    % [phi_basis, theta_basis] = analytical_basis(geo.f_center, 1, 0.5);
    % nematic_ = phi_basis;
    % nematic_ = nematic_ ./ vecnorm(nematic_, 2, 2);
    % nematic = v2c(V2v(geo, nematic_)).^2;
else
    load(dir + sprintf("geo%d.mat", start), ...
     "M", "P", "velocity", "pressure", "nematic_", "p", "o", "r"); 
    geo = Geometry(M, P);
    nematic = v2c(V2v(geo, nematic_)).^2;
end

%% main loop (not supposed to be modified)
for t = (start + 1):p.T
    %%% incremental potential minimization
    [~, K, ~, div, KTK, DTD] = geo.evolving_operators();
    P0 = P(:); 
    % P(:) = P0 + p.dt * velocity;
    fa = sign(p.alpha) * double(ttt(K, sptensor(veronese(geo, nematic)), [1, 3, 5], [1, 2, 3]));
    eps_f = Inf; eps_d = Inf; j = 0;
    %%% energy wrapper (geo is updated in the loop)
    energy = @(pressure, geo) ...
        (geo.V(:) - P0)' * (p.kappa_inv * (KTK + o.k * DTD)) * (geo.V(:) - P0) ...
        - p.dt * (geo.V(:) - P0)' * (p.kappa_inv * fa(:) + div' * pressure) ...
        + p.dt * geo.willmore_energy(1);
    % energy = @(pressure, geo) + p.dt * geo.willmore_energy(1);
    use_precond = true;
    while ((eps_f > o.tol_f) || (eps_d > o.tol_d)) && (j < o.max_iter)
        H = preconditioner(geo, o.metric);
        % [~, notPD] = chol(H)
        % if notPD
        %     warning("not PD");
        % end
        %%% energy 
        E = energy(pressure, geo);
        %%% update bending force
        fb = geo.bending_force(1);
        % b = - p.dt * (- fb(:));
        b = - 2 * (p.kappa_inv * (KTK + o.k * DTD)) * (P(:) - P0) - p.dt * (- p.kappa_inv * fa(:) - fb(:) - div' * pressure); 
        % %%% gradient descent/ascent
        % P(:) = P0 + (H + 2 * o.h * KTK + o.k * DTD) \ (H * (P(:) - P0) + p.dt * o.h * (fa(:) + fb(:) + div' * pressure));

        % %%% compute gradient direction dP (full step)
        if use_precond
            % P1 = P(:); % 
            % P2 = P0 + (H + 2 * o.h * KTK + o.k * DTD) \ (H * (P(:) - P0) + p.dt * o.h * (fa(:) + fb(:) + div' * pressure));
            % dP = P2 - P1; 
            dP = H \ b;
            % dP = b;
        else
            dP = b;
        end

        %%%
        accepted = false;
        alpha = o.h;
        for ls_iter = 1:ls.maxiter %Line search
            geo_new = Geometry(M, reshape(P(:) + alpha * dP, [], 3));
            dE = E - energy(pressure, geo_new);
            % alpha * b' * dP
            % dP == b
            % alpha * norm(b)^2
            if (dE >= ls.c * alpha * b' * dP)
                P = geo_new.V;
                geo = geo_new;
                accepted = true;
                use_precond = true;
                break;
            else
                alpha = ls.tau * alpha;
            end
        end
        if ~accepted
            % Gradient descent fallback (steepest descent in force direction)
            % P(:) = P(:) + 1e-4 * norm(P(:)) / norm(b) * b;
            % geo = Geometry(M, P);      
            use_precond = false;
            % warning('Line search failed to accept, switched to gradient step\n');
        end

        expan = div * (P(:) - P0);
        pressure = pressure  - o.eta * expan;
        % pressure = (geo.mass2 + o.e * geo.f_lap) \ (geo.mass2 * (pressure - o.eta * expan));
        %%% measure residual
        % eps_f = norm_f(b, geo.v_area); eps_d = 0;
        eps_f = norm(dP); eps_d = norm(expan);
        %  norm_d(expan, geo.f_area);
        if verbose fprintf("t = %d, j = %d, ls_iter = %d,  E = %0.10g, eps_f = %0.4g, eps_d = %0.4g \n", t, j, ls_iter, E, eps_f, eps_d); end

        %%% update geometry
        % geo = Geometry(M, P);
        j = j + 1;
    end
    %%% diffusion
    % nematic = diffuse(geo, nematic, p.dt / (1e-8 + abs(p.alpha)));
    % nematic = (geo.mass2 + p.dt / (1e-8 + abs(p.alpha)) * geo.face_bochner_laplacian(2)) \ (geo.mass2 * nematic); 
    %%% advection
    nematic = advect(M, reshape(P0, [], 3), reshape(P, [], 3), nematic, 2);
    %%% normalize
    nematic = nematic ./ abs(nematic);
    %%% save data
    [P, velocity] = rm_rigid(P, (P(:) - P0) / p.dt, geo.v_area);
    nematic_ = v2V(geo, c2v(nematic.^(1/2)));
    save(dir + sprintf("geo%d.mat", t), "M", "P", "velocity", "pressure", "nematic_", "fa", "fb", "p", "o", "r");
    fprintf("Save geo%d.mat at j =%d, eps_f = %0.4g, eps_d = %0.4g \n", t, j, eps_f, eps_d);
    %%% update geometry, remesh if needed
    geo = Geometry(M, P);
    if ~geo.is_delaunay(5e-2) && hasRemesher
        geo_pre = geo; 
        fprintf("Remeshing. t = %d \n", t);
        [M, P] = remeshing(int32(M), P, int32([]), r.edge_length, int32(r.n_iter)); M = cast(M, "double");
        geo = Geometry(M, P);
        [velocity, pressure, nematic] = map_data(geo, geo_pre, velocity, pressure, nematic);
    end               
end

%% helper functions
function H = preconditioner(geo, metric)
    % preconditioner for incremental potential minimization
    mass0 = spdiags(geo.v_area, 0, geo.mesh.n_v, geo.mesh.n_v);
    mass0_inv = spdiags(1./geo.v_area, 0, geo.mesh.n_v, geo.mesh.n_v);
    switch  metric 
        % case "bih"
        %     H = geo.lap * mass0_inv * geo.lap + mass0;
        %     H = blkdiag(H, H, H);
        case "bih"
            H = mass0_inv * geo.lap * mass0_inv * geo.lap;
            H = blkdiag(H, H, H) + speye(3*geo.mesh.n_v);
        % case "lap"
        %     H = geo.lap + mass0;
        %     H = blkdiag(H, H, H);
        case "lap"
            H = mass0_inv * geo.lap;
            H = blkdiag(H, H, H) + speye(3*geo.mesh.n_v);
        case "mass"
            H = blkdiag(mass0, mass0, mass0);
        case "l2"
            H = speye(3*geo.mesh.n_v);
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

function q_bra_hat = initialize_nematic(face, vertex, rank)
    % initialize nematic field by rayleigh quotient
    k = 2; %% k-atic
    geo = Geometry(face, vertex);
    L = geo.face_bochner_laplacian(k);
    [V, ~] = eigs(L, geo.mass2, rank, "smallestabs");
    q_bra_hat = V(:, rank) ./ vecnorm(V(:, rank), 2, 2);
end

function q_bra_hat = initialize_nematic2(face, vertex)
end

function [velocity, pressure, q_bra_hat] = map_data(geo, geo_pre, velocity_pre, pressure_pre, q_bra_hat_pre)
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
    %%% interpolate face data - nematic
    QQ_pre = veronese(geo_pre, q_bra_hat_pre);
    [QQ_pre_v, neighbor] = geo_pre.mesh.face_to_vertex(reshape(QQ_pre, geo_pre.mesh.n_f, []));
    QQ_pre_v = QQ_pre_v ./ neighbor;
    QQ = interpolate(geo_pre.F, face, uv, QQ_pre_v);
    q_bra_hat = veronese_inv(geo, QQ);
end

function QQ = veronese(geo, q_bra_hat)
    % veronese map from complex to tensor
    Q = v2V(geo, c2v(q_bra_hat.^(1/2)));
    QQ = Q(:, :, ones(1,3)) .* permute(Q(:,:,ones(1,3)), [1,3,2]);
end

function q_bra_hat = veronese_inv(geo, QQ)
    % inverse veronese map from tensor to complex
    eigen = @(i) eigs(reshape(QQ(i, :), 3, 3), 1 ,'largestabs');
    [Q, ~] = arrayfun(eigen, 1:geo.mesh.n_f, 'UniformOutput', false);
    Q = cell2mat(Q)';
    q = V2v(geo, Q);
    q_bra_hat = v2c(q ./ vecnorm(q, 2, 2)).^2;
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

function Fq_bra_hat = advect(mesh, vertex1, vertex2, q_bra_hat, k)
    % Lie advection of nematic field
    q_hat = q_bra_hat.^(1/k);
    F = pushforward(mesh, vertex1, vertex2);
    q = c2v(q_hat);
    Fq = squeeze(pagemtimes(permute(F, [2, 3, 1]), ...
                           permute(q, [2, 3, 1])))';
    Fq = Fq ./ vecnorm(Fq, 2, 2); % normalize
    Fq_bra_hat = v2c(Fq).^k;
end

function nematic = diffuse(geo, nematic, dt)
    QQ = V2Q(v2V(geo, c2v(nematic.^(1/2))));
    QQ = reshape(QQ, size(QQ, 1), []);
    QQ = (geo.mass2 + dt * geo.f_lap) \ (geo.mass2 * QQ);
    QQ = realign(geo, QQ);
    QQ = reshape(QQ, size(QQ, 1), 3, 3);
    nematic = v2c(V2v(geo, Q2V(geo, QQ))).^2;
end

function q = V2v(geo, Q)
    % R3 realization to local chart
    q = [dot(geo.f_basis_u, Q, 2), dot(geo.f_basis_v, Q, 2)];
end

function Q = v2V(geo, q)
    % R3 realization of tangent vector
    Q = geo.f_basis_u .* q(:, 1) + geo.f_basis_v .* q(:, 2);
end

function q = c2v(q_hat)
    % complex to vector
    q = [real(q_hat), imag(q_hat)];
end
function q_hat = v2c(q)
    % vector to complex
    q_hat = q(:,1) + 1i * q(:,2);
end

function QQ = realign(geo, QQ)
    % make QQ rank one 
    QQ = reshape(QQ, size(QQ, 1), 3, 3);
    QQ = V2Q(Q2V(geo, QQ));
    QQ = reshape(QQ, size(QQ, 1), []);
end

function Q = Q2V(geo, QQ)
    % map from QQ = s(Q'Q - I/3) to Q

    % QQ = rm_trace(QQ);
    eigen = @(i) eigs(reshape(QQ(i, :), 3, 3), 1 ,'largestabs');
    [Q, ~] = arrayfun(eigen, 1:size(QQ, 1), 'UniformOutput', false);
    Q = cell2mat(Q)';
    Q = Q - dot(Q, geo.f_normal, 2) .* geo.f_normal;
    Q = Q ./ vecnorm(Q, 2, 2);
end

function QQ = V2Q(Q)
    QQ = Q(:, :, ones(1,3)) .* permute(Q(:,:,ones(1,3)), [1,3,2]);
    % QQ = rm_trace(QQ);
end

function [P, theta, phi] = torus_project(P, R, r)
    % project points to torus, R is the major radius, r is the minor radius
    % phi is the angle around the torus, theta is the angle around the minor circle

    theta = atan2(P(:, 2), sqrt(P(:, 1) .* P(:, 1) + P(:, 3) .* P(:, 3)) - R);
    phi = atan2(P(:, 3), P(:, 1));
    P(:, 1) = (R + r .* cos(theta)) .* cos(phi);
    P(:, 3) = (R + r .* cos(theta)) .* sin(phi);
    P(:, 2) = r * sin(theta);
end

function [phi_basis, theta_basis] = analytical_basis(P, R, r)
    % analytically obtain the basis of a torus

    [P, theta, phi] = torus_project(P, R, r);
    phi_basis = [(R + r .* cos(theta)) .* (- sin(phi)), zeros(size(P, 1), 1), (R + r .* cos(theta)) .* cos(phi) ];
    theta_basis = [r .* (-sin(theta)) .* cos(phi), r * cos(theta), r .* (-sin(theta)) .* sin(phi)];
    % phi_basis = phi_basis ./ vecnorm(phi_basis, 2, 2);
    % theta_basis = theta_basis ./ vecnorm(theta_basis, 2, 2);
end
