% verbose = false;
%---  directory
% dir = "./data/tests/rod9/"; 
addpath("evolve_nematic/")
[status, msg, msgID] = mkdir(dir); 
%---  check remeshing 
hasRemesher = (exist('remeshing', 'file') ~= 0);
if ~hasRemesher
    warning('Remesher is not installed and has been disabled.');
end

%---  system init
%---  0 for new simulation, otherwise continue from geo"start".mat
% start = 0; 
if start == 0
    %---  geometry
    [P, M] = subdivided_sphere(4);
    geo = Geometry(M, P);
    %---  parameters
    p.dt = 0.05; % time step
    p.k = 10; % bulk viscosity 
    % p.kappa = 50; % rigidity [>0, 50]
    % p.alpha = Inf; % activity \pm [0.5, Inf]: extensile (+), contractile (-) 
    p.T = 700; % total time (frames)
    p.area0 = geo.area; % area 
    %---  optimizer
    o.h = 1; % optimizer step
    o.eta = 1; % multiplier step
    o.tol_f = 1e-4; % force tolerance 
    o.tol_d = 1e-3; % constraint tolerance
    o.max_iter = 10000; % maximum iterations
    %---  line search 
    lsr.c       = 1e-4; % small required decrease in E (monotone E decrease)
    lsr.tau     = 0.5; % backtracking factor
    lsr.max_iter = 20; % max backtracking iterations
    %---  remesh
    r.edge_length = mean(geo.he_length); % target edge length
    r.n_iter = 50; % remeshing iterations
    %---  initialize
    velocity = zeros(size(P, 1) * 3, 1);
    lambda = 0.;
    fa = zeros(size(P, 1) * 3, 1);
    fb = geo.bending_force(1);
    %---  initialize nematics
    nematic = initialize_nematic(M, P, 1);
    for i = 1:20
        nematic = diffuse(geo, nematic, 0.1);
    end
    nematic_ = v2V(geo, c2v(nematic.^(1/2)));
    save(dir + "geo0.mat", "M", "P", "velocity", "lambda", "nematic_", "fa", "fb", "p", "o", "r");
    fprintf("Save geo0.mat \n");
else
    load(sprintf("%sgeo%d.mat", dir, start), ...
     "M", "P", "velocity", "lambda", "nematic_", "p", "o", "r"); 
    geo = Geometry(M, P);
    nematic = v2c(V2v(geo, nematic_)).^2;
    p.T = 700;
    %---  line search 
    lsr.c       = 1e-4; % small required decrease in E (monotone E decrease)
    lsr.tau     = 0.5; % backtracking factor
    lsr.max_iter = 20; % max backtracking iterations
end

%% main loop (not supposed to be modified)
for t = (start + 1):p.T
    %---  initialization
    [~, K, ~, div, KTK, DTD] = geo.evolving_operators();
    P0 = P(:); 
    fa = sign(p.alpha) * double(ttt(K, sptensor(veronese(geo, nematic)), [1, 3, 5], [1, 2, 3]));
    eps_f = Inf; eps_d = Inf; j = 0;
    alpha_mem   = o.h; 

    %---  energy functions
    lagrangian = @(P_, willmore_, lambda_, area_) (P_ - P0)' * (KTK + p.k * DTD) * (P_ - P0) ...
                + p.kappa * p.dt * willmore_ ...
                - p.dt * (P_ - P0)' * fa(:) ...
                - p.dt * (area_ - p.area0) * lambda_;

    %---  approx Hessian
    mass0 = spdiags(geo.v_area, 0, geo.mesh.n_v, geo.mesh.n_v);
    mass0_inv = spdiags(1./geo.v_area, 0, geo.mesh.n_v, geo.mesh.n_v);
    bih = geo.lap * mass0_inv * geo.lap;
    bih = blkdiag(bih, bih, bih);
    Hess = 2 * (KTK + p.k * DTD) ...
        + 0.5 * p.kappa * p.dt * bih ...
        + 1e-3 * r.edge_length^(-2) * blkdiag(mass0, mass0, mass0);

    %---  schur complement 1: guess pressure
    %{
        schur = @(x) p.dt * div * (Hess \ (div' * x)); 
        rhs = div * P0 ...
            - p.dt *  p.kappa_inv * div * (Hess \ fa(:)) ...
                - 2 * p.kappa_inv * div * (Hess \ ((KTK + p.k * DTD) * P0));
                [pressure,~] = pcg(schur, rhs, 1e-1, 50); 
    %}
    %---  schur complement 1: guess position
    %{
        rhs = p.dt * p.kappa_inv * fa(:) ...
            + 2 * p.kappa_inv * (KTK + p.k * DTD) * P0 ...
                + p.dt * div' * pressure;
                [P(:),~] = pcg(Hess, rhs, 1e-5, 200);  
        geo = Geometry(M, P);
    %}

    %--- guess position with velocity 
    P(:) = P0 + p.dt * velocity;
    geo = Geometry(M, P);

    while ((eps_f > o.tol_f) || (eps_d > o.tol_d)) && (j < o.max_iter)
        %---  gradient descent 1: energy
        E = lagrangian(P(:), geo.willmore_energy(1), lambda, geo.area);
        %---  gradient descent 2: gradient 
        fb = geo.bending_force(1);
        twoHn = geo.lap * P;
        b = - 2 * (KTK + p.k * DTD) * (P(:) - P0) ...
            + p.kappa * p.dt * fb(:) ...
            + p.dt * fa(:) ...
            + p.dt * lambda * twoHn(:);
    
        %---  gradient descent 3: precondition
        %{
            [dP,FLAG,RELRES,ITER] = pcg(Hess,b,1e-3, 100);
            if FLAG ~= 0
                fprintf("FLAG: %d, residual: %0.4g, iter: %d \n", FLAG, RELRES,ITER);
            end
        %}
        dP = Hess \ b;

        %---  force norm 
        quad = dP' * b;
        if quad <= 0 
            message = sprintf('quad = %0.4g should be positive definite', quad);
            error(message);
        else
            eps_f = sqrt(quad);
        end
        
        %---  gradient descent 4: Armijo backtracking with lambda fixed  
        accepted = false;
        alpha = min(alpha_mem, 1 / max(1, eps_f));   % adjust alpha based on force norm 
        for ls_iter = 1:lsr.max_iter
            P_new = reshape(P(:) + alpha * dP, [], 3);
            geo_new = Geometry(M, P_new);
            dE = E - lagrangian(P_new(:), geo_new.willmore_energy(1), lambda, geo_new.area);
            pred = alpha * quad;
            if dE >= lsr.c * pred
                P = P_new;
                geo = geo_new;
                accepted = true;
                alpha_mem   = min(alpha/lsr.tau, o.h);  % grow next try
                break;
            else          
                alpha = lsr.tau * alpha;  % shrink
            end
        end
        if ~accepted
            error("backtracking fails");
        end
        %---  gradient ascent 1: update pressure
        %{
            expan = (div * (P(:) - P0)) ./ geo.f_area;
            eps_d = sqrt(mean(expan.^2));  % RMS, area-weighted
            pressure = pressure - o.eta * expan;
        %}
        %---  gradient ascent 2: update tension
        darea = geo.area - p.area0;
        eps_d = abs(darea);
        lambda = lambda - o.eta * darea;
        %---  report
        if verbose 
            fprintf("t = %d, j = %d, ls_iter = %d,  E = %0.10g, eps_f = %0.4g, eps_d = %0.4g, darea = %0.4g \n", t, j, ls_iter, E, eps_f, eps_d, darea); 
        end
        %---  index
        j = j + 1;
    end
    %--- stop integration if not converging
    if j >= o.max_iter
        warning("Terminating at t = %d because j reached o.max_iter = %d. eps_f = %0.4g, eps_d = %0.4g", ...
                 t, o.max_iter, eps_f, eps_d);
        break;  % breaks out of the time loop
    end
    %---  diffusion
    if abs(p.alpha) < Inf
        nematic = diffuse(geo, nematic, p.dt / (eps + abs(p.alpha)));
        % nematic = (geo.mass2 + p.dt / (1e-8 + abs(p.alpha)) * geo.face_bochner_laplacian(2)) \ (geo.mass2 * nematic); 
    end
    %---  advection
    nematic = advect(M, reshape(P0, [], 3), reshape(P, [], 3), nematic, 2);
    %---  normalize
    nematic = nematic ./ abs(nematic);
    %---  save data
    [P, velocity] = rm_rigid(P, (P(:) - P0) / p.dt, geo.v_area);
    geo.V = P;
    nematic_ = v2V(geo, c2v(nematic.^(1/2)));
    save(dir + sprintf("geo%d.mat", t), "M", "P", "velocity", "lambda", "nematic_", "fa", "fb", "p", "o", "r");
    fprintf("Save geo%d.mat at j = %d, eps_f = %0.4g, eps_d = %0.4g \n", t, j, eps_f, eps_d);
    %---  remesh if needed
    if ~geo.is_delaunay(1e-2) && hasRemesher
        geo_pre = geo; 
        fprintf("Remeshing. t = %d \n", t);
        [M, P] = remeshing(int32(M), P, int32([]), r.edge_length, int32(r.n_iter)); M = cast(M, "double");
        geo = Geometry(M, P);
        [velocity, nematic] = map_data(geo, geo_pre, velocity, nematic);
    end               
end

%% helper functions

function q_bra_hat = initialize_nematic(face, vertex, rank)
    % initialize nematic field by rayleigh quotient
    k = 2; %% k-atic
    geo = Geometry(face, vertex);
    L = geo.face_bochner_laplacian(k);
    [V, ~] = eigs(L, geo.mass2, rank, "smallestabs");
    q_bra_hat = V(:, rank) ./ vecnorm(V(:, rank), 2, 2);
end

function [velocity, q_bra_hat] = map_data(geo, geo_pre, velocity_pre, q_bra_hat_pre)
    % interpolate data from previous geometry to current geometry
    kdtree = KDTreeSearcher(geo_pre.f_center);
    %---  interpolate vertex data - velocity
    [face, uv, count, fail] = project(geo_pre.V, geo_pre.F, geo.V, kdtree, 6);
    if fail
        error("projection failed.");
    end
    velocity = interpolate(geo_pre.F, face, uv, reshape(velocity_pre, [], 3));
    velocity = velocity(:);
    %---  interpolate face data
    [face, uv, count, fail] = project(geo_pre.V, geo_pre.F, geo.f_center, kdtree, 6);
    if fail
        error("projection failed.");
    end
    %---  interpolate face data - pressure
    %{
        [pressure_pre_v, neighbor] = geo_pre.mesh.face_to_vertex(pressure_pre);
        pressure_pre_v = pressure_pre_v ./ neighbor;
        pressure = interpolate(geo_pre.F, face, uv, pressure_pre_v);
    %}
    %---  interpolate face data - nematic
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