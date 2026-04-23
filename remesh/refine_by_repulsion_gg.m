function [M_refined, P_refined] = refine_by_repulsion_gg(M, P, p_in)
    % p_in: struct with fields pnd, kappa, beta, co (parameters)
    % M: Nx3 vertices
    % P: Mx3 faces

    % Step 1: Compute bending force (fb), repulsion force (fr), pressure force (fp)
    geo_00 = Geometry(M,P);
    pie = p_in.pnd*(16*pi*p_in.kappa);
    
    fb = geo_00.bending_force(p_in.kappa);
    dphi = TPE_grad_nsquared_gpt(geo_00, p_in.co);
    fr = p_in.beta * dphi;

    fp = pie * (geo_00.v_area .* geo_00.v_normal);

    % Step 2: Compute per-vertex norms of forces
    norm_matrix = zeros(size(fb,1), 3);
    norm_matrix(:,1) = vecnorm(fr, 2, 2);
    norm_matrix(:,2) = vecnorm(fp, 2, 2);
    norm_matrix(:,3) = vecnorm(fb, 2, 2);

    % Step 3: Compute per-face RGB force magnitude sums
    rgb_forces = zeros(size(P,1), 3);
    for f = 1:size(P,1)
        verts = P(f,:);
        rgb_forces(f,:) = sum(norm_matrix(verts, :), 1);
    end

    % Normalize each rgb vector by its max component
    rgb_forces = rgb_forces ./ max(rgb_forces, [], 2);

    % Step 4: Select faces where repulsive force dominates (normalized to 1)
    ftrlist = find(rgb_forces(:,1) == 1);

    % Step 5: Refine those faces using refine_sub
    if ~isempty(ftrlist)
        [M_refined, P_refined] = refine_sub_gg(M, P, ftrlist);
    else
        M_refined = M;
        P_refined = P;
    end
end
