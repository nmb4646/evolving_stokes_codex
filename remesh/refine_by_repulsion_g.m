function geo_out = refine_by_repulsion_g(geo_in, p_in)
    geo_out = geo_in;
    pie = p_in.pnd*(16*pi*p_in.kappa);
    
    fb = geo_out.bending_force(p_in.kappa);
    dphi = TPE_grad_truncated(geo_out, p_in.co,.2);
    fr = p_in.beta * dphi;

    fp = pie * (geo_out.v_area .* geo_out.v_normal);

    % compute per-vertex norms
    norm_matrix = zeros(size(fb));
    for i = 1:length(fb)
        norm_matrix(i,1) = norm(fr(i,:));
        norm_matrix(i,2) = norm(fp(i,:));
        norm_matrix(i,3) = norm(fb(i,:));
    end

    % per-face "rgb" force magnitude
    rgb_forces = zeros(size(geo_out.F));
    for f = 1:size(geo_out.F,1)
        for v = geo_out.F(f,:)
            rgb_forces(f,:) = rgb_forces(f,:) + norm_matrix(v,:);
        end
        rgb_forces(f,:) = rgb_forces(f,:) / max(rgb_forces(f,:));
    end

    % select faces where repulsive force dominates
    ftrlist = find(rgb_forces(:,1) == 1);

    if ~isempty(ftrlist)
        geo_out = refine_faces_g(geo_out, ftrlist); % batch refinement
    end
end