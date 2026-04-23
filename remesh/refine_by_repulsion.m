function [multigrid_out] = refine_by_repulsion(geo_in,p_in) %in - Nc X Nc, can be a geometry

    pie = p_in.pnd*(16*pi*p_in.kappa);
    
    %%% update bending force
    fb = geo_in.bending_force(p_in.kappa);
    
    %%% repulsive force
    dphi = TPE_grad_truncated(geo_in,3,.2);
    fr = 250*p_in.beta*dphi;
    
    %%% functional gradient
    weighted_normals = geo_in.v_area.*geo_in.v_normal; %L^2
    fp = pie * weighted_normals;
    
    norm_matrix = 0*fb;
    rgb_forces = double(0*geo_in.F);
    for i = 1:length(fb)
        norm_matrix(i,1) = norm(fr(i,:));
        norm_matrix(i,2) = norm(fp(i,:));
        norm_matrix(i,3) = norm(fb(i,:)); 
    end
    
    for f = 1:length(geo_in.F)
        for vind = geo_in.F(f,:)
            rgb_forces(f,:) = rgb_forces(f,:) + norm_matrix(vind,:);
        end
        rgb_forces(f,:) = rgb_forces(f,:)/max(rgb_forces(f,:));
    end
    
    ftrlist = 1:length(geo_in.F);
    ftrlist = sort(ftrlist(rgb_forces(:,1)>.1),'descend')';
    
    if ~isempty(ftrlist)
        multigrid_out =refine_sub(geo_in,ftrlist);
    else

        multigrid_out.M = geo_in.F; 
        multigrid_out.P = geo_in.V; 
        multigrid_out.n = geo_in.v_normal; 
        multigrid_out.a = geo_in.v_area;
    end

end