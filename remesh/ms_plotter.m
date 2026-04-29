close all;clc;

pnd =  0;
beta = (pnd/.9)*(-4e-8);
view_azi = 154; view_ele =21; alph = .9; %view_azi=-35;view_ele=-68;
edge = 2.5;
timetime=0;
gamy=-1;
Da = 1000;
Sd = 1;
rotate = true;
zoom=false;
show_multigrid=false;
cross_section=false; csx=.5;
color_mode = "curvature"; % options: "permeation_flux", "curvature", "uniform"
local_average_colors = true;
color_smooth_iters = 7;
color_smooth_blend = 1.0;
interpolate_vertex_colors = true;


%streamline bulk viz
show_bulk_streamlines = false;
bulk_grid_n = 40;
bulk_exclusion_scale = 1.1;
bulk_streamline_color = [0.1, 0.1, 0.1];
bulk_seed_n = 8;
bulk_seed_mode = "offset_surface"; % options: "offset_surface", "yz_plane"
bulk_seed_x = 1;
bulk_offset_scale = 1.8;
bulk_vertex_stride = 6;
bulk_seed_exclusion_scale = 0.6;
bulk_direction_arrows = true;
bulk_arrow_scale = 0.18;
bulk_arrow_stride = 10;
bulk_stream_step = 0.06;
bulk_stream_maxvert = 2000;

run_tag = sprintf('gamy_%+.2f_Da_%+.2f_Sd_%+.2f', gamy, Da, Sd);
run_tag = strrep(run_tag, '+', 'p');
run_tag = strrep(run_tag, '-', 'm');
name = "./data/" + run_tag + ".gif";
directory = "./data/ms_batch_data/";
delete(name)
gif_frame_count = 0;
size_x = 1000; size_y = 900;
set(gcf, 'Position',  [200, 200, size_x, size_y], ...
    'Color', 'w', ...
    'InvertHardcopy', 'off')
gif_delay = 1 / 14;
gif_frame_path = "./data/.ms_plotter_frame.png";
if isfile(gif_frame_path)
    delete(gif_frame_path)
end

%Determine initial face properties to track if the mesh effects the dynamics
geoj = load(directory + run_tag + sprintf("/geo%d.mat",1));
geoj_Geom = Geometry(geoj.M,geoj.P);

% Area change
area0 = geoj_Geom.area;

%determine last frame
folder = directory + run_tag; %change this to your actual path
files = dir(fullfile(folder, 'geo*.mat'));
tf = int32(length(files));
ti = 1; ts = int32(max(1, ceil(double(tf)/1000)));
% OVERRIDES



% Static coloration
geoj = load(directory + run_tag + sprintf("/geo%d.mat",tf)); 
geoj_Geom = Geometry(geoj.M,geoj.P);mcs = geoj_Geom.v_mean_curvature;
% disp(geoj_Geom.volume)
% disp(geoj_Geom.willmore_energy(1)-8*pi)
A0_range = mcs;
if local_average_colors
    A0_range = smooth_vertex_field(A0_range, geoj.M, color_smooth_iters, color_smooth_blend);
end
A0_min = min(A0_range); A0_max = max(A0_range);
if color_mode == "permeation_flux"
    [~, ~, ~, ~, KTK_final, DTD_final] = geoj_Geom.evolving_operators();
    twoHn_final = geoj_Geom.lap * geoj.P;
    fb_final = geoj_Geom.bending_force(1);
    fs_final = reshape(geoj.lambda * twoHn_final, size(geoj.P));
    fv_final = reshape((- 2 * (KTK_final + geoj.p.k * DTD_final) * geoj.velocity / geoj.p.dt), size(geoj.P));
    flux_range = geoj.p.pnd + dot(fb_final + fs_final + fv_final, geoj_Geom.v_normal, 2);
    if local_average_colors
        flux_range = smooth_vertex_field(flux_range, geoj.M, color_smooth_iters, color_smooth_blend);
    end
    flux_abs_max = max(abs(flux_range));
    if flux_abs_max == 0
        flux_abs_max = 1;
    end
end

% Static velo vector scale
tvelocity = reshape(geoj.velocity,[length(geoj.velocity)/3,3]) - (dot(reshape(geoj.velocity,[length(geoj.velocity)/3,3])',geoj_Geom.v_normal'))'.*geoj_Geom.v_normal;
tvelocity = tvelocity(:);
tv_max = max(tvelocity)*5;
tv_color = 'k';

%Uniform coloration
uni_color = [.6,.6,.6];
pi_nondim = pnd;
display_mode = 3;
mintf = 1;
customtime = [ti:10*ts:5320,5320:ts:tf];
for n = ti:ts:tf
    mintf = min(mintf,geoj.p.dt);
    hold off;
    ax = gca;
    set(ax, 'Units', 'normalized', 'Position', [0.03, 0.03, 0.94, 0.94], ...
        'ActivePositionProperty', 'position');
    geoj = load(directory + run_tag + sprintf("/geo%d.mat",n));
    L = length(geoj.P(:,1));
    geoj_Geom = Geometry(geoj.M,geoj.P);
    V = geoj_Geom.mesh.n_v;
    stretch_factor = 1;
    
    if true
        fi = trisurf(geoj.M,geoj.P(:,1),geoj.P(:,2),geoj.P(:,3), ...
            'FaceColor',uni_color,'Edgecolor',[.3,.3,.3], ...
            'FaceAlpha', alph); hold on;
    end
    %Display Settings

    %Face colors
    mcs = geoj_Geom.v_gaussian_curvature;
    A0 = mcs;
    if local_average_colors
        A0 = smooth_vertex_field(A0, geoj.M, color_smooth_iters, color_smooth_blend);
    end
    
    %Color faces
    if color_mode == "curvature"
        if interpolate_vertex_colors
            set(fi,'FaceColor','interp',...
               'FaceVertexCData',A0,...
               'CDataMapping','scaled');
            clim([A0_min, A0_max]);
        else
            A0_face = mean(A0(geoj_M_safe(geoj.M)), 2);
            set(fi,'FaceColor','flat',...
               'FaceVertexCData',(A0_face - A0_min)/(max(A0_max - A0_min, eps))*255,...
               'CDataMapping','direct');
        end
        %[v_ids,f_id] = get_surf_id(geoj.M,geoj.P);

        %set(fi,'FaceColor','flat',...
         %  'FaceVertexCData',255*(f_id)/2,...
          % 'CDataMapping','direct');
    elseif color_mode == "permeation_flux"
        [~, ~, ~, ~, KTK, DTD] = geoj_Geom.evolving_operators();
        twoHn = geoj_Geom.lap * geoj.P;
        fb = geoj_Geom.bending_force(1);
        fs = reshape(geoj.lambda * twoHn, size(geoj.P));
        fv = reshape((- 2 * (KTK + geoj.p.k * DTD) * geoj.velocity / geoj.p.dt), size(geoj.P));
        qperm = geoj.p.pnd + dot(fb + fs + fv, geoj_Geom.v_normal, 2);
        if local_average_colors
            qperm = smooth_vertex_field(qperm, geoj.M, color_smooth_iters, color_smooth_blend);
        end
        if interpolate_vertex_colors
            set(fi,'FaceColor','interp',...
               'FaceVertexCData',qperm,...
               'CDataMapping','scaled');
        else
            face_flux = mean(qperm(geoj_M_safe(geoj.M)), 2);
            set(fi,'FaceColor','flat',...
               'FaceVertexCData',face_flux,...
               'CDataMapping','scaled');
        end
        clim([-flux_abs_max, flux_abs_max]);
        colormap(turbo(256));
        colorbar;
    else
        set(fi,'FaceColor',uni_color);
    end

    if false
        set(gcf,'Renderer','opengl')
        [v_ids,f_id] = get_surf_id(geoj.M,geoj.P);
        % Surface 1
        idx1 = (f_id == 1);
       h1 = trisurf(geoj.M(idx1,:), geoj.P(:,1), geoj.P(:,2), geoj.P(:,3), ...
                'FaceColor','flat','FaceVertexCData',(A0(idx1,:) - A0_min)/(A0_max - A0_min)*255, 'FaceAlpha',.9, 'EdgeColor',[.1,.1,.1],'CDataMapping','direct'); hold on;

        % Surface 2
        idx2 = (f_id == 2);
        h2 = trisurf(geoj.M(idx2,:), geoj.P(:,1), geoj.P(:,2), geoj.P(:,3), ...
                'FaceColor','flat','FaceVertexCData',(A0(idx2,:) - A0_min)/(A0_max - A0_min)*255, 'FaceAlpha',1.0, 'EdgeColor',[.5,.5,.5],'CDataMapping','direct');
    end



    %lighting gouraud
% % lt1 = camlight(-149,-23);
% % lt2 = camlight("right","local");
% 
% %camlight("top")
% %light

% lt3 = camlight('headlight');
% material dull
% h2.SpecularStrength = .3;
% h2.AmbientStrength = .5;



        
    %Display velo?
    if true
        qq = quiver3(geoj.P(:,1),geoj.P(:,2),geoj.P(:,3),geoj.velocity(1:L),geoj.velocity(L+1:L*2),geoj.velocity(L*2+1:end),4,'k');
        %qq.Color = [.9, 0.9, 0.9];
    end

    %Display bulk streamlines seeded from a uniform y-z plane?
    if show_bulk_streamlines
        [~, ~, ~, ~, KTK, DTD] = geoj_Geom.evolving_operators();
        twoHn = geoj_Geom.lap * geoj.P;
        fb = geoj_Geom.bending_force(1);
        fs = reshape(geoj.lambda * twoHn, size(geoj.P));
        fv = reshape((- 2 * (KTK + geoj.p.k * DTD) * geoj.velocity / geoj.p.dt), size(geoj.P));
        f_mem = fb + fs + fv;

        x_grid = linspace(-stretch_factor*edge, stretch_factor*edge, bulk_grid_n);
        y_grid = linspace(-edge, edge, bulk_grid_n);
        z_grid = linspace(-edge, edge, bulk_grid_n);
        [Xg, Yg, Zg] = meshgrid(x_grid, y_grid, z_grid);
        targets = [Xg(:), Yg(:), Zg(:)];

        field_options = struct('near_scale', 3.0, 'very_near_scale', 1.2, 'chunk_size', 96);
        u_bulk = stokeslet_SLP_field(targets, geoj.P, geoj.M, f_mem, field_options) + shear_flow(targets, gamy);
        Ux = reshape(u_bulk(:,1), size(Xg));
        Uy = reshape(u_bulk(:,2), size(Yg));
        Uz = reshape(u_bulk(:,3), size(Zg));

        min_edge = min(geoj_Geom.he_length);
        exclusion_radius2 = (bulk_exclusion_scale * min_edge)^2;
        min_dist2 = inf(size(targets,1), 1);
        chunk_size = 256;
        for first = 1:chunk_size:size(geoj.P, 1)
            last = min(first + chunk_size - 1, size(geoj.P, 1));
            vertex_chunk = geoj.P(first:last, :);
            dx = targets(:,1) - reshape(vertex_chunk(:,1), 1, []);
            dy = targets(:,2) - reshape(vertex_chunk(:,2), 1, []);
            dz = targets(:,3) - reshape(vertex_chunk(:,3), 1, []);
            dist2_chunk = dx.^2 + dy.^2 + dz.^2;
            min_dist2 = min(min_dist2, min(dist2_chunk, [], 2));
        end
        mask = reshape(min_dist2 < exclusion_radius2, size(Xg));
        Ux(mask) = NaN;
        Uy(mask) = NaN;
        Uz(mask) = NaN;

        if bulk_seed_mode == "offset_surface"
            seed_idx = 1:bulk_vertex_stride:size(geoj.P, 1);
            offset_distance = max(bulk_offset_scale, bulk_exclusion_scale + 0.25) * min_edge;
            seed_targets = [
                geoj.P(seed_idx,:) + offset_distance * geoj_Geom.v_normal(seed_idx,:);
                geoj.P(seed_idx,:) - offset_distance * geoj_Geom.v_normal(seed_idx,:)
            ];
        else
            seed_y = linspace(min(y_grid), max(y_grid), bulk_seed_n);
            seed_z = linspace(min(z_grid), max(z_grid), bulk_seed_n);
            [seed_Y, seed_Z] = meshgrid(seed_y, seed_z);
            seed_X = bulk_seed_x * ones(size(seed_Y));
            seed_targets = [seed_X(:), seed_Y(:), seed_Z(:)];
        end
        min_seed_dist2 = inf(size(seed_targets,1), 1);
        for first = 1:chunk_size:size(geoj.P, 1)
            last = min(first + chunk_size - 1, size(geoj.P, 1));
            vertex_chunk = geoj.P(first:last, :);
            dx = seed_targets(:,1) - reshape(vertex_chunk(:,1), 1, []);
            dy = seed_targets(:,2) - reshape(vertex_chunk(:,2), 1, []);
            dz = seed_targets(:,3) - reshape(vertex_chunk(:,3), 1, []);
            dist2_chunk = dx.^2 + dy.^2 + dz.^2;
            min_seed_dist2 = min(min_seed_dist2, min(dist2_chunk, [], 2));
        end
        seed_exclusion_radius2 = (bulk_seed_exclusion_scale * min_edge)^2;
        seed_mask = min_seed_dist2 >= seed_exclusion_radius2;
        seed_x_valid = seed_targets(seed_mask,1);
        seed_y_valid = seed_targets(seed_mask,2);
        seed_z_valid = seed_targets(seed_mask,3);
        launch_ok = ...
            ~isnan(interp3(Xg, Yg, Zg, Ux, seed_x_valid, seed_y_valid, seed_z_valid, 'linear')) & ...
            ~isnan(interp3(Xg, Yg, Zg, Uy, seed_x_valid, seed_y_valid, seed_z_valid, 'linear')) & ...
            ~isnan(interp3(Xg, Yg, Zg, Uz, seed_x_valid, seed_y_valid, seed_z_valid, 'linear'));
        seed_x_valid = seed_x_valid(launch_ok);
        seed_y_valid = seed_y_valid(launch_ok);
        seed_z_valid = seed_z_valid(launch_ok);
        streams_fwd = stream3(Xg, Yg, Zg, Ux, Uy, Uz, ...
            seed_x_valid, seed_y_valid, seed_z_valid, [bulk_stream_step, bulk_stream_maxvert]);
        streams_bwd = stream3(Xg, Yg, Zg, -Ux, -Uy, -Uz, ...
            seed_x_valid, seed_y_valid, seed_z_valid, [bulk_stream_step, bulk_stream_maxvert]);
        for pass = 1:2
            if pass == 1
                stream_set = streams_fwd;
            else
                stream_set = streams_bwd;
            end
            for si = 1:numel(stream_set)
                stream = stream_set{si};
                if isempty(stream)
                    continue
                end
                plot3(stream(:,1), stream(:,2), stream(:,3), ...
                    'Color', bulk_streamline_color, 'LineWidth', 1.0);
                if bulk_direction_arrows && pass == 1 && size(stream,1) >= 2
                    idx = 1:bulk_arrow_stride:(size(stream,1) - 1);
                    anchor = stream(idx, :);
                    tangent = stream(idx + 1, :) - stream(idx, :);
                    tangent_norm = vecnorm(tangent, 2, 2);
                    keep = tangent_norm > 0;
                    if any(keep)
                        anchor = anchor(keep, :);
                        tangent = tangent(keep, :) ./ tangent_norm(keep);
                        quiver3(anchor(:,1), anchor(:,2), anchor(:,3), ...
                            bulk_arrow_scale * tangent(:,1), ...
                            bulk_arrow_scale * tangent(:,2), ...
                            bulk_arrow_scale * tangent(:,3), ...
                            0, 'Color', bulk_streamline_color, 'LineWidth', 0.8, 'MaxHeadSize', 2.0);
                    end
                end
            end
        end
    end

    %Display tangent velo?
    if false
        tvelocity = reshape(geoj.velocity,[length(geoj.velocity)/3,3]) - (dot(reshape(geoj.velocity,[length(geoj.velocity)/3,3])',geoj_Geom.v_normal'))'.*geoj_Geom.v_normal;
        tvelocity = tvelocity(:)/(tv_max);
        quiver3(geoj.P(:,1),geoj.P(:,2),geoj.P(:,3),tvelocity(1:L),tvelocity(L+1:L*2),tvelocity(L*2+1:end),5,tv_color)
    end

    %Display bending force?
    if false
        fb = geoj_Geom.bending_force;
        quiver3(geoj.P(:,1),geoj.P(:,2),geoj.P(:,3),fb(:,1),fb(:,2),fb(:,3),tv_color)
    end
        %Display osmotic force?
    if false
        weighted_normals = geoj_Geom.v_area.*geoj_Geom.v_normal;
        quiver3(geoj.P(:,1),geoj.P(:,2),geoj.P(:,3),weighted_normals(:,1),weighted_normals(:,2),weighted_normals(:,3),tv_color)
    end

    %Display repulsive force?
    if false
        dphi = TPE_grad_nsquared_gpt(geoj_Geom,geoj.p.co);
        uu = reshape(geoj.velocity,[V,3]);
        DtD = zeros([V,3,3]);
        fr = zeros([V,3]);
        for vi = 1:V
            DtD(vi,:,:) = dphi(vi,:)'*dphi(vi,:);
            fr(vi,:) = -reshape(dphi(vi,:)'*dphi(vi,:),[3,3])*uu(vi,:)';
        end
        quiver3(geoj.P(:,1),geoj.P(:,2),geoj.P(:,3),dphi(:,1),dphi(:,2),dphi(:,3),4,Color=[.9,.4,.4])
    end
    
    % Coloring



    %alpha(alph)
    timetime=timetime + ts*geoj.p.dt;
    xlim([-stretch_factor*edge,stretch_factor*edge]); ylim([-edge,edge]); zlim([-edge,edge]);
    %title(sprintf("o.h = %5f o.eta = %5f",geoj.o.h,geoj.o.eta))

    view(view_azi,view_ele)
    set(gca,'Color',[.9,.9,.9])
    set(gcf, 'Position',  [100, 100, size_x, size_y], ...
        'Color', 'w', ...
        'InvertHardcopy', 'off')
    axis manual;
    axis vis3d;
    pbaspect([1 1 1]);
    if true
        axis off;
        grid off;
    end
    
    %%VIEW OVERRIDE
    %ylim([-0.5 -.3]);zlim([-.6 .75])
    %view(0,0)

    if cross_section
        ylim([-2 csx]);%zlim([-.6 .75])
        view(-180,0)
    end
    drawnow;
    frame = getframe(gcf);
    frame_rgb = frame2im(frame);
    [frame_ind, map] = rgb2ind(frame_rgb, 256);
    if gif_frame_count == 0
        imwrite(frame_ind, map, name, 'gif', ...
            'LoopCount', inf, ...
            'DelayTime', gif_delay, ...
            'DisposalMethod', 'restoreBG');
    else
        imwrite(frame_ind, map, name, 'gif', ...
            'WriteMode', 'append', ...
            'DelayTime', gif_delay, ...
            'DisposalMethod', 'restoreBG');
    end
    gif_frame_count = gif_frame_count + 1;
    %disp((geoj_Geom.area - area0)/area0);
    %disp(geoj.p.dt)



end

%Rotation of final state
if rotate
    rspeed=0;
    for n = 1:15
        drawnow;
        frame = getframe(gcf);
        frame_rgb = frame2im(frame);
        [frame_ind, map] = rgb2ind(frame_rgb, 256);
        if gif_frame_count == 0
            imwrite(frame_ind, map, name, 'gif', ...
                'LoopCount', inf, ...
                'DelayTime', gif_delay, ...
                'DisposalMethod', 'restoreBG');
        else
            imwrite(frame_ind, map, name, 'gif', ...
                'WriteMode', 'append', ...
                'DelayTime', gif_delay, ...
                'DisposalMethod', 'restoreBG');
        end
        gif_frame_count = gif_frame_count + 1;
    end
end

function field_out = smooth_vertex_field(field_in, F, n_iters, blend)
    if n_iters <= 0 || blend <= 0
        field_out = field_in;
        return
    end

    adjacency = vertex_adjacency_matrix(F);
    degree = sum(adjacency, 2);
    field_out = field_in(:);

    for iter = 1:n_iters
        neighbor_sum = adjacency * field_out;
        neighbor_avg = field_out;
        has_neighbors = degree > 0;
        neighbor_avg(has_neighbors) = neighbor_sum(has_neighbors) ./ degree(has_neighbors);
        field_out = (field_out + blend * neighbor_avg) ./ (1 + blend);
    end
end

function adjacency = vertex_adjacency_matrix(F)
    n_v = max(F(:));
    edges = [F(:, [1 2]); F(:, [2 3]); F(:, [3 1])];
    edges = sort(edges, 2);
    edges = unique(edges, 'rows');
    adjacency = sparse(edges(:,1), edges(:,2), 1, n_v, n_v);
    adjacency = adjacency + adjacency.';
end

function F_out = geoj_M_safe(F_in)
    F_out = F_in;
end
