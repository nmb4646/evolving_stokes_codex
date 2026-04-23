clear;clc;
load("./data/qs_batch_data/pnd_2.43/geo1.mat")
% geo_pre = Geometry(M,P);
% 
% P_max = max(max(M(1:4476)));
% v_ids = zeros([length(P),1]);
% v_ids(1:P_max)=1;
% f_ids = zeros([length(M),1]);
% for f = 1:length(M)
%     points = M(f,:);
%     f_ids(f) = mean(v_ids(points));
% end
% save("./data/qs_batch_data/pnd_2.43/geo1.mat", "M", "P", "velocity", "lambda", "o", "r","p","v_ids","f_ids");
% v_ids_old = v_ids;
% 
% 
% for n=3:600
% load(sprintf("./data/qs_batch_data/pnd_2.43/geo%d.mat",n))
% v_ids = v_ids_old;
% geo = Geometry(M,P);
% kdtree = KDTreeSearcher(geo_pre.f_center);
% [face, uv, count, fail] = project(geo_pre.V, geo_pre.F, geo.V, kdtree, 16);
% 
% 
% v_ids_new = interpolate(geo_pre.F, face, uv, v_ids);
% f_ids = zeros([length(M),1]);
% %v_to_f
% for f = 1:length(M)
%     points = M(f,:);
%     f_ids(f) = mean(v_ids_new(points));
% end
% 
% v_ids_old = v_ids_new;
% geo_pre = geo;
% 
% save(sprintf("./data/qs_batch_data/pnd_2.43/geo%d.mat",n), "M", "P", "velocity", "lambda", "o", "r","p","v_ids","f_ids");
% end

%trisurf(M,P(:,1),P(:,2),P(:,3),'FaceVertexCData',f_ids_new)


load("./data/qs_batch_data/pnd_2.43/geo600.mat")

Nv = size(P,1);

% --- Build vertex adjacency ---
I = [M(:,1); M(:,2); M(:,3)];
J = [M(:,2); M(:,3); M(:,1)];

A = sparse([I; J], [J; I], true, Nv, Nv);

% --- Connected components on vertices ---
v_id = conncomp(graph(A))';   % Nv x 1

% --- Face IDs (inherit from any vertex) ---
f_id = v_id(M(:,1));          % Nf x 1

trisurf(M,P(:,1),P(:,2),P(:,3),'FaceVertexCData',f_id)


function velocity = map_id(geo, geo_pre, id_pre)
    % interpolate data from previous geometry to current geometry
    kdtree = KDTreeSearcher(geo_pre.f_center);
    %%% interpolate vertex data - velocity
    [face, uv, count, fail] = project(geo_pre.V, geo_pre.F, geo.V, kdtree, 6);
    if fail
        error("projection failed.");
        disp("projection failed.")
    end
    velocity = interpolate(geo_pre.F, face, uv, reshape(velocity_pre, [], 3));
    velocity = velocity(:);
end