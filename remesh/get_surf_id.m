
function [v_id_out, f_id_out] = get_surf_id(M_in,P_in)

Nv = size(P_in,1);

% --- Build vertex adjacency ---
I = [M_in(:,1); M_in(:,2); M_in(:,3)];
J = [M_in(:,2); M_in(:,3); M_in(:,1)];

A = sparse([I; J], [J; I], true, Nv, Nv);

% --- Connected components on vertices ---
v_id_out = conncomp(graph(A))';   % Nv x 1

% --- Face IDs (inherit from any vertex) ---
f_id_out = v_id_out(M_in(:,1));          % Nf x 1

%trisurf(M,P(:,1),P(:,2),P(:,3),'FaceVertexCData',f_id)

end