function[geo_out] = create_inner_layer(geo_in)
M = geo_in.F; P = geo_in.V;
M2 = M;
%M2(:,2) = M(:,1); M2(:,1) = M(:,2); %permute indices to flip normal direction

P2 = P - .5*geo_in.v_normal;
M3 = [M;M2+length(P)];
P3 = [P;P2];
geo_out = Geometry(M3,P3);

end
