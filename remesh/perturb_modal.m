function[geo_out] = perturb_modal(geo_in,mode,factor)

[phi,theta,r] = cart2sph(geo_in.V(:,1),geo_in.V(:,2),geo_in.V(:,3));

di = [phi,-theta + pi/2];
%Decompose into harmonics

L_max = 10;
Yn = getSH(L_max,di,'real');
clm = Yn\r;
r_rec = zeros;
for n = 1:length(clm)
    r_rec = r_rec + clm(n)*Yn(:,n);
end
norm(abs(r_rec) -r);

l_idx = mode*(mode+1)+1;
[P_out(:,1),P_out(:,2),P_out(:,3)] = sph2cart(phi,theta,r + factor*Yn(:,l_idx));

geo_out = Geometry(geo_in.F,P_out);

end