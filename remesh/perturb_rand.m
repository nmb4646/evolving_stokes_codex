function[geo_out] = perturb_rand(geo_in,factor)


P_out = geo_in.V + factor*randn(length(geo_in.V),1).*(geo_in.v_normal);












geo_out = Geometry(geo_in.F,P_out);

end