function [v] = reduced_volume(geo_in)
v = geo_in.volume/((4/3 * pi * sqrt(geo_in.area/(4*pi))^3));
end