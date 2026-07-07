function[a] = excess_area(geo_in)
a = geo_in.area/((3*geo_in.volume/(4*pi))^(2/3)) - 4*pi;
end