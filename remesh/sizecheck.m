[P, M] = subdivided_sphere(25);
P(:,3) = 1.1*P(:,3);
geo=Geometry(M,P);
disp(mean(geo.he_length))