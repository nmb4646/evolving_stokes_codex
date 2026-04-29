close;clear;clc
[P1, M] = subdivided_sphere(12);
geo1 = Geometry(M, P1);

P2 = .9*P1; %+ .0001*randn(size(P1));
geo2 = Geometry(M, P2);

trisurf(M,P2(:,1),P2(:,2),P2(:,3))

norm(dot(geo1.bending_force(1)./geo1.v_area,P2-P1,2))

geo2.willmore_energy(1)/geo2.area-geo1.willmore_energy(1)/geo1.area
