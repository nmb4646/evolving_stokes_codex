close all;

[P,M] = pearledVesicleMesh(.3,35,35);

save("budded_vesicle.mat", "M", "P")


figure
trisurf(M, P(:,1), P(:,2), P(:,3), ...
    'EdgeColor','none', 'FaceColor',[0.5 0.8 0.4]);
axis equal off
camlight
lighting gouraud