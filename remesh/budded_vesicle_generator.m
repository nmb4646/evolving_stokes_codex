close all; clear;

[P,M] = pearledVesicleMesh(.3,20,20);
%[M,P] = initial_dumbbell(.35);
save("budded_vesicle.mat", "M", "P")


figure
trisurf(M, P(:,1), P(:,2), P(:,3), ...
    'EdgeColor',[0,0,0], 'FaceColor',[0.5 0.8 0.4]);
axis equal off
camlight
lighting gouraud