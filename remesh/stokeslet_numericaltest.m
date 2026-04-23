clear;clc;close all;

%%% Testing the numerical implementation of the stokeslet operator

n_vals = 3:17; dx = [];trans_error=[]; normal_error = []; rot_error=[]; comb_error=[]; lin_error=[]; cubenorm_error =[];
for n = n_vals
    [P, M] = subdivided_sphere(n);%P = rfpza(P);
    N = length(P);
    geo=Geometry(M,P); 

    % trans test field: f = constant, u= 
    f = zeros(size(P));
    F = 23/7;    
    for k = 1:length(P)
    f(k,:) = [0,0,F];end

    slpcache = stokeslet_SLP_triangle_setup(M);

    u = stokeslet_SLP_triangle(P,M,f,slpcache);

    dx = [dx; mean(geo.he_length)]; 
    trans_error = [trans_error;max(vecnorm(u - (2/3)*f,2,2))/(F*1.5)];
    N = length(P);

    %normal test field
    f = geo.v_normal;

    u = stokeslet_SLP_triangle(P,M,f,slpcache);

    relerr = max(vecnorm(u, 2, 2));
    normal_error = [normal_error; relerr];

    %rotational test field
   
    ez = 0*P; ez(:,3) = 1;
    f = F*cross(ez,P,2);

    u = stokeslet_SLP_triangle(P,M,f,slpcache);
    rot_error = [rot_error;max(vecnorm(u - (1/3)*f,2,2))/((1/3)*F)];


    % combined field
    v=(1/3)*f;
    f = f + 3*geo.v_normal;
    slpcache = stokeslet_SLP_triangle_setup(M);
    u = stokeslet_SLP_triangle(P,M,f,slpcache);
    comb_error = [comb_error;max(vecnorm(u - v,2,2))/((1/3)*F)];

    % linear field
    A = [1000 3 1; 2 2 2; 2 .001 -1];
    f = P * A.';

    v = (1/6)*P*(A-A').' + (1/5)*(P *((A+A')/2).' - P.*(trace(A)/3));

    slpcache = stokeslet_SLP_triangle_setup(M);
    u = stokeslet_SLP_triangle(P,M,f,slpcache);
    lin_error = [lin_error;max(vecnorm(u - v,2,2))/max(vecnorm(v,2,2))];

    % normal field on a cube
    
    [P,M] = subdivided_weirdshape(n); geo=Geometry(M,P);

    f = geo.v_normal;
    slpcache = stokeslet_SLP_triangle_setup(M);
    u = stokeslet_SLP_triangle(P,M,f,slpcache);

    relerr = max(vecnorm(u, 2, 2));
    cubenorm_error = [cubenorm_error; relerr];


    
end


p_trans = polyfit(log(dx),log(trans_error),1);
p_normal = polyfit(log(dx),log(normal_error),1);
p_rot = polyfit(log(dx),log(rot_error),1);
p_comb = polyfit(log(dx),log(comb_error),1);
p_lin = polyfit(log(dx),log(lin_error),1);
p_cube = polyfit(log(dx),log(cubenorm_error),1);

% Second analytical test: rigid body normalational traction gives rigid body
% motion 

figure;
trisurf(M,P(:,1),P(:,2),P(:,3)); hold on;
quiver3(P(:,1),P(:,2),P(:,3),f(:,1),f(:,2),f(:,3))
quiver3(P(:,1),P(:,2),P(:,3),u(:,1),u(:,2),u(:,3))

figure;
loglog(dx,trans_error,LineWidth=2,DisplayName="uniform force");hold on;
loglog(dx,normal_error,LineWidth=2,DisplayName="normal force");
loglog(dx,rot_error,LineWidth=2,DisplayName="rotational force");
loglog(dx,comb_error,LineWidth=2,DisplayName="rot+norm force");
loglog(dx,lin_error,LineWidth=2,DisplayName="general force f = Ax");
loglog(dx,cubenorm_error,LineWidth=2,DisplayName="normal force on cube");
xlabel("dx",FontSize=17); ylabel("\epsilon",FontSize=25)
set(gca,"YScale","log")
title("Analytical test cases")


legend;

