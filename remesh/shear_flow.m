function [u_out] = shear_flow(P_in,gamy)
u_out = zeros(size(P_in));   % Nx3
u_out(:,1) = gamy * P_in(:,3);
end