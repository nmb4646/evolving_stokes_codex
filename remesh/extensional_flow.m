function U = extensional_flow(P, gamy)
%EXTENSIONAL_FLOW  Incompressible uniaxial extensional flow.
%
%   U = extensional_flow(P, gamy)
%
%   P    : N x 3 point matrix, columns are [x y z]
%   gamy : nondimensional extension rate
%   U    : N x 3 velocity field
%
%   Flow:
%       u = gamy * x
%       v = -0.5 * gamy * y
%       w = -0.5 * gamy * z
%
%   The extension axis is x, so this is appropriate for a dumbbell
%   initially aligned with the x-axis.

    if size(P,2) ~= 3
        error('P must be N x 3.');
    end

    U = zeros(size(P));
    U(:,1) =  gamy      * P(:,1);
    U(:,2) = -0.5*gamy * P(:,2);
    U(:,3) = -0.5*gamy * P(:,3);
end