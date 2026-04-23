function u_out = stokeslet_SLP(P, v_area, f)
% P      : Nx3 points
% v_area : Nx1 areas
% f      : Nx3 traction
%
% u      : Nx3 velocity (single-layer potential)

N = size(P,1);
u_out = zeros(N,3);

coef = 1 / (8*pi);

for i = 1:N
    % vector from sources to target i
    r = P(i,:) - P;           % Nx3
    r2 = sum(r.^2, 2);        % Nx1
    rnorm = sqrt(r2);         % Nx1

    % avoid singularity (ignore i=j)
    mask = (rnorm > 0);

    r = r(mask,:);
    rnorm = rnorm(mask);
    r2 = r2(mask);
    
    fj = f(mask,:);
    Aj = v_area(mask);

    % scalar factors
    inv_r = 1 ./ rnorm;
    inv_r3 = 1 ./ (rnorm.^3);

    % term1: (I / r) * f
    term1 = fj .* inv_r;   % Nx3

    % term2: (r r^T / r^3) f = r * (dot(r,f)/r^3)
    dot_rf = sum(r .* fj, 2);              % Nx1
    term2 = r .* (dot_rf .* inv_r3);       % Nx3

    % combine and weight by area
    contrib = (term1 + term2) .* Aj;

    % sum contributions
    u_out(i,:) = coef * sum(contrib, 1);
end
end