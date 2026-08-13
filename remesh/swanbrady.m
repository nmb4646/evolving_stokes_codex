function C = swanbrady(a, b, delta)
%SWANBRADY One-stresslet-reflection orbital/spin frequency coefficient.
%
%   C = SWANBRADY(a, b, delta) returns the dimensionless coefficient C in
%
%       omega_o = C .* Omega,
%
%   for two identical spherical particles rotating with equal angular
%   velocity normal to a no-slip wall. The geometry is
%
%       d = 2*a + delta,      h = a + b,
%
%   where a is the particle radius, b is the particle-wall surface gap,
%   and delta is the particle-particle surface gap. DELTA may be a row or
%   column vector; C has the same size and orientation as DELTA. All three
%   length inputs must use the same units.
%
%   This is the first stresslet-reflection reduction of the wall-bounded
%   Swan-Brady mobility, including the corresponding torque-to-spin
%   correction. It does not include particle-particle or particle-wall
%   lubrication corrections and should not be used in the small-gap limit.
%
%   Example:
%       a = 0.5;                 % micrometers
%       b = 0.05;                % micrometers
%       delta = [0.2 0.3 0.4];   % micrometers
%       C = swanbrady(a, b, delta)
%
%   The middle entry is approximately 0.0690958749.

    validateattributes(a, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'positive'}, mfilename, 'a', 1);
    validateattributes(b, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'nonnegative'}, mfilename, 'b', 2);
    validateattributes(delta, {'numeric'}, ...
        {'real', 'finite', 'vector', 'nonnegative'}, mfilename, 'delta', 3);

    % Dimensionless center separation and center height.
    D = 2 + delta ./ a;
    H = 1 + b ./ a;

    % Distance from one sphere center to the image of the other.
    R = sqrt(D.^2 + 4 .* H.^2);
    xi = D ./ R;
    z = 2 .* H ./ R;

    % Direct torque-to-translation and torque-to-rotation mobilities,
    % normalized by 6*pi*eta*a^2 and 6*pi*eta*a^3, respectively.
    A = (3/4) .* (D.^(-2) - D ./ R.^3);

    B = 3/4 ...
        - 3 ./ (32 .* H.^3) ...
        - 3 ./ (8 .* D.^3) ...
        + (3/8) .* (1 - 3 .* z.^2) ./ R.^3;

    % Torque-to-stresslet couplings for the xy and yz stresslet modes.
    p = -(9 .* sqrt(2) ./ 8) .* (D.^(-3) - D.^2 ./ R.^5);
    q = -(9 .* sqrt(2) ./ 4) .* D .* H ./ R.^5;

    % Self wall-reflection coupling from the yz stresslet mode back to
    % tangential translation.
    r = (3 .* sqrt(2) ./ 160) .* ...
        (15 ./ H.^2 - 12 ./ H.^4 + 5 ./ H.^6);

    % Pair wall-reflection terms returning the induced stresslets to the
    % tangential velocity.
    Phi = (15/4) .* z.^2 ./ R.^2 ...
        + (4 - 20 .* z.^2) ./ R.^4 ...
        - 5 .* (1 - 7 .* z.^2) ./ R.^6;

    s = -(6 .* sqrt(2) ./ 5) ./ D.^4 ...
        + (3 .* sqrt(2) ./ 10) .* xi .* Phi;

    t = (3 .* sqrt(2) ./ 10) .* z .* ...
        (4 ./ R.^4 - 10 ./ R.^6 + Phi);

    % One induced-stresslet reflection. The factor 10/9 follows from
    % (M_ES,0)^(-1) = (20*pi*eta*a^3/3) I in the orthonormal stresslet basis.
    A1 = A + (10/9) .* (p .* s + q .* (r - t));
    B1 = B - (10/9) .* (p.^2 + q.^2);

    C = (2 ./ D) .* (A1 ./ B1);
end
