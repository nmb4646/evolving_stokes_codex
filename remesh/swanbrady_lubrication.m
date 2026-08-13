function C = swanbrady(a, b, delta)
%SWANBRADY Orbital/spin coefficient with sphere-wall lubrication.
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
%   The far-field mobility is the first stresslet-reflection reduction of
%   the wall-bounded Swan-Brady mobility. A sphere-wall lubrication
%   correction is then added at resistance level as
%
%       R = R_far + R_wall - R_wall_far.
%
%   R_wall is the all-gap fit for translation parallel to a wall given by
%   Dunstan et al. (2012), Appendix A, Eq. (A1). It recovers the
%   Goldman-Cox-Brenner near-contact resistance
%
%       R_wall/(6*pi*eta*a) ~ -(8/15)*log(b/a) + 0.9588.
%
%   This is the leading wall-lubrication correction in the symmetry-reduced
%   orbital mode, not a full inversion of every stresslet-level wall
%   resistance block. There is deliberately NO particle-particle
%   lubrication correction. Rotation is normal to the wall, so its
%   single-wall resistance does not enter C when Omega is prescribed; it
%   changes the required torque, not the force-free ratio omega_o/Omega.
%
%   References:
%     J. W. Swan and J. F. Brady, Phys. Fluids 22, 103301 (2010), Eq. (43).
%     J. Dunstan et al., Phys. Biol. 9, 066003 (2012), Appendix A, Eq. (A1).
%     A. J. Goldman et al., Chem. Eng. Sci. 22, 637-651 (1967).
%
%   Example:
%       a = 0.5;                 % micrometers
%       b = 0.05;                % micrometers
%       delta = [0.2 0.3 0.4];   % micrometers
%       C = swanbrady(a, b, delta)
%
%   The middle entry is approximately 0.0663964341.

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

    % Antisymmetric tangential-translation mobility. This mode produces
    % the mutual orbit: the two spheres translate in opposite directions.
    % Tself is the single-sphere Swan-Brady parallel mobility; Tpair is the
    % direct pair mobility (unbounded RPY plus its no-slip-wall image).
    Tself = 1 - (1/16) .* ...
        (9 ./ H - 2 ./ H.^3 + 1 ./ H.^5);

    Tpair = 3 ./ (4 .* D) + 1 ./ (2 .* D.^3) ...
        - (1/4) .* ( ...
            3 .* (1 + z.^2 ./ 2) ./ R ...
            + 2 .* (1 - 3 .* z.^2) ./ R.^3 ...
            - 2 .* (1 - 5 .* z.^2) ./ R.^5);

    T0 = Tself - Tpair;
    T1 = T0 - (10/9) .* (s.^2 + (r - t).^2);

    % At exact wall contact, tangential resistance diverges logarithmically,
    % so force-free orbital translation vanishes for finite prescribed spin.
    if b == 0
        C = zeros(size(delta));
        return
    end

    % Positive magnitude of the all-gap sphere-wall resistance for
    % translation parallel to the wall, normalized by 6*pi*eta*a.
    % log1p is used to evaluate log((H-1)/H) accurately at large H.
    logGap = log1p(-1 ./ H);
    zetaWall = -( ...
        (2.5295 - 1.9963 .* H) .* logGap ...
        - 2.9963 ...
        + 0.9689 ./ H ...
        + 0.5993 ./ H.^2 ...
        + 0.4691 ./ H.^3);

    % The same one-stresslet reduction for an isolated sphere supplies the
    % wall contribution already present in the far-field model. Subtracting
    % its inverse avoids double-counting that wall resistance.
    T_wall_far = Tself - (10/9) .* r.^2;
    deltaZeta = zetaWall - 1 ./ T_wall_far;

    % Add deltaZeta to the translational entry of the inverse 2-by-2
    % mobility in the antisymmetric-translation/symmetric-rotation mode.
    % No particle-particle resistance or lubrication term is added here.
    mobilityDet = T1 .* B1 - A1.^2;
    C = (2 ./ D) .* A1 ./ ...
        (B1 + deltaZeta .* mobilityDet);
end
