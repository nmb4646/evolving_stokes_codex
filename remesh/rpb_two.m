function C = rpb_two(a, b, delta)
%RPB_TWO Rotne-Prager-Blake model for two rigid two-sphere peanuts.
%
%   C = RPB_TWO(a,b,delta) returns the dimensionless coefficient C in
%
%       omega_o = C .* Omega,
%
%   for two identical rigid peanuts near a no-slip plane wall.  Each
%   peanut consists of two equal, tangent spheres of radius a aligned
%   normal to the wall.  The two peanut axes are separated horizontally by
%
%       d = 2*a + delta.
%
%   The lower sphere has wall gap b, so the lower and upper sphere-center
%   heights are a+b and 3*a+b.  Both peanuts are prescribed the same spin
%   Omega about +z and are otherwise force-free in the tangential (y)
%   direction.  DELTA may be a row or column vector; C has the same shape.
%
%   Hydrodynamics:
%     * generalized Rotne-Prager-Blake force/torque mobility for all four
%       spheres (Uy/Omega_z coupled to Fy/Lz), including the no-slip wall;
%     * rigid-body constraints applied to the two lobes of each peanut;
%     * an excess sphere-wall lubrication resistance on each LOWER lobe;
%     * deliberately NO particle-particle lubrication correction.
%
%   The RPB block is obtained by applying the translational Faxen operator
%   to Blake's half-space Green function.  The lower-lobe wall correction
%   replaces the isolated-sphere RPB parallel-wall resistance by the
%   all-gap fit of Dunstan et al. (2012), Appendix A, Eq. (A1).
%
%   References:
%     J. W. Swan and J. F. Brady, Phys. Fluids 19, 113306 (2007).
%     J. Dunstan et al., Phys. Biol. 9, 066003 (2012), Appendix A.
%
%   Example:
%       a = 0.5;
%       b = 0.05;
%       delta = [0.2 0.3 0.4];
%       C = rpb_two(a,b,delta)

    validateattributes(a, {'numeric'}, ...
        {'real','finite','scalar','positive'}, mfilename, 'a', 1);
    validateattributes(b, {'numeric'}, ...
        {'real','finite','scalar','nonnegative'}, mfilename, 'b', 2);
    validateattributes(delta, {'numeric'}, ...
        {'real','finite','vector','nonnegative'}, mfilename, 'delta', 3);

    % Scale all lengths by a and set eta=1.  With the corresponding force,
    % torque, translation and rotation scales, C is independent of units.
    Dvec = 2 + delta(:) ./ a;
    Hb = 1 + b ./ a;
    Ht = 3 + b ./ a;

    % Exact contact with the stationary wall gives divergent tangential
    % lubrication resistance, hence zero force-free orbital translation.
    if b == 0
        C = zeros(size(delta));
        return
    end

    % Replace the lower-sphere RPB self resistance by the fitted all-gap
    % single-wall resistance.  This is an EXCESS resistance: the wall part
    % already contained in RPB is subtracted to avoid double counting.
    Tself = 1 - (1/16) * (9/Hb - 2/Hb^3 + 1/Hb^5);
    zetaFar = 1 / Tself;

    logGap = log1p(-1/Hb);
    zetaWall = -( ...
        (2.5295 - 1.9963*Hb) * logGap ...
        - 2.9963 ...
        + 0.9689/Hb ...
        + 0.5993/Hb^2 ...
        + 0.4691/Hb^3);
    deltaZeta = zetaWall - zetaFar;

    % Maps body velocities [U1 Omega1 U2 Omega2] to the eight lobe DOFs
    % [Uy_1b Oz_1b Uy_1t Oz_1t Uy_2b Oz_2b Uy_2t Oz_2t].
    K = [1 0 0 0;
         0 1 0 0;
         1 0 0 0;
         0 1 0 0;
         0 0 1 0;
         0 0 0 1;
         0 0 1 0;
         0 0 0 1];

    Cvec = zeros(size(Dvec));

    for k = 1:numel(Dvec)
        D = Dvec(k);

        % Lobe order: peanut 1 bottom/top, peanut 2 bottom/top.
        x = [-D/2; -D/2; D/2; D/2];
        z = [Hb; Ht; Hb; Ht];

        M = zeros(8,8);

        for ii = 1:4
            rows = (2*ii-1):(2*ii);
            for jj = 1:4
                cols = (2*jj-1):(2*jj);
                q = x(ii) - x(jj);

                if ii == jj
                    % Isolated sphere plus its RPB wall reflection.
                    block = diag([1/(6*pi), 1/(8*pi)]) ...
                        + rpb_wall(q, z(ii), z(jj));
                else
                    % Free-space RPY pair plus the Blake-wall reflection.
                    block = rpb_direct(q, z(ii), z(jj)) ...
                        + rpb_wall(q, z(ii), z(jj));
                end

                M(rows,cols) = block;
            end
        end

        % Reciprocity makes M symmetric analytically; remove only floating
        % point antisymmetry before converting mobility to resistance.
        M = 0.5 * (M + M.');
        Rbead = M \ eye(8);

        % Wall lubrication only on the lower lobe of each peanut.  The
        % resistance scale for a tangential sphere is 6*pi*eta*a.
        Rbead(1,1) = Rbead(1,1) + 6*pi*deltaZeta;
        Rbead(5,5) = Rbead(5,5) + 6*pi*deltaZeta;

        % Project bead resistance to the two rigid bodies.
        Rbody = K.' * Rbead * K;

        % Each peanut is force-free in y.  Prescribe Omega1=Omega2=1;
        % linearity then gives U/Omega directly.
        t = [1 3];
        r = [2 4];
        U = -Rbody(t,t) \ (Rbody(t,r) * [1;1]);

        % For centers at x=+-D/2, omega_o=(U2-U1)/D.  Lengths and
        % translations are both scaled by a, so this is C=omega_o/Omega.
        Cvec(k) = (U(2) - U(1)) / D;
    end

    C = reshape(Cvec, size(delta));
end


function block = rpb_direct(q, z, w)
%RPB_DIRECT Free-space generalized RPY block for r >= 2a, with a=eta=1.
% Rows are [Uy, Omega_z], columns are [Fy, Lz].

    dz = z - w;
    r2 = q^2 + dz^2;
    r = sqrt(r2);

    block = [ ...
        (1/r + 2/(3*r^3))/(8*pi),              q/(8*pi*r^3); ...
                       -q/(8*pi*r^3), (3*dz^2/r2 - 1)/(16*pi*r^3)];
end


function block = rpb_wall(q, z, w)
%RPB_WALL Wall-reflected generalized RPB block, with a=eta=1.
% q is receiver x minus source x; z and w are receiver/source heights.

    p = z + w;
    r2 = q^2 + p^2;

    % Faxen-corrected yy component of the Blake image system.  Written in
    % this compact form, it is algebraically identical to applying the
    % receiver and source translational Faxen operators to Blake's tensor.
    n = -3*r2^3 ...
        - 6*w*z*r2^2 ...
        + 2*(2*p^2 - q^2)*r2 ...
        + 2*(q^2 - 4*p^2);

    block = [ ...
        n/(24*pi*r2^(7/2)),                  -q/(8*pi*r2^(3/2)); ...
          q/(8*pi*r2^(3/2)), (q^2 - 2*p^2)/(16*pi*r2^(5/2))];
end
