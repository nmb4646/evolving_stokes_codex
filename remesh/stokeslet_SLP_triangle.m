function u_out = stokeslet_SLP_triangle(P, F, f, cache)
% Triangle-quadrature single-layer Stokes potential evaluated at vertices.
%
% P     : n_v x 3 vertex positions
% F     : n_f x 3 triangle connectivity
% f     : n_v x 3 traction values stored at vertices
% cache : output of stokeslet_SLP_triangle_setup(F)

if nargin < 4 || isempty(cache)
    cache = stokeslet_SLP_triangle_setup(F);
end

coef = 1 / (8 * pi);
n_v = size(P, 1);
n_f = size(F, 1);

v1 = P(F(:, 1), :);
v2 = P(F(:, 2), :);
v3 = P(F(:, 3), :);
f1 = f(F(:, 1), :);
f2 = f(F(:, 2), :);
f3 = f(F(:, 3), :);

face_cross = cross(v2 - v1, v3 - v1, 2);
face_area = 0.5 * vecnorm(face_cross, 2, 2);

bary = cache.face_bary;
src_pos = bary(1) * v1 + bary(2) * v2 + bary(3) * v3;
f_q = bary(1) * f1 + bary(2) * f2 + bary(3) * f3;
src_force = face_area .* f_q;

u_out = zeros(n_v, 3);
for first = 1:cache.chunk_size:n_v
    last = min(first + cache.chunk_size - 1, n_v);
    target = P(first:last, :);
    u_out(first:last, :) = coef * accumulate_sources(target, src_pos, src_force);
end

for i = 1:n_v
    incident = cache.incident{i};
    if isempty(incident)
        continue
    end

    x = P(i, :);
    correction = zeros(1, 3);
    for j = 1:size(incident, 1)
        face_idx = incident(j, 1);
        local_idx = incident(j, 2);
        verts = P(F(face_idx, :), :);
        traction = f(F(face_idx, :), :);
        correction = correction ...
            + face_self_duffy(x, verts, traction, local_idx, cache) ...
            - face_centroid_rule(x, verts, traction);
    end
    u_out(i, :) = u_out(i, :) + coef * correction;
end
end

function u = accumulate_sources(target, src_pos, src_force)
r = reshape(target, [], 1, 3) - reshape(src_pos, 1, [], 3);
r2 = sum(r .^ 2, 3);
rnorm = sqrt(r2);
inv_r = 1 ./ rnorm;
inv_r3 = 1 ./ (r2 .* rnorm);

dot_rf = r(:, :, 1) .* reshape(src_force(:, 1), 1, []) ...
    + r(:, :, 2) .* reshape(src_force(:, 2), 1, []) ...
    + r(:, :, 3) .* reshape(src_force(:, 3), 1, []);

u = zeros(size(target, 1), 3);
for d = 1:3
    fd = reshape(src_force(:, d), 1, []);
    u(:, d) = sum(inv_r .* fd + r(:, :, d) .* dot_rf .* inv_r3, 2);
end
end

function u = face_centroid_rule(x, verts, traction)
e1 = verts(2, :) - verts(1, :);
e2 = verts(3, :) - verts(1, :);
area = 0.5 * norm(cross(e1, e2));
bary = [1/3, 1/3, 1/3];
y = bary(1) * verts(1, :) + bary(2) * verts(2, :) + bary(3) * verts(3, :);
fq = bary(1) * traction(1, :) + bary(2) * traction(2, :) + bary(3) * traction(3, :);
r = x - y;
r2 = dot(r, r);
rnorm = sqrt(r2);
inv_r = 1 / rnorm;
inv_r3 = 1 / (r2 * rnorm);
u = area * (fq * inv_r + r * (dot(r, fq) * inv_r3));
end

function u = face_self_duffy(x, verts, traction, local_idx, cache)
order = [local_idx, mod(local_idx, 3) + 1, mod(local_idx + 1, 3) + 1];
verts_local = verts(order, :);
traction_local = traction(order, :);

a = verts_local(1, :);
e1 = verts_local(2, :) - a;
e2 = verts_local(3, :) - a;
jac_scale = norm(cross(e1, e2));

u = zeros(1, 3);
for q = 1:numel(cache.duffy_w)
    u_param = cache.duffy_u(q);
    v_param = cache.duffy_v(q);
    bary = [1 - u_param, u_param * (1 - v_param), u_param * v_param];
    y = bary(1) * verts_local(1, :) + bary(2) * verts_local(2, :) + bary(3) * verts_local(3, :);
    fq = bary(1) * traction_local(1, :) + bary(2) * traction_local(2, :) + bary(3) * traction_local(3, :);

    r = x - y;
    r2 = dot(r, r);
    rnorm = sqrt(r2);
    inv_r = 1 / rnorm;
    inv_r3 = 1 / (r2 * rnorm);
    jacobian = jac_scale * u_param;

    u = u + cache.duffy_w(q) * jacobian * (fq * inv_r + r * (dot(r, fq) * inv_r3));
end
end
