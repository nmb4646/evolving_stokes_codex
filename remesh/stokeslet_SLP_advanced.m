function u_out = stokeslet_SLP_advanced(P, F, f, cache)
% Higher-order triangle-quadrature single-layer Stokes potential.
%
% P     : n_v x 3 vertex positions
% F     : n_f x 3 triangle connectivity
% f     : n_v x 3 traction values stored at vertices
% cache : output of stokeslet_SLP_advanced_setup(F)
%
% Note to future edits:
% This routine is opt-in only. Do not replace calls to
% stokeslet_SLP_triangle with this function unless specifically instructed.

if nargin < 4 || isempty(cache) || ~isstruct(cache) || ~isfield(cache, 'is_advanced')
    cache = stokeslet_SLP_advanced_setup(F);
end

coef = 1 / (8 * pi);
n_v = size(P, 1);
n_f = size(F, 1);
n_q = size(cache.regular_bary, 1);

v1 = P(F(:, 1), :);
v2 = P(F(:, 2), :);
v3 = P(F(:, 3), :);
f1 = f(F(:, 1), :);
f2 = f(F(:, 2), :);
f3 = f(F(:, 3), :);

face_cross = cross(v2 - v1, v3 - v1, 2);
face_area = 0.5 * vecnorm(face_cross, 2, 2);
face_size = sqrt(face_area);
centroid = (v1 + v2 + v3) / 3;

src_pos = zeros(n_f * n_q, 3);
src_force = zeros(n_f * n_q, 3);
for q = 1:n_q
    rows = (q - 1) * n_f + (1:n_f);
    bary = cache.regular_bary(q, :);
    weights = cache.regular_weights(q) * face_area;
    src_pos(rows, :) = bary(1) * v1 + bary(2) * v2 + bary(3) * v3;
    f_q = bary(1) * f1 + bary(2) * f2 + bary(3) * f3;
    src_force(rows, :) = weights .* f_q;
end

u_out = zeros(n_v, 3);
for first = 1:cache.chunk_size:n_v
    last = min(first + cache.chunk_size - 1, n_v);
    target = P(first:last, :);
    u_chunk = accumulate_sources(target, src_pos, src_force);

    r_centroid = reshape(target, [], 1, 3) - reshape(centroid, 1, [], 3);
    dist_centroid = sqrt(sum(r_centroid .^ 2, 3));
    near_mask = dist_centroid < reshape(cache.near_scale * face_size, 1, []);

    for local_idx = 1:size(target, 1)
        global_idx = first + local_idx - 1;
        incident = cache.incident{global_idx};
        if isempty(incident)
            incident_faces = zeros(0, 1);
            incident_local = zeros(0, 1);
        else
            incident_faces = incident(:, 1);
            incident_local = incident(:, 2);
        end

        near_faces = find(near_mask(local_idx, :));
        candidate_faces = unique([near_faces(:); incident_faces(:)])';
        if isempty(candidate_faces)
            continue
        end

        x = P(global_idx, :);
        correction = zeros(1, 3);
        for jj = 1:numel(candidate_faces)
            face_idx = candidate_faces(jj);
            verts = [v1(face_idx, :); v2(face_idx, :); v3(face_idx, :)];
            traction = [f1(face_idx, :); f2(face_idx, :); f3(face_idx, :)];
            regular_term = face_regular_rule(x, verts, traction, cache.regular_bary, cache.regular_weights);

            incident_match = find(incident_faces == face_idx, 1);
            if ~isempty(incident_match)
                refined_term = face_self_duffy(x, verts, traction, incident_local(incident_match), cache);
            else
                panel_dist = point_triangle_distance(x, verts);
                depth = choose_adaptive_depth(panel_dist, face_size(face_idx), cache);
                if depth <= 0
                    continue
                end
                refined_term = panel_adaptive(x, verts, traction, cache.regular_bary, cache.regular_weights, depth);
            end

            correction = correction + refined_term - regular_term;
        end

        u_chunk(local_idx, :) = u_chunk(local_idx, :) + correction;
    end

    u_out(first:last, :) = coef * u_chunk;
end
end

function depth = choose_adaptive_depth(dist, face_size, cache)
scaled_dist = dist / max(face_size, eps);
if scaled_dist >= cache.near_scale
    depth = 0;
elseif scaled_dist >= 1.5
    depth = 1;
elseif scaled_dist >= cache.very_near_scale
    depth = min(2, cache.max_adaptive_depth);
else
    depth = cache.max_adaptive_depth;
end
end

function u = accumulate_sources(target, src_pos, src_force)
r = reshape(target, [], 1, 3) - reshape(src_pos, 1, [], 3);
r2 = sum(r .^ 2, 3);
rnorm = sqrt(r2);

tol = 1e-14;
mask = rnorm < tol;
rnorm(mask) = Inf;
r2(mask) = Inf;

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

function u = face_regular_rule(x, verts, traction, bary, weights)
area = 0.5 * norm(cross(verts(2, :) - verts(1, :), verts(3, :) - verts(1, :)));
u = zeros(1, 3);
for q = 1:size(bary, 1)
    b = bary(q, :);
    y = b(1) * verts(1, :) + b(2) * verts(2, :) + b(3) * verts(3, :);
    fq = b(1) * traction(1, :) + b(2) * traction(2, :) + b(3) * traction(3, :);
    u = u + stokeslet_kernel(x, y, weights(q) * area * fq);
end
end

function u = panel_adaptive(x, verts, traction, bary, weights, depth)
if depth <= 0
    u = face_regular_rule(x, verts, traction, bary, weights);
    return
end

m12 = 0.5 * (verts(1, :) + verts(2, :));
m23 = 0.5 * (verts(2, :) + verts(3, :));
m31 = 0.5 * (verts(3, :) + verts(1, :));
t12 = 0.5 * (traction(1, :) + traction(2, :));
t23 = 0.5 * (traction(2, :) + traction(3, :));
t31 = 0.5 * (traction(3, :) + traction(1, :));

sub_verts = {
    [verts(1, :); m12; m31], [traction(1, :); t12; t31];
    [m12; verts(2, :); m23], [t12; traction(2, :); t23];
    [m31; m23; verts(3, :)], [t31; t23; traction(3, :)];
    [m12; m23; m31], [t12; t23; t31]
};

u = zeros(1, 3);
for k = 1:4
    u = u + panel_adaptive(x, sub_verts{k, 1}, sub_verts{k, 2}, bary, weights, depth - 1);
end
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

function u = stokeslet_kernel(x, y, weighted_force)
r = x - y;
r2 = dot(r, r);
rnorm = sqrt(r2);
if rnorm < 1e-14
    u = zeros(1, 3);
    return
end
inv_r = 1 / rnorm;
inv_r3 = 1 / (r2 * rnorm);
u = weighted_force * inv_r + r * (dot(r, weighted_force) * inv_r3);
end

function dist = point_triangle_distance(x, verts)
a = verts(1, :);
b = verts(2, :);
c = verts(3, :);

ab = b - a;
ac = c - a;
ap = x - a;
d1 = dot(ab, ap);
d2 = dot(ac, ap);
if d1 <= 0 && d2 <= 0
    dist = norm(ap);
    return
end

bp = x - b;
d3 = dot(ab, bp);
d4 = dot(ac, bp);
if d3 >= 0 && d4 <= d3
    dist = norm(bp);
    return
end

vc = d1 * d4 - d3 * d2;
if vc <= 0 && d1 >= 0 && d3 <= 0
    v = d1 / (d1 - d3);
    proj = a + v * ab;
    dist = norm(x - proj);
    return
end

cp = x - c;
d5 = dot(ab, cp);
d6 = dot(ac, cp);
if d6 >= 0 && d5 <= d6
    dist = norm(cp);
    return
end

vb = d5 * d2 - d1 * d6;
if vb <= 0 && d2 >= 0 && d6 <= 0
    w = d2 / (d2 - d6);
    proj = a + w * ac;
    dist = norm(x - proj);
    return
end

va = d3 * d6 - d5 * d4;
if va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0
    edge = c - b;
    w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
    proj = b + w * edge;
    dist = norm(x - proj);
    return
end

normal = cross(ab, ac);
dist = abs(dot(ap, normal)) / norm(normal);
end
