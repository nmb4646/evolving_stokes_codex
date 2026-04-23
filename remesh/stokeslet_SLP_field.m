function u_out = stokeslet_SLP_field(X, P, F, f, options)
% Off-surface Stokes single-layer potential with near-panel refinement.
%
% X : n_t x 3 target points
% P : n_v x 3 membrane vertices
% F : n_f x 3 triangle connectivity
% f : n_v x 3 traction values at vertices
%
% options.near_scale      : threshold in face-size units for refined quadrature
% options.very_near_scale : threshold in face-size units for adaptive subdivision
% options.chunk_size      : target chunk size

if nargin < 5
    options = struct();
end
if ~isfield(options, 'near_scale')
    options.near_scale = 2.5;
end
if ~isfield(options, 'very_near_scale')
    options.very_near_scale = 0.9;
end
if ~isfield(options, 'chunk_size')
    options.chunk_size = 128;
end

coef = 1 / (8 * pi);
n_t = size(X, 1);

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
centroid_force = face_area .* ((f1 + f2 + f3) / 3);

[tri_bary, tri_weight] = dunavant_rule_7();

u_out = zeros(n_t, 3);
for first = 1:options.chunk_size:n_t
    last = min(first + options.chunk_size - 1, n_t);
    target = X(first:last, :);
    u_chunk = accumulate_centroid(target, centroid, centroid_force);

    r = reshape(target, [], 1, 3) - reshape(centroid, 1, [], 3);
    dist = sqrt(sum(r.^2, 3));

    near_mask = dist < reshape(options.near_scale * face_size, 1, []);
    near_rows = any(near_mask, 2);
    near_ids = find(near_rows);
    for ii = 1:numel(near_ids)
        local_idx = near_ids(ii);
        global_idx = first + local_idx - 1;
        face_ids = find(near_mask(local_idx, :));
        x = X(global_idx, :);
        for jj = 1:numel(face_ids)
            face_idx = face_ids(jj);
            verts = [v1(face_idx, :); v2(face_idx, :); v3(face_idx, :)];
            traction = [f1(face_idx, :); f2(face_idx, :); f3(face_idx, :)];
            centroid_term = panel_centroid(x, verts, traction);
            if dist(local_idx, face_idx) < options.very_near_scale * face_size(face_idx)
                refined_term = panel_adaptive(x, verts, traction, tri_bary, tri_weight, 1);
            else
                refined_term = panel_triangle_quadrature(x, verts, traction, tri_bary, tri_weight);
            end
            u_chunk(local_idx, :) = u_chunk(local_idx, :) - centroid_term + refined_term;
        end
    end

    u_out(first:last, :) = coef * u_chunk;
end
end

function u = accumulate_centroid(target, src_pos, src_force)
r = reshape(target, [], 1, 3) - reshape(src_pos, 1, [], 3);
r2 = sum(r.^2, 3);
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

function u = panel_centroid(x, verts, traction)
area = 0.5 * norm(cross(verts(2,:) - verts(1,:), verts(3,:) - verts(1,:)));
y = (verts(1,:) + verts(2,:) + verts(3,:)) / 3;
fq = (traction(1,:) + traction(2,:) + traction(3,:)) / 3;
u = stokeslet_kernel(x, y, area * fq);
end

function u = panel_triangle_quadrature(x, verts, traction, bary, weights)
area = 0.5 * norm(cross(verts(2,:) - verts(1,:), verts(3,:) - verts(1,:)));
u = zeros(1, 3);
for q = 1:size(bary, 1)
    b = bary(q, :);
    y = b(1) * verts(1,:) + b(2) * verts(2,:) + b(3) * verts(3,:);
    fq = b(1) * traction(1,:) + b(2) * traction(2,:) + b(3) * traction(3,:);
    u = u + stokeslet_kernel(x, y, weights(q) * area * fq);
end
end

function u = panel_adaptive(x, verts, traction, bary, weights, depth)
if depth <= 0
    u = panel_triangle_quadrature(x, verts, traction, bary, weights);
    return
end

m12 = 0.5 * (verts(1,:) + verts(2,:));
m23 = 0.5 * (verts(2,:) + verts(3,:));
m31 = 0.5 * (verts(3,:) + verts(1,:));
t12 = 0.5 * (traction(1,:) + traction(2,:));
t23 = 0.5 * (traction(2,:) + traction(3,:));
t31 = 0.5 * (traction(3,:) + traction(1,:));

sub_verts = {
    [verts(1,:); m12; m31], [traction(1,:); t12; t31];
    [m12; verts(2,:); m23], [t12; traction(2,:); t23];
    [m31; m23; verts(3,:)], [t31; t23; traction(3,:)];
    [m12; m23; m31], [t12; t23; t31]
};

u = zeros(1, 3);
for k = 1:4
    u = u + panel_adaptive(x, sub_verts{k,1}, sub_verts{k,2}, bary, weights, depth - 1);
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

function [bary, weights] = dunavant_rule_7()
bary = [
    1/3, 1/3, 1/3;
    0.059715871789770, 0.470142064105115, 0.470142064105115;
    0.470142064105115, 0.059715871789770, 0.470142064105115;
    0.470142064105115, 0.470142064105115, 0.059715871789770;
    0.797426985353087, 0.101286507323456, 0.101286507323456;
    0.101286507323456, 0.797426985353087, 0.101286507323456;
    0.101286507323456, 0.101286507323456, 0.797426985353087
];
weights = [
    0.225000000000000;
    0.132394152788506;
    0.132394152788506;
    0.132394152788506;
    0.125939180544827;
    0.125939180544827;
    0.125939180544827
];
end
