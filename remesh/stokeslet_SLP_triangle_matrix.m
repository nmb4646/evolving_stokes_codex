function S = stokeslet_SLP_triangle_matrix(P, F, cache)
% Explicit matrix for stokeslet_SLP_triangle.
%
% P     : n_v x 3 vertex positions
% F     : n_f x 3 triangle connectivity
% cache : output of stokeslet_SLP_triangle_setup(F)
%
% S maps f(:) to u(:), where f and u are n_v x 3 arrays.

if nargin < 3 || isempty(cache)
    cache = stokeslet_SLP_triangle_setup(F);
end

coef = 1 / (8 * pi);
n_v = size(P, 1);
n_f = size(F, 1);

v1 = P(F(:, 1), :);
v2 = P(F(:, 2), :);
v3 = P(F(:, 3), :);

face_cross = cross(v2 - v1, v3 - v1, 2);
face_area = 0.5 * vecnorm(face_cross, 2, 2);

bary = cache.face_bary;
src_pos = bary(1) * v1 + bary(2) * v2 + bary(3) * v3;

max_entries = 27 * n_v * n_f + 81 * n_f;
rows = zeros(max_entries, 1);
cols = zeros(max_entries, 1);
vals = zeros(max_entries, 1);
cursor = 0;

for first = 1:cache.chunk_size:n_v
    last = min(first + cache.chunk_size - 1, n_v);
    target_ids = (first:last).';
    target = P(target_ids, :);
    n_t = numel(target_ids);

    r = reshape(target, [], 1, 3) - reshape(src_pos, 1, [], 3);
    r2 = sum(r .^ 2, 3);
    rnorm = sqrt(r2);
    inv_r = 1 ./ rnorm;
    inv_r3 = 1 ./ (r2 .* rnorm);

    target_grid = repmat(target_ids, 1, n_f);
    face_weight = reshape(coef * face_area, 1, []);

    for target_dim = 1:3
        row_grid = target_grid + (target_dim - 1) * n_v;
        r_target_dim = r(:, :, target_dim);
        for source_dim = 1:3
            G_component = (target_dim == source_dim) * inv_r ...
                + r_target_dim .* r(:, :, source_dim) .* inv_r3;
            weighted_component = face_weight .* G_component;
            for local_idx = 1:3
                col_grid = repmat(F(:, local_idx).' + (source_dim - 1) * n_v, n_t, 1);
                entry_vals = bary(local_idx) * weighted_component;
                next = cursor + n_t * n_f;
                rows(cursor + 1:next) = row_grid(:);
                cols(cursor + 1:next) = col_grid(:);
                vals(cursor + 1:next) = entry_vals(:);
                cursor = next;
            end
        end
    end
end

for target_idx = 1:n_v
    incident = cache.incident{target_idx};
    if isempty(incident)
        continue
    end

    x = P(target_idx, :);
    for j = 1:size(incident, 1)
        face_idx = incident(j, 1);
        local_idx = incident(j, 2);
        verts = P(F(face_idx, :), :);
        source_vertices = F(face_idx, :);

        centroid_blocks = face_centroid_rule_blocks(x, verts);
        self_blocks = face_self_duffy_blocks(x, verts, local_idx, cache);
        correction_blocks = self_blocks - centroid_blocks;

        for local_source = 1:3
            source_idx = source_vertices(local_source);
            [block_rows, block_cols, block_vals] = block_entries( ...
                target_idx, source_idx, coef * correction_blocks(:, :, local_source), n_v);
            next = cursor + numel(block_vals);
            rows(cursor + 1:next) = block_rows;
            cols(cursor + 1:next) = block_cols;
            vals(cursor + 1:next) = block_vals;
            cursor = next;
        end
    end
end

S = sparse(rows(1:cursor), cols(1:cursor), vals(1:cursor), 3 * n_v, 3 * n_v);
end

function G = stokeslet_block(x, y)
r = x - y;
r2 = dot(r, r);
rnorm = sqrt(r2);
G = eye(3) / rnorm + (r(:) * r(:).') / (r2 * rnorm);
end

function blocks = face_centroid_rule_blocks(x, verts)
e1 = verts(2, :) - verts(1, :);
e2 = verts(3, :) - verts(1, :);
area = 0.5 * norm(cross(e1, e2));
bary = [1/3, 1/3, 1/3];
y = bary(1) * verts(1, :) + bary(2) * verts(2, :) + bary(3) * verts(3, :);
G = area * stokeslet_block(x, y);

blocks = zeros(3, 3, 3);
for local_idx = 1:3
    blocks(:, :, local_idx) = bary(local_idx) * G;
end
end

function blocks = face_self_duffy_blocks(x, verts, local_idx, cache)
order = [local_idx, mod(local_idx, 3) + 1, mod(local_idx + 1, 3) + 1];
verts_local = verts(order, :);

a = verts_local(1, :);
e1 = verts_local(2, :) - a;
e2 = verts_local(3, :) - a;
jac_scale = norm(cross(e1, e2));

blocks_local = zeros(3, 3, 3);
for q = 1:numel(cache.duffy_w)
    u_param = cache.duffy_u(q);
    v_param = cache.duffy_v(q);
    bary = [1 - u_param, u_param * (1 - v_param), u_param * v_param];
    y = bary(1) * verts_local(1, :) + bary(2) * verts_local(2, :) + bary(3) * verts_local(3, :);
    jacobian = jac_scale * u_param;
    G = cache.duffy_w(q) * jacobian * stokeslet_block(x, y);

    for local_idx = 1:3
        blocks_local(:, :, local_idx) = blocks_local(:, :, local_idx) + bary(local_idx) * G;
    end
end

blocks = zeros(3, 3, 3);
for local_idx = 1:3
    blocks(:, :, order(local_idx)) = blocks_local(:, :, local_idx);
end
end

function [rows, cols, vals] = block_entries(target_idx, source_idx, block, n_v)
rows = zeros(9, 1);
cols = zeros(9, 1);
vals = zeros(9, 1);
cursor = 0;
for target_dim = 1:3
    row = target_idx + (target_dim - 1) * n_v;
    for source_dim = 1:3
        cursor = cursor + 1;
        rows(cursor) = row;
        cols(cursor) = source_idx + (source_dim - 1) * n_v;
        vals(cursor) = block(target_dim, source_dim);
    end
end
end
