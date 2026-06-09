function A = stokeslet_SLP_galerkin_matrix(P, F, options)
% P1 Galerkin matrix for the Stokes single-layer potential.
%
% P       : n_v x 3 vertex positions
% F       : n_f x 3 triangle connectivity
% options : optional struct
%           .rule       "dunavant7" (default) or "centroid"
%           .tol        distance tolerance for coincident quadrature points
%           .symmetrize true by default
%
% A maps nodal tractions to weak/tested velocities:
%
%   A_ij = int_Gamma int_Gamma phi_i(x) K(x,y) phi_j(y) dA_y dA_x
%
% with 3 x 3 vector blocks and component-major vector ordering, matching
% stokeslet_SLP_triangle_matrix.

if nargin < 3 || isempty(options)
    options = struct();
end

if ~isfield(options, "rule") || isempty(options.rule)
    options.rule = "dunavant7";
end
if ~isfield(options, "tol") || isempty(options.tol)
    options.tol = 1e-14;
end
if ~isfield(options, "symmetrize") || isempty(options.symmetrize)
    options.symmetrize = true;
end

coef = 1 / (8 * pi);
n_v = size(P, 1);
n_f = size(F, 1);

[bary, q_weights] = triangle_rule(options.rule);
n_q = size(bary, 1);

v1 = P(F(:, 1), :);
v2 = P(F(:, 2), :);
v3 = P(F(:, 3), :);

face_cross = cross(v2 - v1, v3 - v1, 2);
face_area = 0.5 * vecnorm(face_cross, 2, 2);

q_pos = zeros(n_f, n_q, 3);
q_weight = zeros(n_f, n_q);
for q = 1:n_q
    b = bary(q, :);
    q_pos(:, q, :) = b(1) * v1 + b(2) * v2 + b(3) * v3;
    q_weight(:, q) = q_weights(q) * face_area;
end

max_entries = 81 * n_f * n_f;
rows = zeros(max_entries, 1);
cols = zeros(max_entries, 1);
vals = zeros(max_entries, 1);
cursor = 0;

for target_face = 1:n_f
    target_vertices = F(target_face, :);

    for source_face = target_face:n_f
        source_vertices = F(source_face, :);
        local_blocks = face_pair_blocks( ...
            reshape(q_pos(target_face, :, :), n_q, 3), q_weight(target_face, :), ...
            reshape(q_pos(source_face, :, :), n_q, 3), q_weight(source_face, :), ...
            bary, coef, options.tol);

        [block_rows, block_cols, block_vals] = local_block_entries( ...
            target_vertices, source_vertices, local_blocks, n_v);
        next = cursor + numel(block_vals);
        rows(cursor + 1:next) = block_rows;
        cols(cursor + 1:next) = block_cols;
        vals(cursor + 1:next) = block_vals;
        cursor = next;

        if source_face ~= target_face
            [block_rows, block_cols, block_vals] = local_block_entries_transpose( ...
                target_vertices, source_vertices, local_blocks, n_v);
            next = cursor + numel(block_vals);
            rows(cursor + 1:next) = block_rows;
            cols(cursor + 1:next) = block_cols;
            vals(cursor + 1:next) = block_vals;
            cursor = next;
        end
    end
end

A = sparse(rows(1:cursor), cols(1:cursor), vals(1:cursor), 3 * n_v, 3 * n_v);

if options.symmetrize
    A = 0.5 * (A + A.');
end
end

function blocks = face_pair_blocks(x_q, x_w, y_q, y_w, bary, coef, tol)
blocks = zeros(3, 3, 3, 3);
n_q = size(bary, 1);

for qx = 1:n_q
    x = x_q(qx, :);
    phix = bary(qx, :);

    for qy = 1:n_q
        y = y_q(qy, :);
        phiy = bary(qy, :);
        r = x - y;
        r2 = dot(r, r);
        rnorm = sqrt(r2);
        if rnorm < tol
            continue
        end

        weight = coef * x_w(qx) * y_w(qy);
        G = weight * (eye(3) / rnorm + (r(:) * r(:).') / (r2 * rnorm));

        for target_local = 1:3
            for source_local = 1:3
                blocks(:, :, target_local, source_local) = ...
                    blocks(:, :, target_local, source_local) ...
                    + phix(target_local) * phiy(source_local) * G;
            end
        end
    end
end
end

function [rows, cols, vals] = local_block_entries(target_vertices, source_vertices, blocks, n_v)
rows = zeros(81, 1);
cols = zeros(81, 1);
vals = zeros(81, 1);
cursor = 0;

for target_local = 1:3
    target_idx = target_vertices(target_local);
    for source_local = 1:3
        source_idx = source_vertices(source_local);
        block = blocks(:, :, target_local, source_local);
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
end
end

function [rows, cols, vals] = local_block_entries_transpose(target_vertices, source_vertices, blocks, n_v)
rows = zeros(81, 1);
cols = zeros(81, 1);
vals = zeros(81, 1);
cursor = 0;

for target_local = 1:3
    target_idx = target_vertices(target_local);
    for source_local = 1:3
        source_idx = source_vertices(source_local);
        block = blocks(:, :, target_local, source_local).';
        for source_dim = 1:3
            row = source_idx + (source_dim - 1) * n_v;
            for target_dim = 1:3
                cursor = cursor + 1;
                rows(cursor) = row;
                cols(cursor) = target_idx + (target_dim - 1) * n_v;
                vals(cursor) = block(source_dim, target_dim);
            end
        end
    end
end
end

function [bary, weights] = triangle_rule(rule)
switch string(rule)
    case "centroid"
        bary = [1/3, 1/3, 1/3];
        weights = 1;
    case "dunavant7"
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
    otherwise
        error("Unsupported triangle quadrature rule '%s'.", string(rule));
end
end
