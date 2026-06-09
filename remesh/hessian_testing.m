clc; clear; close all;

p.k = 1000;
p.dt = 0.1;
p.Sd = 1;
p.Da = 1;

[P, M] = subdivided_sphere(7);
geo = Geometry(M, P);

n_v = geo.mesh.n_v;
I3 = speye(3 * n_v);

mass0 = spdiags(geo.v_area, 0, n_v, n_v);
M3 = blkdiag(mass0, mass0, mass0);

mass0_inv = spdiags(1 ./ geo.v_area, 0, n_v, n_v);
bih = geo.lap * mass0_inv * geo.lap;
bih = blkdiag(bih, bih, bih);

[~, ~, ~, ~, KTK, DTD] = geo.evolving_operators();
surface_hess = 2 * (KTK + p.k * DTD) + 0.5 * p.dt * bih;

slp_cache = stokeslet_SLP_triangle_setup(M);
S = stokeslet_SLP_triangle_matrix(P, M, slp_cache);
%S = 0.5 * (M3 * S + S.' * M3);

N = normal_projection_matrix(geo.v_normal);
NN = M3 * N;

H = [
    surface_hess,              -I3;
    -I3,     p.dt^2 * p.Sd * (S + p.Da * NN)
];




fprintf("size(H) = %d x %d\n", size(H, 1), size(H, 2));
fprintf("relative symmetry error = %.3e\n", norm(H - H.', "fro") / norm(H, "fro"));

all_eigs = sort(eig(full(S)));
small_eigs = all_eigs(1:min(10, numel(all_eigs)));
fprintf("smallest eigenvalue = %.16e\n", all_eigs(1));

disp("smallest eigenvalues:");
disp(small_eigs);

[~, chol_flag] = chol(H + 1e-12 * speye(size(H)));
fprintf("chol flag after 1e-12 shift = %d\n", chol_flag);

function N = normal_projection_matrix(v_normal)
    n_v = size(v_normal, 1);
    rows = zeros(9 * n_v, 1);
    cols = zeros(9 * n_v, 1);
    vals = zeros(9 * n_v, 1);
    cursor = 0;

    for i = 1:n_v
        ni = v_normal(i, :).';
        block = ni * ni.';
        for row_dim = 1:3
            row = i + (row_dim - 1) * n_v;
            for col_dim = 1:3
                cursor = cursor + 1;
                rows(cursor) = row;
                cols(cursor) = i + (col_dim - 1) * n_v;
                vals(cursor) = block(row_dim, col_dim);
            end
        end
    end

    N = sparse(rows, cols, vals, 3 * n_v, 3 * n_v);
end