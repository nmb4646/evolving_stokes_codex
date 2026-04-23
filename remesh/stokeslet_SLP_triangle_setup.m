function cache = stokeslet_SLP_triangle_setup(F)
% Precompute topology-only data for the triangle-based Stokes SLP.

n_f = size(F, 1);
n_v = max(F(:));

cache.face_bary = [1/3, 1/3, 1/3];
cache.face_weights = 1;
cache.chunk_size = 128;

[duffy_nodes, duffy_weights] = unit_interval_gauss(4);
[U, V] = ndgrid(duffy_nodes, duffy_nodes);
[WU, WV] = ndgrid(duffy_weights, duffy_weights);
cache.duffy_u = U(:);
cache.duffy_v = V(:);
cache.duffy_w = WU(:) .* WV(:);

face_ids = repelem((1:n_f)', 3, 1);
local_ids = repmat((1:3)', n_f, 1);
incident_raw = accumarray(F(:), (1:numel(face_ids))', [n_v, 1], ...
    @(idx) { [face_ids(idx), local_ids(idx)] });

cache.incident = incident_raw;
end

function [nodes, weights] = unit_interval_gauss(n)
switch n
    case 2
        nodes = [-1; 1] / sqrt(3);
        weights = [1; 1];
    case 3
        nodes = [-sqrt(3/5); 0; sqrt(3/5)];
        weights = [5; 8; 5] / 9;
    case 4
        nodes = [
            -0.8611363115940526;
            -0.3399810435848563;
             0.3399810435848563;
             0.8611363115940526
        ];
        weights = [
            0.3478548451374538;
            0.6521451548625461;
            0.6521451548625461;
            0.3478548451374538
        ];
    otherwise
        error("Unsupported Gauss order %d.", n);
end

nodes = 0.5 * (nodes + 1);
weights = 0.5 * weights;
end
