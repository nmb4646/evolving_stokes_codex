function cache = stokeslet_SLP_advanced_setup(F)
% Topology-only cache for the advanced triangle-based Stokes SLP.
%
% Note to future edits:
% This cache exists to support the opt-in advanced SLP only. Keep the
% original triangle setup as the default unless specifically instructed.

n_f = size(F, 1);
n_v = max(F(:));

cache.is_advanced = true;
cache.chunk_size = 96;
cache.near_scale = 2.5;
cache.very_near_scale = 0.9;
cache.max_adaptive_depth = 3;
cache.face_bary = [1/3, 1/3, 1/3];
cache.face_weights = 1;

[cache.regular_bary, cache.regular_weights] = dunavant_rule_7();

[duffy_nodes, duffy_weights] = unit_interval_gauss(8);
[U, V] = ndgrid(duffy_nodes, duffy_nodes);
[WU, WV] = ndgrid(duffy_weights, duffy_weights);
cache.duffy_u = U(:);
cache.duffy_v = V(:);
cache.duffy_w = WU(:) .* WV(:);

face_ids = repelem((1:n_f)', 3, 1);
local_ids = repmat((1:3)', n_f, 1);
cache.incident = accumarray(F(:), (1:numel(face_ids))', [n_v, 1], ...
    @(idx) { [face_ids(idx), local_ids(idx)] });
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

function [nodes, weights] = unit_interval_gauss(n)
switch n
    case 8
        nodes = [
            -0.9602898564975363;
            -0.7966664774136267;
            -0.5255324099163290;
            -0.1834346424956498;
             0.1834346424956498;
             0.5255324099163290;
             0.7966664774136267;
             0.9602898564975363
        ];
        weights = [
            0.1012285362903763;
            0.2223810344533745;
            0.3137066458778873;
            0.3626837833783620;
            0.3626837833783620;
            0.3137066458778873;
            0.2223810344533745;
            0.1012285362903763
        ];
    otherwise
        error("Unsupported Gauss order %d.", n);
end

nodes = 0.5 * (nodes + 1);
weights = 0.5 * weights;
end
