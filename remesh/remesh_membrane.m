function [M2, P2, info] = remesh_membrane(M, P, base_h, n_iter, p)
%REMESH_MEMBRANE Thin solver-facing wrapper around remesh_surface.

    if nargin < 5 || isempty(p)
        p = struct();
    end
    backend = string(field_or(p, "remesh_backend", "legacy"));
    options = field_or(p, "remesh_options", struct());
    options.backend = backend;
    options.legacy_iterations = n_iter;

    geo = Geometry(double(M), double(P));
    if lower(backend) == "mmgs"
        [target_h, target_info] = compute_target_edge_length(geo, base_h, options);
    else
        target_h = repmat(double(base_h), size(P, 1), 1);
        target_info = struct("adaptive", false, "base_h", base_h, ...
            "h_min", base_h, "h_mean", base_h, "h_max", base_h);
    end
    [P2, M2, info] = remesh_surface(P, M, target_h, options);
    info.target_size = target_info;
end

function value = field_or(data, field, fallback)
    if isfield(data, field) && ~isempty(data.(field))
        value = data.(field);
    else
        value = fallback;
    end
end
