function [P, geo, info] = correct_remesh_constraints(geo, target_area, target_volume, p)
%CORRECT_REMESH_CONSTRAINTS Select legacy or mass-weighted correction.

    if nargin < 4 || isempty(p)
        p = struct();
    end
    mode = lower(string(field_or(p, "constraint_projection", "legacy")));
    options = field_or(p, "constraint_projection_options", struct());
    if mode == "legacy" || mode == "unweighted"
        [P, geo, info] = newton_correct_volume(geo, target_area, target_volume, options);
        info.method = "legacy";
    elseif mode == "mass" || mode == "mass_weighted" || mode == "weighted"
        targets = struct("area", target_area, "volume", target_volume);
        options.weighting = "mass";
        [P, geo, info] = project_surface_constraints(geo.V, geo.F, targets, options);
        if ~info.converged && field_or(options, "warn_on_failure", true)
            warning("correct_remesh_constraints:ProjectionNotConverged", ...
                ["Mass-weighted area/volume projection did not fully converge " ...
                "(%s): relA %.3e, relV %.3e."], ...
                info.reason, info.rel_area_error, info.rel_volume_error);
        end
    elseif mode == "none" || mode == "off"
        P = geo.V;
        info = struct("method", "none", "converged", true, "iterations", 0, ...
            "rel_area_error", abs(geo.area - target_area) / abs(target_area), ...
            "rel_volume_error", abs(geo.volume - target_volume) / abs(target_volume));
    else
        error("Unknown constraint_projection '%s'.", mode);
    end
end

function value = field_or(data, field, fallback)
    if isfield(data, field) && ~isempty(data.(field))
        value = data.(field);
    else
        value = fallback;
    end
end
