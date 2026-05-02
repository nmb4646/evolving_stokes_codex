function [vertex, velocity] = rm_rigid_patched(vertex, velocity, v_area, opts)
    % Remove selected rigid motion from a geometry update.
    %
    % Modes:
    %   "translation"  subtract area-weighted translational velocity and
    %                  recenter vertex.
    %   "all"          additionally remove finite rigid rotation from vertex
    %                  using a weighted Procrustes alignment to opts.P0.
    %
    % For mode "all", pass opts as a struct with:
    %   opts.mode = "all";
    %   opts.P0   = previous-step vertices, n_v by 3
    %   opts.dt   = timestep
    %
    % This keeps the saved displacement and velocity in the same frame:
    %   velocity = (vertex_aligned - centered_P0) / dt.

    if nargin < 4 || isempty(opts)
        opts = "all";
    end

    if ischar(opts) || isstring(opts)
        mode = string(opts);
        P0 = [];
        dt = [];
    elseif isstruct(opts)
        if isfield(opts, 'mode')
            mode = string(opts.mode);
        else
            mode = "all";
        end
        if isfield(opts, 'P0')
            P0 = opts.P0;
        else
            P0 = [];
        end
        if isfield(opts, 'dt')
            dt = opts.dt;
        else
            dt = [];
        end
    else
        error("opts must be a string or struct.");
    end

    n_v = size(vertex, 1);
    velocity = reshape(velocity, n_v, 3);
    v_area = v_area(:);
    mass = sum(v_area);

    com = sum(vertex .* v_area, 1) / mass;
    u_L = sum(velocity .* v_area, 1) / mass;

    vertex = vertex - com;
    velocity = velocity - u_L;

    switch mode
        case {"translation", "trans", "none"}
            velocity = velocity(:);

        case {"all", "rigid", "translation_rotation"}
            if isempty(P0) || isempty(dt)
                [vertex, velocity] = remove_instantaneous_rotation(vertex, velocity, v_area);
                velocity = velocity(:);
                return
            end

            P0 = reshape(P0, n_v, 3);
            com0 = sum(P0 .* v_area, 1) / mass;
            P0 = P0 - com0;

            R = weighted_alignment_rotation(vertex, P0, v_area);
            vertex = vertex * R;
            velocity = (vertex - P0) / dt;
            velocity = velocity(:);

        otherwise
            error("Unknown rigid-motion mode '%s'. Use 'translation' or 'all'.", mode);
    end
end

function [vertex, velocity] = remove_instantaneous_rotation(vertex, velocity, v_area)
    n_v = size(vertex, 1);

    Rsq_id = eye(3);
    Rsq_id = permute(Rsq_id(:, :, ones(n_v, 1)), [3, 1, 2]);
    Rsq_id = Rsq_id .* sum(vertex .* vertex, 2);
    VVT = vertex(:, :, ones(3, 1)) .* permute(vertex(:, :, ones(3, 1)), [1, 3, 2]);
    moment = sum((Rsq_id - VVT) .* v_area, 1);

    ang_momentum = sum(cross(vertex, velocity, 2) .* v_area, 1);
    w = (squeeze(moment) \ ang_momentum')';
    velocity = velocity - cross(repmat(w, n_v, 1), vertex, 2);
end

function R = weighted_alignment_rotation(X, Y, weights)
    % Row-vector convention: X * R best aligns with Y.
    H = X' * (weights .* Y);
    [U, ~, V] = svd(H);
    R = U * V';
    if det(R) < 0
        V(:, end) = -V(:, end);
        R = U * V';
    end
end
