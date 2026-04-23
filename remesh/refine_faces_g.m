function geo_out = refine_faces_g(geo_in, ftr_list)
    geo_out = geo_in;
    F = geo_out.F;

    edge_map = containers.Map('KeyType','char','ValueType','int32');
    new_faces = [];

    for k = 1:length(ftr_list)
        ftr = ftr_list(k);
        tri = F(ftr, :);  % original face to refine
        v1 = tri(1); v2 = tri(2); v3 = tri(3);

        % Midpoints of each edge (with reuse)
        a = get_or_create_midpoint(v1, v2);
        b = get_or_create_midpoint(v2, v3);
        c = get_or_create_midpoint(v3, v1);

        % Maintain CCW orientation based on original face [v1, v2, v3]
        % Orientation-preserving subdivision:
        %   - triangle [v1, a, c]
        %   - triangle [v2, b, a]
        %   - triangle [v3, c, b]
        %   - triangle [a, b, c]

        new_faces = [new_faces;
            v1, a, c;
            v2, b, a;
            v3, c, b;
            a , b, c];
    end

    % Remove refined faces and append new ones
    F(ftr_list, :) = [];
    F = [F; new_faces];
    geo_out.F = F;

    % midpoint caching with consistent orientation
    function idx = get_or_create_midpoint(i, j)
        key = edge_key(i, j);
        if isKey(edge_map, key)
            idx = edge_map(key);
        else
            vi = geo_out.V(i, :);
            vj = geo_out.V(j, :);
            vmid = 0.5 * (vi + vj);
            geo_out.V(end+1, :) = vmid;
            idx = size(geo_out.V, 1);
            edge_map(key) = idx;
        end
    end
end

function key = edge_key(i, j)
    % Edge key that's independent of direction, for symmetry
    if i > j
        [i, j] = deal(j, i);
    end
    key = sprintf('%d_%d', i, j);
end
