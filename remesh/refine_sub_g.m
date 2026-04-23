function geo_out = refine_sub_g(geo_in, ftr)
    geo_out = geo_in;
    vtr = geo_out.F(ftr, :);
    geo_out.F(ftr, :) = [];

    edge_map = containers.Map('KeyType','char','ValueType','int32');
    new_vertices = [];
    
    % Helper function: get or create midpoint vertex
    function vidx = get_or_create_midpoint(i, j)
        key = edge_key(i, j);
        if isKey(edge_map, key)
            vidx = edge_map(key);
        else
            vi = geo_out.V(i, :);
            vj = geo_out.V(j, :);
            vmid = 0.5 * (vi + vj);
            geo_out.V(end+1, :) = vmid;
            vidx = size(geo_out.V, 1);
            edge_map(key) = vidx;
        end
    end

    % Create new midpoints (with caching)
    a = get_or_create_midpoint(vtr(1), vtr(2));
    b = get_or_create_midpoint(vtr(2), vtr(3));
    c = get_or_create_midpoint(vtr(3), vtr(1));

    % Create 4 new faces
    geo_out.F(end+1, :) = [vtr(1), a, c];
    geo_out.F(end+1, :) = [vtr(2), b, a];
    geo_out.F(end+1, :) = [vtr(3), c, b];
    geo_out.F(end+1, :) = [a, b, c];
end

% Utility function to generate consistent key for an edge
function key = edge_key(i, j)
    if i > j
        [i, j] = deal(j, i);
    end
    key = sprintf('%d_%d', i, j);
end