function data = interpolate(F, face, uv, vertex_data) 
    % linearly interpolate vertex data based on (face, barycentric coord).
    % 
    % Inputs:
    %   F: face matrix
    %   face: N by 1, face index
    %   uv: N by 2, point uv coordinate
    %   vertex_data: #V by M, where M is the dimension of the data type
    % Outputs:
    %   data: N by M, point uv coordinate
    vd1 = vertex_data(F(face, 1), :);
    vd2 = vertex_data(F(face, 2), :);
    vd3 = vertex_data(F(face, 3), :);
    data = vd1 .* uv(:, 1) + vd2 .* uv(:, 2) + vd3 .* (1 - uv(:, 1) - uv(:, 2));
end