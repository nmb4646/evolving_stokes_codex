close all; clc;

% Visualize the current on-surface SLP quadrature scheme for one target
% vertex:
% - non-incident faces use one centroid quadrature point
% - incident faces use the 4x4 Duffy-transformed square rule

subdivisions = 5;
target_vertex_id = [];
show_all_centroids = true;
centroid_alpha = 0.18;
duffy_size = 42;
centroid_size = 8;
target_size = 140;

[P, F] = subdivided_sphere(subdivisions);

% Apply a smooth deterministic deformation so the mesh is not perfectly
% spherical.
r = vecnorm(P, 2, 2);
theta = atan2(P(:,2), P(:,1));
phi = acos(max(-1, min(1, P(:,3) ./ r)));
shape_factor = 1 + 0.22 * sin(2 * phi) .* cos(3 * theta) + 0.10 * P(:,3);
P = P .* shape_factor;

geo = Geometry(F, P);

if isempty(target_vertex_id)
    [~, target_vertex_id] = max(P(:,3));
end

incident_mask = any(F == target_vertex_id, 2);
incident_faces = find(incident_mask);
regular_faces = find(~incident_mask);

centroids = (P(F(:,1),:) + P(F(:,2),:) + P(F(:,3),:)) / 3;
regular_centroids = centroids(regular_faces, :);

[duffy_u, duffy_v] = duffy_nodes_4x4();
duffy_points = zeros(numel(incident_faces) * numel(duffy_u), 3);
duffy_face_id = zeros(numel(incident_faces) * numel(duffy_u), 1);

cursor = 1;
for k = 1:numel(incident_faces)
    face_id = incident_faces(k);
    local_idx = find(F(face_id,:) == target_vertex_id, 1);
    order = [local_idx, mod(local_idx, 3) + 1, mod(local_idx + 1, 3) + 1];
    verts = P(F(face_id, order), :);

    for q = 1:numel(duffy_u)
        u = duffy_u(q);
        v = duffy_v(q);
        bary = [1 - u, u * (1 - v), u * v];
        duffy_points(cursor, :) = bary(1) * verts(1,:) + bary(2) * verts(2,:) + bary(3) * verts(3,:);
        duffy_face_id(cursor) = face_id;
        cursor = cursor + 1;
    end
end

figure('Color', 'w');
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for panel = 1:2
    nexttile;
    trisurf(F, P(:,1), P(:,2), P(:,3), ...
        'FaceColor', [0.82, 0.85, 0.89], ...
        'EdgeColor', [0.55, 0.58, 0.63], ...
        'FaceAlpha', 0.92);
    hold on;
    axis equal;

    trisurf(F(incident_faces,:), P(:,1), P(:,2), P(:,3), ...
        'FaceColor', [0.99, 0.79, 0.30], ...
        'EdgeColor', [0.75, 0.48, 0.10], ...
        'FaceAlpha', 0.65);

    if show_all_centroids
        scatter3(regular_centroids(:,1), regular_centroids(:,2), regular_centroids(:,3), ...
            centroid_size, [0.22, 0.48, 0.88], 'filled', ...
            'MarkerFaceAlpha', 0.75, 'MarkerEdgeAlpha', centroid_alpha);
    end

    scatter3(duffy_points(:,1), duffy_points(:,2), duffy_points(:,3), ...
        duffy_size, [0.90, 0.22, 0.15], 'filled', 'MarkerEdgeColor', 'k');

    scatter3(P(target_vertex_id,1), P(target_vertex_id,2), P(target_vertex_id,3), ...
        target_size, [0.05, 0.05, 0.05], 'filled', 'MarkerEdgeColor', [1, 1, 1], 'LineWidth', 1.0);

    target_normal = geo.v_normal(target_vertex_id, :);
    quiver3(P(target_vertex_id,1), P(target_vertex_id,2), P(target_vertex_id,3), ...
        0.35 * target_normal(1), 0.35 * target_normal(2), 0.35 * target_normal(3), ...
        0, 'Color', [0.05, 0.05, 0.05], 'LineWidth', 1.4, 'MaxHeadSize', 2.5);

    title(sprintf('Quadrature For Target Vertex %d', target_vertex_id));
    xlabel('x'); ylabel('y'); zlabel('z');
    view(145, 24);
    grid on;
    camlight headlight;
    lighting gouraud;
end

legend({
    'Mesh', ...
    'Incident faces (Duffy-treated)', ...
    'Regular-face centroids', ...
    'Duffy quadrature points', ...
    'Target vertex'
    }, 'Location', 'southoutside', 'Orientation', 'horizontal');

fprintf('Target vertex: %d\n', target_vertex_id);
fprintf('Incident faces: %d\n', numel(incident_faces));
fprintf('Regular faces: %d\n', numel(regular_faces));
fprintf('Regular quadrature points: %d (1 centroid per non-incident face)\n', size(regular_centroids, 1));
fprintf('Duffy quadrature points: %d (%d per incident face)\n', size(duffy_points, 1), numel(duffy_u));

function [u_all, v_all] = duffy_nodes_4x4()
nodes = [
    -0.8611363115940526;
    -0.3399810435848563;
     0.3399810435848563;
     0.8611363115940526
];
nodes = 0.5 * (nodes + 1);
[U, V] = ndgrid(nodes, nodes);
u_all = U(:);
v_all = V(:);
end
