function [M_refined, P_refined] = refine_sub_gg(M, P, ftr)
    % M: Nx3 array of vertex coordinates
    % P: Mx3 array of triangle indices (faces)
    % ftr: refinement parameter or index (e.g., which edges/faces to refine)
    %
    % Returns:
    % M_refined, P_refined: refined mesh preserving orientation and no duplicates

    % Example: simple midpoint refinement on edges associated with ftr

    % 1. Find edges to refine based on ftr
    % (For demo, assume ftr indicates a set of face indices to refine)
    faces_to_refine = ftr;

    % 2. Extract all edges of these faces
    edges = [P(faces_to_refine, [1 2]);
             P(faces_to_refine, [2 3]);
             P(faces_to_refine, [3 1])];

    % Sort edges so [min max] to identify unique edges
    edges = sort(edges, 2);

    % 3. Find unique edges
    [unique_edges, ~, ic] = unique(edges, 'rows');

    % 4. Compute midpoints of unique edges
    midpoints = (M(unique_edges(:,1), :) + M(unique_edges(:,2), :)) / 2;

    % 5. Append midpoints to vertex list
    M_refined = [M; midpoints];

    % 6. Create mapping from old edges to new midpoint indices
    midpoint_indices = (size(M,1)+1):(size(M,1)+size(midpoints,1));

    % 7. For each face to refine, split into 4 smaller faces
    P_refined = P;
    new_faces = [];

    for i = 1:length(faces_to_refine)
        f_idx = faces_to_refine(i);
        verts = P(f_idx, :);

        % Find midpoints for edges of this face
        e1 = sort([verts(1), verts(2)]);
        e2 = sort([verts(2), verts(3)]);
        e3 = sort([verts(3), verts(1)]);

        % Find indices of midpoints in unique_edges
        mp1 = midpoint_indices(find(ismember(unique_edges, e1, 'rows')));
        mp2 = midpoint_indices(find(ismember(unique_edges, e2, 'rows')));
        mp3 = midpoint_indices(find(ismember(unique_edges, e3, 'rows')));

        % Build 4 new faces preserving orientation
        % Original vertices: v1, v2, v3
        v1 = verts(1);
        v2 = verts(2);
        v3 = verts(3);

        % New faces:
        % Triangle 1: v1, mp1, mp3
        % Triangle 2: mp1, v2, mp2
        % Triangle 3: mp3, mp2, v3
        % Triangle 4: mp1, mp2, mp3

        new_faces = [new_faces;
                     v1, mp1, mp3;
                     mp1, v2, mp2;
                     mp3, mp2, v3;
                     mp1, mp2, mp3];
    end

    % 8. Replace the refined faces with new ones
    % Remove old faces being refined
    P_refined(faces_to_refine, :) = [];

    % Append new refined faces
    P_refined = [P_refined; new_faces];
end
