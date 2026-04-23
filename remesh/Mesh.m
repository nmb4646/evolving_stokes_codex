classdef Mesh
    properties
        F {mustBeNumeric};
        n_v {mustBeNumeric};
        n_f {mustBeNumeric};
        n_he {mustBeNumeric};
        
        v_face {mustBeNumeric};
        v_n_he {mustBeNumeric};
        v_he {mustBeNumeric};

        f_he {mustBeNumeric};
        
        he_src {mustBeNumeric};
        he_dst {mustBeNumeric};
        he_face {mustBeNumeric};
        he_next {mustBeNumeric};
        he_prev {mustBeNumeric};
        he_flip {mustBeNumeric};
    end
    methods
        function obj = Mesh(F, V)
            obj.F = F;
            obj.n_v = size(V, 1);
            obj.n_f = size(F, 1); % number of faces
            obj.n_he = 3 * obj.n_f; % number of half-edges
            
            % Identity map
            I_f = (1:obj.n_f)'; % I_f: F -> F identity map
            I_he = (1:obj.n_he)'; % I_he: H -> H identity map
            
            % half-edge maps
            obj.he_src = [F(:,1); F(:,2); F(:,3)]; % src: H -> V source map
            obj.he_dst = [F(:,2); F(:,3); F(:,1)]; % dst: H -> V destination map
            obj.he_face = [I_f; I_f; I_f]; % face: H -> F face map

            % face maps
            obj.f_he = reshape(I_he, [obj.n_f, 3]);
            
            % halfedge traversal maps
            I_he_m = reshape(I_he, [], 3);
            he_next_m = [I_he_m(:,2), I_he_m(:,3), I_he_m(:,1)];
            obj.he_next = he_next_m(:); % he_next: H -> H next map
            obj.he_prev = obj.he_next(obj.he_next);
            adj = sparse(obj.he_src, obj.he_dst, I_he, obj.n_v, obj.n_v); % adj: (src, dst) -> H adjacency map
            obj.he_flip = full(adj(sub2ind([obj.n_v, obj.n_v], obj.he_dst, obj.he_src))); % he_flip: H -> H flip map
            
            % vertex halfedge iterator
            [obj.v_n_he, obj.v_he] = obj.vertex_halfedges();
        end

        function [v_n_he, v_he] = vertex_halfedges(obj)
            % number of halfedges at each vertex
            v_n_he = accumarray(obj.he_src, 1); 
            % first halfedge (by index) at each vertex
            v_he1 = accumarray(obj.he_src, (1:obj.n_he)', [obj.n_v, 1], @min); 
            % list of halfedges at each vertex 
            v_he = zeros(obj.n_v, max(v_n_he)); v_he(:, 1) = v_he1; %
            for i = 2:max(v_n_he)
                v_he(:, i) = obj.he_flip(obj.he_prev(v_he(:, i - 1)));
            end
        end

        function [sum, n_neighbor] = face_to_vertex(obj, data_f)
            sum = zeros(obj.n_v, size(data_f, 2));
            n_neighbor = accumarray(obj.F(:), 1, [obj.n_v, 1]);
            for i = 1:size(data_f, 2)
                sum(:, i) = accumarray(obj.F(:), [data_f(:, i); data_f(:, i); data_f(:, i)], [obj.n_v, 1]);
            end
        end
        
        function [sum, n_neighbor] = halfedge_to_vertex(obj, data_he)
            sum = zeros(obj.n_v, size(data_he, 2));
            n_neighbor = accumarray(obj.he_src(:), 1, [obj.n_v, 1]);
            for i = 1:size(data_he, 2)
                sum(:, i) = accumarray(obj.he_src(:), data_he(:, i), [obj.n_v, 1]);
            end
        end
    end
end