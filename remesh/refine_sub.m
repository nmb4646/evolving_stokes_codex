function [multigrid_out] = refine_sub(geo_in,ftrlist_in)

    M_out = geo_in.F; P_out = geo_in.V; n_out = geo_in.v_normal; a_out = geo_in.v_area;
    for ftr = ftrlist_in(:)'
        vtr = M_out(ftr,:); M_out(ftr,:) = [];
        L = length(P_out);
        % find the 3 midpoints
        % for each midpoint, if it already exists, use the existing one;
        % otherwise, create a new one!
        v12 = .5*(P_out(vtr(1),:)+P_out(vtr(2),:));
        v23 = .5*(P_out(vtr(2),:)+P_out(vtr(3),:));
        v31 = .5*(P_out(vtr(3),:)+P_out(vtr(1),:));

        if ~any(P_out==v12,'all')
            
            P_out = [P_out; .5*(P_out(vtr(1),:)+P_out(vtr(2),:))]; a = length(P_out);
            n_out = [n_out; .5*(n_out(vtr(1),:)+n_out(vtr(2),:))];
            a_out = [a_out; .5*(a_out(vtr(1),:)+a_out(vtr(2),:))]; 
            %Last line is heuristic; this should be much much more complex
            %but it only works to increase the repulsive force on refined
            %areas so I can leave it for now i think
        else
            a = find(ismember(P_out, v12, 'rows'));
        end

        if ~any(P_out==v23,'all')
            P_out = [P_out; .5*(P_out(vtr(2),:)+P_out(vtr(3),:))]; c = length(P_out);
            n_out = [n_out; .5*(n_out(vtr(2),:)+n_out(vtr(3),:))];
            a_out = [a_out; .5*(a_out(vtr(2),:)+a_out(vtr(3),:))];
        else
            c = find(ismember(P_out, v23, 'rows'));
        end

        if ~any(P_out==v31,'all')
            P_out = [P_out; .5*(P_out(vtr(3),:)+P_out(vtr(1),:))]; b = length(P_out);
            n_out = [n_out; .5*(n_out(vtr(3),:)+n_out(vtr(1),:))];
            a_out = [a_out; .5*(a_out(vtr(3),:)+a_out(vtr(1),:))];
        else
            b = find(ismember(P_out, v31, 'rows'));
        end

        M_out = [M_out; [vtr(1),a,b]];
        M_out = [M_out; [vtr(2),a,c]];
        M_out = [M_out; [vtr(3),b,c]];
        M_out = [M_out; [a,b,c]];
    end
    multigrid_out.M = M_out;
    multigrid_out.P = P_out;
    multigrid_out.n = n_out;
    multigrid_out.a = a_out;

end