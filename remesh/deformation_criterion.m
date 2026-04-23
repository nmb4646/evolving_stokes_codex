function out = deformation_criterion(geo_in)
    if min(geo_in.corner_angle) <.74
        out=true;
    else
        out=false;
    end
end

