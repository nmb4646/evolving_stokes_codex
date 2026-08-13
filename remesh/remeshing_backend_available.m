function available = remeshing_backend_available(p)
%REMESHING_BACKEND_AVAILABLE Check the selected remesher without running it.

    backend = "legacy";
    if nargin >= 1 && isfield(p, "remesh_backend") && ~isempty(p.remesh_backend)
        backend = lower(string(p.remesh_backend));
    end
    script_dir = fileparts(mfilename("fullpath"));
    if backend == "legacy"
        addpath(fullfile(script_dir, "isoremesh"));
        available = exist("remeshing", "file") ~= 0;
    elseif backend == "mmgs"
        addpath(fullfile(script_dir, "mmgs"));
        available = exist("mmgs_remesh_mex", "file") == 3;
    else
        error("Unknown remesh_backend '%s'.", backend);
    end
end
