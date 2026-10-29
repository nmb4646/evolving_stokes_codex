function output_path = curvature_flow_panels(run_folder, output_path, requested_ids, python_executable)
%CURVATURE_FLOW_PANELS Shaded curvature surface with curved velocity arrows.
%   curvature_flow_panels(RUN_FOLDER, OUTPUT_PATH, [35 85 135 185])
%   uses the PyVista renderer in curvature_flow_render.py. Its project-local
%   environment is .venv-curvature-flow. Override with MEMBRANE_RENDER_PYTHON
%   or the optional fourth argument. Both displayed fields use four averaging
%   passes. Curved arrows are integrated and projected onto the display mesh.

if nargin < 1 || isempty(run_folder)
    error('Supply a folder containing geoN.mat files.');
end
if nargin < 2 || isempty(output_path)
    output_path = fullfile(run_folder, 'curvature_flow_panels.png');
end
if nargin < 3 || isempty(requested_ids)
    saved = dir(fullfile(run_folder, 'geo*.mat'));
    ids = [];
    for i = 1:numel(saved)
        token = regexp(saved(i).name, '^geo(\d+)\.mat$', 'tokens', 'once');
        if ~isempty(token)
            ids(end+1) = str2double(token{1}); %#ok<AGROW>
        end
    end
    ids = unique(sort(ids));
    if numel(ids) < 4
        error('Expected at least four saved geoN.mat frames.');
    end
    if numel(ids) > 4
        ids = ids(2:end); % The initial state often stores zero velocity.
    end
    requested_ids = ids(round(linspace(1, numel(ids), 4)));
end
requested_ids = requested_ids(:).';
if numel(requested_ids) ~= 4 || any(~isfinite(requested_ids)) || ...
        any(requested_ids < 0 | requested_ids ~= fix(requested_ids))
    error('requested_ids must contain four nonnegative integer frame numbers.');
end

root = fileparts(mfilename('fullpath'));
if nargin < 4 || isempty(python_executable)
    python_executable = getenv('MEMBRANE_RENDER_PYTHON');
    if isempty(python_executable)
        python_executable = fullfile(root, '.venv-curvature-flow', 'bin', 'python');
    end
end
if ~isfile(python_executable)
    error(['Renderer Python was not found. Install the environment described ' ...
        'in curvature_flow_README.md or supply the Python path as argument 4.']);
end
output_path = char(output_path);
[out_dir, name, ext] = fileparts(output_path);
if ~strcmpi(ext, '.png')
    error('OUTPUT_PATH must end in .png.');
end
if isempty(out_dir)
    out_dir = pwd;
end
command_parts = {shell_quote(python_executable), ...
    shell_quote(fullfile(root, 'curvature_flow_render.py')), ...
    shell_quote(run_folder), '--output', shell_quote(out_dir), ...
    '--image-name', shell_quote([name ext]), '--frames', ...
    sprintf('%d ', requested_ids)};
status = system(strjoin(command_parts, ' '), '-echo');
if status ~= 0
    error('Surface renderer failed; see the diagnostic above.');
end
if ~isfile(output_path)
    error('The renderer did not create the requested output image.');
end
end

function quoted = shell_quote(value)
% POSIX quoting, including filenames that contain apostrophes.
quote = char(39);
quoted = [quote strrep(char(value), quote, ...
    [quote char(34) quote char(34) quote]) quote];
end
