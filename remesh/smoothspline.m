function [xSmooth, ySmooth] = smoothing(x, y, nPoints)
% smoothInterp1 - Smooth a 2D line using shape-preserving cubic interpolation
%
% Syntax:
%   [xSmooth, ySmooth] = smoothInterp1(x, y)
%   [xSmooth, ySmooth] = smoothInterp1(x, y, nPoints)
%
% Inputs:
%   x       - x data (vector)
%   y       - y data (vector)
%   nPoints - number of points in the smoothed curve (default: 300)
%
% Outputs:
%   xSmooth - fine x values
%   ySmooth - smoothed y values

    if nargin < 3
        nPoints = 300; % default number of points
    end
    
    % Create a fine grid along x
    xSmooth = linspace(min(x), max(x), nPoints);
    
    % Interpolate using shape-preserving cubic (pchip)
    ySmooth = interp1(x, y, xSmooth, 'pchip');
end
