function [data_out] = smoothing(data, nPoints)
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

    x = data(:,1);
    y = data(:,2);
    
    % Create a fine grid along x
    xSmooth = linspace(min(x), max(x), nPoints);
    
    
    % Interpolate using shape-preserving cubic (pchip)
    ySmooth = interp1(x, y, xSmooth, 'pchip');

    data_out = [xSmooth',ySmooth'];
end
