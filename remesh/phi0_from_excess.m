function [phi0_out] = phi0_from_excess(excess_in)
%PHI0_FROM_EXCESS Interpolate initial tilt phi0 from excess area.
%   phi0_out = phi0_from_excess(excess_in) uses tabulated data with excess
%   area in column 1 and phi0 in column 2. The output has the same shape as
%   excess_in.

kantsler_data = [.09, .25
    0.12358916478555304, 0.22176308539944906
0.18623024830699775, 0.20695592286501377
0.25225733634311515, 0.19765840220385675
0.33860045146726864, 0.1825068870523416
0.43510158013544015, 0.17389807162534435
0.5468397291196389, 0.16460055096418733
0.6670428893905191, 0.1553030303030303
0.7990970654627539, 0.14738292011019283
0.9446952595936794, 0.1391184573002755
1.0987584650112867, 0.13292011019283748
1.2714446952595937, 0.12637741046831957
1.5, .12];

[excess_data, order] = sort(kantsler_data(:, 1));
phi0_data = kantsler_data(order, 2);

excess_shape = size(excess_in);
excess_query = excess_in(:);

outside = excess_query < excess_data(1) | excess_query > excess_data(end);
if any(outside)
    warning("phi0_from_excess:OutOfRange", ...
        "Some excess-area inputs are outside the tabulated range [%.6g, %.6g]; returning NaN there.", ...
        excess_data(1), excess_data(end));
end

phi0_out = interp1(excess_data, phi0_data, excess_query, "pchip", NaN);
phi0_out = reshape(phi0_out, excess_shape);

end
