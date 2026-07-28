clear; close; clc;

%%% Parameters
figure_position = [600, 200, 900, 700];
figure_color = 'w';

% Current-work simulation data to extract.
broussinos_Sd = 1e-6;
broussinos_Da = 0;
broussinos_gamy_values = [8e-7,8e-6,8e-5]; % Scalar or array, e.g. [4e-7, 8e-7, 1.2e-6].
broussinos_v_values = [0.988940, 0.97734, 0.96597, 0.954817, 0.9438, ...
    0.93314, 0.9226, 0.91227, 0.90213, 0.89217, 0.88240, 0.87280];
broussinos_tf = "last"; % "last" or a numeric frame id, e.g. 200.

% Axis data options.
kantsler_year = 2011; % Options: 2009 or 2011.
x_axis_mode = "Lam"; % Options: "excess_area" or "Lam".
tilt_scaled_by_pi = true; % true plots psi/pi; false plots psi in radians.

x_limits = [0.05, 1.5]; % Excess-area limits; auto-converted if x_axis_mode = "Lam".
y_limits = [0.09, 0.25]; % psi/pi limits; auto-converted if tilt_scaled_by_pi = false.
axis_line_width = 1.25;
axis_font_size = 18;
axis_font_name = 'Helvetica';
axis_font_weight = 'normal';

curve_line_width = 1.5;
scatter_marker_size = 100;
scatter_line_width = 2;

shaqfeh_color = 'k';
misbah_color = 'k';
kantsler_edge_color = 'k';
broussinos_edge_color = 'r';
broussinos_color_map = "lines";
broussinos_color_order = [.8 .1 .1;.1 .8 .1; .1 .1 .8]; % Optional nGamy x 3 RGB matrix. Leave [] to use broussinos_color_map.

if x_axis_mode == "excess_area"
    title_text = "Membrane shape vs tilt angle, \chi = 8";
else
    title_text = "Membrane shape vs tilt angle";
end
title_font_size = 30;
title_font_weight = 'normal';
x_label_text = "";
y_label_text = "";
label_font_size = 22;
label_font_weight = 'normal';

legend_font_size = 20;
legend_location = 'northeast';

if ~ismember(kantsler_year, [2009, 2011])
    error("kantsler_year must be either 2009 or 2011.");
end
if kantsler_year == 2011 ...
        && ~ismember(lower(string(x_axis_mode)), ["lam", "lambda"])
    error(["kantsler_year = 2011 requires x_axis_mode = ""Lam"" " ...
        "because its x data are already in Lambda coordinates."]);
end

shaqfeh_data = [0.09142212189616253, 0.21728650137741048
0.1574492099322799, 0.20695592286501377
0.21839729119638826, 0.1993801652892562
0.3081264108352144, 0.19042699724517906
0.3876975169300226, 0.18353994490358128
0.48250564334085777, 0.17665289256198347
0.577313769751693, 0.17079889807162535
0.6992099322799097, 0.16391184573002754
0.8278781038374717, 0.15736914600550964
0.9633182844243792, 0.15151515151515152
1.1241534988713318, 0.14462809917355374
1.327313769751693, 0.13705234159779614];

kantsler_data=[0.12358916478555304, 0.22176308539944906
0.18623024830699775, 0.20695592286501377
0.25225733634311515, 0.19765840220385675
0.33860045146726864, 0.1825068870523416
0.43510158013544015, 0.17389807162534435
0.5468397291196389, 0.16460055096418733
0.6670428893905191, 0.1553030303030303
0.7990970654627539, 0.14738292011019283
0.9446952595936794, 0.1391184573002755
1.0987584650112867, 0.13292011019283748
1.2714446952595937, 0.12637741046831957];

kantsler_data_2011=[0.185501066098081, 0.7474397590361446
0.15565031982942432, 0.7326807228915663
0.2046908315565032, 0.7305722891566265
0.18763326226012794, 0.7200301204819277
0.2601279317697228, 0.7189759036144578
0.2302771855010661, 0.7063253012048193
0.1812366737739872, 0.7052710843373494
0.2046908315565032, 0.6926204819277109
0.21321961620469082, 0.6768072289156627
0.5692963752665245, 0.5408132530120482
0.5948827292110874, 0.5460843373493975
0.4520255863539446, 0.560843373493976
0.4605543710021322, 0.5650602409638554
0.48187633262260127, 0.6367469879518073
0.4371002132196162, 0.6240963855421687
0.43496801705756927, 0.6019578313253012
0.4605543710021322, 0.5987951807228916
0.23880597014925373, 0.7073795180722892
0.5373134328358209, 0.48704819277108435
0.5565031982942431, 0.5123493975903615
0.5778251599147122, 0.5028614457831325
0.5970149253731343, 0.5060240963855421
0.6204690831556503, 0.49653614457831324
0.6055437100213219, 0.5334337349397591
0.5863539445628998, 0.5302710843373494
0.6567164179104478, 0.5186746987951807
0.6311300639658849, 0.5176204819277108
0.6439232409381663, 0.507078313253012
0.6439232409381663, 0.4786144578313253
0.6481876332622601, 0.46069277108433737
0.6268656716417911, 0.4585843373493976
0.7249466950959488, 0.5239457831325302
0.7398720682302772, 0.5123493975903615
0.6950959488272921, 0.5102409638554217
0.7078891257995735, 0.5007530120481928
0.7164179104477612, 0.49231927710843376
0.6673773987206822, 0.49337349397590363
0.6716417910447761, 0.4796686746987952
0.6886993603411514, 0.47545180722891567
0.6908315565031983, 0.49021084337349397
0.720682302771855, 0.4786144578313253
0.7078891257995735, 0.4628012048192771
0.7569296375266524, 0.4859939759036145
0.7654584221748401, 0.47545180722891567
0.7228144989339019, 0.43644578313253013
0.7356076759061834, 0.4628012048192771
0.7547974413646056, 0.4543674698795181
0.7462686567164178, 0.4438253012048193
0.7825159914712153, 0.45647590361445783
0.7739872068230277, 0.44066265060240967
0.7825159914712153, 0.4280120481927711
0.7569296375266524, 0.4311746987951807
0.7334754797441364, 0.411144578313253
0.7718550106609808, 0.40798192771084335
0.8166311300639658, 0.4164156626506024
0.7889125799573561, 0.4164156626506024
0.8251599147121536, 0.37846385542168676
0.8400852878464818, 0.3689759036144578
0.7782515991471215, 0.36159638554216866
0.8422174840085288, 0.39533132530120485
0.8720682302771855, 0.38162650602409637
0.9211087420042644, 0.37003012048192774
0.9509594882729211, 0.36370481927710846
0.9211087420042644, 0.34683734939759037
0.9339019189765458, 0.33629518072289155
0.9509594882729211, 0.3405120481927711
0.9808102345415778, 0.31837349397590364
0.9957356076759062, 0.32786144578313253
1.0063965884861408, 0.33629518072289155
0.9808102345415778, 0.35
0.906183368869936, 0.376355421686747
0.906183368869936, 0.38795180722891565
0.8742004264392323, 0.39849397590361446
0.8912579957356077, 0.39006024096385544
0.8550106609808102, 0.411144578313253
0.8315565031982942, 0.40798192771084335
0.8742004264392323, 0.4543674698795181
0.9275053304904051, 0.42484939759036144
0.9360341151385927, 0.4058734939759036
0.9552238805970149, 0.3963855421686747
1.0042643923240937, 0.39322289156626505
1.019189765458422, 0.39006024096385544
0.9872068230277186, 0.3868975903614458
0.9722814498933902, 0.39533132530120485
1.049040511727079, 0.37846385542168676
1.072494669509595, 0.36370481927710846
1.021321961620469, 0.36370481927710846
1.0106609808102345, 0.3774096385542169
1.0383795309168444, 0.3605421686746988
1.049040511727079, 0.36686746987951807
1.068230277185501, 0.3457831325301205
1.0426439232409381, 0.34262048192771083
1.0255863539445629, 0.34789156626506024
0.976545842217484, 0.3689759036144578
0.9381663113006397, 0.38795180722891565
0.9722814498933902, 0.3795180722891566
0.8933901918976546, 0.41325301204819276
0.8486140724946695, 0.4259036144578313
0.8891257995735607, 0.429066265060241
0.8720682302771855, 0.43328313253012046
0.8614072494669509, 0.44171686746987954
0.8081023454157782, 0.4628012048192771
0.8251599147121536, 0.45647590361445783
0.8187633262260128, 0.45015060240963856
0.8294243070362474, 0.43960843373493974
0.8102345415778252, 0.4375
0.7931769722814499, 0.4480421686746988
0.812366737739872, 0.42484939759036144
0.8017057569296375, 0.4322289156626506
0.8336886993603412, 0.41852409638554217
0.7910447761194029, 0.4058734939759036
0.9125799573560768, 0.3974397590361446
0.8869936034115138, 0.40165662650602413
0.8678038379530917, 0.41430722891566263
0.9125799573560768, 0.411144578313253
0.9445628997867803, 0.3753012048192771
0.929637526652452, 0.38162650602409637
0.9872068230277186, 0.36370481927710846
0.997867803837953, 0.3605421686746988
1.0042643923240937, 0.3521084337349398
0.9637526652452025, 0.35843373493975905
0.9466950959488273, 0.3457831325301205
1.072494669509595, 0.3173192771084337
0.720682302771855, 0.44066265060240967
0.744136460554371, 0.4227409638554217
0.6780383795309168, 0.4680722891566265
0.6929637526652452, 0.463855421686747
0.7142857142857143, 0.4743975903614458
0.742004264392324, 0.4733433734939759
0.7569296375266524, 0.46490963855421685
0.6460554371002132, 0.49126506024096384
1.091684434968017, 0.3510542168674699
1.0575692963752665, 0.3552710843373494
1.1002132196162047, 0.3436746987951807
1.115138592750533, 0.3373493975903614
1.1471215351812367, 0.32786144578313253
1.164179104477612, 0.315210843373494
1.185501066098081, 0.31204819277108437
1.2302771855010661, 0.30783132530120483
1.2643923240938166, 0.2941265060240964
1.2494669509594882, 0.288855421686747
1.2132196162046909, 0.3036144578313253
1.2302771855010661, 0.30783132530120483
1.236673773987207, 0.3088855421686747
1.3176972281449892, 0.2783132530120482
1.300639658848614, 0.26671686746987955
1.2643923240938166, 0.2572289156626506
1.2857142857142856, 0.26144578313253014
1.3091684434968016, 0.2635542168674699
1.4093816631130063, 0.22771084337349398
1.4371002132196162, 0.23614457831325303
1.5245202558635393, 0.20768072289156628
1.1279317697228144, 0.31099397590361444
1.115138592750533, 0.3246987951807229
1.0980810234541578, 0.315210843373494
1.0703624733475479, 0.3289156626506024
1.0533049040511726, 0.3341867469879518
1.091684434968017, 0.33207831325301207
1.1449893390191896, 0.3173192771084337
1.1897654584221748, 0.2825301204819277
1.2260127931769722, 0.2783132530120482
1.183368869936034, 0.2993975903614458
1.2196162046908314, 0.28990963855421686
1.1215351812366738, 0.28674698795180725];

misbah_data=[0.08803611738148984, 0.21694214876033058
0.1879232505643341, 0.20075757575757577
0.27257336343115124, 0.19008264462809918
0.4621896162528217, 0.17011019283746556
0.6992099322799097, 0.14910468319559228
0.9260722347629796, 0.13085399449035812
1.1580135440180586, 0.11225895316804407
1.3205417607223475, 0.09951790633608816];

[broussinos_data, broussinos_gamy_by_point] = extract_broussinos_data( ...
    broussinos_Sd, broussinos_Da, broussinos_gamy_values, ...
    broussinos_v_values, broussinos_tf);

[shaqfeh_data, x_label_text, y_label_text, x_limits, y_limits] = apply_axis_options( ...
    shaqfeh_data, x_axis_mode, tilt_scaled_by_pi, x_label_text, y_label_text, ...
    x_limits, y_limits);
[kantsler_data, kantsler_marker, kantsler_filled, kantsler_display_name] = ...
    select_kantsler_data(kantsler_year, kantsler_data, kantsler_data_2011, ...
    x_axis_mode, tilt_scaled_by_pi, x_label_text, y_label_text, ...
    x_limits, y_limits);
[misbah_data, ~, ~, ~, ~] = apply_axis_options( ...
    misbah_data, x_axis_mode, tilt_scaled_by_pi, x_label_text, y_label_text, ...
    x_limits, y_limits);
[broussinos_data, ~, ~, ~, ~] = apply_axis_options( ...
    broussinos_data, x_axis_mode, tilt_scaled_by_pi, x_label_text, y_label_text, ...
    x_limits, y_limits);

fig = figure("Color", figure_color, "Position", figure_position);
ax = axes(fig);
hold(ax, "on");

plot(ax, shaqfeh_data(:,1), shaqfeh_data(:,2), ...
    Color=shaqfeh_color, LineWidth=curve_line_width, ...
    DisplayName="Zhao & Shaqfeh simulations");
plot(ax, misbah_data(:,1), misbah_data(:,2), ...
    Color=misbah_color, LineStyle='--', LineWidth=curve_line_width, ...
    DisplayName="Misbah perturbation theory");
if kantsler_filled
    scatter(ax, kantsler_data(:,1), kantsler_data(:,2), scatter_marker_size, ...
        kantsler_marker, LineWidth=scatter_line_width, ...
        MarkerEdgeColor='k', MarkerFaceColor='k', ...
        DisplayName=kantsler_display_name);
else
    scatter(ax, kantsler_data(:,1), kantsler_data(:,2), scatter_marker_size, ...
        kantsler_marker, LineWidth=scatter_line_width, ...
        MarkerEdgeColor=kantsler_edge_color, ...
        DisplayName=kantsler_display_name);
end
plot_broussinos_scatter(ax, broussinos_data, broussinos_gamy_by_point, ...
    broussinos_gamy_values, broussinos_color_order, broussinos_color_map, ...
    broussinos_edge_color, scatter_marker_size, scatter_line_width);

title(ax, title_text, Interpreter='tex', FontSize=title_font_size, ...
    FontWeight=title_font_weight);
xlabel(ax, x_label_text, FontSize=label_font_size, FontWeight=label_font_weight);
ylabel(ax, y_label_text, Interpreter='tex', FontSize=label_font_size, ...
    FontWeight=label_font_weight);
xlim(ax, x_limits);
ylim(ax, y_limits);
set(ax, LineWidth=axis_line_width, FontSize=axis_font_size, ...
    FontName=axis_font_name, FontWeight=axis_font_weight);
box on;
legend(ax, FontSize=legend_font_size, Location=legend_location);


if x_axis_mode == "excess_area"
    exportgraphics(fig, "./data/figures/figure2a.png", Resolution=300);
else
    exportgraphics(fig, "./data/figures/figure2b.png", Resolution=300);
end

function [data, marker, filled, display_name] = select_kantsler_data( ...
        year, data_2009, data_2011, x_axis_mode, tilt_scaled_by_pi, ...
        x_label_text, y_label_text, x_limits, y_limits)
    validateattributes(year, {'numeric'}, ...
        {'scalar', 'integer'}, mfilename, "kantsler_year");

    switch year
        case 2009
            data = apply_axis_options(data_2009, x_axis_mode, ...
                tilt_scaled_by_pi, x_label_text, y_label_text, ...
                x_limits, y_limits);
            marker = 'o';
            filled = false;
            display_name = "Kantsler & Steinberg measurements";
        case 2011
            if ~ismember(lower(string(x_axis_mode)), ["lam", "lambda"])
                error(["kantsler_year = 2011 requires x_axis_mode = ""Lam"" " ...
                    "because its x data are already in Lambda coordinates."]);
            end
            data = data_2011;
            if tilt_scaled_by_pi
                data(:, 2) = data(:, 2) / pi;
            end
            marker = 's';
            filled = true;
            display_name = "Kantsler & Steinberg 2011 measurements";
        otherwise
            error("kantsler_year must be either 2009 or 2011.");
    end
end

function [data, x_label_text, y_label_text, x_limits, y_limits] = apply_axis_options( ...
        data, x_axis_mode, tilt_scaled_by_pi, x_label_text, y_label_text, ...
        x_limits, y_limits)
    x_axis_mode = lower(string(x_axis_mode));
    tilt_scaled_by_pi = logical(tilt_scaled_by_pi);

    switch x_axis_mode
        case {"excess_area", "excess", "ea"}
            default_x_label = "Membrane excess area";
        case {"lam", "lambda"}
            data(:, 1) = Lam_from_ea(data(:, 1));
            if ~isempty(x_limits)
                x_limits = sort(Lam_from_ea(x_limits));
            end
            default_x_label = "\Lambda";
        otherwise
            error('Unknown x_axis_mode "%s". Use "excess_area" or "Lam".', x_axis_mode);
    end

    if tilt_scaled_by_pi
        default_y_label = "Tilt angle \psi/\pi";
    else
        data(:, 2) = pi * data(:, 2);
        if ~isempty(y_limits)
            y_limits = pi * y_limits;
        end
        default_y_label = "Tilt angle \psi";
    end

    if strlength(string(x_label_text)) == 0
        x_label_text = default_x_label;
    end
    if strlength(string(y_label_text)) == 0
        y_label_text = default_y_label;
    end
end

function [broussinos_data, gamy_by_point] = extract_broussinos_data(Sd, Da, gamy_values, v_values, tf)
    script_dir = fileparts(mfilename("fullpath"));
    if strlength(script_dir) == 0
        script_dir = pwd;
    end
    remesh_dir = find_remesh_dir(script_dir);
    addpath(remesh_dir);
    data_dir = find_data_dir(script_dir, remesh_dir);

    gamy_values = gamy_values(:).';
    v_values = v_values(:).';
    broussinos_data = NaN(numel(gamy_values) * numel(v_values), 2);
    gamy_by_point = NaN(numel(gamy_values) * numel(v_values), 1);

    row = 0;
    for g = 1:numel(gamy_values)
        gamy = gamy_values(g);
        for i = 1:numel(v_values)
            row = row + 1;
            v = v_values(i);
            run_tag = make_run_tag(Sd, Da, gamy, v);
            folder = fullfile(data_dir, run_tag);
            frames = list_geo_frames(folder);
            if isempty(frames)
                warning("No geo*.mat files found for gamy = %.6g, v = %.6g in %s", ...
                    gamy, v, folder);
                continue
            end

            frame_id = select_frame(frames, tf);
            if isempty(frame_id)
                warning("Requested frame %s not found for gamy = %.6g, v = %.6g in %s", ...
                    string(tf), gamy, v, folder);
                continue
            end

            data = load(fullfile(folder, sprintf("geo%d.mat", frame_id)), "M", "P");
            geo = Geometry(data.M, data.P);
            tilt = vesicleTiltDeformation(data.P, data.M, 1, 3, "ray_intersection");
            broussinos_data(row, :) = [excess_area(geo), tilt.psi_over_pi];
            gamy_by_point(row) = gamy;
        end
    end

    valid = all(isfinite(broussinos_data), 2) & isfinite(gamy_by_point);
    broussinos_data = broussinos_data(valid, :);
    gamy_by_point = gamy_by_point(valid);
end

function plot_broussinos_scatter(ax, data, gamy_by_point, gamy_values, color_order, color_map, ...
        single_edge_color, marker_size, line_width)
    if isempty(data)
        warning("No Broussinos/current-work data was extracted.");
        return
    end

    gamy_values = gamy_values(:).';
    if isempty(color_order)
        if isscalar(gamy_values)
            color_order = single_edge_color;
        else
            color_order = feval(char(color_map), numel(gamy_values));
        end
    end

    for g = 1:numel(gamy_values)
        gamy = gamy_values(g);
        mask = abs(gamy_by_point - gamy) <= max(eps(abs(gamy)), eps);
        if ~any(mask)
            continue
        end

        if ischar(color_order) || isstring(color_order)
            marker_color = color_order;
        else
            color_idx = mod(g - 1, size(color_order, 1)) + 1;
            marker_color = color_order(color_idx, :);
        end

        if isscalar(gamy_values)
            display_name = "Current work";
        else
            display_name = sprintf("Current work, \\chi= %.3g", gamy/(1e-6));
        end

        scatter(ax, data(mask, 1), data(mask, 2), marker_size, ...
            LineWidth=line_width, MarkerEdgeColor=marker_color, ...
            DisplayName=display_name);
    end
end

function frame_id = select_frame(frames, tf)
    frame_id = [];
    if isstring(tf) || ischar(tf)
        tf = string(tf);
        if lower(tf) == "last"
            frame_id = max(frames);
            return
        end
        requested = str2double(tf);
    else
        requested = double(tf);
    end

    if ismember(requested, frames)
        frame_id = requested;
    end
end

function frames = list_geo_frames(folder)
    files = dir(fullfile(folder, "geo*.mat"));
    frames = zeros(numel(files), 1);
    keep = false(numel(files), 1);
    for i = 1:numel(files)
        token = regexp(files(i).name, "^geo(\d+)\.mat$", "tokens", "once");
        if ~isempty(token)
            frames(i) = str2double(token{1});
            keep(i) = true;
        end
    end
    frames = sort(frames(keep));
end

function run_tag = make_run_tag(Sd, Da, gamy, v)
    run_tag = sprintf("Sd_%.2e_Da_%.2e_gamy_%+.2e_v_%.2e", Sd, Da, gamy, v);
    run_tag = strrep(run_tag, "+", "p");
    run_tag = strrep(run_tag, "-", "m");
end

function remesh_dir = find_remesh_dir(script_dir)
    candidates = string({script_dir, pwd, fullfile(pwd, "remesh")});
    geometry_path = which("Geometry");
    if strlength(geometry_path) > 0
        candidates(end + 1) = string(fileparts(geometry_path));
    end

    for i = 1:numel(candidates)
        candidate = candidates(i);
        if isfile(fullfile(candidate, "Geometry.m"))
            remesh_dir = char(candidate);
            return
        end
    end

    error("Could not locate remesh directory. Run from the repo root or add remesh/ to the MATLAB path.");
end

function data_dir = find_data_dir(script_dir, remesh_dir)
    candidates = string({
        fullfile(remesh_dir, "data", "fs_batch_data"), ...
        fullfile(script_dir, "data", "fs_batch_data"), ...
        fullfile(pwd, "data", "fs_batch_data"), ...
        fullfile(pwd, "remesh", "data", "fs_batch_data")});

    for i = 1:numel(candidates)
        candidate = candidates(i);
        if isfolder(candidate)
            data_dir = char(candidate);
            return
        end
    end

    error("Could not locate data/fs_batch_data. Checked:%s", sprintf("\n  %s", candidates));
end
