close;clear;clc;

Omega = 20;
diptych = true; % True plots Omega = 20 and 30 side by side; false plots Omega only.
h = 1;
b = .1;
a = 1;


for i = 1 % Visualization parameters
visual.single_figure_position = [100, 100, 900, 650]; % [left bottom width height], pixels.
visual.diptych_figure_position = [100, 100, 1500, 650];
visual.figure_color = [1, 1, 1];
visual.tile_spacing = 'compact';
visual.tile_padding = 'compact';
visual.save_figure = true;
visual.output_filename = fullfile('.', 'data', 'figure', 'peanut_data.png');
visual.output_resolution = 300; % Dots per inch.

visual.data_line_style = '-';
visual.data_line_width = 2.5;
visual.data_marker = 'o';
visual.data_marker_size = 7;
visual.green = [0.10, 0.60, 0.20];
visual.red = [0.90, 0.10, 0.10];
visual.blue = [0.05, 0.30, 0.90];
visual.orange = [1.00, 0.45, 0.05];
visual.green_label = '0.5mM';
visual.red_label = '2.5mM';
visual.blue_label = '0.0mM';
visual.orange_label = '0.1mM';

visual.theory_line_style = '--';
visual.theory_line_width = 1.5;
visual.theory_line_color = [0, 0, 0];

visual.lubrication_line_style = '--';
visual.lubrication_line_width = 1.9;
visual.lubrication_line_color = [0.6, 0.6, 0.6];
visual.lubrication_line_color_2 = [0.3, 0.3, 0.3];

visual.x_scale = 'log';
visual.y_scale = 'log';
visual.x_limits = [5e-2, 1e1];
visual.y_min = 0.05;
visual.y_max_omega_factor = 0.5;
visual.x_label = '$\delta$';
visual.y_label = 'f_{orb}';
visual.x_label_interpreter = 'latex';
visual.y_label_interpreter = 'tex';
visual.label_font_size = 20;
visual.axes_font_size = 12;
visual.axes_line_width = 1.0;
visual.show_grid = false;
visual.show_title = true;
visual.title_format = 'f = %g Hz';
visual.title_interpreter = 'tex';
visual.title_font_size = 18;

visual.show_legend = true;
visual.legend_location = 'west';
visual.legend_font_size = 11;
end

if diptych
    omega_values = [20, 30];
    figure_position = visual.diptych_figure_position;
else
    omega_values = Omega;
    figure_position = visual.single_figure_position;
end

fig = figure('Position', figure_position, 'Color', visual.figure_color);
if diptych
    layout = tiledlayout(fig, 1, 2, ...
        'TileSpacing', visual.tile_spacing, 'Padding', visual.tile_padding);
end

for panel_index = 1:numel(omega_values)
Omega = omega_values(panel_index);

for i=1 % Data
if Omega == 20
red_data = [1.3104145601617796, 5.905956112852665
1.3993933265925178, 5.351097178683386
1.4600606673407484, 4.768025078369906
1.5369059656218402, 4.156739811912226
1.5854398382204247, 3.639498432601881
1.6825075834175935, 2.858934169278997
1.698685540950455, 2.7460815047021945
1.80788675429727, 2.6520376175548592
1.852376137512639, 2.238244514106583
1.896865520728008, 1.8244514106583074
2.0020222446916076, 1.6363636363636365
2.1112234580384226, 1.542319749216301
2.2527805864509607, 0.987460815047022
2.41051567239636, 0.7993730407523512
2.5035389282103138, 0.6300940438871474
2.697674418604651, 0.6489028213166145
2.802831142568251, 0.4890282131661442
2.9079878665318506, 0.5454545454545455
3.110212335692619, 0.2821316614420063
3.2072800808897877, 0.329153605015674
3.39737108190091, 0.18808777429467086];

green_data =[1.3993933265925178, 4.3730407523510975
1.500505561172902, 4.269592476489028
1.5935288169868556, 3.6112852664576804
1.621840242669363, 3.376175548589342
1.6663296258847322, 3.084639498432602
1.7350859453993934, 2.605015673981191
1.7755308392315472, 2.2664576802507836
1.8604651162790697, 2.0689655172413794
2.0060667340748233, 1.9561128526645768
2.058645096056623, 1.64576802507837
2.0950455005055613, 1.401253918495298
2.2366026289180994, 1.2695924764890283
2.309403437815976, 1.2413793103448276
2.3700707785642066, 0.8934169278996865
2.394337714863499, 0.7241379310344828
2.487360970677452, 0.7805642633228841
2.693629929221436, 0.5830721003134797
2.903943377148635, 0.3761755485893417
2.9929221435793734, 0.3479623824451411
3.0859453993933266, 0.3761755485893417
3.1951466127401416, 0.27272727272727276];

orange_data = [1.500505561172902, 3.4231974921630095
1.6097067745197169, 3.357366771159875
1.6946410515672397, 2.9905956112852667
1.8038422649140546, 2.7648902821316614
1.9130434782608696, 2.3322884012539187
2.091001011122346, 1.730407523510972
2.1638018200202227, 1.5611285266457682
2.281092012133468, 1.250783699059561
2.3579373104145605, 1.053291536050157
2.4792719919110215, 0.9404388714733543
2.5763397371081904, 0.7523510971786834
2.6895854398382206, 0.6300940438871474
2.806875631951466, 0.5360501567398119
2.916076845298281, 0.3761755485893417
3.01314459049545, 0.4702194357366771
3.2275025278058647, 0.30094043887147337
3.4337714863498485, 0.17868338557993732
3.660262891809909, 0.18808777429467086
3.9797775530839234, 0.17868338557993732
4.2103134479271995, 0.16927899686520376];

blue_data =[2.0748230535894843, 1.542319749216301
2.171890798786653, 1.523510971786834
2.281092012133468, 1.4200626959247649
2.3822042467138527, 1.213166144200627
2.4792719919110215, 0.9968652037617556
2.51567239635996, 0.9404388714733543
2.5803842264914056, 0.9028213166144201
2.6572295247724975, 0.7993730407523512
2.697674418604651, 0.7429467084639498
2.7866531850353895, 0.7429467084639498
2.8877654196157736, 0.5924764890282131
3.0091001011122347, 0.5454545454545455
3.1061678463094036, 0.45141065830721006
3.2072800808897877, 0.4420062695924765
3.308392315470172, 0.3103448275862069
3.4782608695652177, 0.29153605015673983
3.603640040444894, 0.30094043887147337
3.6966632962588473, 0.2821316614420063
3.801820020222447, 0.2163009404388715
3.9514661274014156, 0.20689655172413796
4.024266936299292, 0.17868338557993732
4.287158746208291, 0.15987460815047022
4.550050556117291, 0.14106583072100315
4.687563195146613, 0.10344827586206898];

elseif Omega == 30
    green_data = [1.2879256965944272, 6.666666666666667
1.4035087719298245, 6.87905604719764
1.502579979360165, 6.631268436578171
1.5768833849329205, 5.663716814159292
1.6429308565531475, 4.95575221238938
1.6842105263157894, 4.412979351032448
1.7791537667698658, 3.7640117994100293
1.8369453044375643, 3.3510324483775813
1.9236326109391124, 2.8436578171091447
1.9938080495356036, 2.359882005899705
2.0804953560371517, 1.935103244837758
2.2291021671826625, 1.4985250737463127
2.4231166150670793, 1.0619469026548671
2.6006191950464395, 0.9203539823008849
2.7698658410732713, 0.6253687315634218
2.897832817337461, 0.5309734513274336
3.0257997936016507, 0.4365781710914454
3.0835913312693495, 0.44837758112094395
3.1950464396284826, 0.3185840707964602
3.3023735810113517, 0.24778761061946902
3.4716202270381835, 0.2949852507374631
3.5872033023735805, 0.1887905604719764
3.7647058823529407, 0.1887905604719764
3.913312693498452, 0.10619469026548672
3.9958720330237356, 0.21238938053097345];
    red_data = [1.2961816305469556, 3.5634218289085546
1.4035087719298245, 5.156342182890856
1.502579979360165, 5.852507374631268
1.5975232198142413, 5.4159292035398225
1.6924664602683177, 4.743362831858407
1.7750257997936014, 3.9882005899705013
1.7997936016511866, 3.7168141592920354
1.8823529411764703, 3.410029498525074
1.9896800825593393, 3.150442477876106
2.0309597523219813, 2.9144542772861355
2.0763673890608874, 2.501474926253687
2.0970072239422084, 2.336283185840708
2.183694530443756, 2.0294985250737465
2.245614035087719, 1.6755162241887906
2.303405572755418, 1.2035398230088497
2.39422084623323, 1.0383480825958702
2.4974200206398347, 0.9321533923303835
2.65015479876161, 0.7079646017699115
2.707946336429308, 0.6253687315634218
2.7946336429308563, 0.5899705014749262
2.897832817337461, 0.3775811209439528
3.029927760577915, 0.2949852507374631
3.1248710010319916, 0.35398230088495575
3.1950464396284826, 0.40117994100294985
3.3023735810113517, 0.36578171091445427];

orange_data =[1.498452012383901, 5.2507374631268435
1.6016511867905054, 4.684365781710914
1.7048503611971102, 4.519174041297935
1.8163054695562435, 3.976401179941003
1.9442724458204332, 3.3392330383480826
2.0804953560371517, 2.7492625368731565
2.1424148606811144, 2.359882005899705
2.1960784313725488, 1.9587020648967552
2.2992776057791535, 1.6991150442477876
2.39422084623323, 1.887905604719764
2.4726522187822497, 1.4867256637168142
2.526315789473684, 1.2271386430678466
2.596491228070175, 0.9911504424778761
2.69969040247678, 1.0265486725663717
2.7905056759545923, 1.0147492625368733
2.9060887512899893, 0.7433628318584071
3.087719298245614, 0.6253687315634218];

blue_data =[1.9071207430340555, 2.6902654867256635
1.9938080495356036, 2.5368731563421827
2.0846233230134157, 2.3244837758112094
2.21671826625387, 1.9705014749262537
2.31578947368421, 1.6991150442477876
2.3983488132094943, 1.4867256637168142
2.4932920536635703, 1.4513274336283186
2.592363261093911, 1.2035398230088497
2.707946336429308, 1.0029498525073746
2.8813209494324044, 0.8849557522123894
2.976264189886481, 0.7787610619469026
3.199174406604747, 0.7079646017699115
3.3642930856553144, 0.49557522123893805
3.512899896800825, 0.3421828908554572
3.6986584107327136, 0.3185840707964602
3.7977296181630544, 0.3775811209439528
3.896800825593395, 0.3421828908554572
3.9876160990712073, 0.3421828908554572
4.094943240454076, 0.22418879056047197
4.2063983488132095, 0.22418879056047197
4.293085655314757, 0.2359882005899705
4.379772961816305, 0.24778761061946902
4.48297213622291, 0.15339233038348082
4.635706914344685, 0.15339233038348082
4.771929824561403, 0.10619469026548672
4.866873065015479, 0.11799410029498525];
end
end


red_data(:, 1) = red_data(:, 1) - 1; % Switch to lubrication distance.
green_data(:, 1) = green_data(:, 1) - 1;
blue_data(:, 1) = blue_data(:, 1) - 1;
orange_data(:, 1) = orange_data(:, 1) - 1;
Hz_vs_d = red_data; % Reference dataset for the scaling calculation below.


fit_range = 3:8; % Data-point indices used for the scaling fits.
scaling_fit = polyfit(log(Hz_vs_d(fit_range, 1)), ...
    log(Hz_vs_d(fit_range, 2)), 1);
alpha = scaling_fit(1);
fprintf('Points %d:%d give Hz ~ d^(%.6g).\n', ...
    fit_range(1), fit_range(end), alpha);



ds = logspace(-1,1,100);
Hz_theory = 2*Omega *((ds+1) .^ (-3) - ((ds+1).^2 + 4*h^2) .^ (-3/2));





if diptych
    ax = nexttile(layout);
else
    ax = axes(fig);
end
hold(ax, 'on')

h_green = plot(ax, green_data(:, 1), green_data(:, 2), ...
    'LineStyle', visual.data_line_style, 'Marker', visual.data_marker, ...
    'Color', visual.green, 'MarkerFaceColor', visual.green, ...
    'LineWidth', visual.data_line_width, 'MarkerSize', visual.data_marker_size, ...
    'DisplayName', visual.green_label);
h_red = plot(ax, red_data(:, 1), red_data(:, 2), ...
    'LineStyle', visual.data_line_style, 'Marker', visual.data_marker, ...
    'Color', visual.red, 'MarkerFaceColor', visual.red, ...
    'LineWidth', visual.data_line_width, 'MarkerSize', visual.data_marker_size, ...
    'DisplayName', visual.red_label);
h_blue = plot(ax, blue_data(:, 1), blue_data(:, 2), ...
    'LineStyle', visual.data_line_style, 'Marker', visual.data_marker, ...
    'Color', visual.blue, 'MarkerFaceColor', visual.blue, ...
    'LineWidth', visual.data_line_width, 'MarkerSize', visual.data_marker_size, ...
    'DisplayName', visual.blue_label);
h_orange = plot(ax, orange_data(:, 1), orange_data(:, 2), ...
    'LineStyle', visual.data_line_style, 'Marker', visual.data_marker, ...
    'Color', visual.orange, 'MarkerFaceColor', visual.orange, ...
    'LineWidth', visual.data_line_width, 'MarkerSize', visual.data_marker_size, ...
    'DisplayName', visual.orange_label);
h_theory = plot(ax, ds, Hz_theory, ...
    'LineStyle', visual.theory_line_style, ...
    'LineWidth', visual.theory_line_width, 'Color', visual.theory_line_color, ...
    'DisplayName', sprintf('Far-field model'));



deltas = logspace(-10,.1,200);
Hz_lubrication = (a./(deltas + 2*a)).*(log(deltas/a)./(log(deltas/a)+log(b/a)))*Omega;
h_lubrication = plot(ax, deltas, Hz_lubrication, ...
    'LineStyle', visual.lubrication_line_style, ...
    'LineWidth', visual.lubrication_line_width, ...
    'Color', visual.lubrication_line_color, ...
    'DisplayName', 'Lubrication model');


Hz_lubrication_2 = (a./(deltas + 2*a)).*(2*log(deltas/a)./(2*log(deltas/a)+log(b/a)))*Omega;
h_lubrication_2 = plot(ax, deltas, Hz_lubrication_2, ...
    'LineStyle', visual.lubrication_line_style, ...
    'LineWidth', visual.lubrication_line_width, ...
    'Color', visual.lubrication_line_color_2, ...
    'DisplayName', 'Lubrication model (two spheres)');


% Hz_swanbrady = rpb_two(a,b,deltas)*Omega;
% h_swanbrady = plot(ax, deltas, Hz_swanbrady, ...
%     'LineStyle', visual.lubrication_line_style, ...
%     'LineWidth', visual.lubrication_line_width, ...
%     'Color', [.9,.1,.1], ...
%     'DisplayName', 'Lubrication model (two spheres)');
box on;


set(ax, 'XScale', visual.x_scale, 'YScale', visual.y_scale, ...
    'FontSize', visual.axes_font_size, 'LineWidth', visual.axes_line_width)
xlim(ax, visual.x_limits)
ylim(ax, [visual.y_min, visual.y_max_omega_factor * Omega])
xlabel(ax, visual.x_label, 'Interpreter', visual.x_label_interpreter, ...
    'FontSize', visual.label_font_size)
ylabel(ax, visual.y_label, 'Interpreter', visual.y_label_interpreter, ...
    'FontSize', visual.label_font_size)
if visual.show_title
    title(ax, sprintf(visual.title_format, Omega), ...
        'Interpreter', visual.title_interpreter, ...
        'FontSize', visual.title_font_size)
end
if visual.show_grid
    grid(ax, 'on')
else
    grid(ax, 'off')
end
if visual.show_legend && (~diptych || panel_index == 1)
    legend(ax, [h_blue, h_orange, h_green, h_red, h_theory, h_lubrication, h_lubrication_2], ...
        'Location', visual.legend_location, ...
        'FontSize', visual.legend_font_size)
end
end

if visual.save_figure
    output_directory = fileparts(visual.output_filename);
    if ~isempty(output_directory) && ~isfolder(output_directory)
        mkdir(output_directory)
    end
    exportgraphics(fig, visual.output_filename, ...
        'Resolution', visual.output_resolution)
end
