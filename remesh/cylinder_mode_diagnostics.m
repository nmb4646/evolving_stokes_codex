function cylinder_mode_diagnostics(result, output_dir, cfg)
%CYLINDER_MODE_DIAGNOSTICS Save geometry, spectrum, and fit diagnostics.

if ~cfg.diagnostics.enabled
    return
end
if ~isfolder(output_dir)
    mkdir(output_dir);
end

save_geometry_metrics(result, output_dir, cfg);
save_core_profiles(result, output_dir, cfg);
save_mode_amplitudes(result, output_dir, cfg);
save_growth_dispersion(result, output_dir, cfg);
save_spectra(result, output_dir, cfg);
if cfg.diagnostics.save_individual_fit_plots
    save_individual_fits(result, output_dir, cfg);
end
end

function save_geometry_metrics(result, output_dir, cfg)
    metrics = result.frame_metrics;
    t = metrics.time;
    fig = new_figure(cfg);
    layout = tiledlayout(fig, 3, 2, "TileSpacing", "compact", "Padding", "compact");

    ax = nexttile(layout);
    plot(ax, t, metrics.surface_area / metrics.surface_area(1) - 1, ...
        "LineWidth", cfg.diagnostics.line_width);
    ylabel(ax, "\Delta A/A_0"); style_axes(ax, cfg);

    ax = nexttile(layout);
    plot(ax, t, metrics.volume / metrics.volume(1) - 1, ...
        "LineWidth", cfg.diagnostics.line_width);
    ylabel(ax, "\Delta V/V_0"); style_axes(ax, cfg);

    ax = nexttile(layout);
    yyaxis(ax, "left");
    plot(ax, t, metrics.mean_core_radius, "LineWidth", cfg.diagnostics.line_width);
    ylabel(ax, "R_0");
    yyaxis(ax, "right");
    plot(ax, t, metrics.core_length, "LineWidth", cfg.diagnostics.line_width);
    ylabel(ax, "L_c");
    style_axes(ax, cfg);

    ax = nexttile(layout);
    semilogy(ax, t, max(metrics.cylindricality_error, realmin), ...
        "LineWidth", cfg.diagnostics.line_width, "DisplayName", "cylindricality");
    hold(ax, "on");
    semilogy(ax, t, max(metrics.reconstruction_error, realmin), ...
        "LineWidth", cfg.diagnostics.line_width, "DisplayName", "reconstruction");
    ylabel(ax, "Relative error");
    legend(ax, "Location", "best");
    style_axes(ax, cfg);

    ax = nexttile(layout);
    yyaxis(ax, "left");
    plot(ax, t, metrics.core_area_fraction, ...
        "LineWidth", cfg.diagnostics.line_width, "DisplayName", "area fraction");
    hold(ax, "on");
    plot(ax, t, metrics.core_axial_fraction, "--", ...
        "LineWidth", cfg.diagnostics.line_width, "DisplayName", "axial fraction");
    ylabel(ax, "Core fraction");
    yyaxis(ax, "right");
    plot(ax, t, metrics.axis_angle_change, ...
        "LineWidth", cfg.diagnostics.line_width, "DisplayName", "axis change");
    ylabel(ax, "Axis change (deg)");
    legend(ax, "Location", "best");
    style_axes(ax, cfg);

    ax = nexttile(layout);
    semilogy(ax, t, max(metrics.projection_condition, 1), ...
        "LineWidth", cfg.diagnostics.line_width, "DisplayName", "condition");
    hold(ax, "on");
    semilogy(ax, t, max(metrics.conjugacy_error, realmin), ...
        "LineWidth", cfg.diagnostics.line_width, "DisplayName", "conjugacy error");
    ylabel(ax, "Projection diagnostics");
    legend(ax, "Location", "best");
    style_axes(ax, cfg);

    xlabel_all(layout, "Physical time");
    title(layout, replace(result.series_id, "_", "\_"), ...
        "FontSize", cfg.diagnostics.title_font_size);
    export_figure(fig, fullfile(output_dir, "cylinder_geometry_metrics.png"), cfg);
end

function save_core_profiles(result, output_dir, cfg)
    profiles = result.profiles;
    profiles = profiles(~cellfun(@isempty, profiles));
    if isempty(profiles)
        return
    end

    fig = new_figure(cfg);
    layout = tiledlayout(fig, numel(profiles), 1, ...
        "TileSpacing", "compact", "Padding", "compact");
    colors = lines(numel(profiles));
    for k = 1:numel(profiles)
        profile = profiles{k};
        ax = nexttile(layout);
        scatter(ax, profile.x(~profile.core_mask), profile.r(~profile.core_mask), ...
            6, [0.72, 0.72, 0.72], "filled");
        hold(ax, "on");
        scatter(ax, profile.x(profile.core_mask), profile.r(profile.core_mask), ...
            8, colors(k, :), "filled");
        xline(ax, profile.x_min, "--", "Color", colors(k, :));
        xline(ax, profile.x_max, "--", "Color", colors(k, :));
        ylabel(ax, "r");
        title(ax, sprintf("Frame %d", profile.frame_index));
        style_axes(ax, cfg);
    end
    xlabel(nexttile(layout, numel(profiles)), "Aligned x");
    title(layout, "Detected cylindrical core", ...
        "FontSize", cfg.diagnostics.title_font_size);
    export_figure(fig, fullfile(output_dir, "cylinder_core_profiles.png"), cfg);
end

function save_mode_amplitudes(result, output_dir, cfg)
    selected = select_modes_for_plot(result, cfg);
    if isempty(selected)
        return
    end

    fig = new_figure(cfg);
    ax = axes(fig);
    hold(ax, "on");
    colors = lines(numel(selected));
    for j = 1:numel(selected)
        g = selected(j);
        A = result.group_amplitude(:, g);
        semilogy(ax, result.time, max(A, realmin), ...
            "Color", colors(j, :), ...
            "LineWidth", cfg.diagnostics.line_width, ...
            "DisplayName", mode_label(result.group_m(g), result.group_n(g)));

        idx = result.fit_details(g).indices;
        if ~isempty(idx)
            detail = result.fit_details(g);
            fit_amplitude = exp(detail.intercept + detail.slope * result.time(idx));
            semilogy(ax, result.time(idx), fit_amplitude, "--", ...
                "Color", colors(j, :), ...
                "LineWidth", 1.2 * cfg.diagnostics.line_width, ...
                "HandleVisibility", "off");
        end
    end
    xlabel(ax, "Physical time");
    ylabel(ax, "Grouped amplitude");
    title(ax, "Cylindrical mode amplitudes");
    legend(ax, "Location", "best");
    style_axes(ax, cfg);
    export_figure(fig, fullfile(output_dir, "cylinder_mode_amplitudes.png"), cfg);
end

function save_growth_dispersion(result, output_dir, cfg)
    rates = result.growth_rates;
    valid = isfinite(rates.growth_rate_sigma) & isfinite(rates.mean_qR);
    if ~any(valid)
        return
    end

    fig = new_figure(cfg);
    ax = axes(fig);
    hold(ax, "on");
    unique_m = unique(rates.m(valid));
    colors = lines(numel(unique_m));
    if cfg.diagnostics.growth_x_axis == "n"
        x_values = rates.n;
        x_label = "Axial mode n";
    elseif cfg.diagnostics.growth_x_axis == "qR"
        x_values = rates.mean_qR;
        x_label = "Mean qR";
    else
        error("CylinderMode:InvalidDiagnosticAxis", ...
            "cfg.diagnostics.growth_x_axis must be 'n' or 'qR'.");
    end
    for k = 1:numel(unique_m)
        selected = valid & rates.m == unique_m(k);
        scatter(ax, x_values(selected), rates.growth_rate_sigma(selected), ...
            cfg.diagnostics.marker_size ^ 2, colors(k, :), "filled", ...
            "DisplayName", sprintf("m = %d", unique_m(k)));
    end
    yline(ax, 0, "k:", "HandleVisibility", "off");
    xlabel(ax, x_label);
    ylabel(ax, "Growth rate \sigma");
    title(ax, "Early cylindrical-mode growth rates");
    legend(ax, "Location", "best");
    style_axes(ax, cfg);
    export_figure(fig, fullfile(output_dir, "cylinder_growth_rates.png"), cfg);
end

function save_spectra(result, output_dir, cfg)
    if isempty(result.group_amplitude)
        return
    end
    first = find(result.frame_metrics.frame_valid, 1, "first");
    last = find(result.frame_metrics.frame_valid, 1, "last");
    if isempty(first) || isempty(last)
        return
    end

    fig = new_figure(cfg);
    ax = axes(fig);
    hold(ax, "on");
    for m = unique(result.group_m).'
        selected = result.group_m == m & result.group_n >= 0 ...
            & ~(result.group_m == 0 & result.group_n == 0);
        if ~any(selected)
            continue
        end
        q_first = result.group_qR(first, selected);
        q_last = result.group_qR(last, selected);
        semilogy(ax, q_first, max(result.group_amplitude(first, selected), realmin), ...
            "o-", "LineWidth", cfg.diagnostics.line_width, ...
            "DisplayName", sprintf("m=%d, initial", m));
        semilogy(ax, q_last, max(result.group_amplitude(last, selected), realmin), ...
            "s--", "LineWidth", cfg.diagnostics.line_width, ...
            "DisplayName", sprintf("m=%d, final", m));
    end
    xlabel(ax, "qR");
    ylabel(ax, "Grouped amplitude");
    title(ax, "Initial and final cylindrical spectra");
    legend(ax, "Location", "best");
    style_axes(ax, cfg);
    export_figure(fig, fullfile(output_dir, "cylinder_spectra.png"), cfg);
end

function save_individual_fits(result, output_dir, cfg)
    valid = find(arrayfun(@(fit) ~isempty(fit.indices), result.fit_details));
    valid = valid(~(result.group_m(valid) == 0 & result.group_n(valid) == 0));
    if isempty(valid)
        return
    end
    valid = valid(1:min(numel(valid), cfg.diagnostics.maximum_individual_fit_plots));
    fit_dir = fullfile(output_dir, "cylinder_mode_growth_fits");
    if ~isfolder(fit_dir)
        mkdir(fit_dir);
    end

    for g = valid(:).'
        detail = result.fit_details(g);
        A = result.group_amplitude(:, g);
        fig = new_figure(cfg);
        ax = axes(fig);
        semilogy(ax, result.time, max(A, realmin), "o-", ...
            "MarkerSize", cfg.diagnostics.marker_size, ...
            "LineWidth", cfg.diagnostics.line_width, ...
            "DisplayName", "Amplitude");
        hold(ax, "on");
        idx = detail.indices;
        fit_amplitude = exp(detail.intercept + detail.slope * result.time(idx));
        semilogy(ax, result.time(idx), fit_amplitude, "k--", ...
            "LineWidth", 1.3 * cfg.diagnostics.line_width, ...
            "DisplayName", sprintf("\\sigma = %.5g", detail.slope));
        xlabel(ax, "Physical time");
        ylabel(ax, "Grouped amplitude");
        title(ax, mode_label(result.group_m(g), result.group_n(g)));
        legend(ax, "Location", "best");
        style_axes(ax, cfg);
        filename = sprintf("cylinder_fit_m_%+d_n_%+d.png", ...
            result.group_m(g), result.group_n(g));
        export_figure(fig, fullfile(fit_dir, filename), cfg);
    end
end

function selected = select_modes_for_plot(result, cfg)
    nonbase = ~(result.group_m == 0 & result.group_n == 0);
    magnitude = max(result.group_amplitude, [], 1, "omitnan").';
    magnitude(~nonbase) = -inf;
    fitted = arrayfun(@(fit) ~isempty(fit.indices), result.fit_details);
    priority = magnitude;
    priority(fitted) = priority(fitted) + max(magnitude(isfinite(magnitude)), [], "omitnan");
    [~, order] = sort(priority, "descend");
    order = order(isfinite(priority(order)));
    selected = order(1:min(numel(order), cfg.diagnostics.maximum_modes_in_amplitude_plot));
end

function label = mode_label(m, n)
    label = sprintf("m = %d, n = %+d", m, n);
end

function fig = new_figure(cfg)
    fig = figure( ...
        "Visible", cfg.diagnostics.figure_visible, ...
        "Color", "w", ...
        "Position", cfg.diagnostics.figure_position);
end

function style_axes(ax, cfg)
    set(ax, ...
        "FontSize", cfg.diagnostics.axes_font_size, ...
        "LineWidth", 1.1, ...
        "Box", "on", ...
        "TickDir", "in");
    grid(ax, "on");
end

function xlabel_all(layout, label)
    axes_list = findall(layout.Parent, "Type", "axes");
    for ax = axes_list(:).'
        xlabel(ax, label);
    end
end

function export_figure(fig, filename, cfg)
    exportgraphics(fig, filename, "Resolution", cfg.diagnostics.resolution);
    close(fig);
end
