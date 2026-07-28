function cfg = cylinder_mode_defaults()
%CYLINDER_MODE_DEFAULTS Default settings for cylindrical mode analysis.

this_dir = fileparts(mfilename("fullpath"));

cfg.data_root = fullfile(this_dir, "data", "fs_batch_data");
cfg.output_root = fullfile(this_dir, "data", "cylinder_mode_analysis");
cfg.continue_on_error = true;
cfg.verbose = true;

cfg.frame.first_index = 0;
cfg.frame.last_index = inf;
cfg.frame.stride = 1;
cfg.frame.max_frames = inf;
cfg.frame.time_field_names = ["time", "analysis_time", "physical_time"];

cfg.alignment.axis_mode = "tracked_pca"; % "tracked_pca", "fixed_initial", or "known".
cfg.alignment.known_axis = [];
cfg.alignment.maximum_axis_jump_degrees = 5;
cfg.alignment.smooth_axis_fraction = 0;

cfg.core.method = "auto"; % "auto", "persistent", or "detect_each_frame".
cfg.core.axial_bins = 100;
cfg.core.profile_smoothing_bins = 5;
cfg.core.maximum_abs_radius_slope = 0.20;
cfg.core.maximum_relative_cross_section_std = 0.12;
cfg.core.minimum_radius_fraction = 0.70;
cfg.core.transition_margin_in_radii = 0.75;
cfg.core.minimum_core_length_in_radii = 4.0;
cfg.core.minimum_vertices = 50;
cfg.core.robust_endpoint_quantile = 0.01;
cfg.core.fallback_half_length_fraction = 0.28;
cfg.core.maximum_relative_axial_fraction_change = 0.20;

cfg.window.type = "tukey"; % "tukey" or "none".
cfg.window.alpha = 0.25;

cfg.modes.m_max = 0;
cfg.modes.n_min = 1;
cfg.modes.n_max = 16;
cfg.modes.include_n_zero_nonaxisymmetric = true;
cfg.modes.projection = "weighted_lstsq"; % "weighted_lstsq" or "quadrature_nudft".
cfg.modes.regularization = 1e-12;
cfg.modes.maximum_condition_number = 1e12;

cfg.growth.mode = "auto"; % "auto" or "fixed".
cfg.growth.fixed_start_time = -inf;
cfg.growth.fixed_end_time = inf;
cfg.growth.minimum_points = 6;
cfg.growth.maximum_points = 40;
cfg.growth.maximum_early_frames = 100;
cfg.growth.minimum_amplitude_ratio = 1.20;
cfg.growth.minimum_signal_to_noise = 3.0;
cfg.growth.maximum_dimensionless_amplitude = 0.05;
cfg.growth.minimum_r_squared = 0.90;
cfg.growth.maximum_relative_qR_change = 0.15;
cfg.growth.maximum_endpoint_slope_change = 0.35;
cfg.growth.absolute_noise_floor = 1e-13;
cfg.growth.accept_poor_fits = true;

cfg.quality.maximum_area_relative_drift = 0.02;
cfg.quality.maximum_volume_relative_drift = 0.02;
cfg.quality.maximum_area_relative_jump = 0.05;
cfg.quality.maximum_volume_relative_jump = 0.05;
cfg.quality.maximum_reconstruction_error = 0.35;
cfg.quality.minimum_core_area_fraction = 0.10;
cfg.quality.minimum_effective_vertices = 30;
cfg.quality.minimum_axial_samples_per_wavelength = 10;
cfg.quality.minimum_azimuthal_samples_per_wavelength = 8;

cfg.pairing.enabled = false;
cfg.pairing.roles = strings(0, 1);     % "ordinary", "reference", "plus", or "minus".
cfg.pairing.group_ids = strings(0, 1);
cfg.pairing.epsilon = [];
cfg.pairing.time_tolerance = 1e-8;

cfg.output.save_mat = true;
cfg.output.save_frame_csv = true;
cfg.output.save_mode_csv = true;
cfg.output.save_growth_csv = true;
cfg.output.save_resolved_json = true;
cfg.output.overwrite = true;

cfg.diagnostics.enabled = true;
cfg.diagnostics.figure_visible = "off";
cfg.diagnostics.figure_position = [100, 100, 1200, 760];
cfg.diagnostics.resolution = 220;
cfg.diagnostics.maximum_modes_in_amplitude_plot = 12;
cfg.diagnostics.save_individual_fit_plots = true;
cfg.diagnostics.maximum_individual_fit_plots = 16;
cfg.diagnostics.growth_x_axis = "n"; % "n" or "qR".
cfg.diagnostics.line_width = 1.7;
cfg.diagnostics.marker_size = 7;
cfg.diagnostics.axes_font_size = 13;
cfg.diagnostics.label_font_size = 15;
cfg.diagnostics.title_font_size = 16;
end
