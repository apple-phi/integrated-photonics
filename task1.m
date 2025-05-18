% Identify the modes supported by a
% symmetric silicon-on-insulator (SOI) planar waveguide with a core
% thickness of 300 nm, operating at 1550 nm. Use the cutoff condition for
% planar waveguides
lambda = 1550e-9; % wavelength in meters
n_Si = sb4.get_n('Si', lambda);
n_SiO2 = sb4.get_n('SiO2', lambda);
num_modes_Si = sb4.planar.find_num_modes(n_Si, n_SiO2, 300e-9, lambda);
disp(['Number of supported modes for Si+SiO2 at 300nm: ', num2str(num_modes_Si(1))]);
for i = 1:num_modes_Si(1)
    disp(['Mode ', num2str(i-1), ' is supported.']);
end

% Sweep the waveguide core thickness from 200 nm to 1 μm in steps of
% 100 nm. For each thickness, determine the number of supported modes,
% and plot the number of modes supported by the planar waveguide against
% its core thickness
thicknesses = 200e-9:100e-9:1e-6; % core thicknesses from 200 nm to 1 μm
num_modes_Si = sb4.planar.find_num_modes(n_Si, n_SiO2, thicknesses, lambda);

figure;
plot(thicknesses * 1e6, num_modes_Si, '-o');
xlabel('Core Thickness (μm)');
ylabel('Number of Supported Modes');
title('Supported Modes vs Core Thickness for SOI Waveguide at 1550 nm');
grid on;
legend('1550 nm');

% Repeat the above analysis using different waveguide core materials (e.g.,
% silicon nitride and lithium niobate), i.e., search for the refractive index
% data for these materials at 1550 nm, import the data into MATLAB, and
% plot the supported mode count – core thickness relationship for each
% material. Compare the results to that of the silicon waveguide and
% comment on how core material affects mode support.

n_SiN = sb4.get_n('Si3N4', lambda); % Silicon Nitride
n_LiNbO3 = sb4.get_n('LiNbO3', lambda); % Lithium Niobate
n_GaAs = sb4.get_n('GaAs', lambda); % Gallium Arsenide
num_modes_SiN = sb4.planar.find_num_modes(n_SiN, n_SiO2, thicknesses, lambda);
num_modes_LiNbO3 = sb4.planar.find_num_modes(n_LiNbO3, n_SiO2, thicknesses, lambda);
num_modes_GaAs = sb4.planar.find_num_modes(n_GaAs, n_SiO2, thicknesses, lambda);
figure;
plot(thicknesses * 1e6, num_modes_Si, '-o', ...
    thicknesses * 1e6, num_modes_SiN, '-o', ...
    thicknesses * 1e6, num_modes_LiNbO3, '-x', ...
    thicknesses * 1e6, num_modes_GaAs, '-s');
xlabel('Core Thickness (μm)');
ylabel('Number of Supported Modes');
title('Supported Modes vs Core Thickness for Different Materials at 1550 nm');
grid on;
legend('Si', 'Si3N4', 'LiNbO3', 'GaAs');
% The analysis shows that silicon nitride supports more modes than silicon
% at the same core thickness, while lithium niobate supports fewer modes.
% This is due to the different refractive index contrasts between the core
% and cladding materials. Silicon nitride has a lower refractive index than
% silicon, leading to a smaller effective index contrast, which allows for
% more modes to be supported. In contrast, lithium niobate has a higher
% refractive index than silicon, resulting in a larger effective index
% contrast, but its waveguide structure supports fewer modes due to its
% specific material properties and the way modes are confined in the
% waveguide.


% Repeat the analysis for a silicon waveguide at different wavelengths (e.g.,
% 1310 nm, 1600 nm, etc.). Consider material dispersion—note that the
% refractive indices of materials vary with wavelength. Again, plot "Number
% of Supported Modes vs. Core Thickness" for each wavelength. Compare
% and discuss how operating wavelength affects the number of supported
% modes.

wavelengths = [850e-9, 1310e-9, 1550e-9, 1600e-9, 2000e-9]; % wavelengths in meters
num_modes_all = zeros(length(thicknesses), length(wavelengths));
for j = 1:length(wavelengths)
    lambda = wavelengths(j);
    n_Si = sb4.get_n('Si', lambda);
    n_SiO2 = sb4.get_n('SiO2', lambda);
    num_modes_all(:, j) = sb4.planar.find_num_modes(n_Si, n_SiO2, thicknesses, lambda);
end
figure;
plot(thicknesses * 1e6, num_modes_all, '-o');
xlabel('Core Thickness (μm)');
ylabel('Number of Supported Modes');
title('Supported Modes vs Core Thickness for SOI Waveguide at Different Wavelengths');
grid on;
legend(arrayfun(@(x) sprintf('%.0f nm', x * 1e9), wavelengths, 'UniformOutput', false));


% Discuss how operating wavelength affects the number of supported modes
% The number of supported modes generally decreases with increasing wavelength
% due to the reduced effective index contrast between the core and cladding.
% Longer wavelengths result in a lower spatial frequency of the modes, which
% can lead to fewer modes being supported in the waveguide. This is particularly
% evident in silicon waveguides, where the refractive index contrast is significant.
% As the wavelength increases, the effective index of the modes approaches that of the cladding,
% leading to a reduction in the number of modes that can propagate. This trend is consistent across
% different core thicknesses, with thicker cores supporting more modes at a given wavelength.
% The analysis shows that for a silicon-on-insulator waveguide, the number of supported modes
% decreases as the wavelength increases, highlighting the importance of considering
% wavelength in the design of photonic devices.


% Consider the Si+SiO2 waveguide at 1550 nm.
% For each mode (separately), calculate and plot the effective index as a function of core thickness, starting from where the mode theoretically exists
% Finish continuously tracking each mode as it increases monotonically before you move on to the next mode.
% This will involve calculating the effective index for each mode at various core thicknesses
% Put all the code here, and use a loop to plot each mode on the same figure.
% Define any additional functions you need for the calculations (e.g., dispersion equations, and a function to calculate the effective index).

disp('Starting effective index calculation for Si/SiO2 waveguide at 1550 nm...');

lambda_neff = 1550e-9; % Wavelength in meters
n_core_neff = sb4.get_n('Si', lambda_neff);
n_clad_neff = sb4.get_n('SiO2', lambda_neff);

% Define a fine range of core thicknesses for plotting n_eff
% Start from a small thickness to observe mode cutoffs
thicknesses_plot_neff = linspace(10e-9, 2000e-9, 500); % 10 nm to 1000 nm, 200 points

figure; % New figure for n_eff plots
hold on;
xlabel('Core Thickness (\mum)');
ylabel('Effective Index (n_{eff})');
title(['Effective Index vs. Core Thickness for Si/SiO2 Waveguide at ', num2str(lambda_neff*1e9), ' nm']);
grid on;

% Add horizontal lines for n_clad and n_core
min_thick_um = min(thicknesses_plot_neff) * 1e6;
max_thick_um = max(thicknesses_plot_neff) * 1e6;
plot([min_thick_um max_thick_um], [n_clad_neff n_clad_neff], 'k:', 'LineWidth', 1, 'DisplayName', 'n_{SiO_2}');
plot([min_thick_um max_thick_um], [n_core_neff n_core_neff], 'k--', 'LineWidth', 1, 'DisplayName', 'n_{Si}');

% Get some distinct colors for plotting modes
num_colors_needed = 10; % Default, adjust if many modes are expected
plot_colors = lines(num_colors_needed);

% Determine maximum mode order 'm' to consider for TE and TM polarizations
% This is based on the thickest waveguide in our range
max_modes_at_max_thickness = sb4.planar.find_num_modes(n_core_neff, n_clad_neff, max(thicknesses_plot_neff), lambda_neff);

% Initialize max mode orders to -1 (meaning no modes of that type if count is 0 or array is malformed)
max_m_TE_neff = -1;
max_m_TM_neff = -1;

if ~isempty(max_modes_at_max_thickness)
    % Check for TE modes
    if max_modes_at_max_thickness(1) > 0
        max_m_TE_neff = max_modes_at_max_thickness(1) - 1; % Max mode order m for TE (0-indexed)
    end

    % Check for TM modes, ensuring the second element exists
    if length(max_modes_at_max_thickness) > 1 && max_modes_at_max_thickness(2) > 0
        max_m_TM_neff = max_modes_at_max_thickness(2) - 1; % Max mode order m for TM (0-indexed)
    end
end

% --- TE Modes ---
disp('Calculating and plotting TE modes effective indices...');
for m = 0:max_m_TE_neff
    current_neff_values_TE = NaN(1, length(thicknesses_plot_neff));

    for i = 1:length(thicknesses_plot_neff)
        d_core = thicknesses_plot_neff(i);
        neff = sb4.planar.find_neff(n_core_neff, n_clad_neff, d_core, lambda_neff, 'TE', m);

        if ~isnan(neff) && isreal(neff) && neff > n_clad_neff && neff < n_core_neff
            if i > 1 && ~isnan(current_neff_values_TE(i-1)) && neff < (current_neff_values_TE(i-1) - 1e-5) % Tolerance for numerical noise
                fprintf('Warning: TE mode %d at thickness %.3f nm, n_eff %.4f might violate monotonicity (prev_neff %.4f).\n', ...
                    m, d_core*1e9, neff, current_neff_values_TE(i-1));
            end
            current_neff_values_TE(i) = neff;
        end
    end

    valid_plot_points_TE = ~isnan(current_neff_values_TE) & current_neff_values_TE > n_clad_neff;
    if any(valid_plot_points_TE)
        first_valid_idx_TE = find(valid_plot_points_TE, 1, 'first');
        plot(thicknesses_plot_neff(first_valid_idx_TE:end) * 1e6, ...
            current_neff_values_TE(first_valid_idx_TE:end), ...
            'Color', plot_colors(mod(m, num_colors_needed)+1,:), ...
            'LineWidth', 1.5, ...
            'DisplayName', sprintf('TE_%d', m));
    end
end

% --- TM Modes ---
disp('Calculating and plotting TM modes effective indices...');
for m = 0:max_m_TM_neff
    current_neff_values_TM = NaN(1, length(thicknesses_plot_neff));

    for i = 1:length(thicknesses_plot_neff)
        d_core = thicknesses_plot_neff(i);
        neff = sb4.planar.find_neff(n_core_neff, n_clad_neff, d_core, lambda_neff, 'TM', m);

        if ~isnan(neff) && isreal(neff) && neff > n_clad_neff && neff < n_core_neff
            if i > 1 && ~isnan(current_neff_values_TM(i-1)) && neff < (current_neff_values_TM(i-1) - 1e-5)
                fprintf('Warning: TM mode %d at thickness %.3f nm, n_eff %.4f might violate monotonicity (prev_neff %.4f).\n', ...
                    m, d_core*1e9, neff, current_neff_values_TM(i-1));
            end
            current_neff_values_TM(i) = neff;
        end
    end

    valid_plot_points_TM = ~isnan(current_neff_values_TM) & current_neff_values_TM > n_clad_neff;
    if any(valid_plot_points_TM)
        first_valid_idx_TM = find(valid_plot_points_TM, 1, 'first');
        plot(thicknesses_plot_neff(first_valid_idx_TM:end) * 1e6, ...
            current_neff_values_TM(first_valid_idx_TM:end), ...
            '--', 'Color', plot_colors(mod(m, num_colors_needed)+1,:), ... % Dashed line for TM modes
            'LineWidth', 1.5, ...
            'DisplayName', sprintf('TM_%d', m));
    end
end

% Finalize plot
ylim([n_clad_neff - 0.05, n_core_neff + 0.05]); % Set y-axis limits for better visualization
legend('show', 'Location', 'southeast'); % Display legend
hold off;

disp('Effective index plotting section complete.');

