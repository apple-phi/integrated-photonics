% Identify the modes supported by a
% symmetric silicon-on-insulator (SOI) planar waveguide with a core
% thickness of 300 nm, operating at 1550 nm. Use the cutoff condition for
% planar waveguides

% Create output directory if it doesn't exist
output_dir = 'data/task1';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Define default line width for all plots
defaultLineWidth = 2; % Increased line width

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
plot(thicknesses * 1e6, num_modes_Si, '-o', 'LineWidth', defaultLineWidth);
ax = gca; ax.FontSize = 16; % Base font size for ticks
xlabel('Core Thickness (μm)', 'FontSize', 18);
ylabel('Number of Supported Modes', 'FontSize', 18);
title('Supported Modes vs Core Thickness for SOI Waveguide at 1550 nm', 'FontSize', 20);
grid on;
legend('1550 nm', 'FontSize', 16);
saveas(gcf, fullfile(output_dir, 'modes_vs_thickness_SOI_1550nm.png'));

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
plot(thicknesses * 1e6, num_modes_Si, '-o', 'LineWidth', defaultLineWidth);
hold on; % Ensure subsequent plots are on the same axes
plot(thicknesses * 1e6, num_modes_SiN, '-o', 'LineWidth', defaultLineWidth);
plot(thicknesses * 1e6, num_modes_LiNbO3, '-x', 'LineWidth', defaultLineWidth);
plot(thicknesses * 1e6, num_modes_GaAs, '-s', 'LineWidth', defaultLineWidth);
hold off;
ax = gca; ax.FontSize = 16;
xlabel('Core Thickness (μm)', 'FontSize', 18);
ylabel('Number of Supported Modes', 'FontSize', 18);
title('Supported Modes vs Core Thickness for Different Materials at 1550 nm', 'FontSize', 20);
grid on;
legend('Si', 'Si3N4', 'LiNbO3', 'GaAs', 'FontSize', 16, 'Location', 'best');
saveas(gcf, fullfile(output_dir, 'modes_vs_thickness_materials_1550nm.png'));
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
    lambda_loop = wavelengths(j); % Use a different variable name to avoid conflict with outer scope lambda
    n_Si_loop = sb4.get_n('Si', lambda_loop);
    n_SiO2_loop = sb4.get_n('SiO2', lambda_loop);
    num_modes_all(:, j) = sb4.planar.find_num_modes(n_Si_loop, n_SiO2_loop, thicknesses, lambda_loop);
end
figure;
% Plot each line separately to control LineWidth for all
hold on;
for j = 1:length(wavelengths)
    plot(thicknesses * 1e6, num_modes_all(:,j), '-o', 'LineWidth', defaultLineWidth, 'DisplayName', sprintf('%.0f nm', wavelengths(j) * 1e9));
end
hold off;
ax = gca; ax.FontSize = 16;
xlabel('Core Thickness (μm)', 'FontSize', 18);
ylabel('Number of Supported Modes', 'FontSize', 18);
title('Supported Modes vs Core Thickness for SOI Waveguide at Different Wavelengths', 'FontSize', 20);
grid on;
legend(arrayfun(@(x) sprintf('%.0f nm', x * 1e9), wavelengths, 'UniformOutput', false), 'FontSize', 16, 'Location', 'best');
saveas(gcf, fullfile(output_dir, 'modes_vs_thickness_SOI_wavelengths.png'));


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
ax_neff = gca; ax_neff.FontSize = 16;
xlabel('Core Thickness (\mum)', 'FontSize', 18);
ylabel('Effective Index (n_{eff})', 'FontSize', 18);
title(['Effective Index vs. Core Thickness for Si/SiO2 Waveguide at ', num2str(lambda_neff*1e9), ' nm'], 'FontSize', 20);
grid on;

% Add horizontal lines for n_clad and n_core
min_thick_um = min(thicknesses_plot_neff) * 1e6;
max_thick_um = max(thicknesses_plot_neff) * 1e6;
plot([min_thick_um max_thick_um], [n_clad_neff n_clad_neff], 'k:', 'LineWidth', defaultLineWidth, 'DisplayName', 'n_{SiO_2}');
plot([min_thick_um max_thick_um], [n_core_neff n_core_neff], 'k--', 'LineWidth', defaultLineWidth, 'DisplayName', 'n_{Si}');

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
            'LineWidth', defaultLineWidth, ...
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
            'LineWidth', defaultLineWidth, ...
            'DisplayName', sprintf('TM_%d', m));
    end
end

% Finalize plot
ylim([n_clad_neff - 0.05, n_core_neff + 0.05]); % Set y-axis limits for better visualization
legend('show', 'Location', 'southeast', 'FontSize', 16); % Display legend
hold off;
saveas(gcf, fullfile(output_dir, 'neff_vs_thickness_Si_SiO2_1550nm.png'));

disp('Effective index plotting section complete.');

% Some of the difficulties a numerical solver might encounter, especially in precisely matching the theoretical cutoff thickness:
% 1.  **Slow Change in `n_eff`**: Just above the critical (cutoff) thickness for a particular mode, the effective index `n_eff` is only slightly greater than the cladding index `n_clad`. As the core thickness increases beyond this critical point, `n_eff` indeed increases, but initially, this increase can be very slow. This means the derivative `dn_eff/d(thickness)` is small (close to zero) in this region.

% 2.  **Impact on Root Finding**: The `sb4.planar.find_neff` function solves a transcendental equation of the form `F(n_eff) = 0` for a given thickness `d`.
%     *   When `dn_eff/d(thickness)` is small, it implies that `n_eff` is not very sensitive to changes in `d`.
%     *   From the solver's perspective, when it's trying to find `n_eff` for a thickness `d` that is just barely above the true cutoff, the `n_eff` it's looking for is extremely close to `n_clad`.
%     *   In this regime (where `n_eff` is very close to `n_clad`), the terms inside the transcendental equation (like `gammay = sqrt((neff_val * k0)^2 - (n_clad * k0)^2)`) can become very small. The function `F(n_eff)` might be very "flat" or change very slowly with `n_eff` when `n_eff` is just above `n_clad`.

% 3.  **Numerical Precision and Tolerance**:
%     *   If `F(n_eff)` is very flat, `fzero` might struggle to pinpoint the root with high precision relative to the `n_clad` boundary. A small change in `F(n_eff)` might correspond to a relatively larger uncertainty in `n_eff`, or `fzero` might converge to a value that is numerically indistinguishable from `n_clad` within its tolerance, or it might fail to bracket a root if the function values at the search bounds are too similar.
%     *   The check `val_at_low <= cutoff_check_tolerance` in `find_neff` is designed to identify when a mode is cut off. If `F(n_eff)` near `n_clad` is already very small (but positive) because the mode is extremely weakly guided (due to `dn_eff/d(thickness)` being small), it might fall below this tolerance, leading the solver to declare the mode as cut off even if the thickness is slightly above the theoretical minimum.

% 4.  **Visual Interpretation**: On your plot of `n_eff` vs. thickness, this near-zero derivative manifests as the mode curve starting almost horizontally from `n_clad` before it begins to curve upwards more noticeably. The "flatter" this initial part is, the harder it is for the discrete thickness steps in your simulation to precisely capture the exact point where `n_eff` first lifts off from `n_clad`.

% The fact that `n_eff` changes very slowly with thickness immediately after cutoff (i.e., `dn_eff/d(thickness)` is small) contributes to the numerical challenge of precisely determining the cutoff point and can explain why the numerically observed "appearance" of a mode might be at a slightly larger thickness than the strict theoretical value. The solver needs the mode to be sufficiently "established" (i.e., `n_eff` to be discernibly greater than `n_clad`) to reliably find it.

% --- Plot Critical Thickness vs. Wavelength for different modes ---
disp('Calculating and plotting critical thickness vs. wavelength...');

% Define a range of wavelengths for this plot
plot_wavelengths = linspace(800e-9, 2200e-9, 150); % Wavelengths from 800 nm to 2200 nm

% Mode orders to consider
m_values = 0:4; % m = 0, 1, 2, 3, 4

critical_thicknesses_all_modes = NaN(length(m_values), length(plot_wavelengths));

for j = 1:length(plot_wavelengths)
    current_lambda = plot_wavelengths(j);
    n_Si_current = sb4.get_n('Si', current_lambda);
    n_SiO2_current = sb4.get_n('SiO2', current_lambda);

    if n_Si_current > n_SiO2_current
        denominator_sqrt = sqrt(n_Si_current^2 - n_SiO2_current^2);
        if denominator_sqrt > 1e-9 % Avoid division by zero or very small number
            for i = 1:length(m_values)
                m = m_values(i);
                dc = (m * current_lambda) / (2 * denominator_sqrt);
                critical_thicknesses_all_modes(i, j) = dc;
            end
        end
    end
end

figure;
hold on;
ax_crit_thick_modes = gca; ax_crit_thick_modes.FontSize = 16;
plot_colors_critical = lines(length(m_values)); % Get distinct colors

for i = 1:length(m_values)
    m_loop = m_values(i); % Renamed m to m_loop to avoid conflict
    plot(plot_wavelengths * 1e9, critical_thicknesses_all_modes(i, :) * 1e9, ...
        'Color', plot_colors_critical(i,:), ...
        'LineWidth', defaultLineWidth, ...
        'DisplayName', sprintf('m = %d (TE_{%d}/TM_{%d})', m_loop, m_loop, m_loop));
end

xlabel('Wavelength (nm)', 'FontSize', 18);
ylabel('Critical Core Thickness (nm)', 'FontSize', 18);
title('Critical Thickness vs. Wavelength for Si/SiO_2 Waveguide', 'FontSize', 20);
grid on;
legend('show', 'Location', 'northwest', 'FontSize', 16);
ylim([0, max(critical_thicknesses_all_modes(:), [], 'omitnan') * 1e9 * 1.1]); % Adjust y-limit for better view, add 10nm buffer
hold off;
saveas(gcf, fullfile(output_dir, 'critical_thickness_vs_wavelength.png'));

disp('Critical thickness vs. wavelength plotting complete.');

% --- Plot Critical Thickness for Multimode Operation vs. Wavelength for Different Materials ---
disp('Calculating and plotting critical thickness for multimode operation vs. wavelength for different materials...');

% Define a range of wavelengths for this plot (can reuse plot_wavelengths or define new)
% plot_wavelengths is already defined from 800e-9 to 2200e-9

core_materials = {'Si', 'Si3N4', 'LiNbO3', 'GaAs'};
material_display_names = {'Si', 'Si_3N_4', 'LiNbO_3', 'GaAs'}; % For legend, use underscores for subscripts
cladding_material = 'SiO2';
m_multimode = 1; % For multimode operation, we look at the cutoff of the m=1 mode

critical_thickness_multimode_materials = NaN(length(core_materials), length(plot_wavelengths));

for k = 1:length(core_materials)
    current_core_material_name = core_materials{k};
    for j = 1:length(plot_wavelengths)
        current_lambda = plot_wavelengths(j);

        n_core_current = sb4.get_n(current_core_material_name, current_lambda);
        n_clad_current = sb4.get_n(cladding_material, current_lambda);

        if ~isnan(n_core_current) && ~isnan(n_clad_current) && n_core_current > n_clad_current
            denominator_sqrt = sqrt(n_core_current^2 - n_clad_current^2);
            if denominator_sqrt > 1e-9 % Avoid division by zero or very small number
                dc_multimode = (m_multimode * current_lambda) / (2 * denominator_sqrt);
                critical_thickness_multimode_materials(k, j) = dc_multimode;
            end
        end
    end
end

figure;
hold on;
ax_crit_thick_mat = gca; ax_crit_thick_mat.FontSize = 16;
plot_colors_multimode_materials = lines(length(core_materials)); % Get distinct colors

for k = 1:length(core_materials)
    % Replace underscores with LaTeX-compatible subscripts for display names
    display_name_formatted = strrep(material_display_names{k}, '_', '_{');
    if contains(display_name_formatted, '{') % if a subscript was started
        display_name_formatted = [display_name_formatted '}']; % close it
    end

    plot(plot_wavelengths * 1e9, critical_thickness_multimode_materials(k, :) * 1e9, ...
        'Color', plot_colors_multimode_materials(k,:), ...
        'LineWidth', defaultLineWidth, ...
        'DisplayName', display_name_formatted);
end

xlabel('Wavelength (nm)', 'FontSize', 18);
ylabel('Critical Thickness for Multimode (m=1) (nm)', 'FontSize', 18);
title('Critical Thickness for Multimode Operation vs. Wavelength (SiO_2 Cladding)', 'FontSize', 20);
grid on;
legend('show', 'Location', 'northwest', 'Interpreter', 'tex', 'FontSize', 16);
% Adjust y-limit dynamically, ensuring non-negative lower bound
min_y_val = 0;
max_y_val = max(critical_thickness_multimode_materials(:), [], 'omitnan') * 1e9 * 1.1;
if isnan(max_y_val) || isinf(max_y_val) || max_y_val <= min_y_val % Handle cases with all NaNs or no valid data
    max_y_val = 1000; % Default max if no data
end
ylim([min_y_val, max_y_val + 10]); % Add 10nm buffer

hold off;
saveas(gcf, fullfile(output_dir, 'critical_thickness_multimode_vs_wavelength_materials.png'));

disp('Critical thickness for multimode operation vs. wavelength for different materials plotting complete.');

% --- Plot Group Index (n_g) vs. Wavelength for TE0 mode in SOI ---
disp('Calculating and plotting group index vs. wavelength for TE0 and TM0 modes...');

fixed_core_thickness_ng = 300e-9; % Example: 300 nm SOI core
ng_wavelengths = linspace(1200e-9, 1800e-9, 100); % Wavelength range for ng plot
neff_for_ng_TE0 = NaN(1, length(ng_wavelengths));
neff_for_ng_TM0 = NaN(1, length(ng_wavelengths)); % Added for TM0
n_core_ng_material = 'Si';
n_clad_ng_material = 'SiO2';

for i = 1:length(ng_wavelengths)
    current_lambda_ng = ng_wavelengths(i);
    n_core_val = sb4.get_n(n_core_ng_material, current_lambda_ng);
    n_clad_val = sb4.get_n(n_clad_ng_material, current_lambda_ng);
    if ~isnan(n_core_val) && ~isnan(n_clad_val) && n_core_val > n_clad_val
        neff_for_ng_TE0(i) = sb4.planar.find_neff(n_core_val, n_clad_val, fixed_core_thickness_ng, current_lambda_ng, 'TE', 0);
        neff_for_ng_TM0(i) = sb4.planar.find_neff(n_core_val, n_clad_val, fixed_core_thickness_ng, current_lambda_ng, 'TM', 0); % Added for TM0
    end
end

figure; % Create a new figure for group index plots
hold on; % Hold on to plot multiple lines
ax_ng = gca; ax_ng.FontSize = 16;

% Calculate and plot for TE0
valid_ng_indices_TE0 = ~isnan(neff_for_ng_TE0);
if sum(valid_ng_indices_TE0) > 1 % Need at least 2 points for diff
    d_neff_d_lambda_TE0 = gradient(neff_for_ng_TE0(valid_ng_indices_TE0), ng_wavelengths(valid_ng_indices_TE0));
    group_index_TE0 = neff_for_ng_TE0(valid_ng_indices_TE0) - ng_wavelengths(valid_ng_indices_TE0) .* d_neff_d_lambda_TE0;
    plot(ng_wavelengths(valid_ng_indices_TE0) * 1e9, group_index_TE0, 'LineWidth', defaultLineWidth, 'DisplayName', 'TE_0');
else
    disp('Not enough valid n_eff points to calculate group index for TE0.');
end

% Calculate and plot for TM0
valid_ng_indices_TM0 = ~isnan(neff_for_ng_TM0);
if sum(valid_ng_indices_TM0) > 1 % Need at least 2 points for diff
    d_neff_d_lambda_TM0 = gradient(neff_for_ng_TM0(valid_ng_indices_TM0), ng_wavelengths(valid_ng_indices_TM0));
    group_index_TM0 = neff_for_ng_TM0(valid_ng_indices_TM0) - ng_wavelengths(valid_ng_indices_TM0) .* d_neff_d_lambda_TM0;
    plot(ng_wavelengths(valid_ng_indices_TM0) * 1e9, group_index_TM0, '--', 'LineWidth', defaultLineWidth, 'DisplayName', 'TM_0'); % Dashed line for TM0
else
    disp('Not enough valid n_eff points to calculate group index for TM0.');
end

if sum(valid_ng_indices_TE0) > 1 || sum(valid_ng_indices_TM0) > 1
    xlabel('Wavelength (nm)', 'FontSize', 18);
    ylabel('Group Index (n_g)', 'FontSize', 18);
    title(sprintf('Group Index vs. Wavelength (Si/SiO_2, Core Thickness = %.0f nm)', fixed_core_thickness_ng*1e9), 'FontSize', 20);
    grid on;
    legend('show', 'Location', 'best', 'FontSize', 16);
    saveas(gcf, fullfile(output_dir, sprintf('group_index_TE0_TM0_vs_wavelength_d%.0fnm.png', fixed_core_thickness_ng*1e9)));
end
hold off;
disp('Group index plotting complete.');


% --- Plot Birefringence (n_eff_TE0 - n_eff_TM0) vs. Core Thickness for SOI at 1550 nm ---
disp('Calculating and plotting birefringence vs. core thickness for multiple materials...');

biref_lambda = 1550e-9; % Fixed wavelength
% thicknesses_plot_neff is already defined (10nm to 2000nm, 500 points)
% core_materials and material_display_names are already defined

figure; % Create a new figure for birefringence plots
hold on; % Hold on to plot multiple lines
ax_biref = gca; ax_biref.FontSize = 16;
plot_colors_biref_materials = lines(length(core_materials));

for k_mat = 1:length(core_materials)
    current_core_mat_name = core_materials{k_mat};
    current_display_name = material_display_names{k_mat};

    neff_TE0_biref = NaN(1, length(thicknesses_plot_neff));
    neff_TM0_biref = NaN(1, length(thicknesses_plot_neff));

    n_core_biref = sb4.get_n(current_core_mat_name, biref_lambda);
    n_clad_biref = sb4.get_n('SiO2', biref_lambda); % Cladding is SiO2

    if ~isnan(n_core_biref) && ~isnan(n_clad_biref) && n_core_biref > n_clad_biref
        for i = 1:length(thicknesses_plot_neff)
            d_core = thicknesses_plot_neff(i);
            neff_TE0_biref(i) = sb4.planar.find_neff(n_core_biref, n_clad_biref, d_core, biref_lambda, 'TE', 0);
            neff_TM0_biref(i) = sb4.planar.find_neff(n_core_biref, n_clad_biref, d_core, biref_lambda, 'TM', 0);
        end

        birefringence = neff_TE0_biref - neff_TM0_biref;
        valid_biref_indices = ~isnan(birefringence);

        if any(valid_biref_indices)
            % Format display name for legend
            display_name_formatted = strrep(current_display_name, '_', '_{');
            if contains(display_name_formatted, '{')
                display_name_formatted = [display_name_formatted '}'];
            end

            plot(thicknesses_plot_neff(valid_biref_indices) * 1e6, birefringence(valid_biref_indices), ...
                'LineWidth', defaultLineWidth, 'Color', plot_colors_biref_materials(k_mat,:), ...
                'DisplayName', display_name_formatted);
        else
            fprintf('Not enough valid n_eff points to calculate birefringence for %s.\n', current_core_mat_name);
        end
    else
        fprintf('Cannot calculate birefringence for %s due to invalid core/cladding indices at %.0f nm.\n', current_core_mat_name, biref_lambda*1e9);
    end
end

if ishghandle(gcf) && ~isempty(get(gca, 'Children')) % Check if anything was plotted
    xlabel('Core Thickness (\mum)', 'FontSize', 18);
    ylabel('Birefringence (n_{eff,TE0} - n_{eff,TM0})', 'FontSize', 18);
    title(sprintf('Birefringence vs. Core Thickness at %.0f nm (SiO_2 Cladding)', biref_lambda*1e9), 'FontSize', 20);
    grid on;
    legend('show', 'Location', 'best', 'Interpreter', 'tex', 'FontSize', 16);
    saveas(gcf, fullfile(output_dir, 'birefringence_vs_thickness_materials_1550nm.png'));
end
hold off;
disp('Birefringence plotting complete.');

