% Identify the modes supported by a
% symmetric silicon-on-insulator (SOI) planar waveguide with a core
% thickness of 300 nm, operating at 1550 nm. Use the cutoff condition for
% planar waveguides
lambda = 1550e-9; % wavelength in meters
% Use CSV-based interpolation for refractive indices
n_Si = get_n('Si', lambda);
n_SiO2 = get_n('SiO2', lambda);
disp(['Using refractive index for Si at ', num2str(lambda*1e9), ' nm: ', num2str(n_Si)]);
disp(['Using refractive index for SiO2 at ', num2str(lambda*1e9), ' nm: ', num2str(n_SiO2)]);

num_modes = find_modes(n_Si, n_SiO2, 300e-9, lambda);
for i = 1:num_modes(1)
    disp(['Mode ', num2str(i-1), ' is supported.']);
end

% Sweep the waveguide core thickness from 200 nm to 1 μm in steps of
% 100 nm. For each thickness, determine the number of supported modes,
% and plot the number of modes supported by the planar waveguide against
% its core thickness
thicknesses = 200e-9:100e-9:1e-6; % core thicknesses from 200 nm to 1 μm
num_modes = find_modes(n_Si, n_SiO2, thicknesses, lambda);

% Plot the results and make it pretty
figure;
plot(thicknesses * 1e6, num_modes, '-o');
xlabel('Core Thickness (μm)');
ylabel('Number of Supported Modes');
title('Supported Modes vs Core Thickness for SOI Waveguide at 1550 nm');
grid on;
legend('1550 nm');

% Repeat the analysis for a silicon waveguide at different wavelengths (e.g.,
% 1310 nm, 1600 nm, etc.). Consider material dispersion—note that the
% refractive indices of materials vary with wavelength. Again, plot "Number
% of Supported Modes vs. Core Thickness" for each wavelength. Compare
% and discuss how operating wavelength affects the number of supported
% modes.

wavelengths = [1310e-9, 1550e-9, 1600e-9, 1700e-9]; % wavelengths in meters
num_modes_all = zeros(length(thicknesses), length(wavelengths));
for j = 1:length(wavelengths)
    lambda = wavelengths(j);
    n_Si = get_n('Si', lambda);
    n_SiO2 = get_n('SiO2', lambda);
    disp(['Using refractive index for Si at ', num2str(lambda*1e9), ' nm: ', num2str(n_Si)]);
    disp(['Using refractive index for SiO2 at ', num2str(lambda*1e9), ' nm: ', num2str(n_SiO2)]);
    num_modes_all(:, j) = find_modes(n_Si, n_SiO2, thicknesses, lambda);
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

function num_modes = find_modes(n_core, n_clad, thicknesses, lambda)
% Calculate the number of supported modes for given core and cladding
% refractive indices, core thicknesses, and wavelength.
num_modes = zeros(size(thicknesses)); % preallocate array for number of modes
for i = 1:length(thicknesses)
    t_g = thicknesses(i); % current core thickness
    LHS = 2 * t_g * sqrt(n_core^2 - n_clad^2); % left-hand side of cutoff condition
    num_modes(i) = sum(LHS > (0:5) * lambda); % count supported modes
end
end

function n = get_n(material, lambda)
% Get refractive index n for a material at wavelength lambda (in meters) using CSV data
% The CSV files are expected to be in the 'data' folder and have columns: wl (um), n

% Discover all CSV files in the data folder
csv_files = dir(fullfile('data', '*.csv'));
material_found = false;

for k = 1:length(csv_files)
    fname = csv_files(k).name;
    % Try to infer material from filename (case-insensitive, ignore spaces and dashes)
    fname_clean = lower(regexprep(fname, '[\s\-_.]', ''));
    material_clean = lower(regexprep(material, '[\s\-_.]', ''));
    if contains(fname_clean, material_clean)
        filename = fullfile('data', fname);
        material_found = true;
        break;
    end
end

if ~material_found
    error(['Material file for "', material, '" not found in data folder.']);
end

T = readtable(filename);
wl_um = T.wl; % wavelength in microns
n_vals = T.n;
lambda_um = lambda * 1e6; % convert input wavelength to microns

% Remove any non-finite or NaN values
valid = isfinite(wl_um) & isfinite(n_vals);
wl_um = wl_um(valid);
n_vals = n_vals(valid);

% Ensure sample points are unique for interp1
[wl_um_unique, idx_unique] = unique(wl_um);
n_vals_unique = n_vals(idx_unique);

if isempty(wl_um_unique) || isempty(n_vals_unique)
    error(['No valid data in file: ', filename]);
end

n = interp1(wl_um_unique, n_vals_unique, lambda_um, 'linear', 'extrap');
end
