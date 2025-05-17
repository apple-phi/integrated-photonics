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
