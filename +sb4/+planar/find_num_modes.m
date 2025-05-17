function num_modes = find_num_modes(n_core, n_clad, thicknesses, lambda)
% FIND_MODES Calculate the number of supported modes for given core and cladding
% refractive indices, core thicknesses, and wavelength.
%
% Parameters:
%   n_core - Refractive index of the core
%   n_clad - Refractive index of the cladding
%   thicknesses - Array of core thicknesses (in meters)
%   lambda - Wavelength (in meters)
%
% Returns:
%   num_modes - Array of supported modes for each thickness

num_modes = zeros(size(thicknesses)); % preallocate array for number of modes
for i = 1:length(thicknesses)
    t_g = thicknesses(i); % current core thickness
    LHS = 2 * t_g * sqrt(n_core^2 - n_clad^2); % left-hand side of cutoff condition
    num_modes(i) = sum(LHS > (0:5) * lambda); % count supported modes
end
end
