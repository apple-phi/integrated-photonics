function n = get_n(material, lambda)
% GET_N Get refractive index n for a material at wavelength lambda (in meters) using CSV data
% The CSV files are expected to be in the 'materials' folder and have columns: wl (um), n
%
% Parameters:
%   material - Material name as a string (e.g., 'Si', 'SiO2')
%   lambda - Wavelength in meters
%
% Returns:
%   n - The interpolated refractive index at the specified wavelength

% Discover all CSV files in the materials folder
csv_files = dir(fullfile('materials', '*.csv'));
material_found = false;

for k = 1:length(csv_files)
    fname = csv_files(k).name;
    % Try to infer material from filename (case-insensitive, ignore spaces and dashes)
    fname_clean = lower(regexprep(fname, '[\s\-_.]', ''));
    material_clean = lower(regexprep(material, '[\s\-_.]', ''));
    if contains(fname_clean, material_clean)
        filename = fullfile('materials', fname);
        material_found = true;
        break;
    end
end

if ~material_found
    error(['Material file for "', material, '" not found in materials folder.']);
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
disp(['Using refractive index for ', material, ' at ', num2str(lambda * 1e9), ' nm: ', num2str(n)]);
end
