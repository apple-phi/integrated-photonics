function neff = find_neff(n_core, n_clad, d, lambda, polarization, m)
% find_neff Solves for the effective index of a symmetric planar waveguide.
% Inputs:
%   n_core: Core refractive index
%   n_clad: Cladding refractive index
%   d: Core thickness (meters)
%   lambda: Wavelength (meters)
%   polarization: 'TE' or 'TM' (string)
%   m: Mode order (0, 1, 2, ...)
% Output:
%   neff: Effective refractive index (NaN if mode is not supported or error)

k0 = 2 * pi / lambda;

% Define the transcendental equation F(neff_val) = 0
% F(neff_val) = ky * d - 2 * phi_c - m * pi
trans_eq = @(neff_val) calculate_transcendental_value(neff_val, n_core, n_clad, d, k0, polarization, m);

% Search range for neff. Use very small epsilon for tighter bounds.
epsilon_bound = 1e-14; % Smaller epsilon for bounds
neff_low_bound = n_clad + epsilon_bound;
neff_high_bound = n_core - epsilon_bound;

if neff_low_bound >= neff_high_bound
    neff = NaN; % No valid range for neff (e.g. n_core <= n_clad)
    return;
end

options = optimset('TolX', 1e-12, 'Display', 'off'); % Tighter fzero tolerance

val_at_low = trans_eq(neff_low_bound);
val_at_high = trans_eq(neff_high_bound);

% Protective checks for NaN/Inf from trans_eq evaluation
if isnan(val_at_low) || isnan(val_at_high) || isinf(val_at_low) || isinf(val_at_high)
    neff = NaN;
    return;
end

% Expected behavior for guided modes:
% val_at_low (neff near n_clad) should be > 0 (if V > m*pi)
% val_at_high (neff near n_core) should be < 0
cutoff_check_tolerance = 1e-9; % Tighter tolerance for checking if mode is below cutoff
if val_at_low <= cutoff_check_tolerance
    % Mode is at or below cutoff
    neff = NaN;
elseif sign(val_at_low) * sign(val_at_high) > 0
    % Signs are the same, and val_at_low > 0. This means val_at_high is also > 0.
    % This is unexpected as F(n_core-) should be < 0.
    % This might indicate numerical issues or that n_core is very close to n_clad,
    % or the mode is extremely weakly guided and neff_high_bound isn't close enough to n_core
    % for the asymptotic behavior to dominate.
    % For robustness, could try a guess, but for now, assume no reliable root found.
    % fprintf('Warning: Signs at boundaries are same and positive for m=%d, pol=%s, d=%.2fnm. val_low=%.2e, val_high=%.2e\n', m, polarization, d*1e9, val_at_low, val_at_high);
    neff = NaN;
else % val_at_low > 0 and val_at_high < 0 (or very close to 0), so signs are different. Root is bracketed.
    try
        [neff_sol, ~, exitflag] = fzero(trans_eq, [neff_low_bound, neff_high_bound], options);
        if exitflag > 0
            neff = neff_sol;
        else
            neff = NaN; % fzero did not converge
        end
    catch ME
        % In case fzero itself errors (e.g. complex values if neff goes out of bounds during iteration)
        % fprintf('fzero error for m=%d, pol=%s, d=%.2fnm: %s\n', m, polarization, d*1e9, ME.message);
        neff = NaN;
    end
end

% Final validation of the found neff
if ~isnan(neff)
    % Ensure neff is strictly within (n_clad, n_core) after fzero.
    % Allow a tiny margin for fzero's own tolerance effects if neff is extremely close to boundaries.
    if ~isreal(neff) || neff <= n_clad + epsilon_bound/10 || neff >= n_core - epsilon_bound/10
        % If neff is not strictly within the core/cladding range (plus a tiny bit of margin for numerical noise if it's right at the boundary)
        % or if it's not real, then invalidate it.
        % The check `neff <= n_clad` (without epsilon) would be even stricter.
        % The key is that task1.m will filter for neff > n_clad_neff anyway.
        if neff <= n_clad || neff >= n_core % Stricter check
            neff = NaN;
        end
    end
end
end

function val = calculate_transcendental_value(neff_val, n_core, n_clad, d, k0, polarization, m)
% Helper function to calculate the value of the transcendental equation part:
% ky * d - 2 * phi_c - m * pi

% Ensure neff_val is strictly within (n_clad, n_core) for physical solutions
% Allow for slight numerical errors if neff_val is at the boundary due to fzero's iteration.
% The main check in find_neff handles the strict bounds.
if neff_val < n_clad - 1e-15 || neff_val > n_core + 1e-15 % Looser check here, strict check in main function
    val = NaN;
    return;
end

% Clamp neff_val to be strictly within (n_clad, n_core) for calculations
% to avoid issues with sqrt of negative numbers if fzero overshoots slightly.
if neff_val <= n_clad
    neff_val = n_clad + 1e-16; % Push slightly into the guided region
end
if neff_val >= n_core
    neff_val = n_core - 1e-16; % Push slightly into the guided region
end


beta_sq = (neff_val * k0)^2;
n_core_k0_sq = (n_core * k0)^2;
n_clad_k0_sq = (n_clad * k0)^2;

arg_ky_sq = n_core_k0_sq - beta_sq;
arg_gammay_sq = beta_sq - n_clad_k0_sq;

% ky and gammay must be real for guided modes.
% arg_ky_sq > 0 means neff_val < n_core.
% arg_gammay_sq > 0 means neff_val > n_clad.
if arg_ky_sq <= 0 || arg_gammay_sq <= 0
    % This case should be caught by the initial neff_val check,
    % but as a safeguard for numerical precision issues near boundaries:
    val = NaN;
    return;
end

ky = sqrt(arg_ky_sq);
gammay = sqrt(arg_gammay_sq);

if strcmpi(polarization, 'TE')
    % For TE modes: phi_c = atan(gammay / ky)
    atan_arg = gammay / ky;
elseif strcmpi(polarization, 'TM')
    % For TM modes: phi_c = atan((n_core/n_clad)^2 * gammay / ky)
    atan_arg = (n_core/n_clad)^2 * (gammay / ky);
else
    error('Invalid polarization specified. Must be ''TE'' or ''TM''.');
end

% Handle cases where ky is very small (neff_val approaches n_core)
% atan_arg can become very large. atan(large_number) approaches pi/2.
if ky < 1e-9 % Effectively neff_val is n_core
    phi_c = pi/2;
else
    phi_c = atan(atan_arg);
end

val = ky * d - 2 * phi_c - m * pi;
end
