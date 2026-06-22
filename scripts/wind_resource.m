function atm = wind_resource(site, vw_range, h_ref)
% WIND_RESOURCE  Wind resource estimation for a site.
%
%   atm = wind_resource(site, vw_range, h_ref)
%
%   Inputs:
%     site      - string  : site identifier ('Linosa', 'Affi', ...)
%                 OR
%                 [lat, lon] : decimal degrees -> ERA5 data fetched automatically
%     vw_range  - [m/s] column vector of wind speeds for the Weibull PDF
%     h_ref     - [m]   reference altitude for the Weibull distribution
%
%   Output:
%     atm  struct with fields:
%       .g            gravitational acceleration [m/s^2]
%       .rho          air density [kg/m^3]
%       .vw_range     wind speed range [m/s]
%       .N_vw         number of wind speed bins
%       .wind_profile [height (m), mean wind speed (m/s)] table
%       .alpha        fitted wind shear exponent [-]
%       .h_ref        reference altitude [m]
%       .v_mean       mean wind speed at h_ref [m/s]
%       .k            Weibull shape parameter [-]
%       .A            Weibull scale parameter [m/s]
%       .gw           Weibull PDF evaluated on vw_range
%
%   ERA5 mode notes:
%     Requires Python with cdsapi and netCDF4 installed, and a valid
%     ~/.cdsapirc file (see https://cds.climate.copernicus.eu).
%     The downloaded NetCDF is cached in <project_root>/cache/era5/
%     so subsequent calls with the same coordinates are instantaneous.
%     The first download may take several minutes due to CDS queue.

%% Constants
atm.g        = 9.81;    % [m/s^2]
atm.rho      = 1.225;   % [kg/m^3]
atm.vw_range = vw_range;
atm.N_vw     = length(vw_range);
atm.h_ref    = h_ref;

%% Wind profile: ERA5 (lat/lon) or hardcoded site table (string)
if isnumeric(site) && numel(site) == 2
    atm.wind_profile = fetch_ERA5(site(1), site(2));
else
    atm.wind_profile = site_table(site);
end

%% Fit power-law wind shear: v(h) = C * h^alpha
% Log-linearisation: log(v) = log(C) + alpha*log(h)
h_data = atm.wind_profile(:,1);
v_data = atm.wind_profile(:,2);
X      = [ones(size(h_data)), log(h_data)];
coeffs = X \ log(v_data);

atm.alpha  = coeffs(2);
atm.v_mean = exp(coeffs(1)) * h_ref^atm.alpha;

%% Weibull distribution at h_ref
atm.k  = 2;                                          % Rayleigh (shape)
atm.A  = atm.v_mean / gamma(1 + 1/atm.k);           % scale parameter
atm.gw = atm.k/atm.A * (vw_range/atm.A).^(atm.k-1) .* ...
         exp(-(vw_range/atm.A).^atm.k);              % Weibull PDF

end


% =========================================================================
function profile = fetch_ERA5(lat, lon)
% Call fetch_ERA5_wind.py and return wind_profile table [h, v_mean].
% Results are cached as NetCDF; only the first call downloads from CDS.

script_dir = fileparts(mfilename('fullpath'));
cache_dir  = fullfile(script_dir, '..', 'cache', 'era5');
py_script  = fullfile(script_dir, 'fetch_ERA5_wind.py');

cmd = sprintf('python "%s" %.4f %.4f "%s"', py_script, lat, lon, cache_dir);
fprintf('wind_resource: querying ERA5 for (%.4f, %.4f) ...\n', lat, lon);
[status, out] = system(cmd);

if status ~= 0
    error(['wind_resource: ERA5 fetch failed.\n' ...
           'Make sure Python, cdsapi and netCDF4 are installed and\n' ...
           '~/.cdsapirc is configured.\nSystem output:\n%s'], out);
end

% Parse JSON output: {"heights": [...], "v_mean": [...]}
data    = jsondecode(strtrim(out));
profile = [data.heights(:), data.v_mean(:)];

end


% =========================================================================
function profile = site_table(site)
% Hardcoded wind profiles from Global Wind Atlas.
% Each row: [height (m), mean wind speed (m/s)]

switch lower(site)

    case 'linosa'
        profile = [ 10,  6.14; ...
                    50,  7.39; ...
                   100,  8.00 ];
    case 'ustica'
        profile = [ 10,  4.91; ...
                    50,  7.17; ...
                   100,   7.3];

    case 'affi'
        profile = [  10,  3.33; ...
                     50,  3.91; ...
                    100,  4.29; ...
                    150,  4.50; ...
                    200,  4.67 ];

    otherwise
        error('wind_resource: unknown site "%s". Add a case or pass [lat, lon].', site);
end

end
