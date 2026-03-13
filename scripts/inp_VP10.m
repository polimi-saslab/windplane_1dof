function par = inp_VP10(par)

%% atmospheric values
par.atm.g = 9.81;     % [m/s^2] gravitational acceleration
par.atm.rho = 1.225;  % [kg/m^3] air density
par.atm.vw_range = (3:0.5:20)';   % [m/s] wind speed range
par.atm.N_vw = length(par.atm.vw_range);
par.N_vw = length(par.atm.vw_range);
par.atm.alpha  = 0.14; % wind shear exponent 
par.atm.h_ref = 100;  % [m] reference altitude

par.atm.A = 10.2; % Weibull form parameter [m/s]
par.atm.k = 2.2; % Weibull shape parameter []
par.atm.gw = par.atm.k/par.atm.A *(par.atm.vw_range/par.atm.A).^(par.atm.k-1).*exp(-(par.atm.vw_range/par.atm.A).^par.atm.k); % Weibull wind distribution
 
%% caractheristics windplane
par.b = 10;      % [m] wing span
par.AR = 5;      % [-] aspect ratio
par.A = par.b^2/par.AR; % [m^2] wing area
par.m = 200; % [kg]

%% Aerodynamic caractheristics
airf = load('airfoil_NACA4421_polars'); % naca4421
par.airf = airf.airf;
par.polars_3D.AoA = par.airf.AoA + par.airf.Cl/(pi*par.AR)*180/pi; % from 2D polars to 3D polars  (elliptical wing)
par.polars_3D.CL = par.airf.Cl;
par.polars_3D.CD = par.airf.Cd + par.airf.Cl.^2./(pi*par.AR); % from 2D polars to 3D polars (elliptical wing)

par.lut_CL = fit(par.polars_3D.AoA,par.polars_3D.CL,'smoothingspline'); % save data as splines for the optimization
par.lut_CD = fit(par.polars_3D.AoA,par.polars_3D.CD,'smoothingspline'); % save data as splines for the optimization

par.lut_CL_casadi = casadi.interpolant('LUT','bspline',{par.polars_3D.AoA},par.polars_3D.CL); % save data as splines for the optimization
par.lut_CD_casadi = casadi.interpolant('LUT','bspline',{par.polars_3D.AoA},par.polars_3D.CD); % save data as splines for the optimization

%% Far wake look-up table
par.lut_Ups_double = load('far_wake_Ups.mat'); % far wake Upsilon term, it is computed offline 
par.lut_Ups = casadi.interpolant('LUT','bspline',{par.lut_Ups_double.lambda0_m(1:end,1),par.lut_Ups_double.eta_m(1,1:end)  }, par.lut_Ups_double.Ups(:)); % save data as splines for the optimization

%% Tether
par.Cd_te = 1;   % [-] tether drag coefficient
par.L_te = 150;  % [m] tether length
par.D_te = 8 * 1e-3; % [m] external tether diameter
par.D_te_el = 2 * 1e-3; % [m] diamenter of the electrical component 
par.A_te = pi*par.D_te^2/4; % [m^2] total sectional area
par.A_te_el = pi*par.D_te_el^2/4; % [m^2]  sectional area of the electrical component 
par.A_te_str = par.A_te-par.A_te_el; % [m^2]  sectional area of the structural component
% par.strenght_te = 1e9; 
par.L_creep_te = [-2.4   8.3  -11.2   5.2]'; % curves for tether life extimation due to creep

%% Onboard turbines
par.P_r = 100e3;   % [W] Rated power
par.P_smooth = 0.1; % [-] max fluctuation of power with respect to rated power
par.xi_t = 0.2;  % [-] R_t/(b/2) onboard wind turbines size
par.R_t = par.xi_t* (par.b/2);  % [-] R_t/(b/2) onboard wind turbines size
par.N_t = 4; % number of turbines
par.A_t = par.N_t * pi * par.R_t^2; % [m^2] total onboard turbine area
par.lut_lambda_t = [1.   1.25 1.5  1.75 2.   2.25 2.5  2.75 3.   3.25 3.5  3.75 4.  ]; % Look-Up Table for the onboard turbine CP and CT
par.lut_C_Pt = [ 0.07261323  0.0998083   0.11937062  0.12506522  0.11724799  0.10274738   0.08053356  0.04920224  0.00957951 -0.03979646 -0.09697764 -0.16253743   -0.23916591];
par.lut_C_Tt = [ 0.09418792  0.12793723  0.14634293  0.1484457   0.13515339  0.11623566    0.09070283  0.05751773  0.01836758 -0.02767859 -0.07832964 -0.13380593   -0.1959909 ]; 

par.CTt_fitp = polyfit(par.lut_lambda_t,par.lut_C_Tt,3); % save data as polynomial for the optimization
par.CPt_fitp = polyfit(par.lut_lambda_t,par.lut_C_Pt,3); % save data as polynomial for the optimization

%% Control and dynamic simulation settings
par.tau = linspace(0,1,51)'; % discretization of time vector between 0 and 1

if strcmp(par.dyn,'dyn') % solve dynamics
    par.k_hb = 10;       % number of harmonics for the harmonic balance method
    par.k_c_theta =1;    % number of harmonics for the control inputs
    par.k_c_t = 2;       % number of harmonics for the control inputs
elseif  strcmp(par.dyn,'steady') % simplify problem to look for steady solutions
    par.k_hb = 0;        % number of harmonics for the harmonic balance method
    par.k_c_theta =0;    % number of harmonics for the control inputs
    par.k_c_t = 0;       % number of harmonics for the control inputs
end
par.solver = 'fmincon'; % fmincon || ipopt

end