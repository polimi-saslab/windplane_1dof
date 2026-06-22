function [out,t,vars,c,ceq] = opt_model(x,vars)

%% organize design variables
[vars] = opt_dv(x,vars);

t.tau = vars.tau;                      % [-] non-dimensional time 
vars.u = vars.u_vw .* vars.atm.vw_ref; % [m/s] Fourier coefficients of the tangential velocity
vars.beta = vars.beta_d*pi/180;        % [rad]
vars.gamma = vars.gamma_d * pi/180;    % [rad]

%% Compute turning radius
out.M   = vars.m/(1/2*vars.atm.rho * vars.A * vars.L_te *vars.CLh) ; 
out.Phi = acos(-out.M/2+sqrt(out.M.^2+4)/2);   % [rad] opening angle of the trajectory
out.R0 = vars.L_te * sin(out.Phi);             % [m] trajectory radius

%% Compute time course of the tangential speed and acceleration from Fourier cofficients
t.u = vars.u(1)*ones(size(vars.tau));
out.u_hat = vars.u(1);
t.du_dt = zeros(size(vars.tau));
for k = 1:vars.k_hb  
    t.u       = t.u    +    vars.u(1+k) * sin(k*vars.tau*2*pi) +    vars.u(vars.k_hb+1+k) * cos(k*vars.tau*2*pi);  % [m/s] tangential speed
end
out.T = 2*pi*out.R0/(trapz(vars.tau,t.u)); % [s] Period
t.time = vars.tau*out.T;                   % [s] time 
t.Psi = -cumtrapz(t.time,t.u)/out.R0;      % [rad] azimuth position

for k = 1:vars.k_hb
    t.du_dt = t.du_dt  + k * 2*pi/out.T* vars.u(1+k)  * cos(k*vars.tau*2*pi) - 2*pi/out.T * k *  vars.u(vars.k_hb+1+k)  * sin(k*vars.tau*2*pi);   % [m/s] tangential acceleration
end

%% Compute time course of the pitch angle from Fourier cofficients
t.theta = vars.theta(1)*ones(size(t.Psi)); % [rad]
for k = 1:vars.k_c_theta
    t.theta    = t.theta     +    vars.theta(1+k) * sin(k*vars.tau*2*pi)    +  vars.theta(vars.k_c_theta+1+k) * cos(k*vars.tau*2*pi);
end

%% Compute time course of the onboard turbines tip-speed ratio from Fourier cofficients
t.lambda_t = vars.lambda_t(1)*ones(size(t.Psi)); % [-]
for k = 1:vars.k_c_t
    t.lambda_t = t.lambda_t  +     vars.lambda_t(1+k) * sin(k*vars.tau*2*pi) +    vars.lambda_t(vars.k_c_t+1+k) * cos(k*vars.tau*2*pi); 
end

%% Evaluate wind speed 
t.h = vars.L_te*cos(out.Phi)*sin(vars.beta) - out.R0*cos(vars.beta) .* sin(t.Psi); % [m] altitude as function of azimuth
out.h_min = vars.L_te*cos(out.Phi)*sin(vars.beta) -  out.R0*cos(vars.beta);        % [m] minimum altitude over the loop
t.h = (t.h + sqrt(t.h.^2))/2+1e-3 ; % this is to prevent errors is the altitude goes negative during the optimization process                                     
t.vw = vars.atm.vw_ref * (t.h./vars.atm.h_ref).^vars.atm.alpha;                          % [m/s] norm of the wind velocity
t.vw1 = (sin(vars.beta)*cos(vars.gamma)*cos(t.Psi) + sin(vars.gamma)*sin(t.Psi)).*t.vw;  % [m/s] wind velocity in the cylindrical coordinate system (tangential)
t.vw2 = (-sin(vars.beta)*cos(vars.gamma)*sin(t.Psi) + sin(vars.gamma)*cos(t.Psi)).*t.vw; % [m/s] wind velocity in the cylindrical coordinate system (radial)
t.vw3 = (cos(vars.beta)*cos(vars.gamma)).*t.vw;                                          % [m/s] wind velocity in the cylindrical coordinate system (axial)

if strcmp(vars.dyn,'steady')
    t.vw1 = 0*t.vw1;
    t.vw2 = 0*t.vw2;
    t.vw3 = -1/(2*pi)*trapz(t.Psi,t.vw3) * ones(size(t.vw1));
    vars.atm.g = 0;
end
out.vw3 = trapz(vars.tau,t.vw3); % [m/s]

%% Compute aerodynamic coefficients of the wing 
t.gamma_n = atan((1-vars.af).*t.vw3./(t.vw1+t.u))*180/pi; % [rad] inflow angle in the near wake
t.AoA =  t.theta + t.gamma_n;
if strcmp(vars.jac,'AD')
    t.CL    = vars.lut_CL_casadi(t.AoA);
    t.CD_3D = vars.lut_CD_casadi(t.AoA);
else
    t.CL   = vars.lut_CL(t.AoA);
    t.CD_3D = vars.lut_CD(t.AoA);
end
out.CLh = trapz(vars.tau,t.CL);
out.CD_te  = vars.Cd_te * vars.D_te * vars.L_te/(4*vars.A) ;
t.CD  = out.CD_te + t.CD_3D;

%% Compute onboard turbines coefficients
t.CTt = polyval(vars.CTt_fitp',t.lambda_t);
t.CPt = polyval(vars.CPt_fitp',t.lambda_t);

%%
t.h_tau = (1/2*vars.atm.rho * vars.A * t.CL .* (t.u+t.vw1).* t.vw3.*(1-vars.af)  ... .* cos(t.phi).^2
    - 1/2*vars.atm.rho * vars.A *t.CD .* (t.u+t.vw1).^2 ...
    - 1/2*vars.atm.rho * vars.A_t * t.CTt .* (t.u+t.vw1).^2 ...
    - vars.m * vars.atm.g .* cos(t.Psi).* cos(vars.beta) ...
    - vars.m * t.du_dt)/(1/2*vars.atm.rho * vars.A * vars.atm.vw_ref^2*50 ) ; 

%% Compute additional outputs
t.L   = 1/2 * vars.atm.rho * vars.A * t.CL .*  (t.u+t.vw1).^2; % [N] lift
out.L = trapz(vars.tau,t.L);% [N] mean lift

t.T_te = (t.L - vars.atm.g*vars.m * sin(vars.beta)) * sqrt(1+tan(out.Phi)^2); % [N] tether force
out.T_te = trapz(vars.tau,t.T_te);% [N] mean tether force

t.sigma_te = t.T_te/vars.A_te_str;  % [Pa] tether stress
out.sigma_te = out.T_te/vars.A_te_str;  % [Pa] mean tether stress

t.phi   = - vars.m*vars.atm.g*(sin(vars.beta)*tan(out.Phi)+sin(t.Psi).*cos(vars.beta))./t.L; % [rad] roll angle
t.phi_d = t.phi*180/pi; % [deg] roll angle

%% Compute power
t.P      = 1/2 * vars.atm.rho *vars.A_t * t.CPt .*  (t.u+t.vw1).^3*0.7; % [W] Power
out.P      = trapz(vars.tau,t.P); % [W] Mean power
out.CP      = trapz(vars.tau,t.P); % [W] Mean power

%% Evaluate far wake induction
out.Gamma_0 =  2*vars.b*vars.u(1)*out.CLh/(pi*vars.AR); % tip vortex strenght
out.y_v = vars.b*pi/8;      % [m] tip vortex spanwise position
out.eta_v = out.y_v/out.R0; % [-] tip vortex non dimensional spanwise position
out.a_r = out.Gamma_0/(4*pi*out.y_v.*out.vw3); % [-] vw(1-ar) is the velocity of the vortex rings
out.lambda = trapz(vars.tau,t.u./t.vw3); % [-] mean wing speed ratio
out.lambda0 = out.lambda/(1-out.a_r);% [-] torsional parameter of the helicoidal wake
if strcmp(vars.jac,'FD') % look-up table for the far wake induction
    out.Ups = interp2(vars.lut_Ups_double.eta_m(1,1:end),vars.lut_Ups_double.lambda0_m(1:end,1)',vars.lut_Ups_double.Ups,out.eta_v,out.lambda0);
elseif strcmp(vars.jac,'AD')
    out.Ups = vars.lut_Ups([out.lambda0,out.eta_v]);
end
out.af = out.a_r * out.Ups;  % far wake induction


%% --------------------  Constraints -------------------- 

%% far wake (slack variable)
wake_eq = vars.af - out.af;
%% CL (slack variable)
CL_eq = vars.CLh - out.CLh;
%% compute Fourier coefficients of the residual of the equation of motion
if strcmp(vars.jac,'FD')
    H = zeros(2*vars.k_hb+1,1);
elseif strcmp(vars.jac,'AD')
    H = casadi.MX(2*vars.k_hb+1,1);
end
H(1) = trapz(vars.tau,t.h_tau);
for k = 1:vars.k_hb
    H(1+k)          =  trapz(vars.tau,t.h_tau.*cos(k*vars.tau*2*pi));
    H(1+vars.k_hb+k) = trapz(vars.tau,t.h_tau.*sin(k*vars.tau*2*pi));
end

%% compute Fourier coefficients of power for the constraint of power smoothing
if strcmp(vars.jac,'FD')
    P = zeros(2*vars.k_hb+1,1);
elseif strcmp(vars.jac,'AD')
    P = casadi.MX(2*vars.k_hb+1,1);
end
P(1) = trapz(vars.tau,t.P);
for k = 1:vars.k_hb
    P(1+k)          =  2*trapz(vars.tau,t.P.*cos(k*vars.tau*2*pi));
    P(1+vars.k_hb+k) = 2*trapz(vars.tau,t.P.*sin(k*vars.tau*2*pi));
end
P_fl = norm([P(2:end)]);
Pfl_c = P_fl/vars.P_r - vars.P_smooth;   % power smoothing constraint
P_Pr_c = out.P/(vars.P_r) - 1; % rated power equality constraint

%% Concatenate constraints
if strcmp (vars.region, 'I')
    ceq = [wake_eq;CL_eq;H];
    if strcmp(vars.dyn,'dyn')
        c = [Pfl_c;1-out.h_min/(2*vars.b); 1-t.sigma_te/(20e6); out.sigma_te/(200e6)-1];
    elseif strcmp(vars.dyn,'steady')
        c = [1-out.h_min/(2*vars.b); 1-out.sigma_te/(20e6); out.sigma_te/(200e6)-1];
    end

    % if strcmp(vars.dyn,'dyn')
    %     c = [Pfl_c; 1-t.sigma_te/(20e6); out.sigma_te/(70e6)-1];
    % elseif strcmp(vars.dyn,'steady')
    %     c = [1-out.sigma_te/(20e6); out.sigma_te/(70e6)-1];
    % end
elseif strcmp (vars.region, 'II')
    ceq = [wake_eq;CL_eq;H;P_Pr_c];
    if strcmp(vars.dyn,'dyn')
        c = [Pfl_c;1-out.h_min/(2*vars.b); 1-t.sigma_te/(0.02e9)];
    elseif strcmp(vars.dyn,'steady')
        c = [1-out.h_min/(2*vars.b); 1-out.sigma_te/(0.02e9)];
    end
end

end


