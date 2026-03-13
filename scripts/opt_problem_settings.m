function [opt,par] = opt_problem_settings(par)

par.dv_name = {'af', 'CLh' ,'beta_d','gamma_d','u_vw', 'lambda_t', 'theta'};                                             % design variable names. 
% a_f is the far wake infuction (slack variable)
% CLh is the mean value of the lift coefficient (slack variable)
% beta_d is the elevation angle of the trajectory in degrees
% gamma_d is the yaw angle of the trajectory in degrees
% u_vw are the Fourier coefficients of the normalization of the tangential velocity with the wind speed at the reference altitude u/v_{w,ref}
% lambda_t are the Fourier coefficients of the onboard turbines tip-speed ratio
% theta are the Fourier coefficients of the wing pitch angle

par.dv_dim = [1, 1, 1, 1, 2*par.k_hb+1 ,2*par.k_c_t+1,2*par.k_c_theta+1];                                                % dimensions
opt.x0 = [0.01 1  30 0  6  zeros(1,2*par.k_hb)    3    zeros(1,2*par.k_c_t) 0  zeros(1,2*par.k_c_theta)   ];             % inital values
opt.lb = [0 0.3   10  -20 1  -10*ones(1,2*par.k_hb)  1  -1*ones(1,2*par.k_c_t)  -30   -10*ones(1,2*par.k_c_theta) ];      % lower bounds
opt.ub = [0.2 1.5  50  20 10  10*ones(1,2*par.k_hb)  4.5  1*ones(1,2*par.k_c_t)     10    10*ones(1,2*par.k_c_theta)];   % upper bounds

par.N_dv = length(opt.x0);

opt.A = [];
opt.b = [];
opt.Aeq = [];
opt.beq = [];

end