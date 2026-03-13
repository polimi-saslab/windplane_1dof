clc
clear 
close all
warning off
%% Figures settings
set(groot,'defaultAxesFontSize',20)
set(groot, 'DefaultLineLineWidth', 2);
set(groot, 'DefaultTextInterpreter', 'Latex')
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');
addpath('scripts','wind clusters','yaml');
addpath('C:\Users\trevi\OneDrive - Politecnico di Milano\Desktop\PhD\Working\casadi-3.7.1-windows64-matlab2018b')
import casadi.*
plt.LineStyle = '-';
plt.LineStyle2 = '--';
plt.dim = [0 3 20 11.5];
plt.b = [0, 0.4470, 0.7410];
plt.r = [0.8500, 0.3250, 0.0980];
plt.y = [0.9290, 0.6940, 0.1250];
plt.v = [0.4940 0.1840 0.5560];

%% Inputs
par.dyn = 'dyn';     % dyn || steady
par = inp_VP10(par); % Input parameters

%% Define design variables & bounds
[opt,par] = opt_problem_settings(par);

%% Power curve
p=1;
curve = struct;
par.region='I';

for vw = par.atm.vw_range' %  for each wind speed
    par.atm.vw_ref = vw;

    % solve optimal problem
    flag_P = 0;
    while flag_P == 0
        tic
        par.jac = 'AD';
        if strcmp(par.jac,'AD')
            [res] = opt_call_opt_AD(opt,par);
        elseif strcmp(par.jac,'FD')
            [res] = opt_call_opt_FD(opt,par);
        end
        par.jac = 'FD';
        [out,t,vars,c,ceq] = opt_model(res.xbest,par); % evaluate optimal solution
        opt.x0 = res.xbest;
        toc

        if out.P>par.P_r && strcmp(par.region,'I') % change to region 2 if power larger than rated power
            opt.lb(3) = opt.x0(3);
            par.region = 'II';
        else
            flag_P = 1;
        end
    end

    % save outputs
    curve = out_save_trends_vw(curve,vars,t,out,p);

    disp(strcat('vw = ',num2str(vw), ' m/s'))
    disp(strcat('P = ',num2str(out.P/1e3), ' kW'))
    disp(strcat('--------'))

    p = p+1;
end


%%
if isscalar(par.atm.vw_range) % plots for a given wind speed

    %%
    disp(strcat('sigma  = ', num2str(out.sigma_te/1e6),' MPa'))
    disp(strcat('T  = ', num2str(out.T),' s'))
    disp(strcat('u  = ', num2str(vars.u(1)),' m/s'))
    disp(strcat('beta  = ', num2str(vars.beta*180/pi),' deg'))
    disp(strcat('gamma  = ', num2str(vars.gamma*180/pi),' deg'))

else % plots as a function of wind speed

    CF = trapz(par.atm.vw_range,curve.P.*par.atm.gw')/par.P_r;
    f_creep = trapz(par.atm.vw_range,par.atm.gw'./10.^(polyval(par.L_creep_te,curve.sigma_te/1e9)) );

    disp(['CF = ',num2str(CF)])
    disp(['Life tether = ',num2str(1/f_creep),' years'])

    %%
    out_plots
    
end
