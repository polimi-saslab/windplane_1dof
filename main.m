clc
clear 
close all
warning off
addpath('scripts');
settings
%% Inputs
par.dyn = 'dyn';     % dyn || steady
par = inp_VP3(par); % Input parameters
% par = inp_VP6(par); % Input parameters
% par = inp_VP30(par); % Input parameters
% par = inp_VP10(par); % Input parameters

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

    % out.P = out.P*0.7;

        if out.P>par.P_r && strcmp(par.region,'I') % change to region 2 if power larger than rated power
            % opt.lb(3) = opt.x0(3);
            par.region = 'II';
            opt.x0 = [0.01   1  30    0  5  zeros(1,2*par.k_hb)    3    zeros(1,2*par.k_c_t) 0  zeros(1,2*par.k_c_theta)   ];             % inital values
        else
            flag_P = 1;
        end
    end

    % save outputs
    curve = out_save_trends_vw(curve,vars,t,out,p);

    disp(strcat('vw = ',num2str(vw), ' m/s'))
    disp(strcat('P = ',num2str(out.P/1e3), ' kW'))
    disp(strcat('CP = ',num2str(out.P/(1/2*par.atm.rho*vw^3*pi*par.b^2))))
    % disp(strcat('xiP = ',num2str(out.P/(1/2*par.atm.rho*vw^3*par.A))))
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

    AEP = CF * par.P_r * 8760 / 1e3; % [kWh/year]

    % Competitor turbines CF and AEP at the same site
    comp = competitor_CF(par.atm);

    disp('--- Windplane ---')
    disp(['CF  = ',num2str(CF,'%.3f')])
    disp(['AEP = ',num2str(AEP,'%.0f'),' kWh/year'])
    disp(['Life tether = ',num2str(1/f_creep),' years'])
    disp('--- Competitors ---')
    fields = fieldnames(comp);
    for i = 1:numel(fields)
        t = comp.(fields{i});
        disp([t.name,': CF = ',num2str(t.CF,'%.3f'), ...
              ',  AEP = ',num2str(t.AEP,'%.0f'),' kWh/year'])
    end

    %%
    out_plots
    
end

%%
figure()
plot(par.atm.vw_range,curve.P/1e3)
hold on
plot(par.atm.vw_range,curve.P_max/1e3,'--')
plot(par.atm.vw_range,curve.P_min/1e3,'-.')
plot([par.atm.vw_range(1), par.atm.vw_range(end)],[0 0],'k:',LineWidth=1.2)
plot(par.atm.vw_range,par.atm.gw*20)
ylabel('$P$ (kW)')
xticks(par.atm.vw_range(1:2:end))
xlim([par.atm.vw_range(1), par.atm.vw_range(end)])
xlabel(strcat('$v_w \,$ (h= ',num2str(par.atm.h_ref),'  m) (m/s)'))
% legend('$P_{mean}$','$P_{max}$','$P_{min}$','Location','southeast')
grid on

%%

% %%
% figure()
% plot(t.Psi*180/pi,t.u)
% 
% %%
% figure()
% plot(t.Psi*180/pi,t.L)
% 
% %%
% figure()
% plot(t.Psi*180/pi,t.AoA)
% 
% %%
% figure()
% plot(t.Psi*180/pi,t.theta)


%%
%% 
% figure('units','centimeters','outerposition',plt.dim)
% hold on
% plot(curve.vw_30.t.tau,curve.vw_30.t.u)
% plot(curve.vw_60.t.tau,curve.vw_60.t.u)
% plot(curve.vw_90.t.tau,curve.vw_90.t.u)
% plot(curve.vw_120.t.tau,curve.vw_120.t.u)
% plot(curve.vw_200.t.tau,curve.vw_200.t.u)
% grid on
% ylabel('$u$ (m/s)')
% xlabel('$t/{T}$ (-)')
% % legend('$v_{w}$ = 3 m/s', '$v_{w}$ = 6 m/s', '$v_{w}$ = 9 m/s', '$v_{w}$ = 12 m/s', '$v_{w}$ = 20 m/s'  )
% box on

%%
figure('units','centimeters','outerposition',plt.dim)
hold on
plot(curve.vw_30.t.tau,curve.vw_30.t.P/1e3)
%% 
plot(curve.vw_60.t.tau,curve.vw_60.t.P/1e3)
plot(curve.vw_90.t.tau,curve.vw_90.t.P/1e3)
plot(curve.vw_120.t.tau,curve.vw_120.t.P/1e3)
plot(curve.vw_200.t.tau,curve.vw_200.t.P/1e3)
grid on
ylabel('$P$ (kW)')
xlabel('$t/{T}$ (-)')
%legend('$v_{w}$ = 3 m/s', '$v_{w}$ = 6 m/s', '$v_{w}$ = 9 m/s', '$v_{w}$ = 12 m/s', '$v_{w}$ = 20 m/s'  )
ylim([-25 , 125])
box on


%%
figure('units','centimeters','outerposition',plt.dim)
hold on
plot(curve.vw_30.t.tau,curve.vw_30.t.lambda_t)
plot(curve.vw_60.t.tau,curve.vw_60.t.lambda_t)
plot(curve.vw_90.t.tau,curve.vw_90.t.lambda_t)
plot(curve.vw_120.t.tau,curve.vw_120.t.lambda_t)
plot(curve.vw_200.t.tau,curve.vw_200.t.lambda_t)
grid on
ylabel('$\lambda_t$ (-)')
xlabel('$t/{T}$ (-)')
legend('$v_{w}$ = 3 m/s', '$v_{w}$ = 6 m/s', '$v_{w}$ = 9 m/s', '$v_{w}$ = 12 m/s', '$v_{w}$ = 20 m/s'  )
box on


%%
figure('units','centimeters','outerposition',plt.dim)
hold on
plot(curve.vw_30.t.tau,curve.vw_30.t.CL)
plot(curve.vw_60.t.tau,curve.vw_60.t.CL)
plot(curve.vw_90.t.tau,curve.vw_90.t.CL)
plot(curve.vw_120.t.tau,curve.vw_120.t.CL)
plot(curve.vw_200.t.tau,curve.vw_200.t.CL)
grid on
ylabel('$C_L$ (-)')
xlabel('$t/{T}$ (-)')
box on
% legend('$v_{w}$ = 3 m/s', '$v_{w}$ = 6 m/s', '$v_{w}$ = 9 m/s', '$v_{w}$ = 12 m/s', '$v_{w}$ = 20 m/s'  )


%%
figure('units','centimeters','outerposition',plt.dim)
hold on
plot(curve.vw_30.t.tau,curve.vw_30.t.sigma_te/1e6)
plot(curve.vw_60.t.tau,curve.vw_60.t.sigma_te/1e6)
plot(curve.vw_90.t.tau,curve.vw_90.t.sigma_te/1e6)
plot(curve.vw_120.t.tau,curve.vw_120.t.sigma_te/1e6)
plot(curve.vw_200.t.tau,curve.vw_200.t.sigma_te/1e6)
grid on
ylabel('$\sigma_{te}$ (MPa)')
xlabel('$t/{T}$ (-)')
legend('$v_{w}$ = 3 m/s', '$v_{w}$ = 6 m/s', '$v_{w}$ = 9 m/s', '$v_{w}$ = 12 m/s', '$v_{w}$ = 20 m/s'  )
