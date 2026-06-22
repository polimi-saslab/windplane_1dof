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


%% Time-domain plots for a few representative wind speeds
% Pick 5 indices spread across vw_range instead of hardcoded vw_XX field
% names, so this section works regardless of par.atm.vw_range content.
idx_sample = round(linspace(1, par.atm.N_vw, 5));
idx_sample = unique(idx_sample, 'stable');
vw_sample  = par.atm.vw_range(idx_sample);
fields_vw  = arrayfun(@(v) strcat('vw_', num2str(v*10)), vw_sample, 'UniformOutput', false);
legend_vw  = arrayfun(@(v) sprintf('$v_{w}$ = %g m/s', v), vw_sample, 'UniformOutput', false);

%%
figure('units','centimeters','outerposition',plt.dim)
hold on
for i = 1:numel(fields_vw)
    plot(curve.(fields_vw{i}).t.tau, curve.(fields_vw{i}).t.P/1e3)
end
grid on
ylabel('$P$ (kW)')
xlabel('$t/{T}$ (-)')
legend(legend_vw,'Location','best')
box on


%%
figure('units','centimeters','outerposition',plt.dim)
hold on
for i = 1:numel(fields_vw)
    plot(curve.(fields_vw{i}).t.tau, curve.(fields_vw{i}).t.lambda_t)
end
grid on
ylabel('$\lambda_t$ (-)')
xlabel('$t/{T}$ (-)')
legend(legend_vw,'Location','best')
box on


%%
figure('units','centimeters','outerposition',plt.dim)
hold on
for i = 1:numel(fields_vw)
    plot(curve.(fields_vw{i}).t.tau, curve.(fields_vw{i}).t.CL)
end
grid on
ylabel('$C_L$ (-)')
xlabel('$t/{T}$ (-)')
legend(legend_vw,'Location','best')
box on


%%
figure('units','centimeters','outerposition',plt.dim)
hold on
for i = 1:numel(fields_vw)
    plot(curve.(fields_vw{i}).t.tau, curve.(fields_vw{i}).t.sigma_te/1e6)
end
grid on
ylabel('$\sigma_{te}$ (MPa)')
xlabel('$t/{T}$ (-)')
legend(legend_vw,'Location','best')
