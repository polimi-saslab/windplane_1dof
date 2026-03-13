figure()

nexttile

yyaxis left
plot(par.atm.vw_range,curve.h)
hold on
plot(par.atm.vw_range,curve.h_max,'--')
plot(par.atm.vw_range,curve.h_min,'-.')
ylabel('$h$ (m)')

yyaxis right
plot(par.atm.vw_range,curve.beta)
hold on
plot(par.atm.vw_range,curve.beta+curve.Phi*180/pi,'--')
plot(par.atm.vw_range,curve.beta-curve.Phi*180/pi,'-.')
ylabel('$\beta (^\circ)$')


xlabel('$t/\mathcal{T}$')
xticks(par.atm.vw_range(1:2:end))
xlim([par.atm.vw_range(1), par.atm.vw_range(end)])
xlabel(strcat('$v_w \,$ (h= ',num2str(par.atm.h_ref),'  m) (m/s)'))
grid on


%%
nexttile
plot(par.atm.vw_range,curve.P/1e3)
hold on
plot(par.atm.vw_range,curve.P_max/1e3,'--')
plot(par.atm.vw_range,curve.P_min/1e3,'-.')
plot([par.atm.vw_range(1), par.atm.vw_range(end)],[0 0],'k:',LineWidth=1.2)
ylabel('$P$ (kW)')
xticks(par.atm.vw_range(1:2:end))
xlim([par.atm.vw_range(1), par.atm.vw_range(end)])
xlabel(strcat('$v_w \,$ (h= ',num2str(par.atm.h_ref),'  m) (m/s)'))
% legend('$P_{mean}$','$P_{max}$','$P_{min}$','Location','southeast')
grid on
%%
% figure('units','centimeters','outerposition',plt.dim)
nexttile
hold on
plot(par.atm.vw_range,curve.sigma_te/1e6)
plot(par.atm.vw_range,curve.sigma_te_max/1e6,'--')
plot(par.atm.vw_range,curve.sigma_te_min/1e6,'-.')
ylabel('$\sigma_{te}$ (MPa)')
xticks(par.atm.vw_range(1:2:end))
xlim([par.atm.vw_range(1), par.atm.vw_range(end)])
xlabel(strcat('$v_w \,$ (h= ',num2str(par.atm.h_ref),'  m) (m/s)'))
grid on

%%
nexttile
hold on
plot(par.atm.vw_range,curve.CLh)
plot(par.atm.vw_range,curve.CL_min)
plot(par.atm.vw_range,curve.CL_max)
plot(par.atm.vw_range,curve.CLt,'--')
ylabel('$C_L$')
xticks(par.atm.vw_range(1:2:end))
xlim([par.atm.vw_range(1), par.atm.vw_range(end)])
xlabel(strcat('$v_w \,$ (h= ',num2str(par.atm.h_ref),'  m) (m/s)'))
grid on


%%
nexttile
hold on
plot(par.atm.vw_range,curve.u)
plot(par.atm.vw_range,curve.u_min)
plot(par.atm.vw_range,curve.u_max)
ylabel('$u$ (m/s)')
xticks(par.atm.vw_range(1:2:end))
xlim([par.atm.vw_range(1), par.atm.vw_range(end)])
xlabel(strcat('$v_w \,$ (h= ',num2str(par.atm.h_ref),'  m) (m/s)'))
grid on

%%
nexttile
hold on
plot(par.atm.vw_range,curve.theta)
plot(par.atm.vw_range,curve.theta_min)
plot(par.atm.vw_range,curve.theta_max)
ylabel('$\theta (^\circ)$')
xticks(par.atm.vw_range(1:2:end))
xlim([par.atm.vw_range(1), par.atm.vw_range(end)])
xlabel(strcat('$v_w \,$ (h= ',num2str(par.atm.h_ref),'  m) (m/s)'))
grid on

%%
nexttile
hold on
plot(par.atm.vw_range,curve.lambda_t)
plot(par.atm.vw_range,curve.lambda_t_min)
plot(par.atm.vw_range,curve.lambda_t_max)
ylabel('$\lambda_t$')
xticks(par.atm.vw_range(1:2:end))
xlim([par.atm.vw_range(1), par.atm.vw_range(end)])
xlabel(strcat('$v_w \,$ (h= ',num2str(par.atm.h_ref),'  m) (m/s)'))
grid on


%%
nexttile
hold on
plot(par.atm.vw_range,curve.af_dv)
xticks(par.atm.vw_range(1:2:end))
xlim([par.atm.vw_range(1), par.atm.vw_range(end)])
xlabel(strcat('$v_w \,$ (h= ',num2str(par.atm.h_ref),'  m) (m/s)'))
ylabel('$a_f$')
grid on

