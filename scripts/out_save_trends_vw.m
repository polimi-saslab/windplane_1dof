function curve = out_save_trends_vw(curve,vars,t,out,p)

vw_txt = strcat('vw_', num2str(vars.atm.vw_ref*10));

curve.(vw_txt).t = t;

curve.P(p) = out.P;
curve.P_min(p) = min(t.P);
curve.P_max(p) = max(t.P);

curve.sigma_te(p) = out.sigma_te;
curve.sigma_te_min(p) = min(t.sigma_te);
curve.sigma_te_max(p) = max(t.sigma_te);

curve.u(p) = vars.u(1);
curve.u_min(p) = min(t.u);
curve.u_max(p) = max(t.u);

curve.h(p) = mean(t.h);
curve.h_min(p) = min(t.h);
curve.h_max(p) = max(t.h);

curve.urel(p) = trapz(t.tau,t.u+t.vw1);
curve.urel_min(p) = min(t.u+t.vw1);
curve.urel_max(p) = max(t.u+t.vw1);

curve.theta(p) = vars.theta(1);
curve.theta_min(p) = min(t.theta);
curve.theta_max(p) = max(t.theta);

curve.lambda_t(p) = vars.lambda_t(1);
curve.lambda_t_min(p) = min(t.lambda_t);
curve.lambda_t_max(p) = max(t.lambda_t);

curve.CTt(p) = trapz(t.tau,t.CTt);
curve.CTt_min(p) = min(t.CTt);
curve.CTt_max(p) = max(t.CTt);

curve.CPt(p) = trapz(t.tau,t.CPt);
curve.CPt_min(p) = min(t.CPt);
curve.CPt_max(p) = max(t.CPt);

curve.T(p) = out.T;

curve.T_te(p) = out.T_te;
curve.T_te_min(p) = min(t.T_te);
curve.T_te_max(p) = max(t.T_te);

curve.beta(p,1) = vars.beta_d;
curve.gamma_d(p,1) = vars.gamma_d;

curve.CLh(p,1) = vars.CLh;
curve.CLt(p,1) = trapz(t.tau,t.CL);
curve.CL_min(p,1) = min(t.CL);
curve.CL_max(p,1) = max(t.CL);

curve.af_dv(p) = vars.af;
curve.Phi(p,1) = out.Phi;

% curve.Delta = par.m * par.atm.g /( out.T_te * cosd(vars.beta_d) * (1+tand(vars.beta_d)^2));
curve.Delta(p) = vars.m * vars.atm.g /( out.L * cosd(vars.beta_d) * (1+tand(vars.beta_d)^2) );
curve.eta(p) = (1-  curve.Delta(p) * tand(vars.beta_d))^3;
curve.P_eta(p) = curve.eta(p) * curve.P(p) ;

end