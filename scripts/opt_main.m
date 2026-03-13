function [f,c,ceq] = opt_main(x,par)

[out,~,~,c,ceq] = opt_model(x,par);
if strcmp (par.region, 'I')
    f = - out.P/par.P_r ;
elseif strcmp (par.region, 'II')
    f = out.sigma_te/1e9;
end

end