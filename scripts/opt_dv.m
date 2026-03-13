function [par] = opt_dv(x,par)

%% Save dimensional design variables in organize structure
p=1;
for q = 1:length(par.dv_name)
    par.(par.dv_name{q}) = x(p:p+par.dv_dim(q)-1);
    p=p+par.dv_dim(q);
end

end