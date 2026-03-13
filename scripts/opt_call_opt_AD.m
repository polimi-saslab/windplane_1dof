function [res] = opt_call_opt_AD(opt,par)
% create optimization problem
%%
x = casadi.MX.sym('x',length(opt.x0)); %
[f,c,ceq] = opt_main(x,par);

if strcmp(par.solver,'ipopt')

    g = [ceq;c];
    lbg = [ 0*ones(size(ceq)); -inf*ones(size(c))];
    ubg = [ 0*ones(size(ceq));  zeros(size(c))];

    nlp = struct('x', x, 'f', f , 'g',  g);

    opts = struct;
    opts.ipopt.max_iter = 200;
    opts.ipopt.acceptable_obj_change_tol = 1e-4;
    opts.ipopt.acceptable_dual_inf_tol = 1e6;

    % Create IPOPT solver object
    solver = casadi.nlpsol('solver', 'ipopt', nlp ,opts);

    % Solve the NLP
    res = solver('x0' , opt.x0,... % solution guess
        'lbx', opt.lb,...           % lower bound on x
        'ubx', opt.ub,...           % upper bound on x
        'lbg', lbg,...           % lower bound on g
        'ubg', ubg);             % upper bound on g

    %%
    res.xbest = full(res.x);
    res.fval = full(res.f);

elseif  strcmp(par.solver,'fmincon')

    g = [ceq;c];
    lbg = [0*ones(size(ceq)); -inf*ones(size(c))];
    ubg = [ 0*ones(size(ceq));  zeros(size(c))];

    % Function to compute objective and its gradient
    f = casadi.Function('f',{x},{f,gradient(f,x)});

    % Function to compute constraint vector, its Jacobian, and bounds
    g = casadi.Function('g',{x},{g,jacobian(g,x),lbg,ubg});

    options = optimoptions('fmincon',...
        'Display','off',...         % iter final
        'Algorithm', 'sqp' , ...     % interior-point sqp active-set
        'MaxIterations', 3000,...
        'MaxFunctionEvaluations',1e6,...
        'SpecifyObjectiveGradient',true,...
        'SpecifyConstraintGradient',true,...        
        'StepTolerance',1e-12,...
        'ConstraintTolerance',1e-12);


    [res.xbest,res.f,res.exitflag,res.output,res.lambda,res.grad,res.hessian] = fmincon(@(x) opt_obj_casadi(f,x),opt.x0, opt.A, opt.b, opt.Aeq, opt.beq, opt.lb, opt.ub, @(x) opt_nonlin_casadi(g,x),options);

    disp(strcat('exit flag       =      ', num2str(res.exitflag)))
end

end

