function [min] = opt_call_opt_FD(opt,par)

%% Necessary for the 
xLast = []; % Last place opt_obj was called
myf = []; % Use for objective at xLast
myc = []; % Use for nonlinear inequality constraint
myceq = []; % Use for nonlinear equality constraint


%%
opts = optimoptions(@fmincon, ...
    'Diagnostics','off',...%                                        'GradObj','off', ...        % Gradients of the objective fun. provided
    'GradConstr', 'off', ...    % Gradients of the constraints provided
    'FiniteDifferenceType','central',... % forward central    'FiniteDifferenceStepSize',opt.dx, ...
    'Display','iter',...         % iter final
    'Algorithm', 'sqp', ...     % optimization algorithm (default: 'interior-point') sqp active-set
    'MaxIterations', 2000,...           
    'MaxFunctionEvaluations', 1e8);%...

[min.xbest,min.fval,min.exitflag,min.output,min.lambda, min.grad , min.hessian ] = fmincon(@(x) objfun(x,par), opt.x0, opt.A, opt.b, opt.Aeq, opt.beq, opt.lb, opt.ub, @(x) constr(x,par), opts);


%% Define functions to avoid double evaluation of the objective function when checking constraints
    function y = objfun(x,par)
        if ~isequal(x,xLast) % Check if computation is necessary
            [myf,myc,myceq] = opt_main(x,par);
            xLast = x;
        end
        % Now compute objective function
        y = myf;
    end

    function [c,ceq] = constr(x,par)
        if ~isequal(x,xLast) % Check if computation is necessary
            [myf,myc,myceq] = opt_main(x,par);
            xLast = x;
        end
        % Now compute constraint function
        c = myc; % In this case, the computation is trivial
        ceq = myceq;
    end

end