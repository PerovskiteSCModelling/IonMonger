function sol_init = precondition(sol_init,params,vectors,matrices,options)
% This function can be used to precondition the cell at a particular
% applied voltage (not equal to the cell's built-in voltage) by computing
% consistent, steady-state initial conditions. The inputs are a vector of 
% the original initial conditions, structures for the parameters, vectors
% and matrices required by the solvers and the options structure for
% ode15s. The output is the new vector of initial conditions.

% Parameter input
[psi, G, Verbose] = struct2array(params, {'psi','G','Verbose'});

% Temporarily turn the following warning into an error
warnon = warning('error','MATLAB:ode15s:IntegrationTolNotMet');

% Compute the Jacobian, mass matrix and initial slope and add to options
% Note that the mass matrix is adjusted so ions are as mobile as electrons
if exist('AnJac','file')
    options.Jacobian = @(t,u) AnJac(t,u,params,vectors,matrices);
else
    options.JPattern = Jac(params);
end
options.Mass = mass_matrix(params,vectors,'precondition');
options.InitialSlope = RHS(0,sol_init,@(t) 0,params,vectors,matrices) ...
    \options.Mass;
options.OutputFcn = []; % Surpress output during preconditioning

% Fix the light in its initial state
params.G = @(x,t) G(x,0);

% Evolve the solution from Vbi to the preconditioning voltage
[~,numsol] = ode15s(@(t,u) RHS(t,u,@(t) psi(0)*t/10,params,vectors,matrices), ...
    [0 10],sol_init,options);
sol_init = numsol(end,:)';

% Reset error message to warning
warning(warnon);

% Define the settings for the call to fsolve
fsoptions = optimoptions('fsolve','MaxIterations',40);
if Verbose, fsoptions.Display = 'iter'; else, fsoptions.Display = 'off'; end

% Use the initial estimate to obtain an approximate steady-state solution
if exist('AnJac.m','file')
    fsoptions.SpecifyObjectiveGradient = true;
    [sol_init,~,exitflag,~] = fsolve(@(u) RHS_AnJac(u,psi, ...
        params,vectors,matrices,'init'),sol_init,fsoptions);
else
    fsoptions.JacobPattern = Jac(params,'init');
    [sol_init,~,exitflag,~] = fsolve(@(u) RHS(0,u,psi, ...
        params,vectors,matrices,'init'),sol_init,fsoptions);
end
if exitflag<1
    warning(['Steady-state initial conditions could not be found to ' ...
        'a high degree of accuracy and may be unphysical.']);
end

% Ensure all the algebraic equations are satisfied as exactly as possible
sol_init = apply_Poisson(sol_init,params,vectors,matrices);

end