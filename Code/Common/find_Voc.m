function [psi, sol_init] = find_Voc(sol_init,psi,params,vectors,matrices,options)
% This function can be used to find the steady-state open-circuit voltage
% (Voc), compute consistent initial conditions and adjust the voltage
% protocol such that it starts from Voc. The inputs are a vector of the
% original initial conditions, structures for the parameters, vectors and
% matrices required by the solvers and the options structure for ode15s.
% The outputs are the new voltage protocol and corresponding vector of
% initial conditions for a cell preconditioned at open-circuit.

% Parameter input
[N, Kn, Kp, G, time, splits, atol, tstar2t, psi2Vap, Verbose] ...
    = struct2array(params, {'N','Kn','Kp','G','time','splits','atol', ...
    'tstar2t','psi2Vap','Verbose'});
dx = vectors.dx;

% Temporarily turn the following warning into an error
warnon = warning('error','MATLAB:ode15s:IntegrationTolNotMet');

% Compute the Jacobian, mass matrix and initial slope and add to options
% Note that the mass matrix is adjusted so ions are as mobile as electrons
options.JPattern = Jac(params);
options.Mass = mass_matrix(params,vectors,'precondition');
options.InitialSlope = RHS(0,sol_init,@(t) 0,params,vectors,matrices) ...
                        \options.Mass;

% Fix the light in its initial state
params.G = @(x,t) G(x,0);

% Define the total current density (nelecting the displacement current) and
% the Event function to find when the current is zero at...
mid = ceil((N+1)/2); % the midpoint of the perovskite layer
J = @(mid,y) Kn./dx(mid).*(y(2*N+3+mid)-y(2*N+2+mid)-(y(2*N+3+mid)+y(2*N+2+mid)).*(y(N+2+mid)-y(N+1+mid))./2) ...
            -Kp./dx(mid).*(y(3*N+4+mid)-y(3*N+3+mid)+(y(3*N+4+mid)+y(3*N+3+mid)).*(y(N+2+mid)-y(N+1+mid))./2);
function [value,isterminal,direction] = Voc_event(t,y,direction)
    value = J(mid,y); % event is when current density equals zero
    isterminal = 1; % stop at open-circuit
    % direction = -1 for an increasing voltage and 1 for decreasing
end

% Try to locate the open-circuit voltage
try % decreasing the voltage by ~1.5V
    if Verbose
        disp('Attempting to find Voc between Vbi and Vbi-1.5V');
    end
    options.Events = @(t,y) Voc_event(t,y,1);
    [~,~,~,ye,~] = ode15s(@(t,u) RHS(t,u,@(t) t,params,vectors,matrices), ...
        [0 75],sol_init,options);
    if isempty(ye), error('No open-circuit events found.'); end
    sol_init = ye(end,:)';  
catch % then try increasing the voltage by ~1.5V
    if Verbose
        disp('Attempting to find Voc between Vbi and Vbi+0.5V');
    end
    options.Events = @(t,y) Voc_event(t,y,-1);
    [~,~,~,ye,~] = ode15s(@(t,u) RHS(t,u,@(t) -t,params,vectors,matrices), ...
        [0 75],sol_init,options);
    if isempty(ye), error('No open-circuit events found.'); end
    sol_init = ye(end,:)';
end    

% Reset error message to warning
warning(warnon);

% Define the settings for the call to fsolve
fsoptions = optimoptions('fsolve','MaxIterations',20);
if Verbose, fsoptions.Display = 'iter'; else, fsoptions.Display = 'off'; end
fsoptions.JacobPattern = Jac(params,'findVoc');

% Use the initial estimate to obtain a steady-state solution
sol_init = fsolve(@(u) RHS(0,u,@(t) 0,params,vectors,matrices,'findVoc'), ...
    sol_init,fsoptions);

% Ensure all the algebraic equations are satisfied as exactly as possible
sol_init = apply_Poisson(sol_init,params,vectors,matrices);

% Extract and output open-circuit voltage from value of potential on left
psioc = sol_init(4*N+5);
fprintf('Found Voc to be %0.5g V \n\n', params.psi2Vap(psioc));

% Adjust the voltage protocol if one exists
if ~isnan(psi(time(end)))
    old_psi = @(t) psi(t);
    T1 = splits(2); psi1 = old_psi(T1);
    psi = @(t) old_psi(t)+psioc.*(t<=T1).*(1-old_psi(t)/psi1);
    if Verbose
        clf(99); figure(99);
        plot(tstar2t(time),psi2Vap(psi(time)));
        xlabel('Time (s)'); ylabel('Applied Voltage (V)');
        title('V(t)');
        drawnow;
    end
end

end