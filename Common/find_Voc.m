function [psi, sol_init] = find_Voc(sol_init,psi,params,vectors,matrices,options)
% This function can be used to find the steady-state open-circuit voltage
% (Voc), compute consistent initial conditions and adjust the voltage
% protocol such that it starts from Voc. The inputs are a vector of the
% original initial conditions, structures for the parameters, vectors and
% matrices required by the solvers and the options structure for ode15s.
% The outputs are the new voltage protocol and corresponding vector of
% initial conditions for a cell preconditioned at open-circuit.

% Parameter input
[N, Kn, Kp, G, time, splits, t2tstar, tstar2t, Vap2psi, psi2Vap, dpf, ...
    pbi, ARp, Verbose, phidisp] ...
    = struct2array(params, {'N','Kn','Kp','G','time','splits','t2tstar', ...
                            'tstar2t','Vap2psi','psi2Vap','dpf','pbi', ...
                            'ARp','Verbose', 'phidisp'});
dx = vectors.dx;

% Temporarily turn the following warning into an error
warnon = warning('error','MATLAB:ode15s:IntegrationTolNotMet');

% Compute the Jacobian, mass matrix and initial slope and add to options
% Note that the mass matrix is adjusted so ions are as mobile as electrons
if exist('AnJac.m','file')
    options.Jacobian = @(t,u) AnJac(t,u,params,vectors,matrices);
else
    options.JPattern = Jac(params);
end
options.Mass = mass_matrix(params,vectors,'precondition');
options.InitialSlope = RHS(0,sol_init,@(t) 0,params,vectors,matrices) ...
                        \options.Mass;

% Fix the light in its initial state
params.G = @(x,t) G(x,0);

% Define the total current density (nelecting the displacement current) and
% the Event function to find when the current is zero at...
mid = ceil((N+1)/2); % the midpoint of the perovskite layer
J = @(mid,y) Kn./dx(mid).*(y(2*N+3+mid)-y(2*N+2+mid) ...
                -(y(2*N+3+mid)+y(2*N+2+mid)).*(y(N+2+mid)-y(N+1+mid))./2) ...
            -Kp./dx(mid).*(y(3*N+4+mid)-y(3*N+3+mid) ...
                +(y(3*N+4+mid)+y(3*N+3+mid)).*(y(N+2+mid)-y(N+1+mid))./2) ...
            +dpf./dx(mid).*(y(1+mid)-y(mid)+(y(1+mid)+y(mid)) ...
                .*(y(N+2+mid)-y(N+1+mid))./2) ...
            -(pbi-2*y(4*N+5))/ARp;
function [value,isterminal,direction] = Voc_event(t,y,direction)
    value = J(mid,y); % event is when current density equals zero
    isterminal = 1; % stop at open-circuit
    % direction = -1 for an increasing voltage and 1 for decreasing
end

% Precondition the cell at 0V
search_time = t2tstar(10);
[~,scan_down] = ode15s(@(t,u) RHS(t,u,@(t) ... Vap2psi(0)*t/search_time, ...
                    Vap2psi(0)*(2-t/search_time)*t/search_time, ...
                    params,vectors,matrices),[0 search_time],sol_init,options);
sol_init = scan_down(end,:)';

if integral(@(x) params.G(x,0),0,1) == 0
    % If there is no illumination, assume that open-circuit is at 0V
    
else
    % Try to locate the Voc between 0 and 2V by slowly increasing the voltage
    if Verbose
        disp('Attempting to find Voc between 0 and 2V');
    end
    options.Events = @(t,y) Voc_event(t,y,-1);
    options.InitialSlope = RHS(0,sol_init,@(t) Vap2psi(0), ...
                            params,vectors,matrices)\options.Mass;
    [~,y,~,ye,~] = ode15s(@(t,u) RHS(t,u,@(t) ...
                    (Vap2psi(0)*(search_time-t)+Vap2psi(2)*t)/search_time, ...
                    params,vectors,matrices),[0 search_time],sol_init,options);
    if isempty(ye)
        error('No open-circuit events found.');
    else
        sol_init = ye(end,:)';
    end
    
    % Try to find a better estimate with stricter tolerances
    try
        if Verbose
            disp('Attempting to find better estimate near the first estimate');
        end
        options.AbsTol = 1e-14;
        options.InitialSlope = [];
        options.OutputFcn = [];
        [~,~,~,ze,~] = ode15s(@(t,u) RHS(t,u,@(t) y(end-1,4*N+5)-phidisp-10*t/search_time, ...
                        params,vectors,matrices),[0 search_time],y(end-1,:)',options);
        if isempty(ze)
            if Verbose
                disp('No better estimate found, continung with original estimate');
            end
        else
            if Verbose
                disp('Better estimate found');
            end
            sol_init = ze(end,:)';
        end
    catch
    end

end            

% Reset error message to warning
warning(warnon);

% Define the settings for the call to fsolve
fsoptions = optimoptions('fsolve','MaxIterations',20);
if Verbose, fsoptions.Display = 'iter'; else, fsoptions.Display = 'off'; end

% Use the initial estimate to obtain an approximate steady-state solution
if exist('AnJac.m','file')
    fsoptions.SpecifyObjectiveGradient = true;
    [sol_init,~,exitflag,~] = fsolve(@(u) RHS_AnJac(u,@(t) 0, ...
        params,vectors,matrices,'findVoc'),sol_init,fsoptions);
else
    fsoptions.JacobPattern = Jac(params,'findVoc');
    [sol_init,~,exitflag,~] = fsolve(@(u) RHS(0,u,@(t) 0, ...
        params,vectors,matrices,'findVoc'),sol_init,fsoptions);
end
if exitflag<1
    warning(['Steady-state initial conditions could not be found to ' ...
        'a high degree of accuracy and may be unphysical.']);
end

% Ensure all the algebraic equations are satisfied as exactly as possible
sol_init = apply_Poisson(sol_init,params,vectors,matrices);

% Extract and output open-circuit voltage from value of potential on left
psioc = sol_init(4*N+5)-phidisp;
fprintf('Found estimate for Voc of %0.5g V \n\n', params.psi2Vap(psioc));

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