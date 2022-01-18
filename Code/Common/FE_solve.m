function dstrbns = FE_solve(params,vectors)
% This function carries out a call to the built-in Matlab solver ode15s and
% returns the solution in dimensionless form. The inputs are the structures
% containing the necessary input parameters and vectors. The output is
% another structure containing the dimensionless solution variables.
% The scripts called by this function are located in the Solver folder.

% Parameter input
[psi, time, rtol, atol, Verbose, findVoc, splits, ...
    OutputFcn, MaxStep, Stats] ... % the options on this line may be []
    = struct2array(params,{'psi','time','rtol','atol','Verbose', ...
    'findVoc','splits','OutputFcn','MaxStep','Stats'});

% Program settings
options = odeset('RelTol',rtol,'AbsTol',atol);
options.Vectorized = 'on';
if isfield(params, 'OutputFcn')
    eval(sprintf('OutputFcn = @(t,y,f) %s(t,y,f,params);',OutputFcn));
end
if Verbose, options.OutputFcn = OutputFcn; end
options.MaxStep = MaxStep;
if Verbose, options.Stats = Stats; else, options.Stats = 'off'; end
options.MassSingular = 'yes';
options.MStateDependence = 'none';

% Select 'flag' according to user-defined voltage protocol
if isnan(psi(time(end))) % if psi ends with NaN, this means open-circuit
    flag = 'open-circuit';
else % the dimensionless applied voltage is a specified function of time
    flag = 'none';
end

% Create matrices for RHS
matrices = create_matrices(params,vectors);

if isfield(params,'input_filename')
    if Verbose, disp(['Using initial distributions from the saved file ' ... 
            params.input_filename]), end
    load(params.input_filename)
    
    % unpack final step of solution
    P = sol.dstrbns.P(end,:);
    phi = sol.dstrbns.phi(end,:);
    n = sol.dstrbns.n(end,:);
    p = sol.dstrbns.p(end,:);
    phiE = sol.dstrbns.phiE(end,:);
    nE = sol.dstrbns.nE(end,:);
    phiH = sol.dstrbns.phiH(end,:);
    pH = sol.dstrbns.pH(end,:);
    
    % interpolate onto new spacial grid and nondimensionalise
    [b,N0,VT,dE,dH,kE,kH] = struct2array(params,{'b','N0','VT','dE','dH',...
        'kE','kH'});
    P = interp1(sol.vectors.x,P,vectors.x*b*1e9)/N0;
    phi = interp1(sol.vectors.x,phi,vectors.x*b*1e9)/VT;
    n = interp1(sol.vectors.x,n,vectors.x*b*1e9)/(kE*dE);
    p = interp1(sol.vectors.x,p,vectors.x*b*1e9)/(kH*dH);
    phiE = interp1(sol.vectors.xE,phiE,vectors.xE*b*1e9)/VT;
    nE = interp1(sol.vectors.xE,nE,vectors.xE*b*1e9)/dE;
    phiH = interp1(sol.vectors.xH,phiH,vectors.xH*b*1e9)/VT;
    pH = interp1(sol.vectors.xH,pH,vectors.xH*b*1e9)/dH;
    
    % eliminate superfluous phi points
    phiE = phiE(1:end-1);
    phiH = phiH(2:end);
    
    sol_start = [P ; phi ; n ; p ; phiE ; nE ; phiH ; pH];
    
    if any(isnan(sol_start)); error(['There was an error in interpolating ' ... 
        'the saved distributions onto the new spacial grid'])
    end
else
    % Compute consistent initial conditions for a cell preconditioned at Vbi
    sol_start = initial_conditions(@(t) 0, params,vectors,matrices);
end

% Initiate cycling through different attempts until solution found
[count, err_count] = deal(0);
while count == err_count
    try
        % Always start the solution procedure from steady state at Vbi
        sol_init = sol_start;
        % Perform a preconditioning step if requested
        if isfield(params,'input_filename')
            % do nothing
        elseif findVoc
            % Precondition the cell at open-circuit
            [psi, sol_init] = find_Voc(sol_init,psi,params,vectors,matrices,options);
        elseif abs(psi(0))>atol
            % Precondition at an applied voltage other than Vbi
            sol_init = precondition(sol_init,params,vectors,matrices,options);
        end
        
        % Compute the Jacobian, mass matrix and initial slope and add to options
        if ~exist('AnJac.m','file')
            options.Jacobian = @(t,u) AnJac(t,u,params,vectors,matrices,flag);
        else
        	options.JPattern = Jac(params,flag);
        end
        options.Mass = mass_matrix(params,vectors,flag);
        options.InitialSlope = RHS(0,sol_init,psi,params,vectors,matrices,flag) ...
                                    \options.Mass;

        %% SOLVE

        % Note that to solve without splits and without the loop, it is simply:
        % [time,numsol] = ode15s(@(t,u) RHS(t,u,psi,params,vectors,matrices,flag), ...
        %     time,sol_init,options);

        % Preallocate and initialise index
        numsol = NaN(length(time),length(sol_init));
        next = 1;

        % Make a call to ode15s for each split and iteratively construct solution
        for i = 1:length(splits)-1
            parttime = [splits(i), time(splits(i)<time & time<splits(i+1)), splits(i+1)];
            [~,partsol] = ode15s(@(t,u) RHS(t,u,psi,params,vectors,matrices,flag), ...
                parttime,sol_init,options);
            lts = length(parttime)-2+any(time==splits(i));
            numsol(next:next+lts,:) = partsol(2-any(time==splits(i)):end,:);
            if Verbose
                fprintf('Split %d/%d complete \n',i,length(splits)-1);
            end
            next = next+lts; sol_init = partsol(end,:);

            % Pass the final gradient dudt to the next call
            options.InitialSlope = ...
                RHS(parttime(end),partsol(end,:)',psi,params,vectors,matrices,flag) ...
                \options.Mass;
        end

    catch ME % if no solution can be obtained, try again up to 2 more times
        err_count = err_count + 1;
        if err_count==3
            disp('Could not obtain solution');
            rethrow(ME); % outputs last error message only
        else
            if Verbose
                disp(['Trying again with stricter tolerances, attempt ' ...
                    num2str(err_count)+1]);
            end

            % Reduce the temporal tolerances for the next attempt
            options.RelTol = options.RelTol/10;
            options.AbsTol = options.AbsTol/10;
        end
    end

    count = count + 1;
end

% Unpack numsol into structure
dstrbns = unpack(numsol,params);

end
