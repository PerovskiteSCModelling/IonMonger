function solution = numericalsolver(params)
% This is the main function which carries out the solution procedure. The
% input is a params structure (created by parameters.m). This script first
% creates another structure containing vectors needed by the solver
% (including the spatial grid) then calls FE_solve.m to obtain a solution,
% re-dimensionalises the solution variables and finally outputs a structure
% containing both the output and the input.


%% Carry out the solution procedure

% Start stopwatch
start_t = tic;

% Parse necessary parameters
if ~isfield(params, 'Verbose'), params.Verbose=true; end
[time, Verbose] ...
    = struct2array(params,{'time','Verbose'});

% Create vectors to parse to the solver (includes the spatial setup)
vectors = create_vectors(params);

% Obtain a solution using a method based on finite elements
if Verbose, disp('Calculating numerical solutions...'); end
dstrbns = FE_solve(params,vectors);

% Calculate the current density and interfacial losses
[J, Jl, Jr] = calculate_currents(params, vectors, dstrbns);


%% Re-dimensionalise all outputs
% including vectors, solution variables, time, voltage and current densities
[b, N0, VT, dE, dH, tstar2t, psi2Vap, Vbi, jay] ...
    = struct2array(params,{'b','N0','VT','dE','dH','tstar2t','psi2Vap', ...
                           'Vbi','jay'});
b = b*1e9; % Change the width b in metres to nanometres
vectors.x    = b*vectors.x;
vectors.dx   = b*vectors.dx;
vectors.xE   = b*vectors.xE;
vectors.dxE  = b*vectors.dxE;
vectors.xH   = b*vectors.xH;
vectors.dxH  = b*vectors.dxH;
dstrbns.P    = N0*dstrbns.P;
dstrbns.phi  = VT*dstrbns.phi;
dstrbns.n    = dE*dstrbns.n;
dstrbns.p    = dH*dstrbns.p;
dstrbns.phiE = VT*dstrbns.phiE;
dstrbns.nE   = dE*dstrbns.nE;
dstrbns.phiH = VT*dstrbns.phiH;
dstrbns.pH   = dH*dstrbns.pH;
time         = tstar2t(time);
V            = psi2Vap(dstrbns.phiE(:,1)/VT);
Vres         = V-Vbi+dstrbns.phiE(:,1)-dstrbns.phiH(:,end);
J            = jay.*J;
Jl           = jay.*Jl;
Jr           = jay.*Jr;


%% Final steps

% Stop stopwatch
timetaken = toc(start_t);

% Package up solution
solution = struct('vectors',vectors, 'params',params, 'dstrbns',dstrbns, ...
        'time',time, 'V',V, 'Vres',Vres, 'J',J, 'Jl',Jl, 'Jr',Jr, ...
        'timetaken',timetaken);

% If completion tasks are specified, perform them
if exist('completion_tasks','file'), completion_tasks(solution); end

end
