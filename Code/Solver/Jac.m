function JJJ = Jac(params,flag)
% This function defines a sparse matrix with zeros where there are always
% zeros in the Jacobian of the RHS and ones indicating entries that need
% to be calculated numerically during the solution process. This matrix is
% passed to ode15s via the JPattern option. The input is a structure
% containing the necessary parameters. The last input, flag, is an optional
% argument which is used to adapt the function, e.g. to simulate
% open-circuit conditions.

% Parameter input
[N, NE, NH] = struct2array(params,{'N','NE','NH'});

% Define the sparse matrix block by block
J11 = gallery('tridiag',N+1,1,1,1);
J13 = spalloc(N+1,N+1,1);
J15 = spalloc(N+1,NE,1);
J17 = spalloc(N+1,NH,1);
J25 = J15; J25(1,NE) = 1;
J27 = J17; J27(N+1,1) = 1;
J51 = spalloc(NE,N+1,1);
J52 = J51; J52(NE,1) = 1;
J55 = gallery('tridiag',NE,1,1,1); J55(1,2) = 0;
J56 = J55; J56(1,1) = 0;
J57 = spalloc(NE,NH,1);
J71 = spalloc(NH,N+1,1);
J72 = J71; J72(1,N+1) = 1;
J75 = spalloc(NH,NE,1);
J77 = gallery('tridiag',NH,1,1,1); J77(NH,NH-1) = 0;
J78 = J77; J78(NH,NH) = 0;

% Combine the blocks to define the sparse matrix
JJJ = sparse([  J11 J11 J13 J13 J15 J15 J17 J17; ... % P equation
                J11 J11 J11 J11 J25 J25 J27 J27; ... % phi equation
                J11 J11 J11 J11 J25 J25 J17 J17; ... % n equation
                J11 J11 J11 J11 J15 J15 J27 J27; ... % p equation
                J51 J52 J52 J51 J55 J56 J57 J57; ... % phiE equation
                J51 J52 J52 J51 J56 J55 J57 J57; ... % nE equation
                J71 J72 J71 J72 J75 J75 J77 J78; ... % phiH equation
                J71 J72 J71 J72 J75 J75 J78 J77; ... % pH equation
                ]);

% Adjust right-hand potential BC to account for any parasitic resistance
JJJ(4*N+2*NE+NH+4,:) = [zeros(1,4*N+4), zeros(1,2*NE), ...
                        zeros(1,NH-2), 1, 1, zeros(1,NH-2), 1, 1];

% Perform any additional step requested by the optional input argument flag
if nargin>1
    if strcmp(flag,'none')
        % Do nothing else
    elseif strcmp(flag,'open-circuit')
        % The following changes are required to model open-circuit conditions
        JJJ(4*N+5,:) = [zeros(1,4*N+4), 1, 1, zeros(1,NE-2), ...
            1, 1, zeros(1,NE-2), zeros(1,2*NH)];
        JJJ(4*N+2*NE+NH+4,:) = [zeros(1,4*N+4), 1, zeros(1,NE-1), ...
            zeros(1,NE), zeros(1,NH-1), 1, zeros(1,NH)];
    elseif strcmp(flag,'init')
        % The following change is required for conservation of ion vacancies
        JJJ(N+1,:) = [ones(1,N+1), zeros(1,3*N+2*NE+2*NH+3)];
    elseif strcmp(flag,'findVoc')
        % The following changes are required to model open-circuit conditions
        JJJ(4*N+5,:) = [zeros(1,4*N+4), 1, 1, zeros(1,NE-2), ...
            1, 1, zeros(1,NE-2), zeros(1,2*NH)];
        JJJ(4*N+2*NE+NH+4,:) = [zeros(1,4*N+4), 1, zeros(1,NE-1), ...
            zeros(1,NE), zeros(1,NH-1), 1, zeros(1,NH)];
        % The following change is required for conservation of ion vacancies
        JJJ(N+1,:) = [ones(1,N+1), zeros(1,3*N+2*NE+2*NH+3)];
    else
        error('The optional input argument is not recognised.');
    end
end

end