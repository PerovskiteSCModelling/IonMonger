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
J16 = spalloc(N+1,NE+1,1);
J17 = spalloc(N+1,NH,1);
J18 = spalloc(N+1,NH+1,1);
J25 = J15; J25(1,NE) = 1;
J26 = J16; J26(1,NE:NE+1) = ones(1,2);
J27 = J17; J27(N+1,1) = 1;
J28 = J18; J28(N+1, 1:2) = ones(1,2);
J38 = spalloc(N+1,NH+1,1); J38(N+1,1) = 1;
J46 = spalloc(N+1,NE+1,1); J46(1,NE+1) = 1;
J51 = spalloc(NE,N+1,1);
J52 = J51; J52(NE,1) = 1;
J55 = gallery('tridiag',NE,1,1,1); J55(1,2) = 0;
J56 = gallery('tridiag',NE+1,1,1,1); J56=J56(1:NE,:); J56(1,:)=zeros(1,NE+1);
J57 = spalloc(NE,NH,1);
J58 = spalloc(NE,NH+1,1);
J61 = spalloc(NE+1,N+1,1);
J62 = J61; J62(NE,1) = 1;
J63 = J61; J63(NE+1,1) = 1;
J65 = gallery('tridiag',NE+1,1,1,1); J65=J65(:,1:NE); J65([1,NE+1],:) = zeros(2,NE);
J66 = gallery('tridiag',NE+1,1,1,1);
J67 = spalloc(NE+1,NH,1);
J68 = spalloc(NE+1,NH+1,1);
J71 = spalloc(NH,N+1,1);
J72 = J71; J72(1,N+1) = 1;
J75 = spalloc(NH,NE,1);
J76 = spalloc(NH,NE+1,1);
J77 = gallery('tridiag',NH,1,1,1); J77(NH,NH-1) = 0;
J78 = gallery('tridiag',NH+1,1,1,1); J78=J78(2:NH+1,:); J78(NH,:)=zeros(1,NH+1);

J81 = spalloc(NH+1,N+1,1);
J82 = J81; J82(2,N+1)=1;
J84 = J81; J84(1,N+1)=1;
J85 = spalloc(NH+1, NE,1);
J86 = spalloc(NH+1, NE+1,1);
J87 = gallery('tridiag',NH+1,1,1,1); J87=J87(:,2:NH+1); J87([1,NH+1],:)=zeros(2,NH);
J88 = gallery('tridiag',NH+1,1,1,1);

% Combine the blocks to define the sparse matrix
JJJ = sparse([  J11 J11 J13 J13 J15 J16 J17 J18; ... % P equation
                J11 J11 J11 J11 J25 J26 J27 J28; ... % phi equation
                J11 J11 J11 J11 J25 J26 J17 J38; ... % n equation
                J11 J11 J11 J11 J15 J46 J27 J28; ... % p equation
                J51 J52 J51 J51 J55 J56 J57 J58; ... % phiE equation
                J61 J62 J63 J61 J65 J66 J67 J68; ... % nE equation
                J71 J72 J71 J71 J75 J76 J77 J78; ... % phiH equation
                J81 J82 J81 J84 J85 J86 J87 J88; ... % pH equation
                ]);

% Adjust right-hand potential BC to account for any parasitic resistance
JJJ(4*N+2*NE+NH+5,:) = [zeros(1,4*N+4), zeros(1,2*NE+1), ...
                        zeros(1,NH-2), 1, 1, zeros(1,NH-1), 1, 1];

% Perform any additional step requested by the optional input argument flag
if nargin>1
    if strcmp(flag,'none')
        % Do nothing else
    elseif strcmp(flag,'open-circuit')
        % The following changes are required to model open-circuit conditions
        JJJ(4*N+5,:) = [zeros(1,4*N+4), 1, 1, zeros(1,NE-2), ...
            1, 1, zeros(1,NE-1), zeros(1,2*NH+1)];
        JJJ(4*N+2*NE+NH+5,:) = [zeros(1,4*N+4), 1, zeros(1,NE-1), ...
            zeros(1,NE+1), zeros(1,NH-1), 1, zeros(1,NH+1)];
    elseif strcmp(flag,'init')
        % The following change is required for conservation of ion vacancies
        JJJ(N+1,:) = [ones(1,N+1), zeros(1,3*N+2*NE+2*NH+5)];
    elseif strcmp(flag,'findVoc')
        % The following changes are required to model open-circuit conditions
        JJJ(4*N+5,:) = [zeros(1,4*N+4), 1, 1, zeros(1,NE-2), ...
            1, 1, zeros(1,NE-1), zeros(1,2*NH+1)];
        JJJ(4*N+2*NE+NH+5,:) = [zeros(1,4*N+4), 1, zeros(1,NE-1), ...
            zeros(1,NE+1), zeros(1,NH-1), 1, zeros(1,NH+1)];
        % The following change is required for conservation of ion vacancies
        JJJ(N+1,:) = [ones(1,N+1), zeros(1,3*N+2*NE+2*NH+5)];
    else
        error('The optional input argument is not recognised.');
    end
end

end
