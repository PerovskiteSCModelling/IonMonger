function M = mass_matrix(params,vectors,flag)
% This function constructs the mass matrix M for the equation M.du/dt=f(u).
% M is a diagonal matrix, where the values on the diagonal are the
% prefactors of the time derivative in each DAE and zero when there is no
% time derivative. As such, Dirichlet BCs and the equations for the
% electric potential correspond to a zero in the mass matrix. The inputs
% are structures containing the necessary parameters and vectors. The last
% input, flag, is an optional argument which is used to adapt the function
% e.g. to simulate open-circuit conditions.

% Parameter input
[sigma, chi, N, NE, NH, kE, kH] = ...
    struct2array(params,{'sigma','chi','N','NE','NH','kE','kH'});
[dx, dxE, dxH] = struct2array(vectors,{'dx','dxE','dxH'});

% Define the mass matrix block by block
M11 = gallery('tridiag',N+1,[dx(1:end-1)/6; dx(end)/6], ...
    [dx(1)/3; (dx(1:end-1)+dx(2:end))/3; dx(end)/3],[dx(1)/6; dx(2:end)/6]);
M12 = spalloc(N+1,N+1,1);
M15 = spalloc(N+1,NE,1);
M17 = spalloc(N+1,NH,1);
M33 = sigma*M11; M33(1,1:2) = kE*M33(1,1:2)+[sigma*dxE(end)/3,0];
M36 = M15; M36(1,end) = sigma*dxE(end)/6;
M44 = sigma*chi*M11; M44(end,N:end) = kH*M44(end,N:end)+[0,sigma*chi*dxH(1)/3];
M48 = M17; M48(end,1) = sigma*chi*dxH(1)/6;
M51 = spalloc(NE,N+1,1);
M55 = spalloc(NE,NE,1);
M57 = spalloc(NE,NH,1);
M63 = M51; M63(end,1) = sigma*dxE(end)/6;
M66 = sigma*gallery('tridiag',NE,dxE(1:end-1)/6, ...
    [0; (dxE(1:end-1)+dxE(2:end))/3],[0; dxE(2:end-1)/6]);
M71 = spalloc(NH,N+1,1);
M75 = spalloc(NH,NE,1);
M77 = spalloc(NH,NH,1);
M84 = M71; M84(1,end) = sigma*chi*dxH(1)/6;
M88 = sigma*chi*gallery('tridiag',NH,[dxH(2:end-1)/6; 0], ...
    [(dxH(1:end-1)+dxH(2:end))/3; 0],dxH(2:end)/6);

% Combine the blocks to define the mass matrix
M = sparse([ ...
    M11 M12 M12 M12 M15 M15 M17 M17; ... % P equation
    M12 M12 M12 M12 M15 M15 M17 M17; ... % phi equation
    M12 M12 M33 M12 M15 M36 M17 M17; ... % n equation
    M12 M12 M12 M44 M15 M15 M17 M48; ... % p equation
    M51 M51 M51 M51 M55 M55 M57 M57; ... % phiE equation
    M51 M51 M63 M51 M55 M66 M57 M57; ... % nE equation
    M71 M71 M71 M71 M75 M75 M77 M77; ... % phiH equation
    M71 M71 M71 M84 M75 M75 M77 M88; ... % pH equation
    ]);

% Perform any additional step requested by the optional input argument flag
if nargin>2
    if strcmp(flag,'none')
        % Do nothing else
    elseif strcmp(flag,'open-circuit')
        % The following changes are required to model open-circuit conditions
        M(4*N+5,4*N+NE+5:4*N+NE+6) = sigma*dxE(1)*[1/3, 1/6];
    elseif strcmp(flag,'precondition')
        % The following changes artificially increase the ion mobility to findVoc
        M(1:N+1,:) = sigma*M(1:N+1,:);
    else
        error('The optional input argument is not recognised.');
    end
end

end