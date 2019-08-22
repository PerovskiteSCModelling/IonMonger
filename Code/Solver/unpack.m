function dstrbns = unpack(numsol,params)
% This function can be used to unpack the raw solution (numsol) into a
% structure with a field for each solution variable. The inputs are the
% ode15s output array and the structure of input parameters.

% Parameter input
[N, NE, NH] = struct2array(params, {'N','NE','NH'});

% Unpack the raw solution into a structure
dstrbns.P   = numsol(:,1:N+1);
dstrbns.phi = numsol(:,N+2:2*N+2);
dstrbns.n   = numsol(:,2*N+3:3*N+3);
dstrbns.p   = numsol(:,3*N+4:4*N+4);
dstrbns.phiE = numsol(:,4*N+5:4*N+NE+4);
dstrbns.nE   = numsol(:,4*N+NE+5:4*N+2*NE+4);
dstrbns.phiH = numsol(:,4*N+2*NE+5:4*N+2*NE+NH+4);
dstrbns.pH   = numsol(:,4*N+2*NE+NH+5:4*N+2*NE+2*NH+4);

end