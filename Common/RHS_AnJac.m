function [F, J] = RHS_AnJac(u,psi0,params,vectors,matrices,flag)
% Function that returns both the objective function values and the Jacobian
% as outputs, for use as the function in a call to fsolve.

if nargin==5
    flag = 'none';
end
    
F = RHS(0,u,psi0,params,vectors,matrices,flag);
J = AnJac(0,u,params,vectors,matrices,flag);

end