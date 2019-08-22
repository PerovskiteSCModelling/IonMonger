function sol_init = initial_conditions(psi0,params,vectors,matrices)
% This function creates a vector containing a steady-state solution to the 
% initial problem. The inputs are structures containing the parameters,
% vectors and matrices needed by the solver.

% Parameter input
[chi, kH, kE, Verbose] = struct2array(params, {'chi','kH','kE','Verbose'});
[x, xE, xH] = struct2array(vectors, {'x','xE','xH'});

% Define uniform profiles for the ion vacancy density and electric potential
P_init    = ones(size(x));
phi_init  = zeros(size(x));

% Compute profiles for the carrier concentrations from quasi-steady BVPs
p_guess = bvpinit(x',@(x) [kH*x+kE*(1-x)/chi; 0*x]);
sol = bvp4c(@(x,p) pode(x,p,params),@(pa,pb) pbcs(pa,pb,params),p_guess);
solx = deval(sol,x'); p_init = solx(1,:)';
n_guess = bvpinit(x',@(x) [chi*kH*x+kE*(1-x); 0*x]);
sol = bvp4c(@(x,p) node(x,p,params),@(pa,pb) nbcs(pa,pb,params),n_guess);
solx = deval(sol,x'); n_init = solx(1,:)';

% Define uniform profiles for the electric potential and carrier
% concentrations across the TLs
phiE_init = zeros(size(xE));
nE_init   = ones(size(xE));
phiH_init = zeros(size(xH));
pH_init   = ones(size(xH));

% Combine the initial conditions into one vector to pass to fsolve
sol_init  = [P_init; phi_init; n_init; p_init; ... % perovskite
             phiE_init; nE_init; ... % electron transport layer
             phiH_init; pH_init]; % hole transport layer
          
% Define the settings for the call to fsolve
fsoptions = optimoptions('fsolve','MaxIterations',20);
if Verbose, fsoptions.Display = 'iter'; else, fsoptions.Display = 'off'; end
fsoptions.JacobPattern = Jac(params,'init');

% Use the initial guess to obtain a steady-state solution
sol_init = fsolve(@(u) RHS(0,u,psi0,params,vectors,matrices,'init'), ...
    sol_init,fsoptions);

% Ensure all the algebraic equations are satisfied as exactly as possible
sol_init = apply_Poisson(sol_init,params,vectors,matrices);

end



% Quasi-steady BVP for the hole concentration
% y(1) = p(x), y(2) = jp(x) = -Kp*dp/dx
function dpdx = pode(x,p,params)
[G, R, Kp, chi] = struct2array(params,{'G','R','Kp','chi'});
dpdx = [-p(2)/Kp; ...
        G(x,0)-R(chi*p(1),p(1),1)];
end
function res = pbcs(pa,pb,params)
[kH, Rl] = struct2array(params,{'kH','Rl'});
res = [pb(1)-kH; ...
       pa(2)+Rl(1,pa(1))];
end

% Quasi-steady BVP for the electron concentration
% y(1) = n(x), y(2) = jn(x) = Kn*dn/dx
function dpdx = node(x,p,params)
[G, R, Kn, chi] = struct2array(params,{'G','R','Kn','chi'});
dpdx = [p(2)/Kn; ...
        -(G(x,0)-R(chi*p(1),p(1),1))];
end
function res = nbcs(pa,pb,params)
[kE, Rr] = struct2array(params,{'kE','Rr'});
res = [pa(1)-kE; ...
       pb(2)-Rr(pa(1),1)];
end